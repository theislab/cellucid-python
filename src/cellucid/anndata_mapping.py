"""
The data-mapping routes of an AnnData server: every array, read on request.

A viewer connected to ``cellucid serve`` or a notebook can re-map the dataset
it shows -- another embedding or other columns of one, another neighbour graph,
gene names from a ``var`` column, expression from a layer or ``raw``, another
latent space, vectors bound to another embedding -- without the server knowing
in advance. This module answers those requests under ``_cellucid/mapping/``:

==========================================================  ===========================
route (relative to the dataset base URL)                    payload
==========================================================  ===========================
``candidates.json``                                         the candidate catalog
``embeddings/{i}/{c0}-{c1}[-{c2}].bin``                     normalized float32, n x d
``vectors/{i}/{vcols}/{e}/{ecols}.bin``                     float32, scaled by embedding
``knn/{p}/connectivity_manifest.json``                      edge-pair manifest
``knn/{p}/connectivity/edges.{src,dst,weights.f64}.bin``    edge payloads
``genes/{axis}/{names}/var_manifest.json``                  compact var manifest
``genes/{axis}/{names}/var/{g}.values.f32``                 one gene, float32
``latent/{i}/obs/{f}.outliers.f32``                         outlier quantiles
==========================================================  ===========================

``{i}``/``{e}`` index the catalog's ``obsm`` list, ``{p}`` its ``obsp`` list,
``{axis}`` is ``X``, ``raw`` or ``layer-{k}`` (``k`` into ``layers``),
``{names}`` is ``index`` or ``column-{k}`` (``k`` into that axis's var columns),
``{f}`` is an obs payload index and ``{g}`` a gene's position on the axis.
Routes address entries by position, never by key, so no URL ever carries a
dataset's vocabulary. Columns are 0-based.

An index the catalog does not have is :class:`MappingNotFound` (HTTP 404). An
entry the catalog lists as unusable, or one whose values fail the load-time
check -- a NaN in an embedding, a graph that is not symmetric, a latent space
over the memory budget -- is :class:`MappingLoadError` (HTTP 422 with
``{"error", "key", "reason"}``). Nothing is validated before it is requested, and no
check is ever replaced by a repair: a non-symmetric matrix is refused, never
symmetrized.

Payloads are built lazily and kept in one byte-bounded LRU. The service is
safe to call from several threads: the cache, the per-entry normalization
factors and the sparse column copies are all guarded by one lock, and building
a payload happens outside it.
"""

from __future__ import annotations

import json
import math
import re
import threading
from collections import OrderedDict
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, cast

import numpy as np
from scipy import sparse

from ._byte_order import little_endian_payload_bytes
from .connectivity_contract import (
    ConnectivityEdgePairs,
    build_connectivity_manifest,
    validate_connectivity_edges,
)
from .continuous_payload_diagnosis import NonFinitePayloadError, diagnose_continuous_payload
from .data_mapping import (
    DataMappingError,
    UniqueNames,
    build_anndata_candidates,
    column_gene_names,
    validate_embedding_columns,
)
from .prepare_data._arrays import (
    _normalize_finite_float32_embedding,
    _require_finite_embedding_source,
    _require_finite_float32_array,
)
from .prepare_data._categories import _json_category_values
from .prepare_data._centroids import _compute_latent_space_quantiles
from .prepare_data._manifest import MANIFEST_FORMAT_VERSION
from .vector_fields import _finite_float32_matrix, scale_vector_field

if TYPE_CHECKING:
    from .anndata_adapter import AnnDataAdapter

#: The prefix every mapping route lives under, relative to the dataset base.
MAPPING_ROUTE_PREFIX = "_cellucid/mapping/"

#: Default bound of the payload cache.
DEFAULT_MAPPING_CACHE_BYTES = 256 * 1024 * 1024

_DECIMAL = r"(?:0|[1-9][0-9]*)"
_INDEX = rf"({_DECIMAL})"
_COLUMNS = rf"({_DECIMAL}(?:-{_DECIMAL}){{0,2}})"
_AXIS = rf"(X|raw|layer-{_DECIMAL})"
_NAMES = rf"(index|column-{_DECIMAL})"
_EDGE_FILES = ("edges.src.bin", "edges.dst.bin", "edges.weights.f64.bin")


class MappingNotFound(LookupError):
    """A mapping route that names no catalog entry (HTTP 404)."""


class MappingLoadError(ValueError):
    """A catalog entry that cannot be served as asked (HTTP 422).

    ``key`` is the AnnData key the refusal is about -- an ``obsm`` or ``obsp``
    key, ``"X"``, ``"raw"``, a layer name or a ``var`` column -- so a viewer can
    mark that one option unusable and keep the rest. ``reason`` is the clause a
    viewer lists beside that key ("NaN at cell 17", "not symmetric"), the same
    words a file read in the browser gives; the message is the full sentence.
    """

    def __init__(self, message: str, *, key: str, reason: str) -> None:
        super().__init__(message)
        self.key = key
        self.reason = reason

    def as_payload(self) -> dict[str, str]:
        return {"error": str(self), "key": self.key, "reason": self.reason}


@dataclass(frozen=True)
class MappingRoute:
    """One parsed mapping route."""

    kind: str
    parts: tuple[Any, ...]


def _columns(text: str) -> tuple[int, ...]:
    return tuple(int(part) for part in text.split("-"))


def parse_mapping_route(path: str) -> MappingRoute | None:
    """Parse the part of a path after ``_cellucid/mapping/``, or return None.

    Only the grammar is checked here; whether an index exists is the
    service's question.
    """
    if path == "candidates.json":
        return MappingRoute("candidates", ())
    match = re.fullmatch(rf"embeddings/{_INDEX}/{_COLUMNS}\.bin", path)
    if match:
        return MappingRoute("embedding", (int(match.group(1)), _columns(match.group(2))))
    match = re.fullmatch(rf"vectors/{_INDEX}/{_COLUMNS}/{_INDEX}/{_COLUMNS}\.bin", path)
    if match:
        return MappingRoute(
            "vectors",
            (
                int(match.group(1)),
                _columns(match.group(2)),
                int(match.group(3)),
                _columns(match.group(4)),
            ),
        )
    match = re.fullmatch(rf"knn/{_INDEX}/connectivity_manifest\.json", path)
    if match:
        return MappingRoute("knn_manifest", (int(match.group(1)),))
    match = re.fullmatch(rf"knn/{_INDEX}/connectivity/([a-z0-9.]+)", path)
    if match and match.group(2) in _EDGE_FILES:
        return MappingRoute("knn_edges", (int(match.group(1)), match.group(2)))
    match = re.fullmatch(rf"genes/{_AXIS}/{_NAMES}/var_manifest\.json", path)
    if match:
        return MappingRoute("var_manifest", (match.group(1), match.group(2)))
    match = re.fullmatch(rf"genes/{_AXIS}/{_NAMES}/var/{_INDEX}\.values\.f32", path)
    if match:
        return MappingRoute("gene", (match.group(1), match.group(2), int(match.group(3))))
    match = re.fullmatch(rf"latent/{_INDEX}/obs/{_INDEX}\.outliers\.f32", path)
    if match:
        return MappingRoute("latent", (int(match.group(1)), int(match.group(2))))
    return None


class _ByteLRU:
    """A least-recently-used cache bounded by the bytes its values hold."""

    def __init__(self, max_bytes: int) -> None:
        self.max_bytes = max_bytes
        self._items: OrderedDict[Any, tuple[Any, int]] = OrderedDict()
        self._bytes = 0

    @property
    def resident_bytes(self) -> int:
        return self._bytes

    def get(self, key: Any) -> Any | None:
        item = self._items.get(key)
        if item is None:
            return None
        self._items.move_to_end(key)
        return item[0]

    def put(self, key: Any, value: Any, size: int) -> None:
        if key in self._items:
            self._bytes -= self._items.pop(key)[1]
        if size > self.max_bytes:
            return
        self._items[key] = (value, size)
        self._bytes += size
        while self._bytes > self.max_bytes:
            _key, (_value, evicted) = self._items.popitem(last=False)
            self._bytes -= evicted

    def clear(self) -> None:
        self._items.clear()
        self._bytes = 0


def _value_defect(values: Any) -> str | None:
    """The first value of a numeric array Cellucid cannot read, as a clause.

    Rows are cells. The words are the browser's ("NaN at cell 17",
    "not representable as float32 at cell 4"), so one array is refused in the
    same terms on every path. Returns None when every value is usable.
    """
    array = np.asarray(values)
    if array.dtype.kind not in {"i", "u", "f"}:
        return "not a numeric array"
    if array.ndim == 1:
        array = array.reshape(-1, 1)
    with np.errstate(over="ignore", invalid="ignore"):
        finite = np.isfinite(array)
        if not finite.all():
            row = int(np.argwhere(~finite)[0][0])
            value = float(array[row][~finite[row]][0])
            kind = "NaN" if math.isnan(value) else ("Infinity" if value > 0 else "-Infinity")
            return f"{kind} at cell {row}"
        narrowed = array.astype(np.float32)
        bad = ~np.isfinite(narrowed) | ((array != 0) & (narrowed == 0))
        if bad.any():
            return f"not representable as float32 at cell {int(np.argwhere(bad)[0][0])}"
    return None


def _size_of(value: Any) -> int:
    """Estimate the bytes one cached value holds."""
    if isinstance(value, bytes | bytearray):
        return len(value)
    if isinstance(value, np.ndarray):
        return int(value.nbytes)
    if isinstance(value, ConnectivityEdgePairs):
        return int(value.sources.nbytes + value.destinations.nbytes + value.weights.nbytes)
    if isinstance(value, UniqueNames):
        return sum(len(name) + 64 for name in value.names) + sum(len(name) for name in value.filled)
    if isinstance(value, tuple | list):
        return sum(_size_of(item) for item in value)
    if isinstance(value, dict):
        return len(json.dumps(value))
    return 64


class AnnDataMapping:
    """Serve the mapping routes of one :class:`~cellucid.AnnDataAdapter`."""

    def __init__(
        self,
        adapter: AnnDataAdapter,
        *,
        max_cache_bytes: int = DEFAULT_MAPPING_CACHE_BYTES,
    ) -> None:
        self._adapter = adapter
        self._lock = threading.Lock()
        self._cache = _ByteLRU(max_cache_bytes)
        self._csc_copies: dict[str, sparse.csc_matrix] = {}
        # Built on first use, then frozen: every route index refers to it, and
        # constructing an adapter must not touch obsp or the layers at all.
        self._frozen_candidates: dict[str, Any] | None = None

    # ------------------------------------------------------------------ catalog

    @property
    def adata(self) -> Any:
        return self._adapter.adata

    @property
    def _candidates(self) -> dict[str, Any]:
        """The frozen candidate catalog every route index refers to."""
        with self._lock:
            if self._frozen_candidates is None:
                self._frozen_candidates = self._build_candidates()
            return self._frozen_candidates

    @property
    def candidates(self) -> dict[str, Any]:
        """The frozen candidate catalog, built on first use."""
        return self._candidates

    @property
    def cache_bytes(self) -> int:
        with self._lock:
            return self._cache.resident_bytes

    def _build_candidates(self) -> dict[str, Any]:
        adapter = self._adapter
        catalog = build_anndata_candidates(
            adapter.adata,
            latent_key=adapter.latent_key,
            knn_key=adapter.knn_key,
        )
        # An array the best guess already loaded and refused is listed with
        # the load-time reason, exactly as the browser lists it, so no viewer
        # offers a choice this server has proven it cannot serve.
        rejected = dict(adapter._guess_rejected)
        catalog["obsm"] = [
            {**item, "usable": False, "reason": rejected[item["key"]]}
            if item["key"] in rejected
            else item
            for item in catalog["obsm"]
        ]
        # The defaults are what this server actually serves, so a viewer that
        # starts from them shows exactly the dataset identity's layout.
        catalog["defaults"] = adapter.get_data_mapping()
        catalog["guessed"] = list(adapter._guessed_slots) + (
            ["expression"] if adapter._expression_guessed else []
        )
        return catalog

    # ---------------------------------------------------------------- lookups

    def _obsm_item(self, index: int) -> dict[str, Any]:
        entries = self._candidates["obsm"]
        if not 0 <= index < len(entries):
            raise MappingNotFound(f"obsm entry {index} does not exist.")
        item = cast(dict[str, Any], entries[index])
        if not item["usable"]:
            raise MappingLoadError(
                f"obsm[{item['key']!r}] {item['reason']}.", key=item["key"], reason=item["reason"]
            )
        return item

    def _obsm_width(self, item: dict[str, Any]) -> int:
        shape = item["shape"]
        return 1 if len(shape) == 1 else int(shape[1])

    def _require_columns(self, item: dict[str, Any], columns: tuple[int, ...]) -> list[int]:
        try:
            return validate_embedding_columns(
                list(columns), dimension=len(columns), width=self._obsm_width(item)
            )
        except DataMappingError as error:
            raise MappingNotFound(f"obsm[{item['key']!r}]: {error}") from error

    def _dense_obsm(self, key: str) -> np.ndarray:
        array = np.asarray(self.adata.obsm[key])
        return array.reshape(-1, 1) if array.ndim == 1 else array

    def _cached(self, key: Any, build: Any) -> Any:
        with self._lock:
            value = self._cache.get(key)
        if value is not None:
            return value
        value = build()
        with self._lock:
            self._cache.put(key, value, _size_of(value))
        return value

    # -------------------------------------------------------------- embeddings

    def _embedding(self, index: int, columns: tuple[int, ...]) -> tuple[np.ndarray, float]:
        """Return one mapped embedding's normalized coordinates and its scale."""
        item = self._obsm_item(index)
        resolved = self._require_columns(item, columns)
        key = item["key"]

        def build() -> tuple[np.ndarray, float]:
            # The load-time check of the data mapping, with the same verdict
            # and text the browser and this server's best guess give ("NaN at
            # cell 17"), so every path refuses one array in one sentence.
            defect = self._adapter._embedding_load_defect(key, resolved)
            if defect is not None:
                raise MappingLoadError(
                    f"Embedding {key!r} cannot be plotted: {defect}.", key=key, reason=defect
                )
            try:
                source = _require_finite_embedding_source(
                    self._dense_obsm(key)[:, resolved],
                    label=f"Embedding {key!r}",
                )
                if not self._adapter.normalize_embeddings:
                    return source, 1.0
                normalized, _center, scale, _range = _normalize_finite_float32_embedding(
                    source, label=f"Embedding {key!r}"
                )
            except (TypeError, ValueError) as error:
                raise MappingLoadError(str(error), key=key, reason=str(error)) from error
            return normalized, float(scale)

        return cast(tuple[np.ndarray, float], self._cached(("embedding", index, tuple(resolved)), build))

    def embedding_bytes(self, index: int, columns: tuple[int, ...]) -> bytes:
        coordinates, _scale = self._embedding(index, columns)
        return cast(
            bytes,
            self._cached(
                ("embedding-bytes", index, columns),
                lambda: little_endian_payload_bytes(coordinates.astype(np.float32)),
            ),
        )

    # ----------------------------------------------------------------- vectors

    def _vector_binding(
        self,
        vector_index: int,
        vector_columns: tuple[int, ...],
        embedding_index: int,
        embedding_columns: tuple[int, ...],
    ) -> tuple[dict[str, Any], list[int], int]:
        """Check that one vector array may be bound to one mapped embedding."""
        vector_item = self._obsm_item(vector_index)
        embedding_item = self._obsm_item(embedding_index)
        resolved_vector = self._require_columns(vector_item, vector_columns)
        self._require_columns(embedding_item, embedding_columns)
        dimension = len(embedding_columns)
        vector_width = self._obsm_width(vector_item)
        same_columns = (
            vector_width == self._obsm_width(embedding_item)
            and tuple(resolved_vector) == embedding_columns
        )
        whole = vector_width == dimension and resolved_vector == list(range(dimension))
        if not (same_columns or whole):
            raise MappingNotFound(
                f"obsm[{vector_item['key']!r}] cannot be bound to obsm[{embedding_item['key']!r}] "
                "at these columns."
            )
        return vector_item, resolved_vector, dimension

    def vector_bytes(
        self,
        vector_index: int,
        vector_columns: tuple[int, ...],
        embedding_index: int,
        embedding_columns: tuple[int, ...],
    ) -> bytes:
        """One vector field read at ``vector_columns``, bound to a mapped embedding.

        The vector array is read at the embedding's own columns when the two
        arrays are equally wide (``velocity_pca`` over ``X_pca``), or whole when
        it is exactly as wide as the embedding's dimension. It is scaled by the
        bound embedding's normalization factor.
        """
        vector_item, resolved_vector, dimension = self._vector_binding(
            vector_index, vector_columns, embedding_index, embedding_columns
        )
        _coordinates, scale = self._embedding(embedding_index, embedding_columns)
        key = vector_item["key"]

        def build() -> bytes:
            selected = self._dense_obsm(key)[:, resolved_vector]
            defect = _value_defect(selected)
            if defect is not None:
                raise MappingLoadError(
                    f"Vector field {key!r} cannot be drawn: {defect}.", key=key, reason=defect
                )
            try:
                values, _dimension = _finite_float32_matrix(
                    self._dense_obsm(key)[:, resolved_vector],
                    label=f"Vector field {key!r}",
                    n_rows=int(self.adata.n_obs),
                    declared_dimension=dimension,
                )
                scaled = scale_vector_field(values, scale_factor=scale, label=f"Vector field {key!r}")
            except (TypeError, ValueError) as error:
                raise MappingLoadError(str(error), key=key, reason=str(error)) from error
            return little_endian_payload_bytes(scaled)

        return cast(
            bytes,
            self._cached(
                ("vectors", vector_index, vector_columns, embedding_index, embedding_columns), build
            ),
        )

    # --------------------------------------------------------------------- knn

    def _obsp_item(self, index: int) -> dict[str, Any]:
        entries = self._candidates["obsp"]
        if not 0 <= index < len(entries):
            raise MappingNotFound(f"obsp entry {index} does not exist.")
        item = cast(dict[str, Any], entries[index])
        if not item["usable"]:
            raise MappingLoadError(
                f"obsp[{item['key']!r}] is {item['reason']}.", key=item["key"], reason=item["reason"]
            )
        return item

    def _edges(self, index: int) -> ConnectivityEdgePairs:
        item = self._obsp_item(index)
        key = item["key"]

        def build() -> ConnectivityEdgePairs:
            # Checked when requested, never in advance and never repaired.
            try:
                return validate_connectivity_edges(self.adata.obsp[key], n_cells=int(self.adata.n_obs))
            except (TypeError, ValueError) as error:
                reason = str(error)
                if "symmetric" in reason:
                    raise MappingLoadError(
                        f"obsp[{key!r}] is not symmetric — choose a connectivities matrix.",
                        key=key,
                        reason="not symmetric — choose a connectivities matrix",
                    ) from error
                raise MappingLoadError(f"obsp[{key!r}]: {reason}", key=key, reason=reason) from error

        return cast(ConnectivityEdgePairs, self._cached(("knn", index), build))

    def connectivity_manifest(self, index: int) -> dict[str, Any]:
        edges = self._edges(index)
        return build_connectivity_manifest(
            n_cells=int(self.adata.n_obs),
            n_edges=edges.n_edges,
            max_neighbors=edges.max_neighbors,
            index_bytes=edges.index_bytes,
            index_dtype=edges.index_dtype,
            compression=None,
        )

    def edge_bytes(self, index: int, filename: str) -> bytes:
        edges = self._edges(index)
        array = {
            "edges.src.bin": edges.sources,
            "edges.dst.bin": edges.destinations,
            "edges.weights.f64.bin": edges.weights,
        }[filename]
        return little_endian_payload_bytes(array)

    def edge_length(self, index: int, filename: str) -> int:
        edges = self._edges(index)
        width = 8 if filename == "edges.weights.f64.bin" else edges.index_bytes
        return edges.n_edges * width

    # ------------------------------------------------------------------- genes

    def _axis(self, axis: str) -> tuple[str, Any, Any, list[dict[str, Any]]]:
        """Resolve ``X``, ``raw`` or ``layer-{k}``: (key, matrix, var frame, var columns)."""
        if axis == "X":
            if self._candidates["X"] is None:
                raise MappingNotFound("This dataset has no X.")
            return "X", self.adata.X, self.adata.var, self._candidates["var"]["columns"]
        if axis == "raw":
            if self._candidates["raw"] is None:
                raise MappingNotFound("This dataset has no raw.")
            raw = self.adata.raw
            return "raw", raw.X, raw.var, self._candidates["raw"]["var"]["columns"]
        index = int(axis.removeprefix("layer-"))
        layers = self._candidates["layers"]
        if not 0 <= index < len(layers):
            raise MappingNotFound(f"layer {index} does not exist.")
        layer = layers[index]
        if not layer["usable"]:
            raise MappingLoadError(
                f"layers[{layer['key']!r}] {layer['reason']}.", key=layer["key"], reason=layer["reason"]
            )
        return layer["key"], self.adata.layers[layer["key"]], self.adata.var, self._candidates["var"]["columns"]

    def _names(self, axis: str, names: str) -> tuple[str | None, Any]:
        """Resolve one axis's gene-name set with the shared make-unique rule."""
        _key, _matrix, var, columns = self._axis(axis)

        def build() -> tuple[str | None, Any]:
            if names == "index":
                source: str | None = None
            else:
                position = int(names.removeprefix("column-"))
                if not 0 <= position < len(columns):
                    raise MappingNotFound(f"var column {position} does not exist.")
                column = columns[position]
                if not column["usable"]:
                    raise MappingLoadError(
                        f"var column {column['key']!r}: {column['reason']}.",
                        key=column["key"],
                        reason=column["reason"],
                    )
                source = column["key"]
            # The one gene-name rule of every path (data_mapping.column_gene_names).
            try:
                return source, column_gene_names(var, source)
            except DataMappingError as error:
                raise MappingLoadError(
                    str(error), key=source or "index", reason=str(error).rstrip(".")
                ) from error

        return cast(tuple[str | None, Any], self._cached(("names", axis, names), build))

    def var_manifest(self, axis: str, names: str) -> dict[str, Any]:
        source, resolved = self._names(axis, names)
        return {
            "_format": MANIFEST_FORMAT_VERSION,
            "n_points": int(self.adata.n_obs),
            "compression": None,
            "quantization": None,
            "_varSchema": {
                "kind": "continuous",
                "pathPattern": "var/{index}.values.f32",
                "ext": "f32",
                "dtype": "float32",
                "quantized": False,
            },
            "fields": [[position, name] for position, name in enumerate(resolved.names)],
            "gene_name_sets": [
                {
                    "source": source,
                    "filled_from_index": resolved.filled_from_index,
                    "renamed": [list(pair) for pair in resolved.renamed],
                    "names": None,
                }
            ],
        }

    def gene_count(self, axis: str) -> int:
        _key, _matrix, var, _columns = self._axis(axis)
        return int(var.shape[0])

    def _column(self, axis_key: str, matrix: Any, gene: int) -> np.ndarray:
        if sparse.issparse(matrix):
            if matrix.format != "csc" and not getattr(self.adata, "isbacked", False):
                with self._lock:
                    copy = self._csc_copies.get(axis_key)
                if copy is None:
                    copy = sparse.csc_matrix(matrix)
                    with self._lock:
                        self._csc_copies[axis_key] = copy
                matrix = copy
            return np.asarray(matrix[:, gene].toarray()).reshape(-1)
        column = matrix[:, gene]
        if sparse.issparse(column):
            return np.asarray(column.toarray()).reshape(-1)
        return np.asarray(column).reshape(-1)

    def gene_bytes(self, axis: str, names: str, gene: int) -> bytes:
        axis_key, matrix, _var, _columns = self._axis(axis)
        _source, resolved = self._names(axis, names)
        if not 0 <= gene < len(resolved.names):
            raise MappingNotFound(f"gene {gene} does not exist on axis {axis}.")
        name = resolved.names[gene]

        def build() -> bytes:
            column = self._column(axis_key, matrix, gene)
            try:
                values = _require_finite_float32_array(column, label=f"Gene {name!r} expression")
            except (TypeError, ValueError) as error:
                diagnosis = diagnose_continuous_payload(column, kind="gene", name=name)
                message = str(NonFinitePayloadError(diagnosis)) if diagnosis else str(error)
                defect = _value_defect(column)
                reason = f"{defect} in gene {name!r}" if defect is not None else message
                raise MappingLoadError(message, key=axis_key, reason=reason) from error
            return little_endian_payload_bytes(values)

        # Keyed by axis and position: the bytes do not depend on the name set.
        return cast(bytes, self._cached(("gene", axis, gene), build))

    # ------------------------------------------------------------------ latent

    def _categorical_field(self, field_index: int) -> str:
        keys = self._adapter.get_obs_keys()
        if not 0 <= field_index < len(keys):
            raise MappingNotFound(f"obs field {field_index} does not exist.")
        field_key = keys[field_index]
        if self._adapter.get_obs_field_kind(field_key) != "category":
            raise MappingNotFound(f"obs field {field_index} is not categorical.")
        return field_key

    def outlier_bytes(self, index: int, field_index: int) -> bytes:
        item = self._obsm_item(index)
        adapter = self._adapter
        field_key = self._categorical_field(field_index)
        latent_key = item["key"]

        # The latent space is already in memory as ``adata.obsm``, and a
        # float32 array is used as it is, without a copy. It is not cached:
        # only the per-cell quantiles are, so a latent space of any size costs
        # its own array plus one float32 per cell.
        def latent() -> np.ndarray:
            array = self._dense_obsm(latent_key)
            defect = _value_defect(array)
            if defect is not None:
                raise MappingLoadError(
                    f"Latent space {latent_key!r} cannot be used: {defect}.",
                    key=latent_key,
                    reason=defect,
                )
            try:
                return _require_finite_float32_array(array, label=f"Latent field {latent_key!r}")
            except (TypeError, ValueError) as error:
                raise MappingLoadError(str(error), key=latent_key, reason=str(error)) from error

        def build() -> bytes:
            values = latent()
            categorical = adapter.adata.obs[field_key].astype("category")
            categories = _json_category_values(categorical.cat.categories, field_key=field_key)
            quantiles = _compute_latent_space_quantiles(
                latent=values,
                codes=categorical.cat.codes.to_numpy(dtype=np.int32),
                categories=categories,
                min_points=adapter.centroid_min_points,
            )
            return little_endian_payload_bytes(quantiles)

        return cast(bytes, self._cached(("outliers", index, field_index), build))

    # ---------------------------------------------------------------- dispatch

    def head_length(self, route: MappingRoute) -> int:
        """The byte length a route answers with, without reading the payload.

        Embedding, vector, gene and latent lengths follow from the catalog.
        A graph's lengths are known once it has been read and checked, so a
        HEAD on an edge file reads it. JSON routes are built.
        """
        n_cells = int(self.adata.n_obs)
        if route.kind == "embedding":
            index, columns = route.parts
            self._require_columns(self._obsm_item(index), columns)
            return n_cells * len(columns) * 4
        if route.kind == "vectors":
            _item, _columns, dimension = self._vector_binding(*route.parts)
            return n_cells * dimension * 4
        if route.kind == "gene":
            axis, names, gene = route.parts
            self._names(axis, names)
            if not 0 <= gene < self.gene_count(axis):
                raise MappingNotFound(f"gene {gene} does not exist on axis {axis}.")
            return n_cells * 4
        if route.kind == "latent":
            index, field_index = route.parts
            self._obsm_item(index)
            self._categorical_field(field_index)
            return n_cells * 4
        if route.kind == "knn_edges":
            index, filename = route.parts
            return self.edge_length(index, filename)
        body, _content_type = self.payload(route)
        return len(body)

    def payload(self, route: MappingRoute) -> tuple[bytes, str]:
        """The bytes and content type of one mapping route."""
        if route.kind == "candidates":
            return json.dumps(self._candidates).encode("utf-8"), "application/json"
        if route.kind == "embedding":
            return self.embedding_bytes(*route.parts), "application/octet-stream"
        if route.kind == "vectors":
            return self.vector_bytes(*route.parts), "application/octet-stream"
        if route.kind == "knn_manifest":
            return (
                json.dumps(self.connectivity_manifest(*route.parts)).encode("utf-8"),
                "application/json",
            )
        if route.kind == "knn_edges":
            return self.edge_bytes(*route.parts), "application/octet-stream"
        if route.kind == "var_manifest":
            return json.dumps(self.var_manifest(*route.parts)).encode("utf-8"), "application/json"
        if route.kind == "gene":
            return self.gene_bytes(*route.parts), "application/octet-stream"
        if route.kind == "latent":
            return self.outlier_bytes(*route.parts), "application/octet-stream"
        raise MappingNotFound(f"Unknown mapping route kind {route.kind!r}.")

    def close(self) -> None:
        with self._lock:
            self._cache.clear()
            self._csc_copies.clear()
