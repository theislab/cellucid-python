"""
AnnData Adapter for Cellucid

Provides a server-side adapter that reads AnnData objects (in memory, backed h5ad,
or zarr stores) and serves data in the same format that the Cellucid web viewer
expects. This allows users to visualize AnnData directly without running
prepare first.

Key features:
- Works with in-memory AnnData, backed h5ad files, and eagerly loaded zarr stores
- Handles sparse matrices (CSR/CSC) transparently with automatic format conversion
- Resolves one embedding per dimension slot: X_umap_1d/X_umap_2d/X_umap_3d, or a
  bare X_umap at the dimension its own column count states, or any obsm array
  named explicitly with ``embedding_keys``; when none of those exist it makes a
  deterministic best guess (see :mod:`cellucid.data_mapping`)
- Reads expression from X, raw, or one layer (``expression_source``)
- Computes centroids and outlier quantiles on-demand with caching
- Lazy loading: gene expression and obs data are loaded on-demand
- No quantization (full float32 precision, gzip compression for network transfer)

Lazy Loading Behavior:
----------------------
- h5ad (backed='r'): Uses HDF5 memory-mapping. Arrays are loaded only when accessed.
  Gene expression columns are fetched individually, minimizing memory usage.
- zarr: ``anndata.read_zarr`` materializes the AnnData arrays in memory.
- in-memory: All data is already in RAM. No lazy loading.

Memory Management:
-----------------
The adapter maintains several caches for performance:
- Embedding cache: Normalized coordinates (one per served embedding)
- Centroid cache: Computed label centroids
- Gene expression LRU cache: Recently accessed gene columns (max 100 genes)
- CSC cache: For CSR matrices, a CSC copy is created for efficient column access

IMPORTANT: Always use the context manager or call close() when done:
    with AnnDataAdapter.from_file(
        "data.h5ad",
        dataset_name="Example",
        dataset_id="example",
    ) as adapter:
        # use adapter
    # resources are automatically released

Usage:
    from cellucid.anndata_adapter import AnnDataAdapter

    # From h5ad file (lazy loading via backed mode)
    with AnnDataAdapter.from_file(
        "/path/to/data.h5ad",
        dataset_name="Example",
        dataset_id="example",
    ) as adapter:
        manifest = adapter.get_obs_manifest()

    # From a zarr store (materialized by anndata.read_zarr)
    adapter = AnnDataAdapter.from_file(
        "/path/to/data.zarr",
        dataset_name="Example",
        dataset_id="example",
    )

    # From in-memory AnnData
    adapter = AnnDataAdapter(
        adata,
        dataset_name="Example",
        dataset_id="example",
    )

    # Use with server
    from cellucid import show_anndata
    show_anndata(
        adata,
        dataset_name="Example",
        dataset_id="example",
    )
"""

from __future__ import annotations

import json
import logging
import math
import re
import sys
import threading
from collections import OrderedDict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, cast

import numpy as np
import pandas as pd
from scipy import sparse

from ._byte_order import little_endian_payload_bytes
from ._compression import deterministic_gzip_compress
from ._contracts import (
    EMBEDDINGS_DIRNAME,
    IDENTITY_VERSION,
    _require_dataset_id,
    _require_dataset_name,
    _require_field_identities,
    _require_native_boolean,
    _require_nonempty_string,
    _require_positive_native_integer,
)
from ._server_base import print_detail
from .connectivity_contract import (
    ConnectivityEdgePairs,
    build_connectivity_manifest,
    validate_connectivity_edges,
)
from .continuous_payload_diagnosis import (
    NonFinitePayloadError,
    diagnose_continuous_payload,
)
from .data_mapping import (
    DIMENSION_SLOTS,
    DIMENSIONAL_UMAP_KEYS,
    VECTOR_KEY_PATTERN,
    DataMappingError,
    EmbeddingKeySpec,
    ObsmEntry,
    UniqueNames,
    best_guess_embedding_slots,
    build_anndata_candidates,
    canonical_embedding_id,
    default_columns,
    default_embedding_slots,
    default_vector_bindings,
    describe_obsm,
    describe_obsm_entry,
    expression_layer,
    make_unique_names,
    nothing_plottable_message,
    obsm_entries,
    obsm_entry_defect,
    slot_dimension,
    unique_names_note,
    validate_data_mapping,
    validate_embedding_columns,
)
from .prepare_data._arrays import (
    _normalize_finite_float32_embedding,
    _require_continuous_obs_values,
    _require_finite_embedding_source,
    _require_finite_float32_array,
)
from .prepare_data._categories import _json_category_values
from .prepare_data._manifest import MANIFEST_FORMAT_VERSION
from .vector_fields import (
    _finite_float32_matrix,
    scale_vector_field,
    validate_vector_fields,
)

if TYPE_CHECKING:
    import anndata

logger = logging.getLogger("cellucid.anndata_adapter")

#: The route, relative to a served dataset, of the data-mapping catalog: the
#: mapping the server serves (``defaults``) and every candidate array.
MAPPING_CANDIDATES_PATH = "_cellucid/mapping/candidates.json"


def _is_empty_gene_name(value: object) -> bool:
    """Whether one gene-name value is empty: None, "", NaN, or a pandas missing value."""
    if value is None or value == "":
        return True
    if isinstance(value, float):
        return math.isnan(value)
    return type(value).__name__ in {"NAType", "NaTType"}


@dataclass(frozen=True)
class ServedEmbedding:
    """One embedding the adapter serves: an ``obsm`` key read at fixed columns.

    ``id`` is the canonical embedding id (``key``, or ``key[c0,c1]`` when only
    some columns are read), ``index`` its position in (dimension, id) order and
    therefore its payload file, ``embeddings/{index}_{dimension}d.bin``.
    """

    id: str
    index: int
    key: str
    columns: tuple[int, ...]
    dimension: int
    width: int

    @property
    def file(self) -> str:
        return f"{EMBEDDINGS_DIRNAME}/{self.index}_{self.dimension}d.bin"


def _require_embedding_keys(value: object) -> dict[str, tuple[str, list[int] | None]]:
    """Validate ``embedding_keys``: ``{"2d": "X_tsne"}`` or ``{"2d": ("X_pca", [0, 1])}``."""
    if value is None:
        return {}
    if isinstance(value, str | bytes) or not isinstance(value, Mapping):
        raise TypeError(
            "embedding_keys must map '1d', '2d', or '3d' to an obsm key or to a "
            "(key, columns) pair, for example {'2d': 'X_tsne'}."
        )
    resolved: dict[str, tuple[str, list[int] | None]] = {}
    for slot, spec in value.items():
        if slot not in DIMENSION_SLOTS:
            raise ValueError(f"embedding_keys key {slot!r} must be '1d', '2d', or '3d'.")
        if isinstance(spec, str):
            key, columns = spec, None
        elif isinstance(spec, tuple) and len(spec) == 2 and isinstance(spec[0], str):
            key = spec[0]
            if isinstance(spec[1], str | bytes) or not isinstance(spec[1], Sequence):
                raise TypeError(f"embedding_keys[{slot!r}] columns must be a sequence of integers.")
            columns = list(spec[1])
        else:
            raise TypeError(
                f"embedding_keys[{slot!r}] must be an obsm key or a (key, columns) pair."
            )
        if not key:
            raise ValueError(f"embedding_keys[{slot!r}] must name a non-empty obsm key.")
        resolved[slot] = (key, columns)
    return resolved


def _require_expression_source(value: object) -> str:
    """Validate ``expression_source``: ``"X"``, ``"raw"`` or ``"layer:<name>"``."""
    if not isinstance(value, str) or not (
        value in ("X", "raw")
        or (value.startswith("layer:") and len(value) > len("layer:"))
    ):
        raise ValueError(
            "expression_source must be 'X', 'raw', or 'layer:<name>', "
            f"got {value!r}."
        )
    return value


def _classify_anndata_path(path: str | Path) -> Literal["h5ad", "zarr"]:
    """Classify one explicitly declared AnnData file or materialized Zarr store."""
    resolved_path = Path(path).resolve()
    if not resolved_path.exists():
        raise FileNotFoundError(f"Path not found: {resolved_path}")

    if resolved_path.is_file():
        if resolved_path.suffix != ".h5ad":
            raise ValueError(f"AnnData files must use the exact .h5ad extension: {resolved_path}")
        return "h5ad"

    if not resolved_path.is_dir():
        raise ValueError(
            f"AnnData path must be an .h5ad file or a materialized Zarr directory: {resolved_path}"
        )

    zarr_v2_group = resolved_path / ".zgroup"
    zarr_v2_attrs = resolved_path / ".zattrs"
    zarr_v3_metadata = resolved_path / "zarr.json"
    has_v2_marker = zarr_v2_group.exists() or zarr_v2_attrs.exists()

    if zarr_v3_metadata.exists():
        raise ValueError(
            f"AnnData Zarr input must use the exact Zarr v2 root contract "
            f"(.zgroup and .zattrs, without zarr.json): {resolved_path}"
        )

    if has_v2_marker:
        if not zarr_v2_group.is_file() or not zarr_v2_attrs.is_file():
            raise ValueError(
                f"Zarr v2 directory must contain both .zgroup and .zattrs files: {resolved_path}"
            )
        try:
            group_metadata = json.loads(zarr_v2_group.read_text(encoding="utf-8"))
            attributes = json.loads(zarr_v2_attrs.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            raise ValueError(
                f"Zarr v2 root metadata must be readable UTF-8 JSON: {resolved_path}"
            ) from error
        if group_metadata != {"zarr_format": 2}:
            raise ValueError(f"Zarr v2 .zgroup must declare exactly zarr_format 2: {resolved_path}")
        if not isinstance(attributes, dict):
            raise ValueError(f"Zarr v2 .zattrs must contain a JSON object: {resolved_path}")
        return "zarr"

    raise ValueError(
        f"AnnData path must be an .h5ad file or a complete Zarr v2 directory: {resolved_path}"
    )


# A payload route is the field's integer index, written exactly as the manifest
# writes it, so '0' resolves and '00', '+0', and ' 0' do not.
_PAYLOAD_INDEX_PATTERN = re.compile(r"0|[1-9][0-9]*", flags=re.ASCII)


def _payload_route_index(component: str) -> int | None:
    """Resolve one canonical decimal payload index, or None when it is not one."""
    if _PAYLOAD_INDEX_PATTERN.fullmatch(component) is None:
        return None
    return int(component)


def _obs_field_kind(series: pd.Series, *, key: str) -> Literal["continuous", "category"]:
    """Classify one obs column, naming the exact escape hatch when it cannot be served.

    Classification rules:
    - Categorical dtype → category
    - Boolean dtype → category
    - Numeric dtype → continuous
    - String/object → category (treated as labels)
    """
    if isinstance(series.dtype, pd.CategoricalDtype) or pd.api.types.is_bool_dtype(series):
        return "category"
    if pd.api.types.is_numeric_dtype(series):
        return "continuous"
    if pd.api.types.is_string_dtype(series) or pd.api.types.is_object_dtype(series):
        return "category"
    raise TypeError(
        f"obs field {key!r} has unsupported dtype {series.dtype!r}. Cellucid serves "
        "categorical, boolean, numeric, string, and object obs columns. Fix: convert "
        f'it, for example adata.obs[{key!r}] = adata.obs[{key!r}].astype(str), or '
        f"leave it out with obs_keys=[...] naming only the fields to serve."
    )


def _resolve_obs_keys(obs: pd.DataFrame, requested: Any) -> list[str]:
    """Resolve the exact obs columns this adapter serves."""
    available = list(obs.columns)
    if requested is None:
        return available
    if isinstance(requested, str | bytes):
        raise TypeError("obs_keys must be a sequence of strings, not a string scalar.")
    resolved = list(requested)
    if not resolved:
        raise ValueError("obs_keys must name at least one obs column.")
    missing = [key for key in resolved if key not in available]
    if missing:
        raise KeyError(
            f"obs_keys contain columns not in obs: {missing}. Available columns: {available}"
        )
    return resolved


def _categorical_storage(
    n_categories: int,
    *,
    field_key: str,
) -> tuple[type[np.uint8] | type[np.uint16], str, int]:
    if n_categories > 65_535:
        raise ValueError(
            f"Categorical field {field_key!r} has {n_categories:,} categories; "
            "the current contract supports at most 65,535"
        )
    if n_categories <= 255:
        return np.uint8, "uint8", 255
    return np.uint16, "uint16", 65_535


#: How much resident memory the served gene columns may occupy, in bytes.
#:
#: Counting entries instead of bytes is what made this dangerous: a gene column
#: is one float32 per cell, so a hundred of them is 29 MB on a 72,000-cell object
#: and **7.3 GB** on an 18,142,044-cell one — and the second is exactly the size
#: of object the cache exists for. A byte budget is the same policy expressed in
#: the unit that actually runs out, and it holds 32 columns of a 2 M-cell object
#: or 3 of an 18 M-cell one.
GENE_EXPRESSION_CACHE_BYTES = 256 * 1024 * 1024


class LRUCache:
    """
    LRU cache over a byte budget, with O(1) operations using OrderedDict.

    WARNING: This class is NOT thread-safe. The AnnDataAdapter is designed
    for single-threaded use within an HTTP server (each request is handled
    sequentially). If you need thread-safe caching, wrap operations in
    external locks or use a thread-safe cache implementation.

    The HTTP server uses a single adapter instance per server, and Python's
    GIL provides some level of safety for simple operations, but concurrent
    access from multiple threads may cause race conditions.

    A single value larger than the whole budget is not an error and is not
    cached: refusing to serve an 18 M-cell gene because it will not fit in a
    cache would be the cache deciding what the server may publish.
    """

    def __init__(self, max_bytes: int = GENE_EXPRESSION_CACHE_BYTES):
        if max_bytes < 1:
            raise ValueError("max_bytes must be at least 1")
        self._cache: OrderedDict[str, Any] = OrderedDict()
        self._sizes: OrderedDict[str, int] = OrderedDict()
        self._max_bytes = max_bytes
        self._resident_bytes = 0

    @staticmethod
    def _value_bytes(value: Any) -> int:
        """How much memory one cached value holds."""
        nbytes = getattr(value, "nbytes", None)
        if isinstance(nbytes, int) and nbytes >= 0:
            return nbytes
        return sys.getsizeof(value)

    @property
    def resident_bytes(self) -> int:
        """How much memory the cache is currently holding."""
        return self._resident_bytes

    def get(self, key: str) -> Any | None:
        """Get item and move to end (most recently used)."""
        if key in self._cache:
            self._cache.move_to_end(key)
            self._sizes.move_to_end(key)
            return self._cache[key]
        return None

    def put(self, key: str, value: Any) -> None:
        """Add item, evicting least-recently-used entries to stay in budget."""
        self._discard(key)
        size = self._value_bytes(value)
        if size > self._max_bytes:
            return
        while self._cache and self._resident_bytes + size > self._max_bytes:
            oldest, _ = self._cache.popitem(last=False)
            self._resident_bytes -= self._sizes.pop(oldest)
        self._cache[key] = value
        self._sizes[key] = size
        self._resident_bytes += size

    def _discard(self, key: str) -> None:
        if key in self._cache:
            del self._cache[key]
            self._resident_bytes -= self._sizes.pop(key)

    def __contains__(self, key: str) -> bool:
        return key in self._cache

    def clear(self) -> None:
        """Clear all cached items and release memory."""
        self._cache.clear()
        self._sizes.clear()
        self._resident_bytes = 0

    def __len__(self) -> int:
        return len(self._cache)


def _to_dense_1d(arr: np.ndarray | sparse.spmatrix) -> np.ndarray:
    """Convert sparse vector to dense numpy array."""
    if sparse.issparse(arr):
        return np.asarray(cast(sparse.spmatrix, arr).toarray()).flatten()
    arr = np.asarray(arr)
    if arr.ndim > 1:
        return arr.flatten()
    return arr


def _to_dense_2d(arr: np.ndarray | sparse.spmatrix) -> np.ndarray:
    """Convert sparse matrix to dense numpy array."""
    if sparse.issparse(arr):
        return np.asarray(cast(sparse.spmatrix, arr).toarray())
    return np.asarray(arr)


class AnnDataAdapter:
    """
    Adapter that wraps AnnData and provides data in Cellucid format.

    This adapter generates all the data that would normally be created by
    prepare, but reads directly from AnnData without creating
    intermediate files. This is slower but more convenient for interactive use.
    """

    def __init__(
        self,
        adata: anndata.AnnData,
        latent_key: str | None = None,
        gene_id_column: str | None = None,
        normalize_embeddings: bool = True,
        centroid_outlier_quantile: float = 0.95,
        centroid_min_points: int = 10,
        *,
        dataset_name: str,
        dataset_id: str,
        obs_keys: Sequence[str] | None = None,
        embedding_keys: Mapping[str, EmbeddingKeySpec] | None = None,
        knn_key: str | None = None,
        expression_source: str = "X",
        vector_field_default: str | None = None,
        serve_vector_fields: bool = False,
        quiet: bool = True,
    ) -> None:
        """
        Initialize the adapter.

        Parameters
        ----------
        adata : AnnData
            AnnData object to adapt. Can be in-memory or backed (h5ad file).
        latent_key : str, optional
            Key in obsm for latent space used for outlier quantile calculation.
            If None, outlier quantiles are explicitly unavailable.
        gene_id_column : str, optional
            Exact column of the served gene axis (``var``, or ``raw.var`` when
            ``expression_source="raw"``) whose values name the genes. If None,
            names come from that axis's index. An empty value takes the
            row's index name, and repeated names are made unique the way
            anndata's ``var_names_make_unique()`` does (``MATR3``, ``MATR3-1``).
        normalize_embeddings : bool
            If True, normalize embeddings to [-1, 1] range (recommended).
        centroid_outlier_quantile : float
            Quantile for outlier removal in centroid computation.
        centroid_min_points : int
            Minimum points per category for centroid computation.
        dataset_name : str
            Exact non-empty human-readable name without surrounding whitespace
            or control characters. Unicode is preserved.
        dataset_id : str
            Portable 1-180 byte ASCII identifier: begin with a letter or
            digit, then use only letters, digits, ``.``, ``_``, or ``-``.
            It cannot end in ``.`` or use a reserved Windows device name.
        obs_keys : sequence of str, optional
            Exact ``obs`` columns to serve, in this order. If None, every
            column is served. Naming the columns is how a dataset carrying an
            unservable column, such as a ``datetime64`` collection date, is
            served without editing ``adata``.
        embedding_keys : mapping, optional
            The embedding to serve in a dimension slot, overriding the default
            for that slot: ``{"2d": "X_tsne"}`` reads the first two columns of
            ``obsm['X_tsne']``, and ``{"2d": ("X_pca", [2, 0])}`` reads the
            0-based columns named. Slots not named follow the default rule:
            ``X_umap_1d``/``X_umap_2d``/``X_umap_3d``, or a bare 1- to 3-column
            ``X_umap``. When neither names any slot, Cellucid makes a
            deterministic best guess from the other ``obsm`` arrays.
        knn_key : str, optional
            The ``obsp`` matrix to serve as the neighbour graph, usually
            ``"connectivities"``. None by default: an edge list costs time and
            memory proportional to the number of stored neighbours, which on a
            large object is the single longest part of startup, and most
            sessions never draw the graph. The matrix is read and checked when
            it is named -- it must be square, finite, non-negative, symmetric,
            and zero on the diagonal -- and a matrix that is not symmetric, such
            as ``distances``, is refused with that reason.
        expression_source : str
            Where gene expression is read: ``"X"`` (the default), ``"raw"``
            (``raw.X``, with genes named from ``raw.var``), or
            ``"layer:<name>"``.
        vector_field_default : str, optional
            Exact field id to select when multiple UMAP vector fields exist.
        serve_vector_fields : bool
            Whether to read, validate, and serve the ``<field>_umap_<n>d`` vector
            arrays in ``obsm``. A key binds to the slot of its dimension only
            when that slot shows ``X_umap_<n>d`` or ``X_umap``. False by
            default, for the same reason as ``knn_key``: every declared field
            is read and checked in full before the server binds, and a session
            that never turns the overlay on pays nothing for it now.
        quiet : bool
            Whether to build without reporting progress. True by default, so
            using this class as a library prints nothing. A server passes its
            own ``quiet`` through, because opening a large object spends minutes
            inside this constructor, and a terminal that says nothing for that
            long is indistinguishable from one that has stopped.
        """
        dataset_name = _require_dataset_name(
            dataset_name,
            label="dataset_name",
        )
        dataset_id = _require_dataset_id(
            dataset_id,
            label="dataset_id",
        )
        self.adata: anndata.AnnData = adata
        if adata.isbacked and adata.file._filemode != "r":
            raise ValueError(
                "Backed AnnData inputs must be opened read-only. Pass the .h5ad path "
                "to Cellucid so it can own the exact read-only backing handle."
            )
        if adata.n_obs <= 0:
            raise ValueError("AnnData datasets must contain at least one observation")
        self.latent_key: str | None = latent_key
        if latent_key is not None:
            if not isinstance(latent_key, str) or not latent_key:
                raise ValueError("latent_key must be None or one non-empty string")
            if latent_key not in adata.obsm:
                raise ValueError(f"latent_key {latent_key!r} was not found in adata.obsm")
        self.gene_id_column: str | None
        if gene_id_column is None:
            self.gene_id_column = None
        else:
            self.gene_id_column = _require_nonempty_string(
                gene_id_column,
                label="gene_id_column",
            )
            if not self.gene_id_column.strip():
                raise ValueError("gene_id_column must be None or a non-blank string.")
        self.normalize_embeddings = _require_native_boolean(
            normalize_embeddings,
            label="normalize_embeddings",
        )
        self.knn_key: str | None
        if knn_key is None:
            self.knn_key = None
        else:
            self.knn_key = _require_nonempty_string(knn_key, label="knn_key")
        self.expression_source: str = _require_expression_source(expression_source)
        self.serve_vector_fields = _require_native_boolean(
            serve_vector_fields,
            label="serve_vector_fields",
        )
        self.quiet = _require_native_boolean(quiet, label="quiet")
        requested_embeddings = _require_embedding_keys(embedding_keys)

        if (
            type(centroid_outlier_quantile) is not float
            or not math.isfinite(centroid_outlier_quantile)
            or not (0.5 < centroid_outlier_quantile < 1.0)
        ):
            raise ValueError(
                "centroid_outlier_quantile must be a finite number strictly between 0.5 and 1.0"
            )
        self.centroid_outlier_quantile: float = centroid_outlier_quantile

        self.centroid_min_points = _require_positive_native_integer(
            centroid_min_points,
            label="centroid_min_points",
        )

        self.dataset_name: str = dataset_name
        self.dataset_id: str = dataset_id

        # The expression matrix and the gene axis it belongs to. raw carries its
        # own var, so reading expression from raw also names genes from raw.var.
        self._expression_matrix: Any = self._resolve_expression_matrix()
        gene_axis = self._gene_axis_frame()
        self._gene_names: UniqueNames = self._resolve_gene_names(gene_axis)
        gene_ids = self._gene_names.names
        if self._expression_matrix is not None:
            # Only names a gene axis actually serves reach the viewer, so only
            # served names are held to the display rule.
            gene_ids = _require_field_identities(
                gene_ids,
                label=(
                    "gene axis index"
                    if self.gene_id_column is None
                    else f"gene axis column {self.gene_id_column!r}"
                ),
            )

        # Public source attribution never includes the local filesystem path.
        self._source_type: str = "memory"

        # Caches for computed data, keyed by served embedding id.
        self._embedding_cache: dict[str, np.ndarray] = {}
        # Per-embedding normalization info (center + scale_factor). Used to
        # scale vector fields into the same normalized space as the points.
        self._embedding_norm: dict[str, dict[str, Any]] = {}

        # Optional per-cell vector fields (e.g. velocity, CellRank drift).
        # Structure matches dataset_identity.json: {"default_field": ..., "fields": {...}}.
        self._vector_fields_metadata: dict[str, Any] | None = None
        self._vector_field_obsm_keys: dict[str, dict[str, str]] = {}
        self._vector_field_payload_index_by_id: dict[str, int] = {}
        self._vector_field_cache: dict[tuple[str, str], np.ndarray] = {}
        self._centroid_cache: dict[tuple[str, str], list[dict[str, Any]]] = {}
        self._outlier_quantile_cache: dict[str, np.ndarray] = {}
        self._latent_space: np.ndarray | None = None
        self._gene_ids_cache: list[str] | None = gene_ids
        self._gene_id_to_idx_cache: dict[str, int] | None = {
            gene_id: index for index, gene_id in enumerate(gene_ids)
        }
        # A gene's payload index is its position in the served gene list, which
        # is exactly the index its var manifest entry declares.
        self._gene_payload_index_by_id: dict[str, int] = {
            gene_id: index for index, gene_id in enumerate(gene_ids)
        }

        self._connectivity_cache: ConnectivityEdgePairs | None = None
        self._connectivity_manifest: dict[str, Any] | None = None
        self._has_connectivity: bool = False

        # CSC cache for efficient column access on CSR matrices
        # Note: This doubles memory for sparse matrices, but is necessary for O(1) gene access
        # Set to None initially; created lazily only when needed
        self._X_csc_cache: sparse.csc_matrix | None = None

        # LRU cache for gene expression values using O(1) OrderedDict
        self._gene_expression_cache: LRUCache = LRUCache()

        # Served embeddings, in (dimension, id) order, and each slot's default.
        self._embeddings: list[ServedEmbedding] = []
        self._embedding_defaults: dict[str, str] = {}
        self._guessed_slots: list[str] = []
        # The data-mapping catalog, built on first request and then frozen.
        # The server is threaded, so building it is guarded by a lock.
        self._mapping_candidates: dict[str, Any] | None = None
        self._mapping_candidates_lock = threading.Lock()
        # Track if adapter has been closed
        self._closed: bool = False
        self._close_complete: bool = False

        selected_obs_keys = _require_field_identities(
            _resolve_obs_keys(self.adata.obs, obs_keys),
            label="Observation field",
        )
        self._obs_keys: list[str] = selected_obs_keys
        # obs/ carries both continuous and categorical payloads, so one index
        # space spans both manifest arrays, in the served column order.
        self._obs_payload_index_by_key: dict[str, int] = {
            key: index for index, key in enumerate(selected_obs_keys)
        }
        # Every served column is classified now, so an unservable dtype fails
        # here with its exact escape hatch instead of inside a later request.
        self._report("Obs columns", f"classifying {len(selected_obs_keys):,}")
        for key in selected_obs_keys:
            series = self.adata.obs[key]
            if _obs_field_kind(series, key=key) == "continuous":
                column = series.to_numpy()
                try:
                    _require_continuous_obs_values(
                        column,
                        key=key,
                        n_cells=int(self.adata.n_obs),
                    )
                except ValueError as exc:
                    # This one refuses to start the server at all, so it is the
                    # message with the least room to be vague: it has to name the
                    # column and say how much of it is unusable, or the only way
                    # to find out is to bisect the object by hand.
                    diagnosis = diagnose_continuous_payload(
                        column, kind="field", name=key
                    )
                    if diagnosis is None:
                        raise
                    raise NonFinitePayloadError(diagnosis) from exc

        note = unique_names_note(
            filled_from_index=self._gene_names.filled_from_index,
            renamed=self._gene_names.renamed,
            source=self.gene_id_column,
        )
        if note is not None and self._expression_matrix is not None:
            self._report("Gene names", note)

        # Resolve the embedding of every dimension slot.
        self._report("Embeddings", "resolving obsm keys")
        self._resolve_embeddings(requested_embeddings)
        self._default_dimension: int = max(item.dimension for item in self._embeddings)
        self._report(
            "Embeddings",
            ", ".join(
                f"{item.dimension}D from obsm[{item.key!r}]"
                + (
                    ""
                    if list(item.columns) == default_columns(item.dimension)
                    else f" columns {list(item.columns)}"
                )
                + (" (guessed)" if f"embeddings.{item.dimension}d" in self._guessed_slots else "")
                for item in self._embeddings
            ),
        )

        # Per-cell vector fields are an asked-for capability, like connectivity.
        if self.serve_vector_fields:
            self._report("Vector fields", "scanning obsm")
            self._vector_fields_metadata = self._detect_vector_fields(
                vector_field_default=vector_field_default,
            )
            served_field_count = (
                len(self._vector_fields_metadata["fields"]) if self._vector_fields_metadata else 0
            )
            self._report("Vector fields", f"{served_field_count} served")
        else:
            if vector_field_default is not None:
                raise ValueError(
                    "vector_field_default names the field to select among served "
                    "vector fields, and vector fields are not being served. Pass "
                    "serve_vector_fields=True, or --vector-fields, to serve them, "
                    "or drop vector_field_default."
                )
            # Not scanning obsm is the point: every declared field is otherwise
            # read and checked in full before the server can bind.
            self._vector_fields_metadata = None

        # Connectivity is an asked-for scientific capability. When it is asked
        # for, its exact edge contract is validated and materialized before this
        # adapter can be used, because the manifest the viewer reads first
        # declares the edge count and the neighbor maximum.
        self._initialize_connectivity()
        if self._connectivity_manifest is not None:
            self._report(
                "Connectivity",
                f"{self._connectivity_manifest['n_edges']:,} edges, "
                f"{self._connectivity_manifest['max_neighbors']:,} neighbors at most",
            )

        logger.info(
            f"AnnDataAdapter initialized: {self.n_cells:,} cells, "
            f"{self.n_genes:,} genes, embeddings: {[item.id for item in self._embeddings]}"
        )

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        *,
        latent_key: str | None = None,
        gene_id_column: str | None = None,
        normalize_embeddings: bool = True,
        centroid_outlier_quantile: float = 0.95,
        centroid_min_points: int = 10,
        dataset_name: str,
        dataset_id: str,
        obs_keys: Sequence[str] | None = None,
        embedding_keys: Mapping[str, EmbeddingKeySpec] | None = None,
        knn_key: str | None = None,
        expression_source: str = "X",
        vector_field_default: str | None = None,
        serve_vector_fields: bool = False,
        quiet: bool = True,
    ) -> AnnDataAdapter:
        """
        Create an adapter from an h5ad file or a materialized zarr store.

        Supports both:
        - .h5ad files: HDF5-based, supports true backed mode with memory-mapping
        - .zarr directories: loaded eagerly by ``anndata.read_zarr``

        Parameters
        ----------
        path : str or Path
            Path to h5ad file or zarr directory.
            - For h5ad: path/to/file.h5ad
            - For zarr: path/to/store.zarr (must be a directory)
        latent_key : str, optional
            Explicit key in ``obsm`` for the latent space.
        gene_id_column : str, optional
            Exact column of the served gene axis whose values name the genes.
            If None, names come from that axis's index.
        normalize_embeddings : bool
            Whether to normalize embeddings into the viewer coordinate range.
        centroid_outlier_quantile : float
            Quantile used when computing categorical centroids.
        centroid_min_points : int
            Minimum category size used for centroid computation.
        dataset_name : str
            Exact non-empty human-readable name without surrounding whitespace
            or control characters. Unicode is preserved.
        dataset_id : str
            Portable 1-180 byte ASCII identifier: begin with a letter or
            digit, then use only letters, digits, ``.``, ``_``, or ``-``.
            It cannot end in ``.`` or use a reserved Windows device name.
        obs_keys : sequence of str, optional
            Exact ``obs`` columns to serve, in this order. If None, every
            column is served.
        embedding_keys : mapping, optional
            The embedding of a dimension slot, such as ``{"2d": "X_tsne"}`` or
            ``{"2d": ("X_pca", [0, 1])}``. Slots not named follow the default
            rule, and a best guess applies when no slot resolves.
        knn_key : str, optional
            The ``obsp`` matrix to serve as the neighbour graph. None by
            default, because building the edge list is the longest part of
            opening a large object.
        expression_source : str
            ``"X"`` (the default), ``"raw"``, or ``"layer:<name>"``.
        vector_field_default : str, optional
            Exact field id required when multiple UMAP vector fields exist.
        serve_vector_fields : bool
            Whether to read, validate, and serve the ``<field>_umap_<n>d`` vector
            arrays in ``obsm``. False by default.
        quiet : bool
            Whether to build without reporting progress. True by default.
        Returns
        -------
        AnnDataAdapter
            Adapter instance wrapping the loaded data.

        Raises
        ------
        FileNotFoundError
            If the path does not exist.
        ValueError
            If the path is not a valid h5ad or zarr store.

        Examples
        --------
        >>> # H5AD is always opened in read-only backed mode
        >>> adapter = AnnDataAdapter.from_file(
        ...     "data.h5ad",
        ...     dataset_name="Example",
        ...     dataset_id="example",
        ... )
        """
        import anndata as ad

        dataset_name = _require_dataset_name(
            dataset_name,
            label="dataset_name",
        )
        dataset_id = _require_dataset_id(
            dataset_id,
            label="dataset_id",
        )
        path = Path(path).resolve()

        source_type = _classify_anndata_path(path)
        if source_type == "zarr":
            adata = ad.read_zarr(path)
            logger.info("Loaded zarr store eagerly: %s", path)
        else:
            adata = ad.read_h5ad(path, backed="r")
            logger.info(
                "Loaded h5ad file: %s (is_backed=%s)",
                path,
                bool(adata.isbacked),
            )

        # Create adapter and store non-sensitive source info. A failed adapter
        # construction must not leave a backed file handle open.
        try:
            adapter = cls(
                adata,
                latent_key=latent_key,
                gene_id_column=gene_id_column,
                normalize_embeddings=normalize_embeddings,
                centroid_outlier_quantile=centroid_outlier_quantile,
                centroid_min_points=centroid_min_points,
                dataset_name=dataset_name,
                dataset_id=dataset_id,
                obs_keys=obs_keys,
                embedding_keys=embedding_keys,
                knn_key=knn_key,
                expression_source=expression_source,
                vector_field_default=vector_field_default,
                serve_vector_fields=serve_vector_fields,
                quiet=quiet,
            )
        except BaseException:
            if getattr(adata, "isbacked", False):
                try:
                    adata.file.close()
                except BaseException:
                    logger.exception(
                        "Failed to close the owned AnnData file after adapter construction failed"
                    )
            raise
        adapter._source_type = source_type
        return adapter

    def _check_closed(self) -> None:
        """Raise error if adapter has been closed."""
        if self._closed:
            raise RuntimeError("AnnDataAdapter has been closed. Create a new adapter instance.")

    @property
    def n_cells(self) -> int:
        """Number of cells."""
        self._check_closed()
        return cast(int, self.adata.n_obs)

    @property
    def n_genes(self) -> int:
        """Number of genes on the served gene axis (``var``, or ``raw.var``)."""
        self._check_closed()
        if self.expression_source == "raw":
            return int(self.adata.raw.n_vars)
        return cast(int, self.adata.n_vars)

    @property
    def is_backed(self) -> bool:
        """Whether the live AnnData object reports an active backed store."""
        self._check_closed()
        return cast(bool, self.adata.isbacked)

    # =========================================================================
    # DIMENSION DETECTION
    # =========================================================================

    def _report(self, label: str, value: str) -> None:
        """Report one phase of a build that can run for minutes."""
        if self.quiet:
            return
        print_detail(label, value)

    def _resolve_expression_matrix(self) -> Any:
        """Return the expression matrix ``expression_source`` names, or None for no X."""
        source = self.expression_source
        if source == "X":
            return self.adata.X
        if source == "raw":
            raw = getattr(self.adata, "raw", None)
            if raw is None or raw.X is None:
                raise ValueError(
                    "expression_source='raw' reads adata.raw.X, and this object has no raw."
                )
            return raw.X
        layer = expression_layer(source)
        if layer is None or layer not in self.adata.layers:
            raise ValueError(
                f"expression_source={source!r} names no layer of this object. "
                f"adata.layers holds: {list(self.adata.layers.keys())}."
            )
        matrix = self.adata.layers[layer]
        if tuple(matrix.shape) != (self.adata.n_obs, self.adata.n_vars):
            raise ValueError(
                f"adata.layers[{layer!r}] has shape {tuple(matrix.shape)}; expected "
                f"({self.adata.n_obs}, {self.adata.n_vars})."
            )
        return matrix

    def _gene_axis_frame(self) -> pd.DataFrame:
        """The ``var`` frame of the served gene axis."""
        if self.expression_source == "raw":
            return cast(pd.DataFrame, self.adata.raw.var)
        return cast(pd.DataFrame, self.adata.var)

    def _resolve_gene_names(self, var: pd.DataFrame) -> UniqueNames:
        """Name every gene of the served axis with the shared make-unique rule."""
        index = [str(name) for name in var.index.tolist()]
        if self.gene_id_column is None:
            return make_unique_names(index)
        if self.gene_id_column not in var.columns:
            axis = "raw.var" if self.expression_source == "raw" else "var"
            raise KeyError(
                f"gene_id_column {self.gene_id_column!r} not found in {axis}. "
                f"Available columns: {list(var.columns)}"
            )
        values = var[self.gene_id_column].tolist()
        for position, value in enumerate(values):
            # Names are never coerced: a numeric column is a different
            # identifier space, and stringifying it would invent names.
            if not isinstance(value, str) and not _is_empty_gene_name(value):
                raise TypeError(
                    f"var column {self.gene_id_column!r} identifier at position {position} "
                    f"must be a string, got {type(value).__name__}."
                )
        try:
            return make_unique_names(values, index=index)
        except DataMappingError as error:
            raise ValueError(f"gene_id_column {self.gene_id_column!r}: {error}") from error

    def _obsm_source(self, key: str, columns: Sequence[int]) -> np.ndarray:
        """Read the selected columns of one dense ``obsm`` array, unnormalized."""
        value = self.adata.obsm[key]
        array = np.asarray(value)
        if array.ndim == 1:
            array = array.reshape(-1, 1)
        return _require_finite_embedding_source(
            array[:, list(columns)],
            label=f"Embedding {key!r}",
        )

    def _embedding_load_defect(self, key: str, columns: list[int]) -> str | None:
        """Return why one ``obsm`` array cannot be plotted at ``columns``, or None.

        This is the load-time check of the data mapping: the selected columns
        must be finite, representable as float32, and span a nonzero range.
        """
        array = np.asarray(self.adata.obsm[key])
        if array.ndim == 1:
            array = array.reshape(-1, 1)
        selected = array[:, columns]
        if selected.dtype.kind not in {"i", "u", "f"}:
            return "not a numeric array"
        with np.errstate(over="ignore", invalid="ignore"):
            finite = np.isfinite(selected)
            if not finite.all():
                row = int(np.argwhere(~finite)[0][0])
                value = float(selected[row][~finite[row]][0])
                kind = "NaN" if math.isnan(value) else ("Infinity" if value > 0 else "-Infinity")
                return f"{kind} at cell {row}"
            narrowed = selected.astype(np.float32)
            bad = ~np.isfinite(narrowed) | ((selected != 0) & (narrowed == 0))
            if bad.any():
                row = int(np.argwhere(bad)[0][0])
                return f"not representable as float32 at cell {row}"
        if self.normalize_embeddings:
            spans = selected.astype(np.float64).max(axis=0) - selected.astype(np.float64).min(axis=0)
            if float(spans.max()) <= 0:
                return "no coordinate variation"
        return None

    def _resolve_embeddings(self, requested: dict[str, tuple[str, list[int] | None]]) -> None:
        """Resolve one embedding per dimension slot, then order and index them.

        Explicitly requested slots are validated and must load. Slots not named
        follow the default rule, and its keys must load too: a caller who wrote
        ``X_umap_2d`` meant it. Only when neither fills any slot is the best
        guess made, and a guessed candidate that fails to load is skipped for
        the next one. When nothing loads at all, the error lists what ``obsm``
        holds.
        """
        n_cells = self.n_cells
        described = describe_obsm(self.adata.obsm, n_obs=n_cells)
        entries = obsm_entries(described)
        entry_by_key: dict[str, ObsmEntry] = {entry.key: entry for entry in entries}

        slots: dict[str, dict[str, Any] | None] = default_embedding_slots(entries)
        for slot, (key, columns) in requested.items():
            entry = entry_by_key.get(key)
            if entry is None:
                raise ValueError(
                    f"embedding_keys[{slot!r}] names obsm[{key!r}], which this object does "
                    f"not have. obsm holds: {[item.key for item in entries]}."
                )
            if not entry.usable or entry.width is None:
                raise ValueError(f"embedding_keys[{slot!r}]: obsm[{key!r}] {entry.reason}.")
            dimension = slot_dimension(slot)
            try:
                resolved_columns = validate_embedding_columns(
                    default_columns(dimension) if columns is None else columns,
                    dimension=dimension,
                    width=entry.width,
                )
            except DataMappingError as error:
                raise ValueError(f"embedding_keys[{slot!r}] obsm[{key!r}]: {error}") from error
            slots[slot] = {"key": key, "columns": resolved_columns}

        guessed: list[str] = []
        if all(ref is None for ref in slots.values()):
            guess = best_guess_embedding_slots(
                entries,
                latent_key=self.latent_key,
                load_check=self._embedding_load_defect,
            )
            slots = guess.slots
            guessed = guess.guessed
            if all(ref is None for ref in slots.values()):
                rejected = dict(guess.rejected)
                descriptions = [
                    describe_obsm_entry(
                        key=item["key"],
                        storage=item["storage"],
                        shape=item["shape"],
                        defect=obsm_entry_defect(
                            item, latent_key=self.latent_key, rejected=rejected
                        ),
                    )
                    for item in described
                ]
                raise ValueError(
                    nothing_plottable_message(n_obs=n_cells, descriptions=descriptions)
                )
        else:
            for slot, ref in slots.items():
                if ref is None:
                    continue
                entry = entry_by_key[ref["key"]]
                if not entry.usable or entry.width is None:
                    raise ValueError(f"Embedding {ref['key']!r} {entry.reason}.")
                dimension = slot_dimension(slot)
                if (
                    slot not in requested
                    and ref["key"] in DIMENSIONAL_UMAP_KEYS
                    and entry.width != dimension
                ):
                    # X_umap_<d>d names its own dimension, so a different width
                    # is a malformed declaration, not a wider array to slice.
                    raise ValueError(
                        f"Embedding {ref['key']!r} must have exactly {dimension} columns; "
                        f"got {entry.width}."
                    )
                try:
                    validate_embedding_columns(
                        ref["columns"], dimension=slot_dimension(slot), width=entry.width
                    )
                except DataMappingError as error:
                    raise ValueError(f"Embedding {ref['key']!r}: {error}") from error
                defect = self._embedding_load_defect(ref["key"], ref["columns"])
                if defect is not None:
                    # The full validator names the count and the first positions.
                    self._obsm_source(ref["key"], ref["columns"])
                    raise ValueError(f"Embedding {ref['key']!r} cannot be plotted: {defect}.")

        resolved: list[tuple[int, str, str, tuple[int, ...], int]] = []
        for slot, ref in slots.items():
            if ref is None:
                continue
            width = entry_by_key[ref["key"]].width
            assert width is not None
            embedding_id = canonical_embedding_id(ref["key"], ref["columns"], width=width)
            resolved.append(
                (slot_dimension(slot), embedding_id, ref["key"], tuple(ref["columns"]), width)
            )
        resolved.sort(key=lambda item: (item[0], item[1]))
        self._embeddings = [
            ServedEmbedding(
                id=embedding_id,
                index=index,
                key=key,
                columns=columns,
                dimension=dimension,
                width=width,
            )
            for index, (dimension, embedding_id, key, columns, width) in enumerate(resolved)
        ]
        self._embedding_defaults = {
            f"{item.dimension}d": item.id for item in self._embeddings
        }
        self._guessed_slots = guessed

    @property
    def embeddings(self) -> list[ServedEmbedding]:
        """The served embeddings, in (dimension, id) order."""
        self._check_closed()
        return list(self._embeddings)

    def _served_embedding(self, embedding_id: object) -> ServedEmbedding:
        embedding_id = _require_nonempty_string(embedding_id, label="embedding_id")
        for item in self._embeddings:
            if item.id == embedding_id:
                return item
        raise ValueError(
            f"Embedding {embedding_id!r} is not served. Served: {[item.id for item in self._embeddings]}."
        )

    def get_data_mapping(self) -> dict[str, Any]:
        """The data mapping this adapter serves, as a canonical DataMapping.

        The embeddings are the served slots; the graph, gene names, expression
        matrix and latent space are the ones this adapter was configured with;
        the vector bindings are the ``obsm`` keys it serves, each bound to its
        embedding's ``obsm`` key and columns.
        """
        self._check_closed()
        slots: dict[str, dict[str, Any] | None] = dict.fromkeys(DIMENSION_SLOTS)
        for item in self._embeddings:
            slots[f"{item.dimension}d"] = {"key": item.key, "columns": list(item.columns)}
        by_id = {item.id: item for item in self._embeddings}
        served = sorted(
            (
                (by_id[embedding_id].key, tuple(by_id[embedding_id].columns), obsm_key)
                for field_bindings in self._vector_field_obsm_keys.values()
                for embedding_id, obsm_key in field_bindings.items()
            ),
        )
        bindings = [
            {"key": obsm_key, "embedding": {"key": key, "columns": list(columns)}}
            for key, columns, obsm_key in served
        ]
        return validate_data_mapping(
            {
                "kind": "anndata",
                "embeddings": slots,
                "knn": self.knn_key,
                "geneNames": "index"
                if self.gene_id_column is None
                else f"column:{self.gene_id_column}",
                "expression": None
                if self._expression_matrix is None
                else self.expression_source,
                "latent": self.latent_key,
                "vectorFields": bindings,
            },
            label="served data mapping",
        )

    def get_mapping_candidates(self) -> dict[str, Any]:
        """The data-mapping catalog of the served object.

        Every ``obsm``, ``obsp``, layer and ``var`` column, with shapes, dtypes
        and whether each can be used, read from metadata alone; ``defaults`` is
        the mapping this adapter serves (:meth:`get_data_mapping`), and
        ``guessed`` names the slots that were best-guessed. The catalog is built
        once and then frozen; this returns a copy.
        """
        self._check_closed()
        with self._mapping_candidates_lock:
            if self._mapping_candidates is None:
                catalog = build_anndata_candidates(
                    self.adata,
                    latent_key=self.latent_key,
                    knn_key=self.knn_key,
                )
                catalog["defaults"] = self.get_data_mapping()
                catalog["guessed"] = list(self._guessed_slots)
                self._mapping_candidates = catalog
            return cast(dict[str, Any], json.loads(json.dumps(self._mapping_candidates)))

    def get_embedding_id_for_payload(self, component: str) -> str | None:
        """Resolve one ``{index}_{dimension}d.bin`` route to its served embedding id."""
        self._check_closed()
        match = re.fullmatch(r"(\d+)_([123])d\.bin", component)
        if match is None:
            return None
        index = _payload_route_index(match.group(1))
        if index is None or not 0 <= index < len(self._embeddings):
            return None
        item = self._embeddings[index]
        if item.dimension != int(match.group(2)):
            return None
        return item.id

    def _detect_vector_fields(
        self,
        *,
        vector_field_default: object,
    ) -> dict[str, Any] | None:
        """Bind and validate every ``<field>_umap_<d>d`` vector by the default rule."""
        served_keys = {item.key for item in self._embeddings}
        keys: list[str] = []
        for key in self.adata.obsm:
            if not isinstance(key, str):
                raise TypeError("AnnData obsm keys must be native strings.")
            if not key:
                raise ValueError("AnnData obsm keys must be non-empty strings.")
            # An embedding is never also a vector field: X_umap_2d matches the
            # vector grammar, and so does any served key that happens to.
            if key in DIMENSIONAL_UMAP_KEYS or key in served_keys:
                continue
            keys.append(key)
        slots = {f"{item.dimension}d": {"key": item.key, "columns": list(item.columns)} for item in self._embeddings}
        id_by_slot = {f"{item.dimension}d": item.id for item in self._embeddings}
        bindings = default_vector_bindings(keys, slots)
        values: dict[str, dict[str, Any]] = {}
        obsm_keys: dict[str, dict[str, str]] = {}
        for binding in bindings:
            match = VECTOR_KEY_PATTERN.fullmatch(binding["key"])
            assert match is not None
            field_id = match.group("field")
            embedding_id = id_by_slot[f"{match.group('dimension')}d"]
            values.setdefault(field_id, {})[embedding_id] = self.adata.obsm[binding["key"]]
            obsm_keys.setdefault(field_id, {})[embedding_id] = binding["key"]

        validated = validate_vector_fields(
            values,
            n_cells=self.n_cells,
            embedding_dimensions={item.id: item.dimension for item in self._embeddings},
            vector_field_default=vector_field_default,
            label_for=lambda field_id, embedding_id: (
                f"Vector field {obsm_keys[field_id][embedding_id]!r}"
            ),
        )
        if not validated.fields:
            # Asking for vector fields asks for whatever ``obsm`` declares and
            # the default rule binds, and declaring none is an ordinary result:
            # a key outside the ``<field>_umap_<n>d`` grammar is not a vector
            # field, and a UMAP vector does not bind to a slot showing another
            # layout. A malformed declaration still fails, in validation above.
            return None

        field_ids = _require_field_identities(
            list(validated.fields),
            label="Vector field",
        )
        self._vector_field_obsm_keys = {
            field_id: {embedding_id: obsm_keys[field_id][embedding_id] for embedding_id in by}
            for field_id, by in validated.fields.items()
        }
        self._vector_field_payload_index_by_id = {
            field_id: index for index, field_id in enumerate(field_ids)
        }
        index_by_embedding = {item.id: item.index for item in self._embeddings}
        fields: dict[str, dict[str, Any]] = {}
        for field_id, by_embedding in validated.fields.items():
            payload_index = self._vector_field_payload_index_by_id[field_id]
            # Same key order as the prepared export and as cellucid-r.
            fields[field_id] = {
                "label": field_id,
                "files": {
                    embedding_id: f"vectors/{payload_index}_{index_by_embedding[embedding_id]}.bin"
                    for embedding_id in by_embedding
                },
            }
        return {
            "default_field": validated.default_field,
            "fields": fields,
        }

    # =========================================================================
    # EMBEDDING DATA
    # =========================================================================

    def get_embedding(self, embedding_id: str) -> np.ndarray:
        """
        Get the normalized coordinates of one served embedding.

        Returns a float64 array of shape (n_cells, dimension). The payload is
        float32, but the rounding happens once, in ``get_embedding_binary()``,
        so a centroid measured from these coordinates is not measured from
        already-rounded ones.
        """
        self._check_closed()
        item = self._served_embedding(embedding_id)
        if item.id in self._embedding_cache:
            return self._embedding_cache[item.id]

        if item.key not in self.adata.obsm:
            raise ValueError(f"Embedding source obsm[{item.key!r}] no longer exists.")
        # Precision-preserving, like the prepared export: the extents are
        # measured at the input's own precision and the coordinates are
        # rounded to float32 once, in get_embedding_binary().
        source_shape = tuple(int(extent) for extent in np.shape(self.adata.obsm[item.key]))
        expected_shape = (self.n_cells,) if item.width == 1 and len(source_shape) == 1 else (
            self.n_cells,
            item.width,
        )
        if source_shape != expected_shape:
            raise ValueError(
                f"Embedding {item.key!r} changed after adapter initialization; "
                f"expected shape {expected_shape}, got {source_shape}."
            )
        raw = self._obsm_source(item.key, item.columns)

        if self.normalize_embeddings:
            raw, center, scale_factor, _max_range = _normalize_finite_float32_embedding(
                raw,
                label=f"Embedding {item.key!r}",
            )
            self._embedding_norm[item.id] = {
                "center": center,
                "scale_factor": float(scale_factor),
            }
        else:
            self._embedding_norm[item.id] = {
                "center": np.zeros((item.dimension,), dtype=np.float64),
                "scale_factor": 1.0,
            }

        self._embedding_cache[item.id] = raw
        return raw

    def get_embedding_binary(self, embedding_id: str, compress: bool = False) -> bytes:
        """Get one served embedding as little-endian float32 bytes."""
        compress = _require_native_boolean(compress, label="compress")
        embedding = self.get_embedding(embedding_id)
        data = little_endian_payload_bytes(embedding.astype(np.float32))

        if compress:
            return deterministic_gzip_compress(data, compresslevel=6)
        return data

    def get_vector_field_binary(
        self, field_id: str, embedding_id: str, compress: bool = False
    ) -> bytes:
        """
        Get a per-cell vector field (displacement vectors) as binary float32 data.

        Vector fields are scaled by the normalization factor of the embedding
        they are bound to, so they are in the same normalized space as that
        embedding's points.
        """
        self._check_closed()

        field = _require_nonempty_string(field_id, label="field_id")
        embedding = _require_nonempty_string(embedding_id, label="embedding_id")
        compress = _require_native_boolean(compress, label="compress")

        if self._vector_fields_metadata is None:
            raise ValueError(f"Vector field {field!r} not available.")
        fields = self._vector_fields_metadata["fields"]
        if field not in fields:
            raise ValueError(f"Vector field '{field}' not available")

        source_keys = self._vector_field_obsm_keys[field]
        if embedding not in source_keys:
            raise ValueError(f"Vector field {field!r} is not bound to embedding {embedding!r}.")
        obsm_key = source_keys[embedding]
        if obsm_key not in self.adata.obsm:
            raise ValueError(f"Vector field source {obsm_key!r} no longer exists.")
        item = self._served_embedding(embedding)

        cache_key = (field, embedding)
        if cache_key not in self._vector_field_cache:
            # The same validator the prepared export uses, so the served
            # bytes are the prepared bytes: the caller's own values are scaled
            # and rounded to float32 exactly once, in scale_vector_field().
            raw, _dimension = _finite_float32_matrix(
                self.adata.obsm[obsm_key],
                label=f"Vector field {obsm_key!r}",
                n_rows=self.n_cells,
                declared_dimension=item.dimension,
            )
            self.get_embedding(embedding)
            if embedding not in self._embedding_norm:
                raise RuntimeError(f"Missing normalization state for embedding {embedding!r}.")
            scale_factor = self._embedding_norm[embedding]["scale_factor"]
            self._vector_field_cache[cache_key] = scale_vector_field(
                raw,
                scale_factor=scale_factor,
                label=f"Vector field {obsm_key!r}",
            )

        data = little_endian_payload_bytes(self._vector_field_cache[cache_key])
        if compress:
            return deterministic_gzip_compress(data, compresslevel=6)
        return data

    def get_vector_field_route(self, component: str) -> tuple[str, str] | None:
        """Resolve one ``{fieldIndex}_{embeddingIndex}.bin`` route to (field, embedding)."""
        self._check_closed()
        match = re.fullmatch(r"(\d+)_(\d+)\.bin", component)
        if match is None or self._vector_fields_metadata is None:
            return None
        field_index = _payload_route_index(match.group(1))
        embedding_index = _payload_route_index(match.group(2))
        if field_index is None or embedding_index is None:
            return None
        field_id = next(
            (
                field
                for field, index in self._vector_field_payload_index_by_id.items()
                if index == field_index
            ),
            None,
        )
        if field_id is None or not 0 <= embedding_index < len(self._embeddings):
            return None
        embedding_id = self._embeddings[embedding_index].id
        if embedding_id not in self._vector_field_obsm_keys.get(field_id, {}):
            return None
        return field_id, embedding_id

    # =========================================================================
    # LATENT SPACE FOR OUTLIER COMPUTATION
    # =========================================================================

    def _get_latent_space(self) -> np.ndarray | None:
        """Get latent space for outlier quantile computation."""
        if self._latent_space is not None:
            return self._latent_space

        if self.latent_key is None or self.latent_key not in self.adata.obsm:
            return None

        latent = _require_finite_float32_array(
            self.adata.obsm[self.latent_key],
            label=f"Latent field {self.latent_key!r}",
        )
        if latent.ndim != 2 or latent.shape[0] != self.n_cells:
            raise ValueError(
                f"Latent field {self.latent_key!r} must have shape "
                f"({self.n_cells}, n_dimensions), got {latent.shape}."
            )
        self._latent_space = latent
        return self._latent_space

    # =========================================================================
    # OBS DATA (CELL METADATA)
    # =========================================================================

    def get_obs_keys(self) -> list[str]:
        """Get the exact obs column names this adapter serves."""
        self._check_closed()
        return list(self._obs_keys)

    def _obs_series(self, key: str) -> pd.Series:
        """Return one served obs column, naming the served set when it is absent."""
        if key not in self._obs_payload_index_by_key:
            raise KeyError(f"obs field '{key}' not found. Available fields: {self._obs_keys}")
        return cast(pd.Series, self.adata.obs[key])

    def get_obs_payload_index(self, key: str) -> int:
        """Return the exact payload index one obs key is served under."""
        self._check_closed()
        try:
            return self._obs_payload_index_by_key[key]
        except KeyError as error:
            raise KeyError(f"Observation field {key!r} has no payload index.") from error

    def get_obs_key_for_payload_index(self, component: str) -> str | None:
        """Resolve one payload-index route component to its exact obs key."""
        self._check_closed()
        index = _payload_route_index(component)
        if index is None or not 0 <= index < len(self._obs_keys):
            return None
        return self._obs_keys[index]

    def get_obs_field_kind(self, key: str) -> Literal["continuous", "category"]:
        """
        Determine if an obs field is continuous or categorical.

        Classification rules:
        - Categorical dtype → category
        - Boolean dtype → category
        - Numeric dtype → continuous
        - String/object → category (treated as labels)
        - Empty column → category (safe default)
        """
        self._check_closed()
        return _obs_field_kind(self._obs_series(key), key=key)

    def get_obs_continuous_values(self, key: str, compress: bool = False) -> bytes:
        """
        Get continuous obs field as binary float32 data.

        Values must remain finite after float32 representation.

        Raises
        ------
        KeyError
            If the field is not found in obs.
        """
        self._check_closed()
        compress = _require_native_boolean(compress, label="compress")

        column = self._obs_series(key).to_numpy()
        try:
            values = _require_continuous_obs_values(
                column,
                key=key,
                n_cells=self.n_cells,
            )
        except ValueError as exc:
            diagnosis = diagnose_continuous_payload(
                column, kind="field", name=key
            )
            if diagnosis is None:
                raise
            raise NonFinitePayloadError(diagnosis) from exc

        data = little_endian_payload_bytes(values)

        if compress:
            return deterministic_gzip_compress(data, compresslevel=6)
        return data

    def get_obs_categorical_codes(
        self,
        key: str,
        compress: bool = False,
    ) -> tuple[bytes, list[str | bool | int | float], int]:
        """
        Get categorical obs field as binary codes.

        Returns:
            (binary_codes, category_list, missing_value)

        Categories are assigned codes 0 to n-1.
        Missing values (NaN) are encoded as the missing_value sentinel.

        Raises
        ------
        KeyError
            If the field is not found in obs.
        """
        self._check_closed()
        compress = _require_native_boolean(compress, label="compress")

        s = self._obs_series(key)
        cat = s.astype("category")
        categories = _json_category_values(
            cat.cat.categories,
            field_key=key,
        )
        codes = cat.cat.codes.to_numpy(dtype=np.int32)

        n_categories = len(categories)
        dtype, _dtype_str, missing_value = _categorical_storage(
            n_categories,
            field_key=key,
        )

        # Convert codes (-1 for NaN -> missing_value)
        codes_typed = np.full(self.n_cells, missing_value, dtype=dtype)
        valid_mask = codes >= 0
        codes_typed[valid_mask] = codes[valid_mask].astype(dtype)

        # Log missing value count
        n_missing = (~valid_mask).sum()
        if n_missing > 0:
            logger.debug(f"obs field '{key}': {n_missing} missing values out of {self.n_cells}")

        data = little_endian_payload_bytes(codes_typed)
        if compress:
            data = deterministic_gzip_compress(data, compresslevel=6)

        return data, categories, int(missing_value)

    # =========================================================================
    # CENTROIDS AND OUTLIER QUANTILES
    # =========================================================================

    def _compute_centroids_for_field(
        self,
        key: str,
        embedding_id: str,
    ) -> list[dict]:
        """Compute centroids for a categorical field on one served embedding."""
        cache_key = (key, embedding_id)
        if cache_key in self._centroid_cache:
            return self._centroid_cache[cache_key]

        coords = self.get_embedding(embedding_id)
        s = self._obs_series(key)
        cat = s.astype("category")
        categories = _json_category_values(
            cat.cat.categories,
            field_key=key,
        )
        codes = cat.cat.codes.to_numpy(dtype=np.int32)

        centroids: list[dict[str, Any]] = []
        for code, label in enumerate(categories):
            mask = codes == code
            idx = np.nonzero(mask)[0]
            n = idx.size

            if n < self.centroid_min_points:
                continue

            pts = coords[idx, :]
            center = pts.mean(axis=0)

            # Remove outliers for centroid computation
            if n > self.centroid_min_points:
                dists = np.linalg.norm(pts - center, axis=1)
                thr = float(np.quantile(dists, self.centroid_outlier_quantile))
                inlier_mask = dists <= thr
                n_in = int(inlier_mask.sum())
                if n_in >= self.centroid_min_points:
                    pts_in = pts[inlier_mask, :]
                    center = pts_in.mean(axis=0)
                    used_count = n_in
                else:
                    used_count = n
            else:
                used_count = n

            centroids.append(
                {
                    "category": label,
                    "position": center.astype(float).tolist(),
                    "n_points": int(used_count),
                }
            )

        self._centroid_cache[cache_key] = centroids
        return centroids

    def get_centroids_for_field(self, key: str) -> dict[str, list[dict]]:
        """Get centroids on every served embedding, keyed by embedding id."""
        return {
            item.id: self._compute_centroids_for_field(key, item.id) for item in self._embeddings
        }

    def _compute_outlier_quantiles(self, key: str) -> np.ndarray:
        """Compute per-cell outlier quantiles for a categorical field."""
        if key in self._outlier_quantile_cache:
            return self._outlier_quantile_cache[key]

        latent = self._get_latent_space()
        if latent is None:
            raise ValueError(
                f"Cannot compute outlier quantiles for {key!r}: "
                "the dataset has no declared latent space"
            )

        s = self._obs_series(key)
        cat = s.astype("category")
        categories = _json_category_values(
            cat.cat.categories,
            field_key=key,
        )
        codes = cat.cat.codes.to_numpy(dtype=np.int32)

        quantiles = np.full(self.n_cells, np.nan, dtype=np.float32)

        for code, _label in enumerate(categories):
            mask = codes == code
            idx = np.nonzero(mask)[0]
            n = idx.size

            if n < self.centroid_min_points:
                continue

            pts = latent[idx, :]
            centroid = pts.mean(axis=0)
            dists = np.linalg.norm(pts - centroid, axis=1)
            sorted_dists = np.sort(dists)
            ranks = np.searchsorted(sorted_dists, dists, side="right")
            cell_quantiles = ranks.astype(np.float32) / n
            quantiles[idx] = cell_quantiles

        self._outlier_quantile_cache[key] = quantiles
        return quantiles

    def get_obs_outlier_quantiles(self, key: str, compress: bool = False) -> bytes:
        """Get outlier quantiles as binary float32 data."""
        self._check_closed()
        compress = _require_native_boolean(compress, label="compress")
        quantiles = self._compute_outlier_quantiles(key)
        data = little_endian_payload_bytes(quantiles)

        if compress:
            return deterministic_gzip_compress(data, compresslevel=6)
        return data

    # =========================================================================
    # VAR DATA (GENE EXPRESSION)
    # =========================================================================

    def get_gene_ids(self) -> list[str]:
        """Get the served gene names: the gene axis's names, made unique."""
        self._check_closed()
        if self._gene_ids_cache is None:
            raise RuntimeError("Gene names were not resolved.")
        return self._gene_ids_cache

    def get_gene_payload_index(self, gene_id: str) -> int:
        """Return the exact payload index one gene is served under."""
        self._check_closed()
        try:
            return self._gene_payload_index_by_id[gene_id]
        except KeyError as error:
            raise KeyError(f"Gene {gene_id!r} has no payload index.") from error

    def get_gene_id_for_payload_index(self, component: str) -> str | None:
        """Resolve one payload-index route component to its exact gene identifier."""
        self._check_closed()
        gene_ids = self.get_gene_ids()
        index = _payload_route_index(component)
        if index is None or not 0 <= index < len(gene_ids):
            return None
        return gene_ids[index]

    def _get_gene_idx(self, gene_id: str) -> int:
        """Get gene index with O(1) lookup."""
        # Ensure cache is populated
        if self._gene_id_to_idx_cache is None:
            self.get_gene_ids()
        gene_id_to_idx = self._gene_id_to_idx_cache
        if gene_id_to_idx is None:
            raise RuntimeError("Gene identifier index was not initialized")

        if gene_id not in gene_id_to_idx:
            raise KeyError(f"Gene '{gene_id}' not found in var")

        return gene_id_to_idx[gene_id]

    def _get_gene_column(self, gene_idx: int) -> np.ndarray:
        """
        Extract a single gene column from the expression matrix.

        Handles multiple data formats efficiently:
        - Dense numpy arrays: Direct column access O(1)
        - CSC sparse matrices: Efficient column access O(nnz/n_cols)
        - CSR sparse matrices: Convert to CSC and cache for repeated access
        - Other sparse formats (COO, DOK, etc.): Convert to CSC and cache
        - None: Rejected because no expression payload exists
        - Backed/chunked arrays: Direct slicing (no caching to preserve lazy loading)

        Memory-Speed Tradeoff:
        For in-memory CSR matrices, the CSC cache doubles memory usage but enables
        O(nnz/n_cols) column access instead of O(nnz). This pays off after
        accessing more than a few genes.

        IMPORTANT: For backed h5ad files, we do NOT create a CSC cache because:
        1. It would force loading the entire X matrix into memory
        2. It defeats the purpose of lazy loading for large datasets
        Instead, we accept slower column access to preserve memory efficiency.

        Parameters
        ----------
        gene_idx : int
            Column index of the gene to extract.

        Returns
        -------
        np.ndarray
            1D array of expression values for all cells.
        """
        self._check_closed()

        X = self._expression_matrix

        if X is None:
            raise ValueError("AnnData has no X expression matrix; gene expression is unavailable")

        # Edge case: invalid gene index
        if gene_idx < 0 or gene_idx >= self.n_genes:
            raise IndexError(f"Gene index {gene_idx} out of range [0, {self.n_genes})")

        # CRITICAL: For backed h5ad files, avoid creating CSC cache to preserve
        # lazy loading behavior. The memory cost of CSC cache would defeat the
        # purpose of backed mode for large datasets.
        if self._source_type == "h5ad" and self.is_backed:
            # For backed h5ad, use direct slicing (slower but memory-efficient)
            # Note: backed h5ad matrices support slicing but may be slower
            col = X[:, gene_idx]
            if sparse.issparse(col):
                return np.asarray(col.toarray()).flatten()
            return np.asarray(col).flatten()

        if sparse.issparse(X):
            # For sparse matrix, CSC is efficient for column access
            if sparse.isspmatrix_csc(X):
                return np.asarray(X.getcol(gene_idx).toarray()).flatten()
            elif sparse.isspmatrix_csr(X):
                # CSR is inefficient for column access - convert to CSC and cache
                # This is a one-time cost that pays off with repeated gene queries
                # Only do this for in-memory data (not backed files)
                if self._X_csc_cache is None:
                    logger.info(
                        f"Converting CSR matrix ({X.shape[0]:,}x{X.shape[1]:,}, "
                        f"{X.nnz:,} non-zeros) to CSC for efficient column access"
                    )
                    self._X_csc_cache = X.tocsc()
                return np.asarray(self._X_csc_cache.getcol(gene_idx).toarray()).flatten()
            else:
                # Other sparse format (COO, LIL, DOK, BSR) - convert to CSC
                if self._X_csc_cache is None:
                    logger.info(f"Converting {type(X).__name__} to CSC for column access")
                    self._X_csc_cache = sparse.csc_matrix(X)
                return np.asarray(self._X_csc_cache.getcol(gene_idx).toarray()).flatten()
        else:
            # Dense matrix - direct column access
            return np.asarray(X[:, gene_idx]).flatten()

    def _get_validated_gene_values(self, gene_id: str) -> np.ndarray:
        """Return one exact finite Float32 per-cell gene vector."""
        if self._expression_matrix is None:
            raise ValueError("AnnData has no X expression matrix; gene expression is unavailable")
        values = self._gene_expression_cache.get(gene_id)
        if values is None:
            gene_idx = self._get_gene_idx(gene_id)
            col = self._get_gene_column(gene_idx)
            try:
                values = _require_finite_float32_array(
                    col,
                    label=f"Gene {gene_id!r} expression",
                )
            except ValueError as exc:
                # Counting is only worth a pass over the column once the cheap
                # all-finite check has already said the column is unusable.
                diagnosis = diagnose_continuous_payload(
                    col, kind="gene", name=gene_id
                )
                if diagnosis is None:
                    raise
                raise NonFinitePayloadError(diagnosis) from exc
            if values.ndim != 1 or values.shape[0] != self.n_cells:
                raise ValueError(
                    f"Gene {gene_id!r} expression must have shape ({self.n_cells},), "
                    f"got {values.shape}."
                )
            self._gene_expression_cache.put(gene_id, values)
        return values

    def get_gene_expression(self, gene_id: str, compress: bool = False) -> bytes:
        """Get expression values for a single gene as binary float32."""
        self._check_closed()
        compress = _require_native_boolean(compress, label="compress")
        values = self._get_validated_gene_values(gene_id)
        data = little_endian_payload_bytes(values)
        if compress:
            return deterministic_gzip_compress(data, compresslevel=6)
        return data

    def get_gene_min_max(self, gene_id: str) -> tuple[float, float]:
        """Get min/max values for a gene (for colormap scaling)."""
        values = self._get_validated_gene_values(gene_id)
        return float(values.min()), float(values.max())

    # =========================================================================
    # CONNECTIVITY DATA
    # =========================================================================

    def _initialize_connectivity(self) -> None:
        """Read and check the neighbour graph ``knn_key`` names, atomically."""
        if self.knn_key is None:
            # Not reading obsp at all is the point: on a large object the edge
            # list is the most expensive thing startup could do, and a session
            # that never draws the graph should never pay for it.
            self._has_connectivity = False
            self._connectivity_manifest = None
            return

        key = self.knn_key
        try:
            present = key in self.adata.obsp
        except Exception as exc:
            raise ValueError(f"Could not inspect adata.obsp for {key!r}") from exc
        if not present:
            raise ValueError(
                f"knn_key names obsp[{key!r}], which this object does not have. "
                f"obsp holds: {list(self.adata.obsp.keys())}. Compute the neighbour "
                "graph with sc.pp.neighbors(adata), or serve without a graph."
            )
        self._has_connectivity = True

        self._report("Connectivity", f"reading obsp[{key!r}]")
        try:
            edges = self._compute_connectivity_edges()
        except Exception as exc:
            self._has_connectivity = False
            reason = str(exc)
            if "symmetric" in reason:
                raise ValueError(
                    f"obsp[{key!r}] is not symmetric \u2014 choose a connectivities matrix. "
                    f"The neighbour graph must be undirected: {reason}"
                ) from exc
            raise ValueError(f"obsp[{key!r}] cannot be served as the neighbour graph: {reason}") from exc

        self._connectivity_manifest = build_connectivity_manifest(
            n_cells=self.n_cells,
            n_edges=edges.n_edges,
            max_neighbors=edges.max_neighbors,
            index_bytes=edges.index_bytes,
            index_dtype=edges.index_dtype,
            compression=None,
        )

    def has_connectivity(self) -> bool:
        """Report whether this adapter serves a connectivity graph.

        This answers what the dataset publishes, not what ``adata.obsp`` holds:
        an object whose graph was not asked for reports False, exactly as an
        object with no graph does, because the identity, the manifest, and the
        edge routes all follow this one answer.
        """
        self._check_closed()
        return self._has_connectivity

    def _compute_connectivity_edges(self) -> ConnectivityEdgePairs:
        """Validate and cache one aligned weighted edge payload."""
        self._check_closed()

        if self._connectivity_cache is not None:
            return self._connectivity_cache

        if not self.has_connectivity() or self.knn_key is None:
            raise ValueError("This adapter serves no neighbour graph.")

        self._connectivity_cache = validate_connectivity_edges(
            self.adata.obsp[self.knn_key],
            n_cells=self.n_cells,
        )
        return self._connectivity_cache

    def get_connectivity_edges(
        self,
        compress: bool = False,
    ) -> tuple[bytes, bytes, bytes, int, int]:
        """
        Get aligned endpoint and Float64 weight payloads.

        Returns:
            (sources_binary, destinations_binary, weights_binary,
            n_edges, max_neighbors)
        """
        self._check_closed()
        compress = _require_native_boolean(compress, label="compress")
        edges = self._compute_connectivity_edges()

        sources_data = little_endian_payload_bytes(edges.sources)
        destinations_data = little_endian_payload_bytes(edges.destinations)
        weights_data = little_endian_payload_bytes(edges.weights)

        if compress:
            sources_data = deterministic_gzip_compress(
                sources_data,
                compresslevel=6,
            )
            destinations_data = deterministic_gzip_compress(
                destinations_data,
                compresslevel=6,
            )
            weights_data = deterministic_gzip_compress(
                weights_data,
                compresslevel=6,
            )

        return (
            sources_data,
            destinations_data,
            weights_data,
            edges.n_edges,
            edges.max_neighbors,
        )

    # =========================================================================
    # MANIFEST GENERATION
    # =========================================================================

    def get_dataset_identity(self) -> dict:
        """Generate dataset_identity.json content."""
        obs_manifest = self.get_obs_manifest()
        continuous_fields = obs_manifest["_continuousFields"]
        categorical_fields = obs_manifest["_categoricalFields"]
        obs_fields = [{"key": field[1], "kind": "continuous"} for field in continuous_fields] + [
            {
                "key": field[1],
                "kind": "category",
                "n_categories": len(field[2]),
            }
            for field in categorical_fields
        ]
        n_continuous = len(continuous_fields)
        n_categorical = len(categorical_fields)

        has_conn = self.has_connectivity()
        n_edges = (
            self._connectivity_manifest["n_edges"]
            if self._connectivity_manifest is not None
            else None
        )

        # Build embeddings metadata: one item per served slot.
        embeddings_meta = {
            "default_dimension": self._default_dimension,
            "defaults": dict(self._embedding_defaults),
            "items": {
                item.id: {
                    "label": item.id,
                    "dimension": item.dimension,
                    "file": item.file,
                    "source": {"key": item.key, "columns": list(item.columns)},
                }
                for item in self._embeddings
            },
        }
        identity = {
            "version": IDENTITY_VERSION,
            "id": self.dataset_id,
            "name": self.dataset_name,
            "description": "Loaded directly from AnnData",
            "created_at": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "cellucid_data_version": "anndata_adapter",
            "stats": {
                "n_cells": self.n_cells,
                "n_genes": 0 if self._expression_matrix is None else self.n_genes,
                "n_obs_fields": len(obs_fields),
                "n_categorical_fields": n_categorical,
                "n_continuous_fields": n_continuous,
                "has_connectivity": has_conn,
                "n_edges": n_edges,
            },
            "embeddings": embeddings_meta,
            "obs_fields": obs_fields,
            # No export_settings. That key records the choices a writer made
            # while producing an export; a dataset served straight from an
            # AnnData object is not an export and made none. It is optional in
            # the format, nothing reads it back, and stating it forced a
            # placeholder dtype that no export can contain. Each categorical
            # field declares its exact uint8/uint16 storage in obs_manifest.json
            # according to its category count.
            "source": {
                "name": {
                    "h5ad": "H5AD file",
                    "zarr": "Zarr store",
                    "memory": "In-memory AnnData",
                }[self._source_type],
            },
        }

        if self._vector_fields_metadata:
            identity["vector_fields"] = self._vector_fields_metadata

        return identity

    def get_obs_manifest(self) -> dict:
        """Generate obs_manifest.json content."""
        obs_keys = self.get_obs_keys()

        continuous_fields = []
        categorical_fields = []

        for key in obs_keys:
            payload_index = self.get_obs_payload_index(key)
            kind = self.get_obs_field_kind(key)

            if kind == "continuous":
                continuous_fields.append([payload_index, key])
            else:
                # Categorical
                s = self._obs_series(key)
                cat = s.astype("category")
                categories = _json_category_values(
                    cat.cat.categories,
                    field_key=key,
                )
                n_categories = len(categories)

                _dtype, dtype_str, missing_value = _categorical_storage(
                    n_categories,
                    field_key=key,
                )

                # One centroid list per served embedding, keyed by its id.
                centroids_by_embedding = self.get_centroids_for_field(key)

                categorical_fields.append(
                    [payload_index, key, categories, dtype_str, missing_value, centroids_by_embedding]
                )

        # Build compact manifest format
        # Path patterns must match prepare format for consistency
        # Note: We don't use compression (.gz) in adapter mode - gzip is applied
        # at the HTTP level when client supports Accept-Encoding: gzip
        obs_schemas = {}
        if continuous_fields:
            obs_schemas["continuous"] = {
                "pathPattern": "obs/{index}.values.f32",
                "ext": "f32",
                "dtype": "float32",
                "quantized": False,
            }
        if categorical_fields:
            if self.latent_key is None:
                obs_schemas["categorical"] = {
                    "codesPathPattern": "obs/{index}.codes.{ext}",
                    "outlierPathPattern": None,
                    "outlierExt": None,
                    "outlierDtype": None,
                    "outlierQuantized": False,
                }
            else:
                obs_schemas["categorical"] = {
                    "codesPathPattern": "obs/{index}.codes.{ext}",
                    "outlierPathPattern": "obs/{index}.outliers.f32",
                    "outlierExt": "f32",
                    "outlierDtype": "float32",
                    "outlierQuantized": False,
                }

        return {
            "_format": MANIFEST_FORMAT_VERSION,
            "n_points": self.n_cells,
            "centroid_outlier_quantile": self.centroid_outlier_quantile,
            "latent_key": self.latent_key,
            "compression": None,
            "_obsSchemas": obs_schemas,
            "_continuousFields": continuous_fields,
            "_categoricalFields": categorical_fields,
        }

    def get_var_manifest(self) -> dict:
        """Generate var_manifest.json content."""
        gene_name_set: dict[str, Any]
        if self._expression_matrix is None:
            fields = []
            gene_name_set = {
                "source": self.gene_id_column,
                "filled_from_index": 0,
                "renamed": [],
                "names": None,
            }
        else:
            fields = [
                [self.get_gene_payload_index(gene_id), gene_id] for gene_id in self.get_gene_ids()
            ]
            gene_name_set = {
                "source": self.gene_id_column,
                "filled_from_index": self._gene_names.filled_from_index,
                "renamed": [list(pair) for pair in self._gene_names.renamed],
                "names": None,
            }

        # Path pattern must match prepare format for consistency
        # Note: We don't use compression (.gz) in adapter mode - gzip is applied
        # at the HTTP level when client supports Accept-Encoding: gzip
        var_schema = {
            "kind": "continuous",
            "pathPattern": "var/{index}.values.f32",
            "ext": "f32",
            "dtype": "float32",
            "quantized": False,
        }

        return {
            "_format": MANIFEST_FORMAT_VERSION,
            "n_points": self.n_cells,
            "compression": None,
            "quantization": None,
            "_varSchema": var_schema,
            "fields": fields,
            "gene_name_sets": [gene_name_set],
        }

    def get_connectivity_manifest(self) -> dict | None:
        """Generate connectivity_manifest.json content."""
        self._check_closed()
        if self._connectivity_manifest is None:
            return None
        return dict(self._connectivity_manifest)

    # =========================================================================
    # CLEANUP AND CONTEXT MANAGER
    # =========================================================================

    def close(self) -> None:
        """
        Close the adapter and release all resources.

        This method:
        1. Clears all caches to free memory (embedding, centroid, CSC, gene expression)
        2. Closes the underlying file handle for backed h5ad files
        3. Marks the adapter as closed to prevent further operations

        Safe to call multiple times. Always call this method when done with the
        adapter, or use the context manager::

            with AnnDataAdapter.from_file(
                "data.h5ad",
                dataset_name="Example",
                dataset_id="example",
            ) as adapter:
                # use adapter
            # automatically cleaned up

        Memory Released
        ---------------

        - Embedding cache (normalized coordinates)
        - Centroid cache (computed label centroids)
        - Outlier quantile cache
        - Gene expression LRU cache (up to 100 gene columns)
        - CSC matrix cache (for CSR->CSC converted matrices)
        - Latent space array
        - Gene ID lookup indices
        """
        if self._close_complete:
            return

        # Capture the backed file before marking the adapter closed.
        file_to_close = self.adata.file if self.adata.isbacked else None

        self._closed = True

        # Clear all caches to free memory
        self._embedding_cache.clear()
        self._embedding_norm.clear()
        self._expression_matrix = None
        self._vector_field_cache.clear()
        self._vector_field_obsm_keys.clear()
        self._vector_field_payload_index_by_id.clear()
        self._vector_fields_metadata = None
        self._centroid_cache.clear()
        self._outlier_quantile_cache.clear()
        self._gene_expression_cache.clear()
        self._latent_space = None
        self._connectivity_cache = None
        self._connectivity_manifest = None
        self._has_connectivity = False
        self._gene_ids_cache = None
        self._gene_id_to_idx_cache = None
        self._gene_payload_index_by_id.clear()
        self._obs_payload_index_by_key.clear()

        # Clear CSC cache (can be large for sparse matrices).
        self._X_csc_cache = None

        # Close backed file handle if applicable
        if file_to_close is not None:
            file_to_close.close()
            logger.debug(f"Closed backed file handle for {self.dataset_name}")

        # Clear reference to adata to help garbage collection
        # (do this last since we needed it for is_backed check)
        self.adata = None
        self._close_complete = True

        logger.debug(f"AnnDataAdapter closed: {self.dataset_name}")

    def __enter__(self) -> AnnDataAdapter:
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit - ensures cleanup."""
        self.close()

    def __repr__(self) -> str:
        if self._closed:
            return "AnnDataAdapter(closed)"
        backed_str = " (backed)" if self.is_backed else ""
        source_str = f" from {self._source_type}" if self._source_type != "memory" else ""
        embeddings = ", ".join(item.id for item in self._embeddings)
        return (
            f"AnnDataAdapter({self.n_cells:,} cells, {self.n_genes:,} genes, "
            f"embeddings=[{embeddings}]{backed_str}{source_str})"
        )
