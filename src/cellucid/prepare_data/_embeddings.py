"""
The embeddings an export carries, and the vector fields bound to them.

An export holds any number of embeddings per dimension. Each one is an *item*
with an id, a display label, a dimension and the 0-based source columns it was
read from. Items are ordered by (dimension, id in code-point order); an item's
position in that order is its payload index, so its coordinates are written to
``embeddings/{index}_{dimension}d.bin``. One item per dimension is that
dimension's default, and the highest dimension is the export's initial one.

Every item is normalized on its own with the one rule the format defines, and a
vector field bound to an item is scaled by that item's factor, so its arrows
keep their length relative to the points they are drawn over.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, cast

import numpy as np
from scipy import sparse

from .._contracts import (
    EMBEDDINGS_DIRNAME,
    _require_display_text,
    _require_field_identities,
    _require_nonempty_string,
)
from ..data_mapping import (
    DIMENSION_SLOTS,
    DataMappingError,
    default_columns,
    dimension_slot,
    validate_embedding_columns,
)
from ..vector_fields import ValidatedVectorFields, scale_vector_field, validate_vector_fields
from ._arrays import _normalize_finite_float32_embedding, _require_finite_embedding_source

VECTORS_DIRNAME = "vectors"


@dataclass(frozen=True)
class Embedding:
    """One embedding to export, with the columns it is read from.

    Parameters
    ----------
    values
        A real ``(n_cells, width)`` array, or a plain ``(n_cells,)`` vector for
        a one-column array. Sparse input is densified.
    columns
        The 0-based columns to plot, one per axis, so their count is the
        embedding's dimension (1, 2, or 3). Required when ``values`` has more
        than three columns; otherwise it defaults to every column. In R the
        same argument is 1-based; the export records 0-based columns either way.
    label
        Text the viewer shows for this embedding, such as ``"t-SNE"``. Defaults
        to the embedding's id.
    source_key
        Where the array came from, such as the AnnData ``obsm`` key. It is
        recorded as provenance in ``dataset_identity.json`` and changes
        nothing else. Defaults to None.
    """

    values: Any
    columns: Sequence[int] | None = None
    label: str | None = None
    source_key: str | None = None


@dataclass(frozen=True)
class _ResolvedEmbedding:
    """One validated, normalized embedding item, in export order."""

    id: str
    index: int
    label: str
    dimension: int
    columns: list[int]
    source_key: str | None
    coordinates: np.ndarray
    scale_factor: float
    original_range: float

    def file(self, gzip_suffix: str) -> str:
        """The exact payload path the format derives for this item."""
        return f"{EMBEDDINGS_DIRNAME}/{self.index}_{self.dimension}d.bin{gzip_suffix}"


@dataclass(frozen=True)
class _ResolvedEmbeddings:
    items: list[_ResolvedEmbedding]
    defaults: dict[str, str]
    n_cells: int

    @property
    def default_dimension(self) -> int:
        return max(item.dimension for item in self.items)

    def by_id(self) -> dict[str, _ResolvedEmbedding]:
        return {item.id: item for item in self.items}


def _dense_source(values: object, *, label: str) -> np.ndarray:
    if sparse.issparse(values):
        values = cast(sparse.spmatrix, values).toarray()
    array = np.asarray(values)
    if array.dtype.kind not in {"i", "u", "f"}:
        raise TypeError(f"{label} must contain real numeric values.")
    if array.ndim == 1:
        return array.reshape(-1, 1)
    if array.ndim != 2:
        raise ValueError(f"{label} must be a 2D array or a 1D vector, got shape {array.shape}.")
    return array


def _resolve_one_embedding(embedding_id: str, value: object) -> tuple[np.ndarray, list[int], str, str | None]:
    """Validate one item and return its selected source columns, unnormalized."""
    label = f"Embedding {embedding_id!r}"
    if isinstance(value, Embedding):
        values = value.values
        requested_columns = value.columns
        display_label = value.label
        source_key = value.source_key
    else:
        values = value
        requested_columns = None
        display_label = None
        source_key = None
    if values is None:
        raise ValueError(f"{label} has no values.")
    array = _dense_source(values, label=label)
    width = int(array.shape[1])
    if requested_columns is None:
        if width not in (1, 2, 3):
            raise ValueError(
                f"{label} has {width} columns, and an embedding plots 1, 2, or 3 of them. "
                f"Pass Embedding(values, columns=[...]) naming the 0-based columns to plot, "
                "for example columns=[0, 1]."
            )
        columns = default_columns(width)
    else:
        if isinstance(requested_columns, str | bytes) or not isinstance(
            requested_columns, Sequence
        ):
            raise TypeError(f"{label} columns must be a sequence of 0-based integers.")
        requested = list(requested_columns)
        if not 1 <= len(requested) <= 3:
            raise ValueError(
                f"{label} columns must name 1, 2, or 3 columns, one per axis; got {len(requested)}."
            )
        try:
            columns = validate_embedding_columns(requested, dimension=len(requested), width=width)
        except DataMappingError as error:
            raise ValueError(f"{label}: {error}") from error
    resolved_label = (
        embedding_id
        if display_label is None
        else _require_display_text(display_label, label=f"{label} label", allow_empty=False)
    )
    resolved_source_key = (
        None
        if source_key is None
        else _require_nonempty_string(source_key, label=f"{label} source_key")
    )
    selected = _require_finite_embedding_source(array[:, columns], label=label)
    return selected, columns, resolved_label, resolved_source_key


def _resolve_embeddings(
    embeddings: object,
    embedding_defaults: object,
) -> _ResolvedEmbeddings:
    """Validate every embedding and the per-dimension defaults, then normalize.

    Everything is checked before anything is normalized, and nothing here
    writes: a defect anywhere fails the export before its first output path.
    """
    if embeddings is None:
        raise ValueError(
            "embeddings is required: map at least one embedding id to an array, "
            "for example embeddings={'X_umap': adata.obsm['X_umap']}."
        )
    if isinstance(embeddings, str | bytes) or not isinstance(embeddings, Mapping):
        raise TypeError("embeddings must be a mapping of embedding ids to arrays or Embedding values.")
    if not embeddings:
        raise ValueError("embeddings must contain at least one embedding.")
    ids = _require_field_identities(list(embeddings), label="Embedding")

    selected: dict[str, tuple[np.ndarray, list[int], str, str | None]] = {}
    n_cells: int | None = None
    first_id: str | None = None
    for embedding_id in ids:
        resolved = _resolve_one_embedding(embedding_id, embeddings[embedding_id])
        rows = int(resolved[0].shape[0])
        if n_cells is None:
            n_cells, first_id = rows, embedding_id
        elif rows != n_cells:
            raise ValueError(
                "All embeddings must have the same number of cells: "
                f"{first_id!r} has {n_cells} and {embedding_id!r} has {rows}."
            )
        selected[embedding_id] = resolved
    if n_cells is None or n_cells <= 0:
        raise ValueError("Embeddings must contain at least one cell.")

    order = sorted(ids, key=lambda embedding_id: (len(selected[embedding_id][1]), embedding_id))
    ids_by_dimension: dict[int, list[str]] = {}
    for embedding_id in order:
        ids_by_dimension.setdefault(len(selected[embedding_id][1]), []).append(embedding_id)

    defaults = _resolve_embedding_defaults(embedding_defaults, ids_by_dimension)

    items: list[_ResolvedEmbedding] = []
    for index, embedding_id in enumerate(order):
        source, columns, display_label, source_key = selected[embedding_id]
        normalized, _center, scale_factor, max_range = _normalize_finite_float32_embedding(
            source,
            label=f"Embedding {embedding_id!r}",
        )
        items.append(
            _ResolvedEmbedding(
                id=embedding_id,
                index=index,
                label=display_label,
                dimension=len(columns),
                columns=columns,
                source_key=source_key,
                coordinates=normalized,
                scale_factor=float(scale_factor),
                original_range=float(max_range),
            )
        )
    return _ResolvedEmbeddings(items=items, defaults=defaults, n_cells=n_cells)


def _resolve_embedding_defaults(
    embedding_defaults: object,
    ids_by_dimension: Mapping[int, list[str]],
) -> dict[str, str]:
    """Resolve the one default item of every dimension the export carries."""
    if embedding_defaults is None:
        requested: dict[str, str] = {}
    else:
        if isinstance(embedding_defaults, str | bytes) or not isinstance(
            embedding_defaults, Mapping
        ):
            raise TypeError(
                "embedding_defaults must map '1d', '2d', or '3d' to an embedding id, "
                "for example {'2d': 'X_umap'}."
            )
        requested = {}
        for slot, embedding_id in embedding_defaults.items():
            if slot not in DIMENSION_SLOTS:
                raise ValueError(
                    f"embedding_defaults key {slot!r} must be '1d', '2d', or '3d'."
                )
            dimension = int(slot[0])
            if dimension not in ids_by_dimension:
                raise ValueError(
                    f"embedding_defaults names a {slot} default, and no {dimension}D "
                    "embedding was supplied."
                )
            if not isinstance(embedding_id, str) or embedding_id not in ids_by_dimension[dimension]:
                raise ValueError(
                    f"embedding_defaults[{slot!r}] is {embedding_id!r}, which is not a "
                    f"{dimension}D embedding. The {dimension}D embeddings are: "
                    f"{ids_by_dimension[dimension]}."
                )
            requested[slot] = embedding_id

    defaults: dict[str, str] = {}
    for dimension in sorted(ids_by_dimension):
        slot = dimension_slot(dimension)
        candidates = ids_by_dimension[dimension]
        if slot in requested:
            defaults[slot] = requested[slot]
        elif len(candidates) == 1:
            defaults[slot] = candidates[0]
        else:
            raise ValueError(
                f"{len(candidates)} {dimension}D embeddings were supplied ({candidates}), so "
                f"embedding_defaults must name the one shown first, for example "
                f"embedding_defaults={{{slot!r}: {candidates[0]!r}}}."
            )
    return defaults


@dataclass(frozen=True)
class _ResolvedVectorFields:
    """Validated, scaled vector fields, keyed field id -> item id, in export order."""

    fields: dict[str, dict[str, np.ndarray]]
    default_field: str | None

    def index_by_field(self) -> dict[str, int]:
        return {field_id: index for index, field_id in enumerate(self.fields)}


def _resolve_vector_fields(
    vector_fields: object,
    *,
    embeddings: _ResolvedEmbeddings,
    vector_field_default: object,
) -> _ResolvedVectorFields:
    """Validate every vector field against the item it is bound to, then scale it."""
    by_id = embeddings.by_id()
    validated: ValidatedVectorFields = validate_vector_fields(
        vector_fields,
        n_cells=embeddings.n_cells,
        embedding_dimensions={item.id: item.dimension for item in embeddings.items},
        vector_field_default=vector_field_default,
    )
    if validated.fields:
        _require_field_identities(list(validated.fields), label="Vector field")
    scaled: dict[str, dict[str, np.ndarray]] = {}
    for field_id, by_embedding in validated.fields.items():
        scaled[field_id] = {
            embedding_id: scale_vector_field(
                vectors,
                scale_factor=by_id[embedding_id].scale_factor,
                label=f"Vector field {field_id!r} on embedding {embedding_id!r}",
            )
            for embedding_id, vectors in by_embedding.items()
        }
    return _ResolvedVectorFields(fields=scaled, default_field=validated.default_field)
