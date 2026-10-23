"""
The launch defaults a dataset opens with, checked against what its store holds.

The same choices reach the viewer two ways: a server passes them in its
``/_cellucid/datasets`` entry (``serve_anndata``, ``show_anndata``, ``cellucid
serve``), and :func:`cellucid.prepare` writes them into the file itself as
``uns['cellucid']['defaults']``. Both are checked here, against the names and
widths of what the store holds (metadata only), and produce the same document,
which the web app checks again with ``requireAnnDataLaunchDefaults``. Keys are
the web app's: ``embedding1d``/``embedding2d``/``embedding3d``, ``knn``,
``expression``, ``geneIdColumn``, ``latentKey``, ``vectorFieldDefault``,
``obsKeys``, ``normalizeEmbeddings`` and ``centroidMinPoints``. A choice left at None is not written.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from ._anndata_files import AnnDataContents
from ._contracts import _require_native_boolean, _require_positive_native_integer
from .data_mapping import (
    DIMENSION_SLOTS,
    DataMappingError,
    EmbeddingKeySpec,
    ObsmEntry,
    default_columns,
    enabled_embedding_keys,
    slot_dimension,
    validate_embedding_columns,
    vector_field_keys,
)


def _require_expression_source(value: object) -> str:
    """Validate ``expression_source``: ``"X"``, ``"raw"`` or ``"layer:<name>"``."""
    if not isinstance(value, str) or not (
        value in ("X", "raw") or (value.startswith("layer:") and len(value) > len("layer:"))
    ):
        raise ValueError(
            f"expression_source must be 'X', 'raw', or 'layer:<name>', got {value!r}."
        )
    return value


def _require_embedding_keys(
    value: object, *, label: str
) -> dict[str, tuple[str, list[int] | None]]:
    """Validate ``{"2d": "X_tsne"}`` or ``{"2d": ("X_pca", [0, 1])}``, named ``label``."""
    if value is None:
        return {}
    if isinstance(value, str | bytes) or not isinstance(value, Mapping):
        raise TypeError(
            f"{label} must map '1d', '2d', or '3d' to an obsm key or to a "
            "(key, columns) pair, for example {'2d': 'X_tsne'}."
        )
    resolved: dict[str, tuple[str, list[int] | None]] = {}
    for slot, spec in value.items():
        if slot not in DIMENSION_SLOTS:
            raise ValueError(f"{label} key {slot!r} must be '1d', '2d', or '3d'.")
        if isinstance(spec, str):
            key, columns = spec, None
        elif isinstance(spec, tuple) and len(spec) == 2 and isinstance(spec[0], str):
            key = spec[0]
            if isinstance(spec[1], str | bytes) or not isinstance(spec[1], Sequence):
                raise TypeError(f"{label}[{slot!r}] columns must be a sequence of integers.")
            columns = list(spec[1])
        else:
            raise TypeError(
                f"{label}[{slot!r}] must be an obsm key or a (key, columns) pair."
            )
        if not key:
            raise ValueError(f"{label}[{slot!r}] must name a non-empty obsm key.")
        resolved[slot] = (key, columns)
    return resolved


def _optional_key(value: object, *, label: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not value or not value.strip():
        raise ValueError(f"{label} must be None or one non-blank string.")
    return value


def _require_obsm_key(contents: AnnDataContents, key: str, *, label: str) -> int:
    """The width of ``obsm[key]``, which must be a 2-D array."""
    if key not in contents.obsm:
        raise ValueError(
            f"{label} names obsm[{key!r}], which this store does not have. "
            f"obsm holds: {sorted(contents.obsm)}."
        )
    width = contents.obsm[key]
    if width is None:
        raise ValueError(f"{label} names obsm[{key!r}], which is not a 2-D array.")
    return width


def launch_defaults(
    contents: AnnDataContents,
    *,
    latent_key: str | None,
    gene_id_column: str | None,
    normalize_embeddings: bool | None,
    centroid_min_points: int | None,
    obs_keys: Sequence[str] | None,
    embedding_keys: Mapping[str, EmbeddingKeySpec] | None,
    knn_key: str | None,
    expression_source: str | None,
    vector_field_default: str | None,
    embedding_label: str = "embedding_keys",
    gene_id_label: str = "gene_id_column",
) -> dict[str, Any]:
    """Check every launch choice against the store and return the ``defaults`` document.

    Only metadata is consulted: each named element must exist and have a shape
    the choice can use. A choice left at None is not written: the viewer then
    takes it from the store's own ``uns['cellucid']['defaults']``, else applies
    its default rule. ``embedding_label`` and ``gene_id_label`` are the names
    the caller's arguments go by, for messages.
    """
    defaults: dict[str, Any] = {}

    for slot, (key, columns) in _require_embedding_keys(embedding_keys, label=embedding_label).items():
        dimension = slot_dimension(slot)
        width = _require_obsm_key(contents, key, label=f"{embedding_label}[{slot!r}]")
        try:
            resolved = validate_embedding_columns(
                default_columns(dimension) if columns is None else columns,
                dimension=dimension,
                width=width,
            )
        except DataMappingError as error:
            raise ValueError(f"{embedding_label}[{slot!r}] (obsm[{key!r}]): {error}") from error
        defaults[f"embedding{slot}"] = {"key": key, "columns": resolved}

    knn = _optional_key(knn_key, label="knn_key")
    if knn is not None:
        if knn not in contents.obsp:
            raise ValueError(
                f"knn_key names obsp[{knn!r}], which this store does not have. "
                f"obsp holds: {sorted(contents.obsp)}. Compute the neighbour graph "
                "with sc.pp.neighbors(adata), or leave knn_key unset."
            )
        defaults["knn"] = knn

    gene_axis_columns = contents.var_columns
    if expression_source is not None:
        expression = _require_expression_source(expression_source)
        if expression == "X" and not contents.has_x:
            raise ValueError("expression_source 'X' names X, which this store does not have.")
        if expression == "raw":
            if contents.raw_var_columns is None:
                raise ValueError("expression_source 'raw' names raw.X, which this store does not have.")
            gene_axis_columns = contents.raw_var_columns
        if expression.startswith("layer:") and expression[len("layer:") :] not in contents.layers:
            raise ValueError(
                f"expression_source {expression!r} names a layer this store does not have. "
                f"layers holds: {sorted(contents.layers)}."
            )
        # The viewer names matrices by their place in the store.
        defaults["expression"] = (
            "raw/X"
            if expression == "raw"
            else f"layers/{expression[len('layer:') :]}"
            if expression.startswith("layer:")
            else "X"
        )

    column = _optional_key(gene_id_column, label=gene_id_label)
    if column is not None:
        if column not in gene_axis_columns:
            axis = "raw.var" if expression_source == "raw" else "var"
            raise ValueError(
                f"{gene_id_label} names {axis}[{column!r}], which this store does not have. "
                f"{axis} holds: {sorted(gene_axis_columns)}."
            )
        defaults["geneIdColumn"] = column

    latent = _optional_key(latent_key, label="latent_key")
    if latent is not None:
        _require_obsm_key(contents, latent, label="latent_key")
        defaults["latentKey"] = latent

    field = _optional_key(vector_field_default, label="vector_field_default")
    if field is not None:
        # A field is drawn from a key `<field>_<basis>` as wide as the embedding
        # `X_<basis>` it is named for, when the viewer shows that embedding as
        # the store opens (the web app's `defaultVectorBindings`): only such a
        # field can be the one selected first.
        shown = enabled_embedding_keys(
            [
                ObsmEntry(key=key, width=width, usable=width is not None and width >= 1)
                for key, width in contents.obsm.items()
            ],
            named_slots={
                slot[len("embedding") :]: ref["key"]
                for slot, ref in defaults.items()
                if slot.startswith("embedding")
            },
            latent_key=defaults.get("latentKey"),
        )
        declared: dict[str, bool] = {}
        for key, parsed in vector_field_keys(contents.obsm).items():
            vector_width = contents.obsm[key]
            if vector_width is not None and vector_width == contents.obsm.get(parsed.embedding_key):
                drawn = parsed.embedding_key in shown
                declared[parsed.field] = declared.get(parsed.field, False) or drawn
        if field not in declared:
            raise ValueError(
                f"vector_field_default names field {field!r}, which this store does not "
                f"declare. Its vector fields (obsm '<field>_<basis>' as wide as "
                f"'X_<basis>'): {sorted(declared)}."
            )
        if not declared[field]:
            raise ValueError(
                f"vector_field_default names field {field!r}, whose embeddings "
                f"(obsm 'X_<basis>' of each '{field}_<basis>') the viewer does not show "
                "when the store opens, so the field is not drawn. Show one of them "
                "(embedding_keys), or choose a field drawn on a shown embedding: "
                f"{sorted(name for name, drawn in declared.items() if drawn)}."
            )
        defaults["vectorFieldDefault"] = field

    if obs_keys is not None:
        if isinstance(obs_keys, str | bytes) or not isinstance(obs_keys, Sequence):
            raise TypeError("obs_keys must be None or a sequence of obs column names.")
        names = list(obs_keys)
        if not names:
            raise ValueError("obs_keys must name at least one obs column, or be None.")
        for name in names:
            if not isinstance(name, str) or not name:
                raise ValueError("obs_keys must hold non-empty strings.")
            if name not in contents.obs_columns:
                raise ValueError(
                    f"obs_keys names obs[{name!r}], which this store does not have. "
                    f"obs holds: {sorted(contents.obs_columns)}."
                )
        if len(set(names)) != len(names):
            raise ValueError("obs_keys names a column more than once.")
        defaults["obsKeys"] = names

    if normalize_embeddings is not None:
        defaults["normalizeEmbeddings"] = _require_native_boolean(
            normalize_embeddings, label="normalize_embeddings"
        )
    if centroid_min_points is not None:
        defaults["centroidMinPoints"] = _require_positive_native_integer(
            centroid_min_points, label="centroid_min_points"
        )
    return defaults
