"""
What ``prepare(data=...)`` takes from an AnnData or a MuData.

Extraction fills the same components a caller passes by hand -- ``embeddings``,
``obs``, ``var``, ``gene_expression``, ``connectivities``, ``vector_fields``,
the identity and the launch defaults -- and :func:`cellucid.prepare` checks and
writes them exactly as it checks and writes explicit ones, so an explicit
argument overrides one component and the output is byte for byte what the
same components passed by hand produce. cellucid-r's ``cellucid_prepare(data
= ...)`` follows the same rule.

Nothing here imports anndata or mudata: an object is recognised by its class,
so a MuData is handled without ``mudata`` installed (it could not exist
without it).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from ..data_mapping import vector_field_keys

#: The launch defaults ``uns['cellucid']['defaults']`` may hold, by the
#: ``prepare()`` argument each one fills. ``knn`` and ``expression`` are not
#: here: ``prepare()`` writes ``obsp['connectivities']`` and ``X`` alone, and
#: the viewer's own rule finds them.
_DEFAULT_ARGUMENTS = {
    "latentKey": "latent_key",
    "obsKeys": "obs_keys",
    "geneIdColumn": "var_gene_id_column",
    "vectorFieldDefault": "vector_field_default",
    "normalizeEmbeddings": "normalize_embeddings",
    "centroidMinPoints": "centroid_min_points",
    "modality": "modality",
}


def _class_names(value: object) -> set[tuple[str, str]]:
    return {(cls.__module__.split(".", 1)[0], cls.__name__) for cls in type(value).__mro__}




def is_anndata(value: object) -> bool:
    """Whether ``value`` is an ``anndata.AnnData``."""
    return ("anndata", "AnnData") in _class_names(value)


def native(value: Any) -> Any:
    """A value read from ``uns`` as plain Python: arrays become lists, NumPy scalars numbers."""
    if isinstance(value, Mapping):
        return {str(key): native(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return [native(item) for item in value.tolist()]
    if isinstance(value, list | tuple):
        return [native(item) for item in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.str_):
        return str(value)
    return value


def cellucid_record(uns: Mapping[str, Any], *, what: str) -> dict[str, Any]:
    """The ``prepare()`` arguments ``uns['cellucid']`` fills: identity, citation, source, defaults."""
    record = uns.get("cellucid") if isinstance(uns, Mapping) else None
    if record is None:
        return {}
    if not isinstance(record, Mapping):
        raise TypeError(f"{what}.uns['cellucid'] must be a mapping, got {type(record).__name__}.")
    # Only the keys read below are converted: a prepared dataset's record also
    # holds its outlier quantiles (one value per cell and categorical column),
    # which prepare() recomputes and would otherwise turn into Python lists.
    record = {
        str(key): native(item)
        for key, item in record.items()
        if str(key) in ("id", "name", "description", "citation", "source", "defaults")
    }
    filled: dict[str, Any] = {}
    for key, argument in (("id", "dataset_id"), ("name", "dataset_name"), ("description", "dataset_description")):
        if key in record:
            filled[argument] = record[key]
    for key, text_key, (text_argument, url_argument) in (
        ("citation", "text", ("citation_text", "citation_url")),
        ("source", "name", ("source_name", "source_url")),
    ):
        linked = record.get(key)
        if isinstance(linked, Mapping):
            if text_key in linked:
                filled[text_argument] = linked[text_key]
            if "url" in linked:
                filled[url_argument] = linked["url"]
    defaults = record.get("defaults")
    if isinstance(defaults, Mapping):
        slots: dict[str, Any] = {}
        for slot in ("1d", "2d", "3d"):
            spec = defaults.get(f"embedding{slot}")
            if isinstance(spec, str):
                slots[slot] = spec
            elif isinstance(spec, Mapping) and isinstance(spec.get("key"), str):
                columns = spec.get("columns")
                slots[slot] = spec["key"] if columns is None else (spec["key"], list(columns))
        if slots:
            filled["embedding_defaults"] = slots
        for key, argument in _DEFAULT_ARGUMENTS.items():
            if defaults.get(key) is not None:
                filled[argument] = defaults[key]
    return filled


def is_cell_array(value: object) -> bool:
    """A per-cell numeric array an ``obsm`` key can hold: 1-D or 2-D numbers (not a table)."""
    if isinstance(value, pd.DataFrame):
        return False
    if sparse.issparse(value):
        matrix: Any = value
        return bool(matrix.ndim == 2 and np.dtype(matrix.dtype).kind in "iuf")
    dtype = getattr(value, "dtype", None)
    shape = getattr(value, "shape", None)
    return dtype is not None and shape is not None and np.dtype(dtype).kind in "iuf" and len(shape) in (1, 2)


def _width(value: Any) -> int:
    shape = value.shape
    return 1 if len(shape) == 1 else int(shape[1])


def split_obsm(obsm: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    """An object's ``obsm`` as ``prepare()`` takes it: embeddings, and vector fields.

    Every numeric array is kept, under its own key. One named
    ``<field>_<basis>`` (scVelo's name) as wide as the embedding ``X_<basis>``
    is the vector field ``<field>`` on that embedding; ``prepare()`` writes it
    back under the same key. Tables and non-numeric arrays are not embeddings
    and are left out.
    """
    arrays = {key: value for key, value in obsm.items() if is_cell_array(value)}
    fields: dict[str, dict[str, Any]] = {}
    vector_keys: set[str] = set()
    for key, parsed in vector_field_keys(arrays).items():
        embedding = arrays.get(parsed.embedding_key)
        if embedding is not None and _width(embedding) == _width(arrays[key]):
            fields.setdefault(parsed.field, {})[parsed.embedding_key] = arrays[key]
            vector_keys.add(key)
    embeddings = {key: value for key, value in arrays.items() if key not in vector_keys}
    return embeddings, fields


def anndata_components(adata: Any) -> dict[str, Any]:
    """Every component ``prepare()`` takes, as ``adata`` holds it.

    ``X`` is the expression matrix and ``var`` its genes (neither when ``X``
    is absent); ``obsm`` gives the embeddings and vector fields
    (:func:`split_obsm`); ``obsp['connectivities']`` the neighbour graph;
    ``uns['cellucid']`` the identity and launch defaults.
    """
    embeddings, fields = split_obsm(adata.obsm)
    components: dict[str, Any] = {"embeddings": embeddings, "obs": adata.obs}
    if fields:
        components["vector_fields"] = fields
    if adata.X is not None:
        components["gene_expression"] = adata.X
        components["var"] = adata.var
    if "connectivities" in adata.obsp:
        components["connectivities"] = adata.obsp["connectivities"]
    components.update(cellucid_record(adata.uns, what="data"))
    return components
