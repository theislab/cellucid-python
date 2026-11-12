"""
What ``prepare(data=...)`` takes from an AnnData or a MuData.

Extraction fills the same components a caller passes by hand -- ``embeddings``,
``obs``, ``var``, ``gene_expression``, ``connectivities``, ``vector_fields``,
the identity and the launch defaults -- and :func:`cellucid.prepare` checks and
writes them exactly as it checks and writes explicit ones, so an explicit
argument overrides one component and the output is byte for byte what the
same components passed by hand produce. cellucid-r's ``cellucid_prepare(data
= ...)`` follows the same rule.

Nothing here imports anndata or mudata: an object is recognized by its class,
so a MuData is handled without ``mudata`` installed (it could not exist
without it).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from ..data_mapping import default_graph, vector_field_keys
from ._totals import NAMED_OVER, OVER_SOURCE_PANEL

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
    """The ``prepare()`` arguments ``uns['cellucid']`` fills: identity, citation, source,
    defaults, and the column of the cells' totals over all genes."""
    record = uns.get("cellucid") if isinstance(uns, Mapping) else None
    if record is None:
        return {}
    if not isinstance(record, Mapping):
        raise TypeError(f"{what}.uns['cellucid'] must be a mapping, got {type(record).__name__}.")
    # Only the keys read below are converted: the record of a Cellucid-ready
    # dataset (one prepare() wrote) also holds its outlier quantiles (one value per cell and categorical column),
    # which prepare() recomputes and would otherwise turn into Python lists.
    record = {
        str(key): native(item)
        for key, item in record.items()
        if str(key) in ("id", "name", "description", "citation", "source", "defaults", "totalCounts")
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
    # The caller's column of the cells' totals (over all genes, or over the
    # source's gene panel) stays named (``total_counts=``, with what it
    # covers, ``total_counts_over=``); totals over the stored genes are summed
    # again from the genes then stored (``._totals``).
    totals = record.get("totalCounts")
    if isinstance(totals, Mapping):
        entries = {None: totals} if "field" in totals else totals
        named = {
            modality: entry
            for modality, entry in entries.items()
            if isinstance(entry, Mapping) and entry.get("over") in NAMED_OVER
        }
        if named:
            filled["total_counts"] = (
                named[None]["field"] if None in named else {str(modality): entry["field"] for modality, entry in named.items()}
            )
            panel = [modality for modality, entry in named.items() if entry["over"] == OVER_SOURCE_PANEL]
            if panel:
                filled["total_counts_over"] = (
                    OVER_SOURCE_PANEL if None in named else {str(modality): OVER_SOURCE_PANEL for modality in panel}
                )
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


#: The component holding ``X``'s ``var`` when the matrix chosen is not ``X``
#: (``raw``): the genes of a ``gene_expression=`` given explicitly.
X_VAR_COMPONENT = "_x_var"


def _x_holds_scaled_values(x: Any) -> bool:
    """Whether ``X`` holds scaled values as the default expression rule reads
    it (``expressionRuleInputs``): a dense float matrix whose first row holds a
    negative value -- ``sc.pp.scale`` centers every gene, so about half of any
    cell's values are negative, while counts and normalized values never are.
    A sparse matrix is not scaled. One row read."""
    if sparse.issparse(x) or hasattr(x, "format") or getattr(x, "ndim", None) != 2:
        return False
    if np.dtype(x.dtype).kind != "f" or x.shape[0] == 0 or x.shape[1] == 0:
        return False
    return bool((np.asarray(x[0:1]) < 0).any())


def default_expression_of(adata: Any, *, where: str = "adata") -> tuple[str | None, Any, Any, str | None]:
    """The matrix the viewer opens ``adata`` with (:func:`..data_mapping.default_expression`):
    its name (``"X"``, ``"raw"``, ``"layer:<name>"``, None), the matrix, its
    genes' ``var``, and a note saying why it is not ``X`` (None when it is)."""
    from ..data_mapping import default_expression

    has_x = adata.X is not None
    raw = adata.raw
    chosen = default_expression(
        has_x=has_x,
        has_raw=raw is not None,
        layers=[str(name) for name in adata.layers],
        raw_has_more_genes=has_x and raw is not None and raw.n_vars > adata.n_vars,
        x_scaled=has_x and raw is not None and _x_holds_scaled_values(adata.X),
    )
    keep_x = (
        "to write X instead: gene_expression=adata.X, var=adata.var" if where == "adata"
        else f"to write X instead, set {where}.raw = None"
    )
    if chosen.expression == "raw":
        why = {
            "raw-more-genes": f"raw holds more genes than X ({raw.n_vars:,} against {adata.n_vars:,})",
            "x-scaled": "X holds scaled values (negative ones, from sc.pp.scale)",
            None: "there is no X",
        }[chosen.reason]
        note = (
            f"{'' if where == 'adata' else where + ': '}raw.X written as X, with raw.var as its genes: {why}, and "
            f"Scanpy and the viewer read genes from raw then{'' if not has_x else f' ({keep_x})'}."
        )
        return "raw", raw.X, raw.var, note
    if chosen.expression is not None and chosen.expression.startswith("layer:"):
        name = chosen.expression[len("layer:"):]
        another = (
            "gene_expression=adata.layers[name]" if where == "adata" else f"{where}.X = {where}.layers[name]"
        )
        note = (
            f"{'' if where == 'adata' else where + ': '}layers[{name!r}] written as X: there is no X or raw, and "
            f"the viewer opens this layer (to write another: {another})."
        )
        return chosen.expression, adata.layers[name], adata.var, note
    if chosen.expression == "X":
        return "X", adata.X, adata.var, None
    return None, None, None, None


def anndata_components(adata: Any, notes: list[tuple[str, str]] | None = None) -> dict[str, Any]:
    """Every component ``prepare()`` takes, as ``adata`` holds it.

    The expression matrix and its genes are the ones the viewer opens the
    AnnData with (:func:`default_expression_of`): ``raw`` when it holds more
    genes than ``X`` or ``X`` is scaled, else ``X``; without ``X``, ``raw``
    or a layer; ``obsm`` gives the embeddings and vector fields
    (:func:`split_obsm`); ``obsp['connectivities']`` the neighbor graph;
    ``uns['cellucid']`` the identity and launch defaults. What is not written
    goes in ``notes``, by the component it would have been part of (shown
    when that component comes from ``adata``), as cellucid-r names what it
    does not write: obs and var columns of no AnnData encoding Cellucid
    writes, obsm entries that are not numeric arrays, layers, ``raw`` and
    other obsp graphs.
    """
    from ._inputs import ObsColors, held_obs_colors, storable_columns

    notes = [] if notes is None else notes
    embeddings, fields = split_obsm(adata.obsm)
    obs, note = storable_columns(adata.obs, "obs")
    if note is not None:
        notes.append(("obs", note))
    for key, value in adata.obsm.items():
        if not is_cell_array(value):
            notes.append(("embeddings", f"obsm[{key!r}] holds no numeric coordinates and was not read."))
    components: dict[str, Any] = {
        "embeddings": embeddings,
        "obs": obs,
        # Scanpy's colors of its columns, uns['<column>_colors'].
        "obs_colors": ObsColors(colors=held_obs_colors(adata.uns, list(obs.columns)), explicit=False),
    }
    if fields:
        components["vector_fields"] = fields
    # The matrix the viewer opens this AnnData with (default_expression), so
    # the file prepared from it opens on the same genes and values.
    expression, matrix, genes, chosen_note = default_expression_of(adata)
    if matrix is not None:
        components["gene_expression"] = matrix
        var, note = storable_columns(genes, "raw.var" if expression == "raw" else "var")
        components["var"] = var
        if note is not None:
            notes.append(("var", note))
    if chosen_note is not None:
        notes.append(("gene_expression", chosen_note))
    if expression not in (None, "X") and adata.X is not None:
        # X's genes, for a gene_expression= given in place of the matrix chosen
        # (a layer of X's genes, X itself): an explicit matrix goes with var.
        components[X_VAR_COMPONENT], _ = storable_columns(adata.var, "var")
    left = [name for name in adata.layers if expression != f"layer:{name}"]
    if left:
        notes.append(("gene_expression", f"layers not read: {left} (prepare() writes one matrix; to write a layer: gene_expression=adata.layers[name])."))
    if adata.raw is not None and expression != "raw":
        notes.append(("gene_expression", "raw not read (prepare() writes one matrix; to write raw: gene_expression=adata.raw.X, var=adata.raw.var)."))
    # The graph the viewer's default rule opens with (default_graph):
    # connectivities, else Seurat's <assay>_snn (an AnnData converted from a
    # Seurat object), written as obsp['connectivities'].
    graph = default_graph(adata.obsp.keys())
    other = [key for key in adata.obsp if key != graph]
    if other:
        notes.append(("connectivities", f"obsp {other} not read (prepare() writes one graph, obsp[{graph!r}], as obsp['connectivities'], the viewer's neighbor graph)." if graph is not None else f"obsp {other} not read (no undirected neighbor graph: prepare() writes obsp['connectivities'] or a Seurat '<assay>_snn' graph)."))
    if graph is not None:
        components["connectivities"] = adata.obsp[graph]
    components.update(cellucid_record(adata.uns, what="data"))
    return components
