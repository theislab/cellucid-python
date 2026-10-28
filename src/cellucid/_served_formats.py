"""
Loom and 10x Cell Ranger data as the web app presents them, for the store a
server makes in priority order (:mod:`cellucid._served_store`).

The browser reads a Loom file or a Cell Ranger matrix through its own layout
views, which present it as AnnData (``views/loom.js``, ``views/tenx-h5.js``,
``views/tenx-outs.js`` in the web app). When Python serves such a file by
path, it presents the same AnnData here, by the same rules, so the store in
priority order holds what the browser would show from the file itself --
the same ``obs``, ``var``, embeddings, graphs and matrices, cell for cell --
with its cells in priority order (``tests/served-priority-store.test.mjs``
in the web app opens both and compares them):

======================================  ======================================
Loom                                    presented as
======================================  ======================================
``/matrix`` (genes x cells)             ``X``, cells x genes, read transposed
``/layers/<k>``                         ``layers/<k>`` (``''`` is ``total``)
``/col_attrs/<k>``, one value per cell  ``obs`` columns
``/col_attrs/<k>``, n x k               ``obsm/<k>``
``<b>_X``/``<b>_Y``[/``<b>_Z``], ...    ``obsm/X_<b>`` (``X_embedding``)
compound ``Embedding``, ``Embeddings_*`` ``obsm/X_embedding[_<id>]``
``CellID`` (``obs_names``, ...)         the obs index
``Gene`` (``var_names``, ...)           the var index, the rest ``var``
``/col_graphs/<g>``                     ``obsp/<g>``, CSR
global attributes                       ``uns``
======================================  ======================================

A Cell Ranger matrix (stored by barcode, AnnData's CSR) is ``X`` with the
names ``scanpy.read_10x_h5`` gives; several feature types (or genomes) are
modalities, each a column window of the one matrix. An ``outs/`` folder adds
Cell Ranger's projections (``obsm``) and clusterings (``obs``), matched by
barcode.

Arrays with one row per cell are read where they are (h5py), lazily, by the
served store; metadata with one value per cell or gene (names, attributes,
projections) is read whole, as anndata reads an ``.h5ad``'s ``obs``.
"""

from __future__ import annotations

import functools
import re
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from ._plural import plural

#: The obs and var index attributes a Loom file names its axes by, in order.
LOOM_OBS_INDEX_KEYS = ("CellID", "obs_names", "cell_names")
LOOM_VAR_INDEX_KEYS = ("Gene", "var_names", "gene_names", "Accession")

#: An embedding stored as one column per axis: the base name, and the axis suffixes.
_STACK_PATTERNS = (
    (re.compile(r"^(.*)_([XYZ])$", re.DOTALL), ("X", "Y", "Z")),
    (re.compile(r"^(.*)_([123])$", re.DOTALL), ("1", "2", "3")),
    (re.compile(r"^(.*[A-Za-z])([123])$", re.DOTALL), ("1", "2", "3")),
)

#: The matrix files of an ``outs/`` folder, preferred first.
OUTS_MATRIX_FILES = (
    "filtered_feature_bc_matrix.h5",
    "sample_filtered_feature_bc_matrix.h5",
    "filtered_gene_bc_matrices_h5.h5",
)

_PROJECTION = re.compile(r"(?:^|/)analysis/(umap|tsne|pca)/([^/]+)/projection\.csv$")
_ARC_PROJECTION = re.compile(r"(?:^|/)analysis/dimensionality_reduction/([^/]+)/(umap|tsne)_projection\.csv$")
_CLUSTERS = re.compile(r"(?:^|/)analysis/clustering/([^/]+)/clusters\.csv$")


def no_embedding_message(what: str, where: str = "") -> str:
    """The web app's refusal of data with no embedding (``views/refusal.js``)."""
    return (
        f"{what} has no cell embedding, and Cellucid draws cells at embedding coordinates. "
        + (f"{where} " if where else "")
        + ("Or compute" if where else "Compute")
        + " one and open the result: in Python, sc.pp.normalize_total(adata); sc.pp.log1p(adata); "
        "sc.pp.pca(adata); sc.pp.neighbors(adata); sc.tl.umap(adata), then adata.write_h5ad(path); "
        "in R, Seurat::NormalizeData(), FindVariableFeatures(), ScaleData(), RunPCA() and RunUMAP(), "
        "then cellucid::cellucid_prepare(path, data = obj)."
    )


@dataclass
class Presented:
    """One AnnData as the browser presents a source: metadata read, arrays where they are."""

    obs: pd.DataFrame
    var: pd.DataFrame
    x: Any
    uns: dict[str, Any] = field(default_factory=dict)
    layers: dict[str, Any] = field(default_factory=dict)
    obsm: dict[str, Any] = field(default_factory=dict)
    obsp: dict[str, Any] = field(default_factory=dict)


@dataclass
class PresentedModalities:
    """Several modalities over the same cells (a Cell Ranger matrix of several
    feature types or genomes): the shared ``obs``/``obsm``, each modality's
    ``var`` and ``X``."""

    obs: pd.DataFrame
    obsm: dict[str, Any]
    modalities: dict[str, Presented]


# =============================================================================
# Text and arrays of an HDF5 file
# =============================================================================

_REFERENCE = re.compile(r"&#x([0-9a-fA-F]+);|&#([0-9]+);|&(lt|gt|quot|apos|amp);")
_ENTITIES = {"lt": "<", "gt": ">", "quot": '"', "apos": "'", "amp": "&"}


def _unescape(text: str) -> str:
    """The XML character references loompy 2 writes into ASCII text, decoded."""
    if "&" not in text:
        return text

    def one(match: re.Match[str]) -> str:
        if match.group(1) is not None:
            return chr(int(match.group(1), 16))
        if match.group(2) is not None:
            return chr(int(match.group(2)))
        return _ENTITIES[match.group(3)]

    return _REFERENCE.sub(one, text)


def _kind(dataset: Any) -> str:
    """The web app's kind of an HDF5 dataset: ``string``, ``int``, ``uint``, ``float``, ``bool``, ``compound``, or ``other``."""
    dtype = dataset.dtype
    if dtype.names is not None:
        return "compound"
    import h5py

    if dtype.kind in "SU" or h5py.check_string_dtype(dtype) is not None:
        return "string"
    return {"i": "int", "u": "uint", "f": "float", "b": "bool"}.get(dtype.kind, "other")


def _text(dataset: Any) -> np.ndarray:
    """A string dataset as text: fixed-length ASCII unescaped (Loom 2), variable-length as UTF-8."""
    import h5py

    info = h5py.check_string_dtype(dataset.dtype)
    if dataset.dtype.kind == "S" or (info is not None and info.length is not None):
        raw = dataset[()]
        return np.asarray([_unescape(value.decode("ascii", "replace") if isinstance(value, bytes) else str(value)) for value in np.ravel(raw)], dtype=object).reshape(np.shape(raw))
    return np.asarray(dataset.asstr()[()], dtype=object)


def _values(dataset: Any) -> np.ndarray:
    return _text(dataset) if _kind(dataset) == "string" else np.asarray(dataset[()])


def _lengths(lengths: set[int]) -> str:
    """Attribute lengths as the browser names them: ``"5 or 6 rows"``."""
    return f"{' or '.join(str(n) for n in sorted(lengths))} {plural(max(lengths, default=0), 'row')}"


def _arrays_of(group: Any) -> dict[str, Any]:
    import h5py

    if not isinstance(group, h5py.Group):
        return {}
    return {key: group[key] for key in group if isinstance(group[key], h5py.Dataset)}


def _frame(index_key: str | None, index: np.ndarray | None, length: int, columns: dict[str, np.ndarray]) -> pd.DataFrame:
    """A dataframe of an index (positions as text when the file names none) and columns."""
    names = pd.Index([str(i) for i in range(length)] if index is None else [str(v) for v in index], dtype=object)
    frame = pd.DataFrame({key: values for key, values in columns.items()}, index=names)
    frame.index.name = None if index_key in (None, "_index") else index_key
    return frame


# =============================================================================
# Loom
# =============================================================================


class Transposed:
    """A stored 2-D array read transposed (a Loom matrix, genes x cells, as cells x genes)."""

    def __init__(self, array: Any) -> None:
        self.array = array
        self.shape = (int(array.shape[1]), int(array.shape[0]))
        self.dtype = np.dtype(array.dtype)

    def __getitem__(self, key: Any) -> np.ndarray:
        rows, columns = key if isinstance(key, tuple) else (key, slice(None))
        if isinstance(rows, slice):
            return np.asarray(self.array[columns, rows]).T
        rows = np.asarray(rows, dtype=np.int64)
        if len(rows) == 0:
            return np.zeros((0, len(range(*columns.indices(self.shape[1])))), dtype=self.dtype)
        # The span of the rows asked for, read as one block (a long list of
        # positions is read by HDF5 one point at a time), then taken in order.
        low, high = int(rows.min()), int(rows.max()) + 1
        block = np.asarray(self.array[columns, low:high]).T
        return np.asarray(block[rows - low])


def _coordinate_stacks(columns: dict[str, Any]) -> tuple[list[tuple[str, list[str]]], set[str]]:
    """Embeddings stored one column per axis (``views/loom.js`` ``coordinateStacks``)."""
    stacks: list[tuple[str, list[str]]] = []
    used: set[str] = set()
    for pattern, axes in _STACK_PATTERNS:
        bases: dict[str, dict[str, str]] = {}
        for key, array in columns.items():
            if key in used or _kind(array) not in ("int", "uint", "float"):
                continue
            match = pattern.match(key)
            if match is None:
                continue
            bases.setdefault(match.group(1), {})[match.group(2)] = key
        for base, by_axis in bases.items():
            if axes[0] not in by_axis or axes[1] not in by_axis:
                continue
            keys = [by_axis[axis] for axis in axes[: 3 if axes[2] in by_axis else 2]]
            used.update(keys)
            name = base.strip("_")
            stacks.append(("X_embedding" if name in ("", "X") else (name if name.startswith("X_") else f"X_{name}"), keys))
    return stacks, used


def _stacked(members: list[np.ndarray]) -> np.ndarray:
    dtype = np.float32 if all(np.asarray(m).dtype == np.float32 for m in members) else np.float64
    return np.column_stack([np.asarray(m, dtype=dtype) for m in members])


def _loom_uns(root: Any) -> dict[str, Any]:
    """The global attributes: Loom 3's ``/attrs`` arrays, else the root's own attributes (Loom 2)."""
    import h5py

    uns: dict[str, Any] = {}
    attrs = root.get("attrs")
    if isinstance(attrs, h5py.Group):
        for key, dataset in _arrays_of(attrs).items():
            if len(dataset.shape) <= 1 and _kind(dataset) != "compound":
                uns[key] = _values(dataset)
        return uns
    for key, value in root.attrs.items():
        values = list(np.ravel(value)) if np.ndim(value) else [value]
        if not values:
            continue
        texts = [v.decode("utf-8", "replace") if isinstance(v, bytes) else v for v in values]
        if all(isinstance(v, str) for v in texts):
            decoded = [_unescape(v) for v in texts]
            uns[key] = np.asarray(decoded, dtype=object) if np.ndim(value) else decoded[0]
        elif all(isinstance(v, int | float | np.integer | np.floating) and not isinstance(v, bool | np.bool_) for v in values):
            numbers = np.asarray(values, dtype=np.float64)
            uns[key] = numbers if np.ndim(value) else numbers[0]
    return uns


def _loom_graph(group: Any, n_cells: int, label: str) -> sparse.csr_matrix | None:
    a, b, w = (group.get(name) for name in ("a", "b", "w"))
    if a is None or b is None or w is None or len(a.shape) != 1 or b.shape[0] != a.shape[0] or w.shape[0] != a.shape[0]:
        return None
    rows, columns, weights = np.asarray(a[()], dtype=np.int64), np.asarray(b[()], dtype=np.int64), np.asarray(w[()], dtype=np.float64)
    for values in (rows, columns):
        outside = (values < 0) | (values >= n_cells)
        if outside.any():
            raise ValueError(
                f"Loom graph {label} names cell {int(values[outside][0])}, outside its {n_cells:,} {plural(n_cells, 'cell')}."
            )
    # Each row's entries by column, kept one by one (never summed), as the browser builds them.
    order = np.lexsort((columns, rows))
    indptr = np.concatenate([[0], np.cumsum(np.bincount(rows, minlength=n_cells))]).astype(np.int64 if len(rows) > 2**31 - 1 else np.int32)
    return sparse.csr_matrix((weights[order], columns[order].astype(np.int32), indptr), shape=(n_cells, n_cells))


def loom_presented(root: Any, name: str) -> Presented:
    """A Loom file as the browser presents it (``views/loom.js``)."""
    matrix = root["matrix"]
    col_attrs = _arrays_of(root.get("col_attrs"))
    row_attrs = _arrays_of(root.get("row_attrs"))
    first, second = (int(extent) for extent in matrix.shape)
    col_lengths = {int(array.shape[0]) for array in col_attrs.values()}
    row_lengths = {int(array.shape[0]) for array in row_attrs.values()}

    def fits(lengths: set[int], extent: int) -> bool:
        return not lengths or lengths == {extent}

    transposed = True
    if not (fits(col_lengths, second) and fits(row_lengths, first)):
        if first != second and fits(col_lengths, first) and fits(row_lengths, second):
            transposed = False
        else:
            raise ValueError(
                f"{name} is a Loom file whose cell and gene attributes do not match its {first} × {second} matrix "
                f"(cell attributes have {_lengths(col_lengths)}, gene attributes {_lengths(row_lengths)}). "
                "Write it again with loompy (loompy.create) "
                "or open the data it was made from."
            )
    n_cells, n_genes = (second, first) if transposed else (first, second)

    def as_matrix(array: Any) -> Any:
        if tuple(array.shape) != tuple(matrix.shape):
            return None
        return Transposed(array) if transposed else array

    obs_index_key = next(
        (key for key in LOOM_OBS_INDEX_KEYS if key in col_attrs and len(col_attrs[key].shape) == 1 and _kind(col_attrs[key]) == "string"),
        None,
    )
    one_dimensional = {key: array for key, array in col_attrs.items() if len(array.shape) == 1 and _kind(array) != "compound"}
    stacks, used = _coordinate_stacks({key: array for key, array in one_dimensional.items() if key != obs_index_key})
    obs_columns = {
        key: _values(array)
        for key, array in one_dimensional.items()
        if key != obs_index_key and key not in used and _kind(array) in ("string", "int", "uint", "float", "bool")
    }
    obs = _frame(obs_index_key, None if obs_index_key is None else _text(col_attrs[obs_index_key]), n_cells, obs_columns)

    obsm: dict[str, Any] = {}
    for key, array in col_attrs.items():
        if len(array.shape) == 2 and _kind(array) in ("int", "uint", "float") and array.shape[1] > 0:
            obsm[key] = np.asarray(array[()])
    for key, keys in stacks:
        if key not in obsm:
            obsm[key] = _stacked([col_attrs[member][()] for member in keys])

    def compound(array: Any, x_name: str, y_name: str, key: str) -> None:
        names = array.dtype.names or ()
        if x_name in names and y_name in names and key not in obsm and array.dtype[x_name].shape == () and array.dtype[y_name].shape == ():
            values = array[()]
            obsm[key] = _stacked([values[x_name], values[y_name]])

    embedding = col_attrs.get("Embedding")
    if embedding is not None and _kind(embedding) == "compound" and len(embedding.shape) == 1:
        compound(embedding, "_X", "_Y", "X_embedding")
    embeddings_x, embeddings_y = col_attrs.get("Embeddings_X"), col_attrs.get("Embeddings_Y")
    if embeddings_x is not None and embeddings_y is not None and _kind(embeddings_x) == "compound" and _kind(embeddings_y) == "compound":
        xs, ys = embeddings_x[()], embeddings_y[()]
        for member in embeddings_x.dtype.names or ():
            key = f"X_embedding_{member}"
            if member not in (embeddings_y.dtype.names or ()) or key in obsm:
                continue
            if embeddings_x.dtype[member].shape != () or embeddings_y.dtype[member].shape != ():
                continue
            obsm[key] = _stacked([xs[member], ys[member]])
    if not obsm:
        raise ValueError(
            no_embedding_message(
                f"{name} is a Loom file that",
                "A Loom file keeps an embedding as cell attributes (an n × 2 X_umap, or _X and _Y); this one has none.",
            )
        )

    var_index_key = next(
        (key for key in LOOM_VAR_INDEX_KEYS if key in row_attrs and len(row_attrs[key].shape) == 1 and _kind(row_attrs[key]) == "string"),
        None,
    )
    var_columns = {
        key: _values(array)
        for key, array in row_attrs.items()
        if key != var_index_key and len(array.shape) == 1 and _kind(array) in ("string", "int", "uint", "float", "bool")
    }
    var = _frame(var_index_key, None if var_index_key is None else _text(row_attrs[var_index_key]), n_genes, var_columns)

    layers: dict[str, Any] = {}
    for key, array in _arrays_of(root.get("layers")).items():
        presented = as_matrix(array)
        if presented is not None:
            layers["total" if key == "" else key] = presented

    obsp: dict[str, Any] = {}
    graphs = root.get("col_graphs")
    if graphs is not None:
        for key in graphs:
            if "/" in key:
                continue
            graph = _loom_graph(graphs[key], n_cells, f"/col_graphs/{key}")
            if graph is not None:
                obsp[key] = graph

    return Presented(obs=obs, var=var, x=as_matrix(matrix), uns=_loom_uns(root), layers=layers, obsm=obsm, obsp=obsp)


# =============================================================================
# 10x Cell Ranger
# =============================================================================


@dataclass(frozen=True)
class CellRangerMatrix:
    """One stored Cell Ranger matrix (CSC over features x barcodes, AnnData's CSR over cells x features)."""

    group: Any
    n_cells: int
    n_features: int
    index: np.ndarray
    columns: dict[str, np.ndarray]
    #: Contiguous runs of features of one type (or the genome), in stored order: ``(name, start, count)``.
    windows: list[tuple[str, int, int]]


def _text_attribute(node: Any, key: str) -> str | None:
    value = node.attrs.get(key)
    if value is None:
        return None
    first = np.ravel(value)[0] if np.ndim(value) else value
    return first.decode() if isinstance(first, bytes) else (first if isinstance(first, str) else None)


def _is_matrix_group(group: Any, names_key: str | None) -> bool:
    import h5py

    return (
        isinstance(group, h5py.Group)
        and all(isinstance(group.get(key), h5py.Dataset) for key in ("barcodes", "data", "indices", "indptr"))
        and (names_key is None or group.get(names_key) is not None)
    )


def _feature_windows(types: np.ndarray, name: str) -> list[tuple[str, int, int]]:
    windows: list[list[Any]] = []
    for feature, value in enumerate(types):
        kind = str(value)
        if windows and windows[-1][0] == kind:
            windows[-1][2] += 1
            continue
        if any(window[0] == kind for window in windows):
            raise ValueError(
                f"{name} keeps the features of type '{kind}' in more than one place, so they cannot be read as one "
                "modality. Read it with scanpy (sc.read_10x_h5(path, gex_only=False)) and write it as .h5ad, then "
                "serve that file."
            )
        windows.append([kind, feature, 1])
    return [(kind, start, count) for kind, start, count in windows]


def cell_ranger_matrices(root: Any, name: str) -> list[CellRangerMatrix]:
    """The matrices of a Cell Ranger ``.h5`` (``views/tenx-h5.js``): one (Cell Ranger 3), or one per genome (2)."""
    if _text_attribute(root, "filetype") == "molecule":
        raise ValueError(
            f"{name} is a Cell Ranger molecule table (molecule_info.h5), not a count matrix. "
            "Serve filtered_feature_bc_matrix.h5 from the same outs/ folder, or the whole outs/ folder."
        )
    matrix = root.get("matrix")
    if _is_matrix_group(matrix, "features"):
        features = matrix["features"]
        n_features = int((features.get("id") or features["name"]).shape[0])
        rename = {"id": "gene_ids", "feature_type": "feature_types"}
        keys = ["id", "feature_type", *[key for key in features if key not in ("id", "feature_type")]]
        columns: dict[str, np.ndarray] = {}
        import h5py

        for key in keys:
            array = features.get(key)
            if key in ("name", "_all_tag_keys") or not isinstance(array, h5py.Dataset):
                continue
            if len(array.shape) != 1 or array.shape[0] != n_features or _kind(array) not in ("string", "int", "uint", "float", "bool"):
                continue
            columns[rename.get(key, key)] = _values(array)
        types = features.get("feature_type")
        windows = (
            _feature_windows(_text(types), name)
            if types is not None and types.shape[0] == n_features
            else [("Gene Expression", 0, n_features)]
        )
        return [
            CellRangerMatrix(
                group=matrix,
                n_cells=int(matrix["barcodes"].shape[0]),
                n_features=n_features,
                index=_text(features["name"]),
                columns=columns,
                windows=windows,
            )
        ]
    genomes = [key for key in root if _is_matrix_group(root[key], "genes")]
    found = []
    for genome in genomes:
        group = root[genome]
        names = group["gene_names"] if group.get("gene_names") is not None else group["genes"]
        found.append(
            CellRangerMatrix(
                group=group,
                n_cells=int(group["barcodes"].shape[0]),
                n_features=int(group["genes"].shape[0]),
                index=_text(names),
                columns={"gene_ids": _text(group["genes"])},
                windows=[(genome, 0, int(group["genes"].shape[0]))],
            )
        )
    if not found:
        raise ValueError(f"{name} holds no Cell Ranger count matrix (no /matrix group, and no genome group).")
    if len({matrix.n_cells for matrix in found}) != 1:
        held = ", ".join(f"{m.windows[0][0]}: {m.n_cells}" for m in found)
        raise ValueError(
            f"{name} holds genomes with different cells ({held}). Read one genome with scanpy "
            "(sc.read_10x_h5(path, genome='…')), write it as .h5ad and serve that."
        )
    return found


def _library_rank(directory: str) -> int:
    return 1 if re.search(r"antibody|crispr|multiplexing|peaks|atac|custom", directory, re.IGNORECASE) else 0


def _library_order(left: str, right: str) -> int:
    rank = _library_rank(left) - _library_rank(right)
    return rank if rank else (-1 if left < right else 1 if left > right else 0)


def _by_library(left: tuple[str, str], right: tuple[str, str]) -> int:
    """``(directory, path)`` pairs by their library (``_library_order`` of the directories)."""
    return _library_order(left[0], right[0])


def _js_number(text: str) -> float:
    """``Number(text)`` as JavaScript reads it (NaN when it is no number)."""
    stripped = text.strip()
    if stripped == "":
        return 0.0
    try:
        return float(stripped)
    except ValueError:
        return float("nan")


def _label_order(left: str, right: str) -> int:
    difference = _js_number(left) - _js_number(right)
    if difference == difference and difference != 0:
        return -1 if difference < 0 else 1
    return -1 if left < right else 1 if left > right else 0


def _parse_csv(text: str, label: str) -> tuple[list[str], list[list[str]]]:
    lines = [line for line in re.split(r"\r?\n", text) if line]
    if not lines:
        raise ValueError(f"{label} is empty.")
    header = lines[0].split(",")
    rows = [line.split(",") for line in lines[1:]]
    for index, row in enumerate(rows):
        if len(row) != len(header):
            raise ValueError(
                f"{label} line {index + 2} has {len(row)} {plural(len(row), 'field')} where its header has {len(header)}."
            )
    return header, rows


def outs_matrix_file(files: tuple[str, ...]) -> tuple[str, str] | None:
    """The matrix file of an ``outs/`` folder (the shallowest of the preferred names) and its folder."""
    for wanted in OUTS_MATRIX_FILES:
        found = sorted(
            (path for path in files if path == wanted or path.endswith(f"/{wanted}")),
            key=lambda path: (path.count("/"), path),
        )
        if found:
            return found[0], found[0][: len(found[0]) - len(wanted)]
    return None


def outs_analysis(folder: Path, files: tuple[str, ...], base: str, barcodes: np.ndarray) -> tuple[dict[str, Any], dict[str, Any]]:
    """Cell Ranger's projections and clusterings beside the matrix (``views/tenx-outs.js``): obs columns, obsm."""
    cell_of = {str(barcode): cell for cell, barcode in enumerate(barcodes)}
    n_cells = len(barcodes)
    projections: dict[str, list[tuple[str, str]]] = {}
    clusterings: list[tuple[str, str]] = []
    for path in files:
        if not path.startswith(base):
            continue
        inside = path[len(base) :]
        projection, arc, clusters = _PROJECTION.search(inside), _ARC_PROJECTION.search(inside), _CLUSTERS.search(inside)
        if projection is not None or arc is not None:
            kind, directory = (projection.group(1), projection.group(2)) if projection is not None else (arc.group(2), arc.group(1))  # type: ignore[union-attr]
            projections.setdefault(kind, []).append((directory, path))
        elif clusters is not None:
            clusterings.append((clusters.group(1), path))
    obsm: dict[str, Any] = {}
    for kind, found in projections.items():
        found.sort(key=functools.cmp_to_key(_by_library))
        for position, (directory, path) in enumerate(found):
            header, rows = _parse_csv((folder / path).read_text(), path)
            width = len(header) - 1
            if width < 1:
                continue
            values = np.full((n_cells, width), np.nan, dtype=np.float64)
            for row in rows:
                cell = cell_of.get(row[0])
                if cell is not None:
                    values[cell] = [_js_number(field) for field in row[1:]]
            key = f"X_{kind}" if position == 0 else f"X_{kind}_{re.sub(r'_?[0-9]+_components$', '', directory) or directory}"
            obsm.setdefault(key, values)
    obs: dict[str, Any] = {}
    clusterings.sort(key=functools.cmp_to_key(_by_library))
    for directory, path in clusterings:
        _, rows = _parse_csv((folder / path).read_text(), path)
        labels = sorted(dict.fromkeys(row[1] for row in rows), key=functools.cmp_to_key(_label_order))
        code_of = {label: code for code, label in enumerate(labels)}
        codes = np.full(n_cells, -1, dtype=np.int32)
        for row in rows:
            cell = cell_of.get(row[0])
            if cell is not None:
                codes[cell] = code_of[row[1]]
        obs[directory] = pd.Categorical.from_codes(codes, categories=pd.Index(labels, dtype=object))
    return obs, obsm


def cell_ranger_presented(
    root: Any,
    name: str,
    *,
    analysis: Callable[[np.ndarray], tuple[dict[str, Any], dict[str, Any]]] | None,
) -> Presented | PresentedModalities:
    """A Cell Ranger matrix as the browser presents it: one AnnData, or one modality per feature type or genome.

    ``analysis(barcodes)``: an ``outs/`` folder's obs columns and obsm (None
    for the file alone, which then has no embedding and is refused, saying to
    serve the folder).
    """
    matrices = cell_ranger_matrices(root, name)
    barcodes = _text(matrices[0].group["barcodes"])
    obs_columns, obsm = ({}, {}) if analysis is None else analysis(barcodes)
    if not obsm:
        where = (
            "Its analysis/ folder holds no UMAP or t-SNE projection (Cell Ranger was run with --nosecondary)."
            if analysis is not None
            else "Cell Ranger writes UMAP and t-SNE projections beside it: serve the whole outs/ folder, which holds "
            "this file and analysis/."
        )
        raise ValueError(no_embedding_message(f"{name} is a 10x Cell Ranger count matrix that", where))
    obs = _frame(None, barcodes, len(barcodes), obs_columns)

    def presented(matrix: CellRangerMatrix, start: int, count: int) -> Presented:
        var = _frame(
            None,
            matrix.index[start : start + count],
            count,
            {key: values[start : start + count] for key, values in matrix.columns.items()},
        )
        return Presented(obs=obs, var=var, x=CellRangerX(matrix, start, count), obsm=obsm)

    entries = [(window_name, matrix, start, count) for matrix in matrices for (window_name, start, count) in matrix.windows]
    if len(entries) == 1:
        _, matrix, start, count = entries[0]
        return presented(matrix, start, count)
    return PresentedModalities(
        obs=obs,
        obsm=obsm,
        modalities={window_name: presented(matrix, start, count) for window_name, matrix, start, count in entries},
    )


@dataclass(frozen=True)
class CellRangerX:
    """A Cell Ranger matrix's features ``[start, start + count)`` (a column window: one feature type)."""

    matrix: CellRangerMatrix
    start: int
    count: int

    @property
    def shape(self) -> tuple[int, int]:
        return (self.matrix.n_cells, self.count)
