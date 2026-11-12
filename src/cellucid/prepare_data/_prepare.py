"""
:func:`prepare`: an AnnData, a MuData, or arrays and data frames in; one Cellucid-ready AnnData or MuData out.
"""

from __future__ import annotations

import functools
import os
import warnings
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any, Literal, ParamSpec, TypeVar

import numpy as np
import pandas as pd
from scipy import sparse

from .._anndata_files import AnnDataContents, is_mudata
from .._console import console_print
from .._contracts import (
    _require_dataset_id,
    _require_dataset_name,
    _require_native_boolean,
    _require_prose_text,
)
from .._launch_defaults import launch_defaults
from .._plural import plural
from ..data_mapping import EmbeddingKeySpec, column_gene_names, default_modality
from ..subsample import _require_uint32
from ._data import (
    X_VAR_COMPONENT,
    anndata_components,
    cellucid_record,
    default_expression_of,
    is_anndata,
    is_cell_array,
    split_obsm,
)
from ._inputs import (
    ObsColors,
    held_obs_colors,
    linked_text,
    require_obs_colors,
    resolve_connectivities,
    resolve_embeddings,
    resolve_expression,
    resolve_genes,
    resolve_obs,
    resolve_vector_fields,
    storable_columns,
    stored_expression,
    written_obs_colors,
)
from ._mudata import (
    Modality,
    MuDataLayout,
    _modality_error,
    check_mudata,
    global_var_index,
    modality_obsm_widths,
    placed_rows,
    resolve_modality,
)
from ._outliers import DEFAULT_MIN_POINTS, stored_outlier_quantiles
from ._priority import CellPriority, prefix_levels, stored_order_of, without_positions
from ._store import DEFAULT_ZARR_FORMAT, Elements, OutputTarget, write_store
from ._totals import (
    OVER_ALL_GENES,
    require_counts,
    require_totals_column,
    stored_row_totals,
    total_counts_column,
    total_counts_columns,
    total_counts_over,
    total_counts_over_by_modality,
    with_total_counts,
)

#: The order ``uns['cellucid']['defaults']`` lists its keys in.
_DEFAULTS_ORDER = (
    "embedding1d",
    "embedding2d",
    "embedding3d",
    "modality",
    "geneIdColumn",
    "latentKey",
    "vectorFieldDefault",
    "obsKeys",
    "normalizeEmbeddings",
    "centroidMinPoints",
)


#: The suffix ``out_path`` ends in for each ``format``, and what that writes.
_FORMATS: dict[str, tuple[str, str]] = {
    "h5ad": (".h5ad", "an .h5ad file"),
    "zarr": (".zarr", "a .zarr folder"),
    "zarr.zip": (".zarr.zip", "a .zarr.zip file"),
    "h5mu": (".h5mu", "an .h5mu file"),
}

#: Which ``format`` to choose; every refusal of ``format`` quotes it.
_FORMAT_CHOICES = (
    'format must be "h5ad", "zarr", "zarr.zip" or "h5mu":\n'
    '  "h5ad"      one .h5ad file: for single-file hosts (Zenodo, Figshare, Hugging Face, '
    "S3, a lab server), local files, cellucid serve and Jupyter.\n"
    '  "zarr"      a .zarr folder of small files: for GitHub (no file over 100 MB, no LFS).\n'
    '  "zarr.zip"  the Zarr store in one .zarr.zip file: wherever an .h5ad goes.\n'
    '  "h5mu"      one .h5mu file, for a MuData (data=mdata): wherever an .h5ad goes.'
)

_P = ParamSpec("_P")
_R = TypeVar("_R")


def _requires_format(function: Callable[_P, _R]) -> Callable[_P, _R]:
    """Refuse a call without ``format`` with the choices, not Python's bare message."""

    @functools.wraps(function)
    def checked(*args: _P.args, **kwargs: _P.kwargs) -> _R:
        if "format" not in kwargs:
            raise TypeError(f"prepare() needs format=, the kind of output to write. {_FORMAT_CHOICES}")
        return function(*args, **kwargs)

    return checked


_FORMAT_OF_SUFFIX = {suffix: name for name, (suffix, _) in _FORMATS.items()}


def _known_suffix(name: str) -> str | None:
    """The output suffix ``name`` ends in, the longest that fits."""
    lowered = name.lower()
    for suffix in (".zarr.zip", ".zarr", ".h5ad", ".h5mu"):
        if lowered.endswith(suffix):
            return suffix
    return None


def _resolve_target(
    out_path: object, *, output_format: object, zarr_format: object, force: bool, mudata: bool
) -> OutputTarget:
    """What ``format`` writes, checked against ``out_path``, and whether it may go there."""
    if not (isinstance(output_format, str) and output_format in _FORMATS):
        raise ValueError(f"{_FORMAT_CHOICES}\nGot {output_format!r}.")
    if output_format == "h5mu" and not mudata:
        raise ValueError(
            'format="h5mu" writes a MuData: pass the MuData as data= (read one with '
            'mudata.read_h5mu; pip install "cellucid[mudata]"). An AnnData is written as '
            'format="h5ad", "zarr" or "zarr.zip".'
        )
    if output_format == "h5ad" and mudata:
        raise ValueError(
            'A MuData is written as format="h5mu", "zarr" or "zarr.zip"; an .h5ad holds one '
            "AnnData. To write one modality alone, pass data=mdata.mod['<name>']."
        )
    if not isinstance(out_path, str | os.PathLike):
        raise TypeError("out_path must be a string or os.PathLike path.")
    suffix, what = _FORMATS[output_format]
    given = Path(out_path).expanduser()
    name = given.name
    present = _known_suffix(name)
    if present != suffix or len(name) == len(suffix):
        stem = name if present is None else name[: len(name) - len(present)]
        instead = repr(str(Path(out_path).parent / f"{stem or 'data'}{suffix}"))
        if present is None:
            problem = f"does not end in {suffix}"
        elif present == suffix:
            problem = "is only the suffix"
        else:
            problem = f'ends in {present}, which format="{_FORMAT_OF_SUFFIX[present]}" writes'
        raise ValueError(
            f'format="{output_format}" writes {what}, and out_path {str(out_path)!r} '
            f"{problem}: write {instead} instead."
        )
    kind: Literal["h5ad", "h5mu", "zarr"] = (
        "zarr" if output_format in ("zarr", "zarr.zip") else "h5mu" if output_format == "h5mu" else "h5ad"
    )
    zipped = output_format == "zarr.zip"
    resolved_format: Literal[2, 3] | None
    if kind != "zarr":
        if zarr_format is not None:
            raise ValueError(
                f'zarr_format applies to format="zarr" or "zarr.zip", not format="{output_format}".'
            )
        resolved_format = None
    elif zarr_format is None:
        resolved_format = DEFAULT_ZARR_FORMAT
    elif zarr_format == 2 and type(zarr_format) is int:
        resolved_format = 2
    elif zarr_format == 3 and type(zarr_format) is int:
        resolved_format = 3
    else:
        raise ValueError(f"zarr_format must be 2 or 3, got {zarr_format!r}.")
    path = given.parent.resolve() / name
    if path.exists() or path.is_symlink():
        if not force:
            raise FileExistsError(f"{path} already exists; pass force=True to replace it.")
        if path.is_symlink():
            raise ValueError(f"{path} is a symbolic link; remove it yourself if it should go.")
        if kind == "zarr" and not zipped:
            if not (path.is_dir() and ((path / ".zgroup").is_file() or (path / "zarr.json").is_file())):
                raise ValueError(
                    f"{path} exists and is not a Zarr store, so force=True does not replace "
                    "it; remove it yourself if it should go."
                )
        elif not path.is_file():
            raise ValueError(
                f"{path} is a folder, and format=\"{output_format}\" writes a file; remove it "
                "yourself if it should go."
            )
    return OutputTarget(path=path, kind=kind, zipped=zipped, zarr_format=resolved_format)


def _uns_value(value: Any) -> Any:
    """A JSON-shaped default as anndata stores it: a list becomes an array."""
    if isinstance(value, dict):
        return {key: _uns_value(item) for key, item in value.items()}
    if isinstance(value, list):
        if all(isinstance(item, str) for item in value):
            return np.asarray(value, dtype=object)
        return np.asarray(value, dtype=np.int64)
    return value


@_requires_format
def prepare(
    out_path: str | os.PathLike[str],
    *,
    format: Literal["h5ad", "zarr", "zarr.zip", "h5mu"],
    data: Any = None,
    embeddings: Mapping[str, Any] | None = None,
    dataset_id: str | None = None,
    dataset_name: str | None = None,
    obs: pd.DataFrame | None = None,
    obs_colors: Mapping[str, Sequence[str]] | None = None,
    var: pd.DataFrame | None = None,
    gene_expression: np.ndarray | sparse.spmatrix | sparse.sparray | None = None,
    connectivities: np.ndarray | sparse.spmatrix | sparse.sparray | None = None,
    vector_fields: Mapping[str, Mapping[str, Any]] | None = None,
    modalities: Sequence[str] | None = None,
    gene_identifiers: Sequence[str] | None = None,
    total_counts: str | Mapping[str, str] | None = None,
    total_counts_over: Literal["all genes", "source panel"] | Mapping[str, str] | None = None,
    dataset_description: str | None = None,
    citation_text: str | None = None,
    citation_url: str | None = None,
    source_name: str | None = None,
    source_url: str | None = None,
    embedding_defaults: Mapping[str, EmbeddingKeySpec] | None = None,
    latent_key: str | None = None,
    obs_keys: Sequence[str] | None = None,
    var_gene_id_column: str | None = None,
    vector_field_default: str | None = None,
    normalize_embeddings: bool | None = None,
    centroid_min_points: int | None = None,
    modality: str | None = None,
    priority_seed: int | None = 0,
    zarr_format: Literal[2, 3] | None = None,
    force: bool = False,
) -> Path:
    """Write a Cellucid-ready AnnData (or MuData): a standard file laid out for the browser.

    The output is an ordinary AnnData that ``anndata.read_h5ad`` or
    ``anndata.read_zarr`` reads, written so the Cellucid web app reads it fast
    over HTTP: the expression matrix gene-major (CSC), arrays chunked for range
    requests and compressed with codecs browsers decode natively, and the
    dataset's name, citation and launch defaults in ``uns['cellucid']``. Open
    it with ``cellucid serve``, :func:`~cellucid.show_anndata`, or the web app
    directly. The written contract is the "Cellucid-ready AnnData" page of the
    documentation; cellucid-r's ``cellucid_prepare()`` writes the same.

    The inputs come from ``data`` -- an AnnData or a MuData -- or are passed
    one by one, or both: every component ``data`` holds is its default, and an
    argument given explicitly replaces that one component. The output is
    byte for byte what the same components passed by hand write.

    A MuData is written as a Cellucid-ready MuData (``format="h5mu"``,
    ``"zarr"`` or ``"zarr.zip"``): a standard MuData whose modalities share
    the cells, each modality a Cellucid-ready AnnData (gene-major ``X``,
    chunked and compressed for the browser), the cells' metadata, shared
    embeddings and graph at the root with ``uns['cellucid']``. The written
    contract is the "Cellucid-ready MuData" section of the specification.

    Every input is checked before anything is written (a MuData modality's
    matrix as it is written, one modality at a time), and the output is
    written beside ``out_path`` under a temporary name and moved into place
    only once it is complete.

    Parameters
    ----------
    out_path
        Where to write, ending in the suffix of ``format``: ``.h5ad``,
        ``.zarr``, ``.zarr.zip`` or ``.h5mu``. A leading ``~`` is expanded.
    format
        Required, with no default -- which one depends on where the data goes
        (``"h5mu"`` is the single-file form of a MuData, like ``"h5ad"``;
        ``"zarr"`` and ``"zarr.zip"`` write a MuData as MuData Zarr):

        - ``"zarr"``: a ``.zarr`` folder of small files, for GitHub (a
          repository or its Pages site). GitHub refuses files over 100 MB and
          answers over HTTP/2, 100 requests at a time, so a folder's one
          request per chunk costs nothing there.
        - ``"h5ad"``: one ``.h5ad`` file, for everything else: single-file
          hosts (Zenodo, Figshare, Hugging Face, S3, a lab server), local
          files, ``cellucid serve`` and Jupyter. Its ranges are merged into
          few requests; most such hosts answer over HTTP/1.1, six requests at
          a time, where a folder is several times slower. What opening reads
          comes first in the file (431 KB of a 561,947-cell atlas of 5,103
          genes), which the web app reads with its first request of 512 KiB.
        - ``"zarr.zip"``: the Zarr store in one ``.zarr.zip`` file (entries
          stored, metadata last), wherever an ``.h5ad`` goes.

        Measured on a 561,947-cell atlas, every read the web app makes, in
        total: over GitHub's HTTP/2 (80 ms, 50 MB/s) 8.3 s from the ``.h5ad``,
        8.0 s from the folder, 8.2 s from the ``.zarr.zip``; over HTTP/1.1
        (80 ms, 50 MB/s) 9.0 s, 46.8 s and 9.0 s; from a local file in a
        browser 2.5, 2.6 and 2.7 s. The ``.zarr.zip`` is 13% smaller than the
        ``.h5ad`` (198 against 228 MB). The "Which format?" section of the
        documentation's preparation guide has the full table.
    data
        An ``anndata.AnnData`` or a ``mudata.MuData`` (``pip install
        "cellucid[mudata]"`` to read and build one) that fills every component
        below an explicit argument does not:

        - AnnData: ``obs`` with Scanpy's colors of its columns
          (``uns['<column>_colors']``, as ``obs_colors``); as
          ``gene_expression`` with its ``var``, the matrix the viewer opens
          the AnnData with (Scanpy's rule, ``cellucid.data_mapping.default_expression``):
          ``raw.X`` with ``raw.var`` when ``raw`` holds more genes than
          ``X`` or ``X`` is scaled (``adata.raw = adata`` before the highly
          variable genes were kept and ``sc.pp.scale``), else ``X`` with
          ``var``; without ``X``, ``raw``, else a layer (a message says
          which, and why). A ``gene_expression=`` given in its place goes
          with ``var``; every numeric ``obsm`` array as an embedding
          under its own key, except a ``<field>_<basis>`` array as wide as
          ``X_<basis>``, which is that vector field; ``obsp['connectivities']``;
          and from ``uns['cellucid']`` the identity (``id``, ``name``,
          ``description``, ``citation``, ``source``) and the launch defaults.
          The other matrices (layers, ``raw`` or ``X``), other ``obsp`` matrices, ``obsm`` entries that are
          not numeric arrays, and ``obs``/``var`` columns of no AnnData
          encoding Cellucid writes (a date, a time span) are not written, and
          a message names each one, as cellucid-r names them (the same column
          given as ``obs=`` is refused, saying how to convert it). To leave
          out a component ``data`` holds, remove it from a copy first.
        - MuData (modalities over shared cells, ``axis`` 0): the cells'
          ``mdata.obs``, ``mdata.obsm`` (the shared embeddings and vector
          fields), ``mdata.obsp['connectivities']`` and ``mdata.uns['cellucid']``
          at the root, as above; each modality whole -- the matrix the viewer
          opens it with and its genes (the AnnData rule above: ``raw.X`` with
          ``raw.var`` when ``raw`` holds more genes than ``X`` or ``X`` is
          scaled, the root ``var`` and ``varmap`` then following ``raw``'s
          genes), ``obs``, numeric ``obsm`` arrays and
          ``obsp['connectivities']`` -- as ``mod/<name>``. Scanpy's colors go with the columns they color:
          the root's in ``mdata.uns``, a modality's own in its ``uns``. Columns and arrays mudata copies to the root
          from a modality (``<name>:<column>``) are not written twice: the
          viewer shows every modality's own under that name. The genes are
          each modality's own, so ``gene_expression``, ``var`` and
          ``gene_identifiers`` are refused with a MuData: change
          ``mdata.mod[name]`` first.
    embeddings
        ``{obsm key: (n_cells, k) array}``, at least one (with a MuData, at
        the root or in a modality). Each is stored whole
        as ``obsm[key]`` in float32, every column, at its own coordinates (the
        viewer normalizes them). A latent space for outlier measurement is one
        of them, named by ``latent_key``. A plain ``(n_cells,)`` vector is one
        column; sparse input is densified. A row holding NaN (or an infinity)
        is a cell the embedding does not place (one computed on some of the
        cells): it is written as given, the viewer leaves it out of that
        embedding's views, and in the latent space it gets no outlier
        quantile; at least one cell must have finite values. ``X_umap_1d``/``X_umap_2d``/``X_umap_3d``,
        else the first 1- to 3-column of ``X_umap``, ``X_UMAP``, ``X_wnn.umap``
        and ``X_ref.umap``, are shown by default.
    dataset_id
        The dataset's stable id: 1 to 180 ASCII letters, digits, ``.``, ``_``
        or ``-``, starting with a letter or digit and not ending with ``.``.
        Sessions and community annotations are keyed by it.
    dataset_name
        The name shown for the dataset, also stored as ``uns['title']``.
    obs
        Cell metadata, one row per cell. Its index names the cells (pandas'
        default ``0 .. n-1`` index names them ``"0"`` to ``"n - 1"``);
        repeated names are made unique, as anndata does. Text columns become
        categoricals when a value repeats or is missing (anndata's rule);
        numbers, booleans, categoricals and pandas nullable integers and
        booleans are stored as they are. Infinities are refused.
    obs_colors
        Scanpy's colors of categorical obs columns, ``{column: [color, ...]}``
        (one color a category, in the column's category order), each written
        as ``uns['<column>_colors']``, where Scanpy keeps them and the viewer
        opens the column in them; a MuData modality's own column is named
        ``<modality>:<column>`` and its colors are written in that modality's
        ``uns``. From ``data``: the colors it holds for every column written
        (``uns['<column>_colors']``); one that cannot be written (not a list
        of text, fewer colors than categories, a color that is not hex, a
        column not written as a categorical) is left out, and a message says
        why. Given here, each entry must name a categorical column written
        and give at least one hex color (``"#rrggbb"``, ``"#rgb"`` or
        ``"#rrggbbaa"``, the forms the viewer reads; Scanpy writes hex) a
        category, or it is refused; ``{}`` writes none.
    var
        Gene metadata, one row per column of ``gene_expression``. Its index
        names the genes; with pandas' default index, ``var_gene_id_column``
        names them. Every column is stored.
    gene_expression
        ``(n_cells, n_genes)``, sparse (any scipy format) or dense. Stored as
        ``X``: CSC when sparse, dense otherwise, in float32 -- or float64 when
        the genes it keeps (after ``gene_identifiers``) hold integer counts
        above 2^24, which float32 cannot hold exactly. Values must be finite.
        A stored zero is dropped.
    connectivities
        ``(n_cells, n_cells)`` neighbor graph, stored as
        ``obsp['connectivities']`` (CSR float32). Exactly symmetric, finite,
        non-negative; its diagonal is left out (a self-loop is no edge: Seurat's
        ``RNA_snn`` has 1 there), and a directed graph (``RNA_nn``,
        ``distances``) is refused. From an AnnData ``data=``, the graph the
        viewer's default rule opens with: ``obsp['connectivities']``, else
        Seurat's ``<assay>_snn``.
    vector_fields
        ``{field: {embedding key: (n_cells, k) array}}``, one column per
        column of the embedding. Each array is stored as
        ``obsm['<field>_<basis>']`` for the embedding ``X_<basis>``, scVelo's
        name (``{"velocity": {"X_umap": v}}`` writes ``velocity_umap``), which
        the viewer draws as the field ``<field>`` over that embedding, so the
        embedding's key must start with ``X_``.
    modalities
        With a MuData, the modalities to write, in this order (names of
        ``mdata.mod``); None writes every one, in ``mdata.mod`` order.
        Several are written as a Cellucid-ready MuData of those modalities
        (the root ``var`` and ``varmap`` hold theirs alone). One is written as
        a plain Cellucid-ready AnnData (``format="h5ad"``, ``"zarr"`` or
        ``"zarr.zip"``): the cells' ``obs``, embeddings, graph and
        ``uns['cellucid']`` from the root, and that modality's ``X``, ``var``
        and ``obsm`` arrays, its rows in the MuData's cell order -- so it must
        have measured every cell -- with ``gene_expression``, ``var`` and
        ``gene_identifiers`` then applying to it. cellucid-r's
        ``cellucid_prepare(modalities = )`` takes the same names.
    gene_identifiers
        Keep only these genes, in this order, matched against the names the
        viewer shows (``var_gene_id_column``, else the var index).
    total_counts
        The ``obs`` column holding each cell's total count over **all** its
        genes, computed before any gene was selected (scanpy's
        ``sc.pp.calculate_qc_metrics(adata)`` run before subsetting to highly
        variable genes writes ``total_counts``; Seurat's ``nCount_RNA``).
        The viewer normalizes a gene by it (``x * 10,000 / total``) from this
        one column, instead of reading every gene once to sum them, and the
        totals are of the whole cell where the file stores a gene panel.
        Recorded as ``uns['cellucid']['totalCounts'] = {"field": <column>,
        "over": "all genes"}``; refused unless ``X`` holds counts (every
        value a whole number of at least 0) and every cell's total is finite
        and at least the sum of its stored genes. With a MuData, a mapping
        of modality to column: ``{"rna": "total_counts"}`` (a column of the
        MuData's ``obs``) or ``{"rna": "rna:total_counts"}`` (one of
        ``mdata.mod['rna'].obs``). When None and ``X`` holds counts, each
        cell's sum over the genes stored is written as
        ``obs['total_counts']``, recorded ``"over": "stored genes"`` -- unless
        ``obs`` already has a column of that name, which is kept as given
        (and recorded only when it is exactly those sums).
    total_counts_over
        What the column ``total_counts`` names covers: ``"all genes"`` (when
        None), every gene measured; or ``"source panel"``, every gene of the
        source's gene panel -- a targeted or integration panel, more genes
        than are stored but not the whole transcriptome (the built-in atlas
        samples' 8,192 integration genes). Recorded as the ``over`` of
        ``uns['cellucid']['totalCounts']``, and checked as ``"all genes"``
        is. With a MuData, one value for every modality ``total_counts``
        names, or a mapping of modality to value (``{"rna": "source
        panel"}``; a modality left out covers all genes). Refused without
        ``total_counts``.
    dataset_description
        One paragraph about the dataset.
    citation_text, citation_url
        How to cite the dataset; the address needs the text.
    source_name, source_url
        Where the data comes from; the address needs the name.
    embedding_defaults
        The embedding each dimension opens on: ``{"2d": "X_tsne"}`` (its
        first columns) or ``{"2d": ("X_pca", [2, 0])}`` (0-based columns).
    latent_key
        The embedding (an ``embeddings`` key) the viewer's outlier filter
        measures in: a latent space such as scVI or PCA, never a 2D/3D
        layout. The outlier quantiles of every categorical obs column (each
        cell's distance to its category's center, as the share of its
        category at that distance or closer) are computed here and stored in
        ``uns['cellucid']['outlierQuantiles']``, so the viewer filters
        outliers at once instead of reading the whole latent space to
        compute them.
    obs_keys
        The obs columns the viewer lists, in this order; every column is
        stored either way.
    var_gene_id_column
        The var column whose values name the genes in the viewer. Omitted,
        the var index names them; Data mapping can choose a column later.
    vector_field_default
        The ``vector_fields`` field the overlay opens on; one of its
        embeddings must be shown when the file opens (a chosen or default
        embedding, or a 1-3 column array that is not the latent space).
    normalize_embeddings
        False shows embeddings at their own coordinates.
    centroid_min_points
        The fewest cells a category needs to get outlier quantiles (10 when
        None).
    modality
        With a MuData, the modality the viewer opens on (one of
        ``mdata.mod``); ``var_gene_id_column`` names a column of its ``var``.
        When None the viewer opens on the first gene-expression modality
        (``rna``, ``gex``, ``gene expression``, ``gene_expression``, ``sct``,
        any case), else the first. Refused for an AnnData, which has one
        feature space.
    priority_seed
        The seed of the cells' priority order, an integer in
        ``[0, 2**32 - 1]``. The cells are stored by ascending priority -- the
        ``uniform`` rank of the subsampling specification
        (:func:`~cellucid.subsample_cells`) with this seed -- so every
        subsample Cellucid opens with this seed is the first rows of every
        array, and reads nothing of the cells it leaves out. Each row's
        original position is kept in ``obs['cellucid_original_position']``;
        :func:`~cellucid.restore_original_order` puts the cells back in the
        order given. None stores the cells in the order given, without a
        priority.
    zarr_format
        2 or 3, for ``format="zarr"`` or ``"zarr.zip"``; Zarr v2 when None.
        Refused with ``format="h5ad"`` and ``"h5mu"``.
    force
        Replace an existing output of the same kind. A folder that is not a
        Zarr store is never replaced.

    The launch defaults (``embedding_defaults`` to ``modality``)
    are written to ``uns['cellucid']['defaults']`` only when given; a choice
    left at None stays the viewer's own.

    Returns
    -------
    Path
        The written file or folder.

    Raises
    ------
    TypeError
        A required input is missing (``dataset_id`` or ``dataset_name`` that
        ``data`` does not name, ``embeddings`` without ``data``) or an input
        has the wrong type.
    ValueError
        An input is invalid or the inputs do not fit together (``out_path``
        not ending in the suffix of ``format``, a MuData of several
        modalities with ``format="h5ad"``, an AnnData with ``format="h5mu"``,
        ``zarr_format`` with ``"h5ad"`` or ``"h5mu"``, ...).
    FileExistsError
        ``out_path`` exists and ``force`` is False.
    """
    force = _require_native_boolean(force, label="force")
    if priority_seed is not None:
        priority_seed = _require_uint32(priority_seed, "priority_seed")
    explicit = {
        "embeddings": embeddings,
        "dataset_id": dataset_id,
        "dataset_name": dataset_name,
        "obs": obs,
        "obs_colors": None if obs_colors is None else require_obs_colors(obs_colors),
        "var": var,
        "gene_expression": gene_expression,
        "connectivities": connectivities,
        "vector_fields": vector_fields,
        "gene_identifiers": gene_identifiers,
        "total_counts": total_counts,
        "total_counts_over": total_counts_over,
        "dataset_description": dataset_description,
        "citation_text": citation_text,
        "citation_url": citation_url,
        "source_name": source_name,
        "source_url": source_url,
        "embedding_defaults": embedding_defaults,
        "latent_key": latent_key,
        "obs_keys": obs_keys,
        "var_gene_id_column": var_gene_id_column,
        "vector_field_default": vector_field_default,
        "normalize_embeddings": normalize_embeddings,
        "centroid_min_points": centroid_min_points,
        "modality": modality,
    }
    given = {key: value for key, value in explicit.items() if value is not None}
    if "total_counts" in given and "total_counts_over" not in given:
        # A column named here replaces the data's, and what the data said its
        # column covers goes with that column.
        given["total_counts_over"] = None
    if data is not None and not is_mudata(data) and not is_anndata(data):
        raise TypeError(
            f"data must be an anndata.AnnData or a mudata.MuData, got {type(data).__name__}. "
            "Pass arrays and data frames as embeddings=, obs=, gene_expression= and the rest."
        )
    positions: np.ndarray | None = None
    if data is not None:
        # An object read from a file in a verified priority order holds its
        # cells in that stored order: each row's original position is known,
        # and the way back to the original order and on to the file's is one
        # permutation (CellPriority.of), so a file prepared from it names
        # every cell as that file did -- and one in the same order is written
        # as it is, with no copy of the object.
        stored = stored_order_of(data)
        if stored.notice is not None:
            warnings.warn(stored.notice, UserWarning, stacklevel=3)
        positions = stored.positions
    selected = _selected_modalities(data, modalities)
    notes: list[tuple[str, str]] = []
    multimodal = data is not None and is_mudata(data) and (selected is None or len(selected) > 1)
    target = _resolve_target(
        out_path, output_format=format, zarr_format=zarr_format, force=force, mudata=multimodal
    )
    if multimodal:
        refused = [key for key in ("gene_expression", "var", "gene_identifiers") if key in given]
        if refused:
            raise ValueError(
                f"{refused[0]} describes the genes of one matrix, and a MuData's genes are each "
                "modality's own. Change the modality before prepare, for example "
                "mdata.mod['rna'] = mdata.mod['rna'][:, genes].copy(); mdata.update(), or "
                "write one modality with modalities=['rna']."
            )
        names = check_mudata(data, selected)
        every = [str(name) for name in data.mod]
        arguments = {**_mudata_root_components(data, every, notes), **given, "priority_seed": priority_seed}
        _say_not_written(notes, given)
        _require_identity(arguments, data_given=True)
        write_store(target, _mudata_layout(data, names, arguments, positions))
    else:
        if selected is not None:
            if "modality" in given:
                raise ValueError(
                    "modality names the modality the viewer opens on, and modalities= writes "
                    "one, which is stored as an AnnData with no modalities."
                )
            components = _single_modality_components(data, selected[0], notes)
        else:
            components = anndata_components(data, notes) if data is not None else {}
        x_var = components.pop(X_VAR_COMPONENT, None)
        if x_var is not None and "gene_expression" in given and "var" not in given:
            # An explicit matrix replaces the one chosen (raw) and goes with X's genes.
            components["var"] = x_var
        arguments = {**components, **given, "priority_seed": priority_seed}
        _say_not_written(notes, given)
        _require_identity(arguments, data_given=data is not None)
        if "embeddings" not in given and data is None:
            raise TypeError(
                "prepare() needs embeddings=: at least one {obsm key: (n_cells, k) array}, "
                "for example embeddings={'X_umap': adata.obsm['X_umap']}, or data= an "
                "AnnData or MuData that holds them."
            )
        if "embeddings" not in given and not arguments["embeddings"]:
            raise ValueError(
                "data holds no embedding: its obsm has no numeric array. Compute one, for "
                "example sc.pp.neighbors(adata); sc.tl.umap(adata), or pass embeddings=."
            )
        write_store(target, _anndata_elements(arguments, positions))
    return target.path


def _say_not_written(notes: list[tuple[str, str]], given: Mapping[str, Any]) -> None:
    """Say what of ``data`` is not written, for each component an explicit argument does not replace."""
    for component, note in notes:
        if component not in given:
            console_print(note)


def _require_identity(arguments: dict[str, Any], *, data_given: bool) -> None:
    """``dataset_id`` and ``dataset_name`` are given, or ``data`` names itself."""
    for key, record in (("dataset_id", "id"), ("dataset_name", "name")):
        if arguments.get(key) is None:
            where = (
                f"data's uns['cellucid'] has no {record!r}, so pass {key}="
                if data_given
                else f"pass {key}="
            )
            raise TypeError(f"prepare() needs {key}: {where}.")


def _cellucid_identity(arguments: Mapping[str, Any]) -> dict[str, Any]:
    """``uns['cellucid']`` without its defaults: schema, id, name, description, citation, source."""
    dataset_id = _require_dataset_id(arguments["dataset_id"], label="dataset_id")
    dataset_name = _require_dataset_name(arguments["dataset_name"], label="dataset_name")
    description = arguments.get("dataset_description")
    if description is not None:
        description = _require_prose_text(description, label="dataset_description", allow_empty=True)
    citation = linked_text(
        arguments.get("citation_text"),
        arguments.get("citation_url"),
        text_label="citation_text",
        url_label="citation_url",
        text_key="text",
    )
    source = linked_text(
        arguments.get("source_name"),
        arguments.get("source_url"),
        text_label="source_name",
        url_label="source_url",
        text_key="name",
    )
    cellucid: dict[str, Any] = {"schema": 1, "id": dataset_id, "name": dataset_name}
    if description is not None:
        cellucid["description"] = description
    if citation is not None:
        cellucid["citation"] = citation
    if source is not None:
        cellucid["source"] = source
    return cellucid


def _checked_defaults(contents: AnnDataContents, arguments: Mapping[str, Any]) -> dict[str, Any]:
    """The launch defaults ``uns['cellucid']['defaults']`` holds, checked against ``contents``."""
    embedding_defaults = arguments.get("embedding_defaults")
    defaults = launch_defaults(
        contents,
        latent_key=arguments.get("latent_key"),
        var_gene_id_column=arguments.get("var_gene_id_column"),
        normalize_embeddings=arguments.get("normalize_embeddings"),
        centroid_min_points=arguments.get("centroid_min_points"),
        obs_keys=arguments.get("obs_keys"),
        embedding_defaults=embedding_defaults,
        knn_key=None,
        expression_source=None,
        vector_field_default=arguments.get("vector_field_default"),
        modality=arguments.get("modality"),
    )
    # An embedding default named by its key alone is stored as the key: the
    # viewer then plots its first columns, as it does for any key.
    for slot, spec in (embedding_defaults or {}).items():
        if isinstance(spec, str):
            defaults[f"embedding{slot}"] = spec
    return defaults


def _with_defaults(
    cellucid: dict[str, Any],
    defaults: Mapping[str, Any],
    *,
    obs_columns: dict[str, Any],
    latent: np.ndarray | None,
    total_counts: dict[str, Any] | None,
    priority: dict[str, Any] | None,
) -> dict[str, Any]:
    """``cellucid`` with its launch defaults, for a latent space its outlier
    quantiles, the record of the cells' totals (:mod:`._totals`), and its priority."""
    if defaults:
        cellucid["defaults"] = _uns_value({key: defaults[key] for key in _DEFAULTS_ORDER if key in defaults})
    if latent is not None:
        cellucid["outlierQuantiles"] = stored_outlier_quantiles(
            obs_columns,
            defaults["latentKey"],
            latent,
            defaults.get("centroidMinPoints", DEFAULT_MIN_POINTS),
        )
    if total_counts is not None:
        cellucid["totalCounts"] = total_counts
    if priority is not None:
        cellucid["priority"] = priority
    return cellucid


def _require_vector_field_default(field: object, fields: Sequence[str]) -> None:
    if field is not None and field not in fields:
        if not fields:
            raise ValueError("vector_field_default names a field, and no vector_fields were supplied.")
        raise ValueError(
            f"vector_field_default must name one of vector_fields: {list(fields)}; got {field!r}."
        )


def _anndata_elements(arguments: Mapping[str, Any], given: np.ndarray | None = None) -> Elements:
    """One Cellucid-ready AnnData's elements, every component checked.

    ``given`` is each row's original position when the components are an
    object's in a verified stored order (:func:`prepare`), None when each row
    is its own.
    """
    cellucid = _cellucid_identity(arguments)
    obsm, n_cells = resolve_embeddings(arguments["embeddings"])
    vectors, fields = resolve_vector_fields(arguments.get("vector_fields"), obsm, n_cells)
    obs_index, obs_columns = resolve_obs(arguments.get("obs"), n_cells)
    obs_columns = without_positions(obs_columns, "obs", verified=given is not None)

    gene_expression = arguments.get("gene_expression")
    var_gene_id_column = arguments.get("var_gene_id_column")
    if gene_expression is None:
        supplied = [
            name for name in ("var", "gene_identifiers", "var_gene_id_column") if arguments.get(name) is not None
        ]
        if supplied:
            raise ValueError(
                f"{supplied[0]} describes the genes of gene_expression, and no "
                "gene_expression was supplied."
            )
        expression = None
        genes = None
        var_index = np.asarray([], dtype=object)
        var_columns: dict[str, Any] = {}
    else:
        expression = resolve_expression(gene_expression, n_cells)
        genes, var_index, var_columns = resolve_genes(
            int(expression.shape[1]),
            arguments.get("var"),
            var_gene_id_column,
            arguments.get("gene_identifiers"),
        )

    connectivities = arguments.get("connectivities")
    obsp = {} if connectivities is None else {
        "connectivities": resolve_connectivities(connectivities, n_cells)
    }

    # Every per-cell element in the file's order (a priority order, or the
    # original one for an object in a stored order), before anything is
    # derived from it (the outlier quantiles), so the file holds what its rows give.
    priority = CellPriority.of(n_cells, arguments.get("priority_seed"), given)
    levels: np.ndarray | None = None
    if priority is not None:
        obsm = priority.rows_of_mapping(obsm)
        vectors = priority.rows_of_mapping(vectors)
        obs_index = priority.rows(obs_index)
        obs_columns = priority.obs_columns(obs_columns)
        obsp = {key: priority.graph(graph) for key, graph in obsp.items()}
        if priority.recorded:
            levels = prefix_levels(n_cells)
    order, new_row = (None, None) if priority is None else priority.rows_moved
    # Gene-major in the file's order, in one pass (:mod:`._gene_major`).
    x = None if expression is None else stored_expression(expression, genes=genes, order=order, new_row=new_row, levels=levels)
    del expression
    # Each cell's total count, of the rows as written (:mod:`._totals`).
    obs_columns, totals = with_total_counts(
        obs_columns,
        stored_row_totals(x, n_cells),
        total_counts_column(arguments.get("total_counts")),
        cell_names=obs_index,
        has_matrix=x is not None,
        over=total_counts_over(arguments.get("total_counts_over"), named=arguments.get("total_counts") is not None),
    )

    _require_vector_field_default(arguments.get("vector_field_default"), fields)
    all_obsm = {**obsm, **vectors}
    latent_key = arguments.get("latent_key")
    if latent_key is not None and latent_key not in obsm:
        raise ValueError(f"latent_key must name one of the embeddings: {list(obsm)}; got {latent_key!r}.")
    defaults = _checked_defaults(
        AnnDataContents(
            n_obs=n_cells,
            n_vars=len(var_index),
            obs_columns=frozenset(obs_columns),
            var_columns=frozenset(var_columns),
            obsm={key: int(values.shape[1]) for key, values in all_obsm.items()},
            obsp=frozenset(obsp),
            layers=frozenset(),
            has_x=x is not None,
            raw_var_columns=None,
        ),
        arguments,
    )
    cellucid = _with_defaults(
        cellucid,
        defaults,
        obs_columns=obs_columns,
        latent=None if latent_key is None else obsm[latent_key],
        total_counts=totals,
        priority=None if priority is None or not priority.recorded else priority.record(levels),
    )
    colors, _ = written_obs_colors(arguments.get("obs_colors"), obs_columns)
    return Elements(
        x=x,
        obs_index=obs_index,
        obs=obs_columns,
        var_index=var_index,
        var=var_columns,
        obsm=all_obsm,
        obsp=obsp,
        uns={"title": cellucid["name"], "cellucid": cellucid, **colors},
    )


def _selected_modalities(data: Any, modalities: object) -> list[str] | None:
    """The modalities ``modalities=`` chooses, checked against ``data``; None for every one."""
    if modalities is None:
        return None
    if data is None:
        raise ValueError("modalities chooses among the modalities of data, and no data was supplied.")
    if not is_mudata(data):
        raise ValueError(
            "modalities chooses among a MuData's modalities, and data is an AnnData, which has "
            "one feature space."
        )
    if isinstance(modalities, str | bytes) or not isinstance(modalities, Sequence):
        raise TypeError("modalities must be a sequence of modality names, such as ['rna', 'prot'].")
    chosen = list(modalities)
    if not all(isinstance(name, str) for name in chosen):
        raise TypeError("modalities must hold modality names (strings).")
    if not chosen or len(set(chosen)) != len(chosen):
        raise ValueError("modalities must name at least one modality of data, each once.")
    available = [str(name) for name in data.mod]
    unknown = [name for name in chosen if name not in available]
    if unknown:
        raise ValueError(
            f"modalities names {unknown}, which data does not hold. Its modalities: {available}."
        )
    return chosen


def _single_modality_components(mdata: Any, name: str, notes: list[tuple[str, str]] | None = None) -> dict[str, Any]:
    """One modality of a MuData as an AnnData's components (``modalities=[name]``).

    The root's cells, obs, embeddings, vector fields, graph and identity, and
    the modality's matrix, genes and ``obsm`` arrays, its rows placed in the
    MuData's cell order. A modality that did not measure every cell is
    refused: an AnnData has a row for every cell.
    """
    every = [str(item) for item in mdata.mod]
    check_mudata(mdata, [name])
    modality = mdata.mod[name]
    obsmap = np.asarray(mdata.obsmap[name])
    measured = int(np.count_nonzero(obsmap))
    if measured != int(mdata.n_obs):
        raise ValueError(
            f"Modality {name!r} measured {measured:,} of the {int(mdata.n_obs):,} {plural(int(mdata.n_obs), 'cell')}, and one "
            "modality is written as an AnnData, which has a value for every cell. Write it with "
            f"another modality as a MuData (modalities=[..., {name!r}]), or subset the MuData to "
            "its cells first."
        )
    rows = obsmap.astype(np.int64) - 1
    in_order = bool(np.array_equal(rows, np.arange(rows.size)))

    def placed(values: Any) -> Any:
        return values if in_order else values[rows]

    components = _mudata_root_components(mdata, every, notes)
    own = {key: placed(value) for key, value in modality.obsm.items() if is_cell_array(value)}
    clash = sorted(set(own) & set(components["embeddings"]))
    if clash:
        raise ValueError(
            f"Modality {name!r} and the MuData's root both hold obsm {clash}; an AnnData has one "
            "obsm. Rename one of them first."
        )
    embeddings, fields = split_obsm({**components["embeddings"], **own})
    components["embeddings"] = embeddings
    if fields:
        merged = {**components.get("vector_fields", {})}
        for field, arrays in fields.items():
            merged.setdefault(field, {}).update(arrays)
        components["vector_fields"] = merged
    # The matrix the viewer opens the modality with (default_expression_of):
    # raw when it holds more genes than X or X is scaled, else X.
    where = f"mdata.mod[{name!r}]"
    expression, chosen, genes, chosen_note = default_expression_of(modality, where=where)
    if chosen is not None:
        components["gene_expression"] = placed(chosen)
        var, note = storable_columns(genes, f"{where}.raw.var" if expression == "raw" else f"{where}.var")
        components["var"] = var
        if note is not None and notes is not None:
            notes.append(("var", note))
        if expression not in ("X", None) and modality.X is not None:
            components[X_VAR_COMPONENT], _ = storable_columns(modality.var, f"{where}.var")
    if chosen_note is not None and notes is not None:
        notes.append(("gene_expression", chosen_note))
    if "connectivities" not in components and "connectivities" in modality.obsp:
        graph = modality.obsp["connectivities"]
        components["connectivities"] = graph if in_order else graph[rows][:, rows]
    return components


def _copied_from_modality(key: str, names: Sequence[str], keys_of: Callable[[str], Any]) -> bool:
    """Whether ``key`` is mudata's root copy ``<modality>:<key>`` of a modality's own entry."""
    name, separator, inner = key.partition(":")
    return bool(separator) and name in names and inner in keys_of(name)


def _mudata_root_components(mdata: Any, names: Sequence[str], notes: list[tuple[str, str]] | None = None) -> dict[str, Any]:
    """What a MuData holds at its root, as ``prepare()`` takes it (:func:`._data.split_obsm`).

    ``names`` are every modality of ``mdata``: mudata's root copies of any of
    them (``<name>:<column>``) and its boolean masks are not written, whichever
    modalities are.
    """
    obs = mdata.obs
    copied = [
        column
        for column in obs.columns
        if _copied_from_modality(str(column), names, lambda name: mdata.mod[name].obs.columns)
    ]
    root_obsm = {
        key: value
        for key, value in mdata.obsm.items()
        if not (key in names and np.asarray(value).dtype == bool)
        and not _copied_from_modality(str(key), names, lambda name: mdata.mod[name].obsm)
    }
    embeddings, fields = split_obsm(root_obsm)
    stored_obs, note = storable_columns(obs.drop(columns=copied), "mdata.obs")
    if note is not None and notes is not None:
        notes.append(("obs", note))
    # Scanpy's colors of the root's columns, and of each modality's own (as
    # `<modality>:<column>`), where Scanpy run on the modality keeps them.
    colors = held_obs_colors(mdata.uns, list(stored_obs.columns))
    for name in names:
        colors.update(held_obs_colors(mdata.mod[name].uns, list(mdata.mod[name].obs.columns), prefix=f"{name}:"))
    components: dict[str, Any] = {
        "obs": stored_obs,
        "obs_colors": ObsColors(colors=colors, explicit=False),
        "embeddings": embeddings,
    }
    if fields:
        components["vector_fields"] = fields
    if "connectivities" in mdata.obsp:
        components["connectivities"] = mdata.obsp["connectivities"]
    components.update(cellucid_record(mdata.uns, what="data"))
    return components


def _modality_contents(modality: Modality) -> AnnDataContents:
    """What one modality holds, as the launch defaults are checked against it."""
    return AnnDataContents(
        n_obs=len(modality.obs_index),
        n_vars=len(modality.var_index),
        obs_columns=frozenset(modality.obs),
        var_columns=frozenset(modality.var),
        obsm={key: int(values.shape[1]) for key, values in modality.obsm.items()},
        obsp=frozenset(modality.obsp),
        layers=frozenset(),
        has_x=modality.matrix is not None,
        raw_var_columns=None,
    )


def _written_var_axis(
    mdata: Any, modalities: Mapping[str, Modality]
) -> tuple[np.ndarray, dict[str, Modality]]:
    """The root ``var`` and ``varmap`` of the genes written.

    mudata's root var is every modality's features in modality order, and a
    modality's varmap places them there (1-based, 0 for another modality's).
    Made from what is written when that is not the MuData's own axis: a
    selection of modalities (``modalities=``), or a modality written from
    its ``raw`` (``Modality.genes_from_raw``: raw's genes, not its var's).
    """
    from dataclasses import replace

    names = [
        np.asarray(item.var_index, dtype=object) if item.genes_from_raw
        else np.asarray([str(value) for value in mdata.mod[name].var_names], dtype=object)
        for name, item in modalities.items()
    ]
    var_index = np.concatenate(names) if names else np.asarray([], dtype=object)
    placed: dict[str, Modality] = {}
    offset = 0
    for (name, item), own in zip(modalities.items(), names, strict=True):
        varmap = np.zeros(len(var_index), dtype=np.uint32)
        count = len(own)
        varmap[offset : offset + count] = np.arange(1, count + 1, dtype=np.uint32)
        offset += count
        placed[name] = replace(item, varmap=varmap)
    return var_index, placed


def _modality_totals(
    modality: Modality, named: str | None, root_obs: Mapping[str, Any], over: str
) -> tuple[Modality, dict[str, str] | None]:
    """One modality with its cells' totals in its ``obs`` (:mod:`._totals`),
    and their record; ``named`` the column ``total_counts=`` names for it: a
    column of the MuData's ``obs``, or ``<modality>:<column>`` of its own,
    covering ``over`` (:func:`._totals.total_counts_over`)."""
    from dataclasses import replace

    n_rows = len(modality.obs_index)
    matrix = None if modality.matrix is None else modality.matrix()
    try:
        sums = stored_row_totals(matrix, n_rows)
        del matrix
        own = f"{modality.name}:"
        if named is not None and not named.startswith(own):
            if named not in root_obs:
                raise ValueError(
                    f"total_counts names column {named!r}, which neither the MuData's obs nor "
                    f"mdata.mod[{modality.name!r}].obs ({own}<column>) holds."
                )
            counts = require_counts(sums, has_matrix=modality.matrix is not None)
            # The root column at each of the modality's rows (obsmap is 1-based).
            rows = np.empty(n_rows, dtype=np.int64)
            placed = np.flatnonzero(modality.obsmap)
            rows[modality.obsmap[placed].astype(np.int64) - 1] = placed
            values = root_obs[named]
            if not isinstance(values, np.ndarray):
                raise ValueError(f"total_counts names obs column {named!r}, which is not a column of numbers.")
            require_totals_column(values[rows], named, counts, modality.obs_index)
            return modality, {"field": named, "over": over}
        column = None if named is None else named[len(own) :]
        obs, record = with_total_counts(
            modality.obs,
            sums,
            column,
            cell_names=modality.obs_index,
            field_name=lambda name: f"{own}{name}",
            has_matrix=modality.matrix is not None,
            over=over,
        )
    except (TypeError, ValueError) as error:
        raise _modality_error(modality.name, error) from error
    return replace(modality, obs=obs), record


def _mudata_layout(mdata: Any, names: Sequence[str], arguments: Mapping[str, Any], given: np.ndarray | None = None) -> MuDataLayout:
    """One Cellucid-ready MuData, every component but the modalities' matrices checked.

    ``given`` is each cell's original position when the MuData is in a
    verified stored order (:func:`prepare`), None when each row is its own.
    """
    cellucid = _cellucid_identity(arguments)
    n_cells = int(mdata.n_obs)
    supplied = arguments.get("embeddings") or {}
    if supplied:
        obsm, rows = resolve_embeddings(supplied)
        if rows != n_cells:
            raise ValueError(
                f"The embeddings have {rows:,} {plural(rows, 'row')}, and the MuData has {n_cells:,} {plural(n_cells, 'cell')} "
                "(mdata.n_obs): a root embedding has one row per cell of the MuData."
            )
    else:
        obsm = {}
    vectors, fields = resolve_vector_fields(arguments.get("vector_fields"), obsm, n_cells)
    given_obs = arguments.get("obs")
    if given_obs is not None and len(given_obs) != n_cells:
        raise ValueError(
            f"obs has {len(given_obs):,} {plural(len(given_obs), 'row')}, and the MuData has {n_cells:,} {plural(n_cells, 'cell')} (mdata.n_obs): "
            "obs has one row per cell of the MuData."
        )
    obs_index, obs_columns = resolve_obs(given_obs, n_cells)
    obs_columns = without_positions(obs_columns, "obs", verified=given is not None)
    connectivities = arguments.get("connectivities")
    obsp = {} if connectivities is None else {
        "connectivities": resolve_connectivities(connectivities, n_cells)
    }
    masks = {
        name: np.asarray(mdata.obsm[name], dtype=bool)
        for name in names
        if name in mdata.obsm and np.asarray(mdata.obsm[name]).dtype == bool
    }
    modalities = {name: resolve_modality(mdata, name, stored=given is not None) for name in names}
    named = total_counts_columns(arguments.get("total_counts"))
    over = total_counts_over_by_modality(arguments.get("total_counts_over"), named)
    unknown = sorted(set(named) - set(names))
    if unknown:
        raise ValueError(f"total_counts names {unknown}, which are not modalities written: {list(names)}.")
    totals: dict[str, dict[str, str]] = {}
    for name in names:
        modalities[name], record = _modality_totals(
            modalities[name], named.get(name), obs_columns, over.get(name, OVER_ALL_GENES)
        )
        if record is not None:
            totals[name] = record
    if list(names) == [str(name) for name in mdata.mod] and not any(item.genes_from_raw for item in modalities.values()):
        var_index = global_var_index(mdata)
    else:
        var_index, modalities = _written_var_axis(mdata, modalities)
    # The cells in priority order, each modality's rows by its cells' global
    # priority, so a prefix of the cells is a prefix of every modality.
    priority = CellPriority.of(n_cells, arguments.get("priority_seed"), given)
    levels: np.ndarray | None = None
    if priority is not None:
        obsm = priority.rows_of_mapping(obsm)
        vectors = priority.rows_of_mapping(vectors)
        obs_index = priority.rows(obs_index)
        obs_columns = priority.obs_columns(obs_columns)
        obsp = {key: priority.graph(graph) for key, graph in obsp.items()}
        masks = {name: priority.rows(mask) for name, mask in masks.items()}
        if priority.recorded:
            levels = prefix_levels(n_cells)
        modalities = {name: priority.modality(item, levels) for name, item in modalities.items()}
    if not obsm and not any(modality.obsm for modality in modalities.values()):
        raise ValueError(
            "This MuData holds no embedding: neither mdata.obsm nor any modality's obsm has "
            "a numeric array. Compute one, for example sc.pp.neighbors(mdata.mod['rna']); "
            "sc.tl.umap(mdata.mod['rna']), or pass embeddings=."
        )

    _require_vector_field_default(arguments.get("vector_field_default"), fields)
    all_obsm = {**obsm, **vectors}
    placed_widths = modality_obsm_widths(modalities)
    latent_key = arguments.get("latent_key")
    latent: np.ndarray | None = None
    if latent_key is not None:
        if latent_key in obsm:
            latent = obsm[latent_key]
        elif latent_key in placed_widths:
            # A cell the modality did not measure has a NaN row: it gets no
            # outlier quantile, and the other cells are measured without it.
            name, _, key = latent_key.partition(":")
            modality = modalities[name]
            latent = placed_rows(modality.obsm[key], modality.obsmap)
        else:
            raise ValueError(
                f"latent_key must name one of the embeddings: {[*obsm, *placed_widths]}; got {latent_key!r}."
            )
    chosen = arguments.get("modality")
    opened = modalities.get(chosen) if isinstance(chosen, str) else None
    if opened is None:
        opened = modalities[default_modality(names) or names[0]]
    gene_column = arguments.get("var_gene_id_column")
    opened_var = mdata.mod[opened.name].var
    if gene_column is not None and gene_column in opened_var.columns:
        # The column names the genes of the modality the viewer opens on, by
        # the rule an AnnData's var_gene_id_column is checked with (an absent
        # one is named by the launch-default check below).
        try:
            column_gene_names(opened_var, gene_column)
        except ValueError as error:
            raise _modality_error(opened.name, ValueError(f"var_gene_id_column: {error}")) from error
    defaults = _checked_defaults(
        AnnDataContents(
            n_obs=n_cells,
            n_vars=len(opened.var_index),
            obs_columns=frozenset(
                [
                    *obs_columns,
                    *(f"{name}:{column}" for name, item in modalities.items() for column in item.obs),
                    *(f"{name}:measured" for name in masks),
                ]
            ),
            var_columns=frozenset(opened.var),
            obsm={**{key: int(values.shape[1]) for key, values in all_obsm.items()}, **placed_widths},
            obsp=frozenset(
                [*obsp, *(f"{name}:{key}" for name, item in modalities.items() for key in item.obsp)]
            ),
            layers=frozenset(),
            has_x=True,
            raw_var_columns=None,
            modalities={name: _modality_contents(item) for name, item in modalities.items()},
        ),
        arguments,
    )
    cellucid = _with_defaults(
        cellucid,
        defaults,
        obs_columns=obs_columns,
        latent=latent,
        total_counts=totals or None,
        priority=None if priority is None or not priority.recorded else priority.record(levels),
    )
    from dataclasses import replace

    colors, modality_colors = written_obs_colors(
        arguments.get("obs_colors"), obs_columns, {name: item.obs for name, item in modalities.items()}
    )
    modalities = {name: replace(item, colors=modality_colors[name]) for name, item in modalities.items()}
    return MuDataLayout(
        obs_index=obs_index,
        obs=obs_columns,
        obsm=all_obsm,
        masks=masks,
        obsp=obsp,
        uns={"title": cellucid["name"], "cellucid": cellucid, **colors},
        var_index=var_index,
        modalities=tuple(modalities.values()),
    )
