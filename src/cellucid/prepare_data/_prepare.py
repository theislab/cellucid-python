"""
:func:`prepare`: an AnnData, a MuData, or arrays and data frames in; one Cellucid-ready AnnData or MuData out.
"""

from __future__ import annotations

import functools
import os
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any, Literal, ParamSpec, TypeVar

import numpy as np
import pandas as pd
from scipy import sparse

from .._anndata_files import AnnDataContents, is_mudata
from .._contracts import (
    _require_dataset_id,
    _require_dataset_name,
    _require_display_text,
    _require_native_boolean,
)
from .._launch_defaults import launch_defaults
from .._plural import plural
from ..data_mapping import EmbeddingKeySpec, column_gene_names, default_modality
from ..subsample import _require_uint32
from ._data import (
    anndata_components,
    cellucid_record,
    is_anndata,
    is_cell_array,
    split_obsm,
)
from ._inputs import (
    linked_text,
    resolve_connectivities,
    resolve_embeddings,
    resolve_expression,
    resolve_genes,
    resolve_obs,
    resolve_vector_fields,
    stored_expression,
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
from ._priority import CellPriority, in_original_order, prefix_levels, without_positions
from ._store import DEFAULT_ZARR_FORMAT, Elements, OutputTarget, write_store

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
    var: pd.DataFrame | None = None,
    gene_expression: np.ndarray | sparse.spmatrix | sparse.sparray | None = None,
    connectivities: np.ndarray | sparse.spmatrix | sparse.sparray | None = None,
    vector_fields: Mapping[str, Mapping[str, Any]] | None = None,
    modalities: Sequence[str] | None = None,
    gene_identifiers: Sequence[str] | None = None,
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

        - AnnData: ``obs``; ``X`` as ``gene_expression`` with ``var`` (neither
          when ``X`` is absent); every numeric ``obsm`` array as an embedding
          under its own key, except a ``<field>_<basis>`` array as wide as
          ``X_<basis>``, which is that vector field; ``obsp['connectivities']``;
          and from ``uns['cellucid']`` the identity (``id``, ``name``,
          ``description``, ``citation``, ``source``) and the launch defaults.
          Layers, ``raw`` and other ``obsp`` matrices are not written. To leave
          out a component ``data`` holds, remove it from a copy first.
        - MuData (modalities over shared cells, ``axis`` 0): the cells'
          ``mdata.obs``, ``mdata.obsm`` (the shared embeddings and vector
          fields), ``mdata.obsp['connectivities']`` and ``mdata.uns['cellucid']``
          at the root, as above; each modality whole -- its ``X``, ``var``,
          ``obs``, numeric ``obsm`` arrays and ``obsp['connectivities']`` --
          as ``mod/<name>``. Columns and arrays mudata copies to the root
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
        column; sparse input is densified. ``X_umap_1d``/``X_umap_2d``/``X_umap_3d``,
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
        ``(n_cells, n_cells)`` neighbour graph, stored as
        ``obsp['connectivities']`` (CSR float32). Exactly symmetric, finite,
        non-negative, with a zero diagonal.
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
        cell's distance to its category's centre, as the share of its
        category at that distance or closer) are computed here and stored in
        ``uns['cellucid']['outlierQuantiles']``, so the viewer filters
        outliers at once instead of reading the whole latent space to
        compute them.
    obs_keys
        The obs columns the viewer lists, in this order; every column is
        stored either way.
    var_gene_id_column
        The var column whose values name the genes in the viewer. Omitted,
        the viewer chooses: the var index, or, when the index holds database
        identifiers (such as Ensembl ids) and a var column holds gene symbols
        (such as CELLxGENE's ``feature_name``), that column, the identifiers
        staying searchable beside it.
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
        "var": var,
        "gene_expression": gene_expression,
        "connectivities": connectivities,
        "vector_fields": vector_fields,
        "gene_identifiers": gene_identifiers,
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
    if data is not None and not is_mudata(data) and not is_anndata(data):
        raise TypeError(
            f"data must be an anndata.AnnData or a mudata.MuData, got {type(data).__name__}. "
            "Pass arrays and data frames as embeddings=, obs=, gene_expression= and the rest."
        )
    if data is not None:
        # The cells in their original order: a file prepared from an object
        # read in a stored priority order names every cell as that file did.
        data = in_original_order(data)
    selected = _selected_modalities(data, modalities)
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
        arguments = {**_mudata_root_components(data, every), **given, "priority_seed": priority_seed}
        _require_identity(arguments, data_given=True)
        write_store(target, _mudata_layout(data, names, arguments))
    else:
        if selected is not None:
            if "modality" in given:
                raise ValueError(
                    "modality names the modality the viewer opens on, and modalities= writes "
                    "one, which is stored as an AnnData with no modalities."
                )
            components = _single_modality_components(data, selected[0])
        else:
            components = anndata_components(data) if data is not None else {}
        arguments = {**components, **given, "priority_seed": priority_seed}
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
        write_store(target, _anndata_elements(arguments))
    return target.path


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
        description = _require_display_text(description, label="dataset_description", allow_empty=True)
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
        gene_id_column=arguments.get("var_gene_id_column"),
        normalize_embeddings=arguments.get("normalize_embeddings"),
        centroid_min_points=arguments.get("centroid_min_points"),
        obs_keys=arguments.get("obs_keys"),
        embedding_keys=embedding_defaults,
        knn_key=None,
        expression_source=None,
        vector_field_default=arguments.get("vector_field_default"),
        embedding_label="embedding_defaults",
        gene_id_label="var_gene_id_column",
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
    priority: dict[str, Any] | None,
) -> dict[str, Any]:
    """``cellucid`` with its launch defaults, for a latent space its outlier quantiles, and its priority."""
    if defaults:
        cellucid["defaults"] = _uns_value({key: defaults[key] for key in _DEFAULTS_ORDER if key in defaults})
    if latent is not None:
        cellucid["outlierQuantiles"] = stored_outlier_quantiles(
            obs_columns,
            defaults["latentKey"],
            latent,
            defaults.get("centroidMinPoints", DEFAULT_MIN_POINTS),
        )
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


def _anndata_elements(arguments: Mapping[str, Any]) -> Elements:
    """One Cellucid-ready AnnData's elements, every component checked."""
    cellucid = _cellucid_identity(arguments)
    obsm, n_cells = resolve_embeddings(arguments["embeddings"])
    vectors, fields = resolve_vector_fields(arguments.get("vector_fields"), obsm, n_cells)
    obs_index, obs_columns = resolve_obs(arguments.get("obs"), n_cells)
    obs_columns = without_positions(obs_columns, "obs")

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
        x = None
        var_index = np.asarray([], dtype=object)
        var_columns: dict[str, Any] = {}
    else:
        expression, var_index, var_columns = resolve_genes(
            resolve_expression(gene_expression, n_cells),
            arguments.get("var"),
            var_gene_id_column,
            arguments.get("gene_identifiers"),
        )
        x = stored_expression(expression)
        # The matrix before its values took their stored type is not kept.
        del expression

    connectivities = arguments.get("connectivities")
    obsp = {} if connectivities is None else {
        "connectivities": resolve_connectivities(connectivities, n_cells)
    }

    # Every per-cell element in priority order, before anything is derived
    # from it (the outlier quantiles), so the file holds what its rows give.
    priority = CellPriority.of(n_cells, arguments.get("priority_seed"))
    levels: np.ndarray | None = None
    if priority is not None:
        obsm = priority.rows_of_mapping(obsm)
        vectors = priority.rows_of_mapping(vectors)
        obs_index = priority.rows(obs_index)
        obs_columns = priority.obs_columns(obs_columns)
        obsp = {key: priority.graph(graph) for key, graph in obsp.items()}
        if x is not None:
            x = priority.matrix(x)
        levels = prefix_levels(n_cells)

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
        priority=None if priority is None else priority.record(levels),
    )
    return Elements(
        x=x,
        obs_index=obs_index,
        obs=obs_columns,
        var_index=var_index,
        var=var_columns,
        obsm=all_obsm,
        obsp=obsp,
        uns={"title": cellucid["name"], "cellucid": cellucid},
        prefix_levels=levels,
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


def _single_modality_components(mdata: Any, name: str) -> dict[str, Any]:
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

    components = _mudata_root_components(mdata, every)
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
    components["gene_expression"] = placed(modality.X)
    components["var"] = modality.var
    if "connectivities" not in components and "connectivities" in modality.obsp:
        graph = modality.obsp["connectivities"]
        components["connectivities"] = graph if in_order else graph[rows][:, rows]
    return components


def _copied_from_modality(key: str, names: Sequence[str], keys_of: Callable[[str], Any]) -> bool:
    """Whether ``key`` is mudata's root copy ``<modality>:<key>`` of a modality's own entry."""
    name, separator, inner = key.partition(":")
    return bool(separator) and name in names and inner in keys_of(name)


def _mudata_root_components(mdata: Any, names: Sequence[str]) -> dict[str, Any]:
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
    components: dict[str, Any] = {"obs": obs.drop(columns=copied), "embeddings": embeddings}
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


def _selected_var_axis(
    mdata: Any, modalities: Mapping[str, Modality]
) -> tuple[np.ndarray, dict[str, Modality]]:
    """The root ``var`` and ``varmap`` of the modalities written (``modalities=``).

    mudata's root var is every modality's features in modality order, and a
    modality's varmap places them there (1-based, 0 for another modality's);
    for a selection both are made from the modalities written alone.
    """
    from dataclasses import replace

    names = [np.asarray([str(value) for value in mdata.mod[name].var_names], dtype=object) for name in modalities]
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


def _mudata_layout(mdata: Any, names: Sequence[str], arguments: Mapping[str, Any]) -> MuDataLayout:
    """One Cellucid-ready MuData, every component but the modalities' matrices checked."""
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
    obs_columns = without_positions(obs_columns, "obs")
    connectivities = arguments.get("connectivities")
    obsp = {} if connectivities is None else {
        "connectivities": resolve_connectivities(connectivities, n_cells)
    }
    masks = {
        name: np.asarray(mdata.obsm[name], dtype=bool)
        for name in names
        if name in mdata.obsm and np.asarray(mdata.obsm[name]).dtype == bool
    }
    modalities = {name: resolve_modality(mdata, name) for name in names}
    if list(names) == [str(name) for name in mdata.mod]:
        var_index = global_var_index(mdata)
    else:
        var_index, modalities = _selected_var_axis(mdata, modalities)
    # The cells in priority order, each modality's rows by its cells' global
    # priority, so a prefix of the cells is a prefix of every modality.
    priority = CellPriority.of(n_cells, arguments.get("priority_seed"))
    levels: np.ndarray | None = None
    if priority is not None:
        obsm = priority.rows_of_mapping(obsm)
        vectors = priority.rows_of_mapping(vectors)
        obs_index = priority.rows(obs_index)
        obs_columns = priority.obs_columns(obs_columns)
        obsp = {key: priority.graph(graph) for key, graph in obsp.items()}
        masks = {name: priority.rows(mask) for name, mask in masks.items()}
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
            name, _, key = latent_key.partition(":")
            modality = modalities[name]
            if np.any(modality.obsmap == 0):
                raise ValueError(
                    f"latent_key names {latent_key!r}, and modality {name!r} did not measure "
                    "every cell, so the outlier filter could not measure the others. Name an "
                    "embedding every cell has: one of mdata.obsm, or of a modality that "
                    "measured every cell."
                )
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
        priority=None if priority is None else priority.record(levels),
    )
    return MuDataLayout(
        obs_index=obs_index,
        obs=obs_columns,
        obsm=all_obsm,
        masks=masks,
        obsp=obsp,
        uns={"title": cellucid["name"], "cellucid": cellucid},
        var_index=var_index,
        modalities=tuple(modalities.values()),
    )
