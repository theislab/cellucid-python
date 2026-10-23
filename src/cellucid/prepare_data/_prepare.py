"""
:func:`prepare`: arrays and data frames in, one Cellucid-ready AnnData out.
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

from .._anndata_files import AnnDataContents
from .._contracts import (
    _require_dataset_id,
    _require_dataset_name,
    _require_display_text,
    _require_native_boolean,
)
from .._launch_defaults import launch_defaults
from ..data_mapping import EmbeddingKeySpec
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
from ._outliers import DEFAULT_MIN_POINTS, stored_outlier_quantiles
from ._store import DEFAULT_ZARR_FORMAT, Elements, OutputTarget, write_store

#: The order ``uns['cellucid']['defaults']`` lists its keys in.
_DEFAULTS_ORDER = (
    "embedding1d",
    "embedding2d",
    "embedding3d",
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
}

#: Which ``format`` to choose; every refusal of ``format`` quotes it.
_FORMAT_CHOICES = (
    'format must be "h5ad", "zarr" or "zarr.zip":\n'
    '  "h5ad"      one .h5ad file: for single-file hosts (Zenodo, Figshare, Hugging Face, '
    "S3, a lab server), local files, cellucid serve and Jupyter.\n"
    '  "zarr"      a .zarr folder of small files: for GitHub (no file over 100 MB, no LFS).\n'
    '  "zarr.zip"  the Zarr store in one .zarr.zip file: wherever an .h5ad goes.'
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
    for suffix in (".zarr.zip", ".zarr", ".h5ad"):
        if lowered.endswith(suffix):
            return suffix
    return None


def _resolve_target(
    out_path: object, *, output_format: object, zarr_format: object, force: bool
) -> OutputTarget:
    """What ``format`` writes, checked against ``out_path``, and whether it may go there."""
    if not (isinstance(output_format, str) and output_format in _FORMATS):
        raise ValueError(f"{_FORMAT_CHOICES}\nGot {output_format!r}.")
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
    kind: Literal["h5ad", "zarr"] = "h5ad" if output_format == "h5ad" else "zarr"
    zipped = output_format == "zarr.zip"
    resolved_format: Literal[2, 3] | None
    if kind == "h5ad":
        if zarr_format is not None:
            raise ValueError(
                'zarr_format applies to format="zarr" or "zarr.zip", not format="h5ad".'
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
    format: Literal["h5ad", "zarr", "zarr.zip"],
    embeddings: Mapping[str, Any],
    dataset_id: str,
    dataset_name: str,
    obs: pd.DataFrame | None = None,
    var: pd.DataFrame | None = None,
    gene_expression: np.ndarray | sparse.spmatrix | sparse.sparray | None = None,
    connectivities: np.ndarray | sparse.spmatrix | sparse.sparray | None = None,
    vector_fields: Mapping[str, Mapping[str, Any]] | None = None,
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
    zarr_format: Literal[2, 3] | None = None,
    force: bool = False,
) -> Path:
    """Write a Cellucid-ready AnnData: a standard AnnData laid out for the browser.

    The output is an ordinary AnnData that ``anndata.read_h5ad`` or
    ``anndata.read_zarr`` reads, written so the Cellucid web app reads it fast
    over HTTP: the expression matrix gene-major (CSC), arrays chunked for range
    requests and compressed with codecs browsers decode natively, and the
    dataset's name, citation and launch defaults in ``uns['cellucid']``. Open
    it with ``cellucid serve``, :func:`~cellucid.show_anndata`, or the web app
    directly. The written contract is the "Cellucid-ready AnnData" page of the
    documentation; cellucid-r's ``cellucid_prepare()`` writes the same.

    Every input is checked before anything is written, and the output is
    written beside ``out_path`` under a temporary name and moved into place
    only once it is complete.

    Parameters
    ----------
    out_path
        Where to write, ending in the suffix of ``format``: ``.h5ad``,
        ``.zarr`` or ``.zarr.zip``. A leading ``~`` is expanded.
    format
        Required, with no default -- which one depends on where the data goes:

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
    embeddings
        ``{obsm key: (n_cells, k) array}``, at least one. Each is stored whole
        as ``obsm[key]`` in float32, every column, at its own coordinates (the
        viewer normalizes them). A latent space for outlier measurement is one
        of them, named by ``latent_key``. A plain ``(n_cells,)`` vector is one
        column; sparse input is densified. ``X_umap`` (1 to 3 columns) or
        ``X_umap_1d``/``X_umap_2d``/``X_umap_3d`` are shown by default.
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
        The var column whose values name the genes in the viewer.
    vector_field_default
        The ``vector_fields`` field the overlay opens on; one of its
        embeddings must be shown when the file opens (a chosen or default
        embedding, or a 1-3 column array that is not the latent space).
    normalize_embeddings
        False shows embeddings at their own coordinates.
    centroid_min_points
        The fewest cells a category needs to get outlier quantiles (10 when
        None).
    zarr_format
        2 or 3, for ``format="zarr"`` or ``"zarr.zip"``; Zarr v2 when None.
        Refused with ``format="h5ad"``.
    force
        Replace an existing output of the same kind. A folder that is not a
        Zarr store is never replaced.

    The launch defaults (``embedding_defaults`` to ``centroid_min_points``)
    are written to ``uns['cellucid']['defaults']`` only when given; a choice
    left at None stays the viewer's own.

    Returns
    -------
    Path
        The written file or folder.
    """
    force = _require_native_boolean(force, label="force")
    target = _resolve_target(
        out_path, output_format=format, zarr_format=zarr_format, force=force
    )
    dataset_id = _require_dataset_id(dataset_id, label="dataset_id")
    dataset_name = _require_dataset_name(dataset_name, label="dataset_name")
    description = (
        None
        if dataset_description is None
        else _require_display_text(dataset_description, label="dataset_description", allow_empty=True)
    )
    citation = linked_text(
        citation_text,
        citation_url,
        text_label="citation_text",
        url_label="citation_url",
        text_key="text",
    )
    source = linked_text(
        source_name,
        source_url,
        text_label="source_name",
        url_label="source_url",
        text_key="name",
    )

    obsm, n_cells = resolve_embeddings(embeddings)
    vectors, fields = resolve_vector_fields(vector_fields, obsm, n_cells)
    obs_index, obs_columns = resolve_obs(obs, n_cells)

    if gene_expression is None:
        supplied = [
            name
            for name, value in (
                ("var", var),
                ("gene_identifiers", gene_identifiers),
                ("var_gene_id_column", var_gene_id_column),
            )
            if value is not None
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
            var,
            var_gene_id_column,
            gene_identifiers,
        )
        x = stored_expression(expression)

    obsp = {} if connectivities is None else {
        "connectivities": resolve_connectivities(connectivities, n_cells)
    }

    if vector_field_default is not None and vector_field_default not in fields:
        if not fields:
            raise ValueError("vector_field_default names a field, and no vector_fields were supplied.")
        raise ValueError(
            f"vector_field_default must name one of vector_fields: {fields}; "
            f"got {vector_field_default!r}."
        )
    all_obsm = {**obsm, **vectors}
    if latent_key is not None and latent_key not in obsm:
        raise ValueError(f"latent_key must name one of the embeddings: {list(obsm)}; got {latent_key!r}.")
    defaults = launch_defaults(
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
        latent_key=latent_key,
        gene_id_column=var_gene_id_column,
        normalize_embeddings=normalize_embeddings,
        centroid_min_points=centroid_min_points,
        obs_keys=obs_keys,
        embedding_keys=embedding_defaults,
        knn_key=None,
        expression_source=None,
        vector_field_default=vector_field_default,
        embedding_label="embedding_defaults",
        gene_id_label="var_gene_id_column",
    )
    # An embedding default named by its key alone is stored as the key: the
    # viewer then plots its first columns, as it does for any key.
    for slot, spec in (embedding_defaults or {}).items():
        if isinstance(spec, str):
            defaults[f"embedding{slot}"] = spec

    cellucid: dict[str, Any] = {"schema": 1, "id": dataset_id, "name": dataset_name}
    if description is not None:
        cellucid["description"] = description
    if citation is not None:
        cellucid["citation"] = citation
    if source is not None:
        cellucid["source"] = source
    if defaults:
        cellucid["defaults"] = _uns_value({key: defaults[key] for key in _DEFAULTS_ORDER if key in defaults})
    if latent_key is not None:
        cellucid["outlierQuantiles"] = stored_outlier_quantiles(
            obs_columns,
            latent_key,
            obsm[latent_key],
            defaults.get("centroidMinPoints", DEFAULT_MIN_POINTS),
        )

    write_store(
        target,
        Elements(
            x=x,
            obs_index=obs_index,
            obs=obs_columns,
            var_index=var_index,
            var=var_columns,
            obsm=all_obsm,
            obsp=obsp,
            uns={"title": dataset_name, "cellucid": cellucid},
        ),
    )
    return target.path
