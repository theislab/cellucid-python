"""
The AnnData stores a Cellucid server hands to the browser, byte for byte.

The web app reads ``.h5ad`` files and Zarr stores itself, with HTTP range
requests: only the bytes a view needs cross the wire, the gene index of a
row-major matrix is built once in the browser, and many genes are read in one
batch. A server therefore never computes a payload. It names the store, says
what it is, and serves its bytes:

- an ``.h5ad`` file, or a Zarr store packed into one ``.zip``, is one file read
  with ``Range`` requests;
- a Zarr store folder is read key by key; a key that does not exist is a 404,
  which a Zarr reader takes as "absent" (an unwritten chunk is its fill value).
  A browser cannot list a folder, so the store must carry consolidated
  metadata; when it does not, the server builds it in memory from the store's
  own metadata files and answers it in place of the missing document.

This module classifies a path by the input format registry the web app
reads too (``input_formats.json``: the same names, signatures and folder
markers, byte for byte the web app's ``tests/fixtures/input-formats.json``),
reads the little a server reports and checks
before it binds (the cell and gene counts and the names of the stored
elements -- metadata only, never an array), finds arrays stored with a codec
the browser's reader does not decode, and builds missing consolidated
metadata. (An AnnData or MuData object, and a file not in priority order, are
served by :mod:`cellucid._served_store`.)
"""

from __future__ import annotations

import json
import os
import re
import zipfile
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from functools import cache
from importlib import resources
from pathlib import Path, PurePosixPath
from typing import TYPE_CHECKING, Any, Literal

if TYPE_CHECKING:
    import anndata

AnnDataLayout = Literal["file", "zip", "folder"]

#: The version of ``input_formats.json`` this module reads.
INPUT_FORMATS_SCHEMA = 1


class RefusedFormatError(ValueError):
    """A file of a format the registry recognises but Cellucid never opens.

    Its message is the registry's ``refusal`` (the route that converts the
    file), the words the web app and cellucid-r use for the same file.
    """


@cache
def input_formats() -> dict[str, Any]:
    """The input format registry: the formats Cellucid reads, by id.

    The packaged ``input_formats.json`` is byte for byte the web app's
    ``tests/fixtures/input-formats.json`` (and cellucid-r's copy): the name
    endings of each format and layout, the container signatures and the
    folder markers. ``cellucid serve`` names a store with the same format id
    (``kind``) the browser's reader uses.
    """
    document: dict[str, Any] = json.loads(
        resources.files("cellucid").joinpath("input_formats.json").read_text("utf-8")
    )
    if document.get("schema") != INPUT_FORMATS_SCHEMA:
        raise RuntimeError(
            f"input_formats.json has schema {document.get('schema')!r}; "
            f"this package reads schema {INPUT_FORMATS_SCHEMA}."
        )
    return document


def _format_of_name(name: str, layouts: tuple[str, ...]) -> tuple[dict[str, Any], str] | None:
    """The format and layout whose longest name ending ``name`` ends with."""
    lower = name.lower()
    best: tuple[dict[str, Any], str, str] | None = None
    for entry in input_formats()["formats"]:
        for layout in layouts:
            for ending in entry["names"].get(layout, []):
                if lower.endswith(ending) and (best is None or len(ending) > len(best[2])):
                    best = (entry, layout, ending)
    return None if best is None else (best[0], best[1])


def _layout_of_head(head: bytes) -> tuple[str, str] | None:
    """The layout and container a file's first bytes carry the signature of."""
    for signature in input_formats()["signatures"]:
        magic = bytes.fromhex(signature["hex"])
        for offset in signature["offsets"]:
            if head[offset : offset + len(magic)] == magic:
                return signature["layout"], signature["container"]
    return None


def _format_of_container(container: str, layout: str) -> dict[str, Any]:
    """The format a container holds when it arrives in ``layout`` and nothing else tells.

    The registry's first such format: an AnnData (``h5ad`` for HDF5, AnnData
    Zarr for Zarr).
    """
    entry: dict[str, Any]
    for entry in input_formats()["formats"]:
        if entry["container"] == container and layout in entry["names"]:
            return entry
    raise ValueError(f"No input format is a {container} container in a {layout} layout.")


def _file_endings() -> str:
    endings = [
        names[0]
        for entry in input_formats()["formats"]
        for layout, names in entry["names"].items()
        # A format known by its content alone (MuData Zarr) has no ending.
        if layout != "folder" and names
    ]
    return endings[0] if len(endings) == 1 else f"{', '.join(endings[:-1])} or {endings[-1]}"

_ZARR_V2_ROOT_KEYS = (".zgroup", ".zattrs", ".zmetadata")
_ZARR_V2_METADATA_FILES = frozenset({".zgroup", ".zattrs", ".zarray"})
_ZARR_V3_METADATA_FILE = "zarr.json"


@dataclass(frozen=True)
class AnnDataStore:
    """One dataset on disk, as the server serves it: a file or a folder.

    Attributes
    ----------
    path
        The resolved file or folder.
    kind
        The input format id (``input_formats.json``): ``"h5ad"``, ``"zarr"``,
        ``"h5mu"``, ``"mudata-zarr"``, ``"loom"``, ``"tenx-h5"``,
        ``"tenx-outs"``, ``"rds"``, ``"rda"``, ``"mtx"``, ...
    layout
        ``"file"`` (one file read as it is), ``"zip"`` (a zipped Zarr store)
        or ``"folder"`` (a Zarr store folder, or a folder of files such as a
        MatrixMarket folder or a Cell Ranger ``outs/``).
    zip_root
        For a zipped Zarr store, the folder inside the archive that holds the
        store (often the store's own name), ``""`` when it is the archive root.
    zarr_format
        For a Zarr store folder, 2 or 3.
    ending
        The name ending the format names the file by (``".rds.gz"``), which
        the served route keeps; ``""`` for a folder.
    files
        For a folder of files (not a Zarr store), the paths of its files
        relative to the folder, sorted: a browser cannot list a folder, so the
        server lists them.
    """

    path: Path
    kind: str
    layout: AnnDataLayout
    zip_root: str = ""
    zarr_format: int | None = None
    ending: str = ""
    files: tuple[str, ...] | None = None


#: The formats whose store is an AnnData (anndata's own layout).
ANNDATA_KINDS = frozenset({"h5ad", "zarr"})
#: The formats whose store is a MuData: a root of shared cells and ``mod/``.
MUDATA_KINDS = frozenset({"h5mu", "mudata-zarr"})


def _container_of(store: AnnDataStore) -> str:
    """The container (``hdf5``, ``zarr``, ``rds``, ``text``) a store's format is read with."""
    return str(_format_entry(store.kind)["container"])


def _format_entry(kind: str) -> dict[str, Any]:
    entry: dict[str, Any]
    for entry in input_formats()["formats"]:
        if entry["id"] == kind:
            return entry
    raise ValueError(f"Unknown input format {kind!r}.")


def _zip_store_root(path: Path) -> str | None:
    """Return the folder inside a ``.zip`` that holds a Zarr store, or None.

    The root is the shallowest folder holding a root metadata document.
    """
    try:
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
    except (OSError, zipfile.BadZipFile) as error:
        raise ValueError(f"{path.name} is not a readable ZIP archive: {error}") from error
    roots: list[PurePosixPath] = []
    for name in names:
        member = PurePosixPath(name)
        if member.name in (*_ZARR_V2_ROOT_KEYS, _ZARR_V3_METADATA_FILE):
            roots.append(member.parent)
    if not roots:
        return None
    shallowest = min(roots, key=lambda parent: (len(parent.parts), str(parent)))
    return "" if str(shallowest) == "." else str(shallowest)


def _root_encoding(attributes: Mapping[str, Any]) -> str | None:
    """A root's ``encoding-type``, lower case (mudata 0.4 writes ``MuData``)."""
    value = attributes.get("encoding-type")
    text = value.decode() if isinstance(value, bytes) else value if isinstance(value, str) else None
    return None if text is None else text.lower()


def _zarr_root_encoding(read: Any) -> str | None:
    """The root group's ``encoding-type`` of a Zarr store, from its metadata.

    ``read(key)`` returns the bytes of one store key, or None.
    """
    for key in (_ZARR_V3_METADATA_FILE, ".zattrs"):
        raw = read(key)
        if raw is None:
            continue
        try:
            document = json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError):
            return None
        attributes = document.get("attributes", {}) if key == _ZARR_V3_METADATA_FILE else document
        return _root_encoding(attributes) if isinstance(attributes, dict) else None
    return None


def _zarr_kind(encoding: str | None) -> str:
    """``mudata-zarr`` for a MuData root, else the AnnData Zarr format."""
    return "mudata-zarr" if encoding == "mudata" else _format_of_container("zarr", "folder")["id"]


def _hdf5_kind(path: Path, named: dict[str, Any]) -> str:
    """An HDF5 file's format: what its root says it is, else what its name says.

    anndata and mudata mark their root (``encoding-type``); Loom, Cell Ranger
    and h5Seurat files do not, so for those the name decides (the browser
    then checks the layout itself).
    """
    import h5py

    try:
        with h5py.File(path, "r") as file:
            encoding = _root_encoding(file.attrs)
    except OSError:
        # Not readable as HDF5 here: the browser, which reads it, says why.
        encoding = None
    if encoding == "anndata":
        return "h5ad"
    if encoding == "mudata":
        return "h5mu"
    return str(named["id"]) if named["container"] == "hdf5" else "h5ad"


def folder_format_of_paths(paths: list[str]) -> str | None:
    """The format of a folder from the paths of its files, or None.

    ``paths`` start with the folder's own name (``outs/analysis/...``). The
    web app's ``folderFormatOfPaths`` twin: a folder named as a format names
    one (``atlas.zarr``) is that format; else the first format, in registry
    order, with a marker file anywhere in it (``zarr.json``), or a file whose
    name matches one of its patterns, compared without case (a MatrixMarket
    matrix: ``matrix.mtx.gz``, ``GSE158055_covid19_counts.mtx.gz``,
    ``GSE119945_gene_count.txt.gz``).
    """
    split = [path.replace("\\", "/").split("/") for path in paths]
    for parts in split:
        for part in parts[:-1]:
            named = _format_of_name(part, ("folder",))
            if named is not None:
                return str(named[0]["id"])
    registry = input_formats()
    markers_of = registry["folderMarkers"]
    patterns_of = registry["folderMarkerPatterns"]
    leaves = {parts[-1] for parts in split}
    for entry in registry["formats"]:
        markers = markers_of.get(entry["id"])
        if markers is not None and any(marker in leaves for marker in markers):
            return str(entry["id"])
        patterns = patterns_of.get(entry["id"])
        if patterns is not None and any(
            re.search(pattern, leaf, re.IGNORECASE) is not None for leaf in leaves for pattern in patterns
        ):
            return str(entry["id"])
    return None


def _folder_files(root: Path) -> list[str]:
    """Every file below ``root``, relative to it, sorted; hidden entries skipped."""
    found: list[str] = []
    for directory, folders, files in os.walk(root):
        folders[:] = sorted(folder for folder in folders if not folder.startswith("."))
        relative = Path(directory).relative_to(root)
        for name in files:
            if name.startswith("."):
                continue
            path = Path(directory) / name
            if path.is_file():
                found.append((relative / name).as_posix())
    return sorted(found)


def _folders_named() -> str:
    labels = [
        f"{'an' if entry['label'][0] in 'AEIOU' else 'a'} {entry['label']}"
        for entry in input_formats()["formats"]
        if "folder" in entry["names"]
    ]
    return ", ".join(labels[:-1]) + f" or {labels[-1]}"


def classify_anndata_path(path: str | Path) -> AnnDataStore:
    """Classify one file or folder Cellucid opens, by the input format registry.

    The rules are the input format registry's (``input_formats.json``), as
    the web app applies them: a file is taken when its name ends as a format
    names it (``.h5ad``, ``.zarr.zip``, ``.h5mu``, ``.loom``, ``.h5``,
    ``.rds[.gz]``, ``.rda``/``.RData[.gz]``, ...), and is then what its first
    bytes say (an HDF5 file named ``.zip`` is an ``.h5ad``; an HDF5 root marked
    ``mudata`` is an ``.h5mu``), its name deciding where they carry no mark. A
    folder is a Zarr store by its root metadata (an AnnData, or a MuData by
    its root's ``encoding-type``), else what its files make it: a MatrixMarket
    folder, a Cell Ranger ``outs/``.

    Raises
    ------
    FileNotFoundError
        When nothing exists at ``path``.
    RefusedFormatError
        When ``path`` is a format Cellucid recognises but never opens (an
        h5Seurat file): the message says how to convert it.
    ValueError
        When ``path`` is none of the formats Cellucid reads.
    """
    resolved = Path(path).expanduser().resolve()
    if not resolved.exists():
        raise FileNotFoundError(f"Path not found: {resolved}")
    if resolved.is_file():
        return _classify_file(resolved)
    if not resolved.is_dir():
        raise ValueError(f"A path Cellucid opens must be a file or a folder: {resolved}")
    if (resolved / _ZARR_V3_METADATA_FILE).is_file() or (resolved / ".zgroup").is_file():
        zarr_format = 3 if (resolved / _ZARR_V3_METADATA_FILE).is_file() else 2

        def read(key: str) -> bytes | None:
            candidate = resolved / key
            return candidate.read_bytes() if candidate.is_file() else None

        return AnnDataStore(
            resolved, _zarr_kind(_zarr_root_encoding(read)), "folder", zarr_format=zarr_format
        )
    files = _folder_files(resolved)
    kind = folder_format_of_paths([f"{resolved.name}/{file}" for file in files])
    if kind is None or _format_entry(kind)["container"] == "zarr":
        raise ValueError(
            f"{resolved} is not a folder Cellucid reads: it is not a Zarr store (no "
            f"zarr.json or .zgroup at its root), and holds no file that makes it "
            f"{_folders_named()} (such as matrix.mtx.gz with barcodes.tsv.gz and "
            "features.tsv.gz, or filtered_feature_bc_matrix.h5)."
        )
    return AnnDataStore(resolved, kind, "folder", files=tuple(files))


def _compressed_copy_of(name: str) -> dict[str, Any] | None:
    """The format a gzip-compressed copy (``x.h5ad.gz``) is a copy of, or None.

    A copy of a file read in parts, as the web app's registry tells it
    (``compressedFormatOfName``): a format read as one stream names its
    compressed endings itself (``.rds.gz``).
    """
    if not name.lower().endswith(".gz"):
        return None
    named = _format_of_name(name[:-3], ("file", "zip"))
    return named[0] if named is not None and named[0].get("access") == "ranges" else None


def _classify_file(resolved: Path) -> AnnDataStore:
    named = _format_of_name(resolved.name, ("file", "zip"))
    compressed = _compressed_copy_of(resolved.name) if named is None else None
    if compressed is not None:
        # The browser reads a served file in parts; a gzip stream has none.
        raise ValueError(
            f"{resolved.name} is a gzip-compressed copy of {compressed['label']}: it cannot be "
            f"read in parts. Decompress it (gunzip {resolved.name}) and serve the "
            "decompressed file. The web app opens the compressed copy itself, from a "
            "link or from this computer, reading it whole once."
        )
    if named is None:
        raise ValueError(
            f"A file Cellucid opens must be an {_file_endings()} file: {resolved}"
        )
    entry, layout = named
    ending = _name_ending(resolved.name, entry, layout)
    width = max(
        max(signature["offsets"]) + len(signature["hex"]) // 2
        for signature in input_formats()["signatures"]
    )
    with resolved.open("rb") as file:
        sniffed = _layout_of_head(file.read(width))
    if sniffed is not None and (sniffed[0] != layout or sniffed[1] != entry["container"]):
        layout = sniffed[0]
        entry = _format_of_container(sniffed[1], layout)
    if layout == "zip":
        root = _zip_store_root(resolved)
        if root is None:
            raise ValueError(
                f"{resolved.name} holds no Zarr store (no .zgroup, .zmetadata or "
                "zarr.json). A zipped AnnData store is the store folder "
                "written by adata.write_zarr, zipped."
            )
        prefix = f"{root}/" if root else ""
        with zipfile.ZipFile(resolved) as archive:
            members = set(archive.namelist())

            def read(key: str) -> bytes | None:
                name = f"{prefix}{key}"
                return archive.read(name) if name in members else None

            kind = _zarr_kind(_zarr_root_encoding(read))
        return AnnDataStore(resolved, kind, "zip", zip_root=root, ending=".zarr.zip")
    kind = _hdf5_kind(resolved, entry) if entry["container"] == "hdf5" else str(entry["id"])
    refusal = _format_entry(kind).get("refusal")
    if refusal is not None:
        # A format the registry recognises but Cellucid never opens (h5Seurat):
        # refused here, in the registry's words, rather than served for the
        # browser to refuse.
        raise RefusedFormatError(f"{resolved.name}: {refusal}")
    if entry["container"] == "hdf5":
        return AnnDataStore(resolved, kind, "file", ending=_route_ending(resolved, kind, ending))
    return AnnDataStore(resolved, kind, "file", ending=ending)


def _name_ending(name: str, entry: dict[str, Any], layout: str) -> str:
    lower = name.lower()
    endings = [ending for ending in entry["names"].get(layout, []) if lower.endswith(ending)]
    return max(endings, key=len) if endings else ""


def _route_ending(path: Path, kind: str, named_ending: str) -> str:
    """The ending the served route carries: the file's own when it names ``kind``."""
    entry = _format_entry(kind)
    own = _name_ending(path.name, entry, "file")
    return own or (entry["names"].get("file") or [named_ending])[0]


# =============================================================================
# What a store holds, from metadata alone
# =============================================================================


@dataclass(frozen=True)
class AnnDataContents:
    """The names and sizes of what one AnnData store holds.

    Read from metadata only: no array is read. ``obsm`` maps each key to its
    width (columns) when it is a 2-D array, and to None otherwise (a table, a
    1-D array). For a MuData store the root describes the shared cells
    (``obs``, ``obsm``, ``obsp``) and the default modality's features, and
    ``modalities`` describes each modality.
    """

    n_obs: int
    n_vars: int
    obs_columns: frozenset[str]
    var_columns: frozenset[str]
    obsm: Mapping[str, int | None]
    obsp: frozenset[str]
    layers: frozenset[str]
    has_x: bool
    raw_var_columns: frozenset[str] | None
    #: The modalities of a multimodal store (MuData), by name, each described
    #: as its own AnnData; None for a store of one feature space.
    modalities: Mapping[str, AnnDataContents] | None = None


def _text(value: Any) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def _text_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str | bytes):
        return [_text(value)]
    return [_text(item) for item in value]


def _is_group(node: Any) -> bool:
    return hasattr(node, "keys") and not hasattr(node, "dtype")


def _child(node: Any, key: str) -> Any | None:
    try:
        return node.get(key)
    except (KeyError, TypeError, AttributeError):
        return None


def _dataframe(node: Any, label: str) -> tuple[int, frozenset[str]]:
    """The row count and column names of one stored dataframe.

    A current dataframe is a group whose ``_index`` attribute names its index
    array; the columns are its ``column-order`` plus any column stored but not
    listed there (as SeuratDisk writes), the same columns the browser's reader
    shows. Before anndata 0.7 a dataframe was one compound dataset whose fields
    are its columns.
    """
    if node is None:
        raise ValueError(f"The AnnData store has no {label} table.")
    if not _is_group(node):
        names = tuple(getattr(node.dtype, "names", None) or ())
        return int(node.shape[0]), frozenset(name for name in names if name != "index")
    index_key = _text(node.attrs.get("_index", "_index"))
    index = _child(node, index_key)
    if index is None or not hasattr(index, "shape"):
        raise ValueError(f"The AnnData store's {label} table has no index array {index_key!r}.")
    columns = set(_text_list(node.attrs.get("column-order")))
    columns.update(key for key in node if key not in (index_key, "__categories"))
    return int(index.shape[0]), frozenset(columns)


def _obsm_widths(node: Any) -> dict[str, int | None]:
    if node is None:
        return {}
    if not _is_group(node):
        # Before anndata 0.7: one compound dataset, a field per array.
        dtype = node.dtype
        return {
            name: (int(dtype[name].shape[0]) if len(dtype[name].shape) == 1 else None)
            for name in dtype.names or ()
        }
    widths: dict[str, int | None] = {}
    for key in node:
        child = node[key]
        shape = child.attrs.get("shape") if _is_group(child) else getattr(child, "shape", None)
        widths[key] = int(shape[1]) if shape is not None and len(shape) == 2 else None
    return widths


def _group_keys(node: Any) -> frozenset[str]:
    if node is None or not _is_group(node):
        return frozenset()
    return frozenset(node.keys())


def _contents_of_root(root: Any) -> AnnDataContents:
    n_obs, obs_columns = _dataframe(_child(root, "obs"), "obs")
    n_vars, var_columns = _dataframe(_child(root, "var"), "var")
    raw = _child(root, "raw")
    raw_var = _child(raw, "var") if raw is not None and _is_group(raw) else _child(root, "raw.var")
    raw_x = _child(raw, "X") if raw is not None and _is_group(raw) else _child(root, "raw.X")
    raw_var_columns = _dataframe(raw_var, "raw.var")[1] if raw_var is not None and raw_x is not None else None
    return AnnDataContents(
        n_obs=n_obs,
        n_vars=n_vars,
        obs_columns=obs_columns,
        var_columns=var_columns,
        obsm=_obsm_widths(_child(root, "obsm")),
        obsp=_group_keys(_child(root, "obsp")),
        layers=_group_keys(_child(root, "layers")),
        has_x=_child(root, "X") is not None,
        raw_var_columns=raw_var_columns,
    )


@contextmanager
def _open_root(store: AnnDataStore) -> Iterator[Any]:
    """The store's root group, open for reading metadata."""
    try:
        if _container_of(store) == "hdf5":
            import h5py

            with h5py.File(store.path, "r") as file:
                yield file
            return
        import zarr

        if store.layout == "zip":
            zip_store = zarr.storage.ZipStore(store.path, mode="r")
            try:
                yield zarr.open_group(zip_store, path=store.zip_root, mode="r")
            finally:
                zip_store.close()
            return
        yield zarr.open_group(str(store.path), mode="r")
    except (OSError, KeyError) as error:
        raise ValueError(f"{store.path.name} cannot be read as AnnData: {error}") from error


def _modality_names(mod: Any) -> list[str]:
    """The modality names in their stored order: ``mod-order`` when it lists them all."""
    keys = [key for key in mod if _is_group(mod[key])]
    order = _text_list(mod.attrs.get("mod-order"))
    return order if sorted(order) == sorted(keys) else keys


def _mudata_contents(root: Any) -> AnnDataContents:
    """A MuData store as the browser presents it (``views/mudata.js``).

    The shared cells' ``obs``, ``obsm`` and ``obsp``, each modality's own
    ``obs`` columns as ``<modality>:<column>`` and ``obsm``/``obsp`` keys as
    ``<modality>:<key>``, the per-modality boolean masks of the root ``obsm``
    as ``<modality>:measured`` columns, and the default modality's features
    (``data_mapping.default_modality``) at the root when it measured every cell.
    """
    from .data_mapping import default_modality

    axis = root.attrs.get("axis", 0)
    if int(axis) != 0:
        raise ValueError(
            f"This MuData shares features, not cells (axis {int(axis)}); Cellucid shows "
            "cells. Write it with shared cells (axis=0)."
        )
    n_obs, obs_columns = _dataframe(_child(root, "obs"), "obs")
    mod = _child(root, "mod")
    if mod is None or not _is_group(mod):
        raise ValueError("The MuData store has no mod group of modalities.")
    names = _modality_names(mod)
    if not names:
        raise ValueError("The MuData store holds no modality.")
    modalities = {name: _contents_of_root(mod[name]) for name in names}
    columns = set(obs_columns)
    obsm: dict[str, int | None] = {}
    for key, width in _obsm_widths(_child(root, "obsm")).items():
        if key in modalities:
            columns.add(f"{key}:measured")
        else:
            obsm[key] = width
    obsp = set(_group_keys(_child(root, "obsp")))
    for name, contents in modalities.items():
        columns.update(
            f"{name}:{column}" for column in contents.obs_columns if column not in obs_columns
        )
        obsm.update({f"{name}:{key}": width for key, width in contents.obsm.items()})
        obsp.update(f"{name}:{key}" for key in contents.obsp)
    chosen = modalities[str(default_modality(names))]
    at_root = chosen.n_obs == n_obs
    return AnnDataContents(
        n_obs=n_obs,
        n_vars=chosen.n_vars,
        obs_columns=frozenset(columns),
        var_columns=chosen.var_columns if at_root else frozenset(),
        obsm=obsm,
        obsp=frozenset(obsp),
        layers=chosen.layers if at_root else frozenset(),
        has_x=chosen.has_x and at_root,
        raw_var_columns=chosen.raw_var_columns if at_root else None,
        modalities=modalities,
    )


def read_anndata_contents(store: AnnDataStore) -> AnnDataContents | None:
    """Read the names and sizes of what ``store`` holds, from metadata only.

    None for a format whose layout only the browser reads (Loom, Cell Ranger
    files and folders, MatrixMarket folders, R objects): its contents are
    learned when the browser opens it.
    """
    if store.kind not in ANNDATA_KINDS | MUDATA_KINDS:
        return None
    with _open_root(store) as root:
        return _mudata_contents(root) if store.kind in MUDATA_KINDS else _contents_of_root(root)


def _scalar_text(node: Any) -> str | None:
    if node is None or _is_group(node) or getattr(node, "shape", None) != ():
        return None
    value = node[()]
    if hasattr(value, "item"):
        value = value.item()
    return value.decode() if isinstance(value, bytes) else value if isinstance(value, str) else None


def read_cellucid_identity(store: AnnDataStore) -> tuple[str | None, str | None]:
    """The dataset id and name a store declares in ``uns['cellucid']``, if any.

    ``uns['cellucid']`` is the dataset's own description (``schema``, ``id``,
    ``name``, ``description``, ``citation``, ``source`` and launch
    ``defaults``), which the viewer reads from the store itself. Only the id
    and name are read here: they name the route and key sessions.
    """
    if store.kind not in ANNDATA_KINDS | MUDATA_KINDS:
        return None, None
    with _open_root(store) as root:
        uns = _child(root, "uns")
        cellucid = _child(uns, "cellucid") if uns is not None and _is_group(uns) else None
        if cellucid is None or not _is_group(cellucid):
            return None, None
        return _scalar_text(_child(cellucid, "id")), _scalar_text(_child(cellucid, "name"))


# =============================================================================
# Consolidated metadata of a Zarr store folder
# =============================================================================


def _is_array_folder(directory: str, files: list[str]) -> bool:
    """Whether ``directory`` is a Zarr array, below which only chunks are stored."""
    if ".zarray" in files:
        return True
    if _ZARR_V3_METADATA_FILE not in files:
        return False
    try:
        document = json.loads(Path(directory, _ZARR_V3_METADATA_FILE).read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return False
    return isinstance(document, dict) and document.get("node_type") == "array"


def _metadata_files(root: Path, names: frozenset[str]) -> Iterator[tuple[str, Path]]:
    """Every metadata file below ``root`` whose name is in ``names``, by store key.

    An array's folder is not walked below: it holds only chunks (as files or
    nested chunk folders), which can number in the millions in a large store.
    """
    for directory, subdirectories, files in os.walk(root):
        subdirectories.sort()
        if _is_array_folder(directory, files):
            subdirectories.clear()
        for name in sorted(files):
            if name in names:
                path = Path(directory, name)
                yield path.relative_to(root).as_posix(), path


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(f"Zarr metadata {path} is not readable UTF-8 JSON: {error}") from error


def consolidated_zarr_metadata(store: AnnDataStore) -> tuple[str, bytes] | None:
    """The consolidated metadata a Zarr store folder lacks, or None when it has it.

    Returns the store key to answer and its bytes: ``.zmetadata`` for a Zarr v2
    store, or the root ``zarr.json`` with ``consolidated_metadata`` added for a
    v3 store. A browser cannot list a folder, so this one document is how it
    learns every array and group the store holds.
    """
    if store.layout != "folder":
        raise ValueError("Only a Zarr store folder is served key by key.")
    root = store.path
    if store.zarr_format == 3:
        document = _read_json(root / _ZARR_V3_METADATA_FILE)
        if not isinstance(document, dict):
            raise ValueError(f"{root / _ZARR_V3_METADATA_FILE} is not a JSON object.")
        if isinstance(document.get("consolidated_metadata"), dict):
            return None
        nodes: dict[str, Any] = {}
        for key, path in _metadata_files(root, frozenset({_ZARR_V3_METADATA_FILE})):
            node = key.rsplit("/", 1)[0] if "/" in key else ""
            if node:
                nodes[node] = _read_json(path)
        consolidated = {
            **document,
            "consolidated_metadata": {
                "kind": "inline",
                "must_understand": False,
                "metadata": nodes,
            },
        }
        return _ZARR_V3_METADATA_FILE, json.dumps(consolidated).encode("utf-8")
    if (root / ".zmetadata").is_file():
        return None
    metadata = {key: _read_json(path) for key, path in _metadata_files(root, _ZARR_V2_METADATA_FILES)}
    document = {"zarr_consolidated_format": 1, "metadata": metadata}
    return ".zmetadata", json.dumps(document).encode("utf-8")


# =============================================================================
# What the browser's reader decodes
# =============================================================================

#: HDF5 filters the browser's reader undoes (``formats/hdf5/filters.js``):
#: deflate, shuffle, Fletcher-32, LZF, and hdf5plugin's Blosc, LZ4,
#: bitshuffle and Zstandard.
BROWSER_HDF5_FILTERS: frozenset[int] = frozenset({1, 2, 3, 32000, 32001, 32004, 32008, 32015})

#: Zarr v2 compressors the browser's reader decodes
#: (``formats/zarr/codecs/pipeline.js``).
BROWSER_ZARR_V2_COMPRESSORS: frozenset[str] = frozenset(
    {"gzip", "zlib", "blosc", "zstd", "lz4", "bz2", "lzma", "crc32c", "shuffle", "bitround"}
)

#: Zarr v2 filters it undoes: the string codecs, and filters that keep sizes.
BROWSER_ZARR_V2_FILTERS: frozenset[str] = frozenset({"vlen-utf8", "vlen-bytes", "shuffle", "bitround"})

#: Zarr v3 codecs it decodes, by name.
BROWSER_ZARR_V3_CODECS: frozenset[str] = frozenset(
    {
        "transpose",
        "bytes",
        "endian",
        "vlen-utf8",
        "vlen-bytes",
        "sharding_indexed",
        "gzip",
        "zstd",
        "blosc",
        "crc32c",
        *(
            f"numcodecs.{name}"
            for name in ("gzip", "zlib", "blosc", "zstd", "lz4", "bz2", "lzma", "crc32c", "shuffle", "bitround")
        ),
    }
)

_HDF5_FILTER_NAMES = {4: "szip", 5: "nbit", 6: "scaleoffset", 307: "bzip2", 32001: "Blosc",
                      32004: "LZ4", 32008: "bitshuffle", 32015: "Zstandard", 32026: "Blosc2"}


def _hdf5_undecodable(path: Path) -> list[str]:
    import h5py

    found: list[str] = []

    def visit(name: str, node: Any) -> None:
        if not isinstance(node, h5py.Dataset):
            return
        plist = node.id.get_create_plist()
        for index in range(plist.get_nfilters()):
            filter_id = int(plist.get_filter(index)[0])
            if filter_id not in BROWSER_HDF5_FILTERS:
                label = _HDF5_FILTER_NAMES.get(filter_id, f"filter {filter_id}")
                found.append(f"{name} ({label})")

    with h5py.File(path, "r") as file:
        file.visititems(visit)
    return found


def _zarr_v3_undecodable(codecs: Any) -> list[str]:
    names: list[str] = []
    for codec in codecs if isinstance(codecs, list) else []:
        name = codec if isinstance(codec, str) else codec.get("name") if isinstance(codec, dict) else None
        if not isinstance(name, str) or name not in BROWSER_ZARR_V3_CODECS:
            names.append(str(name))
        elif name == "sharding_indexed" and isinstance(codec, dict):
            configuration = codec.get("configuration") or {}
            names.extend(_zarr_v3_undecodable(configuration.get("codecs")))
            names.extend(_zarr_v3_undecodable(configuration.get("index_codecs")))
    return names


def _zarr_document_undecodable(key: str, document: Any) -> list[str]:
    if not isinstance(document, dict):
        return []
    if key.endswith(".zarray"):
        names = []
        compressor = document.get("compressor")
        if isinstance(compressor, dict) and compressor.get("id") not in BROWSER_ZARR_V2_COMPRESSORS:
            names.append(str(compressor.get("id")))
        for codec in document.get("filters") or []:
            if isinstance(codec, dict) and codec.get("id") not in BROWSER_ZARR_V2_FILTERS:
                names.append(str(codec.get("id")))
        return names
    if document.get("node_type") == "array":
        return _zarr_v3_undecodable(document.get("codecs"))
    return []


def _zarr_metadata_documents(store: AnnDataStore) -> Iterator[tuple[str, Any]]:
    names = _ZARR_V2_METADATA_FILES | {_ZARR_V3_METADATA_FILE}
    if store.layout == "folder":
        for key, path in _metadata_files(store.path, frozenset(names)):
            yield key, _read_json(path)
        return
    prefix = f"{store.zip_root}/" if store.zip_root else ""
    with zipfile.ZipFile(store.path) as archive:
        for name in archive.namelist():
            if name.startswith(prefix) and PurePosixPath(name).name in names:
                yield name[len(prefix) :], json.loads(archive.read(name).decode("utf-8"))


def browser_undecodable(store: AnnDataStore) -> list[str]:
    """Every stored array the browser's reader cannot decode, with what stops it.

    Read from metadata only (HDF5 filter pipelines, Zarr compressors, filters
    and codecs). An empty list means the browser reads the whole store.
    """
    container = _container_of(store)
    if store.files is not None:
        # A folder of files: its HDF5 files (a Cell Ranger outs/ matrix) are
        # read in place; its text files are read once and need no codec.
        return [
            f"{file}: {problem}"
            for file in store.files
            if file.lower().endswith(".h5")
            for problem in _hdf5_undecodable(store.path / file)
        ]
    if container == "hdf5":
        return _hdf5_undecodable(store.path)
    if container != "zarr":
        return []
    found = []
    for key, document in _zarr_metadata_documents(store):
        names = _zarr_document_undecodable(key, document)
        if names:
            found.append(f"{key.rsplit('/', 1)[0] if '/' in key else '/'} ({', '.join(names)})")
    return found


def read_into_memory(store: AnnDataStore) -> anndata.AnnData:
    """Read a whole store into memory with anndata (and its codec plugins)."""
    import anndata as ad

    if _container_of(store) == "hdf5":
        return ad.read_h5ad(store.path)
    if store.layout == "zip":
        import zarr

        zip_store = zarr.storage.ZipStore(store.path, mode="r")
        try:
            return ad.read_zarr(zarr.open_group(zip_store, path=store.zip_root, mode="r"))
        finally:
            zip_store.close()
    return ad.read_zarr(store.path)


# =============================================================================
# An in-memory MuData object, as a temporary file
# =============================================================================


def is_mudata(value: object) -> bool:
    """Whether ``value`` is a ``mudata.MuData``, without importing mudata."""
    return any(
        cls.__name__ == "MuData" and cls.__module__.split(".")[0] == "mudata"
        for cls in type(value).__mro__
    )
