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
the browser's reader does not decode, builds missing consolidated metadata,
and writes an in-memory AnnData object -- or a store the browser cannot
decode, read with anndata and its codec plugins -- to a temporary ``.h5ad``
so it is served exactly like a file.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
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
    """The format a container holds when it arrives in ``layout``."""
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
        if layout != "folder"
    ]
    return " or ".join(endings)

_ZARR_V2_ROOT_KEYS = (".zgroup", ".zattrs", ".zmetadata")
_ZARR_V2_METADATA_FILES = frozenset({".zgroup", ".zattrs", ".zarray"})
_ZARR_V3_METADATA_FILE = "zarr.json"


@dataclass(frozen=True)
class AnnDataStore:
    """One AnnData store on disk, as the server serves it.

    Attributes
    ----------
    path
        The resolved file or folder.
    kind
        The input format id (``input_formats.json``): ``"h5ad"`` or
        ``"zarr"``.
    layout
        ``"file"`` (an ``.h5ad``), ``"zip"`` (a zipped Zarr store) or
        ``"folder"`` (a Zarr store folder).
    zip_root
        For a zipped Zarr store, the folder inside the archive that holds the
        store (often the store's own name), ``""`` when it is the archive root.
    zarr_format
        For a Zarr store folder, 2 or 3.
    """

    path: Path
    kind: str
    layout: AnnDataLayout
    zip_root: str = ""
    zarr_format: int | None = None


def _container_of(store: AnnDataStore) -> str:
    """The container (``hdf5`` or ``zarr``) a store's format is read with."""
    for entry in input_formats()["formats"]:
        if entry["id"] == store.kind:
            return str(entry["container"])
    raise ValueError(f"Unknown input format {store.kind!r}.")


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


def classify_anndata_path(path: str | Path) -> AnnDataStore:
    """Classify one ``.h5ad`` file, zipped Zarr store, or Zarr store folder.

    The rules are the input format registry's (``input_formats.json``), as
    the web app applies them: a file is taken when its name ends as a format
    names it, and is then what its first bytes say (an HDF5 file named
    ``.zip`` is an ``.h5ad``), its name deciding only when they carry no
    signature; a folder is a Zarr store by its root metadata.

    Raises
    ------
    FileNotFoundError
        When nothing exists at ``path``.
    ValueError
        When ``path`` is none of the three.
    """
    resolved = Path(path).expanduser().resolve()
    if not resolved.exists():
        raise FileNotFoundError(f"Path not found: {resolved}")
    if resolved.is_file():
        named = _format_of_name(resolved.name, ("file", "zip"))
        if named is None:
            raise ValueError(
                f"An AnnData file must be an {_file_endings()} file: {resolved}"
            )
        width = max(
            max(signature["offsets"]) + len(signature["hex"]) // 2
            for signature in input_formats()["signatures"]
        )
        with resolved.open("rb") as file:
            sniffed = _layout_of_head(file.read(width))
        entry, layout = named
        if sniffed is not None:
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
            return AnnDataStore(resolved, entry["id"], "zip", zip_root=root)
        return AnnDataStore(resolved, entry["id"], "file")
    if not resolved.is_dir():
        raise ValueError(
            f"AnnData path must be an .h5ad file, a zipped Zarr store or a Zarr store "
            f"folder: {resolved}"
        )
    folder_format = _format_of_container("zarr", "folder")["id"]
    if (resolved / _ZARR_V3_METADATA_FILE).is_file():
        return AnnDataStore(resolved, folder_format, "folder", zarr_format=3)
    if (resolved / ".zgroup").is_file():
        return AnnDataStore(resolved, folder_format, "folder", zarr_format=2)
    raise ValueError(
        f"{resolved} is not a Zarr store: it has no zarr.json (Zarr v3) or .zgroup "
        "(Zarr v2) at its root."
    )


# =============================================================================
# What a store holds, from metadata alone
# =============================================================================


@dataclass(frozen=True)
class AnnDataContents:
    """The names and sizes of what one AnnData store holds.

    Read from metadata only: no array is read. ``obsm`` maps each key to its
    width (columns) when it is a 2-D array, and to None otherwise (a table, a
    1-D array).
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


def read_anndata_contents(store: AnnDataStore) -> AnnDataContents:
    """Read the names and sizes of what ``store`` holds, from metadata only."""
    with _open_root(store) as root:
        return _contents_of_root(root)


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
#: deflate, shuffle, Fletcher-32 and LZF.
BROWSER_HDF5_FILTERS: frozenset[int] = frozenset({1, 2, 3, 32000})

#: Zarr v2 compressors the browser's reader decodes
#: (``formats/zarr/codecs/pipeline.js``).
BROWSER_ZARR_V2_COMPRESSORS: frozenset[str] = frozenset(
    {"gzip", "zlib", "blosc", "zstd", "lz4", "crc32c", "shuffle", "bitround"}
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
            for name in ("gzip", "zlib", "blosc", "zstd", "lz4", "crc32c", "shuffle", "bitround")
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
    if _container_of(store) == "hdf5":
        return _hdf5_undecodable(store.path)
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
# An in-memory AnnData object, as a temporary file
# =============================================================================


def _column_major(matrix: Any) -> Any:
    """A sparse matrix in CSC form, so one gene is one contiguous range.

    The browser reads a gene from a CSC matrix with one range request per
    column; from a CSR matrix it first builds a gene index over every stored
    entry. Anything else is written as it is.
    """
    from scipy import sparse

    if sparse.issparse(matrix) and matrix.format != "csc":
        return sparse.csc_matrix(matrix)
    return matrix


def write_temporary_h5ad(adata: anndata.AnnData, *, dataset_id: str) -> Path:
    """Write ``adata`` to a new temporary ``.h5ad`` and return its path.

    The object is not modified. Sparse expression matrices (``X``, every layer,
    ``raw.X``) are written column-major, uncompressed, so the browser reads one
    gene as one contiguous range. The file lives in its own new temporary
    folder; the caller removes that folder with :func:`remove_temporary_h5ad`.
    """
    import anndata as ad

    if adata.n_obs <= 0:
        raise ValueError("AnnData datasets must contain at least one observation")
    if adata.isbacked:
        # The object as it is now, including what was changed since it was opened.
        adata = adata.to_memory()
    folder = Path(tempfile.mkdtemp(prefix="cellucid-anndata-"))
    path = folder / f"{dataset_id}.h5ad"
    try:
        snapshot = ad.AnnData(
            X=None if adata.X is None else _column_major(adata.X),
            obs=adata.obs,
            var=adata.var,
            uns=dict(adata.uns),
            obsm=dict(adata.obsm),
            varm=dict(adata.varm),
            obsp=dict(adata.obsp),
            varp=dict(adata.varp),
            layers={key: _column_major(value) for key, value in adata.layers.items()},
        )
        if adata.raw is not None:
            snapshot.raw = ad.AnnData(
                X=_column_major(adata.raw.X),
                var=adata.raw.var,
                varm=dict(adata.raw.varm),
            )
        snapshot.write_h5ad(path)
    except BaseException:
        shutil.rmtree(folder, ignore_errors=True)
        raise
    return path


def remove_temporary_h5ad(path: Path) -> None:
    """Remove a file :func:`write_temporary_h5ad` wrote, with its folder."""
    shutil.rmtree(path.parent, ignore_errors=True)
