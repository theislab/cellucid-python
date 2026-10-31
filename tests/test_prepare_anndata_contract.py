"""The Cellucid-ready AnnData ``prepare()`` writes, read back exactly.

``prepare()`` writes a standard AnnData that anndata reads without a warning,
laid out for the web app's range reads: X gene-major, chunked arrays, codecs
browsers decode natively, and ``uns['cellucid']``. The written contract is
``docs/user_guide/python_package/c_data_preparation_api/09_cellucid_ready_anndata_specification.md``;
these tests pin each of its rules on all four outputs (``.h5ad``, Zarr v2 and
v3 folders, ``.zarr.zip``).
"""

from __future__ import annotations

import inspect
import io
import json
import re
import warnings
import zipfile
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
import zarr
from scipy import sparse

from cellucid import prepare, restore_original_order
from cellucid._anndata_files import (
    browser_undecodable,
    classify_anndata_path,
    read_cellucid_identity,
)
from cellucid.anndata_server import AnnDataServer
from cellucid.data_mapping import parse_vector_field_key
from cellucid.prepare_data import _store
from cellucid.prepare_data._outliers import latent_outlier_quantiles

N_CELLS = 40
N_GENES = 12

OUTPUTS = [
    pytest.param("data.h5ad", "h5ad", None, id="h5ad"),
    pytest.param("data.zarr", "zarr", 2, id="zarr-v2"),
    pytest.param("data.zarr", "zarr", 3, id="zarr-v3"),
    pytest.param("data.zarr.zip", "zarr.zip", None, id="zarr-zip"),
]


def _inputs() -> dict[str, Any]:
    rng = np.random.default_rng(7)
    counts = sparse.random(
        N_CELLS, N_GENES, density=0.3, format="csr", random_state=3, dtype=np.float64
    )
    counts.data = np.floor(counts.data * 1000) + 1
    graph = sparse.random(N_CELLS, N_CELLS, density=0.1, format="csr", random_state=5)
    graph = graph + graph.T
    graph.setdiag(0)
    graph.eliminate_zeros()
    obs = pd.DataFrame(
        {
            "cell_type": rng.choice(["T", "B", "NK"], N_CELLS),
            "barcode": [f"AAAC{i:04d}" for i in range(N_CELLS)],
            "donor": pd.Categorical(rng.choice(["d2", "d1"], N_CELLS), categories=["d2", "d1"]),
            "score": rng.normal(size=N_CELLS),
            "n_genes": rng.integers(0, 500, N_CELLS),
            "is_doublet": rng.random(N_CELLS) > 0.8,
            "batch": pd.array([1, None] * (N_CELLS // 2), dtype="Int64"),
        },
        index=[f"cell{i}" for i in range(N_CELLS)],
    )
    obs.loc["cell3", "score"] = np.nan
    var = pd.DataFrame(
        {"symbol": [f"GENE{i}" for i in range(N_GENES)], "highly_variable": [True, False] * 6},
        index=[f"ENSG{i:05d}" for i in range(N_GENES)],
    )
    return {
        "embeddings": {
            "X_umap": rng.normal(size=(N_CELLS, 2)),
            "X_umap_3d": rng.normal(size=(N_CELLS, 3)),
            "X_pca": rng.normal(size=(N_CELLS, 6)),
        },
        "dataset_id": "toy-v1",
        "dataset_name": "Toy blood",
        "obs": obs,
        "var": var,
        "gene_expression": counts,
        "connectivities": graph,
        "vector_fields": {
            "velocity": {
                "X_umap": rng.normal(size=(N_CELLS, 2)),
                "X_umap_3d": rng.normal(size=(N_CELLS, 3)),
                "X_pca": rng.normal(size=(N_CELLS, 6)),
            }
        },
    }


def _write(
    tmp_path: Path, name: str, output_format: str, zarr_format: int | None, **overrides: Any
) -> Path:
    arguments = {**_inputs(), **overrides}
    if zarr_format is not None:
        arguments["zarr_format"] = zarr_format
    return prepare(tmp_path / name, format=output_format, **arguments)


def _read(path: Path) -> ad.AnnData:
    """The written AnnData in the order its cells were given (``restore_original_order``)."""
    return restore_original_order(_stored(path))


def _stored(path: Path) -> ad.AnnData:
    """The written AnnData as stored: its cells in priority order."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        if path.suffix == ".h5ad":
            return ad.read_h5ad(path)
        if path.name.endswith(".zip"):
            store = zarr.storage.ZipStore(path, mode="r")
            try:
                return ad.read_zarr(store)
            finally:
                store.close()
        return ad.read_zarr(path)


# =============================================================================
# Round trip
# =============================================================================


@pytest.mark.parametrize(("name", "output_format", "zarr_format"), OUTPUTS)
def test_anndata_reads_back_every_value_exactly(
    tmp_path: Path, name: str, output_format: str, zarr_format: int | None
) -> None:
    given = _inputs()
    adata = _read(_write(tmp_path, name, output_format, zarr_format))

    assert adata.shape == (N_CELLS, N_GENES)
    assert list(adata.obs_names) == list(given["obs"].index)
    assert list(adata.var_names) == list(given["var"].index)

    # X: CSC, float32, and every integer count exactly as given.
    assert sparse.issparse(adata.X) and adata.X.format == "csc"
    assert adata.X.dtype == np.float32
    assert (given["gene_expression"] != adata.X).nnz == 0

    # obs: repeated text is a categorical in code-point order, unique text a
    # string column, the rest as given.
    obs = adata.obs
    assert isinstance(obs["cell_type"].dtype, pd.CategoricalDtype)
    assert list(obs["cell_type"].cat.categories) == ["B", "NK", "T"]
    assert list(obs["cell_type"]) == list(given["obs"]["cell_type"])
    assert obs["barcode"].dtype == object and list(obs["barcode"]) == list(given["obs"]["barcode"])
    assert list(obs["donor"].cat.categories) == ["d2", "d1"]
    np.testing.assert_array_equal(obs["score"].to_numpy(), given["obs"]["score"].to_numpy())
    np.testing.assert_array_equal(obs["n_genes"].to_numpy(), given["obs"]["n_genes"].to_numpy())
    assert obs["is_doublet"].dtype == bool
    assert obs["batch"].dtype == "Int64" and obs["batch"].isna().sum() == N_CELLS // 2
    assert list(adata.var["symbol"]) == list(given["var"]["symbol"])

    # obsm: float32, every column, at the caller's own coordinates.
    for key, values in given["embeddings"].items():
        assert adata.obsm[key].dtype == np.float32
        np.testing.assert_array_equal(adata.obsm[key], values.astype(np.float32))
    # Each vector field is stored under scVelo's name for it, `<field>_<basis>`
    # for the embedding `X_<basis>`, as wide as that embedding.
    for embedding, values in given["vector_fields"]["velocity"].items():
        key = "velocity_" + embedding.removeprefix("X_")
        assert adata.obsm[key].dtype == np.float32
        np.testing.assert_array_equal(adata.obsm[key], values.astype(np.float32))
        parsed = parse_vector_field_key(key, set(adata.obsm.keys()))
        assert parsed is not None and (parsed.field, parsed.embedding_key) == ("velocity", embedding)
    assert sorted(adata.obsm.keys()) == [
        "X_pca", "X_umap", "X_umap_3d", "velocity_pca", "velocity_umap", "velocity_umap_3d",
    ]

    # obsp: symmetric CSR float32.
    graph = adata.obsp["connectivities"]
    assert graph.format == "csr" and graph.dtype == np.float32
    assert (graph != graph.T).nnz == 0
    np.testing.assert_array_equal(graph.toarray(), given["connectivities"].toarray().astype(np.float32))

    assert adata.uns["title"] == "Toy blood"
    assert adata.uns["cellucid"] == {"schema": 1, "id": "toy-v1", "name": "Toy blood"}


def test_a_dense_expression_matrix_stays_dense_and_gene_major(tmp_path: Path) -> None:
    dense = _inputs()["gene_expression"].toarray()
    path = _write(tmp_path, "data.h5ad", "h5ad", None, gene_expression=dense)
    adata = _read(path)
    assert isinstance(adata.X, np.ndarray) and adata.X.dtype == np.float32
    np.testing.assert_array_equal(adata.X, dense.astype(np.float32))
    with h5py.File(path) as file:
        # A tall strip of genes: one gene is one column of whole chunks.
        assert file["X"].chunks == (N_CELLS, N_GENES)
        # Expression values are not shuffled.
        assert not file["X"].shuffle


def test_counts_float32_cannot_hold_exactly_are_stored_as_float64(tmp_path: Path) -> None:
    counts = _inputs()["gene_expression"].copy()
    counts.data[0] = 2**24 + 1
    adata = _read(_write(tmp_path, "data.h5ad", "h5ad", None, gene_expression=counts))
    assert adata.X.dtype == np.float64
    assert (counts != adata.X).nnz == 0


def test_the_type_of_x_is_decided_on_the_genes_kept(tmp_path: Path) -> None:
    counts = _inputs()["gene_expression"].copy().tocsc()
    counts.data[counts.indptr[0]] = 2**24 + 1  # a count of the first gene alone
    kept = [f"ENSG{i:05d}" for i in range(1, N_GENES)]
    adata = _read(_write(tmp_path, "a.h5ad", "h5ad", None, gene_expression=counts, gene_identifiers=kept))
    assert adata.X.dtype == np.float32
    assert (counts[:, 1:] != adata.X).nnz == 0
    adata = _read(_write(tmp_path, "b.h5ad", "h5ad", None, gene_expression=counts.toarray(), gene_identifiers=kept))
    assert adata.X.dtype == np.float32
    np.testing.assert_array_equal(adata.X, counts[:, 1:].toarray())
    adata = _read(_write(tmp_path, "c.h5ad", "h5ad", None, gene_expression=counts, gene_identifiers=kept[::-1] + ["ENSG00000"]))
    assert adata.X.dtype == np.float64


def test_counts_float32_holds_exactly_stay_float32_and_fractions_are_rounded_once(
    tmp_path: Path,
) -> None:
    counts = _inputs()["gene_expression"].copy().astype(np.int64)
    counts.data[0] = 2**24
    adata = _read(_write(tmp_path, "a.h5ad", "h5ad", None, gene_expression=counts))
    assert adata.X.dtype == np.float32 and (counts != adata.X).nnz == 0

    normalized = _inputs()["gene_expression"].copy()
    normalized.data = np.log1p(normalized.data) / 3.0
    adata = _read(_write(tmp_path, "b.h5ad", "h5ad", None, gene_expression=normalized))
    assert adata.X.dtype == np.float32
    np.testing.assert_array_equal(adata.X.toarray(), normalized.toarray().astype(np.float32))


def test_stored_zeros_are_dropped_and_values_kept(tmp_path: Path) -> None:
    counts = _inputs()["gene_expression"].copy()
    counts.data[:3] = 0
    adata = _read(_write(tmp_path, "data.h5ad", "h5ad", None, gene_expression=counts))
    assert adata.X.nnz == counts.nnz - 3
    np.testing.assert_array_equal(adata.X.toarray(), counts.toarray())


@pytest.mark.parametrize("layout", ["csr", "csc"])
def test_a_backed_anndata_is_prepared_from_disk(tmp_path: Path, layout: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """``adata.X`` of a backed AnnData is anndata's on-disk matrix, laid out by blocks
    (:mod:`cellucid.prepare_data._gene_major`), never read whole."""
    from cellucid.prepare_data import _gene_major

    monkeypatch.setattr(_gene_major, "PASS_ENTRIES", 50)  # several blocks over 40 cells
    inputs = _inputs()
    counts = inputs["gene_expression"].asformat(layout)
    source = tmp_path / "source.h5ad"
    ad.AnnData(X=counts, obs=pd.DataFrame(index=[f"c{i}" for i in range(N_CELLS)]),
               var=pd.DataFrame(index=[f"g{i}" for i in range(N_GENES)]),
               obsm={"X_umap": inputs["embeddings"]["X_umap"]}).write_h5ad(source)
    backed = ad.read_h5ad(source, backed="r")
    path = prepare(tmp_path / "out.h5ad", format="h5ad", data=backed, dataset_id="b", dataset_name="B")
    backed.file.close()
    written = _read(path)
    np.testing.assert_array_equal(written.X.toarray(), counts.toarray())


def test_a_lazily_read_sparse_x_is_prepared_block_by_block(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A dask array of sparse blocks (anndata's lazily read X) is never computed
    whole: it is decoded block of rows by block of rows into the layout's
    scratch folder, and written byte for byte as the same matrix in memory."""
    da = pytest.importorskip("dask.array")
    inputs = _inputs()
    counts = inputs["gene_expression"].tocsr()
    lazy = da.from_array(counts, chunks=(7, N_GENES), asarray=False)
    computed: list[tuple[int, ...]] = []
    real_compute = type(lazy).compute

    def compute(self: Any, **kwargs: Any) -> Any:
        computed.append(tuple(int(extent) for extent in self.shape))
        return real_compute(self, **kwargs)

    monkeypatch.setattr(type(lazy), "compute", compute)
    expected = _write(tmp_path / "memory", "data.zarr", "zarr", 2, gene_expression=counts)
    actual = _write(tmp_path / "lazy", "data.zarr", "zarr", 2, gene_expression=lazy)
    assert computed and all(shape[0] <= 7 for shape in computed)
    files = sorted(str(path.relative_to(expected)) for path in expected.rglob("*") if path.is_file())
    assert files == sorted(str(path.relative_to(actual)) for path in actual.rglob("*") if path.is_file())
    assert [name for name in files if (expected / name).read_bytes() != (actual / name).read_bytes()] == []


@pytest.mark.parametrize("kind", ["csc_matrix", "csc_array"])
def test_a_callers_non_canonical_csc_matrix_is_never_changed(tmp_path: Path, kind: str) -> None:
    """Stored zeros, duplicates and unsorted indices are canonicalised on a copy.

    ``csc_matrix()`` of a CSC input shares its arrays, and scipy canonicalises
    in place: the caller's ``indptr``/``indices``/``data`` were rewritten.
    """
    dense = _inputs()["gene_expression"].toarray()
    csc = sparse.csc_matrix(dense)
    # Column 0: indices reversed, its first value split into a duplicate pair,
    # and a stored zero on a row the column has no value in.
    start, stop = csc.indptr[0], csc.indptr[1]
    present = csc.indices[start:stop]
    absent = int(np.flatnonzero(dense[:, 0] == 0)[0])
    rows = np.concatenate([present[::-1], [present[0], absent]])
    values = np.concatenate([csc.data[start:stop][::-1], [0.0, 0.0]])
    values[stop - start - 1] /= 2  # present[0], now last of the reversed run
    values[stop - start] = values[stop - start - 1]
    indices = np.concatenate([rows, csc.indices[stop:]]).astype(np.int32)
    data = np.concatenate([values, csc.data[stop:]])
    indptr = csc.indptr.copy()
    indptr[1:] += len(rows) - (stop - start)
    constructor = getattr(sparse, kind)
    given = constructor((data, indices, indptr), shape=csc.shape)
    snapshot = (given.data.copy(), given.indices.copy(), given.indptr.copy())
    expected = given.toarray()
    written = _read(_write(tmp_path, "data.h5ad", "h5ad", None, gene_expression=given))
    np.testing.assert_array_equal(written.X.toarray(), expected)
    for before, after in zip(snapshot, (given.data, given.indices, given.indptr), strict=True):
        np.testing.assert_array_equal(after, before)


# =============================================================================
# Layout: chunks, codecs, consolidated metadata
# =============================================================================


def test_h5ad_chunks_and_filters(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # Small chunk targets so the toy arrays span several chunks.
    monkeypatch.setattr(_store, "SLICED_CHUNK_ELEMENTS", 16)
    monkeypatch.setattr(_store, "WHOLE_CHUNK_ELEMENTS", 32)
    monkeypatch.setattr(_store, "CELL_CHUNK_ELEMENTS", 24)
    path = _write(tmp_path, "data.h5ad", "h5ad", None)
    with h5py.File(path) as file:
        assert dict(file.attrs) == {"encoding-type": "anndata", "encoding-version": "0.1.0"}
        assert file["X"].attrs["encoding-type"] == "csc_matrix"
        for key in ("X/data", "X/indices"):
            assert file[key].chunks == (16,)
            assert file[key].compression == "gzip"
            assert file[key].compression_opts == _store.DEFLATE_LEVEL
        # Every array wider than a byte is shuffled, the expression values not.
        assert not file["X/data"].shuffle
        assert file["X/indices"].shuffle and file["X/indptr"].shuffle
        assert file["obsm/X_umap"].shuffle and file["obsp/connectivities/data"].shuffle
        assert not file["obs/is_doublet"].shuffle
        assert file["X/indptr"].chunks == (N_GENES + 1,)
        assert file["X/indices"].dtype == np.int32
        # Every array with one row per cell in chunks of CELL_CHUNK_ELEMENTS
        # values; an embedding chunk holds every column of a block of rows.
        assert file["obs/score"].chunks == (24,)
        assert file["obs/_index"].chunks == (24,)
        assert file["obs/cellucid_original_position"].chunks == (24,)
        assert file["obsm/X_pca"].chunks == (24 // 6, 6)
        assert file["obsm/X_umap"].chunks == (12, 2)
        # The rest (var, the graph) in chunks of WHOLE_CHUNK_ELEMENTS.
        assert file["obsp/connectivities/data"].chunks == (32,)
        # File format 1.10 (superblock 3), and every object's header in the
        # metadata block at the start of the file, before any array's values.
        assert file.id.get_create_plist().get_version()[0] == 3
        headers: list[int] = []
        file.visit(lambda name: headers.append(h5py.h5o.get_info(file[name].id).addr))
        assert max(headers) < _store.H5AD_META_BLOCK_BYTES
        assert min(file["X/data"].id.get_chunk_info(0).byte_offset, file["obs/score"].id.get_chunk_info(0).byte_offset) > max(headers)


class _RecordedFile(io.RawIOBase):
    """A file that records where every read the HDF5 library makes ends."""

    def __init__(self, path: Path) -> None:
        self._file = path.open("rb")
        self.ends: list[int] = []

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def seek(self, offset: int, whence: int = 0) -> int:
        return self._file.seek(offset, whence)

    def tell(self) -> int:
        return self._file.tell()

    def readinto(self, buffer: Any) -> int:
        start = self._file.tell()
        count = self._file.readinto(buffer)
        self.ends.append(start + count)
        return count

    def close(self) -> None:
        self._file.close()
        super().close()


def test_h5ad_opening_reads_lie_before_every_other_value(tmp_path: Path) -> None:
    """What opening reads -- the tree, every attribute's string, ``uns`` but
    its outlier quantiles, the gene names, ``X``'s pointers and every
    categorical obs column's labels -- comes first, so the web app reads it
    with its first request of the file."""
    path = _write(tmp_path, "data.h5ad", "h5ad", None, latent_key="X_pca")
    recorded = _RecordedFile(path)
    try:
        with h5py.File(recorded, "r") as file:
            file.visititems(lambda name, item: [item.attrs[key] for key in item.attrs])
            dict(file.attrs)
            file["var/_index"][...]
            file["X/indptr"][...]
            # Every value of uns but the outlier quantiles: its scalars, and
            # the priority's levels.
            file["uns"].visititems(
                lambda name, item: item[()]
                if isinstance(item, h5py.Dataset) and "outlierQuantiles" not in name
                else None
            )
            labels = [
                f"obs/{key}/categories" for key in file["obs"] if file["obs"][key].attrs["encoding-type"] == "categorical"
            ]
            assert labels, "the toy obs has categorical columns"
            for path in labels:
                file[path][...]
            opening_end = max(recorded.ends)
            bulk: list[int] = []

            def first_chunk(name: str, item: Any) -> None:
                opened = (
                    name in ("var/_index", "X/indptr")
                    or name.endswith("/categories")
                    or (name.startswith("uns/") and "outlierQuantiles" not in name)
                )
                if isinstance(item, h5py.Dataset) and item.chunks is not None and item.id.get_num_chunks() and not opened:
                    bulk.append(item.id.get_chunk_info(0).byte_offset)

            file.visititems(first_chunk)
            assert "outlierQuantiles" in file["uns/cellucid"]
    finally:
        recorded.close()
    assert opening_end <= min(bulk)


def test_zarr_v2_codec_and_consolidated_metadata(tmp_path: Path) -> None:
    path = _write(tmp_path, "data.zarr", "zarr", 2)
    metadata = json.loads((path / ".zmetadata").read_text())["metadata"]
    data = metadata["X/data/.zarray"]
    assert data["compressor"] == {"id": "zlib", "level": _store.DEFLATE_LEVEL}
    # Multi-byte values are byte-shuffled before they are deflated; the
    # expression values, strings and one-byte values are not.
    assert data["filters"] is None
    assert metadata["X/indices/.zarray"]["filters"] == [{"id": "shuffle", "elementsize": 4}]
    assert metadata["X/indptr/.zarray"]["filters"] == [{"id": "shuffle", "elementsize": 4}]
    assert metadata["obsp/connectivities/data/.zarray"]["filters"] == [{"id": "shuffle", "elementsize": 4}]
    assert metadata["obs/score/.zarray"]["filters"] == [{"id": "shuffle", "elementsize": 8}]
    assert metadata["obs/batch/values/.zarray"]["filters"] == [{"id": "shuffle", "elementsize": 8}]
    assert metadata["obs/batch/mask/.zarray"]["filters"] is None
    assert metadata["obs/is_doublet/.zarray"]["filters"] is None
    assert metadata["obs/_index/.zarray"]["filters"] == [{"id": "vlen-utf8"}]
    assert data["chunks"] == [min(_store.SLICED_CHUNK_ELEMENTS, data["shape"][0])]
    assert data["dtype"] == "<f4"
    assert metadata["obsm/X_pca/.zarray"]["chunks"] == [N_CELLS, 6]
    assert metadata["X/.zattrs"]["encoding-type"] == "csc_matrix"
    assert metadata["uns/cellucid/.zattrs"]["encoding-type"] == "dict"


def _embedding_defaults_and_metadata() -> dict[str, Any]:
    """Arguments that put a text, a number, a flag and a list in ``uns``."""
    return {
        "citation_text": "Doe et al.",
        "embedding_defaults": {"3d": ("X_pca", [2, 0, 1])},
        "obs_keys": ["cell_type", "score"],
        "normalize_embeddings": False,
        "centroid_min_points": 3,
    }


def test_every_zarr_v2_array_carries_the_store_codec(tmp_path: Path) -> None:
    # uns scalars and an empty array (the var index of a dataset without
    # genes) included: anndata would give those zarr's default codec.
    full = _write(tmp_path, "full.zarr", "zarr", 2, **_embedding_defaults_and_metadata())
    bare = prepare(
        tmp_path / "bare.zarr",
        format="zarr",
        embeddings={"X_umap": np.zeros((3, 2))},
        dataset_id="bare",
        dataset_name="Bare",
    )
    for path in (full, bare):
        metadata = json.loads((path / ".zmetadata").read_text())["metadata"]
        arrays = {key: value for key, value in metadata.items() if key.endswith(".zarray")}
        assert "var/_index/.zarray" in arrays
        for key, array in arrays.items():
            assert array["compressor"] == {"id": "zlib", "level": _store.DEFLATE_LEVEL}, key
            assert array["fill_value"] in (0, 0.0, False, ""), key
    metadata = json.loads((full / ".zmetadata").read_text())["metadata"]
    assert metadata["uns/cellucid/defaults/obsKeys/.zattrs"]["encoding-type"] == "string-array"
    assert metadata["uns/cellucid/schema/.zattrs"]["encoding-type"] == "numeric-scalar"
    assert metadata["uns/cellucid/schema/.zarray"]["dtype"] == "<i8"


#: The codecs the Zarr v3 specification defines, and ``vlen-utf8``, the
#: codec registered for its variable-length strings
#: (github.com/zarr-developers/zarr-extensions). ``numcodecs.*`` codecs are
#: outside both: another Zarr implementation need not read them.
ZARR_V3_SPECIFICATION_CODECS = frozenset(
    {"bytes", "transpose", "gzip", "blosc", "crc32c", "sharding_indexed", "vlen-utf8"}
)


def _blosc(itemsize: int) -> dict[str, Any]:
    configuration = {"typesize": itemsize, "cname": "zlib", "clevel": _store.DEFLATE_LEVEL}
    return {"name": "blosc", "configuration": {**configuration, "shuffle": "shuffle", "blocksize": 0}}


def test_zarr_v3_codec_and_consolidated_metadata(tmp_path: Path) -> None:
    path = _write(tmp_path, "data.zarr", "zarr", 3, **_embedding_defaults_and_metadata())
    root = json.loads((path / "zarr.json").read_text())
    nodes = root["consolidated_metadata"]["metadata"]
    gzip = {"name": "gzip", "configuration": {"level": _store.DEFLATE_LEVEL}}
    # The expression values, strings and one-byte values: gzip alone.
    assert nodes["X/data"]["data_type"] == "float32"
    assert nodes["X/data"]["codecs"][1:] == [gzip]
    assert nodes["obs/_index"]["codecs"][1:] == [gzip]
    # Every other multi-byte array: Blosc's byte shuffle and zlib.
    assert nodes["X/indices"]["codecs"][1:] == [_blosc(4)]
    assert nodes["uns/cellucid/schema"]["codecs"][1:] == [_blosc(8)]
    for key, node in nodes.items():
        if node["node_type"] == "array":
            assert node["codecs"][1:] in ([gzip], [_blosc(np.dtype(_V3_DTYPES[node["data_type"]]).itemsize)]), key


#: The NumPy dtype of each Zarr v3 data type the writers use.
_V3_DTYPES = {
    "bool": "|b1",
    "int8": "|i1",
    "int16": "<i2",
    "int32": "<i4",
    "int64": "<i8",
    "uint32": "<u4",
    "float32": "<f4",
    "float64": "<f8",
    "string": "|O",
}


@pytest.mark.parametrize("output_format", ["zarr", "zarr.zip"])
def test_zarr_v3_uses_only_codecs_of_the_specification(tmp_path: Path, output_format: str) -> None:
    """Checked on the written metadata, whatever the installed zarr-python
    warns about: a codec outside the specification is an interoperability
    defect even where no warning is raised (zarr-python 3.2 and later raise
    none for ``numcodecs.*``)."""
    name = "data.zarr" if output_format == "zarr" else "data.zarr.zip"
    full = _write(tmp_path, name, output_format, 3, **_embedding_defaults_and_metadata())
    bare = prepare(
        tmp_path / ("bare." + name.split(".", 1)[1]),
        format=output_format,
        zarr_format=3,
        embeddings={"X_umap": np.zeros((3, 2))},
        dataset_id="bare",
        dataset_name="Bare",
    )
    for path in (full, bare):
        if output_format == "zarr":
            documents = [json.loads(entry.read_text()) for entry in path.rglob("zarr.json")]
        else:
            with zipfile.ZipFile(path) as archive:
                documents = [
                    json.loads(archive.read(entry)) for entry in archive.namelist() if entry.endswith("zarr.json")
                ]
        root = next(document for document in documents if "consolidated_metadata" in document)
        documents += list(root["consolidated_metadata"]["metadata"].values())
        arrays = [document for document in documents if document["node_type"] == "array"]
        assert arrays
        for array in arrays:
            names = {codec["name"] for codec in array["codecs"]}
            assert names <= ZARR_V3_SPECIFICATION_CODECS, names


def _assert_blosc_frame(frame: bytes, raw: bytes, typesize: int) -> None:
    """``frame`` is the Blosc 1 frame of ``raw`` that c-blosc writes for the
    writers' ``blosc`` codec (cname zlib, clevel 6, byte shuffle, automatic
    block size) compressing in one thread, and that cellucid-r's
    ``.blosc_frame()`` writes: every header field, block, stream and the
    choice between a zlib stream and a stored one. The zlib streams are
    checked by what they inflate to, not byte for byte, since deflate's output
    differs between zlib builds (a CPython links its own or the system's;
    numcodecs bundles one)."""
    import struct
    import zlib

    def splits(block: int) -> bool:
        return typesize <= 16 and block // typesize >= 128

    nbytes = len(raw)
    block = nbytes
    if nbytes >= 32 * 1024:
        block = 32 * 1024 * 2 * 8  # c-blosc's L1, doubled for zlib, x8 at clevel 6
    if splits(block):
        block = min(max(min(block, 1 << 18) * typesize, 1 << 16), 1 << 20)
    block = min(block, nbytes)
    if block > typesize:
        block = block // typesize * typesize
    split = splits(block)
    flags = 0x01 | (0 if split else 0x10) | (3 << 5)  # byte shuffle, zlib
    version, codec_version, stored_flags, stored_typesize, stored_nbytes, stored_block, total = struct.unpack(
        "<BBBBIII", frame[:16]
    )
    assert (version, codec_version, stored_typesize, stored_nbytes, stored_block, total) == (
        2, 1, typesize, nbytes, block, len(frame)
    )
    if nbytes < 128 or stored_flags == flags | 0x02:
        # Copied whole: too small to compress, or nothing gained.
        assert stored_flags == flags | 0x02 and frame[16:] == raw
        return
    assert stored_flags == flags
    n_blocks = -(-nbytes // block)
    starts = struct.unpack(f"<{n_blocks}I", frame[16 : 16 + 4 * n_blocks])
    at = 16 + 4 * n_blocks
    for index in range(n_blocks):
        assert starts[index] == at
        values = raw[index * block : (index + 1) * block]
        shuffled = np.frombuffer(values, np.uint8).reshape(-1, typesize).T.tobytes()
        streams = typesize if split and len(values) == block else 1
        each = len(values) // streams
        for stream in range(streams):
            (length,) = struct.unpack("<i", frame[at : at + 4])
            data = frame[at + 4 : at + 4 + length]
            at += 4 + length
            plane = shuffled[stream * each : (stream + 1) * each]
            assert data == plane if length == each else zlib.decompress(data) == plane
            assert 0 < length <= each
    assert at == len(frame)


@pytest.mark.parametrize("dtype", ["<i2", "<i4", "<f4", "<i8", "<f8"])
@pytest.mark.parametrize("values", ["ramp", "noise", "sparse"])
@pytest.mark.parametrize("length", [5, 31, 100, 4_000, 20_000, 300_000])
def test_a_blosc_chunk_is_the_documented_frame(tmp_path: Path, dtype: str, values: str, length: int) -> None:
    """Every case of the frame: a whole copy (under 128 bytes, or nothing
    gained), one stream (under 128 values), one plane per byte, planes stored
    as they are, and a float64 chunk of 2^18 values, two blocks of 1 MiB."""
    rng = np.random.default_rng(length)
    if values == "ramp":
        array = (np.arange(length) % 1009).astype(dtype)
    elif values == "noise":
        array = rng.integers(-(2**14), 2**14, length).astype(dtype) if dtype[1] == "i" else rng.standard_normal(length).astype(dtype)
    else:
        array = np.where(rng.random(length) < 0.05, rng.integers(1, 100, length), 0).astype(dtype)
    chunk = min(length, _store.WHOLE_CHUNK_ELEMENTS)
    stored = zarr.create_array(
        str(tmp_path / "a.zarr"),
        shape=array.shape,
        chunks=(chunk,),
        dtype=dtype,
        filters=None,
        compressors=_store._blosc_codec(array.dtype.itemsize),
        zarr_format=3,
    )
    stored[...] = array
    padded = np.zeros(-(-length // chunk) * chunk, dtype=dtype)
    padded[:length] = array
    for index in range(len(padded) // chunk):
        raw = padded[index * chunk : (index + 1) * chunk]
        path = tmp_path / "a.zarr" / "c" / str(index)
        if not raw.any():
            # zarr writes no chunk that holds only the fill value.
            assert not path.exists()
            continue
        _assert_blosc_frame(path.read_bytes(), raw.tobytes(), array.dtype.itemsize)


def test_a_zipped_store_is_the_store_folder_at_the_archive_root_entries_stored(
    tmp_path: Path,
) -> None:
    path = _write(tmp_path, "data.zarr.zip", "zarr.zip", None)
    with zipfile.ZipFile(path) as archive:
        entries = archive.infolist()
    names = [entry.filename for entry in entries]
    assert ".zgroup" in names and ".zmetadata" in names and "X/data/0" in names
    assert all(entry.compress_type == zipfile.ZIP_STORED for entry in entries)
    assert classify_anndata_path(path).zip_root == ""
    # Every chunk first, then the metadata documents, the consolidated
    # metadata last: the end of the archive, read when it opens, describes
    # the store.
    documents = [name.rsplit("/", 1)[-1] in _store.ZARR_METADATA_NAMES for name in names]
    assert documents == sorted(documents)
    assert names[-1] == ".zmetadata"


def test_a_zipped_store_holds_each_arrays_chunks_in_grid_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(_store, "SLICED_CHUNK_ELEMENTS", 2)
    for zarr_format, prefix in ((2, "X/data/"), (3, "X/data/c/")):
        path = _write(tmp_path, f"v{zarr_format}.zarr.zip", "zarr.zip", zarr_format)
        with zipfile.ZipFile(path) as archive:
            chunks = [
                name[len(prefix) :]
                for name in archive.namelist()
                if name.startswith(prefix) and name[len(prefix) :].isdigit()
            ]
        assert len(chunks) > 10
        assert chunks == [str(position) for position in range(len(chunks))]



@pytest.mark.parametrize("zarr_format", [2, 3])
def test_a_zipped_store_depends_on_its_entries_alone(tmp_path: Path, zarr_format: int) -> None:
    """No entry carries the time or permissions the store was written with,
    so the same input always gives the same bytes (as cellucid-r writes it)."""
    first = _write(tmp_path / "first", "data.zarr.zip", "zarr.zip", zarr_format)
    with zipfile.ZipFile(first) as archive:
        entries = archive.infolist()
    assert {entry.date_time for entry in entries} == {(1980, 1, 1, 0, 0, 0)}
    assert {entry.external_attr >> 16 for entry in entries} == {0o100644}
    second = _write(tmp_path / "second", "data.zarr.zip", "zarr.zip", zarr_format)
    assert first.read_bytes() == second.read_bytes()


def test_a_zarr_v3_chunk_carries_no_time(tmp_path: Path) -> None:
    """gzip's header has a time field; numcodecs stamps the time of writing
    in it, which would make every run's chunks differ. Every gzip chunk,
    strings included, holds a zero time; every other chunk is the Blosc
    frame of its values (which carries no time)."""
    path = _write(tmp_path, "data.zarr", "zarr", 3)
    root = zarr.open_group(path, mode="r")
    chunks = [entry for entry in path.rglob("*") if entry.is_file() and "/c/" in entry.as_posix()]
    assert chunks
    framed = 0
    for chunk in chunks:
        data = chunk.read_bytes()
        key = chunk.relative_to(path).as_posix().split("/c/")[0]
        if json.loads((path / key / "zarr.json").read_text())["codecs"][-1]["name"] == "blosc":
            array = root[key]
            coords = [int(part) for part in chunk.relative_to(path / key / "c").parts]
            raw = np.zeros(array.chunks, dtype=array.dtype)
            box = tuple(slice(at * size, (at + 1) * size) for at, size in zip(coords, array.chunks, strict=True))
            values = np.asarray(array[box])
            raw[tuple(slice(0, extent) for extent in values.shape)] = values
            _assert_blosc_frame(data, raw.tobytes(), array.dtype.itemsize)
            framed += 1
            continue
        assert data[:2] == b"\x1f\x8b", chunk
        assert data[4:8] == bytes(4), chunk
    assert framed


@pytest.mark.parametrize(("name", "output_format", "zarr_format"), OUTPUTS)
def test_the_browser_decodes_every_array_and_the_server_names_the_dataset(
    tmp_path: Path, name: str, output_format: str, zarr_format: int | None
) -> None:
    path = _write(
        tmp_path,
        name,
        output_format,
        zarr_format,
        embedding_defaults={"2d": "X_umap"},
        obs_keys=["cell_type", "score"],
    )
    store = classify_anndata_path(path)
    assert browser_undecodable(store) == []
    assert read_cellucid_identity(store) == ("toy-v1", "Toy blood")

    # Served as it is: no name or id given, no launch choice repeated (the
    # viewer reads uns['cellucid']['defaults'] from the store itself).
    server = AnnDataServer(path, quiet=True, serve_web_ui=False)
    assert (server.dataset_id, server.dataset_name) == ("toy-v1", "Toy blood")
    assert server.defaults == {}
    entry = server._datasets_document()["datasets"][0]
    assert entry["id"] == "toy-v1" and entry["anndata"]["defaults"] == {}
    if store.layout == "folder":
        # Consolidated metadata is in the store, so nothing is built in memory.
        assert server._served.generated_metadata is None


# =============================================================================
# uns['cellucid']
# =============================================================================


def test_uns_cellucid_holds_exactly_what_was_given(tmp_path: Path) -> None:
    adata = _read(
        _write(
            tmp_path,
            "data.h5ad",
            "h5ad",
            None,
            dataset_description="",
            citation_text="Doe et al. 2026",
            citation_url="https://doi.org/10.0000/x",
            source_name="CELLxGENE",
            embedding_defaults={"2d": "X_umap", "3d": ("X_pca", [2, 0, 1])},
            latent_key="X_pca",
            obs_keys=["cell_type", "score"],
            var_gene_id_column="symbol",
            vector_field_default="velocity",
            normalize_embeddings=False,
            centroid_min_points=5,
        )
    )
    cellucid = adata.uns["cellucid"]
    defaults = cellucid.pop("defaults")
    outliers = cellucid.pop("outlierQuantiles")
    assert outliers["latentKey"] == "X_pca"
    assert outliers["minPoints"] == 5
    assert cellucid == {
        "schema": 1,
        "id": "toy-v1",
        "name": "Toy blood",
        "description": "",
        "citation": {"text": "Doe et al. 2026", "url": "https://doi.org/10.0000/x"},
        "source": {"name": "CELLxGENE"},
    }
    assert defaults["embedding2d"] == "X_umap"
    assert defaults["embedding3d"]["key"] == "X_pca"
    assert list(defaults["embedding3d"]["columns"]) == [2, 0, 1]
    assert defaults["latentKey"] == "X_pca"
    assert list(defaults["obsKeys"]) == ["cell_type", "score"]
    assert defaults["geneIdColumn"] == "symbol"
    assert defaults["vectorFieldDefault"] == "velocity"
    assert isinstance(defaults["normalizeEmbeddings"], bool | np.bool_)
    assert not defaults["normalizeEmbeddings"]
    assert defaults["centroidMinPoints"] == 5
    # Every key is one the web app's launch-default validator accepts.
    assert set(defaults) <= {
        "embedding1d", "embedding2d", "embedding3d", "knn", "expression", "geneIdColumn",
        "latentKey", "vectorFieldDefault", "obsKeys", "normalizeEmbeddings",
        "centroidMinPoints",
    }


# =============================================================================
# Genes
# =============================================================================


def test_gene_identifiers_keep_those_genes_in_that_order(tmp_path: Path) -> None:
    given = _inputs()
    adata = _read(
        _write(
            tmp_path,
            "data.h5ad",
            "h5ad",
            None,
            gene_identifiers=["GENE5", "GENE0", "GENE9"],
            var_gene_id_column="symbol",
        )
    )
    assert list(adata.var["symbol"]) == ["GENE5", "GENE0", "GENE9"]
    assert list(adata.var_names) == ["ENSG00005", "ENSG00000", "ENSG00009"]
    np.testing.assert_array_equal(
        adata.X.toarray(), given["gene_expression"].toarray()[:, [5, 0, 9]]
    )


def test_genes_are_named_by_the_id_column_when_var_has_no_index(tmp_path: Path) -> None:
    var = _inputs()["var"].reset_index(drop=True)
    adata = _read(_write(tmp_path, "data.h5ad", "h5ad", None, var=var, var_gene_id_column="symbol"))
    assert list(adata.var_names) == list(var["symbol"])


def test_repeated_names_are_made_unique_and_reported(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    given = _inputs()
    obs = given["obs"].copy()
    obs.index = ["a", "a", *obs.index[2:]]
    var = given["var"].copy()
    var.index = ["MATR3", "MATR3", *var.index[2:]]
    adata = _read(_write(tmp_path, "data.h5ad", "h5ad", None, obs=obs, var=var))
    assert list(adata.obs_names[:2]) == ["a", "a-1"]
    assert list(adata.var_names[:2]) == ["MATR3", "MATR3-1"]
    printed = capsys.readouterr().out
    assert "Cell names: 1 duplicate name was made unique, e.g. 'a' → 'a-1'." in printed
    assert "Gene names: 1 duplicate name was made unique, e.g. 'MATR3' → 'MATR3-1'." in printed


def test_without_obs_cells_are_named_by_position(tmp_path: Path) -> None:
    adata = _read(_write(tmp_path, "data.h5ad", "h5ad", None, obs=None))
    assert list(adata.obs_names[:3]) == ["0", "1", "2"]
    assert list(adata.obs.columns) == []


@pytest.mark.parametrize(("name", "output_format", "zarr_format"), OUTPUTS)
def test_only_embeddings_is_a_complete_dataset(
    tmp_path: Path, name: str, output_format: str, zarr_format: int | None
) -> None:
    path = prepare(
        tmp_path / name,
        format=output_format,
        embeddings={"X_umap": np.zeros((5, 2))},
        dataset_id="tiny",
        dataset_name="Tiny",
        **({} if zarr_format is None else {"zarr_format": zarr_format}),
    )
    adata = _read(path)
    assert adata.shape == (5, 0) and adata.X is None
    assert list(adata.obsm) == ["X_umap"]


# =============================================================================
# Refusals: before anything is written
# =============================================================================


REFUSALS = [
    pytest.param({"dataset_id": "has space"}, ValueError, "dataset_id", id="bad-id"),
    pytest.param({"dataset_name": " padded"}, ValueError, "dataset_name", id="bad-name"),
    pytest.param({"citation_url": "https://x.org"}, ValueError, "citation_text is required whenever citation_url", id="url-without-text"),
    pytest.param({"source_name": "S", "source_url": "ftp://x.org"}, ValueError, "absolute http", id="non-http-url"),
    pytest.param({"embeddings": {}}, ValueError, "at least one embedding", id="no-embeddings"),
    pytest.param({"embeddings": {"X/umap": np.zeros((N_CELLS, 2))}}, ValueError, "cannot appear in a file name", id="slash-key"),
    pytest.param({"embeddings": {"X_umap": np.zeros((N_CELLS, 2)), "x_UMAP": np.zeros((N_CELLS, 2))}}, ValueError, "differ only in letter case", id="case-collision"),
    pytest.param({"embeddings": {"X_umap": np.full((N_CELLS, 2), np.nan)}}, ValueError, "every cell's values are missing (NaN or infinite)", id="no-cell-placed"),
    pytest.param({"embeddings": {"X_umap": np.where(np.arange(N_CELLS)[:, None] == 3, [np.nan, 0.0], [1e39, 1.0])}}, ValueError, "within the float32 range", id="placed-cell-beyond-float32"),
    pytest.param({"embeddings": {"X_umap": np.zeros((N_CELLS + 1, 2))}, "vector_fields": None}, ValueError, "obs has 40 rows, but the embeddings have 41", id="row-mismatch"),
    pytest.param({"obs": pd.DataFrame({"_index": range(N_CELLS)})}, ValueError, "'_index' is reserved", id="index-column"),
    pytest.param({"obs": pd.DataFrame({"x": [np.inf] * N_CELLS})}, ValueError, "infinite", id="infinite-obs"),
    pytest.param({"obs": pd.DataFrame({"t": pd.to_datetime(["2026-01-01"] * N_CELLS)})}, TypeError, "no AnnData encoding", id="datetime-obs"),
    pytest.param({"obs": pd.DataFrame({"m": ["a", 1] * (N_CELLS // 2)})}, TypeError, "mixes text", id="mixed-obs"),
    pytest.param({"gene_expression": sparse.csr_matrix(np.full((N_CELLS, N_GENES), np.inf))}, ValueError, "finite numbers", id="infinite-expression"),
    pytest.param({"gene_expression": np.zeros((N_GENES, N_CELLS))}, ValueError, "transpose", id="transposed-expression"),
    pytest.param({"gene_expression": None}, ValueError, "no gene_expression was supplied", id="var-without-expression"),
    pytest.param({"var": None}, ValueError, "genes have no names", id="unnamed-genes"),
    pytest.param({"gene_identifiers": ["NOPE"]}, ValueError, "not found in var", id="unknown-gene"),
    pytest.param({"connectivities": sparse.csr_matrix(np.triu(np.ones((N_CELLS, N_CELLS)), 1))}, ValueError, "exactly symmetric", id="asymmetric-graph"),
    pytest.param({"connectivities": -sparse.csr_matrix(np.ones((N_CELLS, N_CELLS)) - np.eye(N_CELLS))}, ValueError, "non-negative", id="negative-graph"),
    pytest.param({"vector_fields": {"velocity": {"X_pca": np.zeros((N_CELLS, 2))}}}, ValueError, "has 2 columns, and 'X_pca' has 6", id="vectors-narrower-than-embedding"),
    pytest.param({"vector_fields": {"velocity": {"X_tsne": np.zeros((N_CELLS, 2))}}}, ValueError, "'X_tsne' is not one of the embeddings", id="vectors-on-unknown-embedding"),
    pytest.param({"embeddings": {"umap": np.zeros((N_CELLS, 2))}, "vector_fields": {"velocity": {"umap": np.zeros((N_CELLS, 2))}}}, ValueError, "needs an 'X_<basis>' key", id="vectors-on-unprefixed-embedding"),
    pytest.param({"vector_fields": {"X_velocity": {"X_umap": np.zeros((N_CELLS, 2))}}}, ValueError, "obsm['X_velocity_umap'], which the viewer reads as an embedding", id="vector-field-named-like-an-embedding"),
    pytest.param({"embeddings": {"X_umap": np.zeros((N_CELLS, 2)), "X_fwd_umap": np.zeros((N_CELLS, 2))}, "vector_fields": {"T_fwd": {"X_umap": np.zeros((N_CELLS, 2))}}}, ValueError, "reads as the field 'T' on 'X_fwd_umap'", id="vector-key-shadowed-by-a-longer-basis"),
    pytest.param({"vector_fields": {"a_b": {"X_umap": np.zeros((N_CELLS, 2))}, "a": {"X_b_umap": np.zeros((N_CELLS, 2))}}, "embeddings": {"X_umap": np.zeros((N_CELLS, 2)), "X_b_umap": np.zeros((N_CELLS, 2))}}, ValueError, "would both be stored as obsm['a_b_umap']", id="two-fields-one-key"),
    pytest.param({"vector_field_default": "drift"}, ValueError, "must name one of vector_fields", id="unknown-vector-default"),
    pytest.param({"latent_key": "velocity_umap"}, ValueError, "latent_key must name one of the embeddings", id="latent-not-embedding"),
    pytest.param({"embedding_defaults": {"2d": "X_tsne"}}, ValueError, "embedding_defaults['2d'] names obsm['X_tsne']", id="unknown-embedding-default"),
    pytest.param({"embedding_defaults": {"3d": ("X_pca", [0, 0, 1])}}, ValueError, "embedding_defaults['3d']", id="repeated-column"),
    pytest.param({"obs_keys": ["nope"]}, ValueError, "obs_keys names obs['nope']", id="unknown-obs-key"),
    pytest.param({"var_gene_id_column": "nope"}, ValueError, "var_gene_id_column names var['nope']", id="unknown-gene-column"),
    pytest.param({"force": 1}, TypeError, "force", id="non-bool-force"),
]


def test_a_cell_an_embedding_does_not_place_is_written_as_given_and_gets_no_outlier_quantile(tmp_path: Path) -> None:
    """A row holding NaN is a cell the array has no value for (an embedding
    computed on some of the cells): written as given, not drawn by the viewer
    in that embedding, and, in the latent space, given no outlier quantile;
    the other cells are measured without it."""
    rng = np.random.default_rng(3)
    umap = rng.normal(size=(N_CELLS, 2))
    umap[[2, 7]] = np.nan
    latent = rng.normal(size=(N_CELLS, 4))
    latent[5, 1] = np.nan
    obs = pd.DataFrame({"group": pd.Categorical(["a", "b"] * (N_CELLS // 2))})
    adata = _read(
        _write(
            tmp_path, "data.h5ad", "h5ad", None,
            embeddings={"X_umap": umap, "X_latent": latent}, vector_fields=None, obs=obs,
            latent_key="X_latent", centroid_min_points=2,
        )
    )
    written = np.asarray(adata.obsm["X_umap"])
    assert np.isnan(written[[2, 7]]).all()
    assert np.array_equal(written[~np.isnan(umap).any(axis=1)], umap[~np.isnan(umap).any(axis=1)].astype(np.float32))
    quantiles = np.asarray(adata.uns["cellucid"]["outlierQuantiles"]["fields"]["group"])
    expected = latent_outlier_quantiles(latent.astype(np.float32), obs["group"].cat.codes.to_numpy(), 2, 2)
    assert np.isnan(quantiles[5]) and np.array_equal(quantiles, expected, equal_nan=True)


def test_a_self_loop_is_no_edge_seurat_snn_is_written_without_its_diagonal(tmp_path: Path) -> None:
    """Seurat's ``FindNeighbors`` counts each cell among its own neighbours: its
    ``<assay>_snn`` graph has 1 on the diagonal. A self-loop is no edge anyone can
    draw, so the diagonal is left out (as the web app and cellucid-r leave it
    out), never the graph refused."""
    ring = sparse.diags([np.ones(N_CELLS - 1), np.ones(N_CELLS - 1)], [-1, 1], format="csr") * 0.5
    snn = (ring + sparse.identity(N_CELLS, format="csr")).tocsr()
    written = _read(_write(tmp_path, "data.h5ad", "h5ad", None, connectivities=snn)).obsp["connectivities"]
    assert written.diagonal().sum() == 0
    assert (written != ring.astype(np.float32)).nnz == 0


@pytest.mark.parametrize(("overrides", "error", "match"), REFUSALS)
def test_a_defect_is_refused_before_anything_is_written(
    tmp_path: Path, overrides: dict[str, Any], error: type[Exception], match: str
) -> None:
    with pytest.raises(error, match=re.escape(match)):
        _write(tmp_path, "data.h5ad", "h5ad", None, **overrides)
    assert list(tmp_path.iterdir()) == []


#: The four choices and what each is for, as every refusal of ``format`` lists them.
FORMAT_CHOICES = (
    'format must be "h5ad", "zarr", "zarr.zip" or "h5mu":\n'
    '  "h5ad"      one .h5ad file: for single-file hosts (Zenodo, Figshare, Hugging Face, '
    "S3, a lab server), local files, cellucid serve and Jupyter.\n"
    '  "zarr"      a .zarr folder of small files: for GitHub (no file over 100 MB, no LFS).\n'
    '  "zarr.zip"  the Zarr store in one .zarr.zip file: wherever an .h5ad goes.\n'
    '  "h5mu"      one .h5mu file, for a MuData (data=mdata): wherever an .h5ad goes.'
)


def test_format_is_required_and_its_absence_names_the_choices(tmp_path: Path) -> None:
    signature = inspect.signature(prepare)
    assert signature.parameters["format"].default is inspect.Parameter.empty
    assert list(signature.parameters)[:2] == ["out_path", "format"]
    with pytest.raises(TypeError) as raised:
        prepare(tmp_path / "data.h5ad", **_inputs())
    assert str(raised.value) == (
        f"prepare() needs format=, the kind of output to write. {FORMAT_CHOICES}"
    )
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("output_format", ["H5AD", "zarr.ZIP", "zip", ".h5ad", None, 2])
def test_a_format_outside_the_four_is_refused_with_the_choices(
    tmp_path: Path, output_format: object
) -> None:
    with pytest.raises(ValueError) as raised:
        prepare(tmp_path / "data.h5ad", format=output_format, **_inputs())  # type: ignore[arg-type]
    assert str(raised.value) == f"{FORMAT_CHOICES}\nGot {output_format!r}."
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    ("name", "output_format", "message", "instead"),
    [
        pytest.param("data.zarr", "h5ad", 'format="h5ad" writes an .h5ad file, and out_path {path!r} ends in .zarr, which format="zarr" writes: write {instead!r} instead.', "data.h5ad", id="h5ad-to-zarr"),
        pytest.param("data.zarr.zip", "h5ad", 'format="h5ad" writes an .h5ad file, and out_path {path!r} ends in .zarr.zip, which format="zarr.zip" writes: write {instead!r} instead.', "data.h5ad", id="h5ad-to-zip"),
        pytest.param("data.h5ad", "zarr", 'format="zarr" writes a .zarr folder, and out_path {path!r} ends in .h5ad, which format="h5ad" writes: write {instead!r} instead.', "data.zarr", id="zarr-to-h5ad"),
        pytest.param("data.zarr.zip", "zarr", 'format="zarr" writes a .zarr folder, and out_path {path!r} ends in .zarr.zip, which format="zarr.zip" writes: write {instead!r} instead.', "data.zarr", id="zarr-to-zip"),
        pytest.param("data.zarr", "zarr.zip", 'format="zarr.zip" writes a .zarr.zip file, and out_path {path!r} ends in .zarr, which format="zarr" writes: write {instead!r} instead.', "data.zarr.zip", id="zip-to-zarr"),
        pytest.param("data.h5ad", "zarr.zip", 'format="zarr.zip" writes a .zarr.zip file, and out_path {path!r} ends in .h5ad, which format="h5ad" writes: write {instead!r} instead.', "data.zarr.zip", id="zip-to-h5ad"),
        pytest.param("data.loom", "h5ad", 'format="h5ad" writes an .h5ad file, and out_path {path!r} does not end in .h5ad: write {instead!r} instead.', "data.loom.h5ad", id="unknown-suffix"),
        pytest.param(".zarr", "zarr", 'format="zarr" writes a .zarr folder, and out_path {path!r} is only the suffix: write {instead!r} instead.', "data.zarr", id="bare-suffix"),
    ],
)
def test_an_out_path_that_disagrees_with_format_is_refused_with_the_path_to_write(
    tmp_path: Path, name: str, output_format: str, message: str, instead: str
) -> None:
    path = str(tmp_path / name)
    with pytest.raises(ValueError) as raised:
        prepare(path, format=output_format, **_inputs())  # type: ignore[arg-type]
    assert str(raised.value) == message.format(path=path, instead=str(tmp_path / instead))
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    ("name", "output_format", "suffix"),
    [
        pytest.param("Data.H5AD", "h5ad", "Data.H5AD", id="h5ad"),
        pytest.param("Data.ZARR", "zarr", "Data.ZARR", id="zarr"),
        pytest.param("Data.Zarr.Zip", "zarr.zip", "Data.Zarr.Zip", id="zarr-zip"),
    ],
)
def test_the_suffix_matches_format_in_any_letter_case(
    tmp_path: Path, name: str, output_format: str, suffix: str
) -> None:
    assert _write(tmp_path, name, output_format, None).name == suffix


@pytest.mark.parametrize(
    ("name", "output_format", "zarr_format", "match"),
    [
        pytest.param("data.h5ad", "h5ad", 2, 'zarr_format applies to format="zarr" or "zarr.zip", not format="h5ad".', id="zarr-format-for-h5ad"),
        pytest.param("data.zarr", "zarr", 4, "zarr_format must be 2 or 3, got 4.", id="unknown-zarr-format"),
        pytest.param("data.zarr.zip", "zarr.zip", True, "zarr_format must be 2 or 3, got True.", id="bool-zarr-format"),
    ],
)
def test_zarr_format_is_2_or_3_and_only_for_a_zarr_output(
    tmp_path: Path, name: str, output_format: str, zarr_format: Any, match: str
) -> None:
    with pytest.raises(ValueError, match=re.escape(match)):
        _write(tmp_path, name, output_format, zarr_format)
    assert list(tmp_path.iterdir()) == []


def test_an_existing_output_is_replaced_only_with_force(tmp_path: Path) -> None:
    path = _write(tmp_path, "data.zarr", "zarr", 2)
    with pytest.raises(FileExistsError, match="pass force=True"):
        _write(tmp_path, "data.zarr", "zarr", 2)
    _write(tmp_path, "data.zarr", "zarr", 3, force=True)
    assert (path / "zarr.json").is_file() and not (path / ".zgroup").exists()
    assert sorted(item.name for item in tmp_path.iterdir()) == ["data.zarr"]


def test_force_never_replaces_a_folder_that_is_not_a_zarr_store(tmp_path: Path) -> None:
    folder = tmp_path / "results.zarr"
    folder.mkdir()
    (folder / "notes.txt").write_text("keep")
    with pytest.raises(ValueError, match="is not a Zarr store"):
        _write(tmp_path, "results.zarr", "zarr", 2, force=True)
    assert (folder / "notes.txt").read_text() == "keep"
    file_folder = tmp_path / "data.h5ad"
    file_folder.mkdir()
    with pytest.raises(ValueError, match="is a folder"):
        _write(tmp_path, "data.h5ad", "h5ad", None, force=True)


def test_a_failed_write_leaves_the_previous_output_and_no_partial_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write(tmp_path, "data.h5ad", "h5ad", None)
    before = path.read_bytes()

    def fail(self: Any, elements: Any) -> None:
        raise OSError("disk full")

    monkeypatch.setattr(_store._Writer, "write", fail)
    with pytest.raises(OSError, match="disk full"):
        _write(tmp_path, "data.h5ad", "h5ad", None, dataset_name="Other", force=True)
    assert path.read_bytes() == before
    assert sorted(item.name for item in tmp_path.iterdir()) == ["data.h5ad"]


def test_sparse_positions_are_sized_by_their_axis_and_pointers_by_the_values() -> None:
    """Billions of stored values over millions of cells keep 4-byte row positions."""
    from types import SimpleNamespace

    many = SimpleNamespace(format="csc", shape=(2_000_000, 30_000), nnz=4_000_000_000)
    assert _store._index_dtype(many) is np.int64
    assert _store._position_dtype(many) is np.int32
    wide = SimpleNamespace(format="csr", shape=(10, 2**31), nnz=5)
    assert _store._index_dtype(wide) is np.int32
    assert _store._position_dtype(wide) is np.int64


def test_what_an_object_holds_and_prepare_does_not_write_is_named(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """``data=``: a column with no AnnData encoding Cellucid writes (a date)
    is left out with a note, as cellucid-r leaves it out, and so are layers,
    ``raw``, other graphs and obsm tables, each named; the same column given
    as ``obs=`` is refused with the conversion that keeps it."""
    inputs = _inputs()
    adata = ad.AnnData(
        X=inputs["gene_expression"],
        obs=pd.DataFrame(
            {"when": pd.to_datetime(["2026-01-01"] * N_CELLS), "kept": np.arange(N_CELLS)},
            index=[f"c{i}" for i in range(N_CELLS)],
        ),
        var=pd.DataFrame(index=[f"g{i}" for i in range(N_GENES)]),
        obsm={"X_umap": inputs["embeddings"]["X_umap"], "table": pd.DataFrame({"a": np.arange(N_CELLS)}, index=[f"c{i}" for i in range(N_CELLS)])},
        layers={"counts": inputs["gene_expression"]},
        obsp={"distances": sparse.csr_matrix((N_CELLS, N_CELLS))},
    )
    adata.raw = adata.copy()
    path = prepare(tmp_path / "o.zarr", format="zarr", data=adata, dataset_id="o", dataset_name="O")
    said = capsys.readouterr().out
    assert "obs: not stored, no AnnData encoding Cellucid writes: ['when']" in said
    assert "obsm['table'] holds no numeric coordinates and was not read." in said
    assert "layers not read: ['counts']" in said and "raw not read" in said
    assert "obsp ['distances'] not read" in said
    written = ad.read_zarr(path)
    assert "when" not in written.obs and "kept" in written.obs
    # An explicit obs= replaces data's: nothing is said of data's columns, and its date is refused.
    with pytest.raises(TypeError, match="no AnnData encoding"):
        prepare(tmp_path / "e.zarr", format="zarr", data=adata, obs=adata.obs, dataset_id="o", dataset_name="O")
    assert "obs: not stored" not in capsys.readouterr().out


def test_what_a_killed_writer_left_is_removed_by_the_next_one(tmp_path: Path) -> None:
    """A ``prepare()`` killed before it finished leaves its staging beside the
    target (``.<name>.<pid>.partial``); the next one removes it, and the
    target a killed writer had moved aside when the target is there again --
    never a moved-aside copy that is the only one, nor a running writer's."""
    import os

    target = tmp_path / "out.zarr"
    dead = 999999999
    (tmp_path / f".out.zarr.{dead}.partial").mkdir()
    (tmp_path / f".out.zarr.{os.getppid()}.partial").mkdir()  # a running process's
    prepare(target, format="zarr", **_inputs())
    assert not (tmp_path / f".out.zarr.{dead}.partial").exists()
    assert (tmp_path / f".out.zarr.{os.getppid()}.partial").exists()
    (tmp_path / f".out.zarr.{dead}.replaced").mkdir()
    prepare(target, format="zarr", force=True, **_inputs())
    assert not (tmp_path / f".out.zarr.{dead}.replaced").exists()
    # The target gone: the moved-aside copy is the only one, and it stays.
    other = tmp_path / "other.zarr"
    (tmp_path / f".other.zarr.{dead}.replaced").mkdir()
    prepare(other, format="zarr", **_inputs())
    assert (tmp_path / f".other.zarr.{dead}.replaced").exists()
