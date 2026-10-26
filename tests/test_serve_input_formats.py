"""``cellucid serve`` serves every format of the input format registry by path.

The browser reads every format itself (``input_formats.json``), so a server
names the format with the registry's id (``anndata.file.kind``), serves the
bytes -- one file with exact ranges, sent by the kernel; a Zarr store folder
key by key; a folder of files file by file, listed because a browser cannot
list a folder -- and reports cell and gene counts only where it can read them
from metadata (null for an R object, a MatrixMarket folder, a Loom or 10x
file: the browser learns them when it opens the data).
"""

from __future__ import annotations

import gzip
import http.client
import json
import zipfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest

from cellucid._anndata_files import (
    browser_undecodable,
    classify_anndata_path,
    folder_format_of_paths,
    read_anndata_contents,
    read_cellucid_identity,
)
from cellucid.anndata_server import AnnDataServer
from cellucid.cli import _detect_data_format


def _write(path: Path, data: bytes) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return path


def _loom(path: Path) -> Path:
    with h5py.File(path, "w") as file:
        file.create_dataset("matrix", data=np.arange(12, dtype=np.float32).reshape(3, 4))
        file.create_group("row_attrs").create_dataset("Gene", data=np.array([b"a", b"b", b"c"]))
        file.create_group("col_attrs").create_dataset("CellID", data=np.array([b"1", b"2", b"3", b"4"]))
        file.attrs["LOOM_SPEC_VERSION"] = "3.0.0"
    return path


def _tenx_h5(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as file:
        file.attrs["filetype"] = "matrix"
        file.attrs["version"] = 2
        matrix = file.create_group("matrix")
        matrix.create_dataset("barcodes", data=np.array([b"AAAC-1", b"AAAG-1"]))
        matrix.create_dataset("data", data=np.array([1, 2], dtype=np.int32))
        matrix.create_dataset("indices", data=np.array([0, 1], dtype=np.int64))
        matrix.create_dataset("indptr", data=np.array([0, 1, 2], dtype=np.int64))
        matrix.create_dataset("shape", data=np.array([2, 2], dtype=np.int32))
    return path


def _mudata_h5mu(path: Path) -> Path:
    """A MuData store written the way mudata writes one (no mudata needed)."""
    from anndata.io import write_elem

    rng = np.random.default_rng(3)
    obs_names = pd.Index([f"c{i}" for i in range(6)], dtype=object)
    rna = ad.AnnData(
        X=rng.random((6, 4)).astype(np.float32),
        obs=pd.DataFrame({"leiden": pd.Categorical(list("aabbcc"))}, index=obs_names),
        var=pd.DataFrame(index=pd.Index([f"g{i}" for i in range(4)], dtype=object)),
    )
    adt = ad.AnnData(
        X=rng.random((6, 2)).astype(np.float32),
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=pd.Index(["CD3", "CD4"], dtype=object)),
    )
    adt.obsm["X_umap"] = rng.random((6, 2)).astype(np.float32)
    with h5py.File(path, "w", userblock_size=512) as file:
        file.attrs["encoding-type"] = "MuData"
        file.attrs["encoding-version"] = "0.1.0"
        file.attrs["axis"] = 0
        write_elem(file, "obs", pd.DataFrame(index=obs_names))
        write_elem(file, "var", pd.DataFrame(index=pd.Index([*rna.var_names, *adt.var_names], dtype=object)))
        write_elem(file, "obsm", {"X_umap": rng.random((6, 2)).astype(np.float32), "rna": np.ones(6, bool), "prot": np.ones(6, bool)})
        mod = file.create_group("mod")
        write_elem(mod, "prot", adt)
        write_elem(mod, "rna", rna)
        mod.attrs["mod-order"] = ["prot", "rna"]
        write_elem(file, "uns", {"cellucid": {"id": "cite-demo", "name": "CITE demo"}})
    return path


def _mtx_folder(folder: Path, prefix: str = "", *, features: str = "features") -> Path:
    _write(folder / f"{prefix}matrix.mtx.gz", gzip.compress(
        b"%%MatrixMarket matrix coordinate integer general\n2 3 2\n1 1 4\n2 3 1\n"
    ))
    _write(folder / f"{prefix}barcodes.tsv.gz", gzip.compress(b"AAAC-1\nAAAG-1\nAAAT-1\n"))
    _write(folder / f"{prefix}{features}.tsv.gz", gzip.compress(b"ENSG1\tA\tGene Expression\nENSG2\tB\tGene Expression\n"))
    return folder


def _rds(path: Path) -> Path:
    # saveRDS output: a gzip stream of R's XDR serialization. The server never
    # parses it, so the header is enough to name it; GEO gzips it once more.
    return _write(path, gzip.compress(b"X\n\x00\x00\x00\x03" + bytes(64)))


@contextmanager
def _serving(path: Path, **identity: str) -> Iterator[tuple[AnnDataServer, str, int]]:
    server = AnnDataServer(
        path,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        **({"dataset_name": "Served", "dataset_id": "served"} | identity),
    )
    server.start_background()
    assert server._server is not None
    host, port = server._server.server_address[:2]
    try:
        yield server, host, port
    finally:
        server.stop()


def _get(host: str, port: int, path: str, headers: dict[str, str] | None = None) -> tuple[int, dict[str, str], bytes]:
    connection = http.client.HTTPConnection(host, port, timeout=10)
    try:
        connection.request("GET", path, headers=headers or {})
        response = connection.getresponse()
        return response.status, dict(response.getheaders()), response.read()
    finally:
        connection.close()


# =============================================================================
# Classification
# =============================================================================


@pytest.mark.parametrize(
    ("make", "kind", "layout", "ending"),
    [
        (lambda tmp: _loom(tmp / "cells.loom"), "loom", "file", ".loom"),
        (lambda tmp: _tenx_h5(tmp / "filtered_feature_bc_matrix.h5"), "tenx-h5", "file", ".h5"),
        (lambda tmp: _mudata_h5mu(tmp / "cite.h5mu"), "h5mu", "file", ".h5mu"),
        # The root says MuData whatever the name says.
        (lambda tmp: _mudata_h5mu(tmp / "cite.h5"), "h5mu", "file", ".h5mu"),
        (lambda tmp: _rds(tmp / "pbmc.rds"), "rds", "file", ".rds"),
        (lambda tmp: _write(tmp / "GSM1_obj.rds.gz", gzip.compress(gzip.compress(b"X\n"))), "rds", "file", ".rds.gz"),
        (lambda tmp: _rds(tmp / "objects.RData"), "rda", "file", ".rdata"),
        (lambda tmp: _rds(tmp / "objects.rda.gz"), "rda", "file", ".rda.gz"),
    ],
)
def test_every_registry_file_format_is_named_by_its_registry_id(
    tmp_path: Path, make, kind: str, layout: str, ending: str
) -> None:
    path = make(tmp_path)
    store = classify_anndata_path(path)
    assert (store.kind, store.layout, store.ending) == (kind, layout, ending)
    assert _detect_data_format(path) == kind


def test_an_h5ad_named_h5_is_an_h5ad_by_its_root(tmp_path: Path) -> None:
    path = tmp_path / "cells.h5"
    ad.AnnData(X=np.ones((2, 2), np.float32)).write_h5ad(path)
    store = classify_anndata_path(path)
    assert (store.kind, store.ending) == ("h5ad", ".h5ad")


def test_folders_are_named_by_their_files(tmp_path: Path) -> None:
    mtx = classify_anndata_path(_mtx_folder(tmp_path / "filtered_feature_bc_matrix"))
    assert (mtx.kind, mtx.layout) == ("mtx", "folder")
    assert mtx.files == ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")

    # GEO names each sample's files with its accession.
    geo = classify_anndata_path(
        _mtx_folder(tmp_path / "geo", "GSM4567_PBMC_", features="genes")
    )
    assert geo.kind == "mtx"

    outs = tmp_path / "outs"
    _tenx_h5(outs / "filtered_feature_bc_matrix.h5")
    _write(outs / "analysis" / "umap" / "2_components" / "projection.csv", b"Barcode,UMAP-1,UMAP-2\n")
    _write(outs / ".DS_Store", b"")
    store = classify_anndata_path(outs)
    assert (store.kind, store.layout) == ("tenx-outs", "folder")
    assert store.files == (
        "analysis/umap/2_components/projection.csv",
        "filtered_feature_bc_matrix.h5",
    )

    zarr_folder = tmp_path / "store.zarr"
    ad.AnnData(X=np.ones((2, 2), np.float32)).write_zarr(zarr_folder)
    assert classify_anndata_path(zarr_folder).kind == "zarr"
    assert classify_anndata_path(zarr_folder).files is None

    empty = tmp_path / "nothing"
    _write(empty / "notes.txt", b"")
    with pytest.raises(ValueError, match="not a folder Cellucid reads.*matrix.mtx.gz"):
        classify_anndata_path(empty)


def test_the_folder_rule_is_the_web_apps(tmp_path: Path) -> None:
    # Cases of the web app's folderFormatOfPaths (input-formats.js).
    assert folder_format_of_paths(["x.zarr/.zgroup"]) == "zarr"
    assert folder_format_of_paths(["outs/filtered_feature_bc_matrix.h5", "outs/matrix.mtx.gz"]) == "tenx-outs"
    assert folder_format_of_paths(["f/matrix.mtx"]) == "mtx"
    assert folder_format_of_paths(["f/GSM1_matrix.mtx.gz", "f/GSM1_barcodes.tsv.gz"]) == "mtx"
    assert folder_format_of_paths(["f/GSM1-matrix.mtx.gz", "f/GSM1-features.tsv.gz"]) == "mtx"
    # A prefixed marker alone, with nothing else sharing its prefix, is not a folder.
    assert folder_format_of_paths(["f/GSM1_matrix.mtx.gz", "f/GSM2_barcodes.tsv.gz"]) is None
    assert folder_format_of_paths(["f/notes.txt"]) is None


def test_contents_and_identity_are_read_only_where_python_can_read_them(tmp_path: Path) -> None:
    rds = classify_anndata_path(_rds(tmp_path / "pbmc.rds"))
    assert read_anndata_contents(rds) is None
    assert read_cellucid_identity(rds) == (None, None)
    assert browser_undecodable(rds) == []
    loom = classify_anndata_path(_loom(tmp_path / "cells.loom"))
    assert read_anndata_contents(loom) is None

    mudata = classify_anndata_path(_mudata_h5mu(tmp_path / "cite.h5mu"))
    contents = read_anndata_contents(mudata)
    assert contents is not None and contents.modalities is not None
    assert list(contents.modalities) == ["prot", "rna"]
    assert (contents.n_obs, contents.n_vars) == (6, 4)  # rna is the default modality
    assert {"rna:leiden", "rna:measured", "prot:measured"} <= contents.obs_columns
    assert set(contents.obsm) == {"X_umap", "prot:X_umap"}
    assert read_cellucid_identity(mudata) == ("cite-demo", "CITE demo")


def test_a_mudata_that_shares_features_is_refused_with_what_to_do(tmp_path: Path) -> None:
    path = _mudata_h5mu(tmp_path / "axis.h5mu")
    with h5py.File(path, "r+") as file:
        file.attrs["axis"] = 1
    with pytest.raises(ValueError, match=r"shares features, not cells.*axis=0"):
        read_anndata_contents(classify_anndata_path(path))


# =============================================================================
# Serving
# =============================================================================


def test_an_r_object_is_served_as_its_bytes_with_ranges_and_null_counts(tmp_path: Path) -> None:
    path = _write(tmp_path / "GSM1_seurat.rds.gz", gzip.compress(gzip.compress(b"X\n" + bytes(4000))))
    with _serving(path) as (server, host, port):
        listing = json.loads(_get(host, port, "/_cellucid/datasets")[2])
        file = listing["datasets"][0]["anndata"]["file"]
        assert file == {
            "url": "_cellucid/anndata/served.rds.gz",
            "kind": "rds",
            "layout": "file",
            "size": path.stat().st_size,
        }
        health = json.loads(_get(host, port, "/_cellucid/health")[2])
        assert health["format"] == "rds"
        assert (health["n_cells"], health["n_genes"]) == (None, None)
        status, headers, body = _get(host, port, f"/{file['url']}", {"Range": "bytes=0-9"})
        assert status == 206
        assert headers["Content-Range"] == f"bytes 0-9/{path.stat().st_size}"
        assert body == path.read_bytes()[:10]


@pytest.mark.parametrize("prefix", ["", "GSM4567_"])
def test_a_matrix_market_folder_lists_its_files_and_serves_each(tmp_path: Path, prefix: str) -> None:
    folder = _mtx_folder(tmp_path / "sample", prefix)
    with _serving(folder) as (_server, host, port):
        file = json.loads(_get(host, port, "/_cellucid/datasets")[2])["datasets"][0]["anndata"]["file"]
        assert file["url"] == "_cellucid/anndata/served/"
        assert (file["kind"], file["layout"]) == ("mtx", "folder")
        assert file["files"] == [
            {"path": name, "size": (folder / name).stat().st_size}
            for name in sorted(path.name for path in folder.iterdir())
        ]
        for entry in file["files"]:
            status, _headers, body = _get(host, port, f"/{file['url']}{entry['path']}")
            assert status == 200
            assert body == (folder / entry["path"]).read_bytes()
        assert _get(host, port, f"/{file['url']}../served/x")[0] == 404
        assert _get(host, port, f"/{file['url']}missing.tsv.gz")[0] == 404


def test_a_cell_ranger_outs_folder_is_served_with_its_subfolders(tmp_path: Path) -> None:
    outs = tmp_path / "outs"
    _tenx_h5(outs / "filtered_feature_bc_matrix.h5")
    projection = _write(outs / "analysis" / "umap" / "2_components" / "projection.csv", b"Barcode,UMAP-1,UMAP-2\nAAAC-1,0,1\n")
    with _serving(outs) as (_server, host, port):
        file = json.loads(_get(host, port, "/_cellucid/datasets")[2])["datasets"][0]["anndata"]["file"]
        assert file["kind"] == "tenx-outs"
        paths = [entry["path"] for entry in file["files"]]
        assert "analysis/umap/2_components/projection.csv" in paths
        status, _headers, body = _get(host, port, f"/{file['url']}analysis/umap/2_components/projection.csv")
        assert (status, body) == (200, projection.read_bytes())


def test_a_mudata_names_itself_and_is_served_with_its_counts(tmp_path: Path) -> None:
    path = _mudata_h5mu(tmp_path / "cite.h5mu")
    server = AnnDataServer(path, host="127.0.0.1", port=0, quiet=True, serve_web_ui=False)
    assert (server.dataset_id, server.dataset_name) == ("cite-demo", "CITE demo")
    assert server._datasets_document()["datasets"][0]["anndata"]["file"]["url"] == (
        "_cellucid/anndata/cite-demo.h5mu"
    )
    assert (server.server_info["n_cells"], server.server_info["n_genes"]) == (6, 4)
    assert server.server_info["format"] == "h5mu"


def test_a_format_python_cannot_read_must_be_named(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="dataset_name and dataset_id are required"):
        AnnDataServer(_rds(tmp_path / "pbmc.rds"), host="127.0.0.1", port=0, quiet=True, serve_web_ui=False)


def test_a_zipped_mudata_zarr_is_named_by_its_root(tmp_path: Path) -> None:
    archive = tmp_path / "cite.zarr.zip"
    with zipfile.ZipFile(archive, "w") as file:
        file.writestr("cite.zarr/.zgroup", '{"zarr_format":2}')
        file.writestr("cite.zarr/.zattrs", '{"encoding-type":"MuData","axis":0}')
    store = classify_anndata_path(archive)
    assert (store.kind, store.layout, store.zip_root) == ("mudata-zarr", "zip", "cite.zarr")


def test_an_in_memory_mudata_is_served_from_a_temporary_cellucid_ready_h5mu() -> None:
    import mudata
    rng = np.random.default_rng(5)
    names = pd.Index([f"c{i}" for i in range(8)], dtype=object)
    rna = ad.AnnData(X=rng.random((8, 3)).astype(np.float32), obs=pd.DataFrame(index=names))
    rna.obsm["X_umap"] = rng.random((8, 2)).astype(np.float32)
    adt = ad.AnnData(X=rng.random((8, 2)).astype(np.float32), obs=pd.DataFrame(index=names))
    adt.var_names = ["CD3", "CD4"]
    mdata = mudata.MuData({"rna": rna, "prot": adt})
    mdata.obsm["X_umap"] = rng.random((8, 2)).astype(np.float32)
    server = AnnDataServer(
        mdata,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="CITE",
        dataset_id="cite",
        modality="prot",
    )
    try:
        temporary = server.store.path
        assert (server.store.kind, temporary.suffix) == ("h5mu", ".h5mu")
        assert (server.server_info["n_cells"], server.server_info["format"]) == (8, "h5mu")
        assert server.defaults["modality"] == "prot"
        contents = read_anndata_contents(server.store)
        assert contents is not None and set(contents.modalities or {}) == {"rna", "prot"}
    finally:
        server.stop()
    assert not temporary.exists()
