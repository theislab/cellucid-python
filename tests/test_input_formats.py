"""The input format registry, shared with the web app and cellucid-r.

``src/cellucid/input_formats.json`` is byte for byte the web app's
``tests/fixtures/input-formats.json`` (written from
``assets/js/data/input-formats.js``) and cellucid-r's
``inst/input_formats.json``; each repository pins the same
SHA-256. ``classify_anndata_path`` (``cellucid serve``) applies its rules as
the browser does: a file is taken by its name, then is what its first bytes
say; a server names it with the registry's format id.
"""

from __future__ import annotations

import hashlib
import zipfile
from importlib import resources
from pathlib import Path

import pytest

from cellucid._anndata_files import classify_anndata_path, input_formats

#: The same digest is pinned in the web app
#: (tests/input-formats-registry.test.mjs) and in cellucid-r
#: (tests/testthat/test-input-formats.R). Update all three together.
INPUT_FORMATS_SHA256 = "7c1e079c33dcd7d1abb17e4996e6e6cf35f1968327a9097701dd9a1ca860c902"

HDF5_SIGNATURE = bytes.fromhex("894844460d0a1a0a")


def test_the_packaged_registry_is_the_shared_document() -> None:
    data = resources.files("cellucid").joinpath("input_formats.json").read_bytes()
    assert hashlib.sha256(data).hexdigest() == INPUT_FORMATS_SHA256
    document = input_formats()
    assert document["schema"] == 1
    ids = [entry["id"] for entry in document["formats"]]
    assert ids[:2] == ["h5ad", "zarr"]
    assert {"loom", "tenx-h5", "tenx-outs", "h5mu", "mudata-zarr"} <= set(ids)


def _zarr_zip(path: Path) -> Path:
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("cells.zarr/.zgroup", '{"zarr_format":2}')
    return path


def test_a_file_is_taken_by_its_name_and_is_what_its_bytes_say(tmp_path: Path) -> None:
    h5ad = tmp_path / "cells.h5ad"
    h5ad.write_bytes(HDF5_SIGNATURE + bytes(8))
    assert (classify_anndata_path(h5ad).kind, classify_anndata_path(h5ad).layout) == ("h5ad", "file")

    # HDF5 bytes under a .zip name are an .h5ad; ZIP bytes under an .h5ad name a zipped store.
    misnamed = tmp_path / "cells.zip"
    misnamed.write_bytes(HDF5_SIGNATURE + bytes(8))
    assert (classify_anndata_path(misnamed).kind, classify_anndata_path(misnamed).layout) == ("h5ad", "file")
    archive = _zarr_zip(tmp_path / "store.h5ad")
    store = classify_anndata_path(archive)
    assert (store.kind, store.layout, store.zip_root) == ("zarr", "zip", "cells.zarr")

    # A user block before the HDF5 signature is looked past, as the browser does.
    user_block = tmp_path / "block.h5ad"
    user_block.write_bytes(bytes(512) + HDF5_SIGNATURE)
    assert classify_anndata_path(user_block).kind == "h5ad"

    # A name no format has is refused before anything is read.
    other = tmp_path / "cells.txt"
    other.write_bytes(HDF5_SIGNATURE)
    with pytest.raises(ValueError, match=r"must be an \.h5ad, \.zarr\.zip, "):
        classify_anndata_path(other)


def test_an_h5seurat_file_is_refused_in_the_registry_words(tmp_path: Path) -> None:
    import h5py

    from cellucid._anndata_files import RefusedFormatError

    refusal = next(entry for entry in input_formats()["formats"] if entry["id"] == "h5seurat")["refusal"]
    assert "SeuratDisk::LoadH5Seurat" in refusal and "saveRDS" in refusal
    path = tmp_path / "pbmc.h5Seurat"
    with h5py.File(path, "w") as file:
        file.create_group("assays")
        file.create_dataset("cell.names", data=[b"AAAC-1"])
    with pytest.raises(RefusedFormatError) as caught:
        classify_anndata_path(path)
    assert str(caught.value) == f"pbmc.h5Seurat: {refusal}"

    # An AnnData that SeuratDisk's Convert() named .h5Seurat is what its root says.
    converted = tmp_path / "converted.h5Seurat"
    with h5py.File(converted, "w") as file:
        file.attrs["encoding-type"] = "anndata"
    assert classify_anndata_path(converted).kind == "h5ad"


def test_cellucid_serve_refuses_an_h5seurat_file_before_serving(tmp_path: Path, capsys) -> None:
    import h5py

    from cellucid.cli import main

    path = tmp_path / "pbmc.h5Seurat"
    with h5py.File(path, "w") as file:
        file.create_group("assays")
    assert main(["serve", str(path), "--quiet"]) == 1
    error = capsys.readouterr().err
    assert "This is an h5Seurat file (SeuratDisk), which Cellucid does not open." in error
    assert "Unable to detect format" not in error
