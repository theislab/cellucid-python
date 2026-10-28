"""``prepare(data=...)``: an AnnData or a MuData whole, and the Cellucid-ready MuData.

- ``data=adata`` writes byte for byte what the same components passed one by
  one write, and an explicit argument replaces one component;
- ``data=mdata`` writes a Cellucid-ready MuData (``.h5mu``, MuData Zarr folder
  and ``.zarr.zip``) that mudata reads back value for value, and that the web
  app's own reader (``cellucid/assets/js/data/anndata``, run with Node from the
  sibling checkout) reads as mudata does: the modalities, their genes and
  values (NaN where a modality did not measure a cell), the embeddings placed
  among the cells, the modality it opens on;
- every MuData Cellucid cannot show is refused before anything is written,
  with what to do;
- ``modality`` is a launch default, checked like the others;
- a session applies to a MuData as the viewer showed it;
- the R writer's MuData Zarr matches this one (``cross_writer_mudata.py``),
  once cellucid-r ships its fixture.

The written contract is §5 of
``docs/user_guide/python_package/c_data_preparation_api/09_cellucid_ready_anndata_specification.md``.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import warnings
from pathlib import Path
from typing import Any

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pytest
import zarr
from scipy import sparse

from cellucid import prepare, restore_original_order
from cellucid._anndata_files import AnnDataContents
from cellucid._launch_defaults import launch_defaults

REPOSITORY = Path(__file__).resolve().parents[1]
WEB_ROOT = REPOSITORY.parent / "cellucid"
R_CROSS_WRITER = REPOSITORY.parent / "cellucid-r" / "tests" / "testthat" / "fixtures" / "cross_writer_mudata"


# =============================================================================
# AnnData
# =============================================================================


def _adata() -> ad.AnnData:
    rng = np.random.default_rng(11)
    n = 30
    counts = sparse.random(n, 8, density=0.4, format="csr", random_state=2, dtype=np.float64)
    counts.data = np.floor(counts.data * 50) + 1
    graph = sparse.random(n, n, density=0.1, format="csr", random_state=4)
    graph = graph + graph.T
    graph.setdiag(0)
    graph.eliminate_zeros()
    adata = ad.AnnData(
        X=counts,
        obs=pd.DataFrame(
            {"cell_type": rng.choice(["T", "B"], n), "score": rng.normal(size=n)},
            index=[f"cell{i}" for i in range(n)],
        ),
        var=pd.DataFrame({"symbol": [f"S{i}" for i in range(8)]}, index=[f"G{i}" for i in range(8)]),
    )
    adata.obsm["X_umap"] = rng.normal(size=(n, 2))
    adata.obsm["X_pca"] = rng.normal(size=(n, 5))
    adata.obsm["velocity_umap"] = rng.normal(size=(n, 2))
    adata.obsm["protein_table"] = pd.DataFrame({"CD3": rng.normal(size=n)}, index=adata.obs_names)
    adata.obsp["connectivities"] = graph
    adata.obsp["distances"] = graph
    adata.layers["counts"] = counts
    adata.uns["cellucid"] = {
        "id": "toy-anndata",
        "name": "Toy AnnData",
        "description": "Thirty cells.",
        "citation": {"text": "Doe 2026", "url": "https://example.org/doe"},
        "defaults": {"embedding2d": "X_umap", "latentKey": "X_pca", "obsKeys": np.asarray(["cell_type"], dtype=object)},
    }
    return adata


def _files(store: Path) -> dict[str, bytes]:
    return {
        path.relative_to(store).as_posix(): path.read_bytes()
        for path in sorted(store.rglob("*"))
        if path.is_file()
    }


@pytest.mark.parametrize("zarr_format", [2, 3])
def test_data_writes_what_the_same_components_write_byte_for_byte(tmp_path: Path, zarr_format: int) -> None:
    adata = _adata()
    from_data = prepare(tmp_path / "data.zarr", format="zarr", zarr_format=zarr_format, data=adata)
    by_hand = prepare(
        tmp_path / "hand.zarr",
        format="zarr",
        zarr_format=zarr_format,
        embeddings={"X_umap": adata.obsm["X_umap"], "X_pca": adata.obsm["X_pca"]},
        vector_fields={"velocity": {"X_umap": adata.obsm["velocity_umap"]}},
        dataset_id="toy-anndata",
        dataset_name="Toy AnnData",
        dataset_description="Thirty cells.",
        citation_text="Doe 2026",
        citation_url="https://example.org/doe",
        obs=adata.obs,
        var=adata.var,
        gene_expression=adata.X,
        connectivities=adata.obsp["connectivities"],
        embedding_defaults={"2d": "X_umap"},
        latent_key="X_pca",
        obs_keys=["cell_type"],
    )
    assert _files(from_data) == _files(by_hand)
    zipped_data = prepare(tmp_path / "data.zarr.zip", format="zarr.zip", data=adata)
    zipped_hand = prepare(tmp_path / "hand.zarr.zip", format="zarr.zip", data=ad.read_zarr(by_hand) if zarr_format == 2 else adata)
    assert zipped_data.read_bytes() == zipped_hand.read_bytes()


def test_data_takes_every_component_it_holds_and_leaves_out_what_a_cellucid_ready_anndata_has_none_of(
    tmp_path: Path,
) -> None:
    adata = _adata()
    written = restore_original_order(ad.read_h5ad(prepare(tmp_path / "data.h5ad", format="h5ad", data=adata)))
    assert sorted(written.obsm) == ["X_pca", "X_umap", "velocity_umap"]
    assert list(written.obsp) == ["connectivities"]
    assert list(written.layers) == []
    assert written.obs_names.tolist() == adata.obs_names.tolist()
    assert written.var_names.tolist() == adata.var_names.tolist()
    np.testing.assert_array_equal(written.X.toarray(), adata.X.toarray())
    assert written.uns["cellucid"]["defaults"]["latentKey"] == "X_pca"
    assert written.uns["cellucid"]["defaults"]["embedding2d"] == "X_umap"
    assert list(written.uns["cellucid"]["defaults"]["obsKeys"]) == ["cell_type"]
    assert written.uns["cellucid"]["citation"]["url"] == "https://example.org/doe"
    assert "outlierQuantiles" in written.uns["cellucid"]


def test_an_explicit_argument_replaces_one_component(tmp_path: Path) -> None:
    adata = _adata()
    written = ad.read_h5ad(
        prepare(
            tmp_path / "data.h5ad",
            format="h5ad",
            data=adata,
            embeddings={"X_pca": adata.obsm["X_pca"]},
            vector_fields={},
            embedding_defaults={"2d": ("X_pca", [1, 0])},
            gene_expression=adata.layers["counts"][:, :4],
            var=adata.var.iloc[:4],
            dataset_name="Renamed",
        )
    )
    assert list(written.obsm) == ["X_pca"]
    assert written.n_vars == 4
    assert written.uns["title"] == "Renamed"
    assert written.uns["cellucid"]["id"] == "toy-anndata"
    assert dict(written.uns["cellucid"]["defaults"]["embedding2d"])["key"] == "X_pca"


def test_a_prepared_file_prepares_again_with_no_identity(tmp_path: Path) -> None:
    first = prepare(tmp_path / "first.h5ad", format="h5ad", data=_adata())
    again = ad.read_zarr(prepare(tmp_path / "again.zarr", format="zarr", data=ad.read_h5ad(first)))
    assert again.uns["cellucid"]["id"] == "toy-anndata"
    assert again.uns["cellucid"]["defaults"]["latentKey"] == "X_pca"


@pytest.mark.parametrize(
    ("change", "error", "message"),
    [
        pytest.param(lambda a: a.uns.pop("cellucid"), TypeError, "prepare() needs dataset_id: data's uns['cellucid'] has no 'id', so pass dataset_id=.", id="no-identity"),
        pytest.param(lambda a: [a.obsm.pop(key) for key in list(a.obsm)], ValueError, "data holds no embedding: its obsm has no numeric array.", id="no-embedding"),
    ],
)
def test_what_data_lacks_is_named(tmp_path: Path, change: Any, error: type[Exception], message: str) -> None:
    adata = _adata()
    change(adata)
    with pytest.raises(error, match=re.escape(message)):
        prepare(tmp_path / "data.h5ad", format="h5ad", data=adata)
    assert list(tmp_path.iterdir()) == []


def test_without_data_the_identity_and_an_embedding_are_required(tmp_path: Path) -> None:
    with pytest.raises(TypeError, match=re.escape("prepare() needs embeddings=")):
        prepare(tmp_path / "data.h5ad", format="h5ad", dataset_id="x", dataset_name="X")
    with pytest.raises(TypeError, match=re.escape("prepare() needs dataset_name: pass dataset_name=.")):
        prepare(tmp_path / "data.h5ad", format="h5ad", dataset_id="x", embeddings={"X_umap": np.zeros((2, 2))})


@pytest.mark.parametrize("data", [np.zeros((3, 3)), {"X": np.zeros((3, 3))}, "data.h5ad"])
def test_data_is_an_anndata_or_a_mudata(tmp_path: Path, data: Any) -> None:
    with pytest.raises(TypeError, match=re.escape("data must be an anndata.AnnData or a mudata.MuData")):
        prepare(tmp_path / "data.h5ad", format="h5ad", data=data)


def test_an_anndata_is_not_written_as_h5mu_and_has_no_modality(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match=re.escape('format="h5mu" writes a MuData: pass the MuData as data=')):
        prepare(tmp_path / "data.h5mu", format="h5mu", data=_adata())
    with pytest.raises(ValueError, match=re.escape("modality names 'rna', and this data has one feature space")):
        prepare(tmp_path / "data.h5ad", format="h5ad", data=_adata(), modality="rna")
    assert list(tmp_path.iterdir()) == []


# =============================================================================
# The modality launch default
# =============================================================================


def _contents(modalities: list[str] | None = None) -> AnnDataContents:
    one = AnnDataContents(
        n_obs=3,
        n_vars=2,
        obs_columns=frozenset(),
        var_columns=frozenset({"symbol"}),
        obsm={},
        obsp=frozenset(),
        layers=frozenset(),
        has_x=True,
        raw_var_columns=None,
    )
    return AnnDataContents(
        modalities=None if modalities is None else dict.fromkeys(modalities, one),
        n_obs=3,
        n_vars=2,
        obs_columns=frozenset(),
        var_columns=frozenset({"symbol"}),
        obsm={"X_umap": 2},
        obsp=frozenset(),
        layers=frozenset(),
        has_x=True,
        raw_var_columns=None,
    )


def _defaults(contents: AnnDataContents | None = None, **changes: Any) -> dict[str, Any]:
    arguments: dict[str, Any] = {
        "latent_key": None,
        "gene_id_column": None,
        "normalize_embeddings": None,
        "centroid_min_points": None,
        "obs_keys": None,
        "embedding_keys": None,
        "knn_key": None,
        "expression_source": None,
        "vector_field_default": None,
    }
    arguments.update(changes)
    return launch_defaults(_contents() if contents is None else contents, **arguments)


def test_the_modality_default_names_one_of_the_modalities() -> None:
    multimodal = _contents(["rna", "prot"])
    assert _defaults(multimodal, modality="prot") == {"modality": "prot"}
    assert _defaults(multimodal) == {}
    with pytest.raises(ValueError, match=re.escape("modality names 'atac', which this data does not have. Its modalities: ['rna', 'prot'].")):
        _defaults(multimodal, modality="atac")
    with pytest.raises(ValueError, match=re.escape("modality names 'rna', and this data has one feature space")):
        _defaults(modality="rna")
    with pytest.raises(ValueError, match=re.escape("modality must be None or one non-blank string.")):
        _defaults(multimodal, modality="  ")


def test_a_store_python_does_not_read_has_its_choices_checked_as_written() -> None:
    choices: dict[str, Any] = {
        "latent_key": "X_pca",
        "gene_id_column": "symbol",
        "normalize_embeddings": False,
        "centroid_min_points": 4,
        "obs_keys": ["cell_type"],
        "embedding_keys": {"2d": ("X_pca", [3, 1])},
        "knn_key": "connectivities",
        "expression_source": "layer:counts",
        "vector_field_default": "velocity",
        "modality": "rna",
    }
    assert launch_defaults(None, **choices) == {
        "modality": "rna",
        "embedding2d": {"key": "X_pca", "columns": [3, 1]},
        "knn": "connectivities",
        "expression": "layers/counts",
        "geneIdColumn": "symbol",
        "latentKey": "X_pca",
        "vectorFieldDefault": "velocity",
        "obsKeys": ["cell_type"],
        "normalizeEmbeddings": False,
        "centroidMinPoints": 4,
    }
    for change, message in (
        ({"embedding_keys": {"2d": ("X_pca", [1, 1])}}, "embedding_keys['2d'] (obsm['X_pca'])"),
        ({"embedding_keys": {"2d": ("X_pca", [0])}}, "embedding_keys['2d'] (obsm['X_pca'])"),
        ({"expression_source": "counts"}, "expression_source must be 'X', 'raw', or 'layer:<name>'"),
        ({"obs_keys": ["a", "a"]}, "obs_keys names a column more than once."),
        ({"centroid_min_points": 0}, "centroid_min_points"),
        ({"modality": ""}, "modality must be None or one non-blank string."),
    ):
        with pytest.raises((TypeError, ValueError), match=re.escape(message)):
            launch_defaults(None, **{**choices, **change})


# =============================================================================
# MuData
# =============================================================================


def _mudata(*, prot_cells: int = 40, order: tuple[str, ...] = ("prot", "rna")) -> Any:
    rng = np.random.default_rng(5)
    n = 50
    names = [f"c{i:02d}" for i in range(n)]
    counts = sparse.random(n, 20, density=0.3, format="csr", random_state=3, dtype=np.float32)
    counts.data = np.floor(counts.data * 9) + 1
    graph = sparse.random(n, n, density=0.08, format="csr", random_state=6)
    graph = graph + graph.T
    graph.setdiag(0)
    graph.eliminate_zeros()
    rna = ad.AnnData(
        X=counts,
        obs=pd.DataFrame({"leiden": pd.Categorical(rng.choice(["0", "1", "2"], n))}, index=names),
        var=pd.DataFrame({"symbol": [f"SYM{i}" for i in range(20)]}, index=[f"ENSG{i:03d}" for i in range(20)]),
    )
    rna.obsm["X_umap"] = rng.normal(size=(n, 2))
    rna.obsm["X_pca"] = rng.normal(size=(n, 6))
    rna.obsp["connectivities"] = graph
    measured = names[n - prot_cells :][::-1]
    prot = ad.AnnData(
        X=rng.normal(size=(prot_cells, 4)).astype(np.float32),
        obs=pd.DataFrame({"isotype_ratio": rng.random(prot_cells)}, index=measured),
        var=pd.DataFrame(index=["CD3", "CD4", "CD8", "CD19"]),
    )
    prot.obsm["X_umap"] = rng.normal(size=(prot_cells, 2))
    modalities = {"prot": prot, "rna": rna}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mdata = md.MuData({name: modalities[name] for name in order})
    mdata.obsm["X_wnn_umap"] = rng.normal(size=(n, 2))
    mdata.obs["sample"] = pd.Categorical(rng.choice(["s1", "s2"], n))
    mdata.uns["cellucid"] = {"id": "toy-citeseq", "name": "Toy CITE-seq"}
    return mdata


MUDATA_OUTPUTS = [
    pytest.param("m.h5mu", "h5mu", None, id="h5mu"),
    pytest.param("m.zarr", "zarr", 2, id="zarr-v2"),
    pytest.param("m.zarr", "zarr", 3, id="zarr-v3"),
    pytest.param("m.zarr.zip", "zarr.zip", None, id="zarr-zip"),
]


def _read_mudata(path: Path) -> Any:
    """The written MuData in the order its cells were given (``restore_original_order``)."""
    return restore_original_order(_stored_mudata(path))


def _stored_mudata(path: Path) -> Any:
    """The written MuData as stored: its cells in priority order."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        if path.suffix == ".h5mu":
            return md.read_h5mu(path)
        if path.name.endswith(".zip"):
            return md.read_zarr(zarr.storage.ZipStore(path, mode="r"))
        return md.read_zarr(path)


def _write_mudata(tmp_path: Path, name: str, output_format: str, zarr_format: int | None, mdata: Any = None, **arguments: Any) -> Path:
    if zarr_format is not None:
        arguments["zarr_format"] = zarr_format
    return prepare(tmp_path / name, format=output_format, data=_mudata() if mdata is None else mdata, **arguments)


def _dense(matrix: Any) -> np.ndarray:
    return np.asarray(matrix.toarray() if sparse.issparse(matrix) else matrix)


@pytest.mark.parametrize(("name", "output_format", "zarr_format"), MUDATA_OUTPUTS)
def test_mudata_reads_back_every_modality_value_for_value(
    tmp_path: Path, name: str, output_format: str, zarr_format: int | None
) -> None:
    source = _mudata()
    back = _read_mudata(_write_mudata(tmp_path, name, output_format, zarr_format, source, modality="rna"))
    assert list(back.mod) == ["prot", "rna"]
    assert back.axis == 0
    assert back.obs_names.tolist() == source.obs_names.tolist()
    assert back.obs["sample"].tolist() == source.obs["sample"].tolist()
    np.testing.assert_array_equal(back.obsm["X_wnn_umap"], source.obsm["X_wnn_umap"].astype(np.float32))
    for modality in ("prot", "rna"):
        mine, theirs = back.mod[modality], source.mod[modality]
        assert mine.obs_names.tolist() == theirs.obs_names.tolist()
        assert mine.var_names.tolist() == theirs.var_names.tolist()
        np.testing.assert_array_equal(_dense(mine.X), _dense(theirs.X).astype(np.float32))
        np.testing.assert_array_equal(back.obsmap[modality], source.obsmap[modality])
        np.testing.assert_array_equal(back.varmap[modality], source.varmap[modality])
        for key in theirs.obsm:
            np.testing.assert_array_equal(mine.obsm[key], np.asarray(theirs.obsm[key], dtype=np.float32))
    assert back.mod["rna"].obs["leiden"].tolist() == source.mod["rna"].obs["leiden"].tolist()
    assert list(back.mod["rna"].obsp) == ["connectivities"]
    assert back.uns["cellucid"]["defaults"] == {"modality": "rna"}
    assert back.uns["title"] == "Toy CITE-seq"
    assert "cellucid" not in back.mod["rna"].uns


def _root(path: Path) -> Any:
    if path.suffix == ".h5mu":
        import h5py

        return h5py.File(path, "r")
    if path.name.endswith(".zip"):
        return zarr.open_group(zarr.storage.ZipStore(path, mode="r"), mode="r")
    return zarr.open_group(str(path), mode="r")


@pytest.mark.parametrize(("name", "output_format", "zarr_format"), MUDATA_OUTPUTS)
def test_each_modality_is_a_cellucid_ready_anndata_under_a_mudata_root(
    tmp_path: Path, name: str, output_format: str, zarr_format: int | None
) -> None:
    path = _write_mudata(tmp_path, name, output_format, zarr_format)
    root = _root(path)
    assert root.attrs["encoding-type"] == "MuData"
    assert root.attrs["encoding-version"] == "0.1.0"
    assert int(root.attrs["axis"]) == 0
    assert [str(item) for item in root["mod"].attrs["mod-order"]] == ["prot", "rna"]
    assert sorted(root.keys()) == ["mod", "obs", "obsm", "obsmap", "obsp", "uns", "var", "varm", "varmap", "varp"]
    for modality in ("prot", "rna"):
        group = root["mod"][modality]
        assert group.attrs["encoding-type"] == "anndata"
        assert sorted(group["uns"].keys()) == []
        assert root["obsmap"][modality].dtype == np.uint32
    assert root["mod"]["rna"]["X"].attrs["encoding-type"] == "csc_matrix"
    assert root["mod"]["prot"]["X"].dtype == np.float32
    # mudata's masks, which cells each modality measured, as it keeps them.
    assert root["obsm"]["prot"].dtype == bool
    assert root["obsm"]["prot"].shape == (50, 1)
    if output_format == "h5mu":
        assert path.read_bytes()[:6] == b"MuData"
        assert path.read_bytes()[512:520] == b"\x89HDF\r\n\x1a\n"
        root.close()


def test_columns_mudata_copied_to_the_root_are_not_written_twice(tmp_path: Path) -> None:
    mdata = _mudata()
    mdata.obs["rna:leiden"] = mdata.mod["rna"].obs["leiden"].reindex(mdata.obs_names).values
    mdata.obsm["rna:X_umap"] = mdata.mod["rna"].obsm["X_umap"]
    back = _read_mudata(prepare(tmp_path / "m.h5mu", format="h5mu", data=mdata))
    assert "rna:leiden" not in back.obs.columns or back.obs["rna:leiden"].tolist() == mdata.mod["rna"].obs["leiden"].tolist()
    root = _root(tmp_path / "m.h5mu")
    assert "rna:leiden" not in root["obs"].attrs["column-order"]
    assert "rna:X_umap" not in root["obsm"]
    root.close()
    mdata.obs["batch:id"] = "x"
    with pytest.raises(ValueError, match=re.escape("obs column 'batch:id' contains a character")):
        prepare(tmp_path / "again.h5mu", format="h5mu", data=mdata)


def test_a_modality_that_measured_fewer_cells_is_placed_by_obsmap(tmp_path: Path) -> None:
    mdata = _mudata(prot_cells=10)
    back = _read_mudata(prepare(tmp_path / "m.zarr", format="zarr", data=mdata))
    obsmap = np.asarray(back.obsmap["prot"])
    assert int((obsmap == 0).sum()) == 40
    measured = back.obs_names[obsmap > 0]
    rows = obsmap[obsmap > 0] - 1
    assert back.mod["prot"].obs_names[rows].tolist() == measured.tolist()


def test_a_first_modality_that_is_not_gene_expression_is_kept_first_and_opened_by_the_viewer_rule(
    tmp_path: Path,
) -> None:
    from cellucid.data_mapping import default_modality

    mdata = _mudata(order=("prot", "rna"))
    path = prepare(tmp_path / "m.h5mu", format="h5mu", data=mdata, var_gene_id_column="symbol")
    back = _read_mudata(path)
    assert list(back.mod) == ["prot", "rna"]
    assert default_modality(list(back.mod)) == "rna"
    # The gene column is checked on the modality the viewer opens on.
    assert back.uns["cellucid"]["defaults"]["geneIdColumn"] == "symbol"
    with pytest.raises(ValueError, match=re.escape("var_gene_id_column names var['symbol'], which this store does not have.")):
        prepare(tmp_path / "prot.h5mu", format="h5mu", data=mdata, modality="prot", var_gene_id_column="symbol")


def test_a_mudata_gene_name_column_is_checked_as_an_anndata_one_is(tmp_path: Path) -> None:
    """A float column cannot name genes, in a MuData's opened modality as in an AnnData."""
    mdata = _mudata()
    mdata.mod["rna"].var["score"] = np.linspace(0, 1, mdata.mod["rna"].n_vars)
    with pytest.raises(ValueError, match=r"Modality 'rna': var_gene_id_column: .*'score'"):
        prepare(tmp_path / "m.h5mu", format="h5mu", data=mdata, var_gene_id_column="score")
    assert not (tmp_path / "m.h5mu").exists()


def test_a_modality_embedding_and_latent_space_can_be_launch_defaults(tmp_path: Path) -> None:
    mdata = _mudata()
    back = _read_mudata(
        prepare(
            tmp_path / "m.h5mu",
            format="h5mu",
            data=mdata,
            embedding_defaults={"2d": "rna:X_umap"},
            latent_key="rna:X_pca",
        )
    )
    defaults = back.uns["cellucid"]["defaults"]
    assert defaults["embedding2d"] == "rna:X_umap"
    assert defaults["latentKey"] == "rna:X_pca"
    assert "outlierQuantiles" in back.uns["cellucid"]
    with pytest.raises(ValueError, match=re.escape("latent_key names 'prot:X_umap', and modality 'prot' did not measure every cell")):
        prepare(tmp_path / "p.h5mu", format="h5mu", data=mdata, latent_key="prot:X_umap")


def _axis_one(mdata: Any) -> Any:
    rng = np.random.default_rng(1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        shared = md.MuData(
            {
                "a": ad.AnnData(rng.normal(size=(3, 2)), var=pd.DataFrame(index=["f1", "f2"])),
                "b": ad.AnnData(rng.normal(size=(4, 2)), var=pd.DataFrame(index=["f1", "f2"])),
            },
            axis=1,
        )
    return shared


def _empty_modality(mdata: Any) -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return md.MuData({"rna": mdata.mod["rna"], "atac": ad.AnnData(np.zeros((50, 0), dtype=np.float32), obs=pd.DataFrame(index=mdata.mod["rna"].obs_names))})


def _repeated_cells(mdata: Any) -> Any:
    rna = mdata.mod["rna"].copy()
    rna.obs_names = ["same"] * 2 + rna.obs_names[2:].tolist()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return md.MuData({"rna": rna})


def _no_x(mdata: Any) -> Any:
    mdata.mod["prot"].X = None
    return mdata


def _no_embedding(mdata: Any) -> Any:
    for key in list(mdata.obsm):
        del mdata.obsm[key]
    for name in mdata.mod:
        for key in list(mdata.mod[name].obsm):
            del mdata.mod[name].obsm[key]
    return mdata


@pytest.mark.parametrize(
    ("change", "output_format", "message"),
    [
        pytest.param(_axis_one, "h5mu", "This MuData shares features, not cells (axis=1); Cellucid shows cells", id="axis-1"),
        pytest.param(_empty_modality, "h5mu", "Modality 'atac' has no features (50 x 0). Remove it before prepare: del mdata.mod['atac']; mdata.update().", id="empty-modality"),
        pytest.param(_repeated_cells, "h5mu", "repeat ('same'). MuData places each modality's cells among the global ones by name", id="repeated-cell-names"),
        pytest.param(_no_x, "h5mu", "Modality 'prot' has no X. Give it the matrix the viewer reads", id="no-x"),
        pytest.param(_no_embedding, "h5mu", "This MuData holds no embedding: neither mdata.obsm nor any modality's obsm", id="no-embedding"),
        pytest.param(lambda m: m, "h5ad", 'A MuData is written as format="h5mu", "zarr" or "zarr.zip"; an .h5ad holds one AnnData.', id="h5ad"),
    ],
)
def test_a_mudata_cellucid_cannot_show_is_refused_with_what_to_do(
    tmp_path: Path, change: Any, output_format: str, message: str
) -> None:
    mdata = change(_mudata())
    name = "m.h5ad" if output_format == "h5ad" else "m.h5mu"
    with pytest.raises(ValueError, match=re.escape(message)):
        prepare(tmp_path / name, format=output_format, data=mdata)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("argument", ["gene_expression", "var", "gene_identifiers"])
def test_genes_of_one_matrix_are_refused_with_a_mudata(tmp_path: Path, argument: str) -> None:
    value: Any = {"gene_expression": np.zeros((50, 2)), "var": pd.DataFrame(index=["a", "b"]), "gene_identifiers": ["a"]}[argument]
    with pytest.raises(ValueError, match=re.escape(f"{argument} describes the genes of one matrix, and a MuData's genes are each modality's own.")):
        prepare(tmp_path / "m.h5mu", format="h5mu", data=_mudata(), **{argument: value})


def test_an_unknown_modality_is_refused(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match=re.escape("modality names 'atac', which this data does not have. Its modalities: ['prot', 'rna'].")):
        prepare(tmp_path / "m.h5mu", format="h5mu", data=_mudata(), modality="atac")


def test_zarr_format_is_refused_with_h5mu(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match=re.escape('zarr_format applies to format="zarr" or "zarr.zip", not format="h5mu".')):
        prepare(tmp_path / "m.h5mu", format="h5mu", zarr_format=2, data=_mudata())


# =============================================================================
# The web app reads what prepare() wrote
# =============================================================================

_WEB_READ = r"""
const [web, target, obsmKeysJson, obsKeysJson] = process.argv.slice(1);
const { readFileSync, readdirSync, statSync } = await import('node:fs');
const path = await import('node:path');
const { pathToFileURL } = await import('node:url');
const load = rel => import(pathToFileURL(path.join(web, rel)).href);
const { openDatasetStore } = await load('assets/js/data/anndata/open.js');
const { AnnDataReader } = await load('assets/js/data/anndata/reader.js');
console.log = () => {};
const name = path.basename(target);
let source;
if (statSync(target).isDirectory()) {
  const files = new Map();
  const walk = dir => {
    for (const entry of readdirSync(dir)) {
      const full = path.join(dir, entry);
      if (statSync(full).isDirectory()) walk(full);
      else files.set(path.relative(target, full).split(path.sep).join('/'), new File([readFileSync(full)], entry));
    }
  };
  walk(target);
  source = { files, rootName: name };
} else {
  source = { file: new File([readFileSync(target)], name) };
}
const opened = await openDatasetStore({ source }, { decodedChunkCacheBytes: 64 << 20 });
const reader = new AnnDataReader();
await reader.open(opened.root, { format: opened.readerFormat, name: opened.name, decodedChunkBytes: 64 << 20, link: opened.link });
const probe = await reader.probeMapping();
const nan = v => (Number.isNaN(v) ? null : v);
const out = { format: opened.format, nObs: reader.nObs, defaultModality: probe.defaultModality ?? null, modalities: {}, obsm: {}, obs: {} };
for (const item of probe.modalities ?? []) {
  const { index } = await reader.readVarNames(`mod/${item.name}`, null);
  const columns = await reader.getMatrixColumns(`mod/${item.name}/X`, index.map((_, i) => i));
  out.modalities[item.name] = { nCells: item.nCells, varNames: Array.from(index, String), X: columns.map(column => Array.from(column, nan)) };
}
for (const key of JSON.parse(obsmKeysJson)) out.obsm[key] = Array.from(await reader.readObsmColumns(key, [0, 1]), nan);
for (const key of JSON.parse(obsKeysJson)) {
  const field = await reader.getObsField(key);
  out.obs[key] = field.dtype === 'categorical'
    ? { categories: Array.from(field.categories, String), codes: Array.from(field.codes) }
    : { values: Array.from(field.values, nan) };
}
reader.close();
process.stdout.write(JSON.stringify(out));
"""


def _web_reader() -> str:
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is not installed: the web reader cannot run")
    if not (WEB_ROOT / "assets" / "js" / "data" / "anndata" / "open.js").is_file():
        pytest.skip(f"the web app is not checked out beside this repository ({WEB_ROOT})")
    return node


def _placed(values: np.ndarray, obsmap: np.ndarray) -> np.ndarray:
    placed = np.full((len(obsmap), *values.shape[1:]), np.nan, dtype=np.float64)
    present = obsmap > 0
    placed[present] = values[obsmap[present].astype(np.int64) - 1]
    return placed


@pytest.mark.parametrize(("name", "output_format", "zarr_format"), MUDATA_OUTPUTS)
def test_the_web_app_reads_the_prepared_mudata_as_mudata_does(
    tmp_path: Path, name: str, output_format: str, zarr_format: int | None
) -> None:
    node = _web_reader()
    path = _write_mudata(tmp_path, name, output_format, zarr_format)
    # Opened whole, the web app reads the cells as stored (priority order).
    source = _stored_mudata(path)
    completed = subprocess.run(
        [node, "--input-type=module", "-e", _WEB_READ, str(WEB_ROOT), str(path), json.dumps(["X_wnn_umap", "rna:X_umap", "prot:X_umap"]), json.dumps(["sample", "rna:leiden", "prot:isotype_ratio"])],
        capture_output=True,
        text=True,
        check=False,
        timeout=300,
    )
    assert completed.returncode == 0, completed.stderr
    read = json.loads(completed.stdout)
    assert read["format"] == ("h5mu" if output_format == "h5mu" else "mudata-zarr")
    assert read["nObs"] == source.n_obs
    assert read["defaultModality"] == "rna"
    assert list(read["modalities"]) == ["prot", "rna"]
    for modality in ("prot", "rna"):
        item = read["modalities"][modality]
        own = source.mod[modality]
        obsmap = np.asarray(source.obsmap[modality])
        assert item["nCells"] == own.n_obs
        assert item["varNames"] == own.var_names.tolist()
        expected = _placed(_dense(own.X).astype(np.float32), obsmap)
        got = np.asarray([[np.nan if v is None else v for v in column] for column in item["X"]], dtype=np.float64).T
        np.testing.assert_array_equal(got, expected)
    for key, expected in (
        ("X_wnn_umap", np.asarray(source.obsm["X_wnn_umap"], dtype=np.float64)),
        ("rna:X_umap", _placed(np.asarray(source.mod["rna"].obsm["X_umap"]), np.asarray(source.obsmap["rna"]))),
        ("prot:X_umap", _placed(np.asarray(source.mod["prot"].obsm["X_umap"]), np.asarray(source.obsmap["prot"]))),
    ):
        got = np.asarray([np.nan if v is None else v for v in read["obsm"][key]], dtype=np.float64)
        np.testing.assert_array_equal(got, expected.reshape(-1).astype(np.float32).astype(np.float64))
    assert read["obs"]["sample"]["categories"] == source.obs["sample"].cat.categories.tolist()
    leiden = source.mod["rna"].obs["leiden"]
    rows = np.asarray(source.obsmap["rna"]) - 1
    assert [read["obs"]["rna:leiden"]["categories"][code] for code in read["obs"]["rna:leiden"]["codes"]] == leiden.iloc[rows].astype(str).tolist()
    ratio = _placed(source.mod["prot"].obs["isotype_ratio"].to_numpy()[:, None], np.asarray(source.obsmap["prot"])).ravel()
    got_ratio = np.asarray([np.nan if v is None else v for v in read["obs"]["prot:isotype_ratio"]["values"]], dtype=np.float64)
    np.testing.assert_array_equal(got_ratio, ratio.astype(np.float32).astype(np.float64))


# =============================================================================
# A session applies to a MuData
# =============================================================================


def _session_mudata() -> Any:
    rna = ad.AnnData(
        X=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["c0", "c1"]),
        var=pd.DataFrame(index=["G1", "G2", "G3"]),
    )
    rna.obsm["X_umap"] = np.array([[0.0, 1.0], [1.0, 0.0]])
    prot = ad.AnnData(
        X=np.array([[7.0, 8.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["c1"]),
        var=pd.DataFrame(index=["CD3", "CD4"]),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return md.MuData({"rna": rna, "prot": prot})


def _session_mapping(**changes: Any) -> dict[str, Any]:
    mapping: dict[str, Any] = {
        "kind": "anndata",
        "embeddings": {
            "enabled": [{"key": "rna:X_umap", "columns": [0, 1]}],
            "defaults": {"1d": None, "2d": {"key": "rna:X_umap", "columns": [0, 1]}, "3d": None},
        },
        "knn": None,
        "modality": "rna",
        "geneNames": "index",
        "expression": "X",
        "log1p": {},
        "highlyVariable": None,
        "latent": None,
        "vectorFields": [],
    }
    mapping.update(changes)
    return mapping


def _session_bundle(mapping: dict[str, Any], *, var_count: int, fields: list[dict[str, Any]] | None = None) -> Any:
    from test_session_data_mapping_contract import _Bundle

    return _Bundle(
        {
            "datasetId": "exact",
            "name": None,
            "cellCount": 2,
            "varCount": var_count,
            "cellOrder": {"embedding": {"key": "rna:X_umap", "columns": [0, 1]}, "dimension": 2, "digest": "0123456789abcdef"},
            "dataMapping": mapping,
        },
        fields,
    )


def test_a_session_applies_to_a_mudata_as_the_viewer_showed_it() -> None:
    from test_session_data_mapping_contract import _var_alias

    from cellucid import apply_cellucid_session_to_anndata

    mdata = _session_mudata()
    bundle = _session_bundle(_session_mapping(), var_count=3, fields=[_var_alias("G2")])
    output = apply_cellucid_session_to_anndata(bundle, mdata, expected_dataset_id="exact", add_highlights=False)
    assert output is not mdata
    assert output.obs["copied"].tolist() == [2.0, 5.0]
    assert "copied" not in mdata.obs.columns
    assert output.uns["cellucid"]["session"]["dataset_fingerprint"]["dataMapping"]["modality"] == "rna"
    assert list(output.mod) == ["rna", "prot"]


@pytest.mark.parametrize(
    ("mapping", "var_count", "message"),
    [
        pytest.param(_session_mapping(modality="atac"), 3, "modality", id="unknown-modality"),
        pytest.param({key: value for key, value in _session_mapping().items() if key != "modality"}, 3, "modalit", id="no-modality"),
        pytest.param(_session_mapping(modality="prot"), 3, "varCount 3 != adata.n_vars 2", id="gene-count-of-the-modality"),
    ],
)
def test_a_session_that_does_not_match_the_mudata_is_refused(mapping: dict[str, Any], var_count: int, message: str) -> None:
    from cellucid import apply_cellucid_session_to_anndata

    with pytest.raises(ValueError, match=re.escape(message)):
        apply_cellucid_session_to_anndata(_session_bundle(mapping, var_count=var_count), _session_mudata(), expected_dataset_id="exact")


# =============================================================================
# The R writer's MuData Zarr
# =============================================================================


def test_the_cross_writer_harness_matches_a_store_with_itself(tmp_path: Path) -> None:
    import cross_writer_mudata as harness

    mdata = _mudata(prot_cells=12)
    spec = {
        "arguments": {"dataset_id": "cross", "dataset_name": "Cross", "modality": "rna"},
        "obs": _frame_spec(mdata.obs, mdata.obs_names),
        "obsm": {"X_wnn_umap": _matrix_spec(mdata.obsm["X_wnn_umap"])},
        "modalities": [
            {
                "name": name,
                "X": _matrix_spec(mdata.mod[name].X),
                "obs": _frame_spec(mdata.mod[name].obs, mdata.mod[name].obs_names),
                "var": _frame_spec(mdata.mod[name].var, mdata.mod[name].var_names),
                "obsm": {key: _matrix_spec(value) for key, value in mdata.mod[name].obsm.items()},
            }
            for name in mdata.mod
        ],
    }
    described = tmp_path / "input.json"
    described.write_text(json.dumps(spec), encoding="utf-8")
    for zarr_format in (2, 3):
        other = prepare(
            tmp_path / f"other-v{zarr_format}.zarr", format="zarr", zarr_format=zarr_format,
            data=harness.mudata_from_spec(spec), **spec["arguments"],
        )
        assert harness.main([str(other), str(described)]) == 0
    # A difference is named.
    changed = json.loads(described.read_text(encoding="utf-8"))
    changed["arguments"]["dataset_name"] = "Changed"
    (tmp_path / "changed.json").write_text(json.dumps(changed), encoding="utf-8")
    assert harness.compare(tmp_path / "other-v2.zarr", tmp_path / "changed.json") is not None


def _matrix_spec(value: Any) -> dict[str, Any]:
    if sparse.issparse(value):
        csc = sparse.csc_matrix(value)
        return {"kind": "sparse", "shape": list(csc.shape), "p": csc.indptr.tolist(), "i": csc.indices.tolist(), "x": csc.data.astype(np.float64).tolist()}
    array = np.asarray(value, dtype=np.float64)
    return {"kind": "dense", "type": "double", "shape": list(array.shape), "values": array.ravel().tolist()}


def _frame_spec(frame: pd.DataFrame, index: Any) -> dict[str, Any]:
    columns = []
    for name in frame.columns:
        values = frame[name]
        if isinstance(values.dtype, pd.CategoricalDtype):
            columns.append({"name": name, "type": "factor", "values": values.astype(str).tolist(), "levels": values.cat.categories.astype(str).tolist(), "ordered": bool(values.cat.ordered)})
        elif values.dtype == object:
            columns.append({"name": name, "type": "character", "values": values.tolist()})
        else:
            columns.append({"name": name, "type": "double", "values": values.astype(float).tolist()})
    return {"index": [str(name) for name in index], "columns": columns}


def test_the_r_writer_writes_the_same_mudata_zarr() -> None:
    """cellucid-r's ``cellucid_prepare(data = <multi-assay object>)`` against ``prepare(data=mdata)``.

    cellucid-r keeps, in ``tests/testthat/fixtures/cross_writer_mudata/``, each
    case as ``<case>.zarr`` (its MuData Zarr store) beside ``<case>.json`` (the
    same input, described as ``cross_writer_mudata.py`` reads it).
    """
    import cross_writer_mudata as harness

    if not (REPOSITORY.parent / "cellucid-r").is_dir():
        pytest.skip("cellucid-r is not checked out beside this repository")
    cases = sorted(R_CROSS_WRITER.glob("*.json"))
    assert cases, f"cellucid-r holds no MuData cross-writer case ({R_CROSS_WRITER}/<case>.zarr + <case>.json)"
    for described in cases:
        store = described.with_suffix(".zarr")
        assert store.is_dir(), f"{described.name} has no {store.name} beside it"
        assert harness.compare(store, described) is None, described.name


# =============================================================================
# modalities=: which modalities are written (cellucid-r's modalities =)
# =============================================================================


def test_modalities_writes_the_modalities_named_in_that_order(tmp_path: Path) -> None:
    mdata = _mudata(prot_cells=50)
    path = prepare(tmp_path / "m.zarr", format="zarr", data=mdata, modalities=["rna", "prot"])
    back = md.read_zarr(str(path))
    assert list(back.mod) == ["rna", "prot"]
    assert list(back.var_names) == [*mdata.mod["rna"].var_names, *mdata.mod["prot"].var_names]
    np.testing.assert_array_equal(back.varmap["rna"][:20], np.arange(1, 21))
    np.testing.assert_array_equal(back.varmap["prot"][20:], np.arange(1, 5))
    only = prepare(tmp_path / "rna-only.zarr", format="zarr", data=mdata, modalities=["rna"], force=True)
    assert ad.read_zarr(str(only)).var_names.tolist() == mdata.mod["rna"].var_names.tolist()


def test_one_modality_is_written_as_an_anndata_with_the_roots_cells(tmp_path: Path) -> None:
    mdata = _mudata(prot_cells=50, order=("prot", "rna"))
    path = prepare(tmp_path / "rna.h5ad", format="h5ad", data=mdata, modalities=["rna"])
    back = restore_original_order(ad.read_h5ad(path))
    assert back.obs_names.tolist() == mdata.obs_names.tolist()
    rows = np.asarray(mdata.obsmap["rna"]).astype(np.int64) - 1
    np.testing.assert_array_equal(back.X.toarray(), mdata.mod["rna"].X[rows].toarray())
    assert {"X_wnn_umap", "X_umap", "X_pca"} <= set(back.obsm)
    np.testing.assert_allclose(back.obsm["X_umap"], mdata.mod["rna"].obsm["X_umap"][rows], rtol=1e-6)
    assert "connectivities" in back.obsp
    assert back.uns["cellucid"]["id"] == "toy-citeseq"


def test_modalities_is_refused_where_it_cannot_apply(tmp_path: Path) -> None:
    mdata = _mudata(prot_cells=40)
    with pytest.raises(ValueError, match=r"Modality 'prot' measured 40 of the 50 cells"):
        prepare(tmp_path / "p.h5ad", format="h5ad", data=mdata, modalities=["prot"])
    with pytest.raises(ValueError, match=r"modalities names \['atac'\], which data does not hold"):
        prepare(tmp_path / "m.zarr", format="zarr", data=mdata, modalities=["rna", "atac"])
    with pytest.raises(ValueError, match="each once"):
        prepare(tmp_path / "m.zarr", format="zarr", data=mdata, modalities=["rna", "rna"])
    with pytest.raises(TypeError, match="sequence of modality names"):
        prepare(tmp_path / "m.zarr", format="zarr", data=mdata, modalities="rna")
    with pytest.raises(ValueError, match="an AnnData, which has one feature space"):
        prepare(tmp_path / "a.h5ad", format="h5ad", data=mdata.mod["rna"], modalities=["rna"])
    with pytest.raises(ValueError, match="modality names the modality the viewer opens on"):
        prepare(tmp_path / "r.h5ad", format="h5ad", data=mdata, modalities=["rna"], modality="rna")
    assert not any(tmp_path.iterdir())


def test_both_repositories_hold_the_same_mudata_cross_writer_harness() -> None:
    """cellucid-r runs the same harness, byte for byte, against this writer."""
    copy = REPOSITORY.parent / "cellucid-r" / "tests" / "testthat" / "python" / "cross_writer_mudata.py"
    if not copy.is_file():
        pytest.skip("cellucid-r is not checked out beside this repository")
    assert copy.read_bytes() == (Path(__file__).parent / "cross_writer_mudata.py").read_bytes()


def test_an_obs_of_the_wrong_length_names_the_mudatas_cells(tmp_path: Path) -> None:
    mdata = _mudata()
    with pytest.raises(ValueError, match=r"obs has 3 rows, and the MuData has 50 cells \(mdata.n_obs\)"):
        prepare(tmp_path / "m.zarr", format="zarr", data=mdata, obs=mdata.obs.iloc[:3])
