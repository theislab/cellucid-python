"""compact_v1 quantization requires variation, and prepare() must say so up front.

A quantized field is published as ``[key, minValue, maxValue]`` and the viewer
discards a whole manifest whose ``minValue`` is not strictly below its
``maxValue``, so a payload with no variation has no quantized encoding at all.
Subsetting an expression matrix to one lineage routinely leaves genes detected
in no remaining cell, so the rejection has to name every offending payload
before any file is written rather than one per failed export.
"""

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import cellucid.prepare_data as prepare_data
from cellucid import prepare

LINEAGE_CELLS = 60
OTHER_CELLS = 40
N_GENES = 64
# Genes carried by the discarded lineage only; zero everywhere in the subset.
LINEAGE_SILENT_GENES = ("GENE07", "GENE19", "GENE41", "GENE58")


def _gene_ids() -> list[str]:
    return [f"GENE{index:02d}" for index in range(N_GENES)]


def _subset_to_lineage(sparse_expression: bool = False):
    """Build a haematopoietic-style lineage subset of a larger expression matrix."""
    rng = np.random.default_rng(20260730)
    n_cells = LINEAGE_CELLS + OTHER_CELLS
    gene_ids = _gene_ids()

    expression = rng.random((n_cells, N_GENES), dtype=np.float32) + 0.5
    lineage_mask = np.zeros(n_cells, dtype=bool)
    lineage_mask[:LINEAGE_CELLS] = True
    for gene_id in LINEAGE_SILENT_GENES:
        # Detected only in the cells that the subset discards.
        expression[lineage_mask, gene_ids.index(gene_id)] = 0.0

    subset = expression[lineage_mask, :]
    if sparse_expression:
        subset = sparse.csr_matrix(subset)

    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(["HSC"] * 30 + ["MPP"] * 30),
            "total_counts": rng.random(LINEAGE_CELLS).astype(np.float32) + 1.0,
        }
    )
    latent = rng.random((LINEAGE_CELLS, 8)).astype(np.float32)
    embedding = rng.random((LINEAGE_CELLS, 2)).astype(np.float32)
    var = pd.DataFrame(index=pd.Index(gene_ids, name="gene_id"))
    return subset, obs, var, latent, embedding


def _prepare_lineage_subset(out_dir, *, sparse_expression=False, obs=None, latent=None, **options):
    expression, default_obs, var, default_latent, embedding = _subset_to_lineage(
        sparse_expression=sparse_expression
    )
    resolved_obs = default_obs if obs is None else obs
    options.setdefault("var_quantization", 8)
    options.setdefault("obs_continuous_quantization", 8)
    options.setdefault("centroid_min_points", 10)
    prepare(
        latent_space=default_latent if latent is None else latent,
        obs=resolved_obs,
        var=var,
        gene_expression=expression,
        X_umap_2d=embedding,
        out_dir=out_dir,
        obs_keys=list(resolved_obs.columns),
        dataset_id="lineage-subset",
        dataset_name="Lineage subset",
        obs_categorical_dtype="uint8",
        compression=6,
        force=True,
        **options,
    )


def test_every_zero_variance_gene_of_a_lineage_subset_is_named_in_one_failure(tmp_path):
    out_dir = tmp_path / "lineage-subset"

    with pytest.raises(ValueError) as failure:
        _prepare_lineage_subset(out_dir)

    message = str(failure.value)
    assert f"{len(LINEAGE_SILENT_GENES)} gene(s) with var_quantization=8" in message
    for gene_id in LINEAGE_SILENT_GENES:
        assert f"{gene_id!r}: every value is the constant 0.0" in message
    # Every other gene still has variation and must not be blamed.
    for gene_id in _gene_ids():
        if gene_id not in LINEAGE_SILENT_GENES:
            assert f"'{gene_id}':" not in message
    assert "gene_identifiers=" in message
    assert "var_quantization=None" in message
    assert not out_dir.exists()


def test_zero_variance_genes_are_reported_before_any_file_is_written(tmp_path, monkeypatch):
    def _forbidden_write(*args, **kwargs):
        raise AssertionError("prepare() wrote an output payload before validating quantizability")

    monkeypatch.setattr(prepare_data, "_write_binary", _forbidden_write)

    with pytest.raises(ValueError, match="cannot be quantized"):
        _prepare_lineage_subset(tmp_path / "lineage-subset")


def test_sparse_lineage_subset_reports_the_same_zero_variance_genes(tmp_path):
    with pytest.raises(ValueError) as failure:
        _prepare_lineage_subset(tmp_path / "lineage-subset", sparse_expression=True)

    message = str(failure.value)
    for gene_id in LINEAGE_SILENT_GENES:
        assert f"{gene_id!r}: every value is the constant 0.0" in message


def test_dropping_the_zero_variance_genes_publishes_the_remaining_feature_set(tmp_path):
    out_dir = tmp_path / "lineage-subset"
    retained = [gene_id for gene_id in _gene_ids() if gene_id not in LINEAGE_SILENT_GENES]

    _prepare_lineage_subset(out_dir, gene_identifiers=retained)

    manifest = pd.read_json(out_dir / "var_manifest.json", typ="series")
    exported = [field[0] for field in manifest["fields"]]
    assert exported == retained
    for field in manifest["fields"]:
        assert field[1] < field[2]


def test_full_precision_export_keeps_every_zero_variance_gene(tmp_path):
    out_dir = tmp_path / "lineage-subset"

    _prepare_lineage_subset(out_dir, var_quantization=None)

    manifest = pd.read_json(out_dir / "var_manifest.json", typ="series")
    exported = [field[0] for field in manifest["fields"]]
    assert exported == _gene_ids()
    assert all(len(field) == 1 for field in manifest["fields"])


def test_constant_continuous_obs_field_is_reported_with_the_genes(tmp_path):
    _, obs, _, _, _ = _subset_to_lineage()
    obs = obs.assign(total_counts=np.full(LINEAGE_CELLS, 1.5, dtype=np.float32))

    with pytest.raises(ValueError) as failure:
        _prepare_lineage_subset(tmp_path / "lineage-subset", obs=obs)

    message = str(failure.value)
    assert "1 continuous obs field(s) with obs_continuous_quantization=8" in message
    assert "'total_counts': every value is the constant 1.5" in message
    assert "obs_keys=" in message
    assert "obs_continuous_quantization=None" in message
    # The genes are reported in the same failure, not on a later run.
    assert f"{len(LINEAGE_SILENT_GENES)} gene(s) with var_quantization=8" in message
    assert str(len(LINEAGE_SILENT_GENES) + 1) in message.split(" continuous payload(s)")[0]


def test_constant_continuous_obs_field_survives_full_precision_export(tmp_path):
    out_dir = tmp_path / "lineage-subset"
    _, obs, _, _, _ = _subset_to_lineage()
    obs = obs.assign(total_counts=np.full(LINEAGE_CELLS, 1.5, dtype=np.float32))

    _prepare_lineage_subset(
        out_dir,
        obs=obs,
        var_quantization=None,
        obs_continuous_quantization=None,
    )

    manifest = pd.read_json(out_dir / "obs_manifest.json", typ="series")
    assert [field[0] for field in manifest["_continuousFields"]] == ["total_counts"]
    assert all(len(field) == 1 for field in manifest["_continuousFields"])


def test_generated_outlier_quantiles_without_a_qualifying_category_are_reported(tmp_path):
    _, obs, _, _, _ = _subset_to_lineage()
    obs = obs.assign(cell_type=pd.Categorical([f"C{index // 2}" for index in range(LINEAGE_CELLS)]))

    with pytest.raises(ValueError) as failure:
        _prepare_lineage_subset(
            tmp_path / "lineage-subset",
            obs=obs,
            gene_identifiers=[
                gene_id for gene_id in _gene_ids() if gene_id not in LINEAGE_SILENT_GENES
            ],
        )

    message = str(failure.value)
    assert "generated categorical outlier quantile set(s)" in message
    assert "'cell_type': no category holds at least centroid_min_points=10 cells" in message
    assert "centroid_min_points" in message


def test_generated_outlier_quantiles_without_latent_variation_are_reported(tmp_path):
    _, obs, _, _, _ = _subset_to_lineage()
    obs = obs.assign(cell_type=pd.Categorical(["HSC"] * LINEAGE_CELLS))

    with pytest.raises(ValueError) as failure:
        _prepare_lineage_subset(
            tmp_path / "lineage-subset",
            obs=obs,
            latent=np.zeros((LINEAGE_CELLS, 8), dtype=np.float32),
            gene_identifiers=[
                gene_id for gene_id in _gene_ids() if gene_id not in LINEAGE_SILENT_GENES
            ],
        )

    message = str(failure.value)
    assert "'cell_type': every value is the constant 1.0" in message


def test_float64_variation_that_collapses_under_float32_is_reported_as_constant(tmp_path):
    expression, obs, var, latent, embedding = _subset_to_lineage()
    expression = expression.astype(np.float64)
    collapsing = np.full(LINEAGE_CELLS, 1.0, dtype=np.float64)
    collapsing[0] = 1.0 + 1e-17
    expression[:, _gene_ids().index("GENE03")] = collapsing

    with pytest.raises(ValueError) as failure:
        prepare(
            latent_space=latent,
            obs=obs,
            var=var,
            gene_expression=expression,
            X_umap_2d=embedding,
            out_dir=tmp_path / "lineage-subset",
            obs_keys=list(obs.columns),
            dataset_id="lineage-subset",
            dataset_name="Lineage subset",
            obs_categorical_dtype="uint8",
            var_quantization=8,
            obs_continuous_quantization=8,
            centroid_min_points=10,
            compression=6,
            force=True,
        )

    assert "'GENE03': every value is the constant 1.0" in str(failure.value)


def test_nonfinite_gene_values_are_rejected_before_any_file_is_written(tmp_path, monkeypatch):
    expression, obs, var, latent, embedding = _subset_to_lineage()
    expression = expression.copy()
    expression[3, _gene_ids().index("GENE11")] = np.nan

    def _forbidden_write(*args, **kwargs):
        raise AssertionError("prepare() wrote an output payload before validating gene values")

    monkeypatch.setattr(prepare_data, "_write_binary", _forbidden_write)

    with pytest.raises(ValueError, match="Gene 'GENE11' expression must contain only finite values"):
        prepare(
            latent_space=latent,
            obs=obs,
            var=var,
            gene_expression=expression,
            X_umap_2d=embedding,
            out_dir=tmp_path / "lineage-subset",
            obs_keys=list(obs.columns),
            dataset_id="lineage-subset",
            dataset_name="Lineage subset",
            obs_categorical_dtype="uint8",
            var_quantization=8,
            obs_continuous_quantization=8,
            centroid_min_points=10,
            compression=6,
            force=True,
        )
