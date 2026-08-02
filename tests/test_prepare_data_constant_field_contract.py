"""compact_v1 has a case for a constant continuous field, and prepare() uses it.

A gene expressed at one level in every published cell -- very often zero, once
an atlas is subset to one lineage -- is ordinary scientific data, and so is an
obs column a subset flattened. compact_v1 publishes such a field as equal
bounds with every code ``0``, so dequantization returns the exact constant
instead of an approximation, and nothing anywhere divides by
``maxValue - minValue``.

Every assertion below reads the payload back off disk and decodes it with the
viewer's own arithmetic (``cellucid/assets/js/data/data-loaders.js``), so a
writer that emitted a plausible-looking manifest and a wrong payload would fail
here.
"""

import gzip
import json

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid import prepare
from cellucid.prepare_data import _generation as prepare_generation
from cellucid.prepare_data import _quantization as prepare_quantization

LINEAGE_CELLS = 60
OTHER_CELLS = 40
N_GENES = 64
# Genes carried by the discarded lineage only; zero everywhere in the subset.
LINEAGE_SILENT_GENES = ("GENE07", "GENE19", "GENE41", "GENE58")
# A gene detected at one identical nonzero level in every published cell.
CONSTANT_NONZERO_GENE = "GENE22"
CONSTANT_NONZERO_LEVEL = 2.5

MAX_QUANT = {8: 254, 16: 65534}
# The dtypes a reader decodes a payload with: uint8 codes are single bytes and
# have no byte order, uint16 codes carry the format's published little-endian
# order, which is not the host's on every architecture.
DTYPE = {8: np.dtype(np.uint8), 16: np.dtype("<u2")}


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
    expression[lineage_mask, gene_ids.index(CONSTANT_NONZERO_GENE)] = CONSTANT_NONZERO_LEVEL

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


def _read_payload(out_dir, path_pattern, payload_index):
    relative = path_pattern.replace("{index}", str(payload_index))
    raw = (out_dir / relative).read_bytes()
    if relative.endswith(".gz"):
        raw = gzip.decompress(raw)
    return raw


def _viewer_dequantize(codes, minimum, maximum, bits):
    """Decode exactly as ``dequantize()`` in the viewer's data-loaders.js does.

    JavaScript evaluates the expression in float64 and stores it into a
    Float32Array, so the intermediate stays float64 and only the result narrows.
    """
    scale = (maximum - minimum) / MAX_QUANT[bits]
    return (minimum + codes.astype(np.float64) * scale).astype(np.float32)


def _read_var_field(out_dir, gene_id):
    manifest = json.loads((out_dir / "var_manifest.json").read_text(encoding="utf-8"))
    schema = manifest["_varSchema"]
    for entry in manifest["fields"]:
        if entry[1] == gene_id:
            return manifest, schema, entry
    raise AssertionError(f"{gene_id!r} is absent from var_manifest.json")


def _read_obs_continuous_field(out_dir, key):
    manifest = json.loads((out_dir / "obs_manifest.json").read_text(encoding="utf-8"))
    schema = manifest["_obsSchemas"]["continuous"]
    for entry in manifest["_continuousFields"]:
        if entry[1] == key:
            return manifest, schema, entry
    raise AssertionError(f"{key!r} is absent from obs_manifest.json")


def _assert_constant_field_round_trips(out_dir, schema, entry, *, constant, bits, n_values):
    """The published entry and payload must decode back to the exact constant."""
    payload_index, _key, minimum, maximum = entry

    # The manifest declares the constant case: equal bounds, both the constant.
    assert minimum == maximum, "a constant field must publish equal bounds"
    assert np.float32(minimum) == np.float32(constant)
    assert prepare_quantization._is_constant_continuous_range(minimum, maximum)

    codes = np.frombuffer(
        _read_payload(out_dir, schema["pathPattern"], payload_index),
        dtype=DTYPE[bits],
    )
    assert codes.shape == (n_values,)
    # Every code is 0: the writer never derived a scale for this field.
    assert codes.tolist() == [0] * n_values

    # The general dequantization arithmetic stays exact and finite: it divides
    # by maxQuant, never by (maxValue - minValue).
    decoded = _viewer_dequantize(codes, minimum, maximum, bits)
    expected = np.full(n_values, constant, dtype=np.float32)
    np.testing.assert_array_equal(decoded, expected)
    assert np.isfinite(decoded).all()


# ---------------------------------------------------------------------------
# Genes
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("bits", [8, 16])
def test_a_lineage_subset_publishes_its_undetected_genes_as_constants(tmp_path, bits):
    out_dir = tmp_path / "lineage-subset"

    _prepare_lineage_subset(out_dir, var_quantization=bits, obs_continuous_quantization=bits)

    manifest, schema, _ = _read_var_field(out_dir, LINEAGE_SILENT_GENES[0])
    assert [entry[1] for entry in manifest["fields"]] == _gene_ids()
    for gene_id in LINEAGE_SILENT_GENES:
        _, _, entry = _read_var_field(out_dir, gene_id)
        _assert_constant_field_round_trips(
            out_dir,
            schema,
            entry,
            constant=0.0,
            bits=bits,
            n_values=LINEAGE_CELLS,
        )


@pytest.mark.parametrize("bits", [8, 16])
def test_a_constant_nonzero_gene_is_published_and_recovered_exactly(tmp_path, bits):
    out_dir = tmp_path / "lineage-subset"

    _prepare_lineage_subset(out_dir, var_quantization=bits, obs_continuous_quantization=bits)

    _, schema, entry = _read_var_field(out_dir, CONSTANT_NONZERO_GENE)
    _assert_constant_field_round_trips(
        out_dir,
        schema,
        entry,
        constant=CONSTANT_NONZERO_LEVEL,
        bits=bits,
        n_values=LINEAGE_CELLS,
    )


@pytest.mark.parametrize("bits", [8, 16])
def test_a_varying_gene_beside_the_constants_keeps_its_own_bounds(tmp_path, bits):
    out_dir = tmp_path / "lineage-subset"

    _prepare_lineage_subset(out_dir, var_quantization=bits, obs_continuous_quantization=bits)

    varying = [
        gene_id
        for gene_id in _gene_ids()
        if gene_id not in LINEAGE_SILENT_GENES and gene_id != CONSTANT_NONZERO_GENE
    ]
    for gene_id in varying:
        _, schema, entry = _read_var_field(out_dir, gene_id)
        _payload_index, _key, minimum, maximum = entry
        assert minimum < maximum
        codes = np.frombuffer(
            _read_payload(out_dir, schema["pathPattern"], entry[0]),
            dtype=DTYPE[bits],
        )
        # A field the writer took the general path for reaches both terminals.
        assert codes.min() == 0
        assert codes.max() == MAX_QUANT[bits]


def test_full_precision_export_recovers_the_constant_genes_exactly(tmp_path):
    out_dir = tmp_path / "lineage-subset"

    _prepare_lineage_subset(out_dir, var_quantization=None)

    manifest, schema, _ = _read_var_field(out_dir, LINEAGE_SILENT_GENES[0])
    assert all(len(field) == 2 for field in manifest["fields"])
    for gene_id, constant in (
        (LINEAGE_SILENT_GENES[0], 0.0),
        (CONSTANT_NONZERO_GENE, CONSTANT_NONZERO_LEVEL),
    ):
        _, _, entry = _read_var_field(out_dir, gene_id)
        values = np.frombuffer(
            _read_payload(out_dir, schema["pathPattern"], entry[0]),
            dtype="<f4",
        )
        np.testing.assert_array_equal(values, np.full(LINEAGE_CELLS, constant, dtype=np.float32))


@pytest.mark.parametrize("bits", [8, 16])
def test_a_sparse_lineage_subset_publishes_the_same_constant_genes(tmp_path, bits):
    out_dir = tmp_path / "lineage-subset"

    _prepare_lineage_subset(
        out_dir,
        sparse_expression=True,
        var_quantization=bits,
        obs_continuous_quantization=bits,
    )

    for gene_id, constant in (
        (LINEAGE_SILENT_GENES[0], 0.0),
        (CONSTANT_NONZERO_GENE, CONSTANT_NONZERO_LEVEL),
    ):
        _, schema, entry = _read_var_field(out_dir, gene_id)
        _assert_constant_field_round_trips(
            out_dir,
            schema,
            entry,
            constant=constant,
            bits=bits,
            n_values=LINEAGE_CELLS,
        )


def test_float64_variation_that_collapses_under_float32_is_published_as_a_constant(tmp_path):
    """A range that exists only in float64 is one float32 value, which is a constant."""
    out_dir = tmp_path / "lineage-subset"
    expression, obs, var, latent, embedding = _subset_to_lineage()
    expression = expression.astype(np.float64)
    collapsing = np.full(LINEAGE_CELLS, 1.0, dtype=np.float64)
    collapsing[0] = 1.0 + 1e-17
    expression[:, _gene_ids().index("GENE03")] = collapsing

    prepare(
        latent_space=latent,
        obs=obs,
        var=var,
        gene_expression=expression,
        X_umap_2d=embedding,
        out_dir=out_dir,
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

    _, schema, entry = _read_var_field(out_dir, "GENE03")
    _assert_constant_field_round_trips(
        out_dir,
        schema,
        entry,
        constant=1.0,
        bits=8,
        n_values=LINEAGE_CELLS,
    )


# ---------------------------------------------------------------------------
# Continuous obs fields
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("bits", [8, 16])
def test_a_constant_continuous_obs_field_is_published_and_recovered_exactly(tmp_path, bits):
    out_dir = tmp_path / "lineage-subset"
    _, obs, _, _, _ = _subset_to_lineage()
    obs = obs.assign(total_counts=np.full(LINEAGE_CELLS, 1.5, dtype=np.float32))

    _prepare_lineage_subset(
        out_dir,
        obs=obs,
        var_quantization=bits,
        obs_continuous_quantization=bits,
    )

    _, schema, entry = _read_obs_continuous_field(out_dir, "total_counts")
    _assert_constant_field_round_trips(
        out_dir,
        schema,
        entry,
        constant=1.5,
        bits=bits,
        n_values=LINEAGE_CELLS,
    )


def test_full_precision_export_recovers_the_constant_obs_field_exactly(tmp_path):
    out_dir = tmp_path / "lineage-subset"
    _, obs, _, _, _ = _subset_to_lineage()
    obs = obs.assign(total_counts=np.full(LINEAGE_CELLS, 1.5, dtype=np.float32))

    _prepare_lineage_subset(
        out_dir,
        obs=obs,
        var_quantization=None,
        obs_continuous_quantization=None,
    )

    manifest, schema, entry = _read_obs_continuous_field(out_dir, "total_counts")
    assert all(len(field) == 2 for field in manifest["_continuousFields"])
    values = np.frombuffer(
        _read_payload(out_dir, schema["pathPattern"], entry[0]),
        dtype="<f4",
    )
    np.testing.assert_array_equal(values, np.full(LINEAGE_CELLS, 1.5, dtype=np.float32))


# ---------------------------------------------------------------------------
# Generated categorical outlier quantiles
# ---------------------------------------------------------------------------


def test_constant_generated_outlier_quantiles_are_published(tmp_path):
    """One qualifying category yields one quantile value, which is a constant."""
    out_dir = tmp_path / "lineage-subset"
    _, obs, _, _, _ = _subset_to_lineage()
    obs = obs.assign(cell_type=pd.Categorical(["HSC"] * LINEAGE_CELLS))

    _prepare_lineage_subset(
        out_dir,
        obs=obs,
        latent=np.zeros((LINEAGE_CELLS, 8), dtype=np.float32),
    )

    manifest = json.loads((out_dir / "obs_manifest.json").read_text(encoding="utf-8"))
    schema = manifest["_obsSchemas"]["categorical"]
    entry = manifest["_categoricalFields"][0]
    assert entry[1] == "cell_type"
    assert len(entry) == 8
    minimum, maximum = entry[6], entry[7]
    assert minimum == maximum == 1.0
    codes = np.frombuffer(
        _read_payload(out_dir, schema["outlierPathPattern"], entry[0]),
        dtype=np.uint8,
    )
    assert codes.tolist() == [0] * LINEAGE_CELLS
    np.testing.assert_array_equal(
        _viewer_dequantize(codes, minimum, maximum, 8),
        np.full(LINEAGE_CELLS, 1.0, dtype=np.float32),
    )


def test_generated_outlier_quantiles_without_a_qualifying_category_are_reported(tmp_path):
    """A set with no quantile at all has no value to publish, unlike a constant."""
    _, obs, _, _, _ = _subset_to_lineage()
    obs = obs.assign(cell_type=pd.Categorical([f"C{index // 2}" for index in range(LINEAGE_CELLS)]))

    with pytest.raises(ValueError) as failure:
        _prepare_lineage_subset(tmp_path / "lineage-subset", obs=obs)

    message = str(failure.value)
    assert "generated categorical outlier quantile set(s) cannot be encoded" in message
    assert "'cell_type': no category holds at least centroid_min_points=10 cells" in message
    assert "obs_continuous_quantization=None" in message


def test_unencodable_outlier_quantiles_are_reported_before_any_file_is_written(
    tmp_path, monkeypatch
):
    _, obs, _, _, _ = _subset_to_lineage()
    obs = obs.assign(cell_type=pd.Categorical([f"C{index // 2}" for index in range(LINEAGE_CELLS)]))

    def _forbidden_write(*args, **kwargs):
        raise AssertionError("prepare() wrote an output payload before validating encodability")

    monkeypatch.setattr(prepare_generation, "_write_binary", _forbidden_write)

    with pytest.raises(ValueError, match="cannot be encoded"):
        _prepare_lineage_subset(tmp_path / "lineage-subset", obs=obs)


# ---------------------------------------------------------------------------
# Rejections that survive
# ---------------------------------------------------------------------------


def test_nonfinite_gene_values_are_rejected_before_any_file_is_written(tmp_path, monkeypatch):
    expression, obs, var, latent, embedding = _subset_to_lineage()
    expression = expression.copy()
    expression[3, _gene_ids().index("GENE11")] = np.nan

    def _forbidden_write(*args, **kwargs):
        raise AssertionError("prepare() wrote an output payload before validating gene values")

    monkeypatch.setattr(prepare_generation, "_write_binary", _forbidden_write)

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
