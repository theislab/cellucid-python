import numpy as np
import pandas as pd
import pytest

from cellucid.prepare_data import prepare
from cellucid.prepare_data._quantization import (
    _quantize_continuous,
    _quantize_nullable_outlier_quantiles,
)


def _prepare_with_options(out_dir, **options):
    options.setdefault("obs_categorical_dtype", "uint16")
    prepare(
        latent_space=np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            dtype=np.float32,
        ),
        obs=pd.DataFrame({"score": [0.0, 0.5, 1.0]}),
        X_umap_2d=np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            dtype=np.float32,
        ),
        out_dir=out_dir,
        dataset_id="codec-contract",
        dataset_name="Codec contract",
        centroid_min_points=1,
        force=True,
        **options,
    )


@pytest.mark.parametrize(
    ("option", "value", "message"),
    [
        ("compression", 0, "compression"),
        ("compression", -1, "compression"),
        ("compression", 10, "compression"),
        ("compression", 1.0, "compression"),
        ("compression", True, "compression"),
        ("var_quantization", 7, "var_quantization"),
        ("var_quantization", 8.0, "var_quantization"),
        ("var_quantization", True, "var_quantization"),
        (
            "obs_continuous_quantization",
            7,
            "obs_continuous_quantization",
        ),
        (
            "obs_continuous_quantization",
            16.0,
            "obs_continuous_quantization",
        ),
        ("obs_categorical_dtype", "garbage", "obs_categorical_dtype"),
        ("obs_categorical_dtype", 8, "obs_categorical_dtype"),
        ("centroid_outlier_quantile", 0.5, "centroid_outlier_quantile"),
        ("centroid_outlier_quantile", 1.0, "centroid_outlier_quantile"),
        (
            "centroid_outlier_quantile",
            float("nan"),
            "centroid_outlier_quantile",
        ),
        ("centroid_outlier_quantile", True, "centroid_outlier_quantile"),
    ],
)
def test_invalid_codec_options_reject_before_output_mutation(
    tmp_path,
    option,
    value,
    message,
):
    out_dir = tmp_path / "must-not-exist"

    with pytest.raises(ValueError, match=message):
        _prepare_with_options(out_dir, **{option: value})

    assert not out_dir.exists()


@pytest.mark.parametrize("bits", [0, 7, 8.0, 32, True])
def test_quantizer_requires_exact_current_bit_width(bits):
    with pytest.raises(ValueError, match="exactly 8 or 16"):
        _quantize_continuous(
            np.array([0.0, 1.0], dtype=np.float32),
            bits=bits,
            field_name="score",
        )


@pytest.mark.parametrize(
    ("values", "message"),
    [
        (
            np.array([], dtype=np.float32),
            "no finite values",
        ),
        (
            np.array([2.0, np.nan, 4.0], dtype=np.float32),
            "only finite values",
        ),
        (
            np.array([2.0, np.inf, 4.0], dtype=np.float32),
            "only finite values",
        ),
    ],
)
def test_quantizer_rejects_undefined_scientific_ranges(values, message):
    with pytest.raises(ValueError, match=message):
        _quantize_continuous(values, bits=8, field_name="score")


@pytest.mark.parametrize("bits", [8, 16])
@pytest.mark.parametrize("constant", [0.0, -0.0, 3.0, -12.5, 2.0**-149, 3.4028235e38])
def test_quantizer_encodes_a_constant_field_as_equal_bounds_and_zero_codes(bits, constant):
    values = np.full(7, constant, dtype=np.float32)

    quantized, minimum, maximum, scale = _quantize_continuous(
        values,
        bits=bits,
        field_name="constant",
    )

    assert quantized.dtype == (np.uint8 if bits == 8 else np.uint16)
    assert quantized.tolist() == [0] * 7
    assert minimum == float(np.float32(constant))
    assert maximum == minimum
    # No scale exists for a constant field, and the writer never derives one.
    assert scale == 0.0


@pytest.mark.parametrize("bits", [8, 16])
@pytest.mark.parametrize("constant", [0.0, 3.0, -12.5, 2.0**-149, 3.4028235e38])
def test_the_reader_arithmetic_recovers_a_constant_field_exactly(bits, constant):
    """The published bounds and codes must decode back to the exact constant.

    Not within a quantization step: the same value, bit for bit, under both the
    reader's constant branch and the general dequantization formula.
    """
    max_quant = 254 if bits == 8 else 65534
    quantized, minimum, maximum, _ = _quantize_continuous(
        np.full(5, constant, dtype=np.float32),
        bits=bits,
        field_name="constant",
    )

    # The reader's explicit constant branch.
    assert np.float32(minimum) == np.float32(constant)
    # The general formula, which must not divide by (maximum - minimum).
    decoded = (
        minimum + quantized.astype(np.float64) * ((maximum - minimum) / max_quant)
    ).astype(np.float32)
    np.testing.assert_array_equal(decoded, np.full(5, constant, dtype=np.float32))
    assert np.isfinite(decoded).all()


def test_quantizer_preserves_a_noncollapsed_float32_subnormal_range():
    source_values = np.array([2.0**-150, 2.0**-149, 2.0**-148])
    visible_values = source_values.astype(np.float32)

    assert visible_values.tolist() == [0.0, 2.0**-149, 2.0**-148]
    with np.errstate(over="raise", invalid="raise"):
        quantized, minimum, maximum, scale = _quantize_continuous(
            visible_values,
            bits=8,
            field_name="subnormal",
        )

    assert quantized.tolist() == [0, 127, 254]
    assert minimum == 0.0
    assert maximum == 2.0**-148
    assert np.isfinite(scale)
    decoded = (minimum + quantized.astype(np.float64) / 254 * (maximum - minimum)).astype(
        np.float32
    )
    np.testing.assert_array_equal(decoded, visible_values)


@pytest.mark.parametrize(
    ("values", "bits", "expected"),
    [
        ([-99633.28125, 99994.2109375], 8, [0, 254]),
        ([-99922.4375, 99958.34375], 16, [0, 65534]),
    ],
)
def test_quantizer_maps_exact_float32_bounds_to_terminal_codes(values, bits, expected):
    quantized, *_ = _quantize_continuous(
        np.array(values, dtype=np.float32),
        bits=bits,
        field_name="bounds",
    )

    assert quantized.tolist() == expected


def test_outlier_quantizer_preserves_nan_with_an_exact_reserved_marker():
    quantized, minimum, maximum, scale = _quantize_nullable_outlier_quantiles(
        np.array([2.0, np.nan, 4.0], dtype=np.float32),
        bits=8,
        field_name="cluster_outliers",
    )

    assert quantized.tolist() == [0, 255, 254]
    assert minimum == 2.0
    assert maximum == 4.0
    assert scale == 127.0


@pytest.mark.parametrize("invalid_value", [np.inf, -np.inf])
def test_outlier_quantizer_rejects_infinity(invalid_value):
    with pytest.raises(ValueError, match="finite values or NaN"):
        _quantize_nullable_outlier_quantiles(
            np.array([0.25, invalid_value, 0.75], dtype=np.float32),
            bits=8,
            field_name="cluster_outliers",
        )


def test_prepare_reserves_missing_marker_for_generated_outlier_quantiles(tmp_path):
    out_dir = tmp_path / "outlier-missing-marker"
    embedding = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        dtype=np.float32,
    )
    prepare(
        latent_space=np.array(
            [[0.0, 0.0], [1.0, 0.0], [4.0, 0.0], [8.0, 0.0]],
            dtype=np.float32,
        ),
        obs=pd.DataFrame({"cluster": pd.Categorical(["large", "large", "large", "small"])}),
        X_umap_2d=embedding,
        out_dir=out_dir,
        dataset_id="outlier-missing-marker",
        dataset_name="Outlier missing marker",
        centroid_min_points=2,
        obs_continuous_quantization=8,
        obs_categorical_dtype="uint16",
    )

    outlier_codes = np.fromfile(out_dir / "obs/0.outliers.u8", dtype=np.uint8)
    assert outlier_codes[:3].tolist() == [127, 0, 254]
    assert outlier_codes[3] == 255
