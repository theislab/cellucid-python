"""
The continuous-value codec.

Quantization is the half of the export format the viewer has to undo exactly,
so the constant-field case is named once here and both the writer and the
reader take that named branch instead of dividing by a zero range. The
arithmetic is chunked in float64 and nudged at integer boundaries to match what
cellucid-r produces from the same input.
"""

import numpy as np

from ..continuous_payload_diagnosis import describe_non_finite

_QUANTIZATION_ARITHMETIC_CHUNK_SIZE = 1_048_576

def _is_constant_continuous_range(min_val: float, max_val: float) -> bool:
    """Name compact_v1's constant-field case from one payload's bounds.

    A gene detected in no published cell, or an obs column a subset flattened,
    is ordinary scientific data, so the format encodes it rather than rejecting
    it. The case is declared by equal bounds -- the entry stays exactly
    ``[index, key, minValue, maxValue]``, so neither its shape nor its length
    moves -- and every code is written as ``0``. Writer and reader then take
    the same named branch instead of the general arithmetic: the writer never
    divides by ``maxValue - minValue``, and the reader fills ``minValue``
    straight through. The value returns bit-exact, not within a quantization
    step.

    Sole derivation point on the writer side, so the quantizer, the exporter
    and the tests cannot drift apart over what "constant" means.
    """
    return min_val == max_val


def _quantize_continuous(
    values: np.ndarray,
    *,
    bits: int,
    field_name: str,
) -> tuple[np.ndarray, float, float, float]:
    """
    Quantize one finite continuous float32 vector to uint8 or uint16.

    Parameters
    ----------
    values : np.ndarray
        Float32 values to quantize.
    bits : int
        Number of bits for quantization (8 or 16).
    field_name : str
        Name of field for debug messages.

    Returns
    -------
    quantized : np.ndarray
        Quantized values as uint8 or uint16.
    min_val : float
        Minimum value (for dequantization).
    max_val : float
        Maximum value (for dequantization).
    scale : float
        Scale factor (for dequantization), exactly ``0.0`` for a constant field.
    """
    if type(bits) is not int or bits not in (8, 16):
        raise ValueError("Quantization bits must be exactly 8 or 16.")

    if values.size == 0:
        raise ValueError(
            f"Continuous field {field_name!r} cannot be quantized because it has no finite values."
        )
    if not np.isfinite(values).all():
        raise ValueError(
            f"Continuous field {field_name!r} must contain only finite values. "
            + describe_non_finite(values, f"Continuous field {field_name!r}")
        )

    min_val = float(np.min(values))
    max_val = float(np.max(values))

    if bits == 8:
        max_quant = 254  # 255 is the categorical-outlier missing marker.
        dtype: type[np.uint8] | type[np.uint16] = np.uint8
    else:  # 16 bits
        max_quant = 65534  # 65535 is the categorical-outlier missing marker.
        dtype = np.uint16

    if _is_constant_continuous_range(min_val, max_val):
        # The constant case, taken before any scaling exists. The general path
        # divides by (max - min), which is exactly zero here, so a constant
        # field must never reach it. Every code is 0 and both published bounds
        # are the constant itself, which is what CONSTANT_FIELD_ENCODING
        # describes and what the reader's matching constant branch decodes.
        return np.zeros(values.shape, dtype=dtype), min_val, max_val, 0.0

    scale = max_quant / (max_val - min_val)
    quantized = np.empty(values.shape, dtype=dtype)
    for start in range(0, values.shape[0], _QUANTIZATION_ARITHMETIC_CHUNK_SIZE):
        end = start + _QUANTIZATION_ARITHMETIC_CHUNK_SIZE
        # Values already belong to the viewer's Float32 domain, but a valid
        # subnormal range can require a normalization scale beyond Float32.
        # Bounded Float64 chunks avoid both arithmetic overflow and an N-wide
        # Float64 owner.
        normalized = values[start:end].astype(np.float64, copy=True)
        normalized -= min_val
        normalized *= scale
        # Match the R exporter at integer boundaries where binary arithmetic
        # can land infinitesimally below the mathematically exact code.
        normalized += 1e-8
        np.clip(normalized, 0, max_quant, out=normalized)
        quantized[start:end] = normalized

    return quantized, min_val, max_val, scale


def _quantize_nullable_outlier_quantiles(
    values: np.ndarray,
    *,
    bits: int,
    field_name: str,
) -> tuple[np.ndarray, float, float, float]:
    """Quantize generated outlier quantiles, preserving only NaN as missing."""
    if type(bits) is not int or bits not in (8, 16):
        raise ValueError("Quantization bits must be exactly 8 or 16.")
    if np.isinf(values).any():
        raise ValueError(
            f"Outlier quantiles {field_name!r} must contain only finite values or NaN."
        )

    missing_mask = np.isnan(values)
    quantized_finite, min_val, max_val, scale = _quantize_continuous(
        values[~missing_mask],
        bits=bits,
        field_name=field_name,
    )

    if bits == 8:
        dtype: type[np.uint8] | type[np.uint16] = np.uint8
        missing_value = 255
    else:
        dtype = np.uint16
        missing_value = 65535

    quantized = np.full(values.shape, missing_value, dtype=dtype)
    quantized[~missing_mask] = quantized_finite
    return quantized, min_val, max_val, scale
