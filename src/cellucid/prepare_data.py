#!/usr/bin/env python3
"""
Export raw dataframes/arrays to files used by the WebGL viewer.

Includes memory/disk optimization features:
- Quantization for continuous data (var/gene expression and obs continuous)
- Auto dtype selection for categorical obs based on category count
- Gzip compression for all binary files

Instead of AnnData, accepts:
- latent_space: (n_cells, n_dims) numpy/sparse array for outlier quantile calculation
- X_umap_1d / X_umap_2d / X_umap_3d: explicit embeddings (at least one required)
- vector_fields: dict[str, array] of per-cell displacement vectors (optional)
- obs: pandas DataFrame with cell metadata columns
- var: pandas DataFrame with gene/feature metadata
- gene_expression: (n_cells, n_genes) numpy/sparse array for gene expression matrix
- var_gene_id_column: exact column name in var, or None to use var.index
- connectivities: sparse matrix with KNN connectivities
"""

import gzip
import json
import math
import os
import re
import shutil
import tempfile
from collections.abc import Iterable, Iterator, Sequence
from contextlib import contextmanager
from datetime import datetime, timezone
from numbers import Integral, Real
from pathlib import Path
from typing import Any, Literal, cast

import numpy as np
import pandas as pd
import tqdm
from scipy import sparse

from .connectivity_contract import (
    CONNECTIVITY_BINARY_DIRNAME,
    CONNECTIVITY_MANIFEST_FILENAME,
    ConnectivityEdgePairs,
    build_connectivity_manifest,
    validate_connectivity_edges,
)
from .vector_fields import scale_vector_field, validate_vector_fields

DEFAULT_EXPORT_DIR = Path.cwd() / "exports"
DEFAULT_OBS_DIRNAME = "obs"
DEFAULT_VAR_DIRNAME = "var"

# Manifest format version for compact format
MANIFEST_FORMAT_VERSION = "compact_v1"
JsonScalar = str | bool | int | float
_PORTABLE_COMPONENT_PATTERN = re.compile(
    r"^[A-Za-z0-9][A-Za-z0-9._-]{0,179}$",
    flags=re.ASCII,
)
_WINDOWS_RESERVED_COMPONENTS = {
    "CON",
    "PRN",
    "AUX",
    "NUL",
    *(f"COM{index}" for index in range(1, 10)),
    *(f"LPT{index}" for index in range(1, 10)),
}


def _require_portable_filename_component(
    name: object,
    *,
    label: str = "Field key",
) -> str:
    """Require one exact cross-platform ASCII filename component."""
    if not isinstance(name, str):
        raise TypeError(f"{label} must be a native string.")
    if not _PORTABLE_COMPONENT_PATTERN.fullmatch(name):
        raise ValueError(
            f"{label} {name!r} must be 1-180 ASCII bytes, start with an "
            "ASCII letter or digit, and otherwise contain only ASCII "
            "letters, digits, '.', '_', or '-'."
        )
    if name.endswith("."):
        raise ValueError(f"{label} {name!r} must not end with '.'.")
    windows_stem = name.split(".", 1)[0].upper()
    if windows_stem in _WINDOWS_RESERVED_COMPONENTS:
        raise ValueError(f"{label} {name!r} is a reserved Windows filename.")
    return name


def _assert_unique_filename_components(
    keys: list[str],
    *,
    label: str,
) -> list[str]:
    raw_keys: set[str] = set()
    portable_to_raw: dict[str, str] = {}
    portable_keys: list[str] = []
    for key in keys:
        portable_key = _require_portable_filename_component(
            key,
            label=f"{label} key",
        )
        if key in raw_keys:
            raise ValueError(f"{label} key {key!r} is duplicated.")
        collision_key = portable_key.casefold()
        if collision_key in portable_to_raw:
            raise ValueError(
                f"{label} keys {portable_to_raw[collision_key]!r} and "
                f"{key!r} collide case-insensitively at one payload path."
            )
        raw_keys.add(key)
        portable_to_raw[collision_key] = key
        portable_keys.append(portable_key)
    return portable_keys


def _require_string_identifiers(
    values: Sequence[object],
    *,
    label: str,
) -> list[str]:
    identifiers: list[str] = []
    for index, value in enumerate(values):
        if not isinstance(value, str):
            raise TypeError(
                f"{label} identifier at position {index} must be a string, "
                f"got {type(value).__name__}."
            )
        if not value:
            raise ValueError(f"{label} identifier at position {index} must be non-empty.")
        identifiers.append(value)
    return identifiers


def _require_nonempty_string(value: object, *, label: str) -> str:
    """Require explicit text without deriving or coercing an identity."""
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a non-empty string.")
    if not value:
        raise ValueError(f"{label} must be a non-empty string.")
    return value


def _require_native_boolean(value: object, *, label: str) -> bool:
    """Require one native boolean without truth-value coercion."""
    if type(value) is not bool:
        raise TypeError(f"{label} must be exactly True or False.")
    return value


def _require_positive_native_integer(value: object, *, label: str) -> int:
    """Require one positive native integer without numeric coercion."""
    if type(value) is not int:
        raise TypeError(f"{label} must be a native integer.")
    if value < 1:
        raise ValueError(f"{label} must be a positive integer.")
    return value


def _require_optional_native_string(
    value: object,
    *,
    label: str,
    allow_empty: bool,
) -> str | None:
    """Validate optional identity text without omission or string coercion."""
    if value is None:
        return None
    if type(value) is not str:
        raise TypeError(f"{label} must be None or a native string.")
    if not allow_empty and (not value or not value.strip()):
        raise ValueError(f"{label} must be None or a non-empty string.")
    return value


def _require_finite_float32_array(
    values: object,
    *,
    label: str,
) -> np.ndarray:
    """Validate real finite values that remain finite after float32 conversion."""
    dense = cast(sparse.spmatrix, values).toarray() if sparse.issparse(values) else values
    array = np.asarray(dense)
    if array.dtype.kind not in {"i", "u", "f"}:
        raise TypeError(f"{label} must contain real numeric values.")
    if not np.isfinite(array).all():
        raise ValueError(f"{label} must contain only finite values.")
    with np.errstate(over="ignore", invalid="ignore"):
        float32_array = array.astype(np.float32, copy=False)
    if not np.isfinite(float32_array).all():
        raise ValueError(f"{label} contains values outside the finite float32 range.")
    return float32_array


def _require_continuous_obs_values(
    values: object,
    *,
    key: str,
    n_cells: int,
) -> np.ndarray:
    """Require one exact finite float32 observation vector."""
    array = _require_finite_float32_array(
        values,
        label=f"Continuous obs field {key!r}",
    )
    if array.ndim != 1 or array.shape[0] != n_cells:
        raise ValueError(
            f"Continuous obs field {key!r} must have shape ({n_cells},), got {array.shape}."
        )
    return array


def _normalize_finite_float32_embedding(
    embedding: np.ndarray,
    *,
    label: str,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Normalize one validated embedding with a required nonzero range."""
    if embedding.ndim != 2 or embedding.shape[0] == 0:
        raise ValueError(f"{label} must be a non-empty 2D array.")
    working = embedding.astype(np.float64)
    axis_mins = working.min(axis=0)
    axis_maxs = working.max(axis=0)
    max_range = float((axis_maxs - axis_mins).max())
    if max_range <= 0:
        raise ValueError(f"{label} has no coordinate variation and cannot be normalized.")
    center = (axis_mins + axis_maxs) / 2
    scale_factor = 2.0 / max_range
    normalized = ((working - center) * scale_factor).astype(np.float32)
    if not np.isfinite(normalized).all():
        raise ValueError(f"{label} normalization produced non-finite coordinates.")
    return normalized, center.astype(np.float32), scale_factor, max_range


def _json_category_values(
    values: Iterable[object],
    *,
    field_key: str,
) -> list[JsonScalar]:
    """Preserve one exact JSON-scalar identity for each category label."""
    categories: list[str | bool | int | float] = []
    seen: set[tuple[str, object]] = set()
    for raw_value in values:
        value = raw_value.item() if isinstance(raw_value, np.generic) else raw_value
        token: tuple[str, object]
        if isinstance(value, bool):
            token = ("boolean", value)
        elif isinstance(value, str):
            token = ("string", value)
        elif isinstance(value, Integral):
            integer_value = int(value)
            if abs(integer_value) > 9_007_199_254_740_991:
                raise ValueError(
                    f"Categorical field {field_key!r} contains integer label "
                    f"{integer_value!r} outside JavaScript's exact integer range."
                )
            value = integer_value
            token = ("number", value)
        elif isinstance(value, Real) and math.isfinite(value):
            value = float(value)
            token = ("number", value)
        else:
            raise ValueError(f"Categorical field {field_key!r} labels must be finite JSON scalars.")
        if token in seen:
            raise ValueError(
                f"Categorical field {field_key!r} labels must be unique "
                "after exact JSON representation."
            )
        seen.add(token)
        categories.append(value)
    return categories


def _identity_obs_fields_from_compact_manifest(
    manifest: dict,
) -> list[dict]:
    expected_keys = {
        "_format",
        "n_points",
        "centroid_outlier_quantile",
        "latent_key",
        "compression",
        "_obsSchemas",
        "_continuousFields",
        "_categoricalFields",
    }
    if not isinstance(manifest, dict) or set(manifest) != expected_keys:
        raise ValueError("obs manifest must contain exactly the current compact_v1 fields.")
    if manifest["_format"] != MANIFEST_FORMAT_VERSION:
        raise ValueError("obs manifest must use the current compact_v1 format.")

    continuous_fields = manifest["_continuousFields"]
    categorical_fields = manifest["_categoricalFields"]
    if not isinstance(continuous_fields, list):
        raise ValueError("compact_v1 obs manifest _continuousFields must be a list.")
    if not isinstance(categorical_fields, list):
        raise ValueError("compact_v1 obs manifest _categoricalFields must be a list.")

    identity_fields: list[dict] = []
    manifest_keys: list[str] = []
    for field in continuous_fields:
        if (
            not isinstance(field, list)
            or len(field) not in (1, 3)
            or not isinstance(field[0], str)
            or not field[0]
        ):
            raise ValueError(
                "compact_v1 continuous observation fields must be exact "
                "[key] or [key, minValue, maxValue] tuples."
            )
        key = field[0]
        manifest_keys.append(key)
        identity_fields.append({"key": key, "kind": "continuous"})

    for field in categorical_fields:
        if (
            not isinstance(field, list)
            or len(field) not in (5, 7)
            or not isinstance(field[0], str)
            or not field[0]
            or not isinstance(field[1], list)
        ):
            raise ValueError(
                "compact_v1 categorical observation fields must be exact "
                "five- or seven-member tuples with a category array."
            )
        key = field[0]
        manifest_keys.append(key)
        identity_fields.append(
            {
                "key": key,
                "kind": "category",
                "n_categories": len(field[1]),
            }
        )

    _assert_unique_filename_components(
        manifest_keys,
        label="Observation field",
    )
    return identity_fields


def _to_dense(arr: np.ndarray | sparse.spmatrix) -> np.ndarray:
    """Convert sparse matrix to dense numpy array if necessary."""
    if sparse.issparse(arr):
        return np.asarray(cast(sparse.spmatrix, arr).toarray())
    return np.asarray(arr)


def _output_path_is_writable(
    path: Path,
    description: str,
) -> bool:
    """Admit one new staged artifact path."""
    if path.exists():
        raise FileExistsError(f"Cannot write {description}: staged path already exists: {path}")
    return True


def _publish_export_generation(staging_dir: Path, target_dir: Path) -> None:
    """Publish one complete staged directory, restoring the prior one on failure."""
    if not staging_dir.is_dir():
        raise RuntimeError(f"Staged export directory is missing: {staging_dir}")

    prior_dir: Path | None = None
    if target_dir.exists():
        backup_name = tempfile.mkdtemp(
            prefix=f".{target_dir.name}.cellucid-backup-",
            dir=target_dir.parent,
        )
        prior_dir = Path(backup_name)
        prior_dir.rmdir()
        target_dir.rename(prior_dir)

    try:
        staging_dir.rename(target_dir)
    except BaseException:
        if prior_dir is not None:
            try:
                prior_dir.rename(target_dir)
            except BaseException as restore_error:
                raise RuntimeError(
                    f"Failed to publish {target_dir} and failed to restore its prior "
                    f"generation from {prior_dir}."
                ) from restore_error
        raise

    if prior_dir is not None:
        shutil.rmtree(prior_dir)


@contextmanager
def _exclusive_export_generation(target_dir: Path) -> Iterator[None]:
    """Hold one non-waiting cross-process writer lock for an export target."""
    lock_path = target_dir.parent / f".{target_dir.name}.cellucid.lock"
    try:
        descriptor = os.open(
            lock_path,
            os.O_CREAT | os.O_EXCL | os.O_WRONLY,
            0o600,
        )
    except FileExistsError as error:
        raise RuntimeError(f"An export generation is already active for {target_dir}.") from error

    try:
        with os.fdopen(descriptor, "wb") as lock_file:
            lock_file.write(f"{os.getpid()}\n".encode("ascii"))
            lock_file.flush()
            os.fsync(lock_file.fileno())
        yield
    finally:
        lock_path.unlink()


def _write_binary(
    path: Path,
    data: np.ndarray,
    compression: int | None = None,
) -> Path:
    """
    Write binary data, optionally with gzip compression.

    Parameters
    ----------
    path : Path
        Output path. If compression is enabled, '.gz' will be appended.
    data : np.ndarray
        Data to write.
    compression : int or None
        Gzip compression level (1-9). None means no compression.

    Returns
    -------
    Path
        Actual path written (may have .gz suffix).
    """
    if compression is not None and (type(compression) is not int or not 1 <= compression <= 9):
        raise ValueError("compression must be None or one integer from 1 to 9.")

    if compression is not None:
        gz_path = Path(str(path) + ".gz")
        with gzip.open(gz_path, "wb", compresslevel=compression) as f:
            f.write(data.tobytes())
        return gz_path
    else:
        data.tofile(path)
        return path


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
        Scale factor (for dequantization).
    """
    if type(bits) is not int or bits not in (8, 16):
        raise ValueError("Quantization bits must be exactly 8 or 16.")

    if values.size == 0:
        raise ValueError(
            f"Continuous field {field_name!r} cannot be quantized because it has no finite values."
        )
    if not np.isfinite(values).all():
        raise ValueError(f"Continuous field {field_name!r} must contain only finite values.")

    min_val = float(np.min(values))
    max_val = float(np.max(values))

    if max_val == min_val:
        raise ValueError(
            f"Continuous field {field_name!r} cannot be quantized because "
            f"all values are the constant {min_val!r}."
        )

    if bits == 8:
        max_quant = 254  # 255 is the categorical-outlier missing marker.
        dtype: type[np.uint8] | type[np.uint16] = np.uint8
    else:  # 16 bits
        max_quant = 65534  # 65535 is the categorical-outlier missing marker.
        dtype = np.uint16

    scale = max_quant / (max_val - min_val)
    normalized = (values - min_val) * scale
    quantized = np.clip(normalized, 0, max_quant).astype(dtype)

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


def _validate_prepare_codec_options(
    *,
    compression: int | None,
    var_quantization: int | None,
    obs_continuous_quantization: int | None,
    obs_categorical_dtype: object,
    centroid_outlier_quantile: float | None,
) -> None:
    """Validate the sole current compact_v1 codec configuration."""
    if compression is not None and (type(compression) is not int or not 1 <= compression <= 9):
        raise ValueError("compression must be None or one integer from 1 to 9.")

    for name, value in (
        ("var_quantization", var_quantization),
        ("obs_continuous_quantization", obs_continuous_quantization),
    ):
        if value is not None and (type(value) is not int or value not in (8, 16)):
            raise ValueError(f"{name} must be None, 8, or 16.")

    if type(obs_categorical_dtype) is not str or obs_categorical_dtype not in (
        "uint8",
        "uint16",
    ):
        raise ValueError("obs_categorical_dtype must be exactly 'uint8' or 'uint16'.")

    if centroid_outlier_quantile is not None and (
        isinstance(centroid_outlier_quantile, bool)
        or not isinstance(centroid_outlier_quantile, Real)
        or not math.isfinite(centroid_outlier_quantile)
        or centroid_outlier_quantile <= 0.5
        or centroid_outlier_quantile >= 1
    ):
        raise ValueError(
            "centroid_outlier_quantile must be None or one finite number "
            "strictly between 0.5 and 1."
        )


def _select_category_dtype(
    n_categories: int,
) -> tuple[type[np.uint8] | type[np.uint16], int]:
    """
    Select optimal dtype for category codes based on number of categories.

    Parameters
    ----------
    n_categories : int
        Number of unique categories (not counting missing).

    Returns
    -------
    dtype : np.dtype
        Optimal dtype (uint8 or uint16).
    missing_value : int
        Value to use for missing/NaN codes.
    """
    if n_categories > 65_535:
        raise ValueError(
            f"Categorical field has {n_categories:,} categories; "
            "the current contract supports at most 65,535."
        )
    if n_categories <= 255:
        # uint8 has 255 category codes (0-254), with 255 reserved for missing.
        return np.uint8, 255
    # uint16 has 65,535 category codes (0-65534), with 65535 for missing.
    return np.uint16, 65535


def _compute_centroids_for_field(
    coords: np.ndarray,
    codes: np.ndarray,
    categories: list[JsonScalar],
    outlier_quantile: float = 0.95,
    min_points: int = 10,
) -> list[dict]:
    """
    Compute centroids per category with outlier removal (based on embedding coords for display).

    Works with any dimensionality (1D, 2D, 3D, etc.) - the position will have the same
    number of dimensions as the input coords.
    """
    if coords.shape[0] != codes.shape[0]:
        raise ValueError("coords and codes must have the same length.")

    centroids: list[dict[str, object]] = []

    if not (0.5 < outlier_quantile < 1.0):
        raise ValueError("outlier_quantile must be strictly between 0.5 and 1.0")

    for code, label in enumerate(categories):
        mask = codes == code
        idx = np.nonzero(mask)[0]
        n = idx.size
        if n < min_points:
            continue

        pts = coords[idx, :]  # (n, ndim)
        center = pts.mean(axis=0)

        if n > min_points:
            dists = np.linalg.norm(pts - center, axis=1)
            thr = float(np.quantile(dists, outlier_quantile))
            inlier_mask = dists <= thr
            n_in = int(inlier_mask.sum())
            if n_in >= min_points:
                pts_in = pts[inlier_mask, :]
                center = pts_in.mean(axis=0)
                used_count = n_in
            else:
                used_count = n
        else:
            used_count = n

        centroids.append(
            {
                "category": label,
                "position": center.astype(float).tolist(),
                "n_points": int(used_count),
            }
        )

    return centroids


def _compute_centroids_for_all_dimensions(
    embeddings: dict[int, np.ndarray],
    codes: np.ndarray,
    categories: list[JsonScalar],
    outlier_quantile: float = 0.95,
    min_points: int = 10,
) -> dict[int, list[dict]]:
    """
    Compute centroids for each available dimension.

    Returns a dictionary keyed by dimension (1, 2, 3) with centroid lists for each.
    """
    centroids_by_dim = {}
    for dim, coords in embeddings.items():
        centroids_by_dim[dim] = _compute_centroids_for_field(
            coords, codes, categories, outlier_quantile, min_points
        )
    return centroids_by_dim


def _compute_latent_space_quantiles(
    latent: np.ndarray,
    codes: np.ndarray,
    categories: list[JsonScalar],
    min_points: int = 10,
) -> np.ndarray:
    """
    Compute per-cell outlier quantiles based on latent space distances to category centroids.
    """
    n_cells = latent.shape[0]
    quantiles = np.full(n_cells, np.nan, dtype=np.float32)

    for code, _label in enumerate(categories):
        mask = codes == code
        idx = np.nonzero(mask)[0]
        n = idx.size

        if n < min_points:
            continue

        pts = latent[idx, :]
        centroid = pts.mean(axis=0)
        dists = np.linalg.norm(pts - centroid, axis=1)
        sorted_dists = np.sort(dists)
        ranks = np.searchsorted(sorted_dists, dists, side="right")
        cell_quantiles = ranks.astype(np.float32) / n
        quantiles[idx] = cell_quantiles

    return quantiles


def _prepare_generation(
    latent_space: np.ndarray | sparse.spmatrix | None = None,
    obs: pd.DataFrame | None = None,
    var: pd.DataFrame | None = None,
    gene_expression: np.ndarray | sparse.spmatrix | None = None,
    var_gene_id_column: str | None = None,
    gene_identifiers: Sequence[str] | None = None,
    connectivities: np.ndarray | sparse.spmatrix | None = None,
    out_dir: Path | str = DEFAULT_EXPORT_DIR,
    obs_keys: Sequence[str] | None = None,
    centroid_outlier_quantile: float = 0.95,
    centroid_min_points: int = 10,
    # Optimization parameters (all disabled by default)
    var_quantization: int | None = None,
    obs_continuous_quantization: int | None = None,
    compression: int | None = None,
    # Dataset metadata parameters (for dataset_identity.json)
    *,
    _published_out_dir: Path,
    obs_categorical_dtype: Literal["uint8", "uint16"],
    dataset_name: str,
    dataset_id: str,
    dataset_description: str | None = None,
    source_name: str | None = None,
    source_url: str | None = None,
    source_citation: str | None = None,
    # Multi-dimensional embedding parameters (at least one required)
    X_umap_1d: np.ndarray | None = None,
    X_umap_2d: np.ndarray | None = None,
    X_umap_3d: np.ndarray | None = None,
    # Optional per-cell vector fields aligned to the embedding(s)
    # (e.g. scVelo velocity, CellRank drift vectors)
    vector_fields: dict[str, np.ndarray | sparse.spmatrix] | None = None,
    vector_field_default: str | None = None,
) -> None:
    """
    Export raw data arrays to files used by the WebGL viewer.

    Memory/Disk Optimization Options
    --------------------------------
    var_quantization : int or None
        Bits for gene expression quantization (8, 16, or None for full float32).
        8-bit reduces file size by 4x with minimal visual impact for colormapping.
    obs_continuous_quantization : int or None
        Bits for continuous obs field quantization (8, 16, or None for full float32).
    obs_categorical_dtype : 'uint8' or 'uint16'
        - 'uint8': Store up to 255 categories
        - 'uint16': Store up to 65,535 categories
    compression : int or None
        Gzip compression level (1-9). None disables compression.
        Level 6 is a good balance of speed and size. Files get .gz extension.

    Multi-Dimensional Embeddings
    ----------------------------
    At least one dimensional embedding must be provided. The viewer supports
    switching between different dimensionalities of the same data at runtime.
    All embeddings must have the same number of cells (rows) but different
    column counts matching their dimensionality.

    IMPORTANT: Each embedding is normalized independently to fit within the
    [-1, 1] coordinate range. Within each dimension, the same scale factor is
    used for all axes to preserve aspect ratios. This ensures each dimension
    fills the viewing area optimally without requiring manual zoom adjustment.

    X_umap_1d : np.ndarray, optional
        1D embedding coordinates, shape (n_cells, 1). Stored as points_1d.bin.
    X_umap_2d : np.ndarray, optional
        2D embedding coordinates, shape (n_cells, 2). Stored as points_2d.bin.
    X_umap_3d : np.ndarray, optional
        3D embedding coordinates, shape (n_cells, 3). Stored as points_3d.bin.
        This is the primary visualization and is used for centroid computation.

    vector_fields : dict[str, np.ndarray] or None
        Optional per-cell displacement vectors aligned to the embedding space.
        Every key must use the exact dimension-suffixed AnnData ``obsm`` form
        ``<field>_umap_<dim>d`` (for example, ``velocity_umap_2d`` or
        ``T_fwd_umap_3d``).

        Each value must be shaped exactly ``(n_cells, dim)`` and contain finite
        real values representable as float32.
        Vectors are scaled by the same per-dimension normalization scale as points.
        A field's metadata ``default_dimension`` is exactly its highest available
        vector dimension.
    vector_field_default : str or None
        Exact field id to select initially. Required when more than one vector
        field is declared; omitted for a single unambiguous field.

    Standard Parameters
    -------------------
    latent_space : np.ndarray or sparse matrix
        Latent space for outlier quantile calculation, shape (n_cells, n_dims).
    obs : pd.DataFrame
        Cell metadata, shape (n_cells, n_obs_columns).
    var : pd.DataFrame, optional
        Gene/feature metadata. Required if gene_expression is provided.
    gene_expression : np.ndarray or sparse matrix, optional
        Gene expression matrix, shape (n_cells, n_genes).
    var_gene_id_column : str or None
        Exact non-empty column name in var containing string gene identifiers.
        None uses string identifiers from var.index.
    gene_identifiers : sequence of str, optional
        Unique gene identifiers to export. Every identifier must be a non-empty
        string present in var. If None, all genes are exported.
    connectivities : array or sparse matrix, optional
        Exact weighted undirected graph, shape ``(n_cells, n_cells)``. Values
        must be finite and non-negative, the topology and weights must be
        exactly symmetric, and the diagonal must be zero. Sparse inputs must
        not contain duplicate coordinates or stored zero entries.
    out_dir : Path or str
        Output directory (default: exports/ under the current working directory).
    obs_keys : sequence of str or None
        Which obs columns to export. If None, all columns are exported.
    centroid_outlier_quantile : float
        Quantile of distances to keep as inliers when computing centroids.
    centroid_min_points : int
        Minimum number of points in a category to compute a centroid.
    force : bool
        If True, replace files in the export directory. If False, the export
        directory must be absent or empty.

    Dataset Metadata Parameters
    ---------------------------
    dataset_name : str
        Explicit human-readable name for the dataset.
    dataset_description : str, optional
        Description of the dataset.
    dataset_id : str
        Explicit unique identifier for the dataset.
    source_name : str, optional
        Name of the data source (e.g., "HLCA Consortium").
    source_url : str, optional
        URL to the data source.
    source_citation : str, optional
        Citation text for the data source.
    """
    _validate_prepare_codec_options(
        compression=compression,
        var_quantization=var_quantization,
        obs_continuous_quantization=obs_continuous_quantization,
        obs_categorical_dtype=obs_categorical_dtype,
        centroid_outlier_quantile=centroid_outlier_quantile,
    )
    dataset_name = _require_nonempty_string(
        dataset_name,
        label="dataset_name",
    )
    dataset_id = _require_nonempty_string(
        dataset_id,
        label="dataset_id",
    )

    out_dir = Path(out_dir)
    if out_dir.exists():
        if not out_dir.is_dir():
            raise NotADirectoryError(f"Export path exists but is not a directory: {out_dir}")
        if any(out_dir.iterdir()):
            raise FileExistsError(f"Staged export directory must be empty: {out_dir}")
    obs_manifest_filename = "obs_manifest.json"
    obs_binary_dirname = DEFAULT_OBS_DIRNAME
    var_manifest_filename = "var_manifest.json"
    var_binary_dirname = DEFAULT_VAR_DIRNAME
    obs_binary_dir = out_dir / obs_binary_dirname

    def published_path(path: Path) -> Path:
        return _published_out_dir / path.relative_to(out_dir)

    # =========================================================================
    # MULTI-DIMENSIONAL EMBEDDING VALIDATION & PROCESSING
    # =========================================================================
    # Collect all provided embeddings
    embeddings: dict[int, np.ndarray] = {}
    if X_umap_1d is not None:
        embeddings[1] = _require_finite_float32_array(
            X_umap_1d,
            label="X_umap_1d",
        )
    if X_umap_2d is not None:
        embeddings[2] = _require_finite_float32_array(
            X_umap_2d,
            label="X_umap_2d",
        )
    if X_umap_3d is not None:
        embeddings[3] = _require_finite_float32_array(
            X_umap_3d,
            label="X_umap_3d",
        )

    if not embeddings:
        raise ValueError(
            "At least one dimensional embedding must be provided. "
            "Use X_umap_1d, X_umap_2d, or X_umap_3d."
        )

    # Validate each embedding has correct dimensions
    n_cells = None
    for dim, arr in embeddings.items():
        if arr.ndim != 2:
            raise ValueError(f"X_umap_{dim}d must be a 2D array, got shape {arr.shape}.")
        if arr.shape[1] != dim:
            raise ValueError(
                f"X_umap_{dim}d must have exactly {dim} columns, got {arr.shape[1]}. "
                f"Shape is {arr.shape}."
            )
        if n_cells is None:
            n_cells = arr.shape[0]
        elif arr.shape[0] != n_cells:
            raise ValueError(
                f"All embeddings must have the same number of cells. "
                f"First embedding has {n_cells} cells, but X_umap_{dim}d has {arr.shape[0]} cells."
            )
    if n_cells is None or n_cells <= 0:
        raise ValueError("Embeddings must contain at least one cell.")

    connectivity_edges: ConnectivityEdgePairs | None = None
    if connectivities is not None:
        connectivity_edges = validate_connectivity_edges(
            connectivities,
            n_cells=n_cells,
        )

    # =========================================================================
    # NORMALIZE EACH EMBEDDING INDEPENDENTLY TO FIT WITHIN [-1, 1] RANGE
    # =========================================================================
    # Each dimensional embedding (1D, 2D, 3D) is normalized independently so that
    # it fills the viewing area optimally. Within each dimension, we use the same
    # scale factor for all axes to preserve aspect ratios.
    #
    # This ensures that switching between dimensions doesn't require manual zoom
    # adjustments - each dimension will fill the view appropriately.

    normalization_info = {}
    for dim, arr in embeddings.items():
        normalized, center, scale_factor, max_range = _normalize_finite_float32_embedding(
            arr,
            label=f"X_umap_{dim}d",
        )
        embeddings[dim] = normalized

        # Store info for logging
        normalization_info[dim] = {
            "original_range": max_range,
            "center": center.tolist(),
            "scale_factor": scale_factor,
        }

    # Track available dimensions for metadata
    available_dimensions = sorted(embeddings.keys())

    # The highest declared dimension is the one exact initial dimension.
    default_dimension = max(available_dimensions)

    # Print export settings summary
    print("=" * 60)
    print("Export Settings:")
    print(f"  Output directory: {_published_out_dir}")
    print(f"  Compression: {'gzip level ' + str(compression) if compression else 'disabled'}")
    print(
        f"  Var (gene) quantization: {str(var_quantization) + '-bit' if var_quantization else 'disabled (float32)'}"
    )
    print(
        f"  Obs continuous quantization: {str(obs_continuous_quantization) + '-bit' if obs_continuous_quantization else 'disabled (float32)'}"
    )
    print(f"  Obs categorical dtype: {obs_categorical_dtype}")
    print(f"  Available dimensions: {available_dimensions}")
    print(f"  Default dimension: {default_dimension}D")
    print("  Coordinate normalization (per-dimension, aspect-ratio preserved):")
    for dim in sorted(normalization_info.keys()):
        info = normalization_info[dim]
        print(f"    {dim}D: range {info['original_range']:.2f} → [-1, 1]")
    print("=" * 60)

    # Validate and convert latent space
    if latent_space is None:
        raise ValueError("latent_space is required for outlier quantile calculation.")
    latent = _require_finite_float32_array(
        latent_space,
        label="latent_space",
    )
    if latent.ndim != 2 or latent.shape[0] != n_cells:
        raise ValueError(
            f"latent_space must have shape ({n_cells}, n_dimensions), got {latent.shape}."
        )

    # Validate obs
    if obs is None:
        raise ValueError("obs DataFrame is required.")
    if len(obs) != n_cells:
        raise ValueError(f"obs has {len(obs)} rows, but embeddings have {n_cells} cells.")

    # Resolve and validate every observation payload identity before any output
    # path is created.
    if obs_keys is None:
        obs_keys = list(obs.columns)
    else:
        if isinstance(obs_keys, str | bytes):
            raise TypeError("obs_keys must be a sequence of strings, not a string scalar.")
        obs_keys = list(obs_keys)
        missing = [key for key in obs_keys if key not in obs.columns]
        if missing:
            raise KeyError(
                f"obs_keys contain columns not in obs: {missing}. "
                f"Available columns: {list(obs.columns)}"
            )
    safe_obs_keys = _assert_unique_filename_components(
        obs_keys,
        label="Observation field",
    )

    validated_continuous_obs: dict[str, np.ndarray] = {}
    obs_field_summaries: list[dict[str, object]] = []
    for key in obs_keys:
        series = obs[key]
        if isinstance(series.dtype, pd.CategoricalDtype) or pd.api.types.is_bool_dtype(series):
            kind = "category"
        elif pd.api.types.is_numeric_dtype(series):
            kind = "continuous"
        elif pd.api.types.is_string_dtype(series) or pd.api.types.is_object_dtype(series):
            kind = "category"
        else:
            raise TypeError(f"obs field {key!r} has unsupported dtype {series.dtype!r}.")

        if kind == "continuous":
            validated_continuous_obs[key] = _require_continuous_obs_values(
                series.to_numpy(),
                key=key,
                n_cells=n_cells,
            )
            if obs_continuous_quantization is None:
                dtype_str = "float32"
            elif obs_continuous_quantization == 8:
                dtype_str = "uint8"
            else:
                dtype_str = "uint16"
            obs_field_summaries.append(
                {
                    "key": key,
                    "kind": "continuous",
                    "quantized": obs_continuous_quantization is not None,
                    "quantization_bits": int(obs_continuous_quantization)
                    if obs_continuous_quantization is not None
                    else None,
                    "dtype": dtype_str,
                }
            )
            continue

        categorical = series.astype("category")
        categories = _json_category_values(
            categorical.cat.categories,
            field_key=key,
        )
        n_categories = len(categories)
        if obs_categorical_dtype == "uint8":
            if n_categories > 255:
                raise ValueError(
                    f"Field {key!r} has {n_categories} categories, but uint8 supports at most 255."
                )
            dtype_str = "uint8"
        else:
            if n_categories > 65_535:
                raise ValueError(
                    f"Field {key!r} has {n_categories} categories, "
                    "but uint16 supports at most 65,535."
                )
            dtype_str = "uint16"
        obs_field_summaries.append(
            {
                "key": key,
                "kind": "category",
                "category_count": n_categories,
                "codes_dtype": dtype_str,
                "outlier_quantized": obs_continuous_quantization is not None,
                "outlier_quantization_bits": int(obs_continuous_quantization)
                if obs_continuous_quantization is not None
                else None,
            }
        )

    # Resolve every gene payload identity before any earlier dataset artifact can
    # be written.
    genes_to_export: list[str] = []
    gene_id_to_idx: dict[str, int] = {}
    safe_gene_id_by_id: dict[str, str] = {}
    gene_expr_is_sparse = False
    gene_expression_for_export: np.ndarray | sparse.csc_matrix | None = None
    if gene_expression is not None:
        if var is None:
            raise ValueError("var DataFrame must be provided when gene_expression is given.")

        if sparse.issparse(gene_expression):
            gene_expr_is_sparse = True
            gene_expression_for_export = cast(
                sparse.spmatrix,
                gene_expression,
            ).tocsc()
        else:
            gene_expression_for_export = np.asarray(gene_expression)
        if gene_expression_for_export.ndim != 2:
            raise ValueError(
                "gene_expression must be a 2D array or sparse matrix; "
                f"got shape {gene_expression_for_export.shape}."
            )

        n_expr_cells, n_genes = gene_expression_for_export.shape
        if n_expr_cells != n_cells:
            raise ValueError(
                f"gene_expression has {n_expr_cells} cells, but embeddings have {n_cells} cells."
            )
        if len(var) != n_genes:
            raise ValueError(f"var has {len(var)} rows, but gene_expression has {n_genes} genes.")

        if var_gene_id_column is None:
            all_gene_ids = _require_string_identifiers(
                var.index.tolist(),
                label="var index",
            )
        else:
            if type(var_gene_id_column) is not str:
                raise TypeError(
                    "var_gene_id_column must be a native non-empty string or None; "
                    f"got {type(var_gene_id_column).__name__}."
                )
            if not var_gene_id_column:
                raise ValueError("var_gene_id_column must be a native non-empty string or None.")
            if var_gene_id_column not in var.columns:
                raise KeyError(
                    f"var_gene_id_column '{var_gene_id_column}' not found in var. "
                    f"Available columns: {list(var.columns)}"
                )
            all_gene_ids = _require_string_identifiers(
                var[var_gene_id_column].tolist(),
                label=f"var column {var_gene_id_column!r}",
            )
        safe_all_gene_ids = _assert_unique_filename_components(
            all_gene_ids,
            label="Gene",
        )
        gene_id_to_idx = {gene_id: index for index, gene_id in enumerate(all_gene_ids)}
        safe_gene_id_by_id = dict(zip(all_gene_ids, safe_all_gene_ids, strict=True))

        if gene_identifiers is None:
            genes_to_export = all_gene_ids
        else:
            if isinstance(gene_identifiers, str | bytes):
                raise TypeError(
                    "gene_identifiers must be a sequence of strings, not a string scalar."
                )
            genes_to_export = _require_string_identifiers(
                gene_identifiers,
                label="Requested gene",
            )
            _assert_unique_filename_components(
                genes_to_export,
                label="Requested gene",
            )
            missing_genes = [
                gene_id for gene_id in genes_to_export if gene_id not in gene_id_to_idx
            ]
            if missing_genes:
                raise KeyError(
                    f"gene_identifiers contain identifiers not present in var: {missing_genes}."
                )

    validated_vector_fields = validate_vector_fields(
        vector_fields,
        n_cells=n_cells,
        available_dimensions=available_dimensions,
        vector_field_default=vector_field_default,
    )
    scaled_vector_fields: dict[str, dict[int, np.ndarray]] = {}
    vector_fields_identity: dict[str, object] | None = None
    if validated_vector_fields.fields:
        fields_metadata: dict[str, dict[str, object]] = {}
        gzip_suffix = ".gz" if compression else ""
        _assert_unique_filename_components(
            list(validated_vector_fields.fields),
            label="Vector field",
        )
        for field_id, vectors_by_dimension in validated_vector_fields.fields.items():
            dimensions = sorted(vectors_by_dimension)
            scaled_vector_fields[field_id] = {}
            for dimension, vectors in vectors_by_dimension.items():
                scale_factor = normalization_info[dimension]["scale_factor"]
                scaled_vector_fields[field_id][dimension] = scale_vector_field(
                    vectors,
                    scale_factor=scale_factor,
                    label=f"Vector field {field_id!r} {dimension}D",
                )

            fields_metadata[field_id] = {
                "label": field_id,
                "basis": "umap",
                "available_dimensions": dimensions,
                "default_dimension": max(dimensions),
                "files": {
                    f"{dimension}d": (f"vectors/{field_id}_{dimension}d.bin{gzip_suffix}")
                    for dimension in dimensions
                },
            }

        vector_fields_identity = {
            "default_field": validated_vector_fields.default_field,
            "fields": fields_metadata,
        }

    out_dir.mkdir(parents=True, exist_ok=True)
    obs_binary_dir.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # SAVE DIMENSIONAL EMBEDDING FILES
    # =========================================================================
    for dim, arr in embeddings.items():
        dim_filename = f"points_{dim}d.bin"
        dim_path = out_dir / dim_filename
        check_path = Path(str(dim_path) + ".gz") if compression and compression > 0 else dim_path
        if _output_path_is_writable(check_path, check_path.name):
            actual_path = _write_binary(dim_path, arr, compression)
            suffix = " (gzip)" if compression else ""
            print(
                f"✓ Wrote {dim}D positions ({arr.shape[0]:,} cells × {dim} dims) "
                f"to {published_path(actual_path)}{suffix}"
            )

    # =========================================================================
    # SAVE VECTOR FIELDS (OPTIONAL)
    # =========================================================================
    if scaled_vector_fields:
        vectors_dir = out_dir / "vectors"
        vectors_dir.mkdir(parents=True, exist_ok=True)
        for field_id, vectors_by_dimension in scaled_vector_fields.items():
            for dimension, vectors in vectors_by_dimension.items():
                filename = f"{field_id}_{dimension}d.bin"
                path = vectors_dir / filename
                check_path = Path(str(path) + ".gz") if compression and compression > 0 else path
                if _output_path_is_writable(check_path, check_path.name):
                    actual_path = _write_binary(path, vectors, compression)
                    suffix = " (gzip)" if compression else ""
                    print(
                        f"✓ Wrote vector field '{field_id}' {dimension}D "
                        f"({vectors.shape[0]:,} cells × {dimension} comps) "
                        f"to {published_path(actual_path)}{suffix}"
                    )

    # Check if obs manifest already exists
    obs_manifest_path = out_dir / obs_manifest_filename
    if _output_path_is_writable(obs_manifest_path, "obs manifest"):
        # Compact format: separate lists for continuous and categorical fields
        obs_continuous_fields: list = []
        obs_categorical_fields: list = []
        # Track dtype info for schema (will use first encountered)
        continuous_dtype_info: dict = {}
        categorical_dtype_info: dict = {}

        for key, safe_key in zip(obs_keys, safe_obs_keys, strict=True):
            s = obs[key]

            # Decide kind: continuous vs categorical
            if isinstance(s.dtype, pd.CategoricalDtype) or pd.api.types.is_bool_dtype(s):
                kind = "category"
            elif pd.api.types.is_numeric_dtype(s):
                kind = "continuous"
            else:
                kind = "category"

            if kind == "continuous":
                values = validated_continuous_obs[key]

                # Apply quantization if requested
                if obs_continuous_quantization is not None:
                    quantized, min_val, max_val, scale = _quantize_continuous(
                        values, bits=obs_continuous_quantization, field_name=key
                    )

                    if obs_continuous_quantization == 8:
                        dtype_str = "uint8"
                        ext = "u8"
                    else:
                        dtype_str = "uint16"
                        ext = "u16"

                    value_fname = f"{safe_key}.values.{ext}"
                    value_path = obs_binary_dir / value_fname
                    actual_path = _write_binary(value_path, quantized, compression)

                    # Adjust path in manifest if compressed
                    manifest_path = f"{obs_binary_dirname}/{value_fname}"
                    if compression:
                        manifest_path += ".gz"

                    # Compact format: [key, minValue, maxValue]
                    obs_continuous_fields.append([key, min_val, max_val])
                    if not continuous_dtype_info:
                        continuous_dtype_info["ext"] = ext
                        continuous_dtype_info["dtype"] = dtype_str
                        continuous_dtype_info["quantized"] = True
                        continuous_dtype_info["quantizationBits"] = obs_continuous_quantization
                else:
                    # Full precision
                    value_fname = f"{safe_key}.values.f32"
                    value_path = obs_binary_dir / value_fname
                    actual_path = _write_binary(value_path, values, compression)

                    # Compact format: [key]
                    obs_continuous_fields.append([key])
                    if not continuous_dtype_info:
                        continuous_dtype_info["ext"] = "f32"
                        continuous_dtype_info["dtype"] = "float32"
                        continuous_dtype_info["quantized"] = False

            else:
                # Categorical
                cat = s.astype("category")
                categories = _json_category_values(
                    cat.cat.categories,
                    field_key=key,
                )
                codes = cat.cat.codes.to_numpy(dtype=np.int32)  # -1 for NaN

                if codes.shape[0] != n_cells:
                    raise ValueError(
                        f"Length mismatch for obs['{key}']: {codes.shape[0]} vs {n_cells}"
                    )

                n_categories = len(categories)

                dtype: type[np.uint8] | type[np.uint16]
                if obs_categorical_dtype == "uint8":
                    if n_categories > 255:
                        raise ValueError(
                            f"Field '{key}' has {n_categories} categories, "
                            "but uint8 supports at most 255."
                        )
                    dtype, missing_value = np.uint8, 255
                else:  # uint16
                    if n_categories > 65_535:
                        raise ValueError(
                            f"Field '{key}' has {n_categories} categories, "
                            "but uint16 supports at most 65,535."
                        )
                    dtype, missing_value = np.uint16, 65535

                codes_typed = np.full(n_cells, missing_value, dtype=dtype)
                valid_mask = codes >= 0
                codes_typed[valid_mask] = codes[valid_mask].astype(dtype)

                if dtype == np.uint8:
                    codes_fname = f"{safe_key}.codes.u8"
                    dtype_str = "uint8"
                else:
                    codes_fname = f"{safe_key}.codes.u16"
                    dtype_str = "uint16"

                codes_path = obs_binary_dir / codes_fname
                actual_path = _write_binary(codes_path, codes_typed, compression)

                manifest_codes_path = f"{obs_binary_dirname}/{codes_fname}"
                if compression:
                    manifest_codes_path += ".gz"

                # Compute centroids for all available dimensions
                if centroid_outlier_quantile is None:
                    centroids_by_dim = {dim: [] for dim in embeddings}
                else:
                    centroids_by_dim = _compute_centroids_for_all_dimensions(
                        embeddings,
                        codes,
                        categories,
                        outlier_quantile=centroid_outlier_quantile,
                        min_points=centroid_min_points,
                    )

                # Compute per-cell outlier quantiles based on latent space
                outlier_quantiles = _compute_latent_space_quantiles(
                    latent=latent,
                    codes=codes,
                    categories=categories,
                    min_points=centroid_min_points,
                )

                # Quantize outlier quantiles (they're always 0-1)
                if obs_continuous_quantization is not None:
                    oq_quantized, oq_min, oq_max, oq_scale = (
                        _quantize_nullable_outlier_quantiles(
                            outlier_quantiles,
                            bits=obs_continuous_quantization,
                            field_name=f"{key}_outliers",
                        )
                    )

                    if obs_continuous_quantization == 8:
                        oq_dtype_str = "uint8"
                        oq_ext = "u8"
                    else:
                        oq_dtype_str = "uint16"
                        oq_ext = "u16"

                    outlier_fname = f"{safe_key}.outliers.{oq_ext}"
                    outlier_path = obs_binary_dir / outlier_fname
                    _write_binary(outlier_path, oq_quantized, compression)

                    manifest_outlier_path = f"{obs_binary_dirname}/{outlier_fname}"
                    if compression:
                        manifest_outlier_path += ".gz"

                    # Compact format: [key, categories, codesDtype, codesMissingValue, centroidsByDim, outlierMinValue, outlierMaxValue]
                    # centroidsByDim is a dict keyed by dimension: {"1": [...], "2": [...], "3": [...]}
                    centroids_serializable = {
                        str(dim): cents for dim, cents in centroids_by_dim.items()
                    }
                    obs_categorical_fields.append(
                        [
                            key,
                            categories,
                            dtype_str,
                            int(missing_value),
                            centroids_serializable,
                            oq_min,
                            oq_max,
                        ]
                    )
                    if not categorical_dtype_info:
                        categorical_dtype_info["codesExt"] = "u8" if dtype == np.uint8 else "u16"
                        categorical_dtype_info["outlierExt"] = oq_ext
                        categorical_dtype_info["outlierDtype"] = oq_dtype_str
                        categorical_dtype_info["outlierQuantized"] = True
                else:
                    # Full precision outliers
                    outlier_fname = f"{safe_key}.outliers.f32"
                    outlier_path = obs_binary_dir / outlier_fname
                    _write_binary(outlier_path, outlier_quantiles.astype(np.float32), compression)

                    # Compact format: [key, categories, codesDtype, codesMissingValue, centroidsByDim]
                    # centroidsByDim is a dict keyed by dimension: {"1": [...], "2": [...], "3": [...]}
                    centroids_serializable = {
                        str(dim): cents for dim, cents in centroids_by_dim.items()
                    }
                    obs_categorical_fields.append(
                        [key, categories, dtype_str, int(missing_value), centroids_serializable]
                    )
                    if not categorical_dtype_info:
                        categorical_dtype_info["codesExt"] = "u8" if dtype == np.uint8 else "u16"
                        categorical_dtype_info["outlierExt"] = "f32"
                        categorical_dtype_info["outlierDtype"] = "float32"
                        categorical_dtype_info["outlierQuantized"] = False

        # Build compact manifest with schemas
        gz_suffix = ".gz" if compression else ""

        obs_schemas = {}
        if continuous_dtype_info:
            obs_schemas["continuous"] = {
                "pathPattern": f"{obs_binary_dirname}/{{key}}.values.{continuous_dtype_info['ext']}{gz_suffix}",
                "ext": continuous_dtype_info["ext"],
                "dtype": continuous_dtype_info["dtype"],
                "quantized": continuous_dtype_info.get("quantized", False),
            }
            if continuous_dtype_info.get("quantized"):
                obs_schemas["continuous"]["quantizationBits"] = continuous_dtype_info[
                    "quantizationBits"
                ]

        if categorical_dtype_info:
            obs_schemas["categorical"] = {
                "codesPathPattern": f"{obs_binary_dirname}/{{key}}.codes.{{ext}}{gz_suffix}",
                "outlierPathPattern": f"{obs_binary_dirname}/{{key}}.outliers.{categorical_dtype_info['outlierExt']}{gz_suffix}",
                "outlierExt": categorical_dtype_info["outlierExt"],
                "outlierDtype": categorical_dtype_info["outlierDtype"],
                "outlierQuantized": categorical_dtype_info.get("outlierQuantized", False),
            }

        obs_manifest_payload = {
            "_format": MANIFEST_FORMAT_VERSION,
            "n_points": int(n_cells),
            "centroid_outlier_quantile": float(centroid_outlier_quantile)
            if centroid_outlier_quantile is not None
            else None,
            "latent_key": "latent_space",
            "compression": compression if compression else None,
            "_obsSchemas": obs_schemas,
            "_continuousFields": obs_continuous_fields,
            "_categoricalFields": obs_categorical_fields,
        }
        obs_manifest_path.write_text(json.dumps(obs_manifest_payload), encoding="utf-8")

        total_fields = len(obs_continuous_fields) + len(obs_categorical_fields)
        print(
            f"✓ Wrote obs manifest ({total_fields} fields: {len(obs_continuous_fields)} continuous, "
            f"{len(obs_categorical_fields)} categorical) "
            f"to {published_path(obs_manifest_path)} "
            f"with binaries in {obs_binary_dirname}/"
        )

    expected_centroid_quantile = (
        float(centroid_outlier_quantile) if centroid_outlier_quantile is not None else None
    )
    expected_compression = compression if compression else None
    identity_obs_fields = _identity_obs_fields_from_compact_manifest(obs_manifest_payload)
    if (
        type(obs_manifest_payload.get("n_points")) is not int
        or obs_manifest_payload["n_points"] != int(n_cells)
        or obs_manifest_payload.get("centroid_outlier_quantile") != expected_centroid_quantile
        or obs_manifest_payload.get("latent_key") != "latent_space"
        or obs_manifest_payload.get("compression") != expected_compression
    ):
        raise ValueError(
            f"Observation manifest {obs_manifest_path} does not match the "
            "current export settings. Use force=True to replace it."
        )

    expected_identity_obs_fields: list[dict] = []
    for expected_kind in ("continuous", "category"):
        for field_info in obs_field_summaries:
            if field_info["kind"] != expected_kind:
                continue
            entry = {
                "key": field_info["key"],
                "kind": field_info["kind"],
            }
            if expected_kind == "category":
                entry["n_categories"] = field_info["category_count"]
            expected_identity_obs_fields.append(entry)
    if identity_obs_fields != expected_identity_obs_fields:
        raise ValueError(
            f"Observation manifest {obs_manifest_path} fields do not match "
            "the requested observation fields. Use force=True to replace it."
        )

    # Process gene expression if provided
    exported_var_field_count = 0
    if gene_expression is not None:
        if gene_expression_for_export is None:
            raise RuntimeError("Validated gene expression data is unavailable.")

        var_manifest_path = out_dir / var_manifest_filename
        if _output_path_is_writable(var_manifest_path, "var manifest"):
            var_binary_dir = out_dir / var_binary_dirname
            var_binary_dir.mkdir(parents=True, exist_ok=True)

            var_manifest_fields: list[list[Any]] = []

            for gene_id in tqdm.tqdm(genes_to_export, desc="Exporting genes"):
                gene_idx = gene_id_to_idx[gene_id]
                safe_gene_id = safe_gene_id_by_id[gene_id]

                if gene_expr_is_sparse:
                    sparse_expression = cast(
                        sparse.csc_matrix,
                        gene_expression_for_export,
                    )
                    col = sparse_expression.getcol(gene_idx).toarray().flatten()
                else:
                    col = gene_expression_for_export[:, gene_idx]

                gene_values = _require_finite_float32_array(
                    col,
                    label=f"Gene {gene_id!r} expression",
                )

                if gene_values.ndim != 1 or gene_values.shape[0] != n_cells:
                    raise ValueError(
                        f"Gene {gene_id!r} expression must have shape ({n_cells},), "
                        f"got {gene_values.shape}."
                    )

                # Apply quantization if requested
                if var_quantization is not None:
                    quantized, min_val, max_val, scale = _quantize_continuous(
                        gene_values,
                        bits=var_quantization,
                        field_name=gene_id,
                    )

                    if var_quantization == 8:
                        dtype_str = "uint8"
                        ext = "u8"
                    else:
                        dtype_str = "uint16"
                        ext = "u16"

                    value_fname = f"{safe_gene_id}.values.{ext}"
                    value_path = var_binary_dir / value_fname
                    _write_binary(value_path, quantized, compression)

                    manifest_path = f"{var_binary_dirname}/{value_fname}"
                    if compression:
                        manifest_path += ".gz"

                    # Compact format: [key, minValue, maxValue]
                    var_manifest_fields.append([gene_id, min_val, max_val])
                else:
                    # Full precision
                    value_fname = f"{safe_gene_id}.values.f32"
                    value_path = var_binary_dir / value_fname
                    _write_binary(value_path, gene_values, compression)

                    # Compact format: [key] for non-quantized
                    var_manifest_fields.append([gene_id])

            # Build compact manifest with schema
            gz_suffix = ".gz" if compression else ""
            if var_quantization is not None:
                ext = "u8" if var_quantization == 8 else "u16"
                dtype_str = "uint8" if var_quantization == 8 else "uint16"
                var_schema = {
                    "kind": "continuous",
                    "pathPattern": f"{var_binary_dirname}/{{key}}.values.{ext}{gz_suffix}",
                    "ext": ext,
                    "dtype": dtype_str,
                    "quantized": True,
                    "quantizationBits": var_quantization,
                }
            else:
                var_schema = {
                    "kind": "continuous",
                    "pathPattern": f"{var_binary_dirname}/{{key}}.values.f32{gz_suffix}",
                    "ext": "f32",
                    "dtype": "float32",
                    "quantized": False,
                }

            var_manifest_payload = {
                "_format": MANIFEST_FORMAT_VERSION,
                "n_points": int(n_cells),
                "var_gene_id_column": var_gene_id_column,
                "compression": compression if compression else None,
                "quantization": var_quantization,
                "_varSchema": var_schema,
                "fields": var_manifest_fields,
            }
            var_manifest_path.write_text(json.dumps(var_manifest_payload), encoding="utf-8")
            exported_var_field_count = len(var_manifest_fields)

            compression_info = f", gzip level {compression}" if compression else ""
            quant_info = f", {var_quantization}-bit quantized" if var_quantization else ""
            print(
                f"✓ Wrote var manifest ({len(var_manifest_fields)} genes{quant_info}{compression_info}) "
                f"to {published_path(var_manifest_path)}"
            )
    else:
        print("INFO: Gene expression was not requested; no var artifact was emitted.")

    # Process connectivity data if provided
    # GPU-optimized edge format for instanced rendering
    connectivity_meta: dict[str, int | str | None] = {
        "n_edges": None,
        "max_neighbors": None,
        "index_dtype": None,
    }
    if connectivities is not None:
        if connectivity_edges is None:
            raise RuntimeError("Validated connectivity edge pairs are unavailable.")
        connectivity_manifest_path = out_dir / CONNECTIVITY_MANIFEST_FILENAME
        connectivity_meta["n_edges"] = connectivity_edges.n_edges
        connectivity_meta["max_neighbors"] = connectivity_edges.max_neighbors
        connectivity_meta["index_dtype"] = connectivity_edges.index_dtype

        if _output_path_is_writable(
            connectivity_manifest_path,
            "connectivity manifest",
        ):
            connectivity_binary_dir = out_dir / CONNECTIVITY_BINARY_DIRNAME
            connectivity_binary_dir.mkdir(parents=True, exist_ok=True)

            # Write binary files (column-separated for better compression)
            sources_fname = "edges.src.bin"
            dests_fname = "edges.dst.bin"
            weights_fname = "edges.weights.f64.bin"
            sources_path = connectivity_binary_dir / sources_fname
            dests_path = connectivity_binary_dir / dests_fname
            weights_path = connectivity_binary_dir / weights_fname

            _write_binary(sources_path, connectivity_edges.sources, compression)
            _write_binary(
                dests_path,
                connectivity_edges.destinations,
                compression,
            )
            _write_binary(
                weights_path,
                connectivity_edges.weights,
                compression,
            )

            connectivity_manifest_payload = build_connectivity_manifest(
                n_cells=n_cells,
                n_edges=connectivity_edges.n_edges,
                max_neighbors=connectivity_edges.max_neighbors,
                index_bytes=connectivity_edges.index_bytes,
                index_dtype=connectivity_edges.index_dtype,
                compression=compression,
            )

            connectivity_manifest_path.write_text(
                json.dumps(connectivity_manifest_payload), encoding="utf-8"
            )

            print(
                f"✓ Wrote connectivity ({connectivity_edges.n_edges:,} edges, "
                f"max {connectivity_edges.max_neighbors} neighbors/cell, "
                f"{connectivity_edges.index_dtype}) "
                f"to {published_path(connectivity_binary_dir)}"
            )
    else:
        print("INFO: Connectivity was not requested; no connectivity artifact was emitted.")

    # =========================================================================
    # Generate dataset_identity.json (metadata for multi-dataset support)
    # =========================================================================
    identity_path = out_dir / "dataset_identity.json"

    from cellucid import __version__ as cellucid_version

    # Identity order and counts come from the exact emitted compact manifest.
    n_obs_fields = len(identity_obs_fields)
    n_categorical_fields = sum(field["kind"] == "category" for field in identity_obs_fields)
    n_continuous_fields = sum(field["kind"] == "continuous" for field in identity_obs_fields)

    # Build source info if provided
    source_info = None
    if source_name is not None:
        source_info = {"name": source_name}
        if source_url is not None:
            source_info["url"] = source_url
        if source_citation is not None:
            source_info["citation"] = source_citation

    # Build export settings
    export_settings = {
        "compression": compression if compression else None,
        "var_quantization": var_quantization,
        "obs_continuous_quantization": obs_continuous_quantization,
        "obs_categorical_dtype": obs_categorical_dtype,
    }

    # Build embeddings metadata
    gz_suffix = ".gz" if compression else ""
    embeddings_meta: dict[str, Any] = {
        "available_dimensions": available_dimensions,
        "default_dimension": default_dimension,
        "files": {},
    }
    for dim in available_dimensions:
        embeddings_meta["files"][f"{dim}d"] = f"points_{dim}d.bin{gz_suffix}"

    # Build identity payload
    identity_payload = {
        "version": 2,  # Bumped version for multi-dimensional support
        "id": dataset_id,
        "name": dataset_name,
        "description": dataset_description if dataset_description is not None else "",
        "created_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "cellucid_data_version": cellucid_version,
        "stats": {
            "n_cells": int(n_cells),
            "n_genes": int(exported_var_field_count),
            "n_obs_fields": int(n_obs_fields),
            "n_categorical_fields": int(n_categorical_fields),
            "n_continuous_fields": int(n_continuous_fields),
            "has_connectivity": connectivity_meta.get("n_edges") is not None,
            "n_edges": connectivity_meta.get("n_edges"),
        },
        "embeddings": embeddings_meta,
        "obs_fields": identity_obs_fields,
        "export_settings": export_settings,
    }

    if source_info:
        identity_payload["source"] = source_info

    if vector_fields_identity:
        identity_payload["vector_fields"] = vector_fields_identity

    identity_path.write_text(json.dumps(identity_payload, indent=2), encoding="utf-8")
    print(f"✓ Wrote dataset identity to {published_path(identity_path)}")


def prepare(
    latent_space: np.ndarray | sparse.spmatrix | None = None,
    obs: pd.DataFrame | None = None,
    var: pd.DataFrame | None = None,
    gene_expression: np.ndarray | sparse.spmatrix | None = None,
    var_gene_id_column: str | None = None,
    gene_identifiers: Sequence[str] | None = None,
    connectivities: np.ndarray | sparse.spmatrix | None = None,
    out_dir: Path | str = DEFAULT_EXPORT_DIR,
    obs_keys: Sequence[str] | None = None,
    centroid_outlier_quantile: float = 0.95,
    centroid_min_points: int = 10,
    force: bool = False,
    var_quantization: int | None = None,
    obs_continuous_quantization: int | None = None,
    compression: int | None = None,
    *,
    obs_categorical_dtype: Literal["uint8", "uint16"],
    dataset_name: str,
    dataset_id: str,
    dataset_description: str | None = None,
    source_name: str | None = None,
    source_url: str | None = None,
    source_citation: str | None = None,
    X_umap_1d: np.ndarray | None = None,
    X_umap_2d: np.ndarray | None = None,
    X_umap_3d: np.ndarray | None = None,
    vector_fields: dict[str, np.ndarray | sparse.spmatrix] | None = None,
    vector_field_default: str | None = None,
) -> None:
    """Build and atomically publish one complete canonical Cellucid export generation."""
    force = _require_native_boolean(force, label="force")
    centroid_min_points = _require_positive_native_integer(
        centroid_min_points,
        label="centroid_min_points",
    )
    dataset_description = _require_optional_native_string(
        dataset_description,
        label="dataset_description",
        allow_empty=True,
    )
    source_name = _require_optional_native_string(
        source_name,
        label="source_name",
        allow_empty=False,
    )
    source_url = _require_optional_native_string(
        source_url,
        label="source_url",
        allow_empty=False,
    )
    source_citation = _require_optional_native_string(
        source_citation,
        label="source_citation",
        allow_empty=False,
    )
    if source_name is None and (source_url is not None or source_citation is not None):
        raise ValueError(
            "source_name is required whenever source_url or source_citation is provided."
        )

    target_dir = Path(out_dir)
    if target_dir.name == "":
        raise ValueError("out_dir must name a child export directory.")
    target_dir.parent.mkdir(parents=True, exist_ok=True)
    with _exclusive_export_generation(target_dir):
        if target_dir.is_symlink():
            raise ValueError(f"out_dir must not be a symbolic link: {target_dir}")
        if target_dir.exists():
            if not target_dir.is_dir():
                raise NotADirectoryError(f"Export path exists but is not a directory: {target_dir}")
            if not force and any(target_dir.iterdir()):
                raise FileExistsError(
                    f"Refusing to replace non-empty export directory {target_dir}. "
                    "Pass force=True to publish a complete replacement generation."
                )

        staging_dir = Path(
            tempfile.mkdtemp(
                prefix=f".{target_dir.name}.cellucid-stage-",
                dir=target_dir.parent,
            )
        )
        try:
            _prepare_generation(
                latent_space=latent_space,
                obs=obs,
                var=var,
                gene_expression=gene_expression,
                var_gene_id_column=var_gene_id_column,
                gene_identifiers=gene_identifiers,
                connectivities=connectivities,
                out_dir=staging_dir,
                obs_keys=obs_keys,
                centroid_outlier_quantile=centroid_outlier_quantile,
                centroid_min_points=centroid_min_points,
                var_quantization=var_quantization,
                obs_continuous_quantization=obs_continuous_quantization,
                obs_categorical_dtype=obs_categorical_dtype,
                compression=compression,
                _published_out_dir=target_dir,
                dataset_name=dataset_name,
                dataset_id=dataset_id,
                dataset_description=dataset_description,
                source_name=source_name,
                source_url=source_url,
                source_citation=source_citation,
                X_umap_1d=X_umap_1d,
                X_umap_2d=X_umap_2d,
                X_umap_3d=X_umap_3d,
                vector_fields=vector_fields,
                vector_field_default=vector_field_default,
            )
            _publish_export_generation(staging_dir, target_dir)
        except BaseException:
            if staging_dir.exists():
                try:
                    shutil.rmtree(staging_dir)
                except BaseException as cleanup_error:
                    raise RuntimeError(
                        f"Failed to remove rejected staged generation {staging_dir}."
                    ) from cleanup_error
            raise


def generate_datasets_manifest(
    exports_dir: str | Path = DEFAULT_EXPORT_DIR,
    *,
    default_dataset: str,
) -> Path:
    """
    Validate an exports directory and atomically publish its datasets.json manifest.

    This utility helps maintain the datasets.json manifest file that the frontend uses
    to discover available demo datasets. Run this after adding or removing datasets.

    Parameters
    ----------
    exports_dir : Path or str
        Directory containing dataset subdirectories (default: exports/).
    default_dataset : str
        Exact ID of the dataset selected by default.

    Returns
    -------
    Path
        Path to the generated datasets.json file.

    Example
    -------
    >>> from cellucid.prepare_data import generate_datasets_manifest
    >>> generate_datasets_manifest("./exports", default_dataset="my_dataset")
    """
    exports_dir = Path(exports_dir)
    if not exports_dir.exists():
        raise FileNotFoundError(f"Exports directory not found: {exports_dir}")
    if not exports_dir.is_dir():
        raise NotADirectoryError(f"Exports path must be a directory: {exports_dir}")
    default_dataset = _require_nonempty_string(
        default_dataset,
        label="default_dataset",
    )
    if default_dataset != default_dataset.strip() or re.search(r"[\x00-\x1f\x7f]", default_dataset):
        raise ValueError(
            "default_dataset must be exact text without surrounding whitespace "
            "or control characters."
        )

    print(f"Scanning {exports_dir} for datasets...")

    datasets: list[dict[str, object]] = []
    dataset_ids: set[str] = set()
    directory_names: dict[str, str] = {}

    for subdir in sorted(exports_dir.iterdir()):
        if not subdir.is_dir():
            continue

        directory_name = subdir.name
        if (
            len(directory_name.encode("utf-8")) > 180
            or directory_name in {".", ".."}
            or re.fullmatch(r"[A-Za-z0-9_.-]+", directory_name, flags=re.ASCII) is None
            or directory_name.endswith(".")
            or directory_name.split(".", 1)[0].upper() in _WINDOWS_RESERVED_COMPONENTS
        ):
            raise ValueError(
                f"Dataset directory {directory_name!r} is not a portable URL/path component."
            )
        directory_collision_key = directory_name.casefold()
        if directory_collision_key in directory_names:
            raise ValueError(
                f"Dataset directories {directory_names[directory_collision_key]!r} and "
                f"{directory_name!r} collide case-insensitively."
            )
        directory_names[directory_collision_key] = directory_name

        identity_file = subdir / "dataset_identity.json"
        if not identity_file.is_file():
            raise ValueError(
                f"Dataset directory {directory_name!r} has no dataset_identity.json file."
            )

        try:
            identity = json.loads(identity_file.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            raise ValueError(f"{identity_file} must contain readable UTF-8 JSON.") from error
        if not isinstance(identity, dict):
            raise TypeError(f"{identity_file} must contain a JSON object.")
        if type(identity.get("version")) is not int or identity["version"] != 2:
            raise ValueError(f"{identity_file} version must be exactly 2.")

        dataset_id = _require_nonempty_string(
            identity.get("id"),
            label=f"{identity_file} dataset_id",
        )
        if dataset_id != dataset_id.strip() or re.search(r"[\x00-\x1f\x7f]", dataset_id):
            raise ValueError(
                f"{identity_file} dataset_id must be exact text without surrounding "
                "whitespace or control characters."
            )
        if dataset_id in dataset_ids:
            raise ValueError(f"duplicate dataset id {dataset_id!r}.")
        dataset_ids.add(dataset_id)

        dataset_name = _require_nonempty_string(
            identity.get("name"),
            label=f"{identity_file} dataset_name",
        )
        if dataset_name != dataset_name.strip():
            raise ValueError(
                f"{identity_file} dataset_name must not contain surrounding whitespace."
            )

        dataset_entry: dict[str, object] = {
            "id": dataset_id,
            "path": f"{directory_name}/",
            "name": dataset_name,
        }
        if "description" in identity:
            description = identity["description"]
            if not isinstance(description, str):
                raise TypeError(f"{identity_file} description must be a string.")
            dataset_entry["description"] = description

        stats = identity.get("stats")
        if not isinstance(stats, dict):
            raise TypeError(f"{identity_file} stats must be a JSON object.")
        for count_key in ("n_cells", "n_genes"):
            count = stats.get(count_key)
            if type(count) is not int or count < 0 or count > (1 << 53) - 1:
                raise ValueError(
                    f"{identity_file} stats.{count_key} must be a non-negative safe integer."
                )
            dataset_entry[count_key] = count

        datasets.append(dataset_entry)
        print(f"  ✓ Found dataset: {dataset_entry['name']} ({dataset_entry['id']})")

    if not datasets:
        raise ValueError(f"Exports directory contains no datasets: {exports_dir}")
    if default_dataset not in dataset_ids:
        raise ValueError(
            f"default_dataset {default_dataset!r} is missing from the validated datasets."
        )

    manifest = {"version": 1, "default": default_dataset, "datasets": datasets}
    manifest_bytes = json.dumps(manifest, indent=2).encode("utf-8")
    manifest_path = exports_dir / "datasets.json"

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=".datasets.json.",
        suffix=".tmp",
        dir=exports_dir,
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(manifest_bytes)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary_path, manifest_path)
    except BaseException:
        if temporary_path.exists():
            temporary_path.unlink()
        raise

    print(f"✓ Wrote datasets manifest with {len(datasets)} datasets to {manifest_path}")
    print(f"  Default dataset: {default_dataset}")

    return manifest_path
