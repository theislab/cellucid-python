"""
Building one complete export generation inside a staging directory.

This is the writer proper. It resolves and validates every identity -- fields,
genes, vector fields, category labels -- and proves every continuous payload is
encodable before the first output path is created, so a run either fails
naming every defect at once or produces a whole generation. Nothing here
publishes: the staged directory is handed back to the transaction.
"""

import math
from collections.abc import Sequence
from numbers import Real
from pathlib import Path
from typing import Any, Literal, cast

import numpy as np
import pandas as pd
import tqdm
from scipy import sparse

from .._console import console_print
from .._contracts import (
    JsonScalar,
    _require_dataset_id,
    _require_dataset_name,
    _require_field_identities,
    _require_string_identifiers,
    _require_unique_identifiers,
)
from ..connectivity_contract import (
    CONNECTIVITY_BINARY_DIRNAME,
    CONNECTIVITY_MANIFEST_FILENAME,
    ConnectivityEdgePairs,
    build_connectivity_manifest,
    validate_connectivity_edges,
)
from ..vector_fields import scale_vector_field, validate_vector_fields
from ._arrays import (
    _normalize_finite_float32_embedding,
    _require_continuous_obs_values,
    _require_finite_embedding_source,
    _require_finite_float32_array,
)
from ._binary_io import (
    _final_binary_path,
    _gzip_suffix,
    _output_path_is_writable,
    _write_binary,
)
from ._categories import _declared_categorical_storage, _json_category_values
from ._centroids import _compute_centroids_for_all_dimensions, _compute_latent_space_quantiles
from ._manifest import (
    DEFAULT_OBS_DIRNAME,
    DEFAULT_VAR_DIRNAME,
    MANIFEST_FORMAT_VERSION,
    _declared_obs_payload_paths,
    _declared_var_payload_paths,
    _gene_names_from_compact_manifest,
    _identity_obs_fields_from_compact_manifest,
    _require_declared_payloads_on_disk,
    _require_dense_payload_indices,
    _write_manifest,
)
from ._quantization import _quantize_continuous, _quantize_nullable_outlier_quantiles


def _gene_expression_column(
    matrix: np.ndarray | sparse.csc_matrix,
    *,
    is_sparse: bool,
    gene_index: int,
    gene_id: str,
    n_cells: int,
) -> np.ndarray:
    """Materialize exactly the float32 column that one gene would be exported as.

    The pre-flight finiteness scan and the export loop must agree to the last
    bit, because both the rejection and the published bounds depend on float32
    rounding: a source range that only exists in float64 collapses to one
    exported value, which is exported as a constant field. Both callers go
    through this function so neither can drift from the other.
    """
    if is_sparse:
        column = cast(sparse.csc_matrix, matrix).getcol(gene_index).toarray().flatten()
    else:
        column = cast(np.ndarray, matrix)[:, gene_index]

    values = _require_finite_float32_array(
        column,
        label=f"Gene {gene_id!r} expression",
    )
    if values.ndim != 1 or values.shape[0] != n_cells:
        raise ValueError(
            f"Gene {gene_id!r} expression must have shape ({n_cells},), got {values.shape}."
        )
    return values


def _all_missing_outlier_defect(values: np.ndarray, *, centroid_min_points: int) -> str | None:
    """Name why one generated outlier-quantile set has nothing to encode, or None.

    A set whose every quantile is missing has no value anywhere -- not one
    constant, none at all -- so there is no bound to publish and nothing for a
    reader to decode. That is a ``centroid_min_points`` that no category
    reaches, which is a setting the caller can change, and it is the only
    remaining reason a continuous payload cannot be encoded.
    """
    if np.isnan(values).all():
        return (
            "no category holds at least centroid_min_points="
            f"{centroid_min_points} cells, so every generated quantile is missing"
        )
    return None


def _require_encodable_outlier_quantiles(
    *,
    outlier_defects: list[tuple[str, str]],
    obs_continuous_quantization: int | None,
) -> None:
    """Reject every unencodable outlier-quantile set at once, before any output.

    The check is per field, so the export loop would otherwise fail on the
    first offender, after the fields ahead of it were already written. One run
    names every affected field instead of one run per field.
    """
    if not outlier_defects:
        return

    listing = "\n".join(f"    {name!r}: {reason}" for name, reason in outlier_defects)
    raise ValueError(
        f"{len(outlier_defects)} generated categorical outlier quantile set(s) cannot be "
        "encoded: a set with no quantile at all has no value to publish.\n"
        f"  {len(outlier_defects)} set(s) with "
        f"obs_continuous_quantization={obs_continuous_quantization}:\n"
        f"{listing}\n"
        "  Fix: lower centroid_min_points so a category qualifies, drop the field "
        "from obs_keys=, or pass obs_continuous_quantization=None to export the "
        "generated quantiles as full-precision float32."
    )


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


def _prepare_generation(
    latent_space: np.ndarray | sparse.spmatrix | None = None,
    obs: pd.DataFrame | None = None,
    var: pd.DataFrame | None = None,
    gene_expression: np.ndarray | sparse.spmatrix | None = None,
    var_gene_id_column: str | None = None,
    gene_identifiers: Sequence[str] | None = None,
    connectivities: np.ndarray | sparse.spmatrix | None = None,
    obs_keys: Sequence[str] | None = None,
    centroid_outlier_quantile: float | None = 0.95,
    centroid_min_points: int = 10,
    # Optimization parameters (all disabled by default)
    var_quantization: int | None = None,
    obs_continuous_quantization: int | None = None,
    compression: int | None = None,
    # Dataset metadata parameters (for dataset_identity.json)
    *,
    out_dir: Path,
    _published_out_dir: Path,
    obs_categorical_dtype: Literal["uint8", "uint16"],
    dataset_name: str,
    dataset_id: str,
    created_at: str,
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
        Codes and bounds derive from the viewer-visible float32 values. A gene
        whose values are all one float32 value -- including a source range that
        collapses to one under float32 rounding, and including a gene detected
        in no exported cell -- is published as a constant field: equal bounds
        and every code 0, which decodes back to that exact value.
    obs_continuous_quantization : int or None
        Bits for continuous obs field quantization (8, 16, or None for full float32).
        It follows the same viewer-visible float32 domain as var quantization,
        including the constant-field encoding, and it also governs the
        generated categorical outlier quantiles.
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
        1D embedding coordinates, shape (n_cells, 1). A plain ``(n_cells,)``
        vector is accepted as the single column it declares. Stored as
        points_1d.bin.
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
        real values representable as float32. A ``<field>_umap_1d`` key also
        accepts a plain ``(n_cells,)`` vector as the single column it declares.
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
        Cell metadata, shape (n_cells, n_obs_columns). Every string category
        label must be text the viewer can show as written: non-empty, free of
        control characters and of the zero-width characters U+200B, U+2060, and
        U+FEFF, and without leading or trailing whitespace. Two labels that a
        whitespace-collapsing renderer draws identically are rejected too. No
        label is ever trimmed: trimming rewrites an annotation the caller did
        not ask to change and can merge two categories into one.
    var : pd.DataFrame, optional
        Gene/feature metadata. Required if gene_expression is provided.
    gene_expression : np.ndarray or sparse matrix, optional
        Gene expression matrix, shape (n_cells, n_genes).
    var_gene_id_column : str or None
        Exact non-empty column name in var containing string gene names. None
        uses string names from var.index. Whatever this selects is recorded
        faithfully: Cellucid performs no symbol lookup and ships no mapping, so
        the caller decides what a gene is called. Every identifier in the
        selected axis must be a non-empty string and must be distinct, because
        ``gene_identifiers`` addresses var rows by identifier and a repeated one
        names no single row.
    gene_identifiers : sequence of str, optional
        Unique gene identifiers to export. Every identifier must be a non-empty
        string present in var. If None, all genes are exported. Payload files
        are named by integer index, so an identifier is never a path; the
        exported ones must additionally be text the viewer can draw exactly as
        stored, exactly as ``obs_keys`` narrows which observation keys are held
        to that rule. A var row this argument leaves out reaches no manifest and
        is not checked against it.
    connectivities : array or sparse matrix, optional
        Exact weighted undirected graph, shape ``(n_cells, n_cells)``. Values
        must be finite and non-negative, the topology and weights must be
        exactly symmetric, and the diagonal must be zero. Sparse inputs must
        not contain duplicate coordinates or stored zero entries.
    out_dir : Path
        Staging directory that receives the candidate generation. It is always
        supplied by :func:`prepare` and must be absent or empty.
    obs_keys : sequence of str or None
        Which obs columns to export. If None, all columns are exported.
    centroid_outlier_quantile : float or None
        Quantile of distances to keep as inliers when computing centroids. Pass
        None to skip the computation entirely: every categorical field is then
        published with empty centroid lists.
    centroid_min_points : int
        Minimum number of points in a category to compute a centroid.

    Dataset Metadata Parameters
    ---------------------------
    dataset_name : str
        Exact non-empty human-readable name without surrounding whitespace or
        control characters. Unicode is preserved.
    dataset_description : str, optional
        Description of the dataset. May be empty.
    dataset_id : str
        Portable 1-180 byte ASCII identifier: begin with a letter or digit,
        then use only letters, digits, ``.``, ``_``, or ``-``. It cannot end
        in ``.`` or use a reserved Windows device name.
    source_name : str, optional
        Name of the data source (e.g., "HLCA Consortium").
    source_url : str, optional
        URL to the data source.
    source_citation : str, optional
        Citation text for the data source.

    ``dataset_name``, ``dataset_description``, and the three ``source_*``
    strings are shown to the reader verbatim, so each obeys the same rule as a
    string category label: no control character, none of the zero-width
    characters U+200B, U+2060, or U+FEFF, and no leading or trailing whitespace
    of any kind, U+00A0 NO-BREAK SPACE included.
    """
    _validate_prepare_codec_options(
        compression=compression,
        var_quantization=var_quantization,
        obs_continuous_quantization=obs_continuous_quantization,
        obs_categorical_dtype=obs_categorical_dtype,
        centroid_outlier_quantile=centroid_outlier_quantile,
    )
    dataset_name = _require_dataset_name(
        dataset_name,
        label="dataset_name",
    )
    dataset_id = _require_dataset_id(
        dataset_id,
        label="dataset_id",
    )

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
        embeddings[1] = _require_finite_embedding_source(
            X_umap_1d,
            label="X_umap_1d",
        )
    if X_umap_2d is not None:
        embeddings[2] = _require_finite_embedding_source(
            X_umap_2d,
            label="X_umap_2d",
        )
    if X_umap_3d is not None:
        embeddings[3] = _require_finite_embedding_source(
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
        if arr.ndim == 1 and dim == 1:
            # X_umap_1d declares one coordinate per cell, so a plain length-n
            # vector already is the (n, 1) column it declares. AnnData performs
            # exactly this reshape when such an array is put in obsm, so the
            # AnnData path already accepts it.
            arr = arr.reshape(-1, 1)
            embeddings[dim] = arr
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
    console_print("=" * 60)
    console_print("Export Settings:")
    console_print(f"  Output directory: {_published_out_dir}")
    console_print(
        f"  Compression: {'gzip level ' + str(compression) if compression else 'disabled'}"
    )
    console_print(
        f"  Var (gene) quantization: {str(var_quantization) + '-bit' if var_quantization else 'disabled (float32)'}"
    )
    console_print(
        f"  Obs continuous quantization: {str(obs_continuous_quantization) + '-bit' if obs_continuous_quantization else 'disabled (float32)'}"
    )
    console_print(f"  Obs categorical dtype: {obs_categorical_dtype}")
    console_print(f"  Available dimensions: {available_dimensions}")
    console_print(f"  Default dimension: {default_dimension}D")
    console_print("  Coordinate normalization (per-dimension, aspect-ratio preserved):")
    for dim in sorted(normalization_info.keys()):
        info = normalization_info[dim]
        console_print(f"    {dim}D: range {info['original_range']:.2f} → [-1, 1]")
    console_print("=" * 60)

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
    if not isinstance(obs, pd.DataFrame):
        raise TypeError(f"obs must be a pandas DataFrame, got {type(obs).__name__}.")
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
    obs_keys = _require_field_identities(
        obs_keys,
        label="Observation field",
    )
    # obs/ holds both continuous and categorical payloads, so one index space
    # spans both manifest arrays and follows the order the fields were selected.
    obs_payload_index_by_key = {key: index for index, key in enumerate(obs_keys)}

    validated_continuous_obs: dict[str, np.ndarray] = {}
    validated_categorical_obs: dict[str, tuple[list[JsonScalar], np.ndarray]] = {}
    generated_outlier_quantiles: dict[str, np.ndarray] = {}
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
        codes = categorical.cat.codes.to_numpy(dtype=np.int32)  # -1 for NaN
        if codes.shape[0] != n_cells:
            raise ValueError(f"Length mismatch for obs['{key}']: {codes.shape[0]} vs {n_cells}")
        validated_categorical_obs[key] = (categories, codes)
        # Outlier quantiles are generated here, not in the export loop, so the
        # quantizability of this generated payload is known before any write.
        generated_outlier_quantiles[key] = _compute_latent_space_quantiles(
            latent=latent,
            codes=codes,
            categories=categories,
            min_points=centroid_min_points,
        )
        n_categories = len(categories)
        dtype_str = _declared_categorical_storage(
            obs_categorical_dtype,
            n_categories,
            field_key=key,
        ).dtype_str
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
    gene_payload_index_by_id: dict[str, int] = {}
    gene_expr_is_sparse = False
    gene_expression_for_export: np.ndarray | sparse.csc_matrix | None = None
    if gene_expression is not None:
        if var is None:
            raise ValueError("var DataFrame must be provided when gene_expression is given.")
        if not isinstance(var, pd.DataFrame):
            raise TypeError(f"var must be a pandas DataFrame, got {type(var).__name__}.")

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
        # Uniqueness over the whole var, display text over the exported subset.
        # Every var row is addressable through gene_identifiers=, so a repeated
        # identifier anywhere in var makes that lookup ambiguous and is rejected
        # here. Being drawable is a property of a name the viewer shows, so it
        # is checked on genes_to_export below: an unexported gene reaches no
        # manifest, exactly as an obs column left out of obs_keys= does.
        _require_unique_identifiers(all_gene_ids, label="Gene")
        gene_id_to_idx = {gene_id: index for index, gene_id in enumerate(all_gene_ids)}

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
            _require_unique_identifiers(genes_to_export, label="Requested gene")
            missing_genes = [
                gene_id for gene_id in genes_to_export if gene_id not in gene_id_to_idx
            ]
            if missing_genes:
                raise KeyError(
                    f"gene_identifiers contain identifiers not present in var: {missing_genes}."
                )

        genes_to_export = _require_field_identities(
            genes_to_export,
            label="Gene",
        )
        gene_payload_index_by_id = {
            gene_id: index for index, gene_id in enumerate(genes_to_export)
        }

    validated_vector_fields = validate_vector_fields(
        vector_fields,
        n_cells=n_cells,
        available_dimensions=available_dimensions,
        vector_field_default=vector_field_default,
    )
    scaled_vector_fields: dict[str, dict[int, np.ndarray]] = {}
    vector_field_payload_index_by_id: dict[str, int] = {}
    vector_fields_identity: dict[str, object] | None = None
    if validated_vector_fields.fields:
        fields_metadata: dict[str, dict[str, object]] = {}
        gzip_suffix = _gzip_suffix(compression)
        vector_field_ids = _require_field_identities(
            list(validated_vector_fields.fields),
            label="Vector field",
        )
        vector_field_payload_index_by_id = {
            field_id: index for index, field_id in enumerate(vector_field_ids)
        }
        for field_id, vectors_by_dimension in validated_vector_fields.fields.items():
            dimensions = sorted(vectors_by_dimension)
            payload_index = vector_field_payload_index_by_id[field_id]
            scaled_vector_fields[field_id] = {}
            for dimension, vectors in vectors_by_dimension.items():
                scale_factor = normalization_info[dimension]["scale_factor"]
                scaled_vector_fields[field_id][dimension] = scale_vector_field(
                    vectors,
                    scale_factor=scale_factor,
                    label=f"Vector field {field_id!r} {dimension}D",
                )

            # Key order is the one the output format specification prints and
            # the one cellucid-r emits, so both writers produce the same file.
            fields_metadata[field_id] = {
                "label": field_id,
                "available_dimensions": dimensions,
                "default_dimension": max(dimensions),
                "files": {
                    f"{dimension}d": (f"vectors/{payload_index}_{dimension}d.bin{gzip_suffix}")
                    for dimension in dimensions
                },
                "basis": "umap",
            }

        vector_fields_identity = {
            "default_field": validated_vector_fields.default_field,
            "fields": fields_metadata,
        }

    # =========================================================================
    # ENCODABILITY PRE-FLIGHT
    # =========================================================================
    # Every continuous payload is checked here, before the first output path
    # exists, so one run reports every field that cannot be encoded instead of
    # failing on each one in turn partway through the export. A constant field
    # is not one of them: it has its own encoding (equal bounds, every code 0),
    # so the only payload left with nothing to publish is a generated
    # outlier-quantile set in which every quantile is missing.
    outlier_range_defects: list[tuple[str, str]] = []

    if obs_continuous_quantization is not None:
        for key, quantiles in generated_outlier_quantiles.items():
            defect = _all_missing_outlier_defect(
                quantiles,
                centroid_min_points=centroid_min_points,
            )
            if defect is not None:
                outlier_range_defects.append((key, defect))

    if gene_expression_for_export is not None:
        # Materializing every column here is what makes a non-finite value in
        # any gene fail before the first output path is created rather than
        # partway through the export.
        for gene_id in genes_to_export:
            _gene_expression_column(
                gene_expression_for_export,
                is_sparse=gene_expr_is_sparse,
                gene_index=gene_id_to_idx[gene_id],
                gene_id=gene_id,
                n_cells=n_cells,
            )

    _require_encodable_outlier_quantiles(
        outlier_defects=outlier_range_defects,
        obs_continuous_quantization=obs_continuous_quantization,
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    obs_binary_dir.mkdir(parents=True, exist_ok=True)

    # The export root is reconciled against this register once the whole
    # generation exists, so every artifact created directly beneath it is named
    # here as it is created. The point payloads are the one exception: they are
    # declared by path in dataset_identity.json, and it is that declaration --
    # not the name the write chose -- that has to match what is on disk.
    declared_root_files: set[str] = set()
    declared_root_directories: set[str] = {obs_binary_dirname}

    # =========================================================================
    # SAVE DIMENSIONAL EMBEDDING FILES
    # =========================================================================
    # This cast is the one and only float32 rounding a coordinate gets.
    for dim, arr in embeddings.items():
        dim_filename = f"points_{dim}d.bin"
        dim_path = out_dir / dim_filename
        check_path = _final_binary_path(dim_path, compression)
        if _output_path_is_writable(check_path, check_path.name):
            actual_path = _write_binary(dim_path, arr.astype(np.float32), compression)
            suffix = " (gzip)" if compression else ""
            console_print(
                f"✓ Wrote {dim}D positions ({arr.shape[0]:,} cells × {dim} dims) "
                f"to {published_path(actual_path)}{suffix}"
            )

    # =========================================================================
    # SAVE VECTOR FIELDS (OPTIONAL)
    # =========================================================================
    if scaled_vector_fields:
        vectors_dir = out_dir / "vectors"
        vectors_dir.mkdir(parents=True, exist_ok=True)
        declared_root_directories.add("vectors")
        for field_id, vectors_by_dimension in scaled_vector_fields.items():
            payload_index = vector_field_payload_index_by_id[field_id]
            for dimension, vectors in vectors_by_dimension.items():
                filename = f"{payload_index}_{dimension}d.bin"
                path = vectors_dir / filename
                check_path = _final_binary_path(path, compression)
                if _output_path_is_writable(check_path, check_path.name):
                    actual_path = _write_binary(path, vectors, compression)
                    suffix = " (gzip)" if compression else ""
                    console_print(
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

        for key in obs_keys:
            s = obs[key]
            payload_index = obs_payload_index_by_key[key]

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

                    value_path = obs_binary_dir / f"{payload_index}.values.{ext}"
                    _write_binary(value_path, quantized, compression)

                    # Compact format: [index, key, minValue, maxValue]
                    obs_continuous_fields.append([payload_index, key, min_val, max_val])
                    if not continuous_dtype_info:
                        continuous_dtype_info["ext"] = ext
                        continuous_dtype_info["dtype"] = dtype_str
                        continuous_dtype_info["quantized"] = True
                        continuous_dtype_info["quantizationBits"] = obs_continuous_quantization
                else:
                    # Full precision
                    value_path = obs_binary_dir / f"{payload_index}.values.f32"
                    _write_binary(value_path, values, compression)

                    # Compact format: [index, key]
                    obs_continuous_fields.append([payload_index, key])
                    if not continuous_dtype_info:
                        continuous_dtype_info["ext"] = "f32"
                        continuous_dtype_info["dtype"] = "float32"
                        continuous_dtype_info["quantized"] = False

            else:
                # Categorical: identity and generated quantiles were resolved
                # before any output path existed.
                categories, codes = validated_categorical_obs[key]
                n_categories = len(categories)

                storage = _declared_categorical_storage(
                    obs_categorical_dtype,
                    n_categories,
                    field_key=key,
                )
                dtype = storage.dtype
                dtype_str = storage.dtype_str
                missing_value = storage.missing_value

                codes_typed = np.full(n_cells, missing_value, dtype=dtype)
                valid_mask = codes >= 0
                codes_typed[valid_mask] = codes[valid_mask].astype(dtype)

                codes_fname = f"{payload_index}.codes.{storage.codes_extension}"

                codes_path = obs_binary_dir / codes_fname
                _write_binary(codes_path, codes_typed, compression)

                # Compute centroids for all available dimensions
                if centroid_outlier_quantile is None:
                    centroids_by_dim: dict[int, list[dict]] = {dim: [] for dim in embeddings}
                else:
                    centroids_by_dim = _compute_centroids_for_all_dimensions(
                        embeddings,
                        codes,
                        categories,
                        outlier_quantile=centroid_outlier_quantile,
                        min_points=centroid_min_points,
                    )

                outlier_quantiles = generated_outlier_quantiles[key]

                # Quantize outlier quantiles (they're always 0-1)
                if obs_continuous_quantization is not None:
                    oq_quantized, oq_min, oq_max, oq_scale = _quantize_nullable_outlier_quantiles(
                        outlier_quantiles,
                        bits=obs_continuous_quantization,
                        field_name=f"{key}_outliers",
                    )

                    if obs_continuous_quantization == 8:
                        oq_dtype_str = "uint8"
                        oq_ext = "u8"
                    else:
                        oq_dtype_str = "uint16"
                        oq_ext = "u16"

                    outlier_path = obs_binary_dir / f"{payload_index}.outliers.{oq_ext}"
                    _write_binary(outlier_path, oq_quantized, compression)

                    # Compact format: [index, key, categories, codesDtype,
                    # codesMissingValue, centroidsByDim, outlierMinValue, outlierMaxValue]
                    # centroidsByDim is a dict keyed by dimension: {"1": [...], "2": [...], "3": [...]}
                    centroids_serializable = {
                        str(dim): cents for dim, cents in centroids_by_dim.items()
                    }
                    obs_categorical_fields.append(
                        [
                            payload_index,
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
                        categorical_dtype_info["codesExt"] = storage.codes_extension
                        categorical_dtype_info["outlierExt"] = oq_ext
                        categorical_dtype_info["outlierDtype"] = oq_dtype_str
                        categorical_dtype_info["outlierQuantized"] = True
                else:
                    # Full precision outliers
                    outlier_path = obs_binary_dir / f"{payload_index}.outliers.f32"
                    _write_binary(outlier_path, outlier_quantiles.astype(np.float32), compression)

                    # Compact format: [index, key, categories, codesDtype,
                    # codesMissingValue, centroidsByDim]
                    # centroidsByDim is a dict keyed by dimension: {"1": [...], "2": [...], "3": [...]}
                    centroids_serializable = {
                        str(dim): cents for dim, cents in centroids_by_dim.items()
                    }
                    obs_categorical_fields.append(
                        [
                            payload_index,
                            key,
                            categories,
                            dtype_str,
                            int(missing_value),
                            centroids_serializable,
                        ]
                    )
                    if not categorical_dtype_info:
                        categorical_dtype_info["codesExt"] = storage.codes_extension
                        categorical_dtype_info["outlierExt"] = "f32"
                        categorical_dtype_info["outlierDtype"] = "float32"
                        categorical_dtype_info["outlierQuantized"] = False

        # Build compact manifest with schemas
        gz_suffix = _gzip_suffix(compression)

        obs_schemas = {}
        if continuous_dtype_info:
            obs_schemas["continuous"] = {
                "pathPattern": f"{obs_binary_dirname}/{{index}}.values.{continuous_dtype_info['ext']}{gz_suffix}",
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
                "codesPathPattern": f"{obs_binary_dirname}/{{index}}.codes.{{ext}}{gz_suffix}",
                "outlierPathPattern": f"{obs_binary_dirname}/{{index}}.outliers.{categorical_dtype_info['outlierExt']}{gz_suffix}",
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
        _write_manifest(obs_manifest_path, obs_manifest_payload)
        declared_root_files.add(obs_manifest_filename)

        total_fields = len(obs_continuous_fields) + len(obs_categorical_fields)
        console_print(
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
    _require_declared_payloads_on_disk(
        out_dir,
        directory_name=obs_binary_dirname,
        declared=_declared_obs_payload_paths(obs_manifest_payload),
        axis="Observation",
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
            declared_root_directories.add(var_binary_dirname)

            var_manifest_fields: list[list[Any]] = []

            for gene_id in tqdm.tqdm(genes_to_export, desc="Exporting genes"):
                payload_index = gene_payload_index_by_id[gene_id]

                gene_values = _gene_expression_column(
                    gene_expression_for_export,
                    is_sparse=gene_expr_is_sparse,
                    gene_index=gene_id_to_idx[gene_id],
                    gene_id=gene_id,
                    n_cells=n_cells,
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

                    value_path = var_binary_dir / f"{payload_index}.values.{ext}"
                    _write_binary(value_path, quantized, compression)

                    # Compact format: [index, name, minValue, maxValue]
                    var_manifest_fields.append([payload_index, gene_id, min_val, max_val])
                else:
                    # Full precision
                    value_path = var_binary_dir / f"{payload_index}.values.f32"
                    _write_binary(value_path, gene_values, compression)

                    # Compact format: [index, name] for non-quantized
                    var_manifest_fields.append([payload_index, gene_id])

            # Build compact manifest with schema
            gz_suffix = _gzip_suffix(compression)
            if var_quantization is not None:
                ext = "u8" if var_quantization == 8 else "u16"
                dtype_str = "uint8" if var_quantization == 8 else "uint16"
                var_schema = {
                    "kind": "continuous",
                    "pathPattern": f"{var_binary_dirname}/{{index}}.values.{ext}{gz_suffix}",
                    "ext": ext,
                    "dtype": dtype_str,
                    "quantized": True,
                    "quantizationBits": var_quantization,
                }
            else:
                var_schema = {
                    "kind": "continuous",
                    "pathPattern": f"{var_binary_dirname}/{{index}}.values.f32{gz_suffix}",
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
            _write_manifest(var_manifest_path, var_manifest_payload)
            declared_root_files.add(var_manifest_filename)
            exported_var_field_count = len(var_manifest_fields)

            exported_gene_names = _gene_names_from_compact_manifest(var_manifest_payload)
            if exported_gene_names != genes_to_export:
                raise ValueError(
                    f"Gene manifest {var_manifest_path} names do not match the "
                    "requested genes. Use force=True to replace it."
                )
            _require_declared_payloads_on_disk(
                out_dir,
                directory_name=var_binary_dirname,
                declared=_declared_var_payload_paths(var_manifest_payload),
                axis="Gene",
            )

            compression_info = f", gzip level {compression}" if compression else ""
            quant_info = f", {var_quantization}-bit quantized" if var_quantization else ""
            console_print(
                f"✓ Wrote var manifest ({len(var_manifest_fields)} genes{quant_info}{compression_info}) "
                f"to {published_path(var_manifest_path)}"
            )
    else:
        console_print("INFO: Gene expression was not requested; no var artifact was emitted.")

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
            declared_root_directories.add(CONNECTIVITY_BINARY_DIRNAME)

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

            _write_manifest(connectivity_manifest_path, connectivity_manifest_payload)
            declared_root_files.add(CONNECTIVITY_MANIFEST_FILENAME)

            declared_connectivity_paths: set[str] = set()
            for path_key in ("sourcesPath", "destinationsPath", "weightsPath"):
                declared_path = connectivity_manifest_payload[path_key]
                if not isinstance(declared_path, str) or not declared_path:
                    raise RuntimeError(
                        f"connectivity_manifest.json {path_key} must be one payload path."
                    )
                declared_connectivity_paths.add(declared_path)
            _require_declared_payloads_on_disk(
                out_dir,
                directory_name=CONNECTIVITY_BINARY_DIRNAME,
                declared=declared_connectivity_paths,
                axis="Connectivity",
            )

            console_print(
                f"✓ Wrote connectivity ({connectivity_edges.n_edges:,} edges, "
                f"max {connectivity_edges.max_neighbors} neighbors/cell, "
                f"{connectivity_edges.index_dtype}) "
                f"to {published_path(connectivity_binary_dir)}"
            )
    else:
        console_print("INFO: Connectivity was not requested; no connectivity artifact was emitted.")

    if vector_fields_identity is not None:
        declared_vector_paths: set[str] = set()
        for field_metadata in cast(dict[str, Any], vector_fields_identity["fields"]).values():
            declared_vector_paths.update(field_metadata["files"].values())
        _require_dense_payload_indices(
            list(vector_field_payload_index_by_id.values()),
            axis="Vector field",
        )
        _require_declared_payloads_on_disk(
            out_dir,
            directory_name="vectors",
            declared=declared_vector_paths,
            axis="Vector field",
        )

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
    gz_suffix = _gzip_suffix(compression)
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
        "created_at": created_at,
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

    _write_manifest(identity_path, identity_payload, indent=2)
    declared_root_files.add(identity_path.name)
    console_print(f"✓ Wrote dataset identity to {published_path(identity_path)}")

    # Each axis directory has now been reconciled against the manifest that
    # declares it, but the point payloads sit at the export root and are
    # declared in dataset_identity.json, so nothing above compared the
    # coordinate files the viewer is told to fetch against the files that were
    # written. The declared name and the written name are two expressions of
    # the compression setting, and a generation whose points disagree publishes
    # successfully and then fails in the browser with no coordinates. The root
    # is the last surface where a declaration and a file can disagree, so it is
    # reconciled whole: the manifests and the payload directories this call
    # created are named here too, and anything else that reached the root is
    # refused rather than published.
    for dimension_key, declared_point_path in embeddings_meta["files"].items():
        if not isinstance(declared_point_path, str) or not declared_point_path:
            raise RuntimeError(
                f"dataset_identity.json embeddings.files.{dimension_key} must be one payload path."
            )
        declared_root_files.add(declared_point_path)
    _require_declared_payloads_on_disk(
        out_dir,
        directory_name=None,
        declared=declared_root_files,
        axis="Export",
        declared_directories=declared_root_directories,
    )
