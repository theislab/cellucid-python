"""
AnnData Adapter for Cellucid

Provides a server-side adapter that reads AnnData objects (in memory, backed h5ad,
or zarr stores) and serves data in the same format that the Cellucid web viewer
expects. This allows users to visualize AnnData directly without running
prepare first.

Key features:
- Works with in-memory AnnData, backed h5ad files, and eagerly loaded zarr stores
- Handles sparse matrices (CSR/CSC) transparently with automatic format conversion
- Reads UMAP dimensions only from explicit obsm keys (X_umap_1d, X_umap_2d, X_umap_3d)
- Computes centroids and outlier quantiles on-demand with caching
- Lazy loading: gene expression and obs data are loaded on-demand
- No quantization (full float32 precision, gzip compression for network transfer)

Lazy Loading Behavior:
----------------------
- h5ad (backed='r'): Uses HDF5 memory-mapping. Arrays are loaded only when accessed.
  Gene expression columns are fetched individually, minimizing memory usage.
- zarr: ``anndata.read_zarr`` materializes the AnnData arrays in memory.
- in-memory: All data is already in RAM. No lazy loading.

Memory Management:
-----------------
The adapter maintains several caches for performance:
- Embedding cache: Normalized UMAP coordinates (one per dimension)
- Centroid cache: Computed label centroids
- Gene expression LRU cache: Recently accessed gene columns (max 100 genes)
- CSC cache: For CSR matrices, a CSC copy is created for efficient column access

IMPORTANT: Always use the context manager or call close() when done:
    with AnnDataAdapter.from_file(
        "data.h5ad",
        dataset_name="Example",
        dataset_id="example",
    ) as adapter:
        # use adapter
    # resources are automatically released

Usage:
    from cellucid.anndata_adapter import AnnDataAdapter

    # From h5ad file (lazy loading via backed mode)
    with AnnDataAdapter.from_file(
        "/path/to/data.h5ad",
        dataset_name="Example",
        dataset_id="example",
    ) as adapter:
        manifest = adapter.get_obs_manifest()

    # From a zarr store (materialized by anndata.read_zarr)
    adapter = AnnDataAdapter.from_file(
        "/path/to/data.zarr",
        dataset_name="Example",
        dataset_id="example",
    )

    # From in-memory AnnData
    adapter = AnnDataAdapter(
        adata,
        dataset_name="Example",
        dataset_id="example",
    )

    # Use with server
    from cellucid import show_anndata
    show_anndata(
        adata,
        dataset_name="Example",
        dataset_id="example",
    )
"""

from __future__ import annotations

import json
import logging
import math
import re
from collections import OrderedDict
from datetime import UTC, datetime
from numbers import Integral, Real
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, cast

import numpy as np
import pandas as pd
from scipy import sparse

from ._compression import deterministic_gzip_compress
from .connectivity_contract import (
    ConnectivityEdgePairs,
    build_connectivity_manifest,
    validate_connectivity_edges,
)
from .prepare_data import (
    _assert_unique_filename_components,
    _normalize_finite_float32_embedding,
    _require_continuous_obs_values,
    _require_finite_float32_array,
    _require_native_boolean,
    _require_nonempty_string,
    _require_positive_native_integer,
    _require_string_identifiers,
)
from .vector_fields import (
    is_vector_field_declaration_key,
    scale_vector_field,
    validate_vector_fields,
)

if TYPE_CHECKING:
    import anndata

logger = logging.getLogger("cellucid.anndata_adapter")


def _require_dimension(value: object) -> int:
    """Require one native 1D, 2D, or 3D dimension."""
    if type(value) is not int:
        raise TypeError("dim must be a native integer dimension.")
    if value not in (1, 2, 3):
        raise ValueError(f"Invalid dimension {value}. Expected 1, 2, or 3.")
    return value


def _classify_anndata_path(path: str | Path) -> Literal["h5ad", "zarr"]:
    """Classify one explicitly declared AnnData file or materialized Zarr store."""
    resolved_path = Path(path).resolve()
    if not resolved_path.exists():
        raise FileNotFoundError(f"Path not found: {resolved_path}")

    if resolved_path.is_file():
        if resolved_path.suffix != ".h5ad":
            raise ValueError(f"AnnData files must use the exact .h5ad extension: {resolved_path}")
        return "h5ad"

    if not resolved_path.is_dir():
        raise ValueError(
            f"AnnData path must be an .h5ad file or a materialized Zarr directory: {resolved_path}"
        )

    zarr_v2_group = resolved_path / ".zgroup"
    zarr_v2_attrs = resolved_path / ".zattrs"
    zarr_v3_metadata = resolved_path / "zarr.json"
    has_v2_marker = zarr_v2_group.exists() or zarr_v2_attrs.exists()

    if zarr_v3_metadata.exists():
        raise ValueError(
            f"AnnData Zarr input must use the exact Zarr v2 root contract "
            f"(.zgroup and .zattrs, without zarr.json): {resolved_path}"
        )

    if has_v2_marker:
        if not zarr_v2_group.is_file() or not zarr_v2_attrs.is_file():
            raise ValueError(
                f"Zarr v2 directory must contain both .zgroup and .zattrs files: {resolved_path}"
            )
        try:
            group_metadata = json.loads(zarr_v2_group.read_text(encoding="utf-8"))
            attributes = json.loads(zarr_v2_attrs.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            raise ValueError(
                f"Zarr v2 root metadata must be readable UTF-8 JSON: {resolved_path}"
            ) from error
        if group_metadata != {"zarr_format": 2}:
            raise ValueError(f"Zarr v2 .zgroup must declare exactly zarr_format 2: {resolved_path}")
        if not isinstance(attributes, dict):
            raise ValueError(f"Zarr v2 .zattrs must contain a JSON object: {resolved_path}")
        return "zarr"

    raise ValueError(
        f"AnnData path must be an .h5ad file or a complete Zarr v2 directory: {resolved_path}"
    )


_DIRECT_PAYLOAD_UNSAFE_RE = re.compile(r"[^A-Za-z0-9._-]+", flags=re.ASCII)
_DIRECT_PAYLOAD_EDGE_RE = re.compile(r"^[._]+|[._]+$", flags=re.ASCII)


def _direct_payload_component(name: object, *, label: str) -> str:
    """Encode one display identifier as the current direct-server route key."""
    if not isinstance(name, str):
        raise TypeError(f"{label} must be a native string.")
    if not name:
        raise ValueError(f"{label} must be non-empty.")
    component = _DIRECT_PAYLOAD_UNSAFE_RE.sub("_", name)
    component = _DIRECT_PAYLOAD_EDGE_RE.sub("", component)
    if not component:
        raise ValueError(
            f"{label} {name!r} has no ASCII letters, digits, '.', '_', or '-' "
            "for its direct-server payload route."
        )
    return component


def _assert_unique_direct_payload_components(
    keys: list[str],
    *,
    label: str,
) -> list[str]:
    """Require unique raw identifiers and collision-free direct route keys."""
    raw_keys: set[str] = set()
    component_to_raw: dict[str, str] = {}
    components: list[str] = []
    for key in keys:
        if key in raw_keys:
            raise ValueError(f"{label} key {key!r} is duplicated.")
        component = _direct_payload_component(
            key,
            label=f"{label} key",
        )
        previous = component_to_raw.get(component)
        if previous is not None:
            raise ValueError(
                f"{label} keys {previous!r} and {key!r} collide at direct "
                f"payload component {component!r}."
            )
        raw_keys.add(key)
        component_to_raw[component] = key
        components.append(component)
    return components


def _validate_manifest_field_keys(keys: list[str], *, label: str) -> None:
    _assert_unique_direct_payload_components(keys, label=label)


def _json_category_values(values: Any, *, field_key: str) -> list[str | bool | int | float]:
    categories: list[str | bool | int | float] = []
    seen: set[tuple[str, Any]] = set()
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
                    f"{integer_value!r} outside JavaScript's exact integer range"
                )
            value = integer_value
            token = ("number", value)
        elif isinstance(value, Real) and math.isfinite(value):
            value = float(value)
            token = ("number", value)
        else:
            raise ValueError(f"Categorical field {field_key!r} labels must be finite JSON scalars")
        if token in seen:
            raise ValueError(
                f"Categorical field {field_key!r} labels must be unique "
                "after exact JSON representation"
            )
        seen.add(token)
        categories.append(value)
    return categories


def _categorical_storage(
    n_categories: int,
    *,
    field_key: str,
) -> tuple[type[np.uint8] | type[np.uint16], str, int]:
    if n_categories > 65_535:
        raise ValueError(
            f"Categorical field {field_key!r} has {n_categories:,} categories; "
            "the current contract supports at most 65,535"
        )
    if n_categories <= 255:
        return np.uint8, "uint8", 255
    return np.uint16, "uint16", 65_535


class LRUCache:
    """
    Simple LRU cache with O(1) operations using OrderedDict.

    WARNING: This class is NOT thread-safe. The AnnDataAdapter is designed
    for single-threaded use within an HTTP server (each request is handled
    sequentially). If you need thread-safe caching, wrap operations in
    external locks or use a thread-safe cache implementation.

    The HTTP server uses a single adapter instance per server, and Python's
    GIL provides some level of safety for simple operations, but concurrent
    access from multiple threads may cause race conditions.
    """

    def __init__(self, max_size: int = 100):
        if max_size < 1:
            raise ValueError("max_size must be at least 1")
        self._cache: OrderedDict[str, Any] = OrderedDict()
        self._max_size = max_size

    def get(self, key: str) -> Any | None:
        """Get item and move to end (most recently used)."""
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        return None

    def put(self, key: str, value: Any) -> None:
        """Add item, evicting oldest if at capacity."""
        if key in self._cache:
            self._cache.move_to_end(key)
        else:
            if len(self._cache) >= self._max_size:
                self._cache.popitem(last=False)  # Remove oldest
        self._cache[key] = value

    def __contains__(self, key: str) -> bool:
        return key in self._cache

    def clear(self) -> None:
        """Clear all cached items and release memory."""
        self._cache.clear()

    def __len__(self) -> int:
        return len(self._cache)

    @property
    def max_size(self) -> int:
        """Maximum cache size."""
        return self._max_size


def _to_dense_1d(arr: np.ndarray | sparse.spmatrix) -> np.ndarray:
    """Convert sparse vector to dense numpy array."""
    if sparse.issparse(arr):
        return np.asarray(cast(sparse.spmatrix, arr).toarray()).flatten()
    arr = np.asarray(arr)
    if arr.ndim > 1:
        return arr.flatten()
    return arr


def _to_dense_2d(arr: np.ndarray | sparse.spmatrix) -> np.ndarray:
    """Convert sparse matrix to dense numpy array."""
    if sparse.issparse(arr):
        return np.asarray(cast(sparse.spmatrix, arr).toarray())
    return np.asarray(arr)


class AnnDataAdapter:
    """
    Adapter that wraps AnnData and provides data in Cellucid format.

    This adapter generates all the data that would normally be created by
    prepare, but reads directly from AnnData without creating
    intermediate files. This is slower but more convenient for interactive use.
    """

    def __init__(
        self,
        adata: anndata.AnnData,
        latent_key: str | None = None,
        gene_id_column: str | None = None,
        normalize_embeddings: bool = True,
        centroid_outlier_quantile: float = 0.95,
        centroid_min_points: int = 10,
        *,
        dataset_name: str,
        dataset_id: str,
        vector_field_default: str | None = None,
    ) -> None:
        """
        Initialize the adapter.

        Parameters
        ----------
        adata : AnnData
            AnnData object to adapt. Can be in-memory or backed (h5ad file).
        latent_key : str, optional
            Key in obsm for latent space used for outlier quantile calculation.
            If None, outlier quantiles are explicitly unavailable.
        gene_id_column : str, optional
            Exact column in ``var`` containing gene identifiers. If None,
            identifiers come from ``var.index``.
        normalize_embeddings : bool
            If True, normalize embeddings to [-1, 1] range (recommended).
        centroid_outlier_quantile : float
            Quantile for outlier removal in centroid computation.
        centroid_min_points : int
            Minimum points per category for centroid computation.
        dataset_name : str
            Explicit human-readable dataset name.
        dataset_id : str
            Explicit dataset identifier.
        vector_field_default : str, optional
            Exact field id to select when multiple UMAP vector fields exist.
        """
        self.adata: anndata.AnnData = adata
        if adata.isbacked and adata.file._filemode != "r":
            raise ValueError(
                "Backed AnnData inputs must be opened read-only. Pass the .h5ad path "
                "to Cellucid so it can own the exact read-only backing handle."
            )
        if adata.n_obs <= 0:
            raise ValueError("AnnData datasets must contain at least one observation")
        self.latent_key: str | None = latent_key
        if latent_key is not None:
            if not isinstance(latent_key, str) or not latent_key:
                raise ValueError("latent_key must be None or one non-empty string")
            if latent_key not in adata.obsm:
                raise ValueError(f"latent_key {latent_key!r} was not found in adata.obsm")
        self.gene_id_column: str | None
        if gene_id_column is None:
            self.gene_id_column = None
        else:
            self.gene_id_column = _require_nonempty_string(
                gene_id_column,
                label="gene_id_column",
            )
            if not self.gene_id_column.strip():
                raise ValueError("gene_id_column must be None or a non-blank string.")
        self.normalize_embeddings = _require_native_boolean(
            normalize_embeddings,
            label="normalize_embeddings",
        )

        if (
            type(centroid_outlier_quantile) is not float
            or not math.isfinite(centroid_outlier_quantile)
            or not (0.5 < centroid_outlier_quantile < 1.0)
        ):
            raise ValueError(
                "centroid_outlier_quantile must be a finite number strictly between 0.5 and 1.0"
            )
        self.centroid_outlier_quantile: float = centroid_outlier_quantile

        self.centroid_min_points = _require_positive_native_integer(
            centroid_min_points,
            label="centroid_min_points",
        )

        self.dataset_name: str = _require_nonempty_string(
            dataset_name,
            label="dataset_name",
        )
        self.dataset_id: str = _require_nonempty_string(
            dataset_id,
            label="dataset_id",
        )

        if self.gene_id_column is None:
            raw_gene_ids = self.adata.var.index.tolist()
            gene_label = "var index"
        else:
            if self.gene_id_column not in self.adata.var.columns:
                raise KeyError(
                    f"gene_id_column {self.gene_id_column!r} not found in var. "
                    f"Available columns: {list(self.adata.var.columns)}"
                )
            raw_gene_ids = self.adata.var[self.gene_id_column].tolist()
            gene_label = f"var column {self.gene_id_column!r}"
        gene_ids = _require_string_identifiers(
            raw_gene_ids,
            label=gene_label,
        )
        gene_payload_components = _assert_unique_direct_payload_components(
            gene_ids,
            label="Gene",
        )

        # Public source attribution never includes the local filesystem path.
        self._source_type: str = "memory"

        # Caches for computed data
        self._embedding_cache: dict[int, np.ndarray] = {}
        # Per-dimension normalization info for embeddings (center + scale_factor).
        # Used to scale vector fields into the same normalized space as points.
        self._embedding_norm: dict[int, dict[str, Any]] = {}

        # Optional per-cell vector fields (e.g. velocity, CellRank drift).
        # Structure matches dataset_identity.json: {"default_field": ..., "fields": {...}}.
        self._vector_fields_metadata: dict[str, Any] | None = None
        self._vector_field_obsm_keys: dict[str, dict[int, str]] = {}
        self._vector_field_cache: dict[tuple[str, int], np.ndarray] = {}
        self._centroid_cache: dict[str, list[dict[str, Any]]] = {}
        self._outlier_quantile_cache: dict[str, np.ndarray] = {}
        self._latent_space: np.ndarray | None = None
        self._gene_ids_cache: list[str] | None = gene_ids
        self._gene_id_to_idx_cache: dict[str, int] | None = {
            gene_id: index for index, gene_id in enumerate(gene_ids)
        }
        self._gene_payload_component_by_id: dict[str, str] = dict(
            zip(gene_ids, gene_payload_components, strict=True)
        )
        self._gene_id_by_payload_component: dict[str, str] = dict(
            zip(gene_payload_components, gene_ids, strict=True)
        )

        self._connectivity_cache: ConnectivityEdgePairs | None = None
        self._connectivity_manifest: dict[str, Any] | None = None
        self._has_connectivity: bool = False

        # CSC cache for efficient column access on CSR matrices
        # Note: This doubles memory for sparse matrices, but is necessary for O(1) gene access
        # Set to None initially; created lazily only when needed
        self._X_csc_cache: sparse.csc_matrix | None = None

        # LRU cache for gene expression values using O(1) OrderedDict
        self._gene_expression_cache: LRUCache = LRUCache(max_size=100)

        # UMAP embedding key resolution (dimension -> obsm key)
        self._umap_embedding_key_by_dim: dict[int, str] = {}
        # Track if adapter has been closed
        self._closed: bool = False
        self._close_complete: bool = False

        obs_keys = _require_string_identifiers(
            self.adata.obs.columns.tolist(),
            label="obs columns",
        )
        obs_payload_components = _assert_unique_direct_payload_components(
            obs_keys,
            label="Observation field",
        )
        self._obs_payload_component_by_key: dict[str, str] = dict(
            zip(obs_keys, obs_payload_components, strict=True)
        )
        self._obs_key_by_payload_component: dict[str, str] = dict(
            zip(obs_payload_components, obs_keys, strict=True)
        )
        for key in obs_keys:
            series = self.adata.obs[key]
            if (
                not isinstance(series.dtype, pd.CategoricalDtype)
                and not pd.api.types.is_bool_dtype(series)
                and pd.api.types.is_numeric_dtype(series)
            ):
                _require_continuous_obs_values(
                    series.to_numpy(),
                    key=key,
                    n_cells=int(self.adata.n_obs),
                )

        # Auto-detect available dimensions
        self._available_dimensions: list[int] = self._detect_dimensions()
        self._default_dimension: int = self._select_default_dimension()

        # Detect optional per-cell vector fields aligned to UMAP.
        self._vector_fields_metadata = self._detect_vector_fields(
            vector_field_default=vector_field_default,
        )

        # Connectivity is an advertised scientific capability, so validate and
        # materialize its exact edge contract before this adapter can be used.
        self._initialize_connectivity()

        logger.info(
            f"AnnDataAdapter initialized: {self.n_cells:,} cells, "
            f"{self.n_genes:,} genes, dimensions: {self._available_dimensions}"
        )

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        *,
        latent_key: str | None = None,
        gene_id_column: str | None = None,
        normalize_embeddings: bool = True,
        centroid_outlier_quantile: float = 0.95,
        centroid_min_points: int = 10,
        dataset_name: str,
        dataset_id: str,
        vector_field_default: str | None = None,
    ) -> AnnDataAdapter:
        """
        Create an adapter from an h5ad file or a materialized zarr store.

        Supports both:
        - .h5ad files: HDF5-based, supports true backed mode with memory-mapping
        - .zarr directories: loaded eagerly by ``anndata.read_zarr``

        Parameters
        ----------
        path : str or Path
            Path to h5ad file or zarr directory.
            - For h5ad: path/to/file.h5ad
            - For zarr: path/to/store.zarr (must be a directory)
        latent_key : str, optional
            Explicit key in ``obsm`` for the latent space.
        gene_id_column : str, optional
            Exact column in ``var`` containing gene identifiers. If None,
            identifiers come from ``var.index``.
        normalize_embeddings : bool
            Whether to normalize embeddings into the viewer coordinate range.
        centroid_outlier_quantile : float
            Quantile used when computing categorical centroids.
        centroid_min_points : int
            Minimum category size used for centroid computation.
        dataset_name : str
            Explicit human-readable dataset name.
        dataset_id : str
            Explicit stable dataset identifier.
        vector_field_default : str, optional
            Exact field id required when multiple UMAP vector fields exist.
        Returns
        -------
        AnnDataAdapter
            Adapter instance wrapping the loaded data.

        Raises
        ------
        FileNotFoundError
            If the path does not exist.
        ValueError
            If the path is not a valid h5ad or zarr store.

        Examples
        --------
        >>> # H5AD is always opened in read-only backed mode
        >>> adapter = AnnDataAdapter.from_file(
        ...     "data.h5ad",
        ...     dataset_name="Example",
        ...     dataset_id="example",
        ... )
        """
        import anndata as ad

        dataset_name = _require_nonempty_string(
            dataset_name,
            label="dataset_name",
        )
        dataset_id = _require_nonempty_string(
            dataset_id,
            label="dataset_id",
        )
        path = Path(path).resolve()

        source_type = _classify_anndata_path(path)
        if source_type == "zarr":
            adata = ad.read_zarr(path)
            logger.info("Loaded zarr store eagerly: %s", path)
        else:
            adata = ad.read_h5ad(path, backed="r")
            logger.info(
                "Loaded h5ad file: %s (is_backed=%s)",
                path,
                bool(adata.isbacked),
            )

        # Create adapter and store non-sensitive source info. A failed adapter
        # construction must not leave a backed file handle open.
        try:
            adapter = cls(
                adata,
                latent_key=latent_key,
                gene_id_column=gene_id_column,
                normalize_embeddings=normalize_embeddings,
                centroid_outlier_quantile=centroid_outlier_quantile,
                centroid_min_points=centroid_min_points,
                dataset_name=dataset_name,
                dataset_id=dataset_id,
                vector_field_default=vector_field_default,
            )
        except BaseException:
            if getattr(adata, "isbacked", False):
                try:
                    adata.file.close()
                except BaseException:
                    logger.exception(
                        "Failed to close the owned AnnData file after adapter construction failed"
                    )
            raise
        adapter._source_type = source_type
        return adapter

    def _check_closed(self) -> None:
        """Raise error if adapter has been closed."""
        if self._closed:
            raise RuntimeError("AnnDataAdapter has been closed. Create a new adapter instance.")

    @property
    def n_cells(self) -> int:
        """Number of cells."""
        self._check_closed()
        return cast(int, self.adata.n_obs)

    @property
    def n_genes(self) -> int:
        """Number of genes."""
        self._check_closed()
        return cast(int, self.adata.n_vars)

    @property
    def is_backed(self) -> bool:
        """Whether the live AnnData object reports an active backed store."""
        self._check_closed()
        return cast(bool, self.adata.isbacked)

    # =========================================================================
    # DIMENSION DETECTION
    # =========================================================================

    def _detect_dimensions(self) -> list[int]:
        """
        Validate and resolve the exact supported UMAP keys in ``obsm``.

        Returns
        -------
        list[int]
            Sorted list of available dimensions (e.g., [2, 3])

        Raises
        ------
        ValueError
            If a declared embedding is malformed, ambiguous, or absent.
        """
        self._umap_embedding_key_by_dim = {}
        available: list[int] = []

        for dim in [1, 2, 3]:
            key = f"X_umap_{dim}d"
            if key not in self.adata.obsm:
                continue
            embedding = _require_finite_float32_array(
                self.adata.obsm[key],
                label=f"Embedding {key!r}",
            )
            if embedding.ndim != 2:
                raise ValueError(
                    f"Embedding {key!r} must be a 2D array; got shape {embedding.shape}."
                )
            if embedding.shape[0] != self.n_cells:
                raise ValueError(
                    f"Embedding {key!r} must have exactly {self.n_cells} rows; "
                    f"got {embedding.shape[0]}."
                )
            if embedding.shape[1] != dim:
                raise ValueError(
                    f"Embedding {key!r} must have exactly {dim} columns; got {embedding.shape[1]}."
                )
            if self.normalize_embeddings:
                _normalize_finite_float32_embedding(
                    embedding,
                    label=f"Embedding {key!r}",
                )
            available.append(dim)
            self._umap_embedding_key_by_dim[dim] = key

        if not available:
            raise ValueError(
                "No supported UMAP embedding was declared in adata.obsm. "
                "Declare exactly one or more of 'X_umap_1d', 'X_umap_2d', "
                "or 'X_umap_3d'. "
                f"Available obsm keys: {list(self.adata.obsm.keys())}."
            )

        return sorted(available)

    def _select_default_dimension(self) -> int:
        """Select the highest declared dimension."""
        return max(self._available_dimensions)

    def _detect_vector_fields(
        self,
        *,
        vector_field_default: object,
    ) -> dict[str, Any] | None:
        """Validate UMAP vectors; each field defaults to its highest dimension."""
        candidates: dict[object, object] = {}
        for key in self.adata.obsm:
            if not isinstance(key, str):
                raise TypeError("AnnData obsm keys must be native strings.")
            if not key:
                raise ValueError("AnnData obsm keys must be non-empty strings.")
            if key in {
                "X_umap",
                "X_umap_1d",
                "X_umap_2d",
                "X_umap_3d",
            }:
                continue
            if is_vector_field_declaration_key(key):
                candidates[key] = self.adata.obsm[key]

        validated = validate_vector_fields(
            candidates,
            n_cells=self.n_cells,
            available_dimensions=self._available_dimensions,
            vector_field_default=vector_field_default,
        )
        self._vector_field_obsm_keys = validated.source_keys
        if not validated.fields:
            return None

        fields: dict[str, dict[str, Any]] = {}
        _assert_unique_filename_components(
            list(validated.fields),
            label="Vector field",
        )
        for field_id, vectors_by_dimension in validated.fields.items():
            dimensions = sorted(vectors_by_dimension)
            fields[field_id] = {
                "label": field_id,
                "basis": "umap",
                "available_dimensions": dimensions,
                "default_dimension": max(dimensions),
                "files": {
                    f"{dimension}d": (f"vectors/{field_id}_{dimension}d.bin")
                    for dimension in dimensions
                },
            }

        return {
            "default_field": validated.default_field,
            "fields": fields,
        }

    # =========================================================================
    # EMBEDDING DATA
    # =========================================================================

    def get_embedding(self, dim: int) -> np.ndarray:
        """
        Get embedding coordinates for a dimension.

        Returns normalized Float32 array of shape (n_cells, dim).
        """
        self._check_closed()
        dim = _require_dimension(dim)

        if dim not in self._available_dimensions:
            raise ValueError(
                f"Dimension {dim}D not available. Available: {self._available_dimensions}"
            )

        if dim in self._embedding_cache:
            return self._embedding_cache[dim]

        # Get raw embedding
        key = self._umap_embedding_key_by_dim.get(dim)
        if not key or key not in self.adata.obsm:
            raise ValueError(
                f"Could not find embedding for dimension {dim}. "
                f"Expected an embedding resolved for {dim}D. Available obsm keys: {list(self.adata.obsm.keys())}"
            )
        raw = _require_finite_float32_array(
            self.adata.obsm[key],
            label=f"Embedding {key!r}",
        )
        if raw.ndim != 2 or raw.shape != (self.n_cells, dim):
            raise ValueError(
                f"Embedding {key!r} changed after adapter initialization; "
                f"expected shape ({self.n_cells}, {dim}), got {raw.shape}."
            )

        if self.normalize_embeddings:
            raw, center, scale_factor, _max_range = _normalize_finite_float32_embedding(
                raw,
                label=f"Embedding {key!r}",
            )
            self._embedding_norm[dim] = {
                "center": center,
                "scale_factor": float(scale_factor),
            }
        else:
            self._embedding_norm[dim] = {
                "center": np.zeros((dim,), dtype=np.float32),
                "scale_factor": 1.0,
            }

        self._embedding_cache[dim] = raw
        return raw

    def get_embedding_3d(self, dim: int) -> np.ndarray:
        """
        Get embedding padded to 3D for WebGL rendering.

        1D -> (x, 0, 0)
        2D -> (x, y, 0)
        3D -> (x, y, z)
        """
        embedding = self.get_embedding(dim)
        n_cells = embedding.shape[0]

        if dim == 3:
            return embedding
        elif dim == 2:
            result = np.zeros((n_cells, 3), dtype=np.float32)
            result[:, :2] = embedding
            return result
        elif dim == 1:
            result = np.zeros((n_cells, 3), dtype=np.float32)
            result[:, 0] = embedding[:, 0]
            return result
        else:
            raise ValueError(f"Unsupported dimension: {dim}")

    def get_points_binary(self, dim: int, compress: bool = False) -> bytes:
        """Get embedding as binary data (for HTTP response)."""
        compress = _require_native_boolean(compress, label="compress")
        embedding = self.get_embedding(dim)
        data = embedding.astype(np.float32).tobytes()

        if compress:
            return deterministic_gzip_compress(data, compresslevel=6)
        return data

    def get_vector_field_binary(self, field_id: str, dim: int, compress: bool = False) -> bytes:
        """
        Get a per-cell vector field (displacement vectors) as binary float32 data.

        Vector fields are scaled by the SAME per-dimension normalization scale as
        the embedding points, so they are in the same normalized space as the
        `points_{dim}d.bin` responses.
        """
        self._check_closed()

        field = _require_nonempty_string(field_id, label="field_id")
        dimension = _require_dimension(dim)
        compress = _require_native_boolean(compress, label="compress")

        if self._vector_fields_metadata is None:
            raise ValueError(f"Vector field {field!r} not available.")
        fields = self._vector_fields_metadata["fields"]
        if field not in fields:
            raise ValueError(f"Vector field '{field}' not available")

        source_keys = self._vector_field_obsm_keys[field]
        if dimension not in source_keys:
            raise ValueError(f"Vector field {field!r} does not provide {dimension}D data.")
        obsm_key = source_keys[dimension]
        if obsm_key not in self.adata.obsm:
            raise ValueError(f"Vector field source {obsm_key!r} no longer exists.")

        cache_key = (field, dimension)
        if cache_key not in self._vector_field_cache:
            raw = _require_finite_float32_array(
                self.adata.obsm[obsm_key],
                label=f"Vector field {obsm_key!r}",
            )
            expected_shape = (self.n_cells, dimension)
            if raw.ndim != 2 or raw.shape != expected_shape:
                raise ValueError(
                    f"Vector field {obsm_key!r} changed after adapter "
                    f"initialization; expected shape {expected_shape}, "
                    f"got {raw.shape}."
                )

            self.get_embedding(dimension)
            if dimension not in self._embedding_norm:
                raise RuntimeError(f"Missing normalization state for {dimension}D embedding.")
            scale_factor = self._embedding_norm[dimension]["scale_factor"]
            self._vector_field_cache[cache_key] = scale_vector_field(
                raw,
                scale_factor=scale_factor,
                label=f"Vector field {obsm_key!r}",
            )

        data = self._vector_field_cache[cache_key].tobytes()
        if compress:
            return deterministic_gzip_compress(data, compresslevel=6)
        return data

    # =========================================================================
    # LATENT SPACE FOR OUTLIER COMPUTATION
    # =========================================================================

    def _get_latent_space(self) -> np.ndarray | None:
        """Get latent space for outlier quantile computation."""
        if self._latent_space is not None:
            return self._latent_space

        if self.latent_key is None or self.latent_key not in self.adata.obsm:
            return None

        latent = _require_finite_float32_array(
            self.adata.obsm[self.latent_key],
            label=f"Latent field {self.latent_key!r}",
        )
        if latent.ndim != 2 or latent.shape[0] != self.n_cells:
            raise ValueError(
                f"Latent field {self.latent_key!r} must have shape "
                f"({self.n_cells}, n_dimensions), got {latent.shape}."
            )
        self._latent_space = latent
        return self._latent_space

    # =========================================================================
    # OBS DATA (CELL METADATA)
    # =========================================================================

    def get_obs_keys(self) -> list[str]:
        """Get list of obs column names."""
        self._check_closed()
        return list(self.adata.obs.columns)

    def get_obs_payload_component(self, key: str) -> str:
        """Return the exact direct-server route component for one obs key."""
        self._check_closed()
        try:
            return self._obs_payload_component_by_key[key]
        except KeyError as error:
            raise KeyError(f"Observation field {key!r} has no direct payload route.") from error

    def get_obs_key_for_payload_component(self, component: str) -> str | None:
        """Resolve one direct-server route component to its exact obs key."""
        self._check_closed()
        return self._obs_key_by_payload_component.get(component)

    def get_obs_field_kind(self, key: str) -> Literal["continuous", "category"]:
        """
        Determine if an obs field is continuous or categorical.

        Classification rules:
        - Categorical dtype → category
        - Boolean dtype → category
        - Numeric dtype → continuous
        - String/object → category (treated as labels)
        - Empty column → category (safe default)
        """
        self._check_closed()

        if key not in self.adata.obs.columns:
            raise KeyError(
                f"Field '{key}' not found in obs. Available: {list(self.adata.obs.columns)}"
            )

        s = self.adata.obs[key]

        if isinstance(s.dtype, pd.CategoricalDtype) or pd.api.types.is_bool_dtype(s):
            return "category"
        elif pd.api.types.is_numeric_dtype(s):
            return "continuous"
        elif pd.api.types.is_string_dtype(s) or pd.api.types.is_object_dtype(s):
            return "category"
        raise TypeError(f"obs field {key!r} has unsupported dtype {s.dtype!r}.")

    def get_obs_continuous_values(self, key: str, compress: bool = False) -> bytes:
        """
        Get continuous obs field as binary float32 data.

        Values must remain finite after float32 representation.

        Raises
        ------
        KeyError
            If the field is not found in obs.
        """
        self._check_closed()
        compress = _require_native_boolean(compress, label="compress")

        if key not in self.adata.obs.columns:
            raise KeyError(
                f"obs field '{key}' not found. Available fields: {list(self.adata.obs.columns)}"
            )

        values = _require_continuous_obs_values(
            self.adata.obs[key].to_numpy(),
            key=key,
            n_cells=self.n_cells,
        )

        data = cast(bytes, values.tobytes())

        if compress:
            return deterministic_gzip_compress(data, compresslevel=6)
        return data

    def get_obs_categorical_codes(
        self,
        key: str,
        compress: bool = False,
    ) -> tuple[bytes, list[str | bool | int | float], int]:
        """
        Get categorical obs field as binary codes.

        Returns:
            (binary_codes, category_list, missing_value)

        Categories are assigned codes 0 to n-1.
        Missing values (NaN) are encoded as the missing_value sentinel.

        Raises
        ------
        KeyError
            If the field is not found in obs.
        """
        self._check_closed()
        compress = _require_native_boolean(compress, label="compress")

        if key not in self.adata.obs.columns:
            raise KeyError(
                f"obs field '{key}' not found. Available fields: {list(self.adata.obs.columns)}"
            )

        s = self.adata.obs[key]
        cat = s.astype("category")
        categories = _json_category_values(
            cat.cat.categories,
            field_key=key,
        )
        codes = cat.cat.codes.to_numpy(dtype=np.int32)

        n_categories = len(categories)
        dtype, _dtype_str, missing_value = _categorical_storage(
            n_categories,
            field_key=key,
        )

        # Convert codes (-1 for NaN -> missing_value)
        codes_typed = np.full(self.n_cells, missing_value, dtype=dtype)
        valid_mask = codes >= 0
        codes_typed[valid_mask] = codes[valid_mask].astype(dtype)

        # Log missing value count
        n_missing = (~valid_mask).sum()
        if n_missing > 0:
            logger.debug(f"obs field '{key}': {n_missing} missing values out of {self.n_cells}")

        data = codes_typed.tobytes()
        if compress:
            data = deterministic_gzip_compress(data, compresslevel=6)

        return data, categories, int(missing_value)

    # =========================================================================
    # CENTROIDS AND OUTLIER QUANTILES
    # =========================================================================

    def _compute_centroids_for_field(
        self,
        key: str,
        dim: int,
    ) -> list[dict]:
        """Compute centroids for a categorical field at a given dimension."""
        cache_key = f"{key}_{dim}d"
        if cache_key in self._centroid_cache:
            return self._centroid_cache[cache_key]

        coords = self.get_embedding(dim)
        s = self.adata.obs[key]
        cat = s.astype("category")
        categories = _json_category_values(
            cat.cat.categories,
            field_key=key,
        )
        codes = cat.cat.codes.to_numpy(dtype=np.int32)

        centroids: list[dict[str, Any]] = []
        for code, label in enumerate(categories):
            mask = codes == code
            idx = np.nonzero(mask)[0]
            n = idx.size

            if n < self.centroid_min_points:
                continue

            pts = coords[idx, :]
            center = pts.mean(axis=0)

            # Remove outliers for centroid computation
            if n > self.centroid_min_points:
                dists = np.linalg.norm(pts - center, axis=1)
                thr = float(np.quantile(dists, self.centroid_outlier_quantile))
                inlier_mask = dists <= thr
                n_in = int(inlier_mask.sum())
                if n_in >= self.centroid_min_points:
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

        self._centroid_cache[cache_key] = centroids
        return centroids

    def get_centroids_for_field(self, key: str) -> dict[str, list[dict]]:
        """Get centroids for all available dimensions."""
        result = {}
        for dim in self._available_dimensions:
            centroids = self._compute_centroids_for_field(key, dim)
            result[str(dim)] = centroids
        return result

    def _compute_outlier_quantiles(self, key: str) -> np.ndarray:
        """Compute per-cell outlier quantiles for a categorical field."""
        if key in self._outlier_quantile_cache:
            return self._outlier_quantile_cache[key]

        latent = self._get_latent_space()
        if latent is None:
            raise ValueError(
                f"Cannot compute outlier quantiles for {key!r}: "
                "the dataset has no declared latent space"
            )

        s = self.adata.obs[key]
        cat = s.astype("category")
        categories = _json_category_values(
            cat.cat.categories,
            field_key=key,
        )
        codes = cat.cat.codes.to_numpy(dtype=np.int32)

        quantiles = np.full(self.n_cells, np.nan, dtype=np.float32)

        for code, _label in enumerate(categories):
            mask = codes == code
            idx = np.nonzero(mask)[0]
            n = idx.size

            if n < self.centroid_min_points:
                continue

            pts = latent[idx, :]
            centroid = pts.mean(axis=0)
            dists = np.linalg.norm(pts - centroid, axis=1)
            sorted_dists = np.sort(dists)
            ranks = np.searchsorted(sorted_dists, dists, side="right")
            cell_quantiles = ranks.astype(np.float32) / n
            quantiles[idx] = cell_quantiles

        self._outlier_quantile_cache[key] = quantiles
        return quantiles

    def get_obs_outlier_quantiles(self, key: str, compress: bool = False) -> bytes:
        """Get outlier quantiles as binary float32 data."""
        self._check_closed()
        compress = _require_native_boolean(compress, label="compress")
        quantiles = self._compute_outlier_quantiles(key)
        data = quantiles.tobytes()

        if compress:
            return deterministic_gzip_compress(data, compresslevel=6)
        return data

    # =========================================================================
    # VAR DATA (GENE EXPRESSION)
    # =========================================================================

    def get_gene_ids(self) -> list[str]:
        """Get list of gene identifiers."""
        self._check_closed()

        if self._gene_ids_cache is not None:
            return self._gene_ids_cache

        if self.gene_id_column is None:
            values = self.adata.var.index.tolist()
            label = "var index"
        else:
            if self.gene_id_column not in self.adata.var.columns:
                raise KeyError(
                    f"gene_id_column '{self.gene_id_column}' not found in var. "
                    f"Available columns: {list(self.adata.var.columns)}"
                )
            values = self.adata.var[self.gene_id_column].tolist()
            label = f"var column {self.gene_id_column!r}"

        self._gene_ids_cache = _require_string_identifiers(
            values,
            label=label,
        )
        payload_components = _assert_unique_direct_payload_components(
            self._gene_ids_cache,
            label="Gene",
        )

        # Build index for O(1) lookup
        self._gene_id_to_idx_cache = {gid: idx for idx, gid in enumerate(self._gene_ids_cache)}
        self._gene_payload_component_by_id = dict(
            zip(self._gene_ids_cache, payload_components, strict=True)
        )
        self._gene_id_by_payload_component = dict(
            zip(payload_components, self._gene_ids_cache, strict=True)
        )
        return self._gene_ids_cache

    def get_gene_payload_component(self, gene_id: str) -> str:
        """Return the exact direct-server route component for one gene."""
        self._check_closed()
        try:
            return self._gene_payload_component_by_id[gene_id]
        except KeyError as error:
            raise KeyError(f"Gene {gene_id!r} has no direct payload route.") from error

    def get_gene_id_for_payload_component(self, component: str) -> str | None:
        """Resolve one direct-server route component to its exact gene identifier."""
        self._check_closed()
        return self._gene_id_by_payload_component.get(component)

    def _get_gene_idx(self, gene_id: str) -> int:
        """Get gene index with O(1) lookup."""
        # Ensure cache is populated
        if self._gene_id_to_idx_cache is None:
            self.get_gene_ids()
        gene_id_to_idx = self._gene_id_to_idx_cache
        if gene_id_to_idx is None:
            raise RuntimeError("Gene identifier index was not initialized")

        if gene_id not in gene_id_to_idx:
            raise KeyError(f"Gene '{gene_id}' not found in var")

        return gene_id_to_idx[gene_id]

    def _get_gene_column(self, gene_idx: int) -> np.ndarray:
        """
        Extract a single gene column from the expression matrix.

        Handles multiple data formats efficiently:
        - Dense numpy arrays: Direct column access O(1)
        - CSC sparse matrices: Efficient column access O(nnz/n_cols)
        - CSR sparse matrices: Convert to CSC and cache for repeated access
        - Other sparse formats (COO, DOK, etc.): Convert to CSC and cache
        - None: Rejected because no expression payload exists
        - Backed/chunked arrays: Direct slicing (no caching to preserve lazy loading)

        Memory-Speed Tradeoff:
        For in-memory CSR matrices, the CSC cache doubles memory usage but enables
        O(nnz/n_cols) column access instead of O(nnz). This pays off after
        accessing more than a few genes.

        IMPORTANT: For backed h5ad files, we do NOT create a CSC cache because:
        1. It would force loading the entire X matrix into memory
        2. It defeats the purpose of lazy loading for large datasets
        Instead, we accept slower column access to preserve memory efficiency.

        Parameters
        ----------
        gene_idx : int
            Column index of the gene to extract.

        Returns
        -------
        np.ndarray
            1D array of expression values for all cells.
        """
        self._check_closed()

        X = self.adata.X

        if X is None:
            raise ValueError("AnnData has no X expression matrix; gene expression is unavailable")

        # Edge case: invalid gene index
        if gene_idx < 0 or gene_idx >= self.n_genes:
            raise IndexError(f"Gene index {gene_idx} out of range [0, {self.n_genes})")

        # CRITICAL: For backed h5ad files, avoid creating CSC cache to preserve
        # lazy loading behavior. The memory cost of CSC cache would defeat the
        # purpose of backed mode for large datasets.
        if self._source_type == "h5ad" and self.is_backed:
            # For backed h5ad, use direct slicing (slower but memory-efficient)
            # Note: backed h5ad matrices support slicing but may be slower
            col = X[:, gene_idx]
            if sparse.issparse(col):
                return np.asarray(col.toarray()).flatten()
            return np.asarray(col).flatten()

        if sparse.issparse(X):
            # For sparse matrix, CSC is efficient for column access
            if sparse.isspmatrix_csc(X):
                return np.asarray(X.getcol(gene_idx).toarray()).flatten()
            elif sparse.isspmatrix_csr(X):
                # CSR is inefficient for column access - convert to CSC and cache
                # This is a one-time cost that pays off with repeated gene queries
                # Only do this for in-memory data (not backed files)
                if self._X_csc_cache is None:
                    logger.info(
                        f"Converting CSR matrix ({X.shape[0]:,}×{X.shape[1]:,}, "
                        f"{X.nnz:,} non-zeros) to CSC for efficient column access"
                    )
                    self._X_csc_cache = X.tocsc()
                return np.asarray(self._X_csc_cache.getcol(gene_idx).toarray()).flatten()
            else:
                # Other sparse format (COO, LIL, DOK, BSR) - convert to CSC
                if self._X_csc_cache is None:
                    logger.info(f"Converting {type(X).__name__} to CSC for column access")
                    self._X_csc_cache = sparse.csc_matrix(X)
                return np.asarray(self._X_csc_cache.getcol(gene_idx).toarray()).flatten()
        else:
            # Dense matrix - direct column access
            return np.asarray(X[:, gene_idx]).flatten()

    def _get_validated_gene_values(self, gene_id: str) -> np.ndarray:
        """Return one exact finite Float32 per-cell gene vector."""
        if self.adata.X is None:
            raise ValueError("AnnData has no X expression matrix; gene expression is unavailable")
        values = self._gene_expression_cache.get(gene_id)
        if values is None:
            gene_idx = self._get_gene_idx(gene_id)
            col = self._get_gene_column(gene_idx)
            values = _require_finite_float32_array(
                col,
                label=f"Gene {gene_id!r} expression",
            )
            if values.ndim != 1 or values.shape[0] != self.n_cells:
                raise ValueError(
                    f"Gene {gene_id!r} expression must have shape ({self.n_cells},), "
                    f"got {values.shape}."
                )
            self._gene_expression_cache.put(gene_id, values)
        return values

    def get_gene_expression(self, gene_id: str, compress: bool = False) -> bytes:
        """Get expression values for a single gene as binary float32."""
        self._check_closed()
        compress = _require_native_boolean(compress, label="compress")
        values = self._get_validated_gene_values(gene_id)
        data = values.tobytes()
        if compress:
            return deterministic_gzip_compress(data, compresslevel=6)
        return data

    def get_gene_min_max(self, gene_id: str) -> tuple[float, float]:
        """Get min/max values for a gene (for colormap scaling)."""
        values = self._get_validated_gene_values(gene_id)
        return float(values.min()), float(values.max())

    # =========================================================================
    # CONNECTIVITY DATA
    # =========================================================================

    def _initialize_connectivity(self) -> None:
        """Validate the complete optional connectivity capability atomically."""
        try:
            self._has_connectivity = "connectivities" in self.adata.obsp
        except Exception as exc:
            raise ValueError("Could not inspect adata.obsp for 'connectivities'") from exc

        if not self._has_connectivity:
            self._connectivity_manifest = None
            return

        try:
            edges = self._compute_connectivity_edges()
        except Exception as exc:
            raise ValueError(
                "Invalid adata.obsp['connectivities']; the direct server was not "
                f"initialized: {exc}"
            ) from exc

        self._connectivity_manifest = build_connectivity_manifest(
            n_cells=self.n_cells,
            n_edges=edges.n_edges,
            max_neighbors=edges.max_neighbors,
            index_bytes=edges.index_bytes,
            index_dtype=edges.index_dtype,
            compression=None,
        )

    def has_connectivity(self) -> bool:
        """Check if connectivity data is available."""
        self._check_closed()
        return self._has_connectivity

    def _compute_connectivity_edges(self) -> ConnectivityEdgePairs:
        """Validate and cache one aligned weighted edge payload."""
        self._check_closed()

        if self._connectivity_cache is not None:
            return self._connectivity_cache

        if not self.has_connectivity():
            raise ValueError("No connectivity data in adata.obsp['connectivities']")

        self._connectivity_cache = validate_connectivity_edges(
            self.adata.obsp["connectivities"],
            n_cells=self.n_cells,
        )
        return self._connectivity_cache

    def get_connectivity_edges(
        self,
        compress: bool = False,
    ) -> tuple[bytes, bytes, bytes, int, int]:
        """
        Get aligned endpoint and Float64 weight payloads.

        Returns:
            (sources_binary, destinations_binary, weights_binary,
            n_edges, max_neighbors)
        """
        self._check_closed()
        compress = _require_native_boolean(compress, label="compress")
        edges = self._compute_connectivity_edges()

        sources_data = edges.sources.tobytes()
        destinations_data = edges.destinations.tobytes()
        weights_data = edges.weights.tobytes()

        if compress:
            sources_data = deterministic_gzip_compress(
                sources_data,
                compresslevel=6,
            )
            destinations_data = deterministic_gzip_compress(
                destinations_data,
                compresslevel=6,
            )
            weights_data = deterministic_gzip_compress(
                weights_data,
                compresslevel=6,
            )

        return (
            sources_data,
            destinations_data,
            weights_data,
            edges.n_edges,
            edges.max_neighbors,
        )

    # =========================================================================
    # MANIFEST GENERATION
    # =========================================================================

    def get_dataset_identity(self) -> dict:
        """Generate dataset_identity.json content."""
        obs_manifest = self.get_obs_manifest()
        continuous_fields = obs_manifest["_continuousFields"]
        categorical_fields = obs_manifest["_categoricalFields"]
        obs_fields = [{"key": field[0], "kind": "continuous"} for field in continuous_fields] + [
            {
                "key": field[0],
                "kind": "category",
                "n_categories": len(field[1]),
            }
            for field in categorical_fields
        ]
        n_continuous = len(continuous_fields)
        n_categorical = len(categorical_fields)

        has_conn = self.has_connectivity()
        n_edges = (
            self._connectivity_manifest["n_edges"]
            if self._connectivity_manifest is not None
            else None
        )

        # Build embeddings metadata
        embeddings_meta = {
            "available_dimensions": self._available_dimensions,
            "default_dimension": self._default_dimension,
            "files": {f"{dim}d": f"points_{dim}d.bin" for dim in self._available_dimensions},
        }
        identity = {
            "version": 2,
            "id": self.dataset_id,
            "name": self.dataset_name,
            "description": "Loaded directly from AnnData",
            "created_at": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "cellucid_data_version": "anndata_adapter",
            "stats": {
                "n_cells": self.n_cells,
                "n_genes": 0 if self.adata.X is None else self.n_genes,
                "n_obs_fields": len(obs_fields),
                "n_categorical_fields": n_categorical,
                "n_continuous_fields": n_continuous,
                "has_connectivity": has_conn,
                "n_edges": n_edges,
            },
            "embeddings": embeddings_meta,
            "obs_fields": obs_fields,
            "export_settings": {
                # No file compression (data is served dynamically)
                # Network compression (gzip) is applied transparently by the server
                "compression": None,
                # No quantization - full float32 precision for all values
                # This ensures no data loss but increases network transfer size
                "var_quantization": None,
                "obs_continuous_quantization": None,
                # Each categorical field declares its exact uint8/uint16 storage
                # in obs_manifest.json according to its category count.
                "obs_categorical_dtype": "auto",
            },
            "source": {
                "name": {
                    "h5ad": "H5AD file",
                    "zarr": "Zarr store",
                    "memory": "In-memory AnnData",
                }[self._source_type],
            },
        }

        if self._vector_fields_metadata:
            identity["vector_fields"] = self._vector_fields_metadata

        return identity

    def get_obs_manifest(self) -> dict:
        """Generate obs_manifest.json content."""
        obs_keys = self.get_obs_keys()
        _validate_manifest_field_keys(obs_keys, label="Observation field")

        continuous_fields = []
        categorical_fields = []

        for key in obs_keys:
            kind = self.get_obs_field_kind(key)

            if kind == "continuous":
                continuous_fields.append([key])
            else:
                # Categorical
                s = self.adata.obs[key]
                cat = s.astype("category")
                categories = _json_category_values(
                    cat.cat.categories,
                    field_key=key,
                )
                n_categories = len(categories)

                _dtype, dtype_str, missing_value = _categorical_storage(
                    n_categories,
                    field_key=key,
                )

                # Get centroids for all dimensions
                centroids_by_dim = self.get_centroids_for_field(key)

                # Outlier quantile min/max (always 0-1 since it's a quantile)
                categorical_fields.append(
                    [key, categories, dtype_str, missing_value, centroids_by_dim]
                )

        # Build compact manifest format
        # Path patterns must match prepare format for consistency
        # Note: We don't use compression (.gz) in adapter mode - gzip is applied
        # at the HTTP level when client supports Accept-Encoding: gzip
        obs_schemas = {}
        if continuous_fields:
            obs_schemas["continuous"] = {
                "pathPattern": "obs/{key}.values.f32",
                "ext": "f32",
                "dtype": "float32",
                "quantized": False,
            }
        if categorical_fields:
            if self.latent_key is None:
                obs_schemas["categorical"] = {
                    "codesPathPattern": "obs/{key}.codes.{ext}",
                    "outlierPathPattern": None,
                    "outlierExt": None,
                    "outlierDtype": None,
                    "outlierQuantized": False,
                }
            else:
                obs_schemas["categorical"] = {
                    "codesPathPattern": "obs/{key}.codes.{ext}",
                    "outlierPathPattern": "obs/{key}.outliers.f32",
                    "outlierExt": "f32",
                    "outlierDtype": "float32",
                    "outlierQuantized": False,
                }

        return {
            "_format": "compact_v1",
            "n_points": self.n_cells,
            "centroid_outlier_quantile": self.centroid_outlier_quantile,
            "latent_key": self.latent_key,
            "compression": None,
            "_obsSchemas": obs_schemas,
            "_continuousFields": continuous_fields,
            "_categoricalFields": categorical_fields,
        }

    def get_var_manifest(self) -> dict:
        """Generate var_manifest.json content."""
        if self.adata.X is None:
            fields = []
        else:
            gene_ids = self.get_gene_ids()
            _validate_manifest_field_keys(gene_ids, label="Gene")
            fields = [[gid] for gid in gene_ids]

        # Path pattern must match prepare format for consistency
        # Note: We don't use compression (.gz) in adapter mode - gzip is applied
        # at the HTTP level when client supports Accept-Encoding: gzip
        var_schema = {
            "kind": "continuous",
            "pathPattern": "var/{key}.values.f32",
            "ext": "f32",
            "dtype": "float32",
            "quantized": False,
        }

        return {
            "_format": "compact_v1",
            "n_points": self.n_cells,
            "var_gene_id_column": self.gene_id_column,
            "compression": None,
            "quantization": None,
            "_varSchema": var_schema,
            "fields": fields,
        }

    def get_connectivity_manifest(self) -> dict | None:
        """Generate connectivity_manifest.json content."""
        self._check_closed()
        if self._connectivity_manifest is None:
            return None
        return dict(self._connectivity_manifest)

    # =========================================================================
    # CLEANUP AND CONTEXT MANAGER
    # =========================================================================

    def close(self) -> None:
        """
        Close the adapter and release all resources.

        This method:
        1. Clears all caches to free memory (embedding, centroid, CSC, gene expression)
        2. Closes the underlying file handle for backed h5ad files
        3. Marks the adapter as closed to prevent further operations

        Safe to call multiple times. Always call this method when done with the
        adapter, or use the context manager::

            with AnnDataAdapter.from_file(
                "data.h5ad",
                dataset_name="Example",
                dataset_id="example",
            ) as adapter:
                # use adapter
            # automatically cleaned up

        Memory Released
        ---------------

        - Embedding cache (normalized UMAP coordinates)
        - Centroid cache (computed label centroids)
        - Outlier quantile cache
        - Gene expression LRU cache (up to 100 gene columns)
        - CSC matrix cache (for CSR->CSC converted matrices)
        - Latent space array
        - Gene ID lookup indices
        """
        if self._close_complete:
            return

        # Capture the backed file before marking the adapter closed.
        file_to_close = self.adata.file if self.adata.isbacked else None

        self._closed = True

        # Clear all caches to free memory
        self._embedding_cache.clear()
        self._embedding_norm.clear()
        self._vector_field_cache.clear()
        self._vector_field_obsm_keys.clear()
        self._vector_fields_metadata = None
        self._centroid_cache.clear()
        self._outlier_quantile_cache.clear()
        self._gene_expression_cache.clear()
        self._latent_space = None
        self._connectivity_cache = None
        self._connectivity_manifest = None
        self._has_connectivity = False
        self._gene_ids_cache = None
        self._gene_id_to_idx_cache = None
        self._gene_payload_component_by_id.clear()
        self._gene_id_by_payload_component.clear()
        self._obs_payload_component_by_key.clear()
        self._obs_key_by_payload_component.clear()

        # Clear CSC cache (can be large for sparse matrices).
        self._X_csc_cache = None

        # Close backed file handle if applicable
        if file_to_close is not None:
            file_to_close.close()
            logger.debug(f"Closed backed file handle for {self.dataset_name}")

        # Clear reference to adata to help garbage collection
        # (do this last since we needed it for is_backed check)
        self.adata = None
        self._close_complete = True

        logger.debug(f"AnnDataAdapter closed: {self.dataset_name}")

    def __enter__(self) -> AnnDataAdapter:
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit - ensures cleanup."""
        self.close()

    def __repr__(self) -> str:
        if self._closed:
            return "AnnDataAdapter(closed)"
        backed_str = " (backed)" if self.is_backed else ""
        source_str = f" from {self._source_type}" if self._source_type != "memory" else ""
        return (
            f"AnnDataAdapter({self.n_cells:,} cells, {self.n_genes:,} genes, "
            f"dims={self._available_dimensions}{backed_str}{source_str})"
        )
