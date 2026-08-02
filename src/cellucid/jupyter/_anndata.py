"""
The notebook viewer that reads an AnnData object with no export step.

:class:`AnnDataViewer` starts a :class:`cellucid.anndata_server.AnnDataServer`,
which serves the adapter's virtual export instead of files on disk;
:func:`show_anndata` is the one-call form. Everything else -- hooks, commands,
session capture -- is inherited unchanged from :class:`BaseViewer`.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

from .._contracts import _require_dataset_id, _require_dataset_name
from .._server_base import CELLUCID_WEB_URL
from ._base import BaseViewer, logger

if TYPE_CHECKING:
    import anndata

class AnnDataViewer(BaseViewer):
    """
    Interactive viewer for AnnData objects in Jupyter notebooks.

    This viewer serves AnnData directly without requiring prepare.
    It's more convenient for interactive exploration but slower than using
    pre-exported data. Supports bidirectional communication via hooks.

    Supports:
    - In-memory AnnData objects
    - h5ad files (HDF5-based, with lazy loading via backed mode)
    - zarr stores materialized by ``anndata.read_zarr``

    Example:
        >>> viewer = AnnDataViewer(
        ...     adata,
        ...     dataset_name="Example",
        ...     dataset_id="example",
        ... )
        >>> viewer.display()
        >>>
        >>> @viewer.on_selection
        ... def analyze_selection(event):
        ...     subset = adata[event['cells']]
        ...     sc.pl.violin(subset, ['gene1', 'gene2'])
        >>>
        >>> # From h5ad file with lazy loading
        >>> viewer = AnnDataViewer(
        ...     "/path/to/data.h5ad",
        ...     dataset_name="Example",
        ...     dataset_id="example",
        ... )
    """

    def __init__(
        self,
        data: str | Path | anndata.AnnData,
        port: int | None = None,
        height: int = 600,
        auto_open: bool = True,
        *,
        latent_key: str | None = None,
        gene_id_column: str | None = None,
        normalize_embeddings: bool = True,
        centroid_outlier_quantile: float = 0.95,
        centroid_min_points: int = 10,
        dataset_name: str,
        dataset_id: str,
        obs_keys: Sequence[str] | None = None,
        vector_field_default: str | None = None,
        client_server_url: str | None = None,
        web_source_url: str = CELLUCID_WEB_URL,
        web_cache_dir: str | Path | None = None,
        allowed_hosts: Sequence[str] | None = None,
    ) -> None:
        """
        Initialize the AnnData viewer.

        Args:
            data: AnnData object or path to h5ad file or zarr directory. A
                leading ``~`` in a path is expanded to the user's home directory.
            port: Port for the data server (auto-selected if None).
            height: Height of the embedded viewer in pixels.
            auto_open: Automatically display when created.
            latent_key: Explicit key in ``obsm`` for the latent space.
            gene_id_column: Exact column in ``var`` containing gene identifiers.
                If None, identifiers come from ``var.index``.
            normalize_embeddings: Whether to normalize embeddings into the viewer
                coordinate range.
            centroid_outlier_quantile: Quantile used for categorical centroids.
            centroid_min_points: Minimum category size used for centroid computation.
            dataset_name: Exact non-empty human-readable name without
                surrounding whitespace or control characters.
            dataset_id: Portable 1-180 byte ASCII identifier beginning with a
                letter or digit and otherwise using letters, digits, ``.``,
                ``_``, or ``-``; no trailing ``.`` or Windows device name.
            obs_keys: Exact ``obs`` columns to serve, in this order. If None,
                every column is served. Naming the columns leaves out a column
                Cellucid cannot serve, such as a ``datetime64`` collection date.
            vector_field_default: Exact field id required when multiple UMAP
                vector fields exist.
            client_server_url: Exact browser-reachable data server base URL.
            web_source_url: Origin publishing the web asset inventory.
            web_cache_dir: Directory holding the active verified web build.
            allowed_hosts: Extra ``Host`` names the data server answers to.
                Required when ``client_server_url`` points at a reverse proxy
                such as ``jupyter-server-proxy``: name that proxy's host here.
        """
        dataset_name = _require_dataset_name(
            dataset_name,
            label="dataset_name",
        )
        dataset_id = _require_dataset_id(
            dataset_id,
            label="dataset_id",
        )
        super().__init__(
            port=port,
            height=height,
            client_server_url=client_server_url,
            web_source_url=web_source_url,
            web_cache_dir=web_cache_dir,
            allowed_hosts=allowed_hosts,
        )

        try:
            self.data = data

            self._start_server(
                latent_key=latent_key,
                gene_id_column=gene_id_column,
                normalize_embeddings=normalize_embeddings,
                centroid_outlier_quantile=centroid_outlier_quantile,
                centroid_min_points=centroid_min_points,
                dataset_name=dataset_name,
                dataset_id=dataset_id,
                obs_keys=obs_keys,
                vector_field_default=vector_field_default,
            )
            self._activate()

            if auto_open and self._context["in_jupyter"]:
                self.display()
        except BaseException:
            self._rollback_failed_construction()
            raise

    def _start_server(
        self,
        *,
        latent_key: str | None,
        gene_id_column: str | None,
        normalize_embeddings: bool,
        centroid_outlier_quantile: float,
        centroid_min_points: int,
        dataset_name: str,
        dataset_id: str,
        obs_keys: Sequence[str] | None,
        vector_field_default: str | None,
    ) -> None:
        """Start the AnnData server."""
        from ..anndata_server import AnnDataServer

        self._server = AnnDataServer(
            data=self.data,
            port=self.port,
            host="127.0.0.1",
            open_browser=False,
            quiet=True,
            serve_web_ui=True,
            web_source_url=self.web_source_url,
            web_cache_dir=self.web_cache_dir,
            allowed_hosts=sorted(self.allowed_hosts),
            latent_key=latent_key,
            gene_id_column=gene_id_column,
            normalize_embeddings=normalize_embeddings,
            centroid_outlier_quantile=centroid_outlier_quantile,
            centroid_min_points=centroid_min_points,
            dataset_name=dataset_name,
            dataset_id=dataset_id,
            obs_keys=obs_keys,
            vector_field_default=vector_field_default,
        )
        self._server.start_background()
        self.port = self._server.port
        logger.info(f"Started AnnData server at {self._server.url}")

    @property
    def viewer_url(self) -> str:
        """Get the full viewer URL with anndata flag."""
        from urllib.parse import urlencode

        base = self._get_client_server_url()
        query = urlencode(
            {
                "jupyter": "true",
                "viewerId": self._viewer_id,
                "viewerToken": self._viewer_token,
                "anndata": "true",
            }
        )
        return f"{base}/?{query}"

    def _get_pre_display_html(self) -> str | None:
        """Show warning about AnnData mode being slower."""
        return """
        <div style="background: #fff3cd; border: 1px solid #ffc107; border-radius: 4px;
                    padding: 10px; margin-bottom: 10px; color: #856404;">
            <strong>Note:</strong> Loading data directly from AnnData. This is slower than
            using <code>prepare</code>. For production use, consider exporting
            your data first.
        </div>
        """

    def __repr__(self) -> str:
        status = "running" if self._server and self._server.is_running() else "stopped"
        return f"AnnDataViewer(port={self.port}, status={status})"


def show_anndata(
    data: str | Path | anndata.AnnData,
    height: int = 600,
    *,
    latent_key: str | None = None,
    gene_id_column: str | None = None,
    normalize_embeddings: bool = True,
    centroid_outlier_quantile: float = 0.95,
    centroid_min_points: int = 10,
    dataset_name: str,
    dataset_id: str,
    obs_keys: Sequence[str] | None = None,
    vector_field_default: str | None = None,
    client_server_url: str | None = None,
    web_source_url: str = CELLUCID_WEB_URL,
    web_cache_dir: str | Path | None = None,
    allowed_hosts: Sequence[str] | None = None,
) -> AnnDataViewer:
    """
    Quickly display an AnnData object, h5ad file, or zarr store in a notebook.

    This is the easiest way to visualize AnnData directly without running
    prepare first. However, it's slower than using pre-exported data.

    Args:
        data: AnnData object, path to h5ad file, or path to zarr directory. A
            leading ``~`` in a path is expanded to the user's home directory.
        height: Height of the viewer in pixels.
        latent_key: Explicit key in ``obsm`` for the latent space.
        gene_id_column: Exact column in ``var`` containing gene identifiers.
            If None, identifiers come from ``var.index``.
        normalize_embeddings: Whether to normalize embeddings into the viewer
            coordinate range.
        centroid_outlier_quantile: Quantile used for categorical centroids.
        centroid_min_points: Minimum category size used for centroid computation.
        dataset_name: Exact non-empty human-readable name without surrounding
            whitespace or control characters.
        dataset_id: Portable 1-180 byte ASCII identifier beginning with a
            letter or digit and otherwise using letters, digits, ``.``, ``_``,
            or ``-``; no trailing ``.`` or Windows device name.
        obs_keys: Exact ``obs`` columns to serve, in this order. If None, every
            column is served. Naming the columns leaves out a column Cellucid
            cannot serve, such as a ``datetime64`` collection date.
        vector_field_default: Exact field id required when multiple UMAP
            vector fields exist.
        client_server_url: Exact browser-reachable data server base URL.
        web_source_url: Origin publishing the web asset inventory.
        web_cache_dir: Directory holding the active verified web build.
        allowed_hosts: Extra ``Host`` names the data server answers to.
            Required when ``client_server_url`` points at a reverse proxy such
            as ``jupyter-server-proxy``: name that proxy's host here.

    Returns:
        AnnDataViewer instance for interaction via hooks.

    Supported formats:
        - In-memory AnnData objects
        - .h5ad files (HDF5-based, lazy loading via backed mode)
        - .zarr directories materialized by ``anndata.read_zarr``

    Example:
        >>> from cellucid import show_anndata
        >>> viewer = show_anndata(
        ...     adata,
        ...     dataset_name="Example",
        ...     dataset_id="example",
        ... )
        >>>
        >>> @viewer.on_selection
        ... def handle(event):
        ...     subset = adata[event['cells']]
        ...     sc.pl.violin(subset, 'gene1')

        >>> # With custom options
        >>> viewer = show_anndata(
        ...     adata,
        ...     latent_key="X_pca",
        ...     dataset_name="Example",
        ...     dataset_id="example",
        ...     height=800,
        ... )

    Note:
        For production use or sharing, consider using prepare
        to create optimized binary files, then use show() to display them.
    """
    return AnnDataViewer(
        data=data,
        height=height,
        auto_open=True,
        latent_key=latent_key,
        gene_id_column=gene_id_column,
        normalize_embeddings=normalize_embeddings,
        centroid_outlier_quantile=centroid_outlier_quantile,
        centroid_min_points=centroid_min_points,
        dataset_name=dataset_name,
        dataset_id=dataset_id,
        obs_keys=obs_keys,
        vector_field_default=vector_field_default,
        client_server_url=client_server_url,
        web_source_url=web_source_url,
        web_cache_dir=web_cache_dir,
        allowed_hosts=allowed_hosts,
    )
