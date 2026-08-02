"""
Export raw dataframes/arrays to files used by the WebGL viewer.

Includes memory/disk optimization features:
- Quantization for continuous data (var/gene expression and obs continuous)
- Categorical obs codes stored at the caller's declared width, which every
  field's category count is required to fit
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

The module is a package. Its public surface is exactly the two functions below;
everything else is split by responsibility across private submodules:

``_contracts`` (in the parent package)
    Identity, display-text, and output-directory rules, standard library only.
``_manifest``
    The compact_v1 manifest format and the payload paths it declares.
``_arrays`` / ``_categories`` / ``_quantization`` / ``_centroids``
    Numeric input contracts, label identity and storage, the continuous codec,
    and the measured centroids and outlier quantiles.
``_binary_io``
    Deterministic, atomic payload writes.
``_filesystem`` / ``_transaction`` / ``_locking``
    Durable path primitives, the write-ahead publication transaction, and the
    exclusive writer lock.
``_generation``
    The writer that builds one complete generation in a staging directory.
``_export`` / ``_catalog``
    :func:`prepare` and :func:`generate_datasets_manifest`.
"""

from ._catalog import generate_datasets_manifest
from ._export import prepare

__all__ = ["generate_datasets_manifest", "prepare"]
