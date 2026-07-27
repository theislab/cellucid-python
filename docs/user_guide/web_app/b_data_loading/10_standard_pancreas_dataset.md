# Standard Pancreas sample dataset

Cellucid's public sample catalog includes the pancreas dataset used by the
[scVelo pancreas vignette](https://scvelo.readthedocs.io/en/stable/vignettes/Pancreas.html).
It is a complete known-good example for dimensional navigation, categorical and
continuous metadata, gene coloring, connectivity, and RNA-velocity rendering.
The catalog id is exactly `pancreas`; Suo remains the catalog default.
Its published identity contains exactly 3,696 cells and 3,753 genes.

## Load it

1. Open an ordinary Cellucid startup and dismiss the welcome overlay.
2. Open **Session**.
3. Under **Sample datasets**, select
   **Pancreatic endocrinogenesis (scVelo)**.
4. Wait until the loaded-dataset name changes to
   **Pancreatic endocrinogenesis (scVelo)**. The compact Session statistics
   round both counts to **4K**; **Coloring & Filtering** reports the exact
   **Showing all 3,696 points**, while the published identity and provenance
   pin the exact inventory of **3,753 genes**.

The sample opens in its declared default 3D embedding with **Orbit**
navigation. Selecting 1D or 2D changes the dimension-specific default to
**Planar**, unless you explicitly chose another navigation mode. No camera path
is created or played automatically.

## Verify the scientific controls

The public dataset exposes:

- source-owned categorical observations `clusters_coarse` (5 categories) and
  `clusters` (8 categories);
- the derived categorical convenience alias `cell_type` (8 categories), whose
  category order and per-cell codes exactly equal source-owned `clusters`;
- source-owned continuous observations `S_score` and `G2M_score`;
- 3,753 ordered, source-marked, nonconstant genes for on-demand coloring;
- 33,476 exact weighted connectivity edges;
- independently declared 1D, 2D, and 3D embeddings; and
- one dimension-matched vector field, `velocity_umap`, in 1D, 2D, and 3D.

To inspect velocity, open **Visualization**, keep **Render mode** set to
`Points`, find **Vector Field Overlay**, and enable **Show overlay**. The field
selector reports `velocity_umap`. The overlay is a qualitative view of the
projected stochastic velocity; it does not replace scVelo's numerical
diagnostics or uncertainty checks.

```{figure} ../../../_static/screenshots/vector_field_velocity/pancreas-velocity-3d.png
:alt: Cellucid showing the Pancreatic endocrinogenesis (scVelo) sample in its 3D Orbit view, colored by clusters, with the velocity_umap overlay enabled and vector controls visible.
:width: 1440px

The standard Pancreas sample in Build 2026-07-27.1: `clusters` supplies eight
categories, **Show overlay** is enabled, `velocity_umap` is selected, and the
dimension-owned controls report **3D** with **Orbit** navigation.
```

## What is source-owned and what is derived

The source is the commit-pinned E15.5 H5AD from the scVelo notebooks:

```text
https://raw.githubusercontent.com/theislab/scvelo_notebooks/f6cad69fd509be44c8453205309f4c3c3c37ba34/data/Pancreas/endocrinogenesis_day15.h5ad
```

Its exact 2D UMAP, 50-dimensional PCA, cell and gene order, four source-owned
observation fields (`clusters_coarse`, `clusters`, `S_score`, and `G2M_score`),
highly-variable-gene marker, expression values, and weighted connectivity graph
are retained as source-owned inputs. The public `cell_type` field is a derived
convenience alias of `clusters`, with exactly the same category order and
per-cell codes; it is not a fifth source-owned field. The 1D and 3D UMAPs are
independently optimized deterministic embeddings; neither is sliced, padded,
or copied from the 2D coordinates. Velocity is computed on a separate
analysis-only working copy and projected independently onto each exact
embedding. It never replaces the public expression inventory.

The source contains 27,998 ordered genes. The network payload intentionally
ships the exact 3,753 columns that are both marked
`highly_variable_genes == "True"` and nonconstant. This removes unmarked and
constant per-gene requests without changing cell order, selected-gene order, or
the selected expression values.

## Provenance and reproducibility

The canonical provenance is
[`sources/pancreas.json`](https://github.com/theislab/cellucid-datasets/blob/main/sources/pancreas.json)
in `cellucid-datasets`. It pins:

- the source commit, URL, SHA-256, byte size, and biological citation;
- the exact Python and scientific-package environment;
- deterministic embedding and velocity parameters;
- semantic digests for coordinates, vectors, metadata, genes, expression, and
  connectivity;
- the exact Cellucid 0.9.1 producer-source digest; and
- the complete 3,774-file, 2,477,205-byte generation digest
  `268a36b0c62f1e872bc4ebc271283fe09b9036c2b8db2b7bbf2d7d737556c865`.

The source data correspond to
[GSE132188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132188).
For scientific use, cite
[Bastidas-Ponce et al. (2019)](https://doi.org/10.1242/dev.173849) in addition
to the tools used in your analysis.

## Troubleshooting

- If the sample is not listed, the configured exports host is not serving the
  current `datasets.json`; include the footer **Build** value in a bug report.
- If the overlay controls are absent, confirm that the loaded dataset name is
  **Pancreatic endocrinogenesis (scVelo)** and that `Points` is the active
  render mode.
- If **Show overlay** disables after a dimension change, wait for that
  dimension's embedding and `velocity_umap` payload to finish loading before
  retrying.
- If points or flow look biologically surprising, compare the exact dimension
  with the scVelo analysis and inspect velocity confidence before interpreting
  directionality.

For the complete overlay control reference, see
{doc}`../i_vector_field_velocity/index`. For general sample and remote loading
failures, see {doc}`08_troubleshooting_data_loading`.
