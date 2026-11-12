# Analysis mode: Dot plot

**Audience:** wet lab + computational users (marker panels, cluster annotation)
**What you’ll learn:**
- What the color and the size of a dot mean
- How to enter gene sets, and what happens to genes the dataset lacks
- How the values compare to scanpy’s `pl.dotplot` and `pl.matrixplot`

**Prerequisites:**
- A dataset loaded with gene expression
- A categorical field to group the cells by

---

## What a dot plot shows

One dot per gene and group:

- its **color** is the gene’s mean expression in the group’s cells;
- its **size** is the fraction of the group’s cells that express the gene
  (value above the cutoff, 0 by default).

A cell with no value of a gene — a MuData modality that did not measure it
(CITE-seq protein on some cells, a multiome cell without ATAC), a layer holding
`NaN` — is left out of that gene’s mean, fraction and violin, as differential
expression and marker genes leave it out: they are over the group’s cells
**with a value**. A `cells without a value left out` chip says when a gene has
such cells, the hover names a group’s cells and how many hold a value, and the
expanded view’s summary counts them.

A dot plot answers “which groups express these genes, and how uniformly?” —
the usual way to check a marker panel against clusters. **Matrix plot** shows
the same means as colored cells, without the sizes; **Stacked violin** shows
each gene’s whole distribution in each group.

As in scanpy, the dot area grows as the fraction to the power 1.5, scaled so the
largest fraction, rounded up to the next 10%, fills the largest dot; the size
legend shows 20 to 100% of that fraction.

---

## Group by and genes

- **Group by**: the categorical field whose groups are the rows. The default is
  the field the cells are colored by, when it is categorical and not named
  like a sample; otherwise the best field named like a cell type, then another
  annotation, a lineage, then clusters (the ranks of
  {doc}`../d_fields_coloring_legends/02_field_selector_ux`); otherwise the
  first not named like a sample.
- **Genes**: one or more named gene sets. Type or paste genes (commas,
  semicolons, spaces, tabs or new lines all separate), import a `.txt`,
  `.csv`, `.tsv` or `.gmt` file, or take **From markers** (one set per group of
  the shown Marker genes result). Each set is a bracketed, labeled block of
  columns.

Gene names match exactly, else by case when only one gene matches (`cd3e`
becomes `CD3E`). Genes the dataset does not have are left out and listed on a
`not found` chip; they never stop the run. With another **Expression matrix**
the genes are matched against that matrix's own names (`raw.X` has its own), so a
gene only raw has is found there
([Expression matrix](08_analysis_mode_genes_panel.md#expression-matrix)).

```{figure} ../../../_static/screenshots/analysis/dotplot-result.png
:alt: The Dot plot form with Group by clusters and one gene set, Set 1, holding eight pancreatic markers as chips (Krt19, Sox9, Neurog3, Fev, Ins1, Gcg, Sst, Ghrl), the box to add or paste more, Clear, Add set, Import and Export, More options, the Run button, and under it the line the result opens with: 8 genes, 8 groups of clusters and the cell count.
:width: 468px

**Group by** `clusters` and one gene set of eight pancreatic markers, entered
as chips. After **Run** the line under it says what was drawn: genes, groups
and cells.
```

---

## Groups and cells

Under **More options → Groups and cells**:

- **Groups to show**: all groups, or a chosen few.
- **Cells**: all cells, the cells the view shows, or the cells of one highlight
  page. A group with no cell in scope is left out.

---

## Values

Under **More options → Values**:

| Option | Choices | scanpy |
|---|---|---|
| Expression matrix | Mapped (the matrix the data mapping shows); X, a layer or raw for this plot alone | `layer=` / `use_raw=` |
| Color | Mean expression; mean scaled 0–1 per gene; mean as z-score per gene | `standard_scale=None` / `'var'`; z as `scipy.stats.zscore` over the groups |
| Expressed above | a number, 0 by default | `expression_cutoff` |
| Mean of expressing cells only | off / on | `mean_only_expressed` |

The values averaged are those **Data mapping › Expression values** sets for
the whole app (as stored by default; normalized to 10,000 counts per cell
and/or `log1p`: scanpy's dot plot of `normalize_total(target_sum=1e4)` +
`log1p` data with both on), so the dots agree with the coloring; see
[Values tested](08_analysis_mode_genes_panel.md#values-tested). The expanded
view's **Values** row and `about.values` say which, and the color bar names
them when they are not as stored (`Mean expression · log1p(norm.)`), as the
CSV's means column does (`mean · log1p(norm.)`) — the values of the run, so a
stale result keeps its own.

The mean, fraction and mean of expressing cells match `pl.dotplot`’s
`dot_color_df` and `dot_size_df` (the web app’s test suite holds them to it,
including a layer with negative values, where zeros are not expressed and the
mean includes the negatives).

---

## Layout

Under **More options → Layout** (these redraw without reading again):

- **Plot**: dot plot, matrix plot or stacked violin. A stacked violin draws,
  in each cell, the gene's distribution in the group: its values run up the
  cell (40 equal bins over the gene's range in the plotted cells, zero
  included) and the share of the group's cells at each value runs across it,
  filled with the cell's color value. Each gene's violins share one value
  scale and one width scale.
- **Genes down, groups across**: swap the axes.
- **Group order**: as in the field, largest group first, or by a dendrogram of
  the groups’ mean profiles (correlation distance, average linkage).
- **Color map**: Reds, Viridis or Blues.
- **Largest dot**: automatic (fitted to the room each row has) or 12–32 px.

Every gene and group is named, never one name on another: a panel or window too
small for them (many groups, a phone) draws the plot as large as its names need
and scrolls.

---

## Results and export

- **Sidebar**: the plot, a `not found` chip when genes were left out, and
  **Score as signature**, which opens Gene signature with one signature per
  gene set.
- **Expanded view** (**Expand**): the plot at full size, with a summary of the run, and
  the Values and Layout options beside it (changed there, they change the form
  too, and redraw without running).
- **CSV** (in the expanded view): `dot_plot.csv`, one row per group and gene:
  `group, gene, gene_group, mean, fraction_expressing, n_cells, n_measured`
  (unscaled means, the column named with their values when they are not as
  stored, `mean · log1p(norm.)`; `n_measured` the group’s cells with a value of
  the gene, which the mean and the fraction are over).

```{figure} ../../../_static/screenshots/analysis/dotplot-expanded.png
:alt: The expanded Dot plot of the eight pancreatic markers across the eight clusters: one row per cluster, dot color the mean expression and dot size the fraction of cells expressing it, with its color bar and size legend, Export, and a summary table below.
:width: 1440px

The expanded dot plot: one row per cluster and one column per gene under its
set's bracket, color the mean expression and size the share of cells
expressing it. Each marker lights up its own cell type (Ins1 in Beta, Gcg in
Alpha, Sst in Delta, Ghrl in Epsilon, Neurog3 in the endocrine progenitors).
```

Only the listed genes are read; changing a display option never reads again.
Sessions save the options and the shown result.

---

## Options and defaults

<!-- options:dotplot -->
| Label | Default |
|---|---|
| Group by | The coloring field when it is not named like a sample, else the best field named like a cell type, annotation, lineage or cluster, else the first not named like a sample |
| Genes | One empty set |
| Groups to show | All |
| Cells | All cells |
| Expression matrix | The mapped matrix |
| Color | Mean expression |
| Expressed above | 0 |
| Mean of expressing cells only | Off |
| Plot | Dot plot |
| Genes down, groups across | Off |
| Group order | As in the field |
| Color map | Reds |
| Largest dot | Automatic |
<!-- /options:dotplot -->
