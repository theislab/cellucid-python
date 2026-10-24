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
  like a sample; otherwise the first categorical field named like a cell type
  or cluster; otherwise the
  first not named like a sample.
- **Genes**: one or more named gene sets. Type or paste genes (commas,
  semicolons, spaces, tabs or new lines all separate), import a `.txt`,
  `.csv`, `.tsv` or `.gmt` file, or take **From markers** (one set per group of
  the shown Marker genes result). Each set is a bracketed, labelled block of
  columns.

Gene names match exactly, else by case when only one gene matches (`cd3e`
becomes `CD3E`). Genes the dataset does not have are left out and listed on a
`not found` chip; they never stop the run. With another **Expression matrix**
the genes are matched against that matrix's own names (`raw` has its own), so a
gene only raw has is found there
([Expression matrix](08_analysis_mode_genes_panel.md#expression-matrix)).

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

On a phone the plot keeps its size and its frame scrolls.

---

## Results and export

- **Sidebar**: the plot, a `not found` chip when genes were left out, and
  **Score as signature**, which opens Gene signature with one signature per
  gene set.
- **Expanded view** (⤢): the plot at full size, with a summary of the run.
- **CSV** (in the expanded view): `dot_plot.csv`, one row per group and gene:
  `group, gene, gene_group, mean, fraction_expressing, n_cells` (unscaled means).

Only the listed genes are read; changing a display option never reads again.
Sessions save the options and the shown result.

---

## Options and defaults

<!-- options:dotplot -->
| Label | Default |
|---|---|
| Group by | The coloring field when it is not named like a sample, else the first categorical field named like a cell type or cluster, else the first not named like a sample |
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
