# Analysis mode: Gene Signature (Gene Signature Score)

**Audience:** wet lab + computational users (gene program scoring)  
**Time:** 20–40 minutes  
**What you’ll learn:**
- What a gene signature score represents (in Cellucid terms)
- How to enter gene sets: typed, pasted, imported or taken from marker genes
- What happens to genes the dataset does not have
- The five scoring methods, with their formulas
- How rescaling and the plot change the picture, never the scores

**Prerequisites:**
- A dataset loaded
- Gene expression available
- One or more gene lists (“signatures”) to score

---

## What a gene signature is (user-facing)

A gene signature is a curated gene list meant to represent a biological state,
pathway, or cell program.

Cellucid turns each signature into a **per-cell score**:
- one number per cell and gene set,
- compared across the pages you select.

Intuition for wet lab users:
- a high score means “the genes of the program are high in this cell”,
- a low score means “the program is not expressed strongly in these cells”.

---

## The form

| Where | Options |
| --- | --- |
| Essentials | **Gene sets**, then **Pages** (with **Rest of \<page\>**) |
| More options → Scoring | **Method**; **Control genes**, **Expression bins**, **Random seed** (Module score); **Max rank** (UCell) |
| More options → Display | **Plot**, **Rescale**, **Sets** |

Typing never runs anything: press **Run**. The line above Run says what is
missing (`Add genes to a set`, `Select at least one page`, `None of the genes are
in this dataset`). After a run, changing the genes, the pages or a scoring option
marks the result **Stale · Run**; the display options redraw at once.

### Defaults

<!-- options:signature -->
| Label | Default |
|---|---|
| Gene sets | One empty set |
| Method | Mean |
| Control genes | 50 |
| Expression bins | 25 |
| Random seed | 0 |
| Max rank | 1500 |
| Plot | Violin plot |
| Rescale | None |
| Sets | Side by side |
<!-- /options:signature -->

---

## Gene sets

Each set is a row with its name (rename it in place), its genes as chips, and a
field to add more:

- **Type** a gene and choose it from the suggestions, or **paste** a list: names
  may be separated by commas, semicolons, spaces, tabs or new lines, and each is
  kept once.
- **Import** a file: a gene list (`.txt`), a table (`.csv`, `.tsv`) or a GMT file
  (`.gmt`, as MSigDB writes it, one set per line). A table can be one column (one
  set), a set column and a gene column (`set,gene`, or scanpy’s `group,names`), or
  one column per set named by its header.
- **From markers** — one set per group of the last Marker genes result (shown
  once there is one); Marker genes’ **Score as signature** does the same from the
  result.
- **Add set** scores several programs in one run; **×** on a chip removes a gene,
  **Clear** empties a set.

A row’s header reads `12 genes`, or `12 genes · 2 not found`.

### Genes not in the dataset

A name is matched to the dataset’s genes exactly, or else to the one gene that
matches it ignoring case (`Cd3e` → `CD3E`; the chip says what it matched). A name
that matches nothing, or several genes ignoring case, is **not found**:

- its chip is amber, and **Remove not found** drops all of them;
- it is left out of the score — never an error — and the result names it in a
  **n not found** chip;
- a set with no gene of the dataset is not scored and is named in the result.

No alias table is applied. An export can publish fewer genes than its source
measured, so a signature gene can be absent rather than misspelled; see
{doc}`../d_fields_coloring_legends/07_genes_in_the_built_in_samples`.

---

## Methods

Each cell `i` and set `G` (the set’s genes found in the dataset), on the
expression matrix in use (see {doc}`../b_data_loading/13_data_mapping`):

| Method | Score | Reference |
| --- | --- | --- |
| **Mean** (default) | the mean of `x_ig` over `g ∈ G` | — |
| **Median** | the median of `x_ig` (the mean of the two middle values for an even count) | — |
| **Mean z-score** | the mean of `(x_ig − μ_g) / σ_g`, with `μ_g`, `σ_g` (ddof 1) over every cell of the dataset | — |
| **Module score (scanpy)** | the mean of `x_ig` over `G` less the mean over control genes | scanpy `tl.score_genes` (Satija et al. 2015) |
| **Rank-based (UCell)** | `1 − (Σ_{g∈G} min(r_ig, R) − S(S+1)/2) / (S·R − S(S+1)/2)` | Andreatta & Carmona 2021; the UCell package's `u_stat` (UCell 2.x) |

Details:
- **Mean** and **Median** use the cell’s finite values; a cell with none gets no
  score and is left out of the statistics and plots.
- **Mean z-score** leaves out a gene that is constant over every cell (it has no
  z-score) and names it in a **n constant** chip.
- **Module score** bins every gene by the rank of its mean over all cells into
  **Expression bins** bins (25), draws **Control genes** genes (50) from each bin
  a set gene falls in, and removes the set’s own genes from the draw — scanpy
  1.10’s rule. The draw is seeded by **Random seed** (0), so the same seed gives
  the same controls; it uses its own generator, so a control set drawn from a bin
  larger than **Control genes** is not bit-identical to scanpy’s. The expanded
  view names the controls.
- **UCell** ranks every gene in each cell, highest first, ties sharing their
  average rank; `r_ig` is gene `g`’s rank, `S` the set size and `R` **Max rank**
  (1500, at most the dataset's genes), at or past which a rank counts as `R`.
  This is the current UCell package's score (R `u_stat`): 1 when the set's
  genes are the cell's top `S`, 0 when none is ranked within `R`. (The 2021
  paper's form counted such ranks as `R + 1` over `S·R`, so a cell with no set
  gene in range scored `(S − 1) / (2R)` instead of 0.) A set needs
  `R > (S + 1) / 2`. Scores lie in 0–1, comparable across cells and datasets.

What a run reads:
- **Mean**, **Median**: each set gene once, as its non-zero entries.
- **Mean z-score**, **Module score**: also each needed gene’s mean and spread over
  all cells — for Module score, every gene’s, read once per dataset and kept.
- **UCell**: every gene once, since ranks need them all. The cells are split
  across the analysis workers.

:::{note}
Scores are on the scale of the expression matrix in use: log-normalised values
give log-scale scores, counts give count-scale scores. UCell is the method least
sensitive to that choice.
:::

---

## Display

- **Plot**: Violin plot (default), Box plot, or Histogram.
- **Rescale**: **None** (default); **z-score**, `(x − μ) / σ` (population σ); or
  **0–1**, `(x − min) / (max − min)` — per set, over every page it was scored on,
  so pages stay comparable and two sets are never put on one scale. Rescaling
  changes the picture only; the CSV holds the scores as computed.
- **Sets**: **Side by side** (one distribution per page and set, `Page · Set`) or
  **One at a time** (buttons above the plot choose the set).

---

## Outputs and interpretation

- A distribution per page (and set), with chips for genes not found and constant
  genes.
- **Expanded view:** summary statistics per page; with two or more distributions,
  the same tests Detailed uses for a continuous variable (Welch’s t-test and
  Mann–Whitney U for two, with effect sizes); each set’s genes as chips, amber for
  those left out, and Module score’s controls.

Compare effect sizes before p-values, and check whether a signature is a program
or a proxy for cell type composition. The tests are exploratory; see
{doc}`01_analysis_mental_model`.

---

## Export (CSV)

`gene_signature_scores.csv`, one line per scored cell and set:
`page,set,cell_index,score`. `cell_index` is the cell's row in the whole
dataset in its original order (0-based, also with a subsample open), so the
file joins back to `adata.obs`. See {doc}`09_exporting_analysis_results`.

---

## Edge cases and pitfalls

- **Huge signatures (hundreds of genes)**: each gene is one read, so the time grows
  with the gene count over a slow link. The genes (and a Module score's
  controls) are read a batch at a time, as many as the analysis memory holds,
  and each is added to every page's scores and let go, so a large signature on
  a large dataset never holds every gene at once; a gene in two sets is read
  for each.
- **Housekeeping-dominated signatures**: scores track library size and QC rather
  than biology; try UCell, or Module score.
- **Comparing sets**: Mean scores of two sets are on different footings when their
  genes differ in level; Module score and UCell correct for that.

---

## Troubleshooting (Gene Signature)

### Symptom: “n not found”

The names are not the dataset’s. Read a few names in the field selector and match
that vocabulary; **Remove not found** clears them from the set.

### Symptom: Run is disabled with `None of the genes are in this dataset`

No set has a gene of this dataset. Check the vocabulary (symbols vs Ensembl ids,
human vs mouse names) and that the dataset has gene expression (see
{doc}`../b_data_loading/index`).

---

## Interface reference

```{figure} ../../../_static/screenshots/analysis/gene-signature.png
:alt: The Gene signature panel with one gene set holding Ins1, Ins2 and Iapp, both pages selected, and the signature score of each page drawn beneath.
:width: 488px

Three genes, both pages, mean scoring, no rescaling. The two violins are the
score distribution for each page — which is the whole point of the mode: a
signature is only interpretable **relative to** something else.
```

```{figure} ../../../_static/screenshots/analysis/gene-signature-expanded.png
:alt: The expanded Gene signature view of Ins1, Ins2 and Iapp: the score distribution of each page, the options beside it, a summary table per page and the two-group tests with their effect sizes.
:width: 1440px

Expanded Gene Signature. With exactly two pages selected the statistical panel
shows Welch's t-test and Mann-Whitney U side by side, each with its own effect
size — Cohen's *d* and rank-biserial *r*.
```

:::{tip}
If the two tests disagree, the distributions usually differ in shape rather than
in location, or a handful of extreme cells is carrying the mean.
:::

---

## Next steps

- {doc}`08_analysis_mode_genes_panel` (marker discovery across categorical groups)
- {doc}`09_exporting_analysis_results` (export signature results)
- {doc}`10_troubleshooting_analysis` (missing genes, missing expression)
