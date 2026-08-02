# Which genes the built-in samples publish

**Audience:** anyone exploring the built-in sample datasets, especially wet lab
users looking for a specific marker
**Time:** 5–10 minutes
**What you’ll learn:**
- What the **Gene Expression** search is searching in a built-in sample
- Why a gene that is genuinely in the underlying study can be absent from a sample
- How many genes each sample publishes, and where to check that yourself
- What to do when you need a gene a sample does not carry

---

## The short version

The gene list in **Coloring & Filtering → Gene Expression** is exactly the list
of gene names the loaded dataset published, and nothing else. Cellucid adds no
aliases and looks nothing up.

The four human samples — **Suo**, **Garcia**, **He**, and **Kanemaru** — publish
a subset of the genes their prepared source file contains, because that file
identifies genes by Ensembl accession and only the accessions that could be
translated to a gene symbol are published. So a real gene, measured in the real
study, can produce **No gene matches** here. That is a property of how these
four demonstration samples were prepared, not a limit of Cellucid and not a
defect in the studies.

---

## What you see when a gene is not there

Type a name the loaded dataset does not publish and the dropdown replaces the
result list with three lines:

- **No gene matches**, followed by your query in quotes;
- *This dataset publishes N gene names, chosen when it was prepared and possibly
  a subset of the source data. Check for typos too.*, where *N* is the loaded
  dataset’s own published gene count; and
- a **Why a gene may be missing** link, which opens
  {doc}`05_troubleshooting_fields_legends`.

```{figure} ../../../_static/screenshots/web_app/gene-search-no-match.png
:alt: A field row labelled GENE EXPRESSION with CD19 typed into its search box, and an open dropdown below reading No gene matches, then a sentence stating how many gene names this dataset publishes and that it was chosen when the dataset was prepared and is possibly a subset of the source data with a reminder to check for typos, then a link reading Why a gene may be missing.
:width: 480px

The empty state on the Pancreas sample, searching for a human B-cell marker in
a mouse dataset. The number in the middle line is **this** dataset's own
published gene count, so it is the fastest way to tell a typo apart from an
absent gene.
```

For contrast, a query that does match simply lists the matches:

```{figure} ../../../_static/screenshots/web_app/gene-search-dropdown-open.png
:alt: The same field row with the two letters Rp typed into the search box and a dropdown below listing gene names that contain those letters, one per row.
:width: 480px

Matching is a case-insensitive substring test over the published names only.
There is no alias table and no second identifier to try instead.
```

The viewer states that count and stops there on purpose: an export records the
names it published and nothing about what preparation left behind, so the app
cannot know whether your query was dropped, misspelled, or never measured.

Before concluding a gene is absent, note how the search behaves: it is a
**substring** match and it ignores letter case, so `gata3` finds `GATA3` and the
shorter `gata` matches every published name containing it. At most 100 matches
are listed at once, with a `...and N more` line when there are further ones.
Matching in the Analysis panel is stricter — see
[Where else this list is used](#where-else-this-list-is-used).

---

## The four human samples: two reductions, not one

### Reduction 1 — the prepared source file already holds 8,192 genes

Suo, Garcia, He, and Kanemaru were exported from integrated files produced
upstream of this project, and each of those files holds exactly **8,192**
features rather than a whole transcriptome. The selection to 8,192 was made by
that integration pipeline, before any Cellucid step ran, and no record of the
criterion travels with the file: the gene axis arrives carrying accessions and
no annotation columns at all.

That 8,192 is a per-dataset selection, not a fixed panel shared by the four. The
four lists overlap only partially: any two of them share between 46.9% and 68.8%
of their accessions, only 2,385 accessions are common to all four, and the four
together name 15,637 distinct accessions. Each sample was reduced from its own
data.

### Reduction 2 — only genes that resolve to a symbol are published

Cellucid publishes **one name per gene and keeps no second identifier**, so
whatever a dataset publishes is the only string its gene search can match. An
Ensembl accession such as `ENSG00000004948.16` is not a name a reader would type,
so the export step for these four samples translates accessions into gene
symbols and publishes only the ones it can name:

1. An identifier that is already a symbol passes through unchanged.
2. An `ENSG…` accession is translated through one fixed local symbol table — no
   network lookup — trying the versioned accession first and then the
   version-stripped base.
3. A gene is published **only if that translation yields something that is not
   itself an accession**.

Step 3 fails in two ways. The accession may be absent from the symbol table
entirely, or it may be present but mapped to another accession — the table holds
1,627 entries that return an accession — so the lookup succeeds and still yields
no name. Both count as unnamed, and unnamed genes are left out rather than
published under an accession.

### What each sample publishes

| Sample | Catalog id | Genes in the prepared file | Published | Left out |
| --- | --- | --- | --- | --- |
| Developing human immune system (Suo) | `suo` | 8,192 | 5,103 | 3,089 |
| Developing human gonad (Garcia) | `garcia` | 8,192 | 5,754 | 2,438 |
| Developing human lung (He) | `he` | 8,192 | 5,152 | 3,040 |
| Developing human heart (Kanemaru) | `kanemaru` | 8,192 | 3,691 | 4,501 |

Most of the genes left out are absent from the symbol table altogether — 2,067
for Suo, 1,654 for Garcia, 2,091 for He, and 3,802 for Kanemaru. The remainder
are in the table but map to another accession: 1,022, 784, 949, and 699
respectively.

:::{note}
Every published count above is the `stats.n_genes` value in that sample’s own
`exports/<id>/dataset_identity.json`, which is also the number of entries in its
`var_manifest.json` and the number of payload files in its `var/` directory. Each
sample’s `<id>_2_export_cellucid.ipynb` in
[`cellucid-datasets`](https://github.com/theislab/cellucid-datasets) prints the
full breakdown in its committed output. Read those rather than this page if you
are verifying a build; they change when a sample is regenerated.
:::

---

## Pancreas is filtered for a different reason

The **Pancreatic endocrinogenesis (scVelo)** sample is mouse, and its source
already indexes genes by MGI symbol, so no translation step applies to it and
none of the above is why its gene axis is smaller than its source.

Pancreas publishes the 3,753 of its source’s 27,998 genes that are both marked
`highly_variable_genes == "True"` and nonconstant. That is a separate, pinned
selection made in its build script, and it is described in
{doc}`../b_data_loading/10_standard_pancreas_dataset`. Do not read the symbol
rule above onto it.

---

## Where else this list is used

The published gene axis is the same list every part of the app draws on, so the
same absence shows up beyond the search box:

- coloring by a gene, and every legend and filter that follows from it;
- {doc}`../h_analysis/06_analysis_mode_differential_expression_de` and
  {doc}`../h_analysis/08_analysis_mode_genes_panel`, whose results can only ever
  cover published genes; and
- {doc}`../h_analysis/07_analysis_mode_gene_signature`, where a pasted gene list
  is matched **exactly** and case-sensitively against the published names, with
  no aliasing — a stricter rule than the search box’s substring match.

A marker table or a signature score computed on one of these four samples is
therefore computed over its published genes, not over the study’s full gene
axis. That is worth stating when you report a result from a demonstration
dataset.

---

## What to do if you need a gene that is not there

**First, confirm it is really absent.** Search a short, distinctive substring
rather than the full name, and if the dropdown ends in `...and N more`, narrow
the query until the whole match list is visible.

**If it is absent, no in-app setting brings it back.** The gene panel is fixed by
the export; the viewer has nothing to re-enable, and reloading the sample
produces the same list.

**Work from a dataset that carries the gene.** The built-in samples are fixed
demonstration exports, so a full gene axis has to come from data you prepare
yourself:

- Load your own `.h5ad` or `AnnData` through {doc}`../b_data_loading/04_server_tutorial`
  or {doc}`../b_data_loading/05_jupyter_tutorial`, which needs no export step.
- Or prepare an export of your own with
  {doc}`../../python_package/c_data_preparation_api/index`, choosing
  `var_gene_id_column` so genes are named the way your readers will type them.
  {doc}`../../python_package/c_data_preparation_api/05_var_gene_metadata` and
  {doc}`../../r_package/c_data_preparation_api/05_var_gene_metadata` cover that
  choice, including resolving accessions to symbols in your own preprocessing.

**If you want these exact studies at full depth,** start from the upstream
accession. {doc}`../b_data_loading/12_sample_dataset_provenance` names the study,
URL, and citation behind each sample and states what is and is not reproducible
about each one.

---

## Next steps

- {doc}`02_field_selector_ux` — how the gene search, rename, delete, and restore
  behave.
- {doc}`05_troubleshooting_fields_legends` — the symptom-driven guide the in-app
  **Why a gene may be missing** link opens.
- {doc}`../b_data_loading/12_sample_dataset_provenance` — what to cite for each
  sample and how reproducible each one is.
