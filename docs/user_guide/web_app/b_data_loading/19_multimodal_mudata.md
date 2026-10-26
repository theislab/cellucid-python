# Multimodal data: MuData and the Modality select

**Audience:** anyone with CITE-seq, multiome or another experiment that
measured several kinds of features on the same cells  
**Time:** 5 minutes  
**What you’ll learn:** how to open a MuData file, how the **Modality** select
in Data mapping works, and what happens to cells a modality did not measure  
**Prerequisites:** none

A **modality** is one feature space measured on the cells: RNA, surface
proteins (ADT), chromatin peaks. Cellucid shows one modality's features at a
time; the cells, their fields, embeddings and neighbour graph are shared.

---

## Open a MuData file

**Open… › File…** opens an `.h5mu`; **Open… › Folder…** (or a drop) opens a
MuData Zarr folder, and **File…** a MuData `.zarr.zip`. Links work as for any
file ({doc}`14_open_from_link`). Files written by mudata 0.2 and later open,
`.h5mu` and Zarr v2 or v3 alike, read in place.

| In the MuData file | In Cellucid |
|---|---|
| the global `obs` | the fields |
| each modality's own `obs` columns (`rna` › `leiden`) | fields named `<modality>:<column>` (`rna:leiden`), as mudata names a pulled column |
| the global `obsm` (a joint embedding such as `X_wnn_umap`) | embeddings |
| each modality's `obsm` (`rna` › `X_umap`) | embeddings named `rna:X_umap` |
| the global `obsp`, each modality's `obsp` | neighbour graphs (`rna:connectivities`) |
| `mod/<name>` | one modality each, in the file's `mod-order` |
| the masks mudata keeps of which cells each modality measured | fields `<modality>:measured` |

mudata 0.4 and later keep each modality's columns in the modality only; the
`<modality>:<column>` fields show them, so a recent file's clusters appear as
they do in a file written by mudata 0.3. A MuData whose modalities share
features rather than cells (`axis` 1) is refused: Cellucid shows cells.

The dataset opens on its gene-expression modality: the first named `rna`,
`gex`, `gene expression`, `gene_expression` or `sct` (in any case), else the
first modality. A file written by `cellucid.prepare()` or `cellucid_prepare()`
may name another in its launch defaults (`uns['cellucid']['defaults']['modality']`).

---

## The Modality select

In the Session panel, open **Data mapping**, then **Genes**. A dataset of
several modalities shows **Modality** above **Expression**; a dataset of one
feature space never does. Each option names the modality, its feature count,
and how many cells it measured when that is not every cell:

```text
prot · 5 features · 36 of 44 cells
rna · 28 features
```

```{figure} ../../../_static/screenshots/data_loading/data-mapping-modality.png
:alt: The Data mapping group of a CITE-seq .h5mu with Genes open: a Modality select reading "rna · 28 features" above Expression, Gene names, Values and Highly variable genes, and Revert and Apply below.
:width: 488px

**Data mapping › Genes** of a CITE-seq `.h5mu`: **Modality** sits above
**Expression**, **Gene names**, **Values** and **Highly variable genes**,
which all follow the modality chosen. Choose one and click **Apply**.
```

Choose a modality and click **Apply**. Expression, Gene names, Values and
Highly variable genes then list that modality's matrices and columns, set to
its defaults, and gene search, colouring by a gene, Marker genes,
Differential expression and Dot plot read its features. Choosing back the
modality in use restores its choices.

A **cell the modality did not measure** has no value, not zero: coloured by a
protein, it is drawn as missing, and an analysis does not count it as an
expressed zero. Proteins measured on 36 of 44 cells read NaN on the other 8.

The modality is part of the data mapping, so a saved session restores it, and
a session of a multimodal dataset is refused on a file whose modalities
differ. A dataset of one feature space has no modality in its mapping: every
session saved before modalities existed still applies.

---

## The same select for 10x feature types

A 10x Cell Ranger matrix of several feature types (Gene Expression with
Antibody Capture) shows one modality per type, and a Cell Ranger 2 file of
several genomes one per genome ({doc}`18_open_loom_and_10x`). Each is a range
of the one matrix, so switching between them reads nothing twice.

---

## Where to go next

- {doc}`13_data_mapping`: every group of the data mapping.
- {doc}`18_open_loom_and_10x`: Loom files and Cell Ranger output.
- {doc}`07_file_format_expectations`: every format Cellucid reads.
