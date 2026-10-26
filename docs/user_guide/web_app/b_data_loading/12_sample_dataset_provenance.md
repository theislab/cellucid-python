# Sample dataset provenance and reproducibility

Cellucid's built-in sample picker offers five datasets, published in
[`cellucid-datasets`](https://github.com/theislab/cellucid-datasets) as
Cellucid-ready AnnData Zarr stores (`anndata/<id>.zarr/`) listed by one dataset
catalog, `anndata/datasets.json`, which GitHub Pages serves at
`https://theislab.github.io/cellucid-datasets/anndata/`. Any AnnData reader
opens them (`anndata.read_zarr`). Each is built from
a published study, and each names its upstream study, that study's URL, and
the citation to use in its own `uns['cellucid']` (`source` and `citation`),
which the **Dataset** box shows (**Details**: **Source**, **Citation**).

All five carry the same reproducibility guarantee. Each is rebuilt by a pinned
build script that verifies its source file, its environment, the writer's own
code, and everything it writes, and that stops rather than publishes when a
rebuild differs anywhere. What differs between them is where the source comes
from: Pancreas downloads a commit-pinned public file, while the four human
atlases are built from integration outputs the dataset repository does not
hold and may not redistribute. This page says exactly what is pinned and what
that does and does not let you do.

## The catalog

| Sample | Catalog id | Cells | Genes | Upstream study |
| --- | --- | --- | --- | --- |
| Developing human immune system (Suo) | `suo` | 561,947 | 5,103 | [E-MTAB-11343](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11343) |
| Developing human gonad (Garcia) | `garcia` | 219,731 | 5,754 | [E-MTAB-10551](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10551) |
| Developing human lung (He) | `he` | 71,650 | 5,152 | [E-MTAB-11278](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11278) |
| Developing human heart (Kanemaru) | `kanemaru` | 131,636 | 3,691 | [Developing Heart Cell Atlas](https://www.heartcellatlas.org/foetal.html) ([doi:10.1101/2024.04.29.591736](https://doi.org/10.1101/2024.04.29.591736)) |
| Pancreatic endocrinogenesis (scVelo) | `pancreas` | 3,696 | 3,753 | [GSE132188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132188) |

`suo` is the catalog default. The counts are the cell and gene counts each
file holds, as the dataset list and the **Dataset** box show them.

The ArrayExpress accessions are the ones holding each study's single-cell
RNA-seq libraries. For Suo that is E-MTAB-11343; E-MTAB-11341 is the same
study's Visium spatial accession, not the source of this sample. Kanemaru is the
fetal heart of the Developing Heart Cell Atlas (Cranley, Kanemaru et al.,
bioRxiv 2024), whose RNA data have no public archive accession, so it is
identified by its preprint DOI and its data portal; the adult-heart study
E-MTAB-12916 shares none of its donors.

Each file holds its embeddings in `obsm`. The four human atlases each carry a
1D, 2D, and 3D UMAP computed from one shared neighbour graph (`X_umap_1d`,
`X_umap_2d`, `X_umap_3d`), the first two components of the unintegrated PCA
(`X_pca`), and the integrated latent space the graph was built on (`scvi`, or
`scanvi` for Suo), which is the dataset's `latentKey`; Garcia, He, and Kanemaru
also carry the 2D UMAP shipped in their source files (`X_umap`). For each of the
four, the launch defaults are `X_umap_1d`, `X_umap_2d`, and `X_umap_3d`, and
the starting view opens on the 3D UMAP. Pancreas has its own set; see
{doc}`10_standard_pancreas_dataset`.

The gene count above is the whole gene axis the viewer can search, and in every
case it is smaller than the upstream study's. The four human samples publish only
the genes whose Ensembl accession resolved to a symbol, out of files that had
already been reduced to 8,192 features upstream; `pancreas` publishes the highly
variable, nonconstant genes of its 27,998-gene source. Neither reduction is
recoverable in the viewer, so if you are looking for a specific gene, read
{doc}`../d_fields_coloring_legends/07_genes_in_the_built_in_samples` before
concluding it is absent from the study.

## How each sample is built

Two build scripts in `cellucid-datasets/scripts/` produce the five
stores, and both delegate everything that does not depend on the dataset
to one shared module, `pinned_build.py`:

- `build_pancreas.py` builds `pancreas`. It downloads the commit-pinned scVelo
  source file, or reads a local copy passed with `--source`.
- `build_atlas.py --dataset <id>` builds `suo`, `garcia`, `he`, and `kanemaru`
  with one recipe: it reads the integration output, applies the published cell
  subset (Suo only), curates `obs`, computes one neighbour graph on the
  integrated latent space with its 1D, 2D, and 3D UMAPs, single-threaded and
  seeded, names the genes, and writes the Cellucid-ready AnnData with
  `cellucid.prepare(format="zarr")`. It reads its inputs from
  `--source-dir`, which defaults to `$CELLUCID_ANNDATA_DIR` and otherwise to a
  `.anndatas/` directory beside the repository.

Each sample has one build record, `sources/<id>.json`, which is the build's
contract rather than a description of it. The environment every sample is built
in is recorded once, in `sources/build-environment.json`, together with
`requirements-build.lock`.

## What the pins cover

- **The environment** — the interpreter implementation and version, the platform
  and machine, the zlib runtime version, the OpenBLAS build with a single
  thread, and the fixed hash seed and random seed, all from
  `sources/build-environment.json`. `requirements-build.lock` is itself pinned
  by SHA-256, every package it lists must be installed at exactly that version,
  and the installed distributions must be exactly that set plus `cellucid`.
- **The producer** — a SHA-256 over each `cellucid` source file
  `prepare()` executes, combined into one manifest digest, plus a check at build time
  that the writing path reached no module the list omits, so the pin cannot
  quietly stop covering the writer. The installed `cellucid` version must equal
  the recorded producer version.
- **The input** — the source file's byte size and SHA-256, verified against the
  actual file before it is read, and digests of every array the build reads from
  it. The four atlases also pin, by size and SHA-256, the upstream gene
  identifier list and the Ensembl-to-symbol table their gene names come from.
- **The analysis** — the parameters and digests of everything the build
  derives: the neighbour graph and the embeddings, the gene-naming census of the
  four atlases, and the velocity projections of Pancreas.
- **The output** — the `uns['cellucid']` the store carries, the published
  expression matrix (shape, stored values and digest), and the file count,
  total byte size, and a SHA-256 store digest over every file of the store.

The environment and producer checks run before the build reads a source,
touches the network, or writes output. The new store is written into a
staging location beside the published one, read back with anndata and
compared element by element with what was written and with everything
`sources/<id>.json` pins, and only then swapped into place. A rebuild that
differed anywhere stops instead of publishing, and a second build reproduces
the store byte for byte. The dataset repository's contract suites re-measure
every published store against the same store digest and open each one with
anndata, so a store that changes under the repository fails there too.

The store digest covers the dataset's own files, not its starting view
({doc}`../l_sessions_sharing/04_official_sample_states`): the session file's
bytes are pinned by the catalog entry's own `session.sha256` and
`session.bytes` instead, so the view can be republished without moving a
store pin.

## What the pins do not claim

- **That the upstream file is immutable.** The Pancreas build downloads its URL
  unless it is handed a local copy. It is the recorded digest, not the URL, that
  turns a changed or unavailable upstream into a failure rather than a silent
  difference.
- **That anyone can rebuild the four atlases.** Their sources are outputs of an
  scVI/scANVI integration pipeline — derivatives of the upstream studies above,
  not archive downloads — and are not redistributable. The repository
  never holds them; `sources/<id>.json` records exactly which bytes rebuild each
  one, and the build refuses any others. Without those exact files you can verify
  a published store but not regenerate it.
- **That the build runs anywhere.** It refuses to start on a different
  interpreter version, platform, BLAS build, or installed package set, rather
  than publishing something different under the same name.

## What to do with this

- **Citing a sample.** Cite the upstream study, not Cellucid and not the
  dataset repository. Each file carries the study name and URL under
  `uns['cellucid']['source']` and the citation under
  `uns['cellucid']['citation']`. To cite Cellucid
  itself, see {doc}`/citing`.
- **Confirming you hold the published store.** Compare the store you have
  against the output contract in `sources/<id>.json`
  (`anndata_contract.recorded`): its store digest, file count, and total byte
  size. These values change
  whenever a dataset is regenerated, so read them from the file rather than
  copying a digest out of documentation.
- **Rebuilding.** `pancreas` can be rebuilt from the dataset repository alone;
  the four atlases also need their pinned source files. Every rebuild needs the
  environment `sources/build-environment.json` pins.
- **Needing a pinned build of your own.** Write a Cellucid-ready AnnData from
  the upstream study yourself with {doc}`../../python_package/c_data_preparation_api/index`, and
  record the source digest and the environment beside the result as described
  in
  {doc}`../../python_package/b_concepts_mental_models/04_dataset_identity_and_reproducibility`.

## Related pages

- {doc}`10_standard_pancreas_dataset` — what the Pancreas sample contains and
  how to load it.
- {doc}`../d_fields_coloring_legends/07_genes_in_the_built_in_samples` — which
  genes each sample publishes and why the rest are not there.
- {doc}`../l_sessions_sharing/04_official_sample_states` — the default
  session the sample catalog names for each sample.
- {doc}`06_dataset_identity_why_it_matters` — what makes two loads the same
  dataset for sessions, sharing, and community annotation.
- {doc}`11_custom_dataset_repository` — publishing and validating a catalog of
  your own.
