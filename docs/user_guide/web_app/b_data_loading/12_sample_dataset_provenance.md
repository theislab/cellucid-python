# Sample dataset provenance and reproducibility

Cellucid's built-in sample picker offers five prepared datasets, published in
[`cellucid-datasets`](https://github.com/theislab/cellucid-datasets). Each is a
real export of a published study, and each names its upstream study, that
study's URL, and the citation to use in its own
`exports-v3/<id>/dataset_identity.json`, under `source`.

All five carry the same reproducibility guarantee. Each is rebuilt by a pinned
build script that verifies its source file, its environment, the exporter's own
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

`suo` is the catalog default. The counts are the `stats.n_cells` and
`stats.n_genes` values each generation publishes in its own
`dataset_identity.json`.

The ArrayExpress accessions are the ones holding each study's single-cell
RNA-seq libraries. For Suo that is E-MTAB-11343; E-MTAB-11341 is the same
study's Visium spatial accession, not the source of this sample. Kanemaru is the
fetal heart of the Developing Heart Cell Atlas (Cranley, Kanemaru et al.,
bioRxiv 2024), whose RNA data have no public archive accession, so it is
identified by its preprint DOI and its data portal; the adult-heart study
E-MTAB-12916 shares none of its donors.

Each generation also declares its embeddings there, under `embeddings`. The
four human atlases each carry a 1D, 2D, and 3D UMAP computed from one shared
neighbour graph (`X_umap_1d`, `X_umap_2d`, `X_umap_3d`) and a 2D view of the
unintegrated PCA (`X_pca`, labelled **PCA (unintegrated)**); Garcia, He, and
Kanemaru also carry the 2D UMAP shipped in their source files (`X_umap`,
labelled **UMAP (source file)**). For each of the four, the declared defaults
are `X_umap_1d`, `X_umap_2d`, and `X_umap_3d`, and the default dimension is 3D.
Pancreas declares its own set; see {doc}`10_standard_pancreas_dataset`.

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
generations, and both delegate everything that does not depend on the dataset
to one shared module, `pinned_build.py`:

- `build_pancreas.py` builds `pancreas`. It downloads the commit-pinned scVelo
  source file, or reads a local copy passed with `--source`.
- `build_atlas.py --dataset <id>` builds `suo`, `garcia`, `he`, and `kanemaru`
  with one recipe: it reads the integration output, applies the published cell
  subset (Suo only), curates `obs`, computes one neighbour graph on the
  integrated latent space with its 1D, 2D, and 3D UMAPs, single-threaded and
  seeded, names the genes, and exports. Its parameters are the ones the
  preparation notebooks in `cellucid-datasets/notebooks/` recorded; the script
  is their deterministic transcription. It reads its inputs from
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
- **The producer** — a SHA-256 over each of the 23 `cellucid` source files the
  export path executes, combined into one manifest digest, plus a check at build
  time that the export path reached no module the list omits, so the pin cannot
  quietly stop covering the exporter. The installed `cellucid` version must
  equal the recorded producer version, which must equal the
  `cellucid_data_version` the generation declares.
- **The input** — the source file's byte size and SHA-256, verified against the
  actual file before it is read, and digests of every array the build reads from
  it. The four atlases also pin, by size and SHA-256, the upstream gene
  identifier list and the Ensembl-to-symbol table their gene names come from.
- **The analysis** — the parameters and digests of everything the build
  derives: the neighbour graph and the embeddings, the gene-naming census of the
  four atlases, and the velocity projections of Pancreas.
- **The output** — the identity, the manifests, a digest of every coordinate and
  vector payload, and the file count, total byte size, and a SHA-256 generation
  digest over the whole exported directory.

The environment and producer checks run before the build reads a source,
touches the network, or writes output. The new generation is written into a
staging directory beside `exports-v3/<id>/`, validated against everything
`sources/<id>.json` pins, and only then swapped into place. A rebuild that
differed anywhere stops instead of publishing. The dataset repository's
contract suite re-measures every published directory against the same
generation digest, so a generation that changes under the repository fails
there too.

The generation digest covers every file in `exports-v3/<id>/` except the two
publication sidecars, `default.cellucid-session` and `state-snapshots.json`.
Those advertise and hold the reviewed starting view
({doc}`../l_sessions_sharing/04_official_sample_states`); the session file's
bytes are pinned by the catalog entry's own `state_sha256` and `state_bytes`
instead, so the view
can be republished without moving a generation pin.

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
  a published generation but not regenerate it.
- **That the build runs anywhere.** It refuses to start on a different
  interpreter version, platform, BLAS build, or installed package set, rather
  than publishing something different under the same name.

## What to do with this

- **Citing a sample.** Cite the upstream study, not Cellucid and not the
  dataset repository. Each `exports-v3/<id>/dataset_identity.json` carries the
  study name, its URL, and the citation text under `source`. To cite Cellucid
  itself, see {doc}`/citing`.
- **Confirming you hold the published generation.** Compare the directory you
  have against `export_contract` in `sources/<id>.json`: `generation_sha256`,
  `expected_file_count`, and `expected_total_bytes`. These values change
  whenever a dataset is regenerated, so read them from the file rather than
  copying a digest out of documentation.
- **Rebuilding.** `pancreas` can be rebuilt from the dataset repository alone;
  the four atlases also need their pinned source files. Every rebuild needs the
  environment `sources/build-environment.json` pins.
- **Needing a pinned build of your own.** Export from the upstream study
  yourself with {doc}`../../python_package/c_data_preparation_api/index`, and
  record the source digest and the environment beside the result as described
  in
  {doc}`../../python_package/b_concepts_mental_models/04_dataset_identity_and_reproducibility`.

## Related pages

- {doc}`10_standard_pancreas_dataset` — what the Pancreas sample contains and
  how to load it.
- {doc}`../d_fields_coloring_legends/07_genes_in_the_built_in_samples` — which
  genes each sample publishes and why the rest are not there.
- {doc}`../l_sessions_sharing/04_official_sample_states` — the reviewed
  starting view each sample publishes beside its scientific files.
- {doc}`06_dataset_identity_why_it_matters` — what makes two loads the same
  dataset for sessions, sharing, and community annotation.
- {doc}`11_custom_dataset_repository` — publishing and validating a catalog of
  your own.
