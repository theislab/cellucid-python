# Sample dataset provenance and reproducibility

Cellucid's built-in sample picker offers five prepared datasets, published in
[`cellucid-datasets`](https://github.com/theislab/cellucid-datasets). Each is a
real export of a published study, and each names its upstream study, that
study's URL, and the citation to use in its own
`exports/<id>/dataset_identity.json`, under `source`.

They do not all carry the same reproducibility guarantee. One is rebuilt by a
pinned script that verifies its source, its environment, the exporter's own
code, and everything it writes. The other four are published artifacts whose
exact bytes are recorded but whose build is not pinned, and which cannot be
rebuilt without source data this project does not distribute. This page says
which is which, so you do not have to infer it.

## The catalog

| Sample | Catalog id | Cells | Genes | Upstream study |
| --- | --- | --- | --- | --- |
| Developing human immune system (Suo) | `suo` | 561,947 | 5,103 | [E-MTAB-11341](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11341) |
| Developing human gonad (Garcia) | `garcia` | 219,731 | 5,754 | [E-MTAB-10551](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10551) |
| Developing human lung (He) | `he` | 71,650 | 5,152 | [E-MTAB-11278](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11278) |
| Developing human heart (Kanemaru) | `kanemaru` | 131,636 | 3,691 | [E-MTAB-12916](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12916) |
| Pancreatic endocrinogenesis (scVelo) | `pancreas` | 3,696 | 3,753 | [GSE132188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132188) |

`suo` is the catalog default. The counts are the `stats.n_cells` and
`stats.n_genes` values each generation publishes in its own
`dataset_identity.json`.

The gene count above is the whole gene axis the viewer can search, and in every
case it is smaller than the upstream study's. The four human samples publish only
the genes whose Ensembl accession resolved to a symbol, out of files that had
already been reduced to 8,192 features upstream; `pancreas` publishes the highly
variable, nonconstant genes of its 27,998-gene source. Neither reduction is
recoverable in the viewer, so if you are looking for a specific gene, read
{doc}`../d_fields_coloring_legends/07_genes_in_the_built_in_samples` before
concluding it is absent from the study.

## Pancreas is pinned

`pancreas` is produced by `cellucid-datasets/scripts/build_pancreas.py`, and
[`sources/pancreas.json`](https://github.com/theislab/cellucid-datasets/blob/main/sources/pancreas.json)
is that script's build contract rather than a description of it. It pins:

- **the source** — repository, commit, URL, SHA-256, and byte size, verified
  against the actual file before the export runs;
- **the environment** — interpreter version and platform, a hashed
  `requirements-pancreas.lock` whose every pin must match what is installed,
  the OpenBLAS build and its single-thread setting, and the fixed seeds;
- **the producer** — a SHA-256 over each of the 20 `cellucid` modules the
  export path executes, plus a check at build time that the export path reached
  no module the list omits, so the pin cannot quietly stop covering the
  exporter;
- **the inputs** — cell and gene counts, category levels, and digests of the
  source expression, embedding, and connectivity matrices; and
- **the output** — the file count, the total byte size, and a SHA-256
  generation digest over the whole exported directory.

The environment and producer checks run before the build touches the network or
the output directory, the source is verified by digest before the export runs,
and the finished directory is measured against the recorded generation digest
before it is published. A rebuild that differed anywhere stops instead of
publishing. The dataset repository's contract suite re-measures the published
directory against that same digest, so a generation that changes under the
repository fails there too.

Two things the pin does not claim. It does not make the upstream file
immutable: the build downloads that URL unless it is handed a local copy of the
file, and it is the recorded digest, not the URL, that turns a changed or
unavailable upstream into a failure rather than a silent difference. And it is
not a claim that the build runs anywhere: it refuses to start on a different
interpreter version, platform, BLAS build, or installed package set, rather
than publishing something different under the same name.

For what this dataset contains and how to load it, see
{doc}`10_standard_pancreas_dataset`.

## The other four are recorded, not pinned

`suo`, `garcia`, `he`, and `kanemaru` are produced by the notebook pairs in
`cellucid-datasets/notebooks/`: one notebook computes the shared neighbour
graph and the 1D, 2D, and 3D embeddings and writes an intermediate `.h5ad`, and
one exports that intermediate. What is known about the four published
generations is recorded in
[`sources/unpinned-generations.json`](https://github.com/theislab/cellucid-datasets/blob/main/sources/unpinned-generations.json),
which is explicit that it describes artifacts rather than pinning a build.

Recorded, and re-checked by the dataset repository's contract suite:

- the SHA-256 generation digest, file count, and total byte size of each
  published directory, measured exactly the way the Pancreas generation is
  measured;
- the producer version each export states about itself in its
  `dataset_identity.json`; and
- the file count printed by the committed output of the export notebook, which
  ties each published directory to one committed notebook run.

Not established:

- that re-running the notebooks would produce these bytes again — nothing there
  is a build contract;
- that the recorded producer version is true; it is the export's own statement
  about itself, with nothing outside the export corroborating it;
- which bytes each export was built from. The preparation notebooks record the
  source file's size, never a checksum, and a matching size is not a matching
  file; and
- that the named notebooks produced these exact directories. The file-count
  agreement is evidence, not identity — two runs differing only in coordinate
  values would agree on it.

Not recoverable:

- a checksum of any source `.h5ad`. None was taken while these generations were
  built, and none can be reconstructed now, because the export format stores no
  cell barcodes;
- a pinned interpreter, library set, BLAS configuration, or thread count. The
  preparation notebooks write the `anndata`, `scanpy`, and `numpy` versions into
  the intermediate `.h5ad`, which is not part of the repository, and nothing
  records `umap-learn`, which governs the embeddings;
- bit-exact embeddings. The notebooks fix `random_state=0`, but scanpy's
  neighbour search is approximate above 4,096 cells and no thread count is
  pinned, so equal inputs need not give equal coordinates; and
- the sources themselves. The notebooks read from a workspace-level
  `.anndatas/` directory that this project does not distribute, so these four
  can be rebuilt only with access to that data.

These four remain complete, published exports of named studies, with fixed
bytes you can check against the recorded digests. What differs is the sense in
which they are reproducible: you can confirm you hold exactly the published
artifact, but you cannot regenerate it from a recipe in the repository.

## What to do with this

- **Citing a sample.** Cite the upstream study, not Cellucid and not the
  dataset repository. Each `exports/<id>/dataset_identity.json` carries the
  study name, its URL, and the citation text under `source`.
- **Confirming you hold the published generation.** Compare the directory you
  have against the recorded generation digest: `sources/pancreas.json` for
  `pancreas`, `sources/unpinned-generations.json` for the other four. Both
  files hold the current values; they change whenever a dataset is regenerated,
  so read them rather than copying a digest out of documentation.
- **Rebuilding.** Only `pancreas` can be rebuilt from the dataset repository
  alone, and only in the environment its build contract pins.
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
