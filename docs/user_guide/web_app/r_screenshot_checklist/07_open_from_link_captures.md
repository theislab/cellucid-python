# The Session panel, Open…, and the Link or server bar: captures to take

Working material, not published. First prepared 2026-10-22 together with the
scenarios in `docs/_tooling/screenshots/scenarios.mjs` (`buildDataLinkScenarios`
and `buildDataLoadingRefusalScenarios`), and rewritten on 2026-10-23 for the
redesigned Session panel: the **Dataset** box replaced the two-column dataset
table, **Open…** (with **File…** and **Folder…**) replaced the **Files**
buttons (**H5AD**, **Zarr ZIP**, **Zarr folder**), one **Link or server** bar
replaced the **Server** and **Open from link** fields with their **Connect**
buttons, and the **Transfers** list is gone (how much has been read is a line
in the Dataset box). Every scenario below was emitted on 2026-10-25 (see the
status note in {doc}`index`), except the two refusals that would replace
`fail-github-catalog-not-found.png`, which keeps its own scenario. Read
{doc}`00_capture_tooling_and_conventions` first.

Every link in these scenarios points at a loopback host of the capture run (the
docs catalog host serving `cellucid-datasets/anndata/`, or the web repository's permissive-CORS mirror), never the
internet, so the images are reproducible offline. The one server scenario
(`loading-address-server-connected`) starts `cellucid serve` on the Pancreas
sample's store on a free loopback port: it needs `cellucid` on the `PATH` (or
`CELLUCID_CLI` naming it).

## Run

```sh
cd cellucid-python/docs/_tooling/screenshots
node capture.mjs shoot loading-session-panel loading-session-panel-dark \
  loading-dataset-box-details loading-open-menu loading-link-control-typed \
  loading-link-catalog-connected loading-address-server-connected \
  loading-best-guess-notification loading-best-guess-show \
  loading-data-mapping-other-arrays filtering-outlier-no-latent \
  filtering-outlier-slider-95 loading-link-refusal-unrecognised \
  loading-link-refusal-not-found --emit
```

`shoot` prints each `{figure}` block with its `:width:`. Paste it where the
table below says, with the caption given here, **in the same change** — an
orphan PNG fails `python -m pytest` (the `referenced == screenshots` contract).
Then delete the images this set replaces, with their `captures.json`
records, and run `python -m pytest tests/test_docs_current_api_contract.py` and
`sphinx-build -W`.

## The images

Files are under `_static/screenshots/data_loading/` unless the scenario's topic
says otherwise. A file **re-captured in place** keeps its name, so the pages
that show it need no edit: their `:alt:` text and captions already describe the
new panel.

| Scenario | File | Page and place | Caption |
| --- | --- | --- | --- |
| `loading-session-panel` | `data-loading-session-panel.png` (**re-captured in place**) | the eight pages that already show it; nothing to paste | The figure blocks stay; their alt text now names the **Dataset** box, **Open…**, the **Link or server** bar, **Session file** and the collapsed **Data mapping**. |
| `loading-session-panel-dark` | `data-loading-session-panel-dark.png` (new) | `09_screenshots.md`, **Loading and session controls**, after the light panel | The same Session panel in the dark theme. |
| `loading-dataset-box-details` | `dataset-box-details.png` (new) | `01_loading_options_overview.md`, end of **The Dataset box**; and `14_open_from_link.md`, **How much has been read** | The **Dataset** box for an `.h5ad` read in place from a web address: the name, cells and genes, format and host, "… read in … requests", and **Details** open with the file, its format and size, how it is read, and what the dataset holds. |
| `loading-open-menu` | `choose-local-data-control.png` (**re-captured in place**) | `03_browser_file_picker_tutorial.md` and `09_screenshots.md`; nothing to paste | **Open…** under **File or folder**, with its menu open: **File…** (AnnData, Seurat/SCE `.rds`, Loom, 10x `.h5`, `.h5mu`), **Folder…** (Zarr, 10x outs, MatrixMarket). |
| `loading-link-control-typed` | `open-from-link-control.png` (new) | `14_open_from_link.md`, end of **The address bar** | The **Link or server** bar with a Zenodo record typed in. **Open** (or Enter) opens it; the page address then carries `?data=zenodo:23288289`. |
| `loading-link-catalog-connected` | `connect-data-link.png` (new; **replaces** `connect-github-catalog.png`) | `14_open_from_link.md`, after **The address bar**; and in place of `connect-github-catalog.png` in `02_local_demo_tutorial.md`, `09_screenshots.md` and `11_custom_dataset_repository.md` | A dataset catalog folder opened from a link: the dropdown lists its datasets, labelled "Dataset from <host>", and the line under the bar reads **Link**, the address and **Disconnect**, which brings the samples back. |
| `loading-address-server-connected` | `connect-remote-server.png` (**re-captured in place**) | `09_screenshots.md`; nothing to paste | The **Link or server** bar connected to a running `cellucid serve`: the line under it reads **Server**, the address and **Disconnect**. |
| `loading-best-guess-notification` | `best-guess-notification.png` (existing scenario, re-capture) | `13_data_mapping.md`, **Files without a UMAP: the best guess** | The notification a file without a UMAP raises: what was guessed, and **Show**. Data mapping does not open on its own. |
| `loading-best-guess-show` | `best-guess-show-checklist.png` (new) | `13_data_mapping.md`, after the notification | **Show** opens **Data mapping › Embeddings** on the guessed row. |
| `loading-data-mapping-other-arrays` | `data-mapping-other-arrays.png` (new) | `13_data_mapping.md`, **Choose the embeddings views can show** | The Embeddings checklist as one flat list, with **Other arrays** open at its end: an array that is off and wider than three columns shows its width. |
| `filtering-outlier-no-latent` | `filtering/outlier-no-latent.png` (new) | `e_filtering/02_outlier_filtering_per_active_field.md`, **Turning it on: choose a latent space** | Without a latent space the outlier filter is its label, its (i) and one **Choose latent space** button; the button's title names the file's candidates. |
| `filtering-outlier-slider-95` | `filtering/outlier-slider-95.png` (**re-captured in place**) | the four pages that show it; nothing to paste | With a latent space (stored quantiles): the slider alone. The state with a slider and one line (*“First use reads X_pca (… MB).”*) is `outlier-cost-line` in {doc}`02_pages_interactions_fields_filtering`; it has no scenario yet. |
| `loading-link-refusal-unrecognised` | `fail-data-link-unrecognised.png` (new) | `14_open_from_link.md`, section **Why a link can be refused**, before the table | Text Cellucid does not recognise is refused before anything is requested, with the kinds of link it opens. |
| `loading-link-refusal-not-found` | `fail-data-link-not-found.png` (new; **replaces** `fail-github-catalog-not-found.png`) | `14_open_from_link.md`, same section, after the table; and in place of `fail-github-catalog-not-found.png` in `02_local_demo_tutorial.md`, `08_troubleshooting_data_loading.md`, `09_screenshots.md` and `11_custom_dataset_repository.md` | A web address that is not a dataset catalog folder: the notification says it holds no `datasets.json` and what a folder link must point at. The dataset that was open stays open. |

## Captures without a scenario (by hand)

| File | How to reach it | Page and place | Caption |
| --- | --- | --- | --- |
| `data_loading/fail-safari-hidden-files.png` (new) | Safari (WebKit's own folder dialog; Playwright WebKit does not reproduce it): **Open… › Folder…**, choose a Zarr **v2** folder (`cellucid-datasets/anndata/pancreas.zarr`); crop the notification "This browser left out the folder's hidden files … Drag the folder onto the page, or open its .zarr.zip." | `03_browser_file_picker_tutorial.md`, **Hidden metadata files** | Safari's folder dialog leaves out a Zarr v2 store's hidden metadata files; drag the folder onto the page or open its `.zarr.zip` instead. |
| `server/open-served-dataset.png` (**re-capture in place**) | `cellucid serve` on the Pancreas store, open its printed Viewer URL, whole window | `04_server_tutorial.md`; nothing to paste | Its alt text now describes the **Dataset** box and the **Server** line under the **Link or server** bar. |
| `web_app/app-overview-cell-type.png` (re-capture in place, optional) | the Suo sample, sidebar scrolled to **Session**, whole window | `a_orientation/03_quick_tour_60_seconds.md` and the six other pages that show it | Its alt text still describes the earlier panel ("a dataset summary … above sample and local-data controls"); update it with the recapture. |

## Captures for the formats tutorials (2026-10-24, by hand)

The tutorials of the R and formats work have no image yet; each needs one, with
the pointer on the control named (the screenshot rules of
{doc}`00_capture_tooling_and_conventions`). The datasets are public and small.

| File | How to reach it | Page and place | Caption |
| --- | --- | --- | --- |
| `data_loading/open-geo-mtx-sample.png` (new) | `?data=geo:GSE336826/GSM9843410_12hpf_scRNA_matrix.mtx.gz`; once open (about 170 MB), the whole window coloured by `Cluster` | `20_tutorials_open_published_data.md`, end of section 1 | A GEO sample's MatrixMarket files with its UMAP table: 4,278 cells, read once. |
| `data_loading/open-zenodo-seurat-rds.png` (new) | `?data=zenodo:8071138/Mice_Mast_cells.Rds`; **Data mapping › Genes** open, pointer on the **Modality** select listing `integrated`, `RNA`, `decontX` | `20_tutorials_open_published_data.md`, step 4 of section 2 | A Seurat object from Zenodo: its assays are modalities. |
| `data_loading/open-cellranger-outs-folder.png` (new) | **Open… › Folder…** on the 1k PBMC `outs/` the tutorial assembles; the Dataset box and the t-SNE coloured by `graphclust` | `20_tutorials_open_published_data.md`, end of section 3, and `18_open_loom_and_10x.md` | Cell Ranger's own t-SNE and clusters, read beside the matrix. |
| `data_loading/modality-select-citeseq.png` (new) | `?data=gh:theislab/cellucid-demo-custom-datasets/catalog&dataset=synthetic-citeseq` (after the demo repository is pushed); **Data mapping › Genes**, pointer on **Modality** (`rna` · 8 features, `prot` · 4 features) | `19_multimodal_mudata.md`, **The Modality select**, and the Python and R CITE-seq tutorials | The Modality select of a MuData. |
| `r_package/rstudio-viewer-pbmc-small.png` (new) | RStudio: `cellucid_view(pbmc_small)` in the console, the Viewer pane showing the 80 cells; pointer on the console line | `r_package/d_viewing_loading/04_view_in_rstudio_positron.md`, **The one line**, and `r_package/f_notebooks_tutorials/04_explore_a_seurat_object_from_r.md`, step 1 | `pbmc_small` in RStudio's Viewer pane, served from memory. |
| `r_package/positron-viewer-pbmc-small.png` (new) | the same in Positron | `r_package/d_viewing_loading/04_view_in_rstudio_positron.md` | The same in Positron. |
| `r_package/rstudio-selection-back.png` (new) | after a confirmed lasso, the console showing `20 cells selected with the lasso` and `table(picked$groups)` beside the Viewer pane | `r_package/f_notebooks_tutorials/04_explore_a_seurat_object_from_r.md`, step 3 | A lasso in the viewer, the cells in R. |
| `notebooks/citeseq-show-anndata.png` (new) | Jupyter: `cellucid.show_anndata("pbmc_1k_citeseq.h5mu")` from the Python CITE-seq tutorial, the embedded viewer coloured by `CD19_TotalSeqB` with **Modality** on `prot` | `python_package/f_notebooks_tutorials/24_citeseq_mudata_prepare_and_show.md`, step 4 | Proteins of a CITE-seq MuData in a notebook. |

## Removed

- `loading-link-file-transfers` (`transfers-file-in-place.png`): the
  **Transfers** list no longer exists. Its content is the read line of the
  Dataset box, captured by `loading-dataset-box-details`.

## What changed in the tool for this set

- `capture.mjs` passes `urls.sampleOrigin` (the CORS mirror that answers range
  requests) to scenarios, and a scenario's `query` may be a function of `urls`
  (a `?data=` link to a loopback host).
- A startup that names its data exactly (`?data=`) keeps the welcome overlay
  hidden; the tool now waits for the overlay only until startup shows it or is
  ready without it (`keepWelcome` still throws if it never opens).
- `PANCREAS_READY` and the other count waits read "Showing all 3,696 **cells**";
  the app has said "cells" since the panel was renamed, so every Pancreas
  scenario failed at its first step before this fix.
