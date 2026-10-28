# Cellucid in your pipeline: notebook round trip (Python and R)

Seven runnable notebooks, each with a twin in R, take Cellucid from "a viewer
I look at" to a step of the analysis: show the data, work in the viewer, and
bring what you made back onto the object in one call — or the other way, push
what you computed into the viewer. Every notebook opens in Google Colab from
its badge, and runs in JupyterLab, Jupyter Notebook 7 and VS Code as it is.

| Topic | Python (AnnData) | R (Seurat, SingleCellExperiment) |
|---|---|---|
| Show a dataset in any notebook, in one line | {doc}`41_show_a_dataset_in_any_notebook` | {doc}`../../r_package/f_notebooks_tutorials/07_view_in_a_notebook` |
| Pull highlights, labels and scores back | {doc}`42_pull_everything_back` | {doc}`../../r_package/f_notebooks_tutorials/08_pull_everything_back` |
| Push clusters, refine them in the viewer, pull them back | {doc}`43_push_clusters_refine_and_pull` (Scanpy Leiden) | {doc}`../../r_package/f_notebooks_tutorials/09_push_clusters_refine_and_pull` (Seurat `FindClusters`) |
| No notebook: the web app's exports, read in one call | {doc}`44_round_trip_without_a_notebook` | {doc}`../../r_package/f_notebooks_tutorials/10_round_trip_without_a_notebook` |
| Gene sets: GMT in, signature scores out | {doc}`45_gene_sets_round_trip` | {doc}`../../r_package/f_notebooks_tutorials/11_gene_sets_round_trip` |
| A selection as an `.h5ad` subset, re-analysed | {doc}`46_export_an_h5ad_subset_and_reanalyse` | {doc}`../../r_package/f_notebooks_tutorials/12_subset_and_reanalyse` |
| Subsampled views: unknown cells, the exact subsample | {doc}`47_subsampled_views` | {doc}`../../r_package/f_notebooks_tutorials/13_subsampled_views` |

The calls, side by side:

| | Python | R |
|---|---|---|
| show | `viewer = cellucid.show_anndata(adata, dataset_name=…, dataset_id=…)` | `viewer <- cellucid_view(seu)` |
| everything back | `viewer.pull(adata)` | `seu <- cellucid_pull(seu, viewer)` |
| from exported files | `cellucid.pull(adata, "exports/")` | `seu <- cellucid_pull(seu, "exports/")` |
| a column into the viewer | `viewer.push(adata.obs["leiden"])` | `cellucid_push(viewer, seu$seurat_clusters, key = "seurat_clusters")` |
| a gene set into the viewer | `viewer.push_gene_set("Beta", genes)` | `cellucid_push_gene_set(viewer, "B", genes)` |
| the viewer's gene sets | `viewer.gene_sets()` | `cellucid_gene_sets(viewer)` |
| GMT files | `cellucid.write_gmt()`, `cellucid.read_gmt()` | `cellucid_write_gmt()`, `cellucid_read_gmt()` |
| cells of a selection or group | `viewer.subset(adata, which)`, `viewer.selection` | `cellucid_subset(seu, viewer, which)`, `cellucid_selection(viewer)` |
| the viewer's subsample | `cellucid.subsample_cells(adata, fraction=…, seed=…)` | `cellucid_subsample_cells(seu, fraction = …, seed = …)` |

## The two-cell pattern

In a notebook, the viewer is a widget that talks to the kernel through the
notebook's own connection, and a notebook hands a widget's messages to the
kernel **between cells**. So:

1. **One cell shows the viewer.** It appears when the cell finishes.
2. **You work in the viewer.** Each highlight, label, score and gene set is
   sent to the kernel as you make it, while the kernel is idle.
3. **The next cell reads it.** `viewer.pull(adata)` (`cellucid_pull(seu,
   viewer)`), `subset`, `gene_sets` and `selection` read what has already
   arrived, at once — nothing is asked and nothing waits.

A call that would wait for the viewer's answer inside the cell that asked
cannot work there: the answer can only arrive after the cell ends. Such calls
(`wait_for_ready()`, `wait_for_event()`; R's `cellucid_wait_ready()`) stop at
once and say what to do instead. Commands never wait: `push()` and
`set_data_mapping()` (R: `cellucid_push()`, `cellucid_set_data_mapping()`) go
as soon as the viewer can take them, so a notebook run top to bottom works. A session file is asked for in one cell and collected in
the next: `request = viewer.request_session_bundle()`, then
`request.result()` (R: `cellucid_request_session()`, then
`cellucid_session_result()`). Commands to the viewer (`highlight_cells`,
`push`, …) go out at once; `push` needs the viewer to have opened its dataset,
so run it in a cell after the viewer shows the cells.

In RStudio, Positron or an R or Python console (a browser-tab viewer:
`cellucid_view()` outside Jupyter, `cellucid.AnnDataViewer`), the same calls
work the same way, and waiting calls are answered while they wait.

## How these notebooks were checked

Every code cell of every notebook here and of its R twin was run, top to
bottom, against the data it uses: the Python ones in an IPython kernel (the
same `show_anndata()` widget a notebook shows), the R ones in an R session and
the viewer cell in an IRkernel. Where a notebook asks you to work in the
viewer, the check played the viewer's part over the viewer's own channel with
what the web app sends; the exported CSV, GMT and `.h5ad` files were written
by the web app's own exporters. The docs build does not execute notebooks.

```{toctree}
:maxdepth: 1
:hidden:

41_show_a_dataset_in_any_notebook
42_pull_everything_back
43_push_clusters_refine_and_pull
44_round_trip_without_a_notebook
45_gene_sets_round_trip
46_export_an_h5ad_subset_and_reanalyse
47_subsampled_views
```
