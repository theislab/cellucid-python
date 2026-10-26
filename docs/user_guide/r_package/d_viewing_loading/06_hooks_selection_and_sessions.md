# Hooks, Selections and Sessions

**Audience:** analysts who want the viewer and R to work together
**What you'll do:** send commands to an open viewer, react to what you do in
it, and bring a selection or a whole session back to R.

Every viewer has a channel (a WebSocket) to the R session that opened it.
Commands sent before the page has connected wait and go out in order as soon as
it does.

## Commands

```r
viewer <- cellucid_view(seu)

cellucid_color_by(viewer, "seurat_clusters")
cellucid_highlight(viewer, WhichCells(seu, idents = "B"), color = "#377eb8")
cellucid_set_visibility(viewer, seu$nCount_RNA < 500, visible = FALSE)
cellucid_clear_highlights(viewer)
cellucid_reset_camera(viewer)
cellucid_set_data_mapping(viewer, modality = "ADT")   # a Seurat object with an ADT assay
```

Cells are given by name, by 1-based position, or as `TRUE`/`FALSE` per cell.
`cellucid_set_data_mapping()` waits for the viewer to apply the mapping and
stops with the viewer's reason when it refuses it. An object is served with one
matrix per assay (its `data` layer, else `counts`) and no layers, so the
mapping's `expression` stays `"X"`; change the modality, the graph (`knn`) or
the embeddings instead. A file given by path may hold layers
(`expression = "layer:counts"`).

```{figure} ../../../_static/screenshots/r_package/cellucid-view-driven-from-r.png
:alt: Cellucid opened by cellucid_view() on a Seurat object in an R session, after R called cellucid_color_by(viewer, "clusters") and cellucid_highlight() on the Ductal cells: the points are coloured by clusters with the Ductal cells highlighted, and the Dataset box names the object.
:width: 1440px

`cellucid_view()` on a Seurat object of the Pancreas sample, opened in a
browser tab, after R ran `cellucid_color_by(viewer, "clusters")` and
`cellucid_highlight()` on the Ductal cells. The Dataset box reads *Zarr v2
folder · R*: the object is served from memory by this R session.
```

## Events

```r
cellucid_on(viewer, "selection", function(event) {
  message(length(event$cells), " cells selected with the ", event$source)
})
cellucid_on(viewer, "click", function(event) str(event))
```

The events are `selection`, `hover`, `click`, `ready`, `data_mapping`,
`command_error`, `pong`, `debug_snapshot` and `session_bundle`. A handler runs
whenever R is idle; one that fails gives a warning and the others still run.
`cellucid_off(viewer, "hover")` removes a kind's handlers.

Cells are 1-based positions in every event, as the commands take them: a
`selection` gives `cells` and `cell_names`, a `hover` or `click` gives `cell`
and `cell_name`. (Python's events are 0-based, as Python indexes.) `pong` and
`debug_snapshot` answer `cellucid_debug_connection(viewer)`, the report to read
first when a viewer stays blank or does not react; the channel's own
keep-alive pings, every 20 seconds, raise no event.

## A selection back in R

```r
selected <- cellucid_selection(viewer)
selected$cells        # 1-based positions
selected$cell_names   # the object's cell names
sub <- subset(seu, cells = selected$cell_names)
```

A selection of any size arrives at once (millions of cells decode in a fraction
of a second).

## The session

```r
path <- cellucid_session(viewer, file = "pbmc.cellucid-session")
```

`cellucid_session()` asks the viewer for its whole session — highlights,
colours, filters, annotations, the camera — and saves the `.cellucid-session`
file the viewer's **Save session** writes. Open it again in the web app with
the same dataset, under the same dataset id. It stops with what to do when the
viewer is not open or does not answer within `timeout` seconds.

## Writing a session back onto the object

```r
seu <- cellucid_apply_session(seu, viewer)
# or, from a saved file:
seu <- cellucid_apply_session(seu, "pbmc.cellucid-session", expected_dataset_id = "seu")
```

`cellucid_apply_session()` adds what was made in the viewer to the object's
cell metadata (Seurat `meta.data`, `colData` of a SingleCellExperiment or
SummarizedExperiment): one `TRUE`/`FALSE` column per highlight group
(`cellucid_highlight__highlight_1`, …), and each user-defined field — a
categorical field as a factor, a copied column or gene as numbers. The session
itself is recorded in `misc$cellucid$session` (Seurat) or
`metadata(x)$cellucid$session`.

The rules are those of the Python package's
`apply_cellucid_session_to_anndata()`, check for check, and every check runs
before anything is added: the dataset id (a viewer's own, or
`expected_dataset_id` for a file), the cell count, the data mapping the
session was saved with, the columns and genes it uses, and every cell index.
A column the object already has is never overwritten, and a highlight saved
relative to a category is refused when that category no longer holds the same
cells. If `cellucid_view()` was given arguments that change what is shown
(`embeddings`, `obs`, `modalities`, …), pass the same ones after the session.
