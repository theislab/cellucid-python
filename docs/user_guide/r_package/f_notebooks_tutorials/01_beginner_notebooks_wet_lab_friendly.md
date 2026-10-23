# Beginner Tutorial (Wet Lab Friendly)

**Audience:** anyone with cells, coordinates and labels
**Time:** 10 minutes
**Goal:** write one file and look at your cells in the Cellucid web app.

You need two things:

- **coordinates for every cell** — a UMAP or t-SNE, one row per cell;
- **labels for every cell** — a table with one row per cell (cluster, sample,
  condition, …).

## Step 1 — Make (or load) the data

Here the data is simulated: 600 cells in three clusters.

```r
library(cellucid)
set.seed(7)

n <- 600
cluster <- sample(c("Alpha", "Beta", "Gamma"), n, replace = TRUE)
centres <- c(Alpha = 0, Beta = 4, Gamma = 8)
umap <- cbind(
  rnorm(n, mean = centres[cluster]),
  rnorm(n, mean = centres[cluster] / 2)
)

cells <- data.frame(
  cluster = factor(cluster),
  condition = sample(c("control", "treated"), n, replace = TRUE),
  quality = runif(n)
)
```

With your own data, `umap` is a two-column matrix and `cells` a data frame, in
the same cell order.

## Step 2 — Write one file

```r
path <- cellucid_prepare(
  out_path = file.path(tempdir(), "my_first_dataset.zarr.zip"),
  format = "zarr.zip",
  embeddings = list(X_umap = umap),
  dataset_id = "my-first-dataset",
  dataset_name = "My first dataset",
  obs = cells,
  force = TRUE
)
path
file.size(path)
```

- `out_path` is where the file goes; it ends in `.zarr.zip` because `format = "zarr.zip"`
  writes one file.
- `embeddings` names the coordinates: `X_umap` is the name UMAPs usually have.
- `dataset_id` is a short id without spaces; `dataset_name` is the title people
  see.
- `force = TRUE` lets you run the same code again and replace the file.

## Step 3 — Open it

1. Go to <https://www.cellucid.com>.
2. In the **Session** panel, under **Files**, click **Zarr ZIP**.
3. Choose the file printed in Step 2.

Colour the cells by `cluster`, `condition` or `quality` from the field list.

## If something goes wrong

- *"obs has … rows, but the embeddings have … cells"* — the table and the
  coordinates describe different numbers of cells.
- *"… already exists; pass force = TRUE to replace it"* — add `force = TRUE`.
- *"format = \"zarr.zip\" writes one .zarr.zip file, and out_path … does not
  end in \".zarr.zip\""* — the name must end in `.zarr.zip`; the message says
  which path to write instead.

Next: {doc}`02_intermediate_notebooks_computational_workflows` adds gene
expression.
