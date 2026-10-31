# Subsampling

```{eval-rst}
.. currentmodule:: cellucid
```

{func}`~cellucid.subsample_cells` gives the cells a Cellucid subsample keeps,
as rows of `obs`: the cells the web app opens for the same request (Open data,
*Cells to open*), and the cells R's `cellucid_subsample_cells()` gives. Use it
to take in a notebook exactly the cells the viewer showed, or to choose cells
the same way before opening.

```python
import cellucid

kept = cellucid.subsample_cells(adata, fraction=0.1, seed=0)                     # uniform
kept = cellucid.subsample_cells(adata, 50_000, method="stratified", groupby="cell_type")
kept = cellucid.subsample_cells(adata, fraction=0.05, method="spatial", basis="X_umap")
kept = cellucid.subsample_cells(adata.n_obs, fraction=0.1, method="blocks", block_size=2048)
small = adata[kept].copy()
```

- The same data, method, size and seed give the same cells, in all three
  implementations: they follow one specification
  ({doc}`../../../web_app/p_developer_docs/19_subsampling_specification`) and
  one shared conformance fixture.
- `uniform` and `blocks` are nested: with one seed, a smaller subsample is
  inside a larger one.
- Any count from 1 is a selection. The web app opens a subsample of at least
  2 cells, since one cell has no coordinate variation to draw
  ({ref}`fewest-cells-open`).
- `blocks` needs `block_size`: the web app derives it from the file's chunk
  layout and shows it with the subsample (Session › Dataset › Details).
- A size of every cell or more returns every cell.
- The rows of `data` are taken as the cells' original positions, so
  `adata[kept]` holds the cells the viewer opens. An object read from a
  Cellucid-ready file in priority order (`anndata.read_zarr` of `prepare()`'s
  output, any sample dataset) holds them in its stored order: it is refused
  with a `ValueError` until {func}`~cellucid.restore_original_order` puts the
  cells back, as `viewer.pull()`, `viewer.subset()` and applying a session
  refuse it. (Its first `k` rows are the `uniform` subsample of `k` cells
  under the file's seed, `uns['cellucid']['priority']['seed']`.)

```python
adata = cellucid.restore_original_order(anndata.read_zarr("pbmc.zarr"))
kept = cellucid.subsample_cells(adata, fraction=0.1, seed=0)
```

## API reference

```{eval-rst}
.. autofunction:: subsample_cells
```

## See also

- {doc}`../../../web_app/b_data_loading/21_subsampled_loading`: opening a
  share of the cells in the web app.
- R: `cellucid_subsample_cells()`, the same cells.
