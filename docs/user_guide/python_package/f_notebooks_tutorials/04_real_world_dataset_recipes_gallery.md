# Real-world dataset recipes (gallery)

Real-world notebooks showing how to prepare various single-cell datasets for visualization with Cellucid.

```{important}
These notebooks are **templates**. They are intentionally explicit about paths, export directories, and validation steps, but you should expect to adapt:
- dataset file locations
- embedding keys (e.g. `X_umap` vs project-specific names)
- gene ID conventions (symbols vs Ensembl)

If you want the conceptual “how exports work” guide first, read:
- {doc}`21_prepare_exports_with_quantization_and_compression`
```

---

## Dataset Gallery

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Apidip
:link: prepare_apidip
:link-type: doc

3D UMAP export workflow for the Apidip dataset.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Braun
:link: prepare_braun
:link-type: doc

Multi-dimensional UMAP export (1D/2D/3D) for the Braun dataset.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Garcia
:link: prepare_garcia
:link-type: doc

Multi-dimensional UMAP export for the Garcia dataset.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` He
:link: prepare_he
:link-type: doc

Multi-dimensional UMAP export for the He dataset.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Kanemaru
:link: prepare_kanemaru
:link-type: doc

Multi-dimensional UMAP export for the Kanemaru dataset.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Miller
:link: prepare_miller
:link-type: doc

Multi-dimensional UMAP export for the Miller dataset.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Suo
:link: prepare_suo
:link-type: doc

Multi-dimensional UMAP export for the Suo dataset.
:::

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` CellRank (_test)
:link: prepare_cellrank_pseudotime_test
:link-type: doc

Real dataset from CellRank with a pseudotime-kernel drift vector field exported to an `exports/_test` folder.
:::

::::

---

## What You'll Learn

Each example notebook demonstrates:

- **Configuration** - Setting up dataset paths and export directories
- **Loading** - Reading AnnData from h5ad or zarr format
- **Embedding** - Working with 1D, 2D, and 3D UMAP coordinates
- **Export** - Preparing optimized binary files for the web viewer
- **Validation** - Sanity-checking output file sizes and structure

---

## Running the Notebooks

```bash
# Clone the repository
git clone https://github.com/theislab/cellucid-python
cd cellucid-python

# Install with examples dependencies
pip install -e ".[dev]"

# Open the notebooks
jupyter lab docs/user_guide/python_package/f_notebooks_tutorials/
```

```{note}
If a notebook refers to a local path that does not exist on your machine, treat it as a configuration placeholder. The notebooks keep all path configuration in one early cell so you can update it safely.
```

---

```{toctree}
:maxdepth: 1
:hidden:

prepare_apidip
prepare_braun
prepare_garcia
prepare_he
prepare_kanemaru
prepare_miller
prepare_suo
prepare_cellrank_pseudotime_test
```
