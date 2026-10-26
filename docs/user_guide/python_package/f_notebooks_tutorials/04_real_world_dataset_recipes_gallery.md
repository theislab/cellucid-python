# Real-world dataset recipe

The preparation notebook in this section uses a real, checked-in Pancreas
AnnData file. It runs from a fresh clone: there are no lab-server paths, hidden
inputs, or values that you must ask someone to supply.

::::{grid} 1
:gutter: 3

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Prepare Pancreas
:link: prepare_pancreas
:link-type: doc

Turn 2,531 real pancreatic endocrinogenesis cells into a compact Cellucid-ready
AnnData, read it back with anndata, and open the result in Cellucid.
:::

::::

## When this example is useful

Run the notebook when you want to:

- see the exact boundary between an `AnnData` object and a Cellucid-ready one,
- learn which embedding, observation, expression, and connectivity inputs map
  to `prepare(...)`,
- write a small network-efficient file before scaling the same pattern to your
  full study, or
- verify a Python environment with a known public input before debugging your
  own data.

The notebook deliberately keeps 128 highly variable genes so it finishes
quickly. Your production file can keep a larger gene panel or all genes.
For the reasoning behind each argument, continue with
{doc}`21_prepare_a_cellucid_ready_anndata`.

## Run the notebook

```bash
git clone https://github.com/theislab/cellucid-python
cd cellucid-python
python -m pip install -e . jupyterlab
jupyter lab docs/user_guide/python_package/f_notebooks_tutorials/prepare_pancreas.ipynb
```

```{note}
The notebook writes its `.h5ad` to a temporary directory and prints the exact
location. Its final cell stops the viewer and removes that directory. Copy the
file elsewhere before cleanup if you want to keep it.
```

```{toctree}
:maxdepth: 1
:hidden:

prepare_pancreas
```
