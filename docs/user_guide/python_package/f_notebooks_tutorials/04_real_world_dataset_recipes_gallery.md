# Real-world dataset recipe

The preparation notebook in this section uses a real, checked-in Pancreas
AnnData file. It runs from a fresh clone: there are no lab-server paths, hidden
inputs, or values that you must ask someone to supply.

::::{grid} 1
:gutter: 3

:::{grid-item-card} {octicon}`beaker;1.5em;sd-mr-1` Prepare Pancreas
:link: prepare_pancreas
:link-type: doc

Turn 2,531 real pancreatic endocrinogenesis cells into a compact prepared
dataset, validate the manifests, and open the result in Cellucid.
:::

::::

## When this example is useful

Run the notebook when you want to:

- see the exact boundary between an `AnnData` object and a prepared dataset,
- learn which embedding, observation, expression, and connectivity inputs map
  to `prepare(...)`,
- make a small network-efficient export before scaling the same pattern to your
  full study, or
- verify a Python environment with a known public input before debugging your
  own data.

The notebook deliberately exports 128 highly variable genes so it finishes
quickly. Your production export can select a larger gene panel or all genes.
For the reasoning behind each size and precision choice, continue with
{doc}`21_prepare_exports_with_quantization_and_compression`.

## Run the notebook

```bash
git clone https://github.com/theislab/cellucid-python
cd cellucid-python
python -m pip install -e . jupyterlab
jupyter lab docs/user_guide/python_package/f_notebooks_tutorials/prepare_pancreas.ipynb
```

```{note}
The notebook writes its prepared dataset to a temporary directory and prints
the exact location. Its final cell stops the viewer and removes that directory.
Copy the export elsewhere before cleanup if you want to keep it.
```

```{toctree}
:maxdepth: 1
:hidden:

prepare_pancreas
```
