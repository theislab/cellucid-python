# Connectivities (neighbour graph)

**Audience:** computational users
**Time:** 5 minutes
**Goal:** store the neighbour graph the viewer draws as edges

`connectivities` is a cells × cells matrix — usually
`adata.obsp["connectivities"]` from `sc.pp.neighbors` — sparse or dense. It is
stored as `obsp['connectivities']`, a CSR float32 matrix.

```python
cellucid.prepare(
    "data.h5ad",
    format="h5ad",
    embeddings={"X_umap": adata.obsm["X_umap"]},
    dataset_id="pbmc-3k",
    dataset_name="PBMC 3k",
    connectivities=adata.obsp["connectivities"],
)
```

The graph must be:

- square, one row and one column per cell;
- **exactly symmetric** (an undirected graph; scanpy's connectivities are);
- finite, non-negative weights within float32's range;
- zero on the diagonal (no self-edges).

Duplicate entries are summed and stored zeros dropped (a zero is no edge).
Anything else is refused before writing:

```text
Connectivity matrix must be exactly symmetric.
```

A directed or asymmetric matrix (`adata.obsp["distances"]`, a kNN graph before
symmetrization) is not a connectivity graph; symmetrize it first, for example
`(g + g.T) / 2`, and clear its diagonal with `g.setdiag(0)`.

The viewer finds `obsp['connectivities']` by its own rule and reads the graph
only when a view draws edges. Its arrays are read whole, in chunks of 262,144
values ({doc}`09_cellucid_ready_anndata_specification`).

## Next

- {doc}`08_vector_fields_velocity_displacement`
- {doc}`../../web_app/c_core_interactions/index`
- The same topic for R, `cellucid_prepare()`: {doc}`../../r_package/c_data_preparation_api/07_connectivities_knn_graph`
