# Data preparation issues

## Start with the invariant named by the exception

Do not reshape, truncate, reorder, or fill scientific arrays merely to make
`prepare()` succeed. Check these identities directly:

- every embedding and vector field has one row per observation;
- expression rows follow the same observation order;
- expression columns follow the rows of `var`;
- connectivity is square, symmetric, and matches the observation count;
- `obs` and `var` columns hold types AnnData can store;
- dataset and gene identifiers are explicit and stable.

## `prepare()` succeeds but the viewer refuses the file

Read the viewer's notification: a refusal naming `uns['cellucid']` names the
key or value it does not accept. Then confirm the file reads in Python
(`anndata.read_h5ad(path)`, `anndata.read_zarr(path)`, or `mudata.read_h5mu`
/ `mudata.read_zarr` for a MuData written by `prepare(data=mdata)`), and for a hosted file
inspect the first failed request in the browser Network panel; its HTTP status
distinguishes a missing file from unreadable content.

## Values look scientifically wrong

Compare a small set of observation and gene identifiers before and after
`prepare()`. Treat transposition, silent reordering, category mismatch and
non-finite coordinates as data-contract failures. Fractional expression values
are rounded to float32 once, when written; integer counts read back exactly.

See {doc}`../c_data_preparation_api/02_input_requirements_global`,
{doc}`../c_data_preparation_api/09_cellucid_ready_anndata_specification`,
and {doc}`../c_data_preparation_api/11_troubleshooting_prepare`.
