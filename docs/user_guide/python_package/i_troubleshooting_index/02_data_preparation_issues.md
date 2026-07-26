# Data preparation issues

## Start with the invariant named by the exception

Do not reshape, truncate, reorder, or fill scientific arrays merely to make an
export complete. Check these identities directly:

- every embedding and vector field has one row per observation;
- expression rows follow the same observation order;
- expression columns follow the declared gene identifiers;
- connectivity is square and matches the observation count;
- categorical and continuous metadata use their declared types;
- dataset and gene identifiers are explicit and stable.

## Export completes but the viewer rejects the folder

Confirm that the selected directory is the export root and contains
`dataset_identity.json` and `obs_manifest.json`. Then inspect the first failed
request in the browser Network panel; its exact filename and HTTP status
distinguish a missing artifact from malformed content.

## Values look scientifically wrong

Compare a small set of observation and gene identifiers before and after export.
Treat transposition, silent reordering, category-code mismatch, non-finite
coordinates, and unexpected quantization as data-contract failures.

See {doc}`../c_data_preparation_api/02_input_requirements_global`,
{doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`,
and {doc}`../c_data_preparation_api/11_troubleshooting_prepare_export`.
