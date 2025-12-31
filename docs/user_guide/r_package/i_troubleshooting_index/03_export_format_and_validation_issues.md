# Export Format and Validation Issues

Use this page when:
- `cellucid_prepare()` ran, but the export folder looks incomplete, or
- the web app fails with “file not found”, “invalid manifest”, or similar.

## Symptom: required files are missing

Required for any export:
- `dataset_identity.json`
- `obs_manifest.json`
- at least one `points_*d.bin[.gz]`
- `obs/` directory with files

**How to confirm**
```r
list.files(out_dir)
```

**Fix**
- Re-export with `force=TRUE`.
- Export to a fresh `out_dir`.

## Symptom: “I only see `.gz` files but the viewer looks for `.bin`”

**Likely cause**
- confusion about compression.

If you set `compression=...`, the exporter writes:
- `points_2d.bin.gz` (not `points_2d.bin`)

The manifests record the correct suffixes. The viewer should follow manifests.

**Fix**
- Don’t rename files manually.
- Re-export if you modified the folder.

## Symptom: `var_manifest.json` exists but `var/` is empty

**Likely causes**
- export was interrupted mid-way
- permission issues in the output directory

**Fix**
- delete the dataset folder and re-export, or set `force=TRUE` to overwrite

## Symptom: gene files exist but gene search doesn’t find them

**Likely causes**
- gene identifiers in `var` aren’t what you expect (`rownames(var)` missing)
- filename sanitization collisions overwrote files

**How to confirm**
- inspect `var_manifest.json` and compare the gene IDs listed to what you expect.

**Fix**
- set `rownames(var)` or `var_gene_id_column` explicitly
- ensure gene IDs are unique after sanitization

## Symptom: “GitHub connect says `datasets.json not found`”

This is not a `cellucid-r` export bug; it’s a hosting-layout issue.

The GitHub loader expects an **exports root**:

```
exports/
  datasets.json
  <dataset>/
    dataset_identity.json
    ...
```

Fix:
- generate and commit `datasets.json` (see {doc}`../d_viewing_loading/02_host_exports_for_sharing`)
- see the full GitHub loader tutorial: {doc}`../../web_app/b_data_loading/02_local_demo_tutorial`

## Validation tools and deep debugging

- Format specification: {doc}`../c_data_preparation_api/09_output_format_specification_exports_directory`
- Validation checklist + binary inspection: {doc}`../d_viewing_loading/03_validate_exports_and_debug_loading`
