# Web App Loading Issues

**Audience:** everyone  
**Time:** 10 minutes  
**What you'll get:** why the browser will not open a folder that exists

Use this page when your export folder exists but Cellucid won’t load it (or loads it incorrectly).

```{note}
The web app has a comprehensive troubleshooting page for loading here:
{doc}`../../web_app/b_data_loading/08_troubleshooting_data_loading`.
This page is the R-export-focused shortcut.
```

## Symptom: “I can’t select a folder / the picker doesn’t show folders”

**Likely causes**
- You denied directory access.
- An embedded or managed-browser policy blocked directory selection.
- You are on mobile; Cellucid's supported loading workflows are desktop
  workflows.

**Fix**
- Open the standalone app in a supported current desktop browser and click
  **Prepared**.
- If organizational policy forbids directory selection, choose server mode:
  {doc}`../../web_app/b_data_loading/04_server_tutorial`

## Symptom: “Dataset loads forever / spinner never stops”

**Likely causes**
- Extremely large exports (many genes/files)
- Browser struggling to scan the directory tree

**Fixes**
- Export fewer genes (gene panel).
- Use quantization and compression.
- Consider server mode for large datasets.

## Symptom: “Canvas is blank / WebGL context lost”

**Likely causes**
- GPU/WebGL issues
- Dataset too large for GPU memory

**Fix**
- Confirm hardware acceleration and WebGL2 are enabled, then update the browser,
  operating system, and GPU driver when applicable.
- Reduce dataset size (subset cells, reduce layers).
- See system requirements: {doc}`../../web_app/a_orientation/02_system_requirements`

## Symptom: “CORS blocked” when loading from a hosted URL

**Likely cause**
- You hosted the files on a server that does not send the required CORS headers.

**Fix**
- Use the supported GitHub exports workflow ({doc}`../d_viewing_loading/02_host_exports_for_sharing`).
- Or use server mode (Python) which handles headers appropriately.

## Symptom: “It loads but labels/fields are wrong”

**Likely cause**
- Row order mismatch in the export inputs.

**Fix**
- Re-export after aligning all inputs to the same `cells` order.
- See: {doc}`../c_data_preparation_api/02_input_requirements_global`

## Next steps

- Validate the export folder: {doc}`../d_viewing_loading/03_validate_exports_and_debug_loading`
