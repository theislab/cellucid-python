# Screenshots

This page is a **screenshot capture checklist** for the Data Loading section.

It exists so you (or a collaborator) can capture screenshots once, systematically, without hunting through every tutorial page.

---

## How placeholders work

All pages in this section currently use:

- `cellucid-python/docs/_static/screenshots/placeholder-screenshot.svg`

Each placeholder is preceded by an HTML comment with:

- what state to capture,
- what to crop/redact,
- suggested caption/alt text,
- suggested filename conventions.

General guidance lives in:

- `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`

Recommended storage location for this section:

- `cellucid-python/docs/_static/screenshots/data_loading/`

---

## Recommended screenshot set (data loading)

### Loader panel (orientation) (1 screenshot)

Capture the “entry point” where users choose local files vs remote server vs GitHub.

Page:
- `b_data_loading/index` (ID: `data-loading-dataset-connections-panel`)

### Overview page (2 screenshots)

1) Data loading controls highlighted (empty state)
2) Successful load state (dataset loaded + visible legend/field selector)

Page:
- `b_data_loading/01_loading_options_overview` (IDs: `data-loading-overview-loader-panel`, `data-loading-overview-success-state`)

### Browser file picker (2 screenshots)

1) The three file picker buttons (folder / h5ad / zarr)
2) Successful exported-folder load state

Page:
- `b_data_loading/03_browser_file_picker_tutorial` (IDs: `data-loading-file-picker-entry-point`, `data-loading-export-folder-success`)

### GitHub exports workflow (2 screenshots)

1) Validating an export locally (folder picker success)
2) GitHub connection input + “connect” UX

Page:
- `b_data_loading/02_local_demo_tutorial` (IDs: `data-loading-github-validate-local-file-picker`, `data-loading-github-connect`)

### Server mode (1 screenshot)

Capture the terminal banner showing:
- server URL
- suggested viewer URL

Page:
- `b_data_loading/04_server_tutorial` (ID: `data-loading-server-terminal-banner`)

### Jupyter embedding (1 screenshot)

Capture the first successful embedded viewer output cell (with any “AnnData mode is slower” warning visible, if applicable).

Page:
- `b_data_loading/05_jupyter_tutorial` (ID: `data-loading-jupyter-embedded-viewer`)

---

## Highly recommended “failure mode” screenshots (debugging gold)

These save huge amounts of time in onboarding and support:

- GitHub connect error showing `datasets.json not found` (or a clear 404)
- CORS / mixed-content error message when connecting to a remote server
- File picker permission prompt (and what it looks like when denied)
- “WebGL context lost” overlay (see: `a_orientation/02_system_requirements`)

If you are documenting vector fields (velocity/drift), capture those screenshots in the overlay section:
- {doc}`../i_vector_field_velocity/08_screenshots`
