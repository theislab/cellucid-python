# Screenshots and diagrams checklist (hooks)

This page is for documentation authors/maintainers.

Goal: when you later add screenshots, they should be:
- reproducible,
- clearly captioned,
- safe to publish (no private data leakage),
- and consistent across pages.

Primary guide:
- `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`

## Where to store assets (recommended)

Create a topic folder:

```text
cellucid-python/docs/_static/screenshots/jupyter_hooks/
```

Naming rules:
- lowercase + hyphens
- prefix with an order number when step-by-step (`01_...`, `02_...`)

## Minimum required screenshots for this chapter

### “Success path” (must-have)

1. Viewer embedded inside a notebook output cell (orientation).
2. A selection made in the viewer (before → after).
3. A Python cell output showing selection indices received (proof hooks work).
4. A Python command affecting the viewer (e.g., highlight or reset camera).

### “Failure path” (high value)

1. “Viewer UI could not be loaded” error page (offline/no cache).
2. “Notebook proxy required” error page (HTTPS mixed content).
3. Browser devtools → Network filtered to `/_cellucid/events` showing a failure (blocked/404/413).

## Minimum required diagrams

1. Architecture overview showing:
   - Python → iframe via `postMessage`
   - iframe → Python via `POST /_cellucid/events`

Placeholders are acceptable until a proper diagram is produced.

## Copy/paste placeholder block (MyST `{figure}`)

Use this template so the capture instructions live next to the figure:

````md
<!-- SCREENSHOT PLACEHOLDER
ID: <unique-id>
Suggested filename: jupyter_hooks/<filename>.png
Where it appears: <page → section>
Capture:
  - UI location: <notebook output / viewer sidebar / browser devtools / etc.>
  - State prerequisites: <what must already be true>
  - Action to reach state: <exact steps>
Crop:
  - Include: <what must be visible>
  - Exclude: <what must not be visible>
Redact:
  - Remove: <tokens, paths, emails, dataset names if private>
Annotations:
  - Callouts: <2–4 callouts max>
Alt text:
  - <one sentence>
Caption:
  - <one sentence: what the reader should learn>
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: <ALT TEXT>
:width: 100%

<CAPTION>
```
````

## Example placeholder: notebook-embedded viewer

<!-- SCREENSHOT PLACEHOLDER
ID: jupyter-hooks-screenshot-embedded-viewer
Suggested filename: jupyter_hooks/01_notebook-embedded-viewer.png
Where it appears: Jupyter Hooks → Quickstart
Capture:
  - UI location: notebook output cell
  - State prerequisites: viewer loaded successfully
  - Action to reach state: run `viewer = show_anndata(...)`
Crop:
  - Include: the iframe and enough notebook context to orient beginners
  - Exclude: local file paths, notebook filenames, private dataset names
Alt text:
  - Notebook cell output containing an embedded Cellucid viewer iframe.
Caption:
  - Cellucid running inside a notebook: the viewer is an iframe connected to a local data server.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder for a notebook cell output showing an embedded Cellucid viewer.
:width: 100%

Placeholder: embedded viewer inside a notebook output cell.
```

