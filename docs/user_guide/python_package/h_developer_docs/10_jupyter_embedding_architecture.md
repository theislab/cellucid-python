# Jupyter embedding architecture

This page explains how Cellucid embeds the web app inside Jupyter/VSCode/Colab notebooks and how the bidirectional communication works.

If you only need the user-facing quickstart, start with:
- {doc}`../d_viewing_apis/index`
- {doc}`../e_jupyter_hooks/index`

---

## The mental model (one sentence)

A notebook viewer is:

> **a local Python server + an iframe pointed at that server + two communication channels** (postMessage commands and HTTP POST events).

---

## Components and responsibilities

### Viewer classes

Implementation:
- `cellucid-python/src/cellucid/jupyter.py`

Key classes:
- `BaseViewer`: shared notebook embedding + hooks/event handling
- `CellucidViewer`: serves an export folder via `CellucidServer`
- `AnnDataViewer`: serves AnnData via `AnnDataServer`

Convenience functions:
- `show(...)` → `CellucidViewer`
- `show_anndata(...)` → `AnnDataViewer`

### Identity and routing

Each viewer has:
- `viewerId`: routes frontend → Python events to the correct viewer
- `viewerToken`: included in Python → frontend commands (postMessage) as a lightweight authenticity token

Each viewer also typically runs its own server on its own port.

---

## How the iframe URL is constructed

The viewer URL is built from the **client-reachable server URL** plus query params:

- `jupyter=true`
- `viewerId=<id>`
- `viewerToken=<token>`
- `anndata=true` (only for AnnDataViewer)

Key method:
- `BaseViewer._get_client_server_url()`

### Why this is complicated

Notebook environments differ:

- Local classic/JupyterLab: browser can usually reach `http://127.0.0.1:<port>`.
- Remote kernels (JupyterHub, SSH to server, etc.): browser cannot reach kernel loopback.
- HTTPS notebooks: browser blocks HTTP loopback as mixed content.
- Colab: uses a special HTTPS port proxy (`google.colab.kernel.proxyPort`).

Cellucid tries to pick a URL that “actually works” by:
- preferring a proxy URL when direct loopback is unlikely,
- probing `/_cellucid/health` before committing to an iframe src,
- and showing an explicit inline error (iframe `srcdoc`) when a proxy is required.

Config override:
- `CELLUCID_CLIENT_SERVER_URL` can force the client URL (see {doc}`06_configuration_env_vars_and_logging`).

---

## Python → frontend: postMessage command channel

When you call:

```python
viewer.highlight_cells([1, 2, 3], color="#ff0000")
```

Python:
1) injects a small JS snippet into the notebook output,
2) finds the iframe by `viewerId`,
3) sends a `postMessage` to `iframe.contentWindow`,
4) includes `viewerId` + `viewerToken` in the message payload.

This channel is used for:
- highlighting cells
- changing color-by
- visibility toggles
- reset camera
- “freeze” before shutdown
- requesting a session bundle

Developer note:
- postMessage is the right tool for “commands” (immediate UI actions).
- it is not used for event delivery (events go over HTTP POST).

---

## Frontend → Python: HTTP POST event channel (hooks)

The embedded viewer posts events to:
- `POST /_cellucid/events`

The server routes events by `viewerId` to the correct viewer callback.

In Python, `BaseViewer`:
- records the event into `viewer.state` (latest snapshot),
- stores recent events in an internal ring buffer,
- triggers hook callbacks registered via:
  - `@viewer.on_selection`, `@viewer.on_hover`, `@viewer.on_click`, `@viewer.on_ready`

Protocol details are documented here:
{doc}`11_hooks_events_protocol_and_schema`.

---

## Session bundle capture (“no download” notebook workflow)

Notebook users often want to:
- interactively create highlights/annotations in the UI,
- then pull that state back into Python for downstream analysis.

Flow:

1) Python sends a command to the iframe:
   - `{type: "requestSessionBundle", requestId: "..."}`
2) The frontend serializes a `.cellucid-session` bundle in-memory.
3) The frontend uploads bytes to:
   - `POST /_cellucid/session_bundle?viewerId=...&requestId=...`
4) The server:
   - validates the bundle magic header,
   - streams the upload to a temp file,
   - routes a `session_bundle` event back to Python (includes path).
5) Python returns a `CellucidSessionBundle` handle to the caller.

This design keeps:
- memory bounded (streaming upload),
- notebook UX simple (no manual file download/upload),
- and allows reuse of the same session format as the web app.

---

## When things go wrong: the “proxy required” inline message

If the notebook is served from HTTPS/remote and direct loopback will not work, the iframe may show a “proxy required” message.

<!-- SCREENSHOT PLACEHOLDER
ID: notebook-proxy-required-message
Suggested filename: developer/notebook-proxy-required.png
Where it appears: Python Package → Developer Docs → Jupyter embedding architecture
Capture:
  - UI location: notebook output cell (iframe area)
  - State prerequisites: run notebook in an HTTPS environment without server proxy enabled
  - Action to reach state: create viewer, iframe tries proxy, fails health probe, shows srcdoc message
Crop:
  - Include: the full inline message and the recommended fixes
  - Exclude: personal notebook names, dataset identifiers
Alt text:
  - Notebook output showing an inline message stating that a proxy is required to load the viewer from a secure or remote notebook environment.
Caption:
  - When the browser cannot reach the kernel’s loopback server directly, Cellucid shows an inline “proxy required” message with concrete fixes (jupyter-server-proxy or CELLUCID_CLIENT_SERVER_URL).
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the inline notebook proxy-required message.
:width: 100%

Notebook proxy-required message: Cellucid explains how to make the server reachable from a secure/remote notebook origin.
```

---

## Troubleshooting

### Symptom: “The viewer iframe is blank”

Likely causes (ordered):
1) the web UI assets could not be loaded (offline + no cache),
2) mixed-content/proxy URL mismatch,
3) the server is not reachable from the browser (remote kernel),
4) the server died (port conflict, exception).

How to confirm:
- run `viewer.debug_connection()` and read:
  - `server_health`, `viewer_index_probe`, `frontend_roundtrip`.

Fix:
- set up `jupyter-server-proxy` or provide `CELLUCID_CLIENT_SERVER_URL`,
- ensure the UI cache is available (see {doc}`09_server_mode_architecture_endpoints_and_security`),
- check server logs for exceptions.

### Symptom: “Hooks worked once, then stopped”

Common causes:
- the notebook output cell was cleared/collapsed and the iframe was destroyed,
- the viewer was garbage-collected (lost references),
- you re-ran a cell and created a new viewer/server, but are watching the old iframe.

Fix:
- keep a reference to the viewer object,
- call `viewer.stop()` when done,
- create a new viewer after notebook reloads.
