# Jupyter embedding + hooks sessions (gallery)

Notebooks focused on embedding the viewer, Python↔frontend hooks, and session/session-bridge workflows.

If you want the conceptual docs (schemas, environment matrix, troubleshooting), start with:
- {doc}`../e_jupyter_hooks/index`

If you want the “hands-on but explained” tutorials, start with:
- {doc}`23_programmatic_highlighting_and_selection_callbacks`
- {doc}`32_session_persistence_and_restoring_analysis_artifacts`

## Notebook list (runnable)

- `jupyter_hooks_sessions_he_developmental.ipynb`
  - selection callbacks
  - highlighting back from Python
  - debugging “hooks don’t fire”
- `jupyter_session_bridge_he_developmental_complete.ipynb`
  - capture a session bundle from the viewer
  - apply it back onto an `AnnData` object
  - validate what changed (and what did not)

## Quick troubleshooting if hooks don’t fire

1) Make sure the viewer is fully loaded (wait for the “ready” event).
2) Call `viewer.debug_connection()` and check:
   - `server_health` is `ok`
   - `web_ui.cache.index_html_exists` is true
3) If you are on a remote Jupyter server, use a proxy (or SSH tunnel) so the browser can reach the kernel-side port.

```{tip}
The fastest “is the plumbing alive?” check is:

`viewer.debug_connection()`
```

```{toctree}
:maxdepth: 1
:hidden:

jupyter_session_bridge_he_developmental_complete
jupyter_hooks_sessions_he_developmental
```
