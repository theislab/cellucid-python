# Session bundles: capture from the browser (`get_session_bundle`)

Session bundles are the “durable artifact” for notebook reproducibility:
- hooks/events are **live** and transient,
- a `.cellucid-session` bundle is a **saved snapshot** of meaningful UI state.

This page documents the **no-download** workflow:

> Python asks the viewer to build a session bundle, the viewer uploads it back to the local server, and Python returns it as a `CellucidSessionBundle` object.

## At a glance

**Audience**
- Wet lab / beginner: you can usually skip this page unless you are collaborating or exporting results.
- Computational users: this is the “bridge UI state back to Python” tool.
- Developers: read the protocol notes and size limits.

**Primary APIs**
- `bundle = viewer.get_session_bundle(timeout=...)`
- `bundle.save("path/to/file.cellucid-session")`

## What is a “session bundle”?

A session bundle is a single binary file (`*.cellucid-session`) produced by the web app that can contain:
- highlight groups and pages,
- user-defined categorical fields,
- other session state the UI considers meaningful and reproducible.

In Python, it is represented by:
- `cellucid.session_bundle.CellucidSessionBundle`

```{note}
The session bundle format is streaming-friendly: Python can read chunk metadata and decode only what it needs.
```

## Step-by-step: capture a session bundle in a notebook

### Step 1 — Make sure the viewer is ready

```python
viewer.wait_for_ready(timeout=60)
```

### Step 2 — Request the bundle (no browser download)

```python
bundle = viewer.get_session_bundle(timeout=60)
bundle
```

You now have a Python object that points to a temp file on disk:

```python
print(bundle.path)
print(bundle.list_chunk_ids())
```

### Step 3 — Save it somewhere persistent (recommended)

The default returned bundle path is typically a temp file.
If you want to keep it:

```python
bundle.save("./my-session.cellucid-session")
```

## How it works (exact, but readable)

When you call `viewer.get_session_bundle()`:

1. Python ensures the viewer is displayed/ready (best effort).
2. Python registers an expected `(viewerId, requestId)` pair on the local server with a TTL.
3. Python sends a command into the iframe:
   - `{ "type": "requestSessionBundle", "requestId": "..." }`
4. The frontend creates a Blob (`.cellucid-session` bytes).
5. The frontend uploads the bytes back to:
   - `POST /_cellucid/session_bundle?viewerId=...&requestId=...`
6. The server validates + streams the upload to a temp file and emits a `session_bundle` event back to Python.
7. `get_session_bundle()` unblocks and returns `CellucidSessionBundle(Path(path))`.

## Limits and performance

- Hard size cap: **512MB** per upload (guards against zip bombs and memory pressure).
- Upload is streamed to disk in chunks (keeps memory bounded).
- The request is TTL-bound; if you wait too long, the server may reject it as “no pending request”.

## Edge cases

- **Viewer not displayed**: `get_session_bundle()` will try to display the viewer (best effort in Jupyter).
- **Offline / no cached UI**: the iframe may show “viewer UI could not be loaded”, so no session capture can occur.
- **Stale iframe after kernel restart**: viewerId mismatch; re-run the viewer cell.
- **Very large sessions**: you can hit the 512MB cap; reduce stored state (fewer highlight groups, etc.).

## Troubleshooting

### Symptom: it hangs until timeout

Likely causes:
- the iframe can’t reach the server URL (proxy/mixed-content),
- the upload request is blocked by a browser policy or proxy,
- the notebook was re-run and the pending request expired.

How to confirm:
- run `viewer.debug_connection()` and check ping/pong
- open browser devtools → Network and look for `/_cellucid/session_bundle`

Fix:
- follow {doc}`14_troubleshooting_hooks`

### Symptom: “No pending session bundle request”

This means the server did not recognize the `(viewerId, requestId)` pair (expired, stale, or wrong viewer).

Fix:
- call `viewer.get_session_bundle()` again (new requestId),
- ensure you’re interacting with the current iframe (not a stale notebook output).

## Next steps

- Apply session bundles back to AnnData: {doc}`11_session_bundles_apply_session_to_anndata`
- Full troubleshooting: {doc}`14_troubleshooting_hooks`

