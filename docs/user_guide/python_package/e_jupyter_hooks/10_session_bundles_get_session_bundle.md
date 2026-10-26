# Session bundles: capture from the browser (`get_session_bundle`)

Session bundles are the “durable artifact” for notebook reproducibility:
- hooks/events are **live** and transient,
- a `.cellucid-session` bundle is a **saved snapshot** of meaningful UI state.

This page documents the **no-download** workflow:

> Python asks the viewer to build a session bundle, the viewer sends it back on its channel to the local server, and Python returns it as a `CellucidSessionBundle` object.

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

1. Python displays the viewer when needed and waits for its exact ready event;
   display or readiness failure is propagated.
2. Python registers an expected `(viewerId, requestId)` pair on the local server with a TTL.
3. Python sends a command to the viewer on its channel:
   - `{ "type": "requestSessionBundle", "requestId": "..." }`
4. The viewer creates a Blob (`.cellucid-session` bytes).
5. The viewer sends the bytes back on the same channel as one binary frame
   (header `{"type": "session_bundle", "requestId": "..."}`, then the bundle).
6. The server streams the frame to a temp file, validates it, and emits a
   `session_bundle` event to Python.
7. `get_session_bundle()` unblocks and returns `CellucidSessionBundle(Path(path))`.
   If the viewer could not make the bundle, or the bundle is invalid, a
   `command_error` naming `requestSessionBundle` arrives instead and
   `get_session_bundle()` raises `RuntimeError` with its reason at once.

If the network drops while the request is pending, the viewer reopens its
channel and sends the bundle on the new socket; the request still answers.

## Limits and performance

- No fixed size cap: the frame is streamed to disk in chunks (memory stays
  bounded), so the only limit is the free space of the temporary directory.
  Each chunk inside the bundle is bounded by its own declared size when it is
  read.
- The request lasts as long as `get_session_bundle(timeout=...)` waits (one
  hour with `timeout=None`); a bundle that arrives after that is dropped.

## Edge cases

- **Viewer not displayed**: in Jupyter, `get_session_bundle()` displays it and
  propagates any display failure.
- **Viewer generation unavailable**: viewer construction/startup raises, so no
  session request is issued.
- **Stale iframe after kernel restart**: its channel is refused and the view
  freezes; re-run the viewer cell.
- **Very large sessions**: the upload needs as much free space in the temporary
  directory as the bundle holds.

## Troubleshooting

### Symptom: it hangs until timeout

Likely causes:
- the viewer's channel is not connected (proxy/mixed-content; the viewer says
  so),
- the viewer is not displayed anywhere (the request waits for it),
- the notebook was re-run and the pending request expired.

How to confirm:
- run `viewer.debug_connection()` and check `channel_connected` and ping/pong
- open browser devtools → Network, filter `WS`, and look at the
  `/_cellucid/ws` connection's frames

Fix:
- follow {doc}`14_troubleshooting_hooks`

### Symptom: “Dropped a session bundle nobody is waiting for” in the log

A bundle arrived after its `get_session_bundle()` call had stopped waiting
(its timeout passed). The bundle is discarded and the channel stays open.

Fix:
- call `viewer.get_session_bundle()` again, with a longer `timeout=` for a
  large session.

## Next steps

- Apply session bundles back to AnnData: {doc}`11_session_bundles_apply_session_to_anndata`
- Full troubleshooting: {doc}`14_troubleshooting_hooks`
- What a session holds, from the web app's side: {doc}`../../web_app/l_sessions_sharing/02_what_gets_saved_and_restored`
