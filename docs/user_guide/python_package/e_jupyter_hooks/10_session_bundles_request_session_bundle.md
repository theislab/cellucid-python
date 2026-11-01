# Session bundles: capture from the browser (`request_session_bundle`)

Session bundles are the “durable artifact” for notebook reproducibility:
- hooks/events are **live** and transient,
- a `.cellucid-session` bundle is a **saved snapshot** of meaningful UI state.

This page documents the **no-download** workflow:

> Python asks the viewer to build a session bundle, the viewer sends it back on its channel, and Python returns it as a `CellucidSessionBundle` object.

```{tip}
To bring highlights, labels, scores and gene sets onto your AnnData you do
not need a bundle at all: `viewer.pull(adata)` writes everything the viewer
holds in one call ({doc}`../g_api_reference_coverage/api/round_trip`). A
bundle is for a durable, shareable snapshot of the whole session.
```

## At a glance

**Audience**
- Wet lab / beginner: you can usually skip this page unless you are collaborating or exporting results.
- Computational users: this is the “bridge UI state back to Python” tool.
- Developers: read the protocol notes and size limits.

**Primary APIs** -- one way, in a notebook widget and a browser-tab viewer
alike (R: `cellucid_request_session()`, then `cellucid_session_result()`):
- `request = viewer.request_session_bundle()`, then `bundle = request.result()`
  (in a notebook, in the next cell)
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

### Step 1 — Ask for the bundle (no browser download)

A notebook widget hears the viewer only **between cells**, so ask in one cell:

```python
request = viewer.request_session_bundle()
```

### Step 2 — Collect it in the next cell

```python
bundle = request.result()
bundle
```

`request.done` says whether the answer arrived. In a notebook widget the
answer cannot arrive while the cell that asked runs, so `result()` there says
so at once. A browser-tab viewer (`AnnDataViewer`) answers while a cell runs,
so there `result(timeout=60)` may follow in the same cell:

```python
bundle = viewer.request_session_bundle().result(timeout=60)
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

When you call `viewer.request_session_bundle()`:

1. Python registers an expected `(viewerId, requestId)` pair for the viewer's
   channel.
3. Python sends a command to the viewer on its channel:
   - `{ "type": "requestSessionBundle", "requestId": "..." }`
4. The viewer creates a Blob (`.cellucid-session` bytes).
5. The viewer sends the bytes back on the same channel as one binary frame
   (header `{"type": "session_bundle", "requestId": "..."}`, then the bundle);
   in a notebook the widget carries it in pieces over the notebook's connection.
6. Python streams the frame to a temp file, validates it, and emits a
   `session_bundle` event.
7. `request.result()` returns `CellucidSessionBundle(Path(path))`.
   If the viewer could not make the bundle, or the bundle is invalid, a
   `command_error` naming `requestSessionBundle` arrives instead and
   `result()` raises `RuntimeError` with its reason.

If the network drops while the request is pending, the viewer reopens its
channel and sends the bundle on the new socket; the request still answers.

## Limits and performance

- No fixed size cap: the frame is streamed to disk in chunks (memory stays
  bounded), so the only limit is the free space of the temporary directory.
  Each chunk inside the bundle is bounded by its own declared size when it is
  read.
- A request waits until it is collected or cancelled (`request.cancel()`, or
  the viewer stops); `result(timeout=...)` bounds only how long that call
  waits, and a later `result()` still collects it.

## Edge cases

- **Viewer not loaded yet**: the request is delivered once the viewer
  connects; collect it in a later cell.
- **Stale output after a kernel restart**: it keeps its last view and says
  the notebook stopped it; re-run the viewer cell.
- **Very large sessions**: the upload needs as much free space in the temporary
  directory as the bundle holds.

## Troubleshooting

### Symptom: `request.result()` raises at once in a notebook

It says the widget hears the viewer only between cells: run
`request.result()` again in the next cell.

### Symptom: `request.result()` says the bundle has not arrived

Likely causes:
- the viewer has not loaded yet, or is not shown in any output,
- the viewer's channel is not connected.

How to confirm:
- run `viewer.debug_connection()` and check `channel_connected` and
  `connected_outputs` (a browser-tab viewer also reports a ping/pong)

Fix:
- follow {doc}`14_troubleshooting_hooks`

### Symptom: “Dropped a session bundle nobody is waiting for” in the log

A bundle arrived after its request was cancelled (`request.cancel()`, or the
viewer stopped). The
bundle is discarded and the channel stays open.

Fix:
- ask again; in a notebook collect it with `request.result()` in a later cell.

## Next steps

- Apply session bundles back to AnnData: {doc}`11_session_bundles_apply_session_to_anndata`
- Full troubleshooting: {doc}`14_troubleshooting_hooks`
- What a session holds, from the web app's side: {doc}`../../web_app/l_sessions_sharing/02_what_gets_saved_and_restored`
