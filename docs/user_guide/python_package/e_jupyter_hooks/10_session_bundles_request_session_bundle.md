# Sessions: capture from the browser (`request_session`)

A session is the “durable artifact” for notebook reproducibility:
- hooks/events are **live** and transient,
- a `.cellucid-session` file is the viewer's **history**: every step a
  person, an agent or a notebook took, the trackers and branches they made,
  and the state each step shows.

This page documents the **no-download** workflow:

> Python asks the viewer for its session, the viewer sends the file back on its channel, and Python returns it as a {class}`~cellucid.CellucidSession`.

```{tip}
To bring highlights, labels, scores and gene sets onto your AnnData you do
not need a session at all: `viewer.pull(adata)` writes everything the viewer
holds in one call ({doc}`../g_api_reference_coverage/api/round_trip`). A
session is for a durable, shareable record of the whole history, and for
applying the state of any step of it.
```

## At a glance

**Audience**
- Wet lab / beginner: you can usually skip this page unless you are collaborating or exporting results.
- Computational users: this is the “bridge UI state back to Python” tool.
- Developers: read the protocol notes and size limits.

**Primary APIs** -- one way, in a notebook widget and a browser-tab viewer
alike (R: `cellucid_request_session()`, then `cellucid_session_result()`):
- `request = viewer.request_session()`, then `session = request.result()`
  (in a notebook, in the next cell)
- `session.save("path/to/file.cellucid-session")`
- `viewer.load_session(path)`: the other way, a saved session's trackers
  added to the viewer's history

## What is a session?

A session is a single binary file (`*.cellucid-session`, session format 4)
written by the web app. It holds the dataset's history: its trackers (each an
independent tree of steps, with named branches), every step (which action it
applied, with which parameters, who took it, its label), and the values of
the viewer's state each step shows -- labels and categories, user-defined
fields, highlight pages and groups, gene sets, filters, colours, views,
camera, layout. The chat bar's conversations are left out of a session a
notebook asks for.

In Python, it is represented by {class}`~cellucid.CellucidSession`, read and
checked whole when it is made (every chunk, array and value against its
digest, every record against the history model). `session.tracker()` is the
tracker the viewer showed (a {class}`~cellucid.history.Tracker`): its
`steps` table (one row per step: id, branch, label, who took it, when),
`step(id).state()` (the state after any step), `branch(name)`, `turns`,
`notes_table()`, `methods()` and `pseudocode()`; `session.trackers` lists
every tracker -- see {doc}`../g_api_reference_coverage/api/sessions`. A file
of session format 3, saved by an earlier Cellucid, is refused with the same
sentence the web app gives; a file written by a newer Cellucid than the
package is refused with "This session was written by a newer Cellucid (X)
than this package reads (Y); update the cellucid package to read it."

## Step-by-step: capture a session in a notebook

### Step 1 — Ask for the session (no browser download)

A notebook widget hears the viewer only **between cells**, so ask in one cell:

```python
request = viewer.request_session()
```

### Step 2 — Collect it in the next cell

```python
session = request.result()
session.trackers
```

`request.done` says whether the answer arrived. In a notebook widget the
answer cannot arrive while the cell that asked runs, so `result()` there says
so at once. A browser-tab viewer (`AnnDataViewer`) answers while a cell runs,
so there `result(timeout=60)` may follow in the same cell:

```python
session = viewer.request_session().result(timeout=60)
```

You now have a Python object that points to a temp file on disk:

```python
print(session.path)
session.tracker().steps[["id", "branch", "label", "origin"]]
```

### Step 3 — Save it somewhere persistent (recommended)

The returned session's path is a temp file. If you want to keep it:

```python
session.save("./my-session.cellucid-session")
```

## How it works (exact, but readable)

When you call `viewer.request_session()`:

1. Python registers an expected `(viewerId, requestId)` pair for the viewer's
   channel.
2. Python sends a command to the viewer on its channel:
   - `{ "type": "requestSessionBundle", "requestId": "..." }`
3. The viewer writes its history as a `.cellucid-session` file of format 4.
4. The viewer sends the bytes back on the same channel as one binary frame
   (header `{"type": "session_bundle", "requestId": "..."}`, then the file);
   in a notebook the widget carries it in pieces over the notebook's connection.
5. Python streams the frame to a temp file, reads and checks it whole, and
   emits a `session_bundle` event.
6. `request.result()` returns `CellucidSession(Path(path))`.
   If the viewer could not make the session, or the file is invalid, a
   `command_error` naming `requestSessionBundle` and the request arrives
   instead and `result()` raises `RuntimeError` with its reason.

If the network drops while the request is pending, the viewer reopens its
channel and sends the session on the new socket; the request still answers.

## The other way: load a session file into the viewer

`viewer.load_session(path)` sends a `.cellucid-session` file (or a
{class}`~cellucid.CellucidSession`) to the viewer as one binary
`load_session` command. Its trackers are added beside the person's own
(nothing of theirs is replaced), and the viewer shows the file's active
tracker. The session must be of the dataset the viewer shows.

```python
loading = viewer.load_session("my-session.cellucid-session")
```

```python
loading.result()   # {"tracker": {"id": ..., "name": ...}, "head": "s7"}
```

The viewer answers with a `session_loaded` event (and a `history_changed` of
kind `loaded`), or with a `command_error`: `result()` then raises
`RuntimeError` with the viewer's reason (another dataset's session, a damaged
file). As with `request_session()`, a notebook widget hears the answer only
between cells; `loading.done` says whether it arrived and `loading.cancel()`
stops waiting. R: `cellucid_load_session()`, then
`cellucid_load_session_result()`.

## Limits and performance

- No fixed size cap: the frame is streamed to disk in chunks (memory stays
  bounded), so the only limit is the free space of the temporary directory.
  Each chunk inside the file is inflated to exactly its declared size when it
  is read.
- A request waits until it is collected or cancelled (`request.cancel()`, or
  the viewer stops); `result(timeout=...)` bounds only how long that call
  waits, and a later `result()` still collects it.

## Edge cases

- **Viewer not loaded yet**: the request is delivered once the viewer
  connects; collect it in a later cell.
- **Another dataset shown**: while the person shows a dataset the notebook
  did not open in the viewer, the viewer refuses to send that dataset's
  session; the request fails with the reason.
- **Stale output after a kernel restart**: it keeps its last view and says
  the notebook stopped it; re-run the viewer cell.
- **Very large sessions**: the upload needs as much free space in the
  temporary directory as the file holds.

## Troubleshooting

### Symptom: `request.result()` raises at once in a notebook

It says the widget hears the viewer only between cells: run
`request.result()` again in the next cell.

### Symptom: `request.result()` says the session has not arrived

Likely causes:
- the viewer has not loaded yet, or is not shown in any output,
- the viewer's channel is not connected.

How to confirm:
- run `viewer.debug_connection()` and check `channel_connected` and
  `connected_outputs` (a browser-tab viewer also reports a ping/pong)

Fix:
- follow {doc}`14_troubleshooting_hooks`

### Symptom: “Dropped a session bundle nobody is waiting for” in the log

A session arrived after its request was cancelled (`request.cancel()`, or the
viewer stopped). The file is discarded and the channel stays open.

Fix:
- ask again; in a notebook collect it with `request.result()` in a later cell.

## Next steps

- Apply a session back to AnnData: {doc}`11_session_bundles_apply_session_to_anndata`
- Full troubleshooting: {doc}`14_troubleshooting_hooks`
- What a session holds, from the web app's side: {doc}`../../web_app/l_sessions_sharing/02_what_gets_saved_and_restored`
