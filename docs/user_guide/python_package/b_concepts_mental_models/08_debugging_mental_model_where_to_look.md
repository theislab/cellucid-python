# Debugging mental model (where to look first)

**Audience:** everyone (beginners: follow the checklist; experts: jump to the relevant symptom)  
**Time:** 10–30 minutes  
**Goal:** have a systematic, low-guesswork workflow for debugging viewer/server/hooks/session issues.

---

## The debugging rule: separate “UI problems” from “data and connection problems”

Most Cellucid failures fall into one of these buckets:

1) **The viewer UI does not appear** (the notebook does not show widgets, the
   browser cannot load the app; for a browser tab, the exact viewer generation
   cannot be established or served)
2) **The browser can’t reach a browser-tab server** (port/network/tunnel/proxy;
   a notebook widget has no server to reach)
3) **The dataset can’t be interpreted** (no plottable embedding, an undecodable array, stale Zarr metadata)
4) **Hooks/events aren’t being delivered** (a cell is still running, the channel is not open, a hook raised, or a different viewer object)
5) **Session capture/apply fails** (asked and collected in the same cell, not ready, exact dataset fingerprint mismatch)

The fastest way to debug is to test these buckets in order.

---

## Step 0 (do this first): `viewer.debug_connection()`

If you have a `viewer` object, run (in a new cell):

```python
report = viewer.debug_connection()
report
```

For a notebook widget (what `show_anndata()` returns) it reports, without
waiting for anything:
- `web_source_url` — where the browser loads the app from,
- `displayed` and `connected_outputs` — the outputs whose page reached the kernel,
- `channel_connected` — whether the viewer opened its channel,
- `dataset` — the id and format the widget serves,
- `state` — `ready` (None until the dataset opened) and the last event,
- `recent_events` — counts by event type; `command_errors` lists the commands
  the viewer refused (a non-empty list means messaging works and the
  *argument* was wrong).

For a browser-tab viewer (`AnnDataViewer`) it also probes its server:
- `server_health` / `server_health_error` — if this fails, the server isn’t
  reachable (port/network problem),
- `anndata_store_probes` — every listed dataset id maps to its store's URL and
  the answer to a `HEAD` and a one-byte `Range` read (a readable store answers
  `206` with `Content-Range`), or to that probe's explicit error,
- `viewer_index_probe` / `viewer_index_probe_error` — whether the active
  verified generation is being served,
- `frontend_roundtrip` — a ping/pong over the channel; `ok=False` means
  Python→frontend messaging is broken or the page isn't open,
- `frontend_debug_snapshot` — the page's URL, origin, server URL and
  connection state,
- `client_server_url` — the exact HTTP(S) base the browser must reach.

```{tip}
When asking for help, paste the `report` (redacting private paths/dataset names) instead of describing symptoms loosely.
```

---

## Step 1: confirm the viewer UI appears

In a notebook, the output shows the widget, which frames the app from
`web_source_url`:

- **an empty output, or a note that the notebook cannot show widgets**:
  install `anywidget` where the notebook *server* runs and reload the page, or
  use a front end that shows widgets
  ({doc}`../e_jupyter_hooks/03_supported_environments_matrix`);
- **“The page at … has not answered as a Cellucid viewer”**: the browser cannot
  load the app (offline, blocked address); allow it, or host a copy and pass
  `web_source_url=`.

A browser-tab viewer serves one exact verified web generation itself; check
`viewer.server_url + "/index.html"`. Its construction establishes that
generation before binding the server and raises on failure. See
{doc}`06_privacy_security_and_offline_vs_online`.

---

## Step 2: confirm a browser-tab server is alive and reachable

Skip this step for a notebook widget: it has no server.

```python
import urllib.request, json
json.loads(urllib.request.urlopen(viewer.server_url + "/_cellucid/health").read())
```

Expected: a JSON payload like `{"status": "ok", ...}`. If this fails, the
server may not be running or the port may be in use.

If the server is alive but the browser cannot reach it, the kernel runs on
another machine. In a notebook, use `show_anndata()` instead: it needs no
address. Otherwise:

### Option A: SSH port forwarding (robust, works everywhere)

1) Start the viewer on a fixed port:

   ```python
   from cellucid import AnnDataViewer

   viewer = AnnDataViewer(
       "data.h5ad",
       port=8765,
       dataset_name="Example",
       dataset_id="example",
   )
   print(viewer.viewer_url)
   ```

2) On your laptop:

   ```bash
   ssh -N -L 8765:127.0.0.1:8765 <user>@<remote-host>
   ```

Now the browser hits your laptop’s localhost port, which tunnels to the remote
kernel. Full guide:
{doc}`../d_viewing_apis/12_remote_servers_ssh_tunneling_and_cloud`.

### Option B: a reverse proxy

If a proxy (a hub's proxy route) carries the server to the browser, two things
must be passed explicitly, because Cellucid detects neither:

1. `client_server_url=` — the proxy’s browser-reachable base URL.
2. `allowed_hosts=` — the proxy front end’s own host name.

The second is required. The proxy forwards the browser’s `Host` header
unchanged, so your loopback-bound Cellucid server receives the hub’s host
(`hub.example.org`) instead of `127.0.0.1`. Cellucid answers only authorities
that name it, and refuses everything else with **421 Misdirected Request** —
that is what stops DNS rebinding (see
{doc}`06_privacy_security_and_offline_vs_online`). Naming the proxy host is how
you tell Cellucid that this particular foreign `Host` is yours.

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer(
    "data.h5ad",
    dataset_name="Example",
    dataset_id="example",
    client_server_url="https://hub.example.org/user/me/proxy/8765",
    allowed_hosts=["hub.example.org"],
)
```

From the command line:

```bash
cellucid serve /path/to/data --allowed-host hub.example.org
```

Each entry is one bare host name — a DNS name or an IP address, with no port,
no scheme, no path and no wildcard — and it matches whatever port the proxy
publishes. Repeat `--allowed-host` (or add list entries) for each additional
name. An unsupported entry is refused at startup with an explicit message
rather than quietly ignored. If this step is missing, every request through
the proxy returns **421 Misdirected Request** and the server log records
`Refused a request whose Host header does not name this server`. The proxy
must also pass WebSockets (the viewer channel).

---

## Step 3: debug hooks/events (frontend → Python)

### Minimal diagnostic hook

```python
@viewer.on_message
def debug(event):
    print(event)
```

If you see nothing when interacting in the UI:
- a cell is still running: a notebook widget's events arrive when the kernel
  is idle, all at once,
- the viewer's channel is not open
  (`viewer.debug_connection()["channel_connected"]` is `False`), or
- you are looking at the output of a different viewer object than the one
  that registered the hooks.

If a hook of yours raises, the exception is not raised in your cell: hooks run
on the viewer channel's thread, so it is logged by the `cellucid.viewer_channel`
logger ("A viewer event hook failed", with the traceback) and the channel stays
open. The event is recorded in `viewer.state` and handed to
`wait_for_event(...)` before any hook runs, and every other hook still runs
for it. The record is logged at
`ERROR` level, which Python writes to stderr even when logging is not
configured (in a notebook it may land in the kernel's log rather than the
cell); a handler on that logger collects it where you want it.

### Browser-side confirmation (browser tab)

Open browser devtools → Network, filter `WS`, and select the `/_cellucid/ws`
connection:
- it should be `101 Switching Protocols`, with `hello` and `welcome` as its
  first frames;
- every selection, hover or click adds a frame.

If the connection was closed with code 4401, the page is older than the
kernel (a restart): the viewer froze and says Python stopped. Recreate the
viewer; do not copy stale viewer URLs across sessions.

---

## Step 4: debug Python → frontend commands

If you send commands like:

```python
viewer.highlight_cells([1, 2, 3])
viewer.set_color_by("cell_type")
```

and nothing happens:
- the viewer may not have connected yet (commands wait and are delivered in
  order once it does),
- the output may have been cleared, or belong to an older viewer,
- the viewer may have refused it: check
  `viewer.debug_connection()["recent_events"]["command_errors"]`.

If nothing is refused and the channel is connected, show the viewer again
(end a cell with `viewer`) or recreate it.

---

## Step 5: debug session capture

### In a notebook

A widget hears the viewer only between cells: ask in one cell, collect in the
next.

```python
request = viewer.request_session_bundle()
```

```python
bundle = request.result()
```

If `request.result()` says the bundle has not arrived (in a notebook, only a
later cell can collect it):
- the viewer never reached ready (it reports ready once the notebook's own
  dataset is on screen; a different dataset chosen in the viewer before that —
  a file, a sample, **None** — replaces it, and that viewer then never reports
  ready),
- the viewer is not shown in any output, or its channel is not connected.

To bring results onto AnnData you need no bundle: `viewer.pull(adata)`.

### In a browser tab

`viewer.request_session_bundle().result(timeout=60)` waits. A “Dropped a
session bundle nobody is waiting for” log means the bundle arrived after its
request was cancelled; request again.

For more: {doc}`05_sessions_to_anndata_bridge`

---

## Step 6: debug “the dataset looks wrong” (identity + ordering)

If selection indices don’t match the cells you expect in Python, the usual culprit is:
- dataset row order changed,
- you are analyzing a different AnnData object than the one the viewer is showing.

Checks:
- confirm which AnnData was passed to `show_anndata(...)`,
- if applying sessions, inspect `bundle.dataset_fingerprint` and `expected_dataset_id`.

See: {doc}`04_dataset_identity_and_reproducibility`

---

## Troubleshooting index (symptom → fix)

### Symptom: “Empty viewer area”

Likely causes:
- the notebook does not show widgets, or the browser cannot load the app
  (Step 1),
- a browser-tab server's generation could not be established, or the server
  is not reachable (Step 2).

Fix:
- run `viewer.debug_connection()` and follow the first failing entry.

### Symptom: “Viewer loads but hooks don’t work”

Likely causes:
- a cell is still running (a widget's events arrive between cells),
- the viewer channel did not open, or closed,
- an output left over from before a kernel restart.

Fix:
- let the cell finish; check `channel_connected`,
- recreate the viewer.

### Symptom: “`apply_cellucid_session_to_anndata(...)` raised an error”

Likely causes:
- the required `expected_dataset_id`, cell count, or gene count does not match
  the bundle fingerprint,
- an output column already exists,
- the session does not contain a requested highlight or user-defined-field
  chunk.

Fix:
- read the raised exception; fingerprint mismatches and column collisions are
  rejected before mutation,
- pass the exact dataset ID for the AnnData object,
- inspect `bundle.list_chunk_ids()` before selecting which current chunks to
  apply.

---

## Next steps

- Privacy/offline model: {doc}`06_privacy_security_and_offline_vs_online`
- Hooks docs: {doc}`../e_jupyter_hooks/index`
