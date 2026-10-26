# Debugging mental model (where to look first)

**Audience:** everyone (beginners: follow the checklist; experts: jump to the relevant symptom)  
**Time:** 10–30 minutes  
**Goal:** have a systematic, low-guesswork workflow for debugging viewer/server/hooks/session issues.

---

## The debugging rule: separate “UI problems” from “server problems”

Most Cellucid failures fall into one of these buckets:

1) **The exact viewer generation cannot be established or served**
2) **The browser can’t reach the Python server** (port/network/tunnel/proxy)
3) **The dataset can’t be interpreted** (no plottable embedding, an undecodable array, stale Zarr metadata)
4) **Hooks/events aren’t being delivered** (the viewer channel is not open, a hook raised, or a different viewer object)
5) **Session capture/apply fails** (not ready, requestId gating, exact dataset fingerprint mismatch)

The fastest way to debug is to test these buckets in order.

---

## Step 0 (do this first): `viewer.debug_connection()`

If you have a `viewer` object, run:

```python
report = viewer.debug_connection()
report
```

This returns a structured dict that checks:
- server probes (`/_cellucid/health`, `/_cellucid/info`, `/_cellucid/datasets`)
- the served store of every listed dataset, keyed by its exact id
- whether the verified viewer UI (`/index.html`) is being served
- whether the viewer's channel (its WebSocket) is connected (`channel_connected`)
- a Python → frontend roundtrip (ping/pong)
- a frontend debug snapshot (iframe URL/origin/user agent)
- recent accepted-event counts grouped by exact event type

### How to read the report (high signal fields)

- `server_health` / `server_health_error`
  - If this fails: the server isn’t reachable (port/network problem).
- `anndata_store_probes`
  - Every listed dataset id maps to its store's URL and the answer to a `HEAD`
    and a one-byte `Range` read (a readable store answers `206` with
    `Content-Range`), or to that probe's explicit error.
- `viewer_index_probe` / `viewer_index_probe_error`
  - If this fails: the active verified generation is not being served.
- `frontend_roundtrip`
  - If `ok=False`: Python→frontend messaging is broken or the iframe isn’t alive.
- `frontend_debug_snapshot`
  - Confirms the live iframe URL, origin, server URL, connection state, parent
    origin, and browser identification.
- `recent_events`
  - Confirms which accepted event types reached Python.
  - `command_errors` lists the commands the viewer refused. A non-empty list
    means messaging works and the *argument* was wrong.
- `client_server_url`
  - This must be the exact HTTP(S) base that the browser can reach.

```{tip}
When asking for help, paste the `report` (redacting private paths/dataset names) instead of describing symptoms loosely.
```

---

## Step 1: confirm the server is alive (no browser required)

In Python, check:

```python
print(viewer.server_url)
```

Then probe the health endpoint:

```python
import urllib.request, json
json.loads(urllib.request.urlopen(viewer.server_url + "/_cellucid/health").read())
```

Expected:
- a JSON payload like `{"status": "ok", ...}`

If this fails:
- the server may not be running,
- the port may be blocked/in use,
- you might be on a remote kernel whose loopback isn’t reachable from your browser (see Step 3).

---

## Step 2: confirm the viewer UI is being served

Cellucid notebook/server mode serves one exact verified web generation from the
Python server.

Check:
- `viewer.server_url + "/index.html"`

Viewer construction establishes that generation before binding the server and
raises on failure. If the published generation is removed while a server is
running, asset requests return an explicit service error. See
{doc}`06_privacy_security_and_offline_vs_online`.

---

## Step 3: confirm the browser can reach the server (remote notebooks)

The most common “it works locally but not on JupyterHub/HPC” issue:
- the browser cannot reach the kernel’s loopback server (`127.0.0.1:<port>`), or
- the notebook page is HTTPS and blocks HTTP loopback (mixed content).

Fix patterns:

### Option A (recommended): `jupyter-server-proxy`

If installed and enabled, the notebook can proxy the viewer through the
notebook server origin. Two things must be passed explicitly, because Cellucid
detects neither:

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
rather than quietly ignored.

If this step is missing, every request through the proxy returns
**421 Misdirected Request** and the server log records
`Refused a request whose Host header does not name this server`.

### Option B: SSH port forwarding (robust, works everywhere)

If your kernel runs on a remote machine:

1) Start the viewer on a fixed port:

   ```python
   from cellucid import AnnDataViewer

   viewer = AnnDataViewer(
       "data.h5ad",
       port=8765,
       dataset_name="Example",
       dataset_id="example",
       client_server_url="http://127.0.0.1:8765",
   )
   print(viewer.viewer_url)
   ```

2) On your laptop:

   ```bash
   ssh -N -L 8765:127.0.0.1:8765 <user>@<remote-host>
   ```

Now the browser hits your laptop’s localhost port, which tunnels to the remote kernel.

See the detailed guide: {doc}`../../web_app/b_data_loading/05_jupyter_tutorial`

---

## Step 4: debug hooks/events (frontend → Python)

### Minimal diagnostic hook

```python
@viewer.on_message
def debug(event):
    print(event)
```

If you see nothing when interacting in the UI:
- the viewer's channel (`/_cellucid/ws`) is not open
  (`viewer.debug_connection()["channel_connected"]` is `False`), or
- you are looking at an iframe of a different viewer object than the one that
  registered the hooks.

If a hook of yours raises, the exception is not raised in your cell: hooks run
on the viewer channel's thread, so it is logged by the `cellucid.viewer_channel`
logger ("A viewer event hook failed", with the traceback) and the channel stays
open. The event is recorded in `viewer.state` and handed to
`wait_for_event(...)` before any hook runs, and every other hook still runs
for it, so one failing hook costs neither (nor `get_session_bundle()`, which
waits for its event). The record is logged at
`ERROR` level, which Python writes to stderr even when logging is not
configured (in a notebook it may land in the kernel's log rather than the
cell); a handler on that logger collects it where you want it.

### Browser-side confirmation

Open browser devtools → Network, filter `WS`, and select the `/_cellucid/ws`
connection:
- it should be `101 Switching Protocols`, with `hello` and `welcome` as its
  first frames;
- every selection, hover or click adds a frame.

If there is no connection, or it closed:
- a proxy in between does not pass WebSockets (the viewer says so),
- or the server URL is unreachable.

If the connection was closed with code 4401, the iframe is older than the
kernel (a restart): the viewer froze and says Python stopped.

Mitigation:
- recreate the viewer,
- avoid copying stale viewer URLs across sessions.

---

## Step 5: debug Python → frontend commands

If you send commands like:

```python
viewer.highlight_cells([1, 2, 3])
viewer.set_color_by("cell_type")
```

and nothing happens:
- the iframe may not be alive,
- the notebook output may have been cleared/collapsed,
- the viewer may not be ready yet.

Use `viewer.debug_connection()` and check:
- `frontend_roundtrip`

If `frontend_roundtrip.ok` is false:
- call `viewer.display()` again,
- or recreate the viewer.

---

## Step 6: debug session capture (`get_session_bundle`)

Common failure modes:

### Timeout waiting for bundle

Likely causes:
- viewer never reached ready (it reports ready once the notebook's own
  dataset is on screen; a different dataset chosen in the viewer before that —
  a file, a sample, **None** — replaces it, and that viewer then never reports
  ready),
- iframe not displayed,
- the viewer's channel is not connected, due to proxy/mixed-content issues.

How to confirm:
- `viewer.debug_connection()` (look at `state.ready`, `frontend_roundtrip`,
  `frontend_debug_snapshot`, and `recent_events`)
- browser network tab (filter `WS`) for the `requestSessionBundle` command and
  the binary frame that answers it on `/_cellucid/ws`

Fix:
- call `viewer.wait_for_ready(timeout=60)` first,
- ensure the viewer is displayed,
- fix remote notebook connectivity (Step 3).

### “Dropped a session bundle nobody is waiting for” (log)

Meaning:
- a bundle arrived after its `get_session_bundle()` call stopped waiting (its
  timeout passed).

Fix:
- request again from the current viewer instance, with a longer `timeout=`.

For more: {doc}`05_sessions_to_anndata_bridge`

---

## Step 7: debug “the dataset looks wrong” (identity + ordering)

If selection indices don’t match the cells you expect in Python, the usual culprit is:
- dataset row order changed,
- you are analyzing a different AnnData object than the one the viewer is showing.

Checks:
- confirm which AnnData was passed to `show_anndata(...)`,
- if applying sessions, inspect `bundle.dataset_fingerprint` and `expected_dataset_id`.

See: {doc}`04_dataset_identity_and_reproducibility`

---

## Troubleshooting index (symptom → fix)

### Symptom: “Blank viewer area”

Likely causes:
- the configured source generation could not be fetched, verified, or
  published during viewer construction,
- mixed-content blocked (HTTPS notebook),
- server not reachable.

Fix:
- run `viewer.debug_connection()` and follow the first failing probe.

### Symptom: “Viewer loads but hooks don’t work”

Likely causes:
- the viewer channel (`/_cellucid/ws`) did not open, or closed,
- an iframe left over from before a kernel restart.

Fix:
- check browser Network tab,
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
