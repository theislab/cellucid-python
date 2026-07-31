# Privacy, security, and network requirements

**Audience:** everyone (especially anyone handling sensitive data)  
**Time:** 15–25 minutes  
**What you’ll learn:**
- When your data stays local vs when it can be exposed
- What the Python server is doing with `cellucid.com` (hosted UI assets)
- What source access the viewer requires in firewalled environments

---

## Mental model (one sentence)

Cellucid is a web app UI, but with `cellucid-python` your **dataset is served
by your Python process**. The main privacy/security decisions are where you
host that server and whether the Python runtime may reach the configured
viewer source.

---

## “Does my dataset get uploaded to cellucid.com?”

### Notebook/server modes (most `cellucid-python` usage)

By default:
- your dataset is served by a local Python server (usually `127.0.0.1:<port>`),
- the exact verified viewer generation is served from that same origin,
- the browser loads data from your local server.

Your dataset is **not uploaded** anywhere automatically.

### Web app file picker mode (no Python server)

If you open the Cellucid web app and use the browser file picker to select an export folder:
- the browser reads files locally,
- nothing is uploaded unless you explicitly host/share the export folder yourself.

Start here: {doc}`../../web_app/b_data_loading/index`

---

## Security model of the Python server (what can be exposed)

### Default: localhost-only

`CellucidServer` and notebook viewers bind to `127.0.0.1` by default.

That means:
- only the same machine can reach the server,
- people on your network cannot access it unless you deliberately bind to a public interface or tunnel it.

### If you bind to a non-localhost host (be careful)

If you run a server on:
- `0.0.0.0` (all interfaces), or
- a public IP / hostname,

then **anyone who can reach that host/port can fetch your dataset files**.

Practical implications:
- this is not an authenticated server by default,
- treat it like “public to whoever can reach it”.

Safe sharing patterns:
- SSH port forwarding (keeps the server bound to localhost on the remote machine)
- VPN + firewall rules
- hosting behind an authenticated reverse proxy (advanced)

---

## The `Host` rule (why a foreign host name is refused)

Binding to `127.0.0.1` stops other machines. It does **not** stop a web page you
happen to have open in another tab.

Any site can publish a DNS record for its own name with a very short lifetime
and then re-point that name at `127.0.0.1`. This is **DNS rebinding**. After the
name flips, the browser believes that site and your Cellucid server are the same
origin: it sends no `Origin` header, applies no CORS check, and hands the page
whatever the server returns — which is your dataset.

So the bind address is not the last line of defence. The `Host` header is:

- Cellucid reads the `Host` on every request, before any route, file read or
  event delivery runs.
- It accepts exactly one well-formed `Host` naming `localhost` or a loopback IP
  literal, on the port the server actually bound. A literal is safe because DNS
  cannot rebind one, and RFC 6761 reserves the name `localhost` for loopback.
- Anything else gets **421 Misdirected Request**, with a fixed body that never
  echoes the attacker’s value. The refused header is written to your server log
  only.
- The rule is applied whenever the bound address is loopback. A deliberate
  non-loopback bind (`--host 0.0.0.0`) is an explicit request for network
  exposure whose legitimate names Cellucid cannot know, so it is not enforced
  there by default.

### Reverse proxies must be declared

A reverse proxy — `jupyter-server-proxy` above all — forwards the browser’s
`Host` verbatim, so your loopback server receives the proxy’s host name and
refuses it. That request is indistinguishable from a rebound one: both arrive on
loopback naming a foreign authority, and a rebound page is same-origin, so it
can forge `X-Forwarded-For`, `X-Forwarded-Host` and friends at will. **There is
no header Cellucid could trust to tell them apart**, which is why the proxy is
never auto-detected — auto-detection would simply re-open the hole.

Name it instead:

```bash
cellucid serve /path/to/data --allowed-host hub.example.org
```

```python
from cellucid import serve

serve("/path/to/data", allowed_hosts=["hub.example.org"])
```

Rules for an entry, all enforced when the server is constructed:

- one bare host name — a DNS name or an IP address literal,
- no port, no scheme, no path, no credentials, no brackets, no zone id,
- **no wildcards**: a wildcard would re-admit rebinding under every name it
  covers, so each name is written out in full,
- matched case-insensitively on **any** port, because the proxy’s front-end port
  is unrelated to the loopback port bound here,
- anything unsupported is refused by name with an actionable error, never
  silently dropped.

Declaring names is also a statement that you know the complete set of
authorities this server should answer to, so it turns the `Host` check on for a
non-loopback bind as well. Loopback authorities on the bound port keep working
either way, so a local request to the exact URL Cellucid printed is never broken
by the declaration.

Leaving `allowed_hosts` unset — the default — leaves the loopback-only rule
exactly as described above.

---

## CORS and why it exists

The Python server sets CORS headers so the viewer can load data in different deployment modes.

Two common scenarios:
- Viewer UI served from `https://www.cellucid.com` but data served from your server (cross-origin).
- Viewer UI served from your local server (same-origin; still safe to keep consistent headers).

If you see browser errors mentioning CORS, it usually means:
- the browser could not reach the server at all (network/proxy),
- or you are mixing HTTPS and HTTP in a way the browser blocks (mixed content).

Start debugging at: {doc}`08_debugging_mental_model_where_to_look`.

---

## Viewer-source access and the local generation

### Why the viewer UI is established at startup

In notebook/server modes, Cellucid serves the web app UI from your Python server so:
- the viewer and dataset share the same origin (avoids mixed-content),
- your notebook embed is more reliable across environments.

Before binding its HTTP server, Cellucid:

1. fetches `cellucid-web-assets.json` from the configured source,
2. downloads every declared object into a staging directory,
3. verifies the complete file set, MIME types, byte lengths, SHA-256 hashes,
   and build identity, and
4. atomically publishes that generation locally.

### What this means without source access

- `cellucid serve`, `show(...)`, and `show_anndata(...)` require the configured
  source at each startup when they serve the viewer UI.
- Source or inventory failure stops startup. A previous local generation is
  kept intact but is not served as a substitute.
- `cellucid serve ... --no-web-ui` is the explicit data-endpoint-only mode and
  does not require viewer-source access.

### Where the cache lives

The default cache is a process-independent temporary directory
(platform-dependent). Select another location explicitly with
`--web-cache-dir PATH` or the `web_cache_dir=` Python argument.

In Python:

```python
from cellucid import get_web_cache_dir
get_web_cache_dir()
```

### Establish or verify the UI generation intentionally

Establish the complete configured source generation:

```python
from cellucid import AnnDataViewer

viewer = AnnDataViewer(
    "data.h5ad",
    auto_open=False,
    dataset_name="My study",
    dataset_id="my-study-v1",
)
viewer.ensure_web_ui_cached(force=True, show_progress=True)
viewer.stop()
```

Verify the selected existing generation without network access or mutation:

```python
viewer.ensure_web_ui_cached(force=False, show_progress=True)
```

Or clear the selected generation:

```python
from cellucid import clear_web_cache
clear_web_cache()
```

```{important}
`force=False` is verification-only. Viewer/server startup always establishes
the configured source generation and reports source failure directly.
```

---

## Privacy implications of session bundles

A `.cellucid-session` file is not “just UI state” — it can encode:
- highlight membership lists (which cells you selected),
- user-defined labels and categories,
- potentially other feature-specific state as the format evolves.

Treat session bundles as data artifacts:
- don’t share them publicly if the underlying dataset is sensitive,
- version them alongside the matching dataset export folder.

---

## Troubleshooting (privacy/network/security)

### Symptom: viewer startup cannot establish the UI generation

Likely causes:
- the configured source is unreachable or blocked,
- an inventory or declared object fails exact verification, or
- the selected cache directory is not writable.

Fix:
- allow the Python runtime to reach the configured `web_source_url`,
- correct the exact response error reported by Cellucid, and
- pass a writable `web_cache_dir` or `--web-cache-dir`.

### Symptom: every request returns `421 Misdirected Request`

Likely cause:
- the browser reaches Cellucid through a reverse proxy (usually
  `jupyter-server-proxy`), which forwards its own `Host`, and that name has not
  been declared.

Confirm:
- the server log contains
  `Refused a request whose Host header does not name this server`.

Fix:
- pass the proxy’s host name as `allowed_hosts=["hub.example.org"]` or
  `--allowed-host hub.example.org`, and keep passing the proxy’s browser base as
  `client_server_url=`.

See {doc}`08_debugging_mental_model_where_to_look` for the full remote-notebook
checklist.

### Symptom: “Someone else on my network can see my dataset”

Likely causes:
- you bound the server to `0.0.0.0` or a public interface,
- your machine firewall allows inbound traffic to the port.

Fix:
- bind to `127.0.0.1`,
- use SSH tunneling for remote access,
- add firewall rules.

---

## Next steps

- Performance/scaling: {doc}`07_performance_mental_model_and_scaling`
- Debugging checklist: {doc}`08_debugging_mental_model_where_to_look`
