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
