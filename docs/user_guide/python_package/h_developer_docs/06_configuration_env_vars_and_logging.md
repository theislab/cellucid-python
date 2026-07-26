# Configuration, environment variables, and logging

This page lists the exact **configuration knobs** available in
`cellucid-python`: CLI flags, function arguments, defaults, and logging
behavior.

---

## Configuration layers (highest priority wins)

In general, configuration flows like this:

1) **Function arguments** (e.g. `serve(..., host=..., port=...)`)
2) **CLI flags** (e.g. `cellucid serve ... --host ... --port ...`)
3) **Defaults** in code (`DEFAULT_HOST`, `DEFAULT_PORT`)

---

## Environment variables

Cellucid has no environment-variable configuration layer. Values such as host,
port, web source, cache directory, and browser-facing URL must come from the
documented CLI flags or function/class arguments.

`VSCODE_PID` and `JPY_PARENT_PID` may be read only to label the current
notebook context in diagnostics; they do not configure routing. Variables such
as `CELLUCID_CLIENT_SERVER_URL` are ignored.

---

## Defaults

Defined in `cellucid-python/src/cellucid/_server_base.py`:

- `DEFAULT_HOST = "127.0.0.1"` (local only by default)
- `DEFAULT_PORT = 8765`

For remote access, explicitly bind:
- `--host 0.0.0.0` (and use SSH tunneling or a firewall)

---

## Viewer and notebook URL arguments

### `web_cache_dir`

Selects the directory where one verified viewer generation is published.
Configure it explicitly:

```bash
cellucid serve /path/to/data --web-cache-dir /path/to/cache
```

Python APIs accept `web_cache_dir=...`. Use `cellucid.get_web_cache_dir()` to
inspect the default temporary location and `cellucid.clear_web_cache()` to
remove its current generation.

### `client_server_url`

Overrides the URL that the **browser** should use to reach the Python server.

This is primarily for remote notebook environments where:
- the kernel runs on a remote VM,
- and `http://127.0.0.1:<port>` is not reachable from your browser.

If you have an HTTPS reverse proxy for the server (or a stable tunnel URL),
pass the exact browser-reachable base URL:

```python
from cellucid import show_anndata
viewer = show_anndata(
    adata,
    dataset_name="Example",
    dataset_id="example",
    client_server_url="https://notebooks.example/proxy/8765",
)
```

The value must be an absolute HTTP or HTTPS URL with no credentials, query,
fragment, trailing slash, or surrounding whitespace.

See details in: {doc}`10_jupyter_embedding_architecture`

---

## Logging

### CLI logging

The CLI uses Python logging with two user-facing knobs:

- `--verbose`: `DEBUG` logs (more detail, stack traces)
- `--quiet`: suppresses most informational prints

### Library logging

When you import `cellucid` as a library, logging is “normal Python logging”:
- by default you may see little/no output,
- configure logging in your application/notebook if you want debug output.

Example:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

---

## CORS / origin policy (important for security)

Servers add CORS headers that allow:
- `localhost` / `127.0.0.1`
- the canonical hosted viewer origin (`https://www.cellucid.com`)
- `https://cellucid.com` and `https://*.cellucid.com`

This is implemented in:
- `cellucid-python/src/cellucid/_server_base.py` (`CORSMixin._get_allowed_origin`)

If you run a server on a public host, do not assume CORS is a full security boundary.
Treat the server as trusted/local unless you have a hard threat model and additional protections (see {doc}`17_security_privacy_and_networking`).

---

## Troubleshooting

### Symptom: “Notebook says mixed content is blocked”

Cause:
- the notebook is served over HTTPS,
- but the viewer tries to load `http://127.0.0.1:<port>`.

Fix options (ordered):
1) Use `jupyter-server-proxy` (recommended).
2) Use a tunnel/reverse proxy and pass its exact base as
   `client_server_url=`.

See: {doc}`10_jupyter_embedding_architecture`

### Symptom: viewer startup cannot establish the UI generation

Cause:
- the source inventory or a declared object could not be fetched or verified,
- or the selected generation directory cannot be written.

Fix:
- confirm the configured `web_source_url` is reachable,
- pass a writable `web_cache_dir`, and
- use `cellucid.clear_web_cache()` only when you intentionally want to remove
  the selected generation before establishing it again.
