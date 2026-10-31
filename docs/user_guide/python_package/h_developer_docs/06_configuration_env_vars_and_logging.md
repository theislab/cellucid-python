# Configuration, environment variables, and logging

This page lists the exact **configuration knobs** available in
`cellucid-python`: CLI flags, function arguments, defaults, and logging
behavior.

---

## Configuration layers (highest priority wins)

There are two layers, and the first given wins:

1) **Function arguments** (e.g. `serve_anndata(..., host=..., port=...)`), or
   the **CLI flags** that `cellucid serve` passes to the same function
   (`--host ... --port ...`)
2) **Defaults** in code (`DEFAULT_HOST`, `DEFAULT_PORT`)

For a dataset's launch choices (embeddings, graph, expression, modality, ...)
an explicit argument or flag wins over the store's own
`uns['cellucid']['defaults']`, which wins over the viewer's default rule.

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

A browser-tab viewer (`AnnDataViewer` with `port=None`) binds
`127.0.0.1:8765` — the one local port Google Drive sign-in accepts — or, when another
viewer holds it, the next free port, and prints a note saying so. A port given
explicitly is used or refused, never replaced. `serve_anndata`,
`AnnDataServer` and `cellucid serve` always bind exactly the port given.

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

Overrides the URL that the **browser** should use to reach a browser-tab
viewer's server ({class}`~cellucid.AnnDataViewer`), for a page opened through
a reverse proxy or a tunnel whose address differs from the server's own:

```python
from cellucid import AnnDataViewer
viewer = AnnDataViewer(
    "data.h5ad",
    dataset_name="Example",
    dataset_id="example",
    client_server_url="https://notebooks.example/proxy/8765",
    allowed_hosts=["notebooks.example"],
)
```

The value must be an absolute HTTP or HTTPS URL with no credentials, query,
fragment, trailing slash, or surrounding whitespace. A notebook viewer
(`show_anndata`, `cellucid.Widget`) has no server address and needs none.

## Logging

### CLI logging

The CLI uses Python logging with two user-facing knobs:

- `--verbose`: `DEBUG` logs (more detail, stack traces)
- `--quiet`: suppresses most informational prints

### Library logging

When you import `cellucid` as a library, logging is “normal Python logging”:
- by default you may see little/no output,
- configure logging in your application/notebook if you want debug output.

The package logs under these logger names:

| Logger | What it reports |
|---|---|
| `cellucid.server` | refused requests (a `Host` that does not name the server); clients that went away mid-response, at `DEBUG` |
| `cellucid.anndata_server` | per-request lines (at `DEBUG`), request failures, background-server failures |
| `cellucid.viewer_channel` | refused or closed viewer channels and failing notebook hooks |
| `cellucid.jupyter` | notebook viewer start, stop and channel registration |

Example:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

---

## CORS / origin policy (important for security)

Servers add CORS headers that allow:
- `localhost` / `127.0.0.1`
- `https://www.cellucid.com` and `https://cellucid.com` (the hosted app's origins, exactly)

This is implemented in:
- `cellucid-python/src/cellucid/_server_base.py` (`CORSMixin._get_allowed_origin`)

If you run a server on a public host, do not assume CORS is a full security boundary.
Treat the server as trusted/local unless you have a hard threat model and additional protections (see {doc}`17_security_privacy_and_networking`).

---

## Troubleshooting

### Symptom: “A browser-tab viewer cannot be reached (mixed content, remote kernel)”

Cause:
- the page is served over HTTPS, or the kernel runs on another machine,
- but the viewer's server is `http://127.0.0.1:<port>` on the kernel's machine.

Fix options (ordered):
1) In a notebook, use `show_anndata()`: a widget needs no address.
2) Use an SSH tunnel, or a reverse proxy whose exact base you pass as
   `client_server_url=` and whose host you pass as `allowed_hosts=[...]`
   (Cellucid does not select a proxy itself).

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
