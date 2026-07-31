# Security, privacy, and networking

This page documents the practical security and privacy considerations for `cellucid-python`.

It is written for:
- people serving sensitive datasets,
- people running Cellucid on shared servers,
- and developers modifying server/hook/session behavior.

---

## Threat model (practical)

Cellucid Python is primarily designed for:
- local use (`127.0.0.1`),
- trusted environments (your laptop, your lab server with SSH access),
- and notebook workflows where the browser and kernel are “yours”.

If you expose a Cellucid server to an untrusted network, you must assume:
- others can fetch dataset files,
- others can trigger event endpoints,
- and you may leak private metadata via URLs/logs.

---

## Network surfaces

### 1) Dataset servers (exported + AnnData)

Binding:
- default host is `127.0.0.1` (local only)
- binding to `0.0.0.0` exposes the server to your network

Recommendation:
- prefer SSH tunneling instead of exposing ports directly.

`Host` validation (defeats DNS rebinding):
- every request is checked before any route, file read or event delivery
- accepted: exactly one well-formed `Host` naming `localhost` or a loopback IP
  literal, on the port actually bound
- refused: everything else, with `421 Misdirected Request` and a body that never
  echoes the value; the header goes to the server log only
- enforced whenever the bound address is loopback

Config:
- `--allowed-host HOST` (repeatable) or `allowed_hosts=[...]` declares extra
  exact host names, which is required behind a reverse proxy such as
  `jupyter-server-proxy` because the proxy forwards the browser's `Host`
  verbatim and a proxied request is indistinguishable from a rebound one
- one bare DNS name or IP literal per entry: no port, no scheme, no wildcard;
  matched on any port
- declaring names also applies the `Host` check to a non-loopback bind
- see {doc}`../b_concepts_mental_models/06_privacy_security_and_offline_vs_online`

### 2) Exact viewer-generation delivery

The server may download the viewer UI from:
- `https://www.cellucid.com`

It caches the files locally.

Privacy note:
- the cached assets are UI code, not your dataset
- but the cache directory location may appear in logs/debug reports; treat it as potentially sensitive.

Config:
- `--web-cache-dir` or `web_cache_dir=`

### 3) Hooks event endpoint

- `POST /_cellucid/events`
- used for frontend → Python hooks

Important:
- events are routed by `viewerId` but are not authenticated beyond that.

Implication:
- if others can reach your server, they can trigger your hook callbacks (which may run arbitrary Python code you wrote in handlers).

### 4) Session bundle upload endpoint

- `POST /_cellucid/session_bundle?viewerId=...&viewerToken=...&requestId=...`
- used for notebook “no download” session capture

Guards:
- requires a pending request registered by Python (`viewerId`, `requestId`)
- validates a magic header (`CELLUCID_SESSION\\n`)
- enforces a max upload size (currently 512MB)
- streams to disk to avoid memory blowups

Still:
- do not expose it publicly for sensitive workflows.

---

## CORS policy (what it does and does not do)

The server sends CORS headers that allow:
- `localhost` / `127.0.0.1`
- `https://www.cellucid.com`
- `https://cellucid.com` and `https://*.cellucid.com`

This helps the browser decide whether JS from one origin can read responses from another.

It does not:
- stop direct HTTP clients,
- prevent same-origin attackers,
- or act as authentication.

---

## Privacy guidance for bug reports and screenshots

Before sharing logs/screenshots:

- replace private dataset names with an explicitly synthetic identifier such as
  `pbmc_demo`
- remove private hostnames and user paths
- remove tokens or secrets
- avoid posting raw `.cellucid-session` bundles from sensitive datasets

When adding screenshots to docs, follow:
- {doc}`15_docs_development_and_style_guide`

---

## Session bundles are untrusted input (developer note)

`.cellucid-session` bundles can come from users, browsers, and shared files.

Python-side parsing and application takes a defensive stance:
- bounds checks on indices/codes
- a required exact dataset ID and exact cell/variable-count fingerprint check
  that raises before mutation on any mismatch
- chunk size guards for gzip decompression

Related code:
- `cellucid-python/src/cellucid/session_bundle.py`
- `cellucid-python/src/cellucid/anndata_session.py`

If you extend session parsing, keep the “untrusted input” stance.

---

## Troubleshooting

### Symptom: “CORS blocked” in browser console

Likely causes:
- you are using a non-allowed origin (custom domain),
- your notebook proxy is rewriting origins,
- you are trying to load the viewer from a file:// URL.

Fix:
- use the server’s `/` URL (same origin),
- use notebook proxy mechanisms (jupyter-server-proxy),
- or serve the viewer from an allowed origin.

### Symptom: “I need to share a server with collaborators”

Recommended pattern:
- run the server on a trusted host,
- use SSH tunneling (`ssh -L <port>:localhost:<port> user@host`),
- and share instructions rather than exposing the port publicly.

If you truly need public hosting, treat it as a separate deployment project with an explicit security review.
