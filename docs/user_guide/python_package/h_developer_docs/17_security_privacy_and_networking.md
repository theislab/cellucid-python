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
- anyone who learns a browser-tab viewer's URL (it carries the viewer's token) can
  open that viewer's channel and trigger your hook callbacks,
- and you may leak private metadata via URLs/logs.

---

## Network surfaces

### 1) The dataset server (`AnnDataServer`)

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
  exact host names, which is required behind a reverse proxy (a hub's proxy
  route) because the proxy forwards the browser's `Host`
  verbatim and a proxied request is indistinguishable from a rebound one
- one bare DNS name or IP literal per entry: no port, no scheme, no wildcard;
  matched on any port
- declaring names also applies the `Host` check to a non-loopback bind
- see {doc}`../b_concepts_mental_models/06_privacy_security_and_offline_vs_online`

### 2) Exact viewer-generation delivery

At every start, before it binds, the server downloads the viewer UI declared
by the inventory `cellucid-web-assets.json` from the web source:
- `https://www.cellucid.com` (default; `--web-source-url` / `web_source_url=`)

It byte-verifies the files and keeps them locally (`--no-web-ui` /
`serve_web_ui=False` skips this and serves the data routes only).

Privacy note:
- the cached assets are UI code, not your dataset
- but the cache directory location may appear in logs/debug reports; treat it as potentially sensitive.

Config:
- `--web-cache-dir` or `web_cache_dir=`

### 3) The viewer channel

- `GET /_cellucid/ws`, a WebSocket per notebook viewer, carrying Python's
  commands and the viewer's events (selections and session bundles as binary
  frames)

Important:
- the upgrade passes the `Host` check and an `Origin` check (the server's own
  host, a loopback name, or a host in `allowed_hosts`); WebSockets are not
  covered by CORS, so this is the cross-site check;
- the first frame must name a registered viewer and its exact `viewerToken`,
  compared with `hmac.compare_digest`; anything else closes the channel
  (`4401`) before a frame is delivered.

Implication:
- the token is generated per viewer. A browser-tab viewer carries it in its
  URL, so anyone who can read that URL and reach the server from an accepted
  origin can open the channel and trigger your hook callbacks, which run the
  Python you wrote in your handlers. A notebook widget never puts it in an
  address, a widget trait or a log: it posts it to the app's exact origin in a
  message.

### 4) Wire capability endpoint

- `GET /_cellucid/protocol`
- used by the web build to find the viewer channel and learn which events this
  `cellucid` accepts before emitting one

No credentials:
- the body is the same for every viewer, every dataset and every request — it
  describes the installed `cellucid` version's protocol, not your data — so
  there is nothing here for a token to protect, and it creates no state.

Implication:
- reaching this route tells someone which `cellucid` protocol generation you
  are running. So does `/_cellucid/health`, which already returns the version.

### 5) Session bundles

- arrive on the viewer channel as one binary frame, for notebook “no
  download” session capture

Guards:
- only on an authenticated channel, and only for a request Python registered
  and is still waiting for (`viewerId`, `requestId`)
- validates a magic header and the complete bundle structure
- streams to disk, so memory stays bounded; the temporary directory's free
  space is the only size bound

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
- `cellucid-python/src/cellucid/anndata_session/`

If you extend session parsing, keep the “untrusted input” stance.

---

## Troubleshooting

### Symptom: “CORS blocked” in browser console

Likely causes:
- you are using a non-allowed origin (custom domain),
- a proxy in front of a browser-tab server is rewriting origins,
- you are trying to load the viewer from a file:// URL.

Fix:
- use the server’s `/` URL (same origin),
- in a notebook, use `show_anndata()` (no server, no CORS),
- or serve the viewer from an allowed origin.

### Symptom: “I need to share a server with collaborators”

Recommended pattern:
- run the server on a trusted host,
- use SSH tunneling (`ssh -L <port>:localhost:<port> user@host`),
- and share instructions rather than exposing the port publicly.

If you truly need public hosting, treat it as a separate deployment project with an explicit security review.
