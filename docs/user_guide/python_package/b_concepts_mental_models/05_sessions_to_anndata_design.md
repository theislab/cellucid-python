# Sessions → AnnData (Design notes)

**Audience:** developer / power user  
**Time:** 20–30 minutes  
**Goal:** understand the protocol and guardrails behind `viewer.get_session_bundle()` and `apply_cellucid_session_to_anndata(...)`.

```{note}
This page is a “deep path” companion to {doc}`05_sessions_to_anndata_bridge`.
Most users should start there.
```

---

## What’s implemented (today)

### 1) Jupyter: “no download” session capture

In a notebook, Python can request the current session bundle from the frontend and receive it as a `CellucidSessionBundle` handle.

High-level flow:
1) Python sends `{"type": "requestSessionBundle", "requestId": ...}` into the viewer iframe (postMessage).
2) The web app serializes the current session into raw `.cellucid-session` bytes.
3) The web app uploads those bytes back to the local Python server:
   - `POST /_cellucid/session_bundle?viewerId=...&requestId=...`
4) The server validates + streams the upload to a temp file (bounded memory).
5) The server routes a `type="session_bundle"` event back to Python so `viewer.get_session_bundle()` can return.

Why this matters:
- Users never have to find/download a file from the browser.
- The artifact is the same `.cellucid-session` format used by the web app.

### 2) Python: streaming-friendly `.cellucid-session` reader

`CellucidSessionBundle`:
- indexes chunk offsets once,
- reads/decompresses chunk payloads on demand.

This keeps memory usage predictable even when sessions include large binary chunks.

### 3) AnnData mutation (safe, minimal, DRY)

`apply_cellucid_session_to_anndata(...)` currently materializes:
- **Highlights** → one boolean `.obs` column per highlight group
- **User-defined categorical fields** → `.obs`/`.var` categorical columns

Plus optional debugging/provenance under:
- `adata.uns["cellucid"]["session"]`

---

## Why it’s implemented this way

### “Bundle first”: one artifact everywhere

Everything is built around the `.cellucid-session` bundle:
- web app Save/Load State
- Jupyter no-download capture
- Python application to AnnData

Benefits:
- reproducibility: the session is a sharable file artifact,
- debuggability: chunks can be inspected independently,
- modularity: features can contribute chunks without a monolithic schema.

### Treat session bundles as untrusted input

Even if the session “came from your own browser”, treat it as external input:
- validate framing (MAGIC header),
- enforce size limits (including gzip “zip bomb” guards),
- bounds-check indices and code arrays,
- skip unknown chunk contributors rather than failing the whole load.

This is critical for:
- robustness (sessions should not brick restore),
- safety (avoid accidental data corruption when applying).

---

## Protocol details (what happens on the wire)

### Frontend → Python events (hooks)

The viewer sends UI events via:
- `POST /_cellucid/events`

The payload includes `viewerId` so the server can route events to the right viewer object.

### Python → frontend commands

Python sends commands via `postMessage` to the iframe. Messages are tagged with:
- `viewerId`
- `viewerToken` (used to authenticate messages even when `targetOrigin='*'` is required by notebook proxies)

### Session bundle upload

Bundle capture uses a gated endpoint:
- `POST /_cellucid/session_bundle?viewerId=...&requestId=...`

Guardrails:
- Python registers an expected `(viewerId, requestId)` pair with a TTL.
- The server accepts exactly one upload for that pair (consumes the pending request).
- A hard size cap exists to prevent accidental huge uploads.

---

## File format (high level)

A `.cellucid-session` file is:

1) MAGIC header: `CELLUCID_SESSION\n`  
2) `manifestByteLength` (u32 LE)  
3) `manifest` JSON bytes (UTF‑8)  
4) repeated chunks:
   - `chunkByteLength` (u32 LE)
   - `chunkBytes` (codec applied, often gzip)

The manifest contains:
- dataset fingerprint
- chunk IDs
- chunk metadata (`kind`, `codec`, size guards)

In Python, you can inspect:

```python
from cellucid import CellucidSessionBundle

bundle = CellucidSessionBundle("my.cellucid-session")
bundle.dataset_fingerprint
bundle.list_chunk_ids()[:20]
```

---

## Identity and mismatch guardrails (critical constraint)

### Cell identity is index-based (today)

Many dataset-dependent chunks align to **row position**:
- highlight membership lists are indices
- user-defined codes are arrays aligned to dataset order

This is why the mismatch policy exists.

### Mismatch policy is intentionally conservative

On mismatch, Python can:
- error (strict pipelines),
- warn + skip (safe default),
- skip silently (rarely recommended; use only when you know what you’re doing).

This prevents the worst failure mode:
> applying “correct-looking” highlights to the wrong biological cells.

---

## Roadmap (what “complete” might look like later)

These are logical extensions (not all implemented yet):
- stable cell IDs (map by `obs_names` so sessions survive reorder/subset),
- apply more field kinds (continuous, booleans, numeric),
- apply rename/delete overlays (opt-in; potentially destructive),
- richer apply reports (per-field diagnostics, partial-apply modes),
- better tooling for inspecting session contents (CLI/summarizers).

---

## Troubleshooting (protocol-level)

### Symptom: “No pending session bundle request”

Meaning:
- the upload arrived without a matching, currently valid requestId (stale tab, TTL expired, wrong viewerId).

Fix:
- request a new bundle again from the *current* viewer instance.

### Symptom: “Invalid session bundle MAGIC header”

Meaning:
- the upload is not a `.cellucid-session` file (corrupted bytes, wrong endpoint).

Fix:
- retry capture; if persistent, check browser console errors and verify you are not behind a proxy that alters request bodies.

---

## Related docs

- User-facing workflow: {doc}`05_sessions_to_anndata_bridge`
- Web app session model: {doc}`../../web_app/l_sessions_sharing/01_session_mental_model`

