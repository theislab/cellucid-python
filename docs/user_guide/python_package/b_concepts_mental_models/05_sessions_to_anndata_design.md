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
1) Python sends `{"type": "requestSessionBundle", "requestId": ...}` to the viewer on its channel (`/_cellucid/ws`).
2) The web app serializes the current session into raw `.cellucid-session` bytes.
3) The web app sends those bytes back on the same channel as one binary frame
   (header `{"type": "session_bundle", "requestId": ...}`, then the bundle).
4) The server streams the frame to a temp file (bounded memory) and validates it.
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
- **User-defined categorical and continuous fields** → cell-aligned `.obs` columns

Plus optional debugging/provenance under:
- `adata.uns["cellucid"]["session"]`

---

## Why it’s implemented this way

### “Bundle first”: one artifact everywhere

Everything is built around the `.cellucid-session` bundle:
- the web app's **Save** and **Load** buttons (**Session file**)
- Jupyter no-download capture
- Python application to AnnData

Benefits:
- reproducibility: the session is a sharable file artifact,
- debuggability: chunks can be inspected independently,
- modularity: features can contribute chunks without a monolithic schema.

### Treat session bundles as untrusted input

Even if the session “came from your own browser”, treat it as external input:
- validate framing (MAGIC header),
- bound every decompression by the size its chunk declares (no gzip “zip
  bomb” can inflate past it),
- bounds-check indices and code arrays,
- validate each declared chunk's exact metadata, dependency order, byte
  framing, codec, and decoded size before consuming it.

This is critical for:
- robustness (sessions should not brick restore),
- safety (avoid accidental data corruption when applying).

---

## Protocol details (what happens on the wire)

### Frontend → Python events (hooks)

Everything travels on the viewer channel: one WebSocket the viewer opens to
`/_cellucid/ws` and authenticates once with its `viewerId` and `viewerToken`
({doc}`../e_jupyter_hooks/05_architecture_viewer_channel`).

The viewer sends UI events as:
- JSON text frames (ready, hover, click, data mapping, command errors)
- binary frames for selections (the selected cell indices as little-endian
  `uint32`, strictly ascending) and session bundles

The channel belongs to one viewer, so events carry no `viewerId`.

### Python → frontend commands

Python sends commands as JSON text frames on the same channel, and a command
that names cells (a highlight, a visibility change) as one binary frame of
`uint32` cell indices. A command sent before the viewer has connected waits
and is sent, in order, when it does.

### Session bundle frames

Guardrails:
- Python registers an expected `(viewerId, requestId)` pair for as long as
  `get_session_bundle(timeout=...)` waits.
- The server accepts exactly one bundle for that pair (consumes the pending
  request); a late one is dropped.
- The frame is streamed to disk, so the temporary directory's free space is
  the only size bound; an invalid bundle makes `get_session_bundle()` raise.

---

## File format (high level)

A `.cellucid-session` file is:

1) MAGIC header: `CELLUCID_SESSION/2\n` (session format 2; any other signature is refused)  
2) `manifestByteLength` (u32 LE)  
3) `manifest` JSON bytes (UTF‑8)  
4) repeated chunks:
   - `chunkByteLength` (u32 LE)
   - `chunkBytes` (codec applied, often gzip)

The manifest contains:
- dataset fingerprint (id, display name, counts, cell order, data mapping)
- the dataset source (sample, link, server, or local file name)
- `requires`: the obs columns and genes the session uses
- chunk IDs
- chunk metadata (`kind`, `codec`, size guards)

In Python, you can inspect:

```python
from cellucid import CellucidSessionBundle

bundle = CellucidSessionBundle("my.cellucid-session")
bundle.dataset_fingerprint
bundle.required_fields
bundle.list_chunk_ids()[:20]
```

---

## Identity and mismatch guardrails (critical constraint)

### Cell identity is index-based (today)

Many dataset-dependent chunks align to **row position**:
- highlight memberships are cell positions (sorted-delta varints, a bitset or
  runs, whichever is smallest), or, for a whole category or most of one, a
  reference to that category's cells checked by count and digest
- user-defined codes are arrays aligned to dataset order (raw, runs, or packed
  1, 2 or 4 bits per code)

This is why exact fingerprint rejection exists.

### Fingerprint rejection is exact

The caller must supply `expected_dataset_id`. Python compares it with the
bundle's `datasetId`, compares `cellCount` with the target AnnData, requires
every key the `dataMapping` names to exist on it, compares `varCount` with
the gene count of the axis that mapping reads, and requires every obs column and
gene listed in the manifest's `requires`. Any mismatch raises `ValueError`
before any mutation is published.

This prevents the worst failure mode:
> applying “correct-looking” highlights to the wrong biological cells.

### `cellOrder` and `dataMapping`: recorded by the viewer, checked here

The manifest fingerprint carries two records beside those three scalars, and both
are required. `cellOrder` is exactly `embedding`, `dimension` (1, 2, or 3), and a
16-character lowercase hexadecimal `digest` of the coordinates of that
embedding. `embedding` is `{"key", "columns"}`: an `obsm` key with the 0-based
columns that were read, one per dimension. It exists because the three scalars
all survive a rewrite of the same dataset from re-sorted input, and index-based
selections would then silently denote different cells.

`dataMapping` records which stored arrays the session was showing: `kind`
(`"anndata"`), `embeddings` (`enabled`, the embeddings that are turned on, and
`defaults`, the `1d`/`2d`/`3d` default of each dimension), `knn`, `geneNames`
(`"index"` or `"column:<name>"`), `expression` (`"X"`, `"raw"`,
`"layer:<name>"`, or null), `log1p`, `highlyVariable`, `latent`, and
`vectorFields`, plus `modality` for a dataset of several modalities (a MuData).

The two readers check them differently:

| | web app | Python |
|---|---|---|
| `cellOrder` and `dataMapping` present | required | required |
| record shapes | validated | validated |
| `dataMapping` keys exist on the dataset | checked | checked against the `AnnData` |
| `cellOrder` digest value | compared against the loaded embedding | never checked |

Python cannot check the digest. It is taken over the normalized coordinates the
viewer drew, while this API applies a session to an `AnnData` whose row order
the caller controls. Enforcing cell-order identity is the viewer's job.

So a bundle this reader accepts may still be refused by the web app. That
direction is expected.

## Troubleshooting (protocol-level)

### Symptom: “Dropped a session bundle nobody is waiting for” (log)

Meaning:
- the bundle arrived after its `get_session_bundle()` call stopped waiting
  (its timeout passed).

Fix:
- request a new bundle again from the *current* viewer instance, with a
  longer `timeout=`.

### Symptom: “This file is not a Cellucid session file of format 2 …”

Meaning:
- the bytes the viewer sent on its channel for a `requestSessionBundle` are
  not a `.cellucid-session` file of the current format (or a file passed to
  `CellucidSessionBundle(...)` is not one). On the channel, the frame is
  dropped and `viewer.get_session_bundle()` raises `RuntimeError` ("The viewer
  did not capture the session bundle: The viewer sent an invalid session
  bundle: …") instead of waiting for its timeout.

Fix:
- make sure the viewer page and the Python package are current versions of
  the same release (reload the viewer after upgrading), then request a new
  bundle.

---

## Related docs

- User-facing workflow: {doc}`05_sessions_to_anndata_bridge`
- Web app session model: {doc}`../../web_app/l_sessions_sharing/01_session_mental_model`
