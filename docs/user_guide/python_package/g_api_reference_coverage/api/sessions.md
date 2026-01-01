# Sessions (`.cellucid-session` bundles)

```{eval-rst}
.. currentmodule:: cellucid
```

Cellucid sessions let you capture *what you did in the viewer* (highlights, user-defined fields, etc.) into a portable artifact:

- File extension: **`.cellucid-session`**
- Use cases:
  - “I selected/annotated cells in the UI; now I want those labels back in Python.”
  - “I want to attach interactive exploration results to a notebook or paper.”
  - “I want a reproducible trail of UI-derived fields.”

This page documents:
- {class}`~cellucid.CellucidSessionBundle` (read/inspect bundles)
- {func}`~cellucid.apply_cellucid_session_to_anndata` (apply bundles to AnnData)

---

## Audience + prerequisites

**Audience**
- Wet lab / beginner: use the “Apply to AnnData” workflow and the troubleshooting section.
- Computational: pay attention to dataset mismatch policies and cell ordering assumptions.
- Developer: read the “Format + safety” section before handling bundles from untrusted sources.

**Prerequisites**
- If applying to AnnData: `pandas` and an `AnnData` object.
- A session bundle, obtained either from the web app UI or from a notebook viewer.

---

## Fast path (apply UI highlights back onto an AnnData)

```python
from cellucid import CellucidSessionBundle

bundle = CellucidSessionBundle("my-session.cellucid-session")

# Apply to AnnData (creates new columns in adata.obs / adata.var)
adata2 = bundle.apply_to_anndata(adata)
```

Now inspect:
- `adata2.obs` for boolean highlight columns and user-defined categorical fields.

---

## Practical path (common workflows)

### 1) How to obtain a session bundle

You have two typical options:

#### A) From the web app UI (download a `.cellucid-session` file)

Use this when you explore in the standalone browser viewer (notebook or not).

<!-- SCREENSHOT PLACEHOLDER
ID: session-export-ui
Where it appears: Sessions (.cellucid-session bundles) → How to obtain a session bundle
Capture:
  - UI location: the Cellucid web app panel/menu where a session export/download is triggered
  - State prerequisites: any non-trivial session state (e.g., at least one highlight group)
  - Action to reach state:
    - create a highlight group (select cells → add to highlight)
    - open the “Export / Session / Save” area
    - click the button that downloads a `.cellucid-session`
Crop:
  - Include: the button/menu item that clearly indicates “download session” + any file extension hint
  - Exclude: private dataset names, account avatars, repo URLs if private
Alt text:
  - Cellucid web app showing the control used to export or download a .cellucid-session bundle.
Caption:
  - Use the session export control to download a `.cellucid-session` bundle containing highlights and other session state.
-->
```{figure} ../../../../_static/screenshots/placeholder-screenshot.svg
:alt: Cellucid web app showing the control used to export or download a .cellucid-session bundle.
:width: 100%

Use the session export control to download a `.cellucid-session` bundle containing highlights and other session state.
```

#### B) From a notebook viewer without downloading (advanced)

Use this when you are running Cellucid inside a notebook and want Python to receive the bundle directly:

```python
viewer.wait_for_ready(timeout=30)
bundle = viewer.get_session_bundle(timeout=60)
```

### 2) Inspect what’s inside a bundle

```python
bundle = CellucidSessionBundle("my-session.cellucid-session")

print(bundle.dataset_fingerprint)
print(bundle.list_chunk_ids()[:10])

# Read decoded payloads (JSON chunks become Python objects)
meta = bundle.decode_chunk("highlights/meta")
```

### 3) Apply to AnnData (what gets created)

Applying a bundle can add:
- **Highlight membership columns** in `adata.obs` (typically boolean)
- **User-defined categorical fields** in `adata.obs` or `adata.var` (depending on source)
- Optional bookkeeping in `adata.uns["cellucid"]` (manifest, fingerprint, and chunk snapshots)

You can choose behavior via policies (see API reference below):
- dataset mismatch handling (`dataset_mismatch=...`)
- column conflict handling (`column_conflict=...`)
- include/exclude deleted user-defined fields

---

## Deep path (format + safety)

### Treat bundles as untrusted input

Session bundles may come from:
- collaborators,
- community annotation workflows,
- or downloaded files of unknown origin.

The bundle reader and AnnData applicator include guards:
- header/magic validation,
- manifest size limits,
- chunk decompression limits (zip-bomb protection),
- bounds checks for indices/codes.

### Format overview (for developers)

The `.cellucid-session` framing is:
1. MAGIC bytes: `CELLUCID_SESSION\n`
2. `manifestByteLength` (u32 little-endian)
3. UTF-8 manifest JSON bytes
4. repeated chunks: `[chunkByteLength (u32 LE), chunkBytes...]`

Each chunk has metadata in the manifest (`kind`, `codec`, expected sizes, ids).

---

## API reference

### `CellucidSessionBundle`

```{eval-rst}
.. autoclass:: CellucidSessionBundle
   :members:
   :show-inheritance:
```

### `apply_cellucid_session_to_anndata`

```{eval-rst}
.. autofunction:: apply_cellucid_session_to_anndata
```

---

## Edge cases (do not skip)

### Dataset mismatch (most common)

A bundle is only meaningfully applicable if it matches the dataset it was created from.
Common mismatch causes:
- you filtered/reordered cells before exporting,
- you regenerated the AnnData with different preprocessing,
- you are applying a session to a subset.

If mismatch is detected, default behavior is to **warn and skip dataset-dependent chunks**.

### Cell ordering assumptions

Highlights and codes are stored as **indices** into the dataset.
If your AnnData row order differs, labels will be assigned to the wrong cells.

### Conflicting column names

Applying a bundle can create columns like `cellucid_highlight__<group_id>`.
If those already exist, choose a conflict policy (`error`, `overwrite`, `suffix`).

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Dataset fingerprint mismatch”
Likely causes:
- You’re applying the bundle to a different AnnData (or a differently ordered one).

Fix:
- Apply to the original AnnData used to create the session.
- If you know what you’re doing, choose `dataset_mismatch="skip"` or `"warn_skip"` to salvage non-dependent chunks.

### Symptom: “Not a .cellucid-session file (invalid MAGIC header)”
Fix:
- Confirm the file is actually a `.cellucid-session` bundle and not (for example) a `.json` export or a renamed file.

### Symptom: “Chunk decompressed too large”
Fix:
- Treat the file as potentially unsafe/corrupted.
- If it is trusted but legitimately large, the bundle writer should include `uncompressedBytes` in metadata; otherwise, re-export the session.

---

### Symptom: “No highlight columns were added”
Likely causes:
- The session contains no highlight data (no highlights were created).
- Dataset mismatch caused dataset-dependent chunks to be skipped.
- You disabled highlight application (`add_highlights=False`).

How to confirm:
- Inspect `bundle.list_chunk_ids()` and look for `highlights/meta`.

Fix:
- Create at least one highlight group in the UI and re-export the session.
- If mismatch is expected, set `dataset_mismatch="skip"` only if you understand the risk.

---

### Symptom: `KeyError: 'highlights/meta'` (or another chunk id)
What it means:
- That chunk does not exist in the bundle (session didn’t include that feature/state).

Fix:
- Treat missing chunks as “feature not used in that session”.
- Guard your code: check `chunk_id in bundle.list_chunk_ids()` before decoding.

---

### Symptom: “Unsupported session chunk codec: …”
Fix:
- Update to a newer `cellucid` version that supports the codec used by the bundle.
- If you control the writer, re-export using supported codecs (`none` or `gzip`).

---

### Symptom: “Column already exists: …”
Likely causes:
- You applied the same session twice, or you already have similarly named columns.

Fix:
- Use `column_conflict="suffix"` (default) to avoid overwriting.
- Or choose `column_conflict="overwrite"` if you intentionally want replacement.

---

### Symptom: “apply_cellucid_session_to_anndata requires pandas”
Fix:
- Install pandas in your environment: `pip install pandas`.

---

## See also

- {doc}`jupyter` for capturing a bundle in-notebook (`viewer.get_session_bundle`)
- {doc}`../../e_jupyter_hooks/index` for deeper hook/event documentation
