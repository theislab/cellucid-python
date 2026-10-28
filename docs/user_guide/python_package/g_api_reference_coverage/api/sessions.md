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
- Computational: pay attention to exact dataset identity and cell ordering assumptions.
- Developer: read the “Format + safety” section before handling bundles from untrusted sources.

**Prerequisites**
- If applying to AnnData: an `AnnData` object. Pandas is part of the Cellucid
  core installation.
- A session bundle, obtained either from the web app UI or from a notebook viewer.

---

## Fast path (apply UI highlights back onto an AnnData)

```python
from cellucid import CellucidSessionBundle

bundle = CellucidSessionBundle("my-session.cellucid-session")

# Apply to AnnData (adds columns to obs and records the session in uns)
adata2 = bundle.apply_to_anndata(
    adata,
    expected_dataset_id="my-study-v1",
)
```

Now inspect:
- `adata2.obs` for boolean highlight columns and user-defined categorical fields.

---

## Practical path (common workflows)

### 1) How to obtain a session bundle

You have two typical options:

#### A) From the web app UI (download a `.cellucid-session` file)

Use this when you explore in the standalone browser viewer (notebook or not).

```{figure} ../../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: The Session panel with the Pancreatic endocrinogenesis sample loaded: the Dataset box with its name, cells and genes, format and host, description and a collapsed Details; Open data with the Sample dataset dropdown, the Open… button and the Link or server bar; Save and Load under Session file; and the collapsed Data mapping group.
:width: 532px

The Session panel gathers every loading path under **Open data** and keeps the **Save** and **Load** buttons under **Session file**, directly below it.
```

#### B) From a notebook viewer without downloading (advanced)

Use this when you are running Cellucid inside a notebook and want Python to receive the bundle directly. A notebook widget hears the viewer only between cells, so ask in one cell:

```python
request = viewer.request_session_bundle()
```

and collect it in the next:

```python
bundle = request.result()
```

The same two steps work in a browser-tab viewer (`AnnDataViewer`), which
answers while a cell runs, so there `result(timeout=60)` may follow in the
same cell. To bring results onto an
AnnData without a bundle, `viewer.pull(adata)` ({doc}`round_trip`).

### 2) Inspect what’s inside a bundle

```python
bundle = CellucidSessionBundle("my-session.cellucid-session")

print(bundle.dataset_fingerprint)   # dataset id, name, counts, cell order, data mapping
print(bundle.dataset_source)        # where it was opened from, or None
print(bundle.required_fields)       # {"obs": [...], "var": [...]} the session uses
print(bundle.list_chunk_ids()[:10])

# Read decoded payloads (JSON chunks become Python objects)
meta = bundle.decode_chunk("highlights/meta")

# Analysis windows are typed JSON: result columns come back as numpy arrays
from cellucid.session_codecs import decode_typed_json
windows = decode_typed_json(bundle.decode_chunk("analysis/windows"))
```

### 3) Apply to AnnData (what gets created)

Applying a bundle can add:
- **Highlight membership columns** in `adata.obs` (typically boolean)
- **User-defined categorical and continuous cell-aligned fields** in
  `adata.obs`
- Optional bookkeeping in `adata.uns["cellucid"]["session"]` (`store_uns=True`,
  the default): the manifest, the dataset fingerprint, the added columns, and
  the `highlights/meta` and `core/field-overlays` chunks

The target may also be a MuData: the session is checked against what the
viewer showed (the cells, and the modality its data mapping names), and the
columns are added to the global `mdata.obs`.

Application requires the exact `expected_dataset_id`. Identity, shape, bounds,
every field the session uses (`requires`: obs columns and genes, all missing
ones named at once), and target-column conflicts are rejected before mutation. You may include or
exclude deleted user-defined fields explicitly.

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

The `.cellucid-session` framing (session format 2) is:
1. MAGIC bytes: `CELLUCID_SESSION/2\n` — the format version is part of the signature
2. `manifestByteLength` (u32 little-endian)
3. UTF-8 manifest JSON bytes
4. repeated chunks: `[chunkByteLength (u32 LE), chunkBytes...]`

The manifest holds exactly `createdAt`, `datasetFingerprint` (`datasetId`,
`name`, `cellCount`, `varCount`, `cellOrder`, `dataMapping`), `source` (where
the dataset was opened from: a sample, a link, a server, or the name of a local
file or folder; `null` for a notebook), `requires` (the obs columns and genes
the session uses) and `chunks`. Each chunk has metadata in the manifest (`kind`,
`codec`, expected sizes, ids).

A file that does not start with `CELLUCID_SESSION/2\n` is refused with the web
app's sentence: *This file is not a Cellucid session file of format 2 (it does not start with CELLUCID_SESSION/2).*

### Compact payloads

| Chunk | Encoding (before the chunk's gzip) |
|---|---|
| `highlights/cells/<groupId>` | byte 0 picks the smallest of: 0 sorted-delta varints, 1 bitset, 2 runs, 3 a bitmask over the cells of the group's category, 4 the whole category (no per-cell bytes) |
| `user-defined/codes/<fieldId>` | byte 0 picks: 0 raw codes, 1 runs, 2 codes packed 1, 2 or 4 bits each |
| `user-defined/values/<fieldId>` | the field's values, little-endian float32, one per cell (`NaN`: no value) |
| `analysis/windows`, `analysis/docked` | typed JSON: settings and results as JSON whose numeric result columns are stored as raw typed arrays |

A category-relative membership (3 or 4) is decoded against
`adata.obs[fieldKey] == categoryName`; the session records that category's cell
count and a digest of its cells, and an AnnData whose category holds other
cells is refused (*The category this selection was saved relative to no longer
holds the same cells.*) rather than silently selecting different cells. The
codecs are in `cellucid.session_codecs` (`decode_cell_set`,
`decode_user_defined_codes`, `decode_typed_json`, `digest_cell_set`).

### Closed contributor/chunk inventory

Application validates every manifest entry before decoding any payload. The
current inventory is:

| Contributor | Exact chunk identity |
|---|---|
| `core-state` | `core/state` |
| `field-overlays` | `core/field-overlays` |
| `dockable-layout` | `ui/dockable-layout` |
| `figure-export` | `ui/figure-export` |
| `highlights-meta` | `highlights/meta` |
| `highlights-cells` | `highlights/cells/<groupId>` |
| `user-defined-codes` | `user-defined/codes/<fieldId>` |
| `user-defined-values` | `user-defined/values/<fieldId>` (a pushed continuous column's float32 values) |
| `gene-sets` | `gene-sets/notebook` (the gene sets a notebook pushed; only when there are some) |
| `analysis-windows` | `analysis/windows` |
| `analysis-docked` | `analysis/docked` |
| `cinematic-camera` | `cinematic/camera` |

An unknown chunk, unknown contributor, or mismatched contributor/chunk pair
raises before payload decoding and before target mutation. AnnData application
materializes only highlight and user-defined field chunks; the other current
chunks remain in the session bundle, while their manifest entries are retained
in `adata.uns` when `store_uns=True`.

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
- you filtered/reordered cells before capturing the session,
- you regenerated the AnnData with different preprocessing,
- you are applying a session to a subset.

Any mismatch raises before the target is mutated.

### Cell ordering assumptions

Highlights and codes are stored as **cell positions** in the dataset (or, for a
whole category, as that category). If your AnnData row order differs, labels
will be assigned to the wrong cells.

A session saved on a **subsample** records which cells were open (the
manifest's `cells`: the specification, the kept rows and their digest). Its
highlights and fields are in the positions of those cells; applying it checks
the digest and that `adata` has the dataset's full cell count, and places
every column at the original rows. A cell that was not open is `<NA>` in a
highlight column and missing in a field's.

### Conflicting column names

Applying a bundle can create columns like `cellucid_highlight__<group_id>`.
If a target column already exists, application raises `ValueError`.

### Backed AnnData

`inplace=False` rejects a backed target before `AnnData.copy()` is called. If
you need an independent result, choose and create it explicitly:

```python
in_memory = backed_adata.to_memory()
result = bundle.apply_to_anndata(
    in_memory,
    expected_dataset_id="my-study-v1",
    inplace=False,
)
```

`inplace=True` is accepted without materializing `X`. It changes the live
backed object's in-memory `obs`/`uns` state and does not write those changes
back to the source H5AD file.

---

## Troubleshooting (symptom → diagnosis → fix)

### Symptom: “Dataset fingerprint mismatch”
Likely causes:
- You’re applying the bundle to a different AnnData (or a differently ordered one).

Fix:
- Apply to the original AnnData used to create the session.
- Pass the exact dataset ID used to create the viewer.

### Symptom: “This file is not a Cellucid session file of format 2 …”
Fix:
- Confirm the file is actually a `.cellucid-session` file saved by this release and not (for example) a `.json` file or a renamed file.

### Symptom: “This session uses the obs column … which this dataset does not have”
The session's `requires` names fields the AnnData lacks: every missing obs
column and gene is listed. Apply the session to the AnnData it was saved on (or
with the same data mapping and gene names).

### Symptom: “The category this selection was saved relative to no longer holds the same cells”
A whole-category highlight is stored as that category. The AnnData's
`obs[fieldKey]` assigns that category to different cells than when the session
was saved; apply the session to the AnnData it was saved on.

### Symptom: “Session chunk uncompressedBytes is smaller than the gzip payload”
What it means:
- A chunk inflates to more bytes than its manifest's `uncompressedBytes`
  declares; decompression stops there (zip-bomb protection).

Fix:
- Treat the file as potentially unsafe or corrupted, and save the session
  again from the viewer.

---

### Symptom: “No highlight columns were added”
Likely causes:
- The session contains no highlight data (no highlights were created).
- Dataset identity or dimensions do not match.
- You disabled highlight application (`add_highlights=False`).

How to confirm:
- Inspect `bundle.list_chunk_ids()` and look for `highlights/meta`.

Fix:
- Create at least one highlight group in the UI and re-export the session.
- Apply only to the exact matching dataset.

---

### Symptom: `KeyError: 'highlights/meta'` (or another chunk id)
What it means:
- That chunk does not exist in the bundle (session didn’t include that feature/state).

Fix:
- Treat missing chunks as “feature not used in that session”.
- Guard your code: check `chunk_id in bundle.list_chunk_ids()` before decoding.

---

### Symptom: “Invalid session manifest (chunk '…' codec)”
What it means:
- A chunk names a codec other than `none` or `gzip`.

Fix:
- Treat the bundle as outside the current contract.
- If you control the writer, create a new bundle using the current `none` or
  `gzip` codec contract.

---

### Symptom: “Unknown current session chunk …” or “requires contributor …”
What it means:
- The manifest contains a chunk outside the closed inventory above, or its
  contributor does not own that exact chunk identity.

Fix:
- Do not edit or partially apply the bundle.
- Capture a new bundle with the matching current web and Python package.

---

### Symptom: “Column already exists: …”
Likely causes:
- You applied the same session twice, or you already have similarly named columns.

Fix:
- Choose non-conflicting `highlights_prefix` and `user_defined_prefix` values,
  or explicitly remove the existing column before applying.

---

## See also

- {doc}`jupyter` for capturing a bundle in-notebook (`viewer.request_session_bundle`)
- {doc}`round_trip` for pulling everything (highlights, labels, scores, gene sets) from a viewer or the web app's exports in one call
- {doc}`../../e_jupyter_hooks/index` for deeper hook/event documentation
