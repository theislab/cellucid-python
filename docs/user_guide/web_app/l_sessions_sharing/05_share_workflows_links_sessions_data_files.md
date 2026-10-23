# Share workflows (links, session bundles and data files)

**Audience:** everyone (especially collaborators working across wet lab + computational)  
**Time:** 20–40 minutes  
**What you’ll learn:**
- The three shareable artifacts: links, session bundles, and dataset files
- Which artifact is appropriate for which collaboration goal
- How to package “the reproducible pair”: dataset file + `.cellucid-session`
- How to preserve the exact dataset identity a session requires

**Prerequisites:**
- You can save a session bundle (see {doc}`03_save_restore_ux`)

---

## The three things you can share (and what they mean)

### 1) Share a link (to the dataset)

What you send:
- a URL (or a GitHub `owner/repo/path`) that lets the recipient load the dataset.

What the recipient gets:
- the dataset, but **not** your exact view state unless you also provide a session bundle.

Best for:
- “Here is the dataset, go explore.”
- demos/tutorials where the dataset is public and hosted.

Not enough for:
- “Please see exactly what I see.”

### 2) Share a session bundle (`.cellucid-session`)

What you send:
- a single file downloaded via **Save** under **Session file**.

What the recipient gets:
- the exact UI configuration you saved (camera, filters, highlights, views, etc.),
- **but only after they load the matching dataset first**.

Best for:
- “Look at this exact selection/camera/filters.”
- review packets (“here are the highlight groups I propose”).

### 3) Share the dataset file

What you send:
- the AnnData file itself: a `.h5ad`, a `.zarr.zip`, or a `.zarr` folder
  (through a shared drive, or zipped as a `.zarr.zip`).

What the recipient gets:
- the data itself, offline-capable.

Best for:
- private datasets,
- offline work,
- long-term reproducibility (you can archive the exact file).

---

## The “portable pair” rule (most important)

If your goal is reproducibility and you are sharing with another human, the safest rule is:

> Share **(dataset file) + (session bundle)** together.

Because:
- the dataset file is the data,
- the session bundle is the state.

If you share only one of them, the recipient will usually be missing something.

---

## Workflow A: send a colleague “exactly what I’m seeing” (recommended)

This is the default 1:1 collaboration workflow.

### Sender steps

1) Confirm the dataset is loaded and looks correct.
2) Save a session: **Session → Session file → Save**.
3) Package the dataset file:
   - If it’s already hosted somewhere your colleague can access, share the link.
   - Otherwise, send the `.h5ad` or `.zarr.zip` file, or place the `.zarr`
     folder in a shared drive.
4) Send:
   - the `.cellucid-session` file, and
   - the dataset access method (file or link).

### Recipient steps

1) Load the dataset first (same method if possible).
2) Load the session file: **Load** → pick the `.cellucid-session` file.
3) Wait for **Session fully restored**. Eager and lazy are internal scheduling
   classes; neither is a separate user-visible success state.

### What success looks like

- The operation reaches **Session fully restored** with no rejection.
- Camera, active field, filters, and snapshot layout match what the sender described.

If the dataset identity differs, Cellucid rejects and rolls back the complete
operation. Go to {doc}`07_versioning_compatibility_and_dataset_identity`.

---

## Workflow B: offline sharing (zip everything)

This is common for private datasets.

### Sender steps

1) Put the dataset file (`.h5ad`, `.zarr.zip`, or the `.zarr` folder) and the
   `.cellucid-session` file into one folder, and zip that folder.
3) Optional: include a short `README.txt` explaining:
   - which dataset this is,
   - which session file to load,
   - what you want the recipient to look at.

### Recipient steps

1) Unzip to a local folder.
2) Open the dataset in Cellucid with **H5AD**, **Zarr ZIP** or **Zarr folder**
   (or drop it on the viewer).
3) Load the session file with **Load**.

---

## Workflow C: “publication package” (reproducibility + figures)

For papers, your goal is usually:
- reproducible provenance, and
- stable artifacts.

Recommended bundle:
- the dataset file (archived, versioned)
- one or more `.cellucid-session` files (key milestones)
- exported figures (SVG/PNG) from the figure export module
- a short `README.md` describing what each session corresponds to

Why:
- figures are the durable output,
- sessions let collaborators reproduce and audit the interactive context.

Related: {doc}`../k_figure_export/index`

---

## Common sharing pitfall: dataset identity rejection

This happens when the recipient loads the dataset via a different identity than the session expects.

Common causes:
- the dataset has no stable id of its own (no `uns['cellucid']['id']`), so
  the id comes from how it was opened — a file name, a link — and differs
  between sender and recipient
  ({doc}`../b_data_loading/06_dataset_identity_why_it_matters`),
- the dataset was written again and its dataset id changed,
- the dataset contents changed but kept the same id. A re-ordering or a
  re-computed embedding *is* caught — the fingerprint digests the displayed
  coordinates — but changed field values, categories or gene identities under
  unchanged counts are not. See
  {doc}`07_versioning_compatibility_and_dataset_identity`.

How to resolve:
- give the dataset a stable id (`prepare(..., dataset_id=...)` writes it into
  `uns['cellucid']`), or open it the same way on both sides, or
- create a new session bundle from the recipient’s dataset load method.

Cellucid never salvages layout or skips only the scientific chunks. The exact
fingerprint and complete current chunk profile succeed together or the whole
restore fails.

Deep dive: {doc}`07_versioning_compatibility_and_dataset_identity`

---

## Next steps

- {doc}`06_collaboration_best_practices` (how teams should organize sessions)
- {doc}`08_security_privacy_and_trust` (what sessions can leak; safe sharing)
- {doc}`10_troubleshooting_sessions` (share/restore failures)
