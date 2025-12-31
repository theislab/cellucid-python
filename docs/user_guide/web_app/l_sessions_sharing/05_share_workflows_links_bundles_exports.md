# Share workflows (links vs bundles vs exports)

**Audience:** everyone (especially collaborators working across wet lab + computational)  
**Time:** 20–40 minutes  
**What you’ll learn:**
- The three shareable artifacts: links, session bundles, and export folders
- Which artifact is appropriate for which collaboration goal
- How to package “the reproducible pair”: dataset export + `.cellucid-session`
- How to avoid the common “dataset mismatch” sharing pitfall

**Prerequisites:**
- You can save a session bundle (see {doc}`03_save_restore_ux`)

---

## The three things you can share (and what they mean)

### 1) Share a link (to the dataset)

What you send:
- a URL (or a GitHub `owner/repo/path`) that lets the recipient load the dataset.

What the recipient gets:
- the dataset, but **not** your exact view state unless you also provide a session bundle or you set up auto-restore.

Best for:
- “Here is the dataset, go explore.”
- demos/tutorials where the dataset is public and hosted.

Not enough for:
- “Please see exactly what I see.”

### 2) Share a session bundle (`.cellucid-session`)

What you send:
- a single file downloaded via **Save State**.

What the recipient gets:
- the exact UI configuration you saved (camera, filters, highlights, views, etc.),
- **but only after they load the matching dataset first**.

Best for:
- “Look at this exact selection/camera/filters.”
- review packets (“here are the highlight groups I propose”).

### 3) Share the dataset export folder

What you send:
- the exported dataset directory (or a zipped copy of it).

What the recipient gets:
- the data itself, offline-capable.

Best for:
- private datasets,
- offline work,
- long-term reproducibility (you can archive the exact export).

---

## The “portable pair” rule (most important)

If your goal is reproducibility and you are sharing with another human, the safest rule is:

> Share **(dataset export folder) + (session bundle)** together.

Because:
- the export folder is the data,
- the session bundle is the state.

If you share only one of them, the recipient will usually be missing something.

---

## Workflow A: send a colleague “exactly what I’m seeing” (recommended)

This is the default 1:1 collaboration workflow.

### Sender steps

1) Confirm the dataset is loaded and looks correct.
2) Save a session: **User data → Session state → Save State**.
3) Package the dataset export folder:
   - If it’s already hosted somewhere your colleague can access, share the link.
   - Otherwise, zip the export folder and send it (or place it in a shared drive).
4) Send:
   - the `.cellucid-session` file, and
   - the dataset access method (folder or link).

### Recipient steps

1) Load the dataset first (same method if possible).
2) Load the session file: **Load State** → pick the `.cellucid-session` file.
3) Wait for eager restore; if the session is large, wait for lazy restore too.

### What success looks like

- No “dataset mismatch” warning.
- Camera, active field, filters, and snapshot layout match what the sender described.

If you get a dataset mismatch warning, go to {doc}`07_versioning_compatibility_and_dataset_identity`.

---

## Workflow B: “open this dataset and it’s already set up” (auto-restore)

This is best for demos and tutorials.

You publish:
- export folder + session file(s) + `state-snapshots.json`.

Then any user who loads the dataset will automatically restore the latest session bundle on startup.

Steps:
1) Place one or more `.cellucid-session` files in the export folder.
2) Add/update `state-snapshots.json` to list them (last entry is restored).
3) Host the folder or distribute it.

Full details: {doc}`04_auto_restore_latest_from_dataset_exports`.

---

## Workflow C: offline sharing (zip everything)

This is common for private datasets.

### Sender steps

1) Zip the dataset export folder (`.zip` or `.tar.gz`).
2) Put the `.cellucid-session` file *inside* the folder before zipping (recommended), or zip it alongside.
3) Optional: include a short `README.txt` explaining:
   - which dataset this is,
   - which session file to load,
   - what you want the recipient to look at.

### Recipient steps

1) Unzip to a local folder.
2) Load the export folder in Cellucid (local folder picker).
3) Load the session file (if it didn’t auto-restore).

---

## Workflow D: “publication package” (reproducibility + figures)

For papers, your goal is usually:
- reproducible provenance, and
- stable artifacts.

Recommended bundle:
- dataset export folder (archived, versioned)
- one or more `.cellucid-session` files (key milestones)
- exported figures (SVG/PNG) from the figure export module
- a short `README.md` describing what each session corresponds to

Why:
- figures are the durable output,
- sessions let collaborators reproduce and audit the interactive context.

Related: {doc}`../k_figure_export/index`

---

## Common sharing pitfall: “it loads but warns dataset mismatch”

This happens when the recipient loads the dataset via a different identity than the session expects.

Common causes:
- the sender loaded from GitHub, the recipient loaded from a local folder (different source type),
- the dataset was re-exported and the dataset id changed,
- the dataset contents changed but kept the same id (dangerous; can silently misapply highlights).

How to resolve:
- align on the same dataset access method (same source type + dataset id), or
- create a new session bundle from the recipient’s dataset load method.

Deep dive: {doc}`07_versioning_compatibility_and_dataset_identity`

---

## Next steps

- {doc}`06_collaboration_best_practices` (how teams should organize sessions)
- {doc}`08_security_privacy_and_trust` (what sessions can leak; safe sharing)
- {doc}`10_troubleshooting_sessions` (share/restore failures)
