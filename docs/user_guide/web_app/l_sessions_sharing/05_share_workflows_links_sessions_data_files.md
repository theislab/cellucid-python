# Share workflows (links, session files and data files)

**Audience:** everyone (especially collaborators working across wet lab + computational)  
**Time:** 20–40 minutes  
**What you'll learn:**
- The three shareable artifacts: links, session files, and dataset files
- Which artifact is appropriate for which collaboration goal
- How to package "the reproducible pair": dataset file + `.cellucid-session`
- How to preserve the exact dataset identity a session requires

**Prerequisites:**
- You can save a session file (see {doc}`03_save_restore_ux`)

---

## The three things you can share (and what they mean)

### 1) Share a link (to the dataset)

What you send:
- a URL (or a GitHub `owner/repo/path`) that lets the recipient load the dataset
  ({doc}`13_share_links_badges_and_embeds`).

What the recipient gets:
- the dataset, but **not** your work on it unless you also send a session file.

Best for:
- "Here is the dataset, go explore."
- demos/tutorials where the dataset is public and hosted.

Not enough for:
- "Please see exactly what I see."

### 2) Share a session file (`.cellucid-session`)

What you send:
- a single file downloaded with **Share…** under **State › Session**.

What the recipient gets:
- what you chose in **Share…**: by default the trackers with their main
  branch and the branch shown, every step of them (including steps you undid),
  and the step you were on, with the views, filters, highlights and analyses
  of each; or the current state only; or the steps only, to replay on their
  own dataset ({doc}`../m_state_history/06_sharing_and_replay`);
- **but only on the matching dataset**. A session saved on a sample, a link
  or a server opens it by itself; for a dataset from your computer they open
  the file first.

Best for:
- "Look at this exact selection, camera and filters."
- review packets ("here are the highlight groups I propose, and how I made
  them").

### 3) Share the dataset file

What you send:
- the dataset file or folder itself, in any format Cellucid reads (an `.h5ad`,
  a `.zarr.zip` or `.zarr` folder, an `.h5mu`, a Seurat or SCE `.rds`, a
  `.loom`, a 10x `.h5`, …), through a shared drive or zipped
  ({doc}`../b_data_loading/07_file_format_expectations`).

What the recipient gets:
- the data itself, offline-capable.

Best for:
- private datasets,
- offline work,
- long-term reproducibility (you can archive the exact file).

---

## The "portable pair" rule (most important)

If your goal is reproducibility and you are sharing with another human, the safest rule is:

> Share **(dataset file or link) + (session file)** together.

Because:
- the dataset file is the data,
- the session file is the history of what you did with it.

If you share only one of them, the recipient will usually be missing something.

---

## Before you send: choose what the file holds

A session file holds history, not only the view on screen: the steps you undid
on the branches you share, their settings and the time each was taken.
**Share…** decides what goes in, and **What's inside** lists it before you save:

- tick only the trackers and branches you mean to share, or share the
  **Current state** only, without how you got there;
- notes and bookmarks are in it by default; conversations with agents and
  Recently deleted only when you tick them; **Times** can be dates only or
  none, and **Dataset** its name only;
- name the tracker and branches so the recipient knows what each line of work
  is ({doc}`../m_state_history/02_branches_trackers_and_comparing`); and
- check the labels you typed (highlight groups, pages, renamed categories)
  against {doc}`08_security_privacy_and_trust`.

---

## Workflow A: send a colleague "exactly what I'm seeing" (recommended)

This is the default 1:1 collaboration workflow.

### Sender steps

1) Confirm the dataset is loaded and looks correct.
2) Save a session: **State → Session → Share… → Save file**.
3) Package the dataset:
   - If it's already hosted somewhere your colleague can access (a sample, a
     link, a server), the session opens it by itself.
   - Otherwise, send the file (an `.h5ad`, `.zarr.zip`, `.rds`, …), or place
     a folder (a `.zarr` store, a 10x `outs/` folder) in a shared drive.
4) Send:
   - the `.cellucid-session` file, and
   - the dataset file, when it was opened from your computer.

### Recipient steps

1) For a dataset from a file, open that file first. A session saved on a
   sample, a link or a server opens its dataset by itself when it is loaded.
2) Load the session file: **State → Session → Load… → Open a session file…** → pick the
   `.cellucid-session` file.
3) Wait for **Session fully restored.**

### What success looks like

- The load ends with **Session fully restored.** and no refusal.
- **State › History** shows the sender's steps, with the step they were on
  marked.
- Camera, active field, filters, and kept views match what the sender
  described.

If the dataset identity differs, Cellucid refuses the whole load and nothing
changes. Go to {doc}`07_versioning_compatibility_and_dataset_identity`.

---

## Workflow B: offline sharing (zip everything)

This is common for private datasets.

### Sender steps

1) Put the dataset file or folder and the `.cellucid-session` file into one
   folder, and zip that folder.
2) Optional: include a short `README.txt` explaining:
   - which dataset this is,
   - which session file to load,
   - what you want the recipient to look at.

### Recipient steps

1) Unzip to a local folder.
2) Open the dataset in Cellucid with **Open…** (**File…** for a file such as an
   `.h5ad`, `.zarr.zip` or `.rds`, **Folder…** for a folder), or drop it on the
   page.
3) Load the session file with **Load… › Open a session file…**.

---

## Workflow C: "publication package" (reproducibility + figures)

For papers, your goal is usually:
- reproducible provenance, and
- stable artifacts.

Recommended package:
- the dataset file (archived, versioned)
- the `.cellucid-session` file (one file holds every step; name its trackers
  and branches by what they are)
- the methods paragraph or the printable report of the steps
  ({doc}`../m_state_history/04_tour_methods_and_report`)
- exported figures (SVG/PNG) from the figure export module
- a short `README.md` describing what each file corresponds to

Why:
- figures are the durable output,
- the session lets collaborators and reviewers reproduce and audit how you got
  there, step by step.

Related: {doc}`../k_figure_export/index`

---

## Common sharing pitfall: dataset identity refusal

This happens when the recipient loads the dataset via a different identity than the session expects.

Common causes:
- the dataset has no stable id of its own (no `uns['cellucid']['id']`), so
  its id comes from the link it was opened from, or, for a file, from its
  format, its file name, its cell and gene counts and its gene names: a
  renamed copy has another id
  ({doc}`../b_data_loading/06_dataset_identity_why_it_matters`),
- the dataset was written again and its dataset id changed,
- the dataset contents changed but kept the same id. A re-ordering or a
  re-computed embedding *is* caught (the fingerprint digests the stored
  coordinates of one embedding), but changed field values, categories or gene
  identities under unchanged counts are not. See
  {doc}`07_versioning_compatibility_and_dataset_identity`.

How to resolve:
- give the dataset a stable id (`prepare(..., dataset_id=...)` writes it into
  `uns['cellucid']`), or open the same file on both sides, or
- create a new session file from the recipient's dataset.

Cellucid never loads part of a session: the dataset fingerprint and the whole
file are checked together, and the load succeeds or nothing changes.

Deep dive: {doc}`07_versioning_compatibility_and_dataset_identity`

---

## Next steps

- {doc}`06_collaboration_best_practices` (how teams should organize sessions)
- {doc}`08_security_privacy_and_trust` (what sessions can leak; safe sharing)
- {doc}`10_troubleshooting_sessions` (share/restore failures)
