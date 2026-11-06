# Session Saving, Restoring, and Sharing

A Cellucid session file (`.cellucid-session`) holds the **history of the
dataset you are looking at**: every tracker, every step of each, the branches,
and the step that is shown. Loading it puts that history back, at the step it
showed, with every step still there to undo to. **Share…** chooses what it
holds: the history of the trackers and branches you tick, the current state
only, or the steps only (a plan to replay on another dataset). It is how you
keep your work outside this browser and how you send it to someone else.

This section is written for **mixed audiences**:

- **Wet lab / non-technical collaborators:** click-by-click "save my work and send it to someone".
- **Computational users:** exact semantics (what is saved, what is left out, and why a load is all or nothing).
- **Power users / demo builders:** the file format, starting views for a catalog, and strict dataset identity.

:::{important}
A session file **does not contain your dataset**.

**Share…** and **Load…** sit in the sidebar's **State** section, under
**Session**. To reopen a session you saved yourself:
1) choose **Load… › Open a session file…** and pick the `.cellucid-session` file;
2) a session saved on a sample, a data link or a Cellucid server opens its
   dataset first; one saved on a file from your computer needs that file
   (in Chrome and Edge, the handle this browser kept of it opens it again
   with one click; elsewhere you choose the file again).

A dataset opened from a catalog (the built-in samples included) can come with
a **starting view** the catalog names: Cellucid puts it back by itself after
the dataset is published and the session's bytes are verified, and the
history starts from it.
:::

:::{tip}
If you're here because "I loaded a session but it looks different", start with:
- {doc}`07_versioning_compatibility_and_dataset_identity` (exact identity and refusal rules)
- {doc}`10_troubleshooting_sessions` (symptom → diagnosis → fix)

The history itself (steps, undo, branches, trackers, the copy kept in this
browser) has its own chapter: {doc}`../m_state_history/index`.
:::

---

## Fast path (choose your goal)

| You want to… | Do this | Why it's the right choice | Start here |
|---|---|---|---|
| Reopen your work later on this computer | Nothing: every step is kept in this browser, and reopening the dataset offers **Continue where you left off?** | The history is saved as you work | {doc}`../m_state_history/03_saved_states_and_autosave` |
| Keep your work as a file, or move it to another computer | **State › Session › Share…** → **Save file**, keep the `.cellucid-session` file | The file holds the history you chose and opens in any browser | {doc}`03_save_restore_ux` |
| Run the same steps on another dataset | Open the other dataset, then **Load… › Replay a session file here…** | Steps are applied again by name; what does not match is flagged, never dropped | {doc}`../m_state_history/06_sharing_and_replay` |
| Start a sample in its reviewed view | Choose the dataset from its catalog; the starting view the catalog names is put back by itself | The data is published first, then the view is integrity-verified | {doc}`04_official_sample_states` |
| Send someone your work | Send the `.cellucid-session` file, and the dataset file when it was opened from your computer | Sessions don't include data; a sample, link or server opens by itself | {doc}`05_share_workflows_links_sessions_data_files` |
| Collaborate on labels with many people | Use **Community annotation** (GitHub-backed) | Sessions are single-user files; annotation is multi-user | {doc}`../j_community_annotation/index` |

---

## What to read (recommended order)

1) {doc}`01_session_mental_model` (what a session file is and how to reason about it)
2) {doc}`02_what_gets_saved_and_restored` (exact list of what is and is not saved)
3) {doc}`03_save_restore_ux` (**Share…** / **Load…** step by step)
4) {doc}`04_official_sample_states` (starting views a catalog names, put back by themselves)
5) {doc}`05_share_workflows_links_sessions_data_files` (how to share with humans)
6) {doc}`10_troubleshooting_sessions` (the big "why did this fail?" map)

From Python and R: read a session file's history -- its trackers, branches,
steps, notes and turns, and the state at any step -- with
`cellucid.read_session()`
({doc}`../../python_package/g_api_reference_coverage/api/sessions`); apply
the labels, highlight groups, fields and gene sets of one step to an AnnData
with `apply_cellucid_session_to_anndata(adata, session, tracker=..., branch=..., step=...)`,
a plan first, then the change, each apply recorded
({doc}`../../python_package/e_jupyter_hooks/11_session_bundles_apply_session_to_anndata`);
ask a notebook viewer for its session
({doc}`../../python_package/e_jupyter_hooks/10_session_bundles_request_session_bundle`),
or push a file into it with `viewer.load_session(path)`.
In R the same is `cellucid_read_session()`, `cellucid_trackers()` and its
companions, `cellucid_apply_session()` with `cellucid_apply_plan()`, and
`cellucid_load_session()`
({doc}`../../r_package/d_viewing_loading/06_hooks_selection_and_sessions`).
Results come back to an AnnData with `viewer.pull(adata)`
({doc}`../../python_package/g_api_reference_coverage/api/round_trip`).

---

## Interface reference

- **State › Session**: **Share…** (chooses what the file holds, then downloads
  the `.cellucid-session` file; {doc}`../m_state_history/06_sharing_and_replay`)
  and **Load…** (every history this browser keeps, **Open a session file…**
  and **Replay a session file here…**; see
  {doc}`../m_state_history/03_saved_states_and_autosave`).
- **State › History**: the steps a session file holds, as a graph with its
  branches and trackers, shown and changed in place
  ({doc}`../m_state_history/01_the_history`,
  {doc}`../m_state_history/02_branches_trackers_and_comparing`).
- **Session** (the panel above **State**): **Open data**, **Subsampling**
  (**Cells to open**, **Genes read for**) and **Data mapping**. A change of the
  cells to open or of the data mapping is a step of the history, so a session
  file keeps it; opening a dataset is not a step.

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-01-session-controls.png
:alt: The State section with the Pancreatic endocrinogenesis sample open: Session with Share… and Load…, and the History with the recording switch, the branch shown and the graph of the steps taken so far, the step shown marked.
:width: 540px

The **State** section holds the session's **Share…** and **Load…** above the
history they write and read.
```

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: The Session panel with the Pancreatic endocrinogenesis sample loaded: the Dataset box with its name, cells and genes, format and host, description and a collapsed Details; Open data with the Sample dataset dropdown, My datasets naming the open sample, the Open… button and the Link or server bar; the Stored line; and the collapsed Subsampling and Data mapping sections.
:width: 532px

The **Session** panel gathers every way a dataset comes in under **Open data**,
with **Subsampling** and **Data mapping** below it. A session is loaded onto
the dataset opened here.
```

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
