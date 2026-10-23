# Session Saving, Restoring, and Sharing

Cellucid sessions are how you **persist and share “what the app looks like right now”**:
camera, active fields, filters, highlights, views/snapshots, and (some) analysis state.

This section is intentionally written for **mixed audiences**:

- **Wet lab / non-technical collaborators:** click-by-click “save my work and send it to someone”.
- **Computational users:** exact semantics (what is saved, what is excluded, and why rejection is atomic).
- **Power users / demo builders:** eager/lazy ordering, large-session performance, and strict dataset identity.

:::{important}
A session bundle (`.cellucid-session`) **does not contain your dataset**.

For a session you saved yourself, you must:
1) load the dataset it was saved on (the same dataset id, cell count, and variable count), then
2) choose **Load** under **Session file** and open the `.cellucid-session` file.

A dataset opened from a catalog (the built-in samples included) can come with
a default session the catalog names: Cellucid applies it automatically after
the dataset is published and the session's bytes are verified.
:::

:::{tip}
If you’re here because “I loaded a session but it looks different”, start with:
- {doc}`07_versioning_compatibility_and_dataset_identity` (exact identity and rejection rules)
- {doc}`10_troubleshooting_sessions` (symptom → diagnosis → fix)
:::

---

## Fast path (choose your goal)

| You want to… | Do this | Why it’s the right choice | Start here |
|---|---|---|---|
| Reopen your work later | **Save** → keep the `.cellucid-session` file | Captures UI state; reproducible “what I was looking at” | {doc}`03_save_restore_ux` |
| Start a sample in its reviewed view | Choose the dataset from its catalog; the default session the catalog names applies automatically | Scientific data publishes first, then the session is integrity-verified | {doc}`04_official_sample_states` |
| Send someone an exact view | Send the **dataset file** + `.cellucid-session` | Sessions don’t include data; file + session is the portable pair | {doc}`05_share_workflows_links_sessions_data_files` |
| Collaborate on labels with many people | Use **Community annotation** (GitHub-backed) | Sessions are single-user artifacts; annotation is multi-user | {doc}`../j_community_annotation/index` |

---

## What to read (recommended order)

1) {doc}`01_session_mental_model` (what a session is and how to reason about it)
2) {doc}`02_what_gets_saved_and_restored` (explicit inclusion/exclusion list)
3) {doc}`03_save_restore_ux` (manual **Save** / **Load**)
4) {doc}`04_official_sample_states` (automatic, integrity-verified default sessions named by a catalog)
5) {doc}`05_share_workflows_links_sessions_data_files` (how to share with humans)
6) {doc}`10_troubleshooting_sessions` (the big “why did this fail?” map)

From Python: read a session from a notebook viewer
({doc}`../../python_package/e_jupyter_hooks/10_session_bundles_get_session_bundle`)
and apply it back to an AnnData
({doc}`../../python_package/e_jupyter_hooks/11_session_bundles_apply_session_to_anndata`).

---

## Interface reference

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: The Cellucid Session panel: the Dataset summary, the Open data group with Sample dataset, Files, Server and Open from link, the Session file group with its Save and Load buttons, and the Data mapping group.
:width: 524px

The Session panel gathers every loading path under **Open data** and keeps the **Save** and **Load** buttons under **Session file**, directly below it.
```

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
