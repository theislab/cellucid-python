# Session Saving, Restoring, and Sharing

Cellucid sessions are how you **persist and share “what the app looks like right now”**:
camera, active fields, filters, highlights, views/snapshots, and (some) analysis state.

This section is intentionally written for **mixed audiences**:

- **Wet lab / non-technical collaborators:** click-by-click “save my work and send it to someone”.
- **Computational users:** exact semantics (what is saved, what isn’t, and why things sometimes restore differently).
- **Power users / demo builders:** progressive restore, large-session performance, dataset identity pitfalls.

:::{important}
A session bundle (`.cellucid-session`) **does not contain your dataset**.

To restore a session you must:
1) load the same dataset (or a dataset Cellucid considers “the same”), then
2) load the `.cellucid-session` file.
:::

:::{tip}
If you’re here because “I loaded a session but it looks different”, start with:
- {doc}`07_versioning_compatibility_and_dataset_identity` (dataset identity + mismatch rules)
- {doc}`10_troubleshooting_sessions` (symptom → diagnosis → fix)
:::

---

## Fast path (choose your goal)

| You want to… | Do this | Why it’s the right choice | Start here |
|---|---|---|---|
| Reopen your work later | **Save State** → keep the `.cellucid-session` file | Captures UI state; reproducible “what I was looking at” | {doc}`03_save_restore_ux` |
| Send someone an exact view | Send **dataset export folder** + `.cellucid-session` | Sessions don’t include data; folder + session is the portable pair | {doc}`05_share_workflows_links_bundles_exports` |
| Collaborate on labels with many people | Use **Community Annotation** (GitHub-backed) | Sessions are single-user artifacts; annotation is multi-user | {doc}`../j_community_annotation/index` |

---

## What to read (recommended order)

1) {doc}`01_session_mental_model` (what a session is and how to reason about it)
2) {doc}`02_what_gets_saved_and_restored` (explicit inclusion/exclusion list)
3) {doc}`03_save_restore_ux` (manual Save State / Load State)
4) {doc}`05_share_workflows_links_bundles_exports` (how to share with humans)
5) {doc}`10_troubleshooting_sessions` (the big “why did this fail?” map)

---

## Interface reference

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: Cellucid Session panel showing sample, local-file, remote-server, GitHub, and session-state controls.
:width: 246px

The Session panel presents each loading path separately and keeps Save State and Load State beside the dataset controls.
```

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
