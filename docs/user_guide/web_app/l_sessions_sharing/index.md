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
1) load the exact matching source route, dataset id, cell count, and variable count, then
2) choose **Load State** and open the `.cellucid-session` file.

Official built-in starting states use a separate, catalog-advertised path:
Cellucid applies the matching static state automatically after it publishes
and verifies the official dataset generation.
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
| Reopen your work later | **Save State** → keep the `.cellucid-session` file | Captures UI state; reproducible “what I was looking at” | {doc}`03_save_restore_ux` |
| Start an official sample in its reviewed view | Choose the built-in sample; its advertised static state applies automatically | Scientific data publishes first, then the state is integrity-verified | {doc}`04_official_sample_states` |
| Send someone an exact view | Send **dataset export folder** + `.cellucid-session` | Sessions don’t include data; folder + session is the portable pair | {doc}`05_share_workflows_links_bundles_exports` |
| Collaborate on labels with many people | Use **Community Annotation** (GitHub-backed) | Sessions are single-user artifacts; annotation is multi-user | {doc}`../j_community_annotation/index` |

---

## What to read (recommended order)

1) {doc}`01_session_mental_model` (what a session is and how to reason about it)
2) {doc}`02_what_gets_saved_and_restored` (explicit inclusion/exclusion list)
3) {doc}`03_save_restore_ux` (manual Save State / Load State)
4) {doc}`04_official_sample_states` (automatic, integrity-verified static states for official samples)
5) {doc}`05_share_workflows_links_bundles_exports` (how to share with humans)
6) {doc}`10_troubleshooting_sessions` (the big “why did this fail?” map)

---

## Interface reference

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: Cellucid Session panel showing sample, local-file, remote-server, GitHub, and session-state controls.
:width: 524px

The Session panel presents each loading path separately and keeps Save State and Load State beside the dataset controls.
```

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
