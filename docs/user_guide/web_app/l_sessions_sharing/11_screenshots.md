# Verified session-state capture

Every figure below was captured from the running application. Use this page as a
picture index: if you are not sure whether what you are looking at is the
expected screen, find it here first.

The captures use a small pinned fixture dataset (120 cells, 6 genes) for the
save/restore round trip, and the published **Pancreas** built-in sample for the
official starting state. Your dataset will differ; the controls and messages
will not.

---

## Where the controls are

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-01-session-controls.png
:alt: The Session accordion in full, with a dataset summary table at the top, the sample picker, local, remote and GitHub loading controls, and the Session state block with Save State and Load State at the bottom.
:width: 524px

The **Session** accordion. Save State and Load State sit at the bottom, below
every control that brings a dataset in — the order is the workflow, because a
session is restored onto a dataset that is already open.
```

---

## Saving

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-02-save-state-action.png
:alt: The Session state block with a mouse pointer over the Save State button.
:width: 488px

**Save State** downloads a file. There is no upload and no server round trip.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-03-save-confirmation.png
:alt: A notification reading "Session saved successfully".
:width: 610px

The one confirmation a save produces.
```

---

## Restoring

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-04-restore-confirmation.png
:alt: A stack of notifications ending with "Session fully restored." and "Session loaded successfully".
:width: 592px

A completed restore reports **Session fully restored.** and then **Session
loaded successfully**. Both are required; either one alone is not a completed
restore.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-07-restored-view.png
:alt: The full application window immediately after a session restore, with the dataset drawn in the viewer and the sidebar showing restored control values.
:width: 1440px

The application immediately after a restore.
```

---

## A control that proves the restore reached the viewer

The four renderer settings are the easiest end-to-end check, because they are
visible in the panel *and* change what the GPU draws.

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-05-renderer-controls-default.png
:alt: The renderer settings block at its defaults, with Shader quality Full, Level-of-Detail unchecked and Frustum culling unchecked.
:width: 480px

Defaults, before anything is saved.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-06-renderer-controls-changed.png
:alt: The renderer settings block with Shader quality Light, Level-of-Detail checked, Force LOD level 4 and Frustum culling checked.
:width: 480px

The state that was saved.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-08-renderer-controls-restored.png
:alt: The renderer settings block after restoring, matching the saved state exactly.
:width: 480px

The same four values after **Load State** into a freshly started app. See
{doc}`03_save_restore_ux` for why this particular check is worth doing once.
```

---

## Official sample starting state

```{figure} ../../../_static/screenshots/sessions_sharing/official-01-sample-state-applied.png
:alt: The Pancreas built-in sample shortly after being chosen, showing a 3-D orbit view coloured by cell type on a light grid background.
:width: 1440px

Choosing a built-in sample applies its reviewed starting state automatically,
after the scientific dataset has been published and the state bytes verified.
Nothing was loaded by hand. See {doc}`04_official_sample_states`.
```

---

## Refusals

A refused session leaves the application exactly as it was. These are the two
you are most likely to meet.

```{figure} ../../../_static/screenshots/sessions_sharing/refusal-01-different-dataset.png
:alt: Notifications reporting that the session was saved on a different dataset than the one open now, naming both the saved and the current cell and gene counts.
:width: 760px

Wrong dataset. The message names both sizes so you can tell which one you meant.
```

```{figure} ../../../_static/screenshots/sessions_sharing/refusal-02-not-a-session-file.png
:alt: Notifications reading "Invalid session file (bad MAGIC header)."
:width: 760px

Wrong file. `Invalid session file (bad MAGIC header).` means the file does not
begin with the session container's signature, so it is not a session bundle at
all — regardless of what it is named.
```

The other identity refusals — a different embedding dimension, a re-ordered
dataset, and a file saved before cell order was recorded — each carry their own
distinct message. They are listed with their causes in
{doc}`10_troubleshooting_sessions` and specified in
{doc}`07_versioning_compatibility_and_dataset_identity`.
