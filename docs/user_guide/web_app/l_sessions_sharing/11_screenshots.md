# Verified session-state capture

Every figure below was captured from the running application. Use this page as a
picture index: if you are not sure whether what you are looking at is the
expected screen, find it here first.

The captures use a small pinned fixture dataset (120 cells, 6 genes) for the
save/restore round trip, and the published **Pancreas** built-in sample for its
catalog default session. Your dataset will differ; the controls and messages
will not.

---

## Where the controls are

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-01-session-controls.png
:alt: The Session accordion showing the Dataset box for the loaded sample, the Open data group with the Sample dataset dropdown, the Open… button and the Link or server bar, the Session file group with Save and Load below them, and the collapsed Data mapping group.
:width: 540px

The **Session** accordion. **Save** and **Load** sit under **Session file**, below
every control that brings a dataset in — the order is the workflow, because a
session is restored onto a dataset: one already open, or the one **Load**
opens first when the session records where it is.
```

---

## Saving

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-02-save-state-action.png
:alt: The Session file group with its (i), Save and Load side by side, and the pointer on Save.
:width: 504px

**Save** downloads a file. There is no upload and no server round trip.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-03-save-confirmation.png
:alt: A notification reading "Session saved successfully" after Save produced a session bundle.
:width: 520px

The one confirmation a save produces.
```

---

## Restoring

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-04-restore-confirmation.png
:alt: The one notification a completed restore publishes: the restore's progress entry, now titled "Loaded" and reading "Session fully restored.".
:width: 520px

A completed restore reports **Session fully restored.** once, on its own
progress entry. Nothing else announces it.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-07-restored-view.png
:alt: The whole application window immediately after a session restore, with the dataset drawn in the viewer and the sidebar showing the restored control values.
:width: 1440px

The application immediately after a restore.
```

---

## A control that proves the restore reached the viewer

The four renderer settings are the easiest end-to-end check, because they are
visible in the panel *and* change what the GPU draws.

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-05-renderer-controls-default.png
:alt: The Cells and Performance controls at their defaults: Shader quality Full (lighting + fog), Point size, the Depth cues sliders, and under Performance Level of detail (LOD) and Frustum culling unchecked and Antialiasing checked.
:width: 496px

Defaults, before anything is saved.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-06-renderer-controls-changed.png
:alt: The same controls changed: Shader quality Light, Level of detail (LOD) checked with Force LOD level at 4, and Frustum culling checked.
:width: 496px

The state that was saved.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-08-renderer-controls-restored.png
:alt: The same controls after a saved session was loaded into a fresh page: Shader quality Light, Level of detail (LOD) checked with Force LOD level 4, and Frustum culling checked, as they were saved.
:width: 496px

The same four values after **Load** into a freshly started app. See
{doc}`03_save_restore_ux` for why this particular check is worth doing once.
```

---

## Catalog default session (official sample)

```{figure} ../../../_static/screenshots/sessions_sharing/official-01-sample-state-applied.png
:alt: The whole application window shortly after the Pancreas built-in sample was chosen, showing the embedding coloured by cell type in a 3-D orbit view on the light Studio backdrop, with the Session accordion open beside it.
:width: 1440px

Choosing a built-in sample applies the default session its catalog names,
after the scientific dataset has been published and the session's bytes verified.
Nothing was loaded by hand. See {doc}`04_official_sample_states`.
```

---

## Refusals

A refused session leaves the application exactly as it was. These are the two
you are most likely to meet.

```{figure} ../../../_static/screenshots/sessions_sharing/refusal-01-different-dataset.png
:alt: A notification titled "Loading failed" reading "Session file: This session was saved on a different dataset than the one that is open now (3,696 cells and 3,753 genes when it was saved, 219,731 cells and 5,754 genes now). Open the dataset the session was saved on, then load the session again."
:width: 760px

Wrong dataset. The message names both sizes so you can tell which one you meant.
```

```{figure} ../../../_static/screenshots/sessions_sharing/refusal-02-not-a-session-file.png
:alt: A Cellucid notification reading "This file is not a Cellucid session file of format 2 (it does not start with CELLUCID_SESSION/2).", with a dismiss button at its right, after a file that is not a session bundle was chosen in Load.
:width: 760px

Wrong file. `This file is not a Cellucid session file of format 2 (it does not start with CELLUCID_SESSION/2).` means the file does not
begin with the session container's signature, so it is not a session bundle at
all — regardless of what it is named.
```

The other identity refusals — an embedding or data-mapping choice the dataset
does not have, a re-ordered dataset, and missing fields or genes — each carry
their own distinct message. They are listed with their causes in
{doc}`10_troubleshooting_sessions` and specified in
{doc}`07_versioning_compatibility_and_dataset_identity`.
