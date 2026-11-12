# Verified session-state capture

Every figure below was captured from the running application. Use this page as a
picture index: if you are not sure whether what you are looking at is the
expected screen, find it here first.

The captures use a small pinned fixture dataset (120 cells, 6 genes) for the
save/load round trip, and the published **Pancreas** built-in sample for its
starting view. Your dataset will differ; the controls and messages will not.

---

## Where the controls are

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-01-session-controls.png
:alt: The State section with the Pancreatic endocrinogenesis sample open: Session with Share… and Load…, and the History with the recording switch, the branch shown and the graph of the steps taken so far, the step shown marked.
:width: 540px

The **State** section holds the session's **Share…** and **Load…** above the
history. **Share…** writes the history you choose to one file; **Load…** puts
a file's history back beside yours, at the step it showed.
```

---

## Saving

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-02-save-state-action.png
:alt: The Session group of the State section with its (i), Share… and Load… side by side, and the pointer on Share….
:width: 504px

**Share…** chooses what the file holds and downloads it. There is no upload and no server round trip.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-03-save-confirmation.png
:alt: A notification reading "Saved Pancreatic endocrinogenesis (scVelo).cellucid-session." after Share… › Save file wrote the session file.
:width: 760px

The one confirmation a share produces.
```

---

## Loading

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-04-restore-confirmation.png
:alt: The one notification a completed load leaves, titled "Loaded "Tracker 1 (shared)" beside yours" and reading "Session fully restored.", with an Undo button.
:width: 676px

A completed load ends in one notice: it says which trackers it added beside
yours, reads **Session fully restored.**, and offers **Undo**.
```

```{figure} ../../../_static/screenshots/sessions_sharing/save-restore-07-restored-view.png
:alt: The whole application window immediately after a session file was loaded, with the dataset drawn in the viewer and the State section showing the loaded history, the step shown marked.
:width: 1440px

The application immediately after a load: the history is back, at the step
it showed.
```

---

## A control that proves the load reached the viewer

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

The same four values after the session was loaded into a freshly started app. See
{doc}`03_save_restore_ux` for why this particular check is worth doing once.
```

---

## Catalog starting view (official sample)

```{figure} ../../../_static/screenshots/sessions_sharing/official-01-sample-state-applied.png
:alt: The whole application window shortly after the Pancreas built-in sample was chosen, showing the embedding colored by cell type in a 3-D orbit view on the Studio (light) background, with the sidebar beside it.
:width: 1440px

Choosing a built-in sample puts back the starting view its catalog names,
after the dataset has been published and the view's bytes verified; the
history starts from it. Nothing was loaded by hand. See
{doc}`04_official_sample_states`.
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
:alt: A Cellucid notification reading "This file is not a Cellucid session file of format 4 (it does not start with CELLUCID_SESSION/4).", with a dismiss button at its right, after a file that is not a session file was chosen in Load.
:width: 760px

Wrong file. `This file is not a Cellucid session file of format 4 (it does not
start with CELLUCID_SESSION/4).` means the file does not begin with the session
file's signature, so it is not a session file at all, regardless of what it is
named.
```

The other refusals (an embedding or data-mapping choice the dataset does not
have, a re-ordered dataset, missing fields or genes, a file of another session
format) each carry their own distinct message. They are listed with their
causes in {doc}`10_troubleshooting_sessions` and specified in
{doc}`07_versioning_compatibility_and_dataset_identity`.
