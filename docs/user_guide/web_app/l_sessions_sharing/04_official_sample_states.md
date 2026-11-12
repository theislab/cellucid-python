# Starting views named by a catalog

Choosing an official sample should open a useful biological view, not make you
recreate the same setup before you can explore. A dataset catalog
(`datasets.json`) can therefore name one small, reviewed **starting view** for
each dataset it lists, in its `session` entry. Cellucid puts it back by itself
when the dataset is opened from that catalog; it is not a second dataset and
it is not a manual **Load…** step. The built-in samples use this, and so can any
catalog you publish ({doc}`../b_data_loading/02_local_demo_tutorial`).

## What happens when you choose a dataset from a catalog

1. Cellucid opens the dataset's AnnData file and loads it.
2. The positions, fields, metadata, and other scientific data become the live
   dataset.
3. Only after the dataset has been published, Cellucid reads the starting
   view the catalog names for it, if it names one.
4. It downloads that one file, stops reading past the byte length the catalog
   declares, refuses a file of any other length, verifies its SHA-256 against
   the catalog, and confirms that the selected dataset has not changed.
5. Cellucid applies the verified static starting view, then starts the
   dataset's history. The starting view is not a step: it is where the
   history begins, its first entry (**Opened**), so **Undo** never goes back
   past it.

Cellucid says nothing when it applies the view: it is part of opening the
sample, not a session you loaded, so the camera and colors simply change.
Only a view that cannot be applied is reported, once, as **Sample view
unavailable** ("The dataset is open, but its published sample view could not
be applied.").

A view of another session format, as can happen for minutes while a release is
deployed, is not read and is no failure: the dataset opens on its default view
with one short line that its starting view isn't available right now.

This order matters: an unverified view file is never treated as scientific
data. If the named file is missing, the wrong length, or fails its digest,
Cellucid reports that after the dataset is available; it does not guess
another filename or substitute a different view.

## What a starting view is

A starting view is a session file of the current format (format 4) with **one
tracker of one root and no step**: the state the dataset opens in, not a
history. Cellucid checks that on the file's manifest, before it reads any step
or value:

- **one tracker**, whose first entry is the step shown and its only one;
- **no source**: it ships beside its dataset, so it does not say where the
  dataset opens from;
- **no cells**: it is captured on every cell and applies to any subsample of
  its dataset. On a subsample the cell-order digest is not compared (the
  catalog pins the file by SHA-256 instead); everything else in the dataset
  fingerprint is;
- and, checked once its one step is read, its one state holds **neither the
  cells to open nor the layout**: which cells are open stays your choice (or
  the link's), and your panels, sidebar and **theme** stay yours, as does
  the viewer background (held in its settings, not put back).

It is read like any other session file: the same signature, manifest and
fingerprint checks, and the same check that the dataset has every field the
view uses. So a starting view must carry the dataset's cell-order record and
data mapping, and it must have been captured in the dataset's default data
mapping, with at most more embeddings turned on.

The same dataset read again on other cells (**Cells to open**) keeps its
history and puts no starting view back.

There is no need to download a starting view and load it: opening the
dataset from its catalog applies it. Loaded by hand (**Load… › Open a session
file…**), it is read like any other session file, and its one tracker,
**Starting view**, is added beside yours. To share a view you made from a
sample, choose the sample, work from its starting view, and use **Share…**.

## How a starting view is made and published

A file from **Share…** is not a starting view: it holds your history (every
tracker and step), the cells it was open on and your layout. Publishing turns
the step it shows into the one state of a starting view and drops the rest.
The sample repository
([`cellucid-datasets`](https://github.com/theislab/cellucid-datasets)) does
this in two scripts, both using the viewer's own code rather than restating
its rules:

1. `scripts/capture_session_preset.mjs` opens the sample from a local copy of
   the catalog in a headless browser (its current starting view is put back,
   exactly as on the site), optionally turns on more embeddings
   (`--enable <obsm key>:<columns>`) and dispatches actions of the viewer
   (`--action '{"type": …, "params": …}'`, repeatable, for example to color
   by a field or set the camera), and saves the history as a session file. You
   can instead open the sample in the app, frame the view by hand and choose
   **Share…** › **Save file**.
2. `scripts/publish_session_preset.mjs` checks the capture against the store
   and publishes it. It reads the store's cell and gene counts and default
   data mapping and recomputes the cell-order digest from the store's
   coordinates, so a capture taken against other data is refused rather than
   relabeled. The viewer's publish step (`trimToPublishedDefault`) then makes
   the step shown the one state, drops every other step, and names the view
   with the identity the catalog serves the dataset under. It refuses a
   capture whose steps made content (highlights, gene sets, a camera path,
   cells hidden one by one, renamed or deleted fields and categories, fields
   you made), naming it, because dropping the steps would lose that content
   silently. The script writes the file and its catalog entry.

If you maintain your own catalog, use those scripts rather than editing a file
by hand, then check the file's new `bytes` and `sha256` in `datasets.json`.

## The catalog entry

A dataset's starting view is its `session` entry in `datasets.json`:

```json
{
  "schema": 2,
  "default": "pancreas",
  "datasets": [
    {
      "id": "pancreas",
      "path": "pancreas.zarr/",
      "session": {
        "path": "pancreas.cellucid-session",
        "bytes": 123456,
        "sha256": "<64 lowercase hexadecimal characters>"
      }
    }
  ]
}
```

- `path` is the file's path inside the catalog folder; it ends in
  `.cellucid-session` and stays below the folder.
- `bytes` is that file's exact length: the viewer stops reading past it and
  refuses a file of any other length, so a published starting view has no
  size limit of its own.
- `sha256` pins the exact bytes of the file, as 64 lowercase hexadecimal
  characters.

All three are required. A dataset without `session` opens without a starting
view, and no request is made for one.

A starting view is published again whenever the reviewed view or the session
format changes, and its `bytes` and `sha256` change with it. Read the current
values out of the catalog rather than copying them from documentation; that
is why the example above shows the field's shape instead of a digest.

The built-in samples' catalog is
[`theislab.github.io/cellucid-datasets/anndata/`](https://theislab.github.io/cellucid-datasets/anndata/datasets.json).
For what each sample is, what to cite, and how far its build is pinned, see
{doc}`../b_data_loading/12_sample_dataset_provenance`.

## What the reviewed view changes

Each built-in sample's starting view:

- colors the view by the `cell_type` observation field;
- shows the dataset's default 3-D embedding in the 3-D **Orbit** camera, at a
  reviewed static position;
- uses the `Points` render mode and the `Full (lighting + fog)` shader;
- has one view, with no kept views, no filters, no highlights and no
  centroid labels; and
- holds an empty camera path with **Autoplay** off.

It does not set the theme or the viewer background: the theme is part of your
layout, which a starting view leaves alone, and the background is yours as the
theme is. The file holds the background it was captured on (**Studio
(light)**), because a state lists every setting, but it is not put back; a
background you never picked follows the theme.

The Pancreas view also turns on the first two components of its PCA as a 2-D
embedding: that PCA is the dataset's latent space, which Cellucid does not
list as an embedding on its own.

```{figure} ../../../_static/screenshots/sessions_sharing/official-01-sample-state-applied.png
:alt: The whole application window shortly after the Pancreas built-in sample was chosen, showing the embedding colored by cell type in a 3-D orbit view on the Studio (light) background, with the sidebar beside it.
:width: 1440px

The Pancreas sample moments after it was chosen. Nothing was clicked beyond
picking the sample: the dataset published first, then its verified starting
view put the view into 3-D **Orbit** with `cell_type` active (the light
theme's default background, **Studio (light)**, is the page's own), and the history started from it. **Share…** and **Load…**
were not involved.
```

You can confirm these from the sidebar: **Visualization** shows the `Points`
render mode, the dataset opens on its 3-D
embedding with the navigation mode set to orbit, and **State › History**
starts at **Opened**.

These are deliberately **static** views. Their camera path is empty, so they
never start a recorded camera path, **Loop playback**, or **Autoplay**. The
camera position changes when the view is applied, but no moving-camera
playback starts.

## Which datasets get a starting view

Only a dataset opened from a catalog (the sample catalog the dataset picker
reads, or a catalog opened from a link in the **Link or server** bar) whose
catalog entry names a `session` gets one. Nothing else is probed:

- a file or folder chosen with **Open…** or dropped from your computer, in
  any format Cellucid reads;
- a single file, store or folder opened from a link;
- a dataset served by `cellucid serve`, `serve_anndata()` or Jupyter; and
- any catalog entry without `session`.

Cellucid does not infer session paths from filenames or look for an unnamed
session. This keeps custom hosting predictable and avoids unrequested network
traffic. A file's own `uns['cellucid']` cannot name a session either: a
`session` key there is ignored
({doc}`../b_data_loading/07_file_format_expectations`).

For those sources, set up the view you want and use the ordinary
{doc}`03_save_restore_ux` workflow. Your own `.cellucid-session` remains a
manual, identity-matched working-state file. See
{doc}`07_versioning_compatibility_and_dataset_identity` for the matching rules.
