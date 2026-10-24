# Default sessions named by a catalog

Choosing an official sample should open a useful biological view, not make you
recreate the same setup before you can explore. A dataset catalog
(`datasets.json`) can therefore name one small, reviewed **default session**
for each dataset it lists. Cellucid applies it automatically when the dataset
is opened from that catalog; it is not a second dataset and it is not a manual
**Load** step. The built-in samples use this, and so can any catalog you
publish ({doc}`../b_data_loading/02_local_demo_tutorial`).

## What happens when you choose a dataset from a catalog

1. Cellucid opens the dataset's AnnData file and loads it.
2. The positions, fields, metadata, and other scientific data become the live
   dataset.
3. Only after the dataset has been published, Cellucid reads the default
   session the catalog names for it — if it names one.
4. It downloads that one file, stops reading past the byte length the catalog
   declares, refuses a file of any other length, verifies its SHA-256 against
   the catalog, and confirms that the selected dataset has not changed.
5. Cellucid applies the verified static starting view.

This order matters: an unverified view file is never treated as scientific
data. If the named session is missing, the wrong length, or fails its digest,
Cellucid reports that failure after the dataset is available; it does not
guess another filename or substitute a different view.

A published default session is a seven-chunk starting state, not a manual user
bundle. Downloading it and choosing **Load** is unsupported and the strict
generic reader rejects it. To create a portable manual state, choose the
dataset, wait for its default view, and then use **Save**.

That is also why a file from **Save** cannot simply be renamed and published. A
save differs from a published default in two ways, and the reader refuses both
while the file itself stays perfectly well-formed:

- **Chunks.** A save writes the whole static chunk registry; a published
  default carries seven of them (everything but the Camera Path). The reader refuses a bundle carrying the extra
  chunks, and one of them, a recorded Camera Path, is refused by name.
- **Source.** A save records where its data came from (a local file, a server,
  a link). A published default is applied to the dataset it is listed beside,
  so it records no source, and the dataset it was captured against must be the
  one the catalog lists. Every session bundle records that dataset, and the
  reader compares it before applying anything.

The publishing step of the sample repository
([`cellucid-datasets`](https://github.com/theislab/cellucid-datasets),
`scripts/publish_session_preset.mjs`) performs both, importing the rules from
the viewer rather than restating them. It refuses to drop a chunk that still
holds content instead of discarding it silently, and it checks the capture's
dataset against the dataset it publishes for, so a capture taken against
different data is refused rather than relabelled. If you maintain your own
catalog, use that script rather than editing a save by hand, then write the
file's new `bytes` and `sha256` into `datasets.json`.

## The catalog entry

A dataset's default session is its `session` entry in `datasets.json`:

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

- `path` is the session file's path inside the catalog folder; it ends in
  `.cellucid-session` and stays below the folder.
- `bytes` is that file's exact length: the viewer stops reading past it and
  refuses a file of any other length, so a published session has no size
  limit of its own.
- `sha256` pins the exact bytes of the file, as 64 lowercase hexadecimal
  characters.

All three are required. A dataset without `session` opens without a default
view, and no request is made for one.

A published session is regenerated whenever the reviewed view or the session
format changes, and its `bytes` and `sha256` change with it. Read the current
values out of the catalog rather than copying them from documentation — that
is why the example above shows the field's shape instead of a digest.

The built-in samples' catalog is
[`theislab.github.io/cellucid-datasets/anndata/`](https://theislab.github.io/cellucid-datasets/anndata/datasets.json).
For what each sample is, what to cite, and how far its build is pinned, see
{doc}`../b_data_loading/12_sample_dataset_provenance`.

## What the reviewed view changes

Each official default session:

- activates the `cell_type` observation field;
- restores a modest, static close-up in the 3-D **Orbit** camera;
- uses the light theme and the **Studio (light)** background; and
- starts without active filters or kept snapshots.

The Pancreas view also enables the first two components of its PCA as a 2D
view: that PCA is the dataset's latent space, which Cellucid does not list as a
view on its own.

```{figure} ../../../_static/screenshots/sessions_sharing/official-01-sample-state-applied.png
:alt: The full application window shortly after choosing the Pancreas built-in sample, showing the dataset coloured by cell type in a 3-D orbit view on a light grid background, with the Session accordion open and its Dataset box filled in.
:width: 1440px

The Pancreas sample moments after it was chosen. Nothing was clicked beyond
picking the sample: the dataset published first, then its verified default
session put the view into 3-D **Orbit** on the light studio background with
`cell_type` active. **Save** and **Load** were not involved.
```

You can confirm each of those four from the sidebar in the figure: the
**Visualization** accordion shows `Light` theme, `Studio (light)` background and
`Points` render mode, and the dataset opens on its 3-D embedding with the
navigation mode set to orbit.

These are deliberately **static-only** presets. They contain no
`cinematic/camera` chunk, so they never restore a recorded Camera Path,
**Loop playback**, or **Autoplay**. The static camera position can change when
the state is applied, but moving-camera playback remains off.

## Which datasets get a default session

Only a dataset opened from a catalog — the sample catalog the dataset picker
reads, or a catalog opened from a link in the **Link or server** bar — whose catalog entry names
a `session` gets one. Nothing else is probed:

- a file or folder chosen with **Open…** or dropped from your computer
  (`.h5ad`, `.zarr.zip`, Zarr folder);
- a single file or Zarr store opened from a link;
- a dataset served by `cellucid serve`, `serve_anndata()` or Jupyter; and
- any catalog entry without `session`.

Cellucid does not infer session paths from filenames or look for an unnamed
session. This keeps custom hosting predictable and avoids unrequested network
traffic. A file's own `uns['cellucid']` cannot name a session either: a
`session` key there is ignored
({doc}`../b_data_loading/07_file_format_expectations`).

For those sources, configure the view you want and use the ordinary
{doc}`03_save_restore_ux` workflow. Your own `.cellucid-session` remains a
manual, identity-matched working-state file. See
{doc}`07_versioning_compatibility_and_dataset_identity` for the matching rules.
