# Official sample starting states

Choosing an official sample should open a useful biological view, not make you
recreate the same setup before you can explore. Each built-in sample therefore
publishes one small, reviewed starting state beside its prepared scientific
files. Cellucid applies that state automatically; it is not a second dataset
and it is not a manual **Load State** step.

## What happens when you choose a sample

1. Cellucid validates and loads the selected prepared dataset.
2. The positions, fields, metadata, and other scientific data become the live
   dataset.
3. Only after the scientific dataset has been published, Cellucid reads the
   state capability advertised for that exact generation.
4. It validates the manifest, downloads the one declared state, verifies the
   state bytes against the catalog's SHA-256, and confirms that the selected
   generation has not changed.
5. Cellucid automatically applies the verified static starting view.

This order matters: an unverified view file is never treated as scientific
data. If an advertised manifest, state, or digest is invalid, Cellucid reports
that failure after the dataset is available; it does not guess another
filename or silently substitute a different preset.

The published `default.cellucid-session` is an internal five-chunk starting
state, not a manual user bundle. Downloading it and choosing **Load State** is
unsupported and the strict generic reader rejects it. To create a portable
manual state, choose the sample, wait for this verified view, and then use
**Save State**.

## The exact publication contract

Every official export generation contains both sidecars:

| Built-in sample | State bundle | Manifest |
| --- | --- | --- |
| Suo | `exports/suo/default.cellucid-session` | `exports/suo/state-snapshots.json` |
| Garcia | `exports/garcia/default.cellucid-session` | `exports/garcia/state-snapshots.json` |
| He | `exports/he/default.cellucid-session` | `exports/he/state-snapshots.json` |
| Kanemaru | `exports/kanemaru/default.cellucid-session` | `exports/kanemaru/state-snapshots.json` |
| Pancreas | `exports/pancreas/default.cellucid-session` | `exports/pancreas/state-snapshots.json` |

The manifest has one exact root, `states`, and one exact entry:

```json
{ "states": ["default.cellucid-session"] }
```

The matching entry in
[`exports/datasets.json`](https://github.com/theislab/cellucid-datasets/blob/main/exports/datasets.json)
advertises `state_manifest` and `state_sha256` as a pair. For example, the Suo
entry declares:

```json
{
  "state_manifest": "state-snapshots.json",
  "state_sha256": "682a065b214eee5e3c7c3bd37dccacf6d39f53de2ac697d0b20d06305f50e578"
}
```

`state_manifest` names the manifest beside that dataset generation.
`state_sha256` pins the exact bytes of `default.cellucid-session`; it is not a
digest of the JSON manifest. Both fields are required for the capability to be
advertised.

A published state is regenerated whenever the reviewed view or the session
format changes, and its `state_sha256` changes with it. Read the current value
out of the catalog rather than copying one from documentation.

You can inspect all five generation directories in the
[`cellucid-datasets` exports](https://github.com/theislab/cellucid-datasets/tree/main/exports).

## What the reviewed view changes

Each official state:

- activates the `cell_type` observation field;
- restores a modest, static close-up in the 3-D **Orbit** camera;
- uses the light theme and grid background; and
- starts without active filters or kept snapshots.

These are deliberately **static-only** presets. They contain no
`cinematic/camera` chunk, so they never restore a recorded Camera Path,
**Loop playback**, or **Autoplay**. The static camera position can change when
the state is applied, but moving-camera playback remains off.

## Which catalogs opt in—and which sources are not probed

Automatic starting states belong to a `local-demo` catalog generation whose
entry explicitly advertises the paired `state_manifest` and `state_sha256`
fields. The official catalog uses this capability. A custom `exportsBaseUrl`
catalog may deliberately opt in by publishing the same paired fields and exact
sidecars.

Custom hosting is therefore not rejected as a class: only explicit catalog
advertisement enables the state requests. Jupyter, GitHub, local-user, and
remote sources are not probed, and neither is any catalog entry without the
complete capability.

- Jupyter (`jupyter`);
- a GitHub repository (`github-repo`);
- a locally selected prepared export (`local-user`);
- a remote prepared-data server (`remote`); and
- any catalog entry without the complete state capability.

Cellucid does not infer sidecar paths from filenames or look for an
unadvertised state. This keeps custom hosting predictable and avoids
unrequested network traffic.

For those zero-probe sources, configure the view you want and use the ordinary
{doc}`03_save_restore_ux` workflow. Your own `.cellucid-session` remains a
manual, identity-matched working-state file. See
{doc}`07_versioning_compatibility_and_dataset_identity` for the matching rules.
