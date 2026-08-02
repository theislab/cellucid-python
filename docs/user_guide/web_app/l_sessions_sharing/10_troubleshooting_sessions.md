# Troubleshooting sessions

**Audience:** everyone reopening or sharing a session

**Time:** 10–45 minutes
**What you’ll learn:**
- how to separate picker, identity, format, cancellation, and resource issues;
- what terminal success looks like; and
- the safest fix for each failure.

## Fast diagnosis

| Symptom | Meaning | First action |
|---|---|---|
| File picker never opens | browser policy, embedded-context restriction, or handler error | use the ordinary browser tab, check Console, and retry from a direct click |
| No dataset is loaded | sessions contain state, not the dataset | load the intended dataset first |
| Dataset identity error | one of source type, dataset id, cell count, or variable count differs | load the exact published generation through the same route |
| Manifest/chunk/profile error | file is not a complete current user bundle | obtain a fresh unchanged **Save State** file |
| Gzip, byte-length, truncation, or trailing-data error | corrupt or dishonest container | transfer/download the original again; do not edit it |
| Progress disappears after Cancel or a newer load | intentional abort | no action; the older restore neither succeeded nor failed |
| **Session fully restored**, but view was already similar | valid complete restore with a subtle visual delta | compare high-signal state below |
| Load is slow | large validated code, highlight, or cache chunks | leave progress open or Cancel; close memory-heavy tabs if needed |

## Know the success boundary

A progress card or an eagerly recognizable camera is not success. Cellucid
reports **Session fully restored** only after:

1. complete framing and manifest validation;
2. exact dataset fingerprint validation;
3. every eager and lazy chunk decode and contributor restore;
4. transaction preparation and commit; and
5. the final UI refresh.

The Session panel then reports **Session loaded successfully**. If any step
fails, earlier mutations roll back and no partial subset is accepted.

## Save State

### The button produces no file

1. Wait until the dataset and controls are ready.
2. Click **Save State** again from the visible Session panel.
3. Check the panel status and browser Console for `Failed to save state`.
4. Check whether your browser or organization blocked downloads from the site
   or embedded Jupyter frame.
5. Open the same application URL in a normal tab and retry if an embedding
   policy blocks downloads.

Chrome, Firefox, Safari/WebKit, and Chromium-based browsers all use the same
current bundle writer. A browser-specific failure is reportable; changing the
file extension or using a different serializer is not a fix.

### Saving takes a long time or the file is large

The likely owners are millions of highlight indices, many categorical code
columns, or analysis-cache artifacts. Save meaningful milestones:

- delete obsolete highlight groups;
- remove unneeded saved views;
- keep only analysis state collaborators need; and
- avoid repeatedly saving near-identical checkpoints.

Do not remove chunks from a saved file. The current inventory is closed and a
manually edited bundle is rejected.

## Load State picker

### Picker does not appear

The picker must be opened by a direct user gesture. Close any covering modal,
click **Load State** once, and inspect Console if nothing happens. Managed
browsers or embedded frames can apply file-input policies; use a standalone tab
with the same application build when permitted.

### File is hidden or cannot be chosen

The file must end exactly in `.cellucid-session`. Unzip an archive first and
remove accidental suffixes such as `.txt`. Renaming an arbitrary file to the
extension does not convert it.

### You canceled the picker

Closing the operating-system picker without choosing a file is a normal no-op.
It creates no restore and should show no red error.

## Complete identity rejection

The fingerprint has exactly:

- `sourceType`;
- `datasetId`;
- `cellCount`;
- `varCount`; and
- `cellOrder`, itself exactly a `dimension` (1, 2, or 3) and a 16-character
  hexadecimal `digest` of the coordinates the session was saved against.

All five must match. A GitHub load and a local folder load are different even
when their files look identical. Fix the route or publication identity, then
retry the unchanged bundle. Cellucid does not restore floating layout alone or
skip filters/highlights.

Read the message before changing anything — each cause has its own:

- **"saved on a different dataset"**, naming both cell and gene counts: one of
  the four scalars differs. Open the intended dataset.
- **"saved while the *N*D view was shown"**: only the cell-order dimension
  differs. Nothing is wrong with the data or the file. Switch back to the
  dimension the message names and load again.
- **"same name and the same number of cells and genes, but its cells are stored
  in a different order"**: the dataset was republished re-sorted under an
  unchanged identity. Load the version the session was saved on, or re-create
  the selections on this one.
- **"saved before Cellucid started recording which cells a selection
  contains"**: the file predates the cell-order record, so its selections can
  never be confirmed. Re-create them and save a new session file. This one is
  not repairable by editing the file.

```{figure} ../../../_static/screenshots/sessions_sharing/refusal-01-different-dataset.png
:alt: Two notifications carrying the same message, that the session was saved on a different dataset than the one open now, naming 120 cells and 6 genes when it was saved against 3,696 cells and 3,753 genes now, and advising to open the dataset the session was saved on.
:width: 760px

The first of the four messages, as it appears. It names **both** sizes — what
the file was saved against and what is open now — so you can tell at a glance
which of the two is the one you meant to open. The message appears twice because
the loading indicator carries it as well as the final result.
```

If the sender changed fields, categories, or gene order while reusing the same
lightweight identity and sizes, publish a new dataset id and create a new
session — those are not covered by the digest. Changed cell order or changed
coordinates are now caught for you. See
{doc}`07_versioning_compatibility_and_dataset_identity`.

## Current-format rejection

The reader accepts one current user-authored profile. Common rejection causes
include:

- extra or missing manifest roots;
- missing, duplicate, aliased, unknown, or reordered singleton chunks;
- missing categorical-code, highlight-membership, or analysis-artifact chunks;
- noncanonical contributor, priority, kind, codec, label, or
  `datasetDependent` metadata;
- mismatched stored or decoded lengths;
- invalid JSON/binary tables; and
- truncated, concatenated, trailing, or malformed gzip data.

Obtain a fresh bundle from **Save State** in the current app. The reader accepts
only the exact schema produced by that current writer.

### You picked a file that is not a session at all

The commonest version of this is not a format problem: it is choosing the wrong
file. Session files are binary and their names are easy to confuse with the
figures and exports saved beside them.

```{figure} ../../../_static/screenshots/sessions_sharing/refusal-02-not-a-session-file.png
:alt: Two notifications reading "Invalid session file (bad MAGIC header)." — one on a loading entry and one as the final result.
:width: 760px

`Invalid session file (bad MAGIC header).` means the file does not begin with
the session container's signature — so it is not a `.cellucid-session` file,
whatever it is named. Check that you picked the right file before investigating
anything else.
```

### Did you download an official `default.cellucid-session`?

That file is an internal, SHA-256-pinned five-chunk sample starting state. It is
applied automatically only through the catalog capability and is intentionally
not a generic manual session.

Choose the official sample, wait for its starting view, then use **Save State**
to create a complete manual file. See {doc}`04_official_sample_states`.

## Cancel and replacement

Pressing **Cancel** on the progress card or starting another manual/automatic
restore aborts the older operation. Its progress is dismissed, its registered
mutations roll back, and it produces neither false success nor a red product
failure.

If you canceled accidentally, wait for any newer operation to finish and start
the desired load again. If an ordinary uncoded `AbortError` appears as a
failure without a user cancel or replacement, collect the diagnostics below;
that is not treated as intentional cancellation.

## “It restored, but looks different”

First require terminal success. Then compare:

1. app build shown in the footer;
2. dataset route, id, cell count, and variable count;
3. live/snapshot layout and active view;
4. camera, dimension, and navigation mode;
5. active obs/var field and category legend;
6. Active filters, including latent-space outlier filtering;
7. highlight pages, group names, counts, and visibility;
8. Camera Path keyframes, Loop, Autoplay, and actual playback state; and
9. analysis windows and cache-backed results.

Different viewport sizes can clamp floating panels so they remain reachable.
That is expected. Missing fields, code columns, highlight membership, or cache
artifacts are not: the strict operation should reject and roll back rather than
report success.

## Progress is slow

The loader yields during large gzip preflight and between lazy chunks so the
browser can process Cancel and a newer request. Very large valid arrays still
take CPU, network, and memory.

- Keep the progress card open if you want the restore to finish.
- Avoid launching another large analysis simultaneously.
- Close unrelated memory-heavy tabs.
- Press **Cancel** if you do not want to wait; do not force-close the page
  unless the browser itself is unresponsive.

If progress never settles and Cancel also cannot run, capture the exact app
build and browser information below.

## Reporting a reproducible bug

Please include:

1. footer build version;
2. browser name/version, operating system, and device architecture;
3. dataset source type, exact dataset id, cell count, and variable count;
4. whether the operation was Save, file Load, URL Load, or automatic official
   sample state;
5. the complete visible error and relevant Console stack;
6. whether Cancel or a newer load was involved;
7. the session file and immutable prepared generation, if sharing is allowed;
   and
8. the smallest reliable reproduction sequence.

Do not post sensitive sessions publicly. See
{doc}`08_security_privacy_and_trust`.

## Next steps

- {doc}`09_edge_cases` (strict expected behavior)
- {doc}`12_reference` (exact bundle and chunk contract)
