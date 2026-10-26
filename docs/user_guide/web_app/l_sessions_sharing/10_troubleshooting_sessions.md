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
| No dataset is loaded | sessions contain state, not the dataset; one saved on a sample, link or server opens its dataset on **Load**, one saved on a local file or in a notebook cannot | open that dataset first, then **Load** |
| Dataset identity error | one of dataset id, cell count, variable count, data mapping, or cell order differs, or the dataset lacks a field or gene the session uses — each with its own message | read the exact message first, then load the exact dataset the session was saved on |
| "not a Cellucid session file of format 2" | the file does not start with the session signature `CELLUCID_SESSION/2` | check that you picked a session file saved by this release |
| Manifest/chunk/profile error | file is not a complete current user bundle | obtain a fresh unchanged file from **Save** |
| Gzip, byte-length, truncation, or trailing-data error | corrupt or dishonest container | transfer/download the original again; do not edit it |
| Progress disappears after Cancel or a newer load | intentional abort | no action; the older restore neither succeeded nor failed |
| **Session fully restored**, but view was already similar | valid complete restore with a subtle visual delta | compare high-signal state below |
| Load is slow | large validated code, highlight, or analysis-result chunks | leave progress open or Cancel; close memory-heavy tabs if needed |

## Know the success boundary

A progress card or an eagerly recognizable camera is not success. Cellucid
reports **Session fully restored** only after:

1. complete framing and manifest validation;
2. exact dataset fingerprint validation;
3. every eager and lazy chunk decode and contributor restore;
4. transaction preparation and commit; and
5. the final UI refresh.

That notification is the restore's one report of success; there is no second
toast or panel status. If any step fails, earlier mutations roll back, no
partial subset is accepted, and the same progress entry ends as **Loading
failed** with the error.

## Save

### The button produces no file

1. Wait until the dataset and controls are ready.
2. Click **Save** again from the visible Session panel.
3. Check for a `Failed to save session` notification and the browser Console for `Failed to save state`.
4. Check whether your browser or organization blocked downloads from the site
   or embedded Jupyter frame.
5. Open the same application URL in a normal tab and retry if an embedding
   policy blocks downloads.

Chrome, Firefox, Safari/WebKit, and Chromium-based browsers all use the same
current bundle writer. A browser-specific failure is reportable; changing the
file extension or using a different serializer is not a fix.

### Saving takes a long time or the file is large

The likely owners are millions of highlight indices, many categorical code
columns, or many analysis results. Save meaningful milestones:

- delete obsolete highlight groups;
- remove unneeded saved views;
- keep only analysis state collaborators need; and
- avoid repeatedly saving near-identical checkpoints.

Do not remove chunks from a saved file. The current inventory is closed and a
manually edited bundle is rejected.

## Load picker

### Picker does not appear

The picker must be opened by a direct user gesture. Close any covering modal,
click **Load** once, and inspect Console if nothing happens. Managed
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

- `datasetId`;
- `name` (shown in messages, never compared);
- `cellCount`;
- `varCount`;
- `cellOrder`, itself exactly the `embedding` the digest was taken over, its
  `dimension` (1, 2, or 3), and a 16-character hexadecimal `digest` of that
  embedding's coordinates; and
- `dataMapping`, the stored embeddings and arrays the dataset was shown with.

All but `name` must match, except that the data mapping may differ when this
dataset can show or use what the session names (its embeddings, graph,
modality, gene names, matrix, …): those are restored. How the dataset was opened does not matter, but a file that states
no id (`uns['cellucid']['id']`) takes one from its file name, so the same data
opened as a renamed file, or from a catalog, can carry a different id. Open the
dataset the message names, then retry the unchanged bundle. Cellucid does not
restore floating layout alone or skip filters/highlights.

Read the message before changing anything — each cause has its own:

- **"saved on a different dataset"**, naming both cell and gene counts: the
  sizes themselves differ. Open the intended dataset.
- **"saved on a different copy of this dataset"**, naming both identities: the
  sizes agree exactly and the publication identity does not. This is one
  dataset reached two ways — most often a state captured against a locally
  served copy of an export the catalog also publishes. The counts are no help
  here, which is why this message names the two identities instead of repeating
  the same two numbers. Open the dataset the state was captured on.
- **"This session shows …, which this dataset does not have"**: the
  session turns on an embedding this dataset lacks (or one it cannot plot,
  with the reason). A session whose embeddings this dataset
  *does* have is not refused: they are restored. Open the dataset version the
  session was saved on.
- **"This session uses … as its neighbour graph, which this dataset does not
  have"** (or modality, gene names, expression matrix, latent space, vector
  fields): the
  session names an array this dataset does not have or cannot use. One it has
  is restored instead of refused.
- **"recorded its cell order on … which this dataset does not have"**: the
  embedding the session was saved in is not part of this dataset.
- **"the same name and the same number of cells and genes, but the cell
  coordinates differ"**: the digest no longer matches. A one-way digest cannot
  say *why*, so the message offers both causes it could be — the cells are
  stored in a different order, or the same cells were exported again from a
  re-computed embedding. Load the version the session was saved on, or re-create
  the selections on this one.
- **"This session uses the observation fields … and the genes …, which the
  dataset that is open now does not have"**: every key the session needs and
  the dataset lacks is named. Open the dataset version the session was saved
  on.
- **"The category this selection was saved relative to no longer holds the
  same cells"**: a highlight made from a category is stored relative to that
  category, and the category's cells changed. Open the dataset version the
  session was saved on.

```{figure} ../../../_static/screenshots/sessions_sharing/refusal-01-different-dataset.png
:alt: A notification titled "Loading failed" reading "Session file: This session was saved on a different dataset than the one that is open now (3,696 cells and 3,753 genes when it was saved, 219,731 cells and 5,754 genes now). Open the dataset the session was saved on, then load the session again."
:width: 760px

The message as it appears: a session saved on the Pancreas sample opened as
`pancreas.h5ad` from this computer, loaded while the Garcia sample is open. (A
session saved on a sample, a link or a server opens its own dataset first
instead, so it is not refused.) It names **both** sizes — what the file was saved
against and what is open now — so you can tell at a glance which of the two is
the one you meant to open.
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

Obtain a fresh bundle from **Save** in the current app. The reader accepts
only the exact schema produced by that current writer.

### You picked a file that is not a session at all

The commonest version of this is not a format problem: it is choosing the wrong
file. Session files are binary and their names are easy to confuse with the
figures and exports saved beside them.

```{figure} ../../../_static/screenshots/sessions_sharing/refusal-02-not-a-session-file.png
:alt: A Cellucid notification reading "This file is not a Cellucid session file of format 2 (it does not start with CELLUCID_SESSION/2).", with a dismiss button at its right, after a file that is not a session bundle was chosen in Load.
:width: 760px

`This file is not a Cellucid session file of format 2 (it does not start with CELLUCID_SESSION/2).` means the file does not begin with
the session container's signature — so it is not a `.cellucid-session` file,
whatever it is named. Check that you picked the right file before investigating
anything else.
```

### Did you download a catalog's default session?

That file is a SHA-256-pinned seven-chunk default session. It is applied
automatically only when its dataset is opened from the catalog that names it,
and is intentionally not a generic manual session.

Choose the dataset from the catalog, wait for its default view, then use **Save**
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
4. every view's embedding, camera, and navigation mode, and whether cameras are linked;
5. active obs/var field and category legend;
6. Filters, including latent-space outlier filtering;
7. highlight pages, group names, counts, and visibility;
8. Camera Path keyframes, Loop, Autoplay, and actual playback state;
9. the analysis panel, analysis windows and their results; and
10. Figure Export settings.

Different viewport sizes can clamp floating panels so they remain reachable.
That is expected. Missing fields, code columns, highlight membership, or
analysis results are not: the strict operation should reject and roll back rather than
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
3. how the dataset was opened, exact dataset id, cell count, and variable count;
4. whether the operation was Save, Load, or an automatic catalog default
   session;
5. the complete visible error and relevant Console stack;
6. whether Cancel or a newer load was involved;
7. the session file and the unchanged dataset file, if sharing is allowed;
   and
8. the smallest reliable reproduction sequence.

Do not post sensitive sessions publicly. See
{doc}`08_security_privacy_and_trust`.

## Next steps

- {doc}`09_edge_cases` (strict expected behavior)
- {doc}`12_reference` (exact bundle and chunk contract)
