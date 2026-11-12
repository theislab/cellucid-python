# Troubleshooting sessions

**Audience:** everyone reopening or sharing a session

**Time:** 10–45 minutes
**What you'll learn:**
- how to separate picker, identity, format, cancellation, and resource issues;
- what a finished load looks like; and
- the safest fix for each failure.

## Fast diagnosis

| Symptom | Meaning | First action |
|---|---|---|
| **Share…** is grayed out | no dataset is open | open a dataset first |
| File picker never opens | browser policy, embedded-context restriction, or handler error | use the ordinary browser tab, check Console, and retry from a direct click |
| No dataset is loaded | sessions contain the history, not the dataset; one saved on a sample, link or server opens its dataset when loaded, one saved on a local file only where this browser kept the file (Chrome and Edge: **Reopen <name>**), one saved in a notebook cannot | click **Reopen <name>** if offered; otherwise open that dataset first, then load the file again |
| Dataset identity error | one of dataset id, cell count, gene count, data mapping, or cell order differs, or the dataset lacks a field or gene the session uses, each with its own message | read the exact message first, then load the exact dataset the session was saved on |
| "saved by an earlier version of Cellucid (session format 3)" | the session was saved before this release, in a format that held one saved state rather than the history | open it in the version of Cellucid that saved it, or re-create it in this version and save it again |
| "saved by a newer version of Cellucid" | the session was saved by a later release | open it in that version of Cellucid |
| "not a Cellucid session file of format 4" | the file does not start with the session signature `CELLUCID_SESSION/4` | check that you picked a session file saved by this release |
| Your file's dataset differs, but you want its steps | a session restores only onto its own dataset | replay its steps by name on the dataset shown: **Load… › Replay a session file here…** ({doc}`../m_state_history/06_sharing_and_replay`) |
| "the file is damaged, so it was not opened; the history shown is unchanged." (what is wrong under **Details**, starting *This session file is damaged:*) | the file breaks the format (a manifest, chunk, step or value that does not check out) | obtain a fresh unchanged file from **Share…** |
| "it does not hold together as a history, so it was not opened; the history shown is unchanged." | a saved state or file whose steps, branches or groups break the history's own rules (the rule and the ids under **Details**) | open another saved state, or report it with its **Details**: a history Cellucid wrote should always hold together |
| "Session restored but for … A part of the state this version of Cellucid does not have (…) is left out" | it was saved by a version with a part of the state this one lacks; the rest is restored | nothing to fix here; open it in the version of Cellucid it was saved with to see that part |
| Gzip, byte-length, truncation, or trailing-data error | corrupt or dishonest file | transfer/download the original again; do not edit it |
| Progress disappears after Cancel or a newer load | intentional stop | no action; the older load neither succeeded nor failed |
| **Session fully restored.**, but the view was already similar | a complete load with a subtle visual change | compare the high-signal state below |
| Load is slow | many steps, large highlight groups, fields you made or analysis results | leave progress open or Cancel; close memory-heavy tabs if needed |
| "Continue where you left off?" did not appear | the dataset has no history kept in this browser holding something of yours, **Start fresh** was chosen for it, or this browser no longer holds it (a private window closed since, cleared site data) | open it from **State › Session › Load…** if it is listed; see {doc}`../m_state_history/03_saved_states_and_autosave` |

## Know the success boundary

A progress card, or a step that is already shown, is not success. Cellucid
reports **Session fully restored.** only after:

1. the session's dataset is open with the cells it was saved on;
2. the signature and manifest check out;
3. the dataset fingerprint and the fields and genes the shown step uses check
   out;
4. the step shown is put back; and
5. every other step's values are read and checked.

That notice, **Session fully restored.**, is the load's one report of
success. If any step fails, what was shown before comes back, the trackers you
had stay, and the load's progress notification ends as **Loading failed**
with the reason.

## Save

### The button produces no file

1. Wait until the dataset and controls are ready; **Share…** stays grayed out
   until a dataset is open.
2. Click **Share…** under **State › Session**, then **Save file** again. When
   no branch is ticked, **Save file** is grayed out and the dialog says
   *Choose a branch to share.*
3. Check for a `Sharing the history failed: …` notification, and the browser
   Console for `[state] Sharing the history failed`.
4. Check whether your browser or organization blocked downloads from the site
   or the notebook output the viewer is shown in.
5. Open the same application URL in a normal tab and retry if an embedding
   policy blocks downloads.

If another dataset is opened, or the data mapping or the open cells change,
while a save is being written, the save fails (the Console says "The dataset
changed while the session was being saved"): save again once the change is
done.

Chrome, Firefox, Safari/WebKit, and Chromium-based browsers all use the same
writer. A browser-specific failure is reportable; changing the file extension
is not a fix.

### Saving takes a long time or the file is large

The likely causes are highlight groups with millions of cells, many fields you
made, or many analysis results, in any step of any tracker. To make the file
smaller:

- delete highlight groups and analyses you no longer need, and the steps that
  made them;
- remove trackers that hold work you do not want to send; and
- avoid saving near-identical files again and again; this browser already
  keeps every step.

Do not remove parts from a saved file. A file edited by hand is refused.

## Load picker

### Picker does not appear

The picker must be opened by a direct user gesture. Close any covering dialog,
click **Load… › Open a session file…** once, and inspect Console if nothing happens. Managed browsers
or embedded frames can apply file-input policies; use a standalone tab with the
same application build when permitted.

### File is hidden or cannot be chosen

The file must end exactly in `.cellucid-session`. Unzip an archive first and
remove accidental suffixes such as `.txt`. Renaming an arbitrary file to the
extension does not convert it.

### You canceled the picker

Closing the operating-system picker without choosing a file is a normal no-op.
It starts no load and should show no red error.

## Complete identity refusal

The fingerprint has exactly:

- `datasetId`;
- `name` (shown in messages, never compared);
- `cellCount`;
- `varCount`;
- `cellOrder`, itself exactly the `embedding` the digest was taken over, its
  `dimension` (1, 2, or 3), and a 16-character hexadecimal `digest` of that
  embedding's stored coordinates, each cell by its position in the dataset's
  own order; and
- `dataMapping`, the stored embeddings and arrays the dataset was shown with.

All but `name` must match, except that the data mapping may differ when this
dataset can show or use what the session names (its embeddings, graph,
modality, gene names, matrix, …): those are put back. How the dataset was
opened does not matter, but a file that states no id (`uns['cellucid']['id']`)
takes one from its format, its file name and its cell and gene counts and gene
names, so the same data opened as a renamed file can carry a different id.
Open the dataset the message names, then retry the unchanged file. Cellucid
does not put back the layout alone or skip filters and highlights.

Read the message before changing anything; each cause has its own:

- **"saved on a different dataset than the one that is open now"**, naming
  both cell and gene counts: the sizes themselves differ. Open the intended
  dataset.
- **"saved on a different dataset. It has the same … cells and … genes"**,
  naming both datasets by the names they show: the sizes agree exactly and the
  identity does not. This is two datasets, or one dataset reached two ways
  under two ids (a renamed copy of a file without its own id, for example).
  Open the dataset the session was saved on.
- **"No dataset is open"**: open the dataset the session was saved on (the
  message names it when it can), then load the session again.
- **"This session shows …, which this dataset does not have"**: the
  session turns on an embedding this dataset lacks (or one it cannot plot,
  with the reason). A session whose embeddings this dataset *does* have is not
  refused: they are put back. Open the dataset version the session was saved
  on.
- **"This session uses … as its neighbor graph, which this dataset does not
  have"** (or modality, gene names, expression matrix, latent space, vector
  fields): the session names an array this dataset does not have or cannot
  use. One it has is put back instead of refused.
- **"recorded its cell order on … which this dataset does not have"**: the
  embedding the session was saved in is not part of this dataset.
- **"the same name and the same number of cells and genes, but the cell
  coordinates differ"**: the digest no longer matches. A one-way digest cannot
  say *why*, so the message offers both causes it could be: the cells are
  listed in a different order, or the same cells were exported again from a
  re-computed embedding. Load the version the session was saved on, or
  re-create the selections on this one.
- **"This session uses the observation fields … and the genes …, which the
  dataset that is open now does not have"**: every key the session's shown
  step needs and the dataset lacks is named. Open the dataset version the
  session was saved on.

```{figure} ../../../_static/screenshots/sessions_sharing/refusal-01-different-dataset.png
:alt: A notification titled "Loading failed" reading "Session file: This session was saved on a different dataset than the one that is open now (3,696 cells and 3,753 genes when it was saved, 219,731 cells and 5,754 genes now). Open the dataset the session was saved on, then load the session again."
:width: 760px

The message as it appears: a session saved on the Pancreas sample opened as
`pancreas.h5ad` from this computer, loaded while the Garcia sample is open. (A
session saved on a sample, a link or a server opens its own dataset first
instead, and so does a local file this browser kept, so it is not refused.) It
names **both** sizes, what the file was saved against and what is open now, so
you can tell at a glance which of the two is the one you meant to open.
```

If the sender changed fields, categories, or gene order while reusing the same
lightweight identity and sizes, publish a new dataset id and save a new
session: those are not covered by the digest. Changed cell order or changed
coordinates are caught for you. See
{doc}`07_versioning_compatibility_and_dataset_identity`.

## Format refusals

The reader accepts session format 4, exactly as the current **Share…** writes
it. A file of another format version is refused naming its version; for format
3 and earlier:

> This session was saved by an earlier version of Cellucid (session format 3), which this version no longer reads: a session file now holds the history of every step (session format 4). It was not opened. Open it in the version of Cellucid that saved it, or re-create the session in this version.

A format-4 file that breaks the format is refused: the notification says the
file is damaged and that the history shown is unchanged, and its **Details**
hold a sentence starting *This session file is damaged:* that names what is
wrong, for example:

- a manifest key missing or extra;
- a chunk listed twice, a chunk no tracker owns, or an eager chunk after a lazy
  one;
- an array that does not match its digest, or a value that does not match its
  hash;
- a step tree that does not hold together (one root, parents and children
  agreeing); or
- a step referring to a value the file lacks.

A file cut short, or with bytes after its last chunk, is refused as an invalid
session file. Obtain a fresh file from **Share…** in the current app.

### You picked a file that is not a session at all

The commonest version of this is not a format problem: it is choosing the wrong
file. Session files are binary and their names are easy to confuse with the
figures and exports saved beside them.

```{figure} ../../../_static/screenshots/sessions_sharing/refusal-02-not-a-session-file.png
:alt: A Cellucid notification reading "This file is not a Cellucid session file of format 4 (it does not start with CELLUCID_SESSION/4).", with a dismiss button at its right, after a file that is not a session file was chosen in Load.
:width: 760px

`This file is not a Cellucid session file of format 4 (it does not start with
CELLUCID_SESSION/4).` means the file does not begin with the session file's
signature, so it is not a `.cellucid-session` file, whatever it is named.
Check that you picked the right file before investigating anything else.
```

### Did you download a catalog's starting view?

That file is a SHA-256-pinned starting view: one tracker of one root, applied
by itself when its dataset is opened from the catalog that names it. You do
not need to load it; loaded by hand, its one tracker, **Starting view**, is
added beside yours.

Choose the dataset from the catalog, wait for its starting view, then use
**Share…** to write your own session file. See {doc}`04_official_sample_states`.

## Cancel and replacement

Pressing **Cancel** on the progress card, or starting another load, stops the
older load. Its progress is dismissed, what was shown before comes back, and it
produces neither false success nor a red failure.

If you canceled accidentally, wait for any newer load to finish and start the
one you wanted again. If a stop appears as a failure without a cancel or a
newer load, collect the diagnostics below.

## "It loaded, but looks different"

First require **Session fully restored.** Then compare:

1. app build shown in the footer;
2. dataset route, id, cell count, and gene count;
3. **State › History**: the tracker in use (the button at the top of the
   branch picker), the branch and the step shown (the sender may have shared
   another branch, the current state only, or a step they did not mean);
4. kept views, their layout, and the focused view;
5. every view's embedding, camera, and navigation mode, and whether cameras
   are linked;
6. active field and category legend;
7. filters, including latent-space outlier filtering and filter by cell count;
8. highlight pages, group names, counts, and visibility;
9. camera path keyframes, **Loop playback**, **Autoplay**, and actual playback
   state;
10. the analysis panel, analysis windows and their results; and
11. Figure Export settings.

Different viewport sizes can clamp floating panels so they remain reachable.
That is expected. Missing fields, highlight cells, or analysis results are
not: the load should refuse and put back what was shown rather than report
success.

## Progress is slow

The loader shows the step that was shown first and reads the rest while it is
shown, and it stays responsive so it can handle **Cancel** and a newer load.
Very large files still take CPU, network, and memory.

- Keep the progress card open if you want the load to finish.
- Avoid launching another large analysis at the same time.
- Close unrelated memory-heavy tabs.
- Press **Cancel** if you do not want to wait; do not force-close the page
  unless the browser itself is unresponsive.

If progress never settles and **Cancel** also cannot run, capture the exact
app build and browser information below.

## Reporting a reproducible bug

Please include:

1. footer build version;
2. browser name/version, operating system, and device architecture;
3. how the dataset was opened, exact dataset id, cell count, and gene count;
4. whether the operation was **Share…**, opening a session file, opening a saved state, a replay, or a
   catalog's starting view;
5. the complete visible error and relevant Console stack;
6. whether **Cancel** or a newer load was involved;
7. the session file and the unchanged dataset file, if sharing is allowed;
   and
8. the smallest reliable reproduction sequence.

Do not post sensitive sessions publicly. See
{doc}`08_security_privacy_and_trust`.

## Next steps

- {doc}`09_edge_cases` (strict expected behavior)
- {doc}`12_reference` (exact file format)
