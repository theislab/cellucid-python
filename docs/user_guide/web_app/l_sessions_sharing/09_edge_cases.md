# Edge cases: sessions

**Audience:** power users and anyone building repeatable workflows

**Time:** 15–30 minutes
**What you'll learn:**
- which unusual outcomes are expected under the strict current reader;
- how cancellation, subsamples, large files and different screens behave; and
- when to save a fresh session.

## Different identity, even with similar files

A session saved on a file opened from GitHub loads on the same file opened
with **Open…** (or a drop) from your computer, from a web address, or from a
server: the route is not part of a dataset's identity, as long as the file
states its own id in `uns['cellucid']['id']`. A file without one gets an id
from how it was opened (its file name, its link), and then the route does
matter. If a publication changes its dataset id, `datasetId` differs, and that
refuses the whole load.

This is intentional. Cellucid will not put back the layout alone or skip the
parts that name cells. Load the dataset identity used by the sender, or save a
fresh session after loading it.

See {doc}`07_versioning_compatibility_and_dataset_identity`.

## Same fingerprint, changed content

The fingerprint checks dataset id, cell count, gene count, the data mapping,
and a `cellOrder` record: a digest of the stored coordinates of the embedding
the live view showed, each cell by its position in the dataset's own order,
plus that embedding and its dimension. It is not a hash of expression values,
observation or gene fields, category inventories, cell identifiers, or the
embeddings that were not displayed.

So the remaining gap is *values*, not order: republishing changed field values,
renamed or re-ordered categories, or a different gene table under the same
identity and the same sizes, while leaving the displayed coordinates untouched,
still passes the guard. Re-ordering cells or recomputing the embedding does
not: both change the digest and are refused with their own message.

Never reuse an identity for changed content. Publish the changed file under a
new immutable id and save a new session.

## The step shown appears before the load is done

A load shows the step the file showed as soon as that step's values are read,
then reads the values of every other step while it is shown. Progress can
therefore stay active after the view has changed; that is not partial
success. **Session fully restored.** appears only once everything is read and
checked, and a new step you take meanwhile waits until then.

If the file's structure is invalid (a value a step names is missing, a
conversation of a tracker it does not hold), the camera, fields, filters,
highlights, windows and playback go back to what was shown before, and the
trackers you had stay. A part that only cannot be put back (a value this
version refuses, one saved by an earlier version of Cellucid) does not undo
the rest: it is left as it was shown and named in **Session restored but for
N parts**, and the history records it as a step of its own
({doc}`../m_state_history/03_saved_states_and_autosave`).
If you press **Cancel** or start another load, the same happens, and the older
progress is dismissed without success or a red failure.

## Loading adds trackers; opening a saved state replaces the history

**Load… › Open a session file…** adds the file's trackers beside yours (a name you already use is
named apart, `"(shared)"`) and leaves your trackers exactly as they were; the
notice's **Undo**, or **Restore** on "Loaded … from <file>" in **⋯ › Recently
deleted…**, moves them to **Recently deleted**. A file of the steps only opens as
a replay, a tracker of its own. Opening one of the saved states in **Load…**
replaces the history shown: with **Autosave in this browser** on, the history you had
stays in its own saved state, and what you opened is kept as another, so
neither is written over. With it off, share a file first if you want to keep
what you had.

## Empty parts replace existing ones

Every part is put back exactly as it was at the step shown, an empty one
included:

- an empty camera path clears the camera path on screen;
- no analysis windows in the file closes the windows that are open; and
- no highlights in the file removes the highlight groups on screen.

A load does not keep stale state from whatever you explored before. Saved
playback settings come back exactly; a camera path saved with **Autoplay** on
starts playing once the load is done.

## Cells to open is a step

Changing **Cells to open** (or **Genes read for**) is a step of the history.
**Undo** reads the dataset again with the cells before, and showing a step of
another subsample reads those cells. The same dataset read again on other cells
keeps its history and puts no starting view back.

A session saved on a subsample records its exact cells. Loading it opens the
dataset with those cells, reading it again if it is open with others, so the
cell-order digest is taken over the same cells it was saved on
({doc}`../b_data_loading/21_subsampled_loading`). A catalog's starting view is
the one exception: it applies to any subsample of its dataset.

## Analysis results in sessions

Differential expression, marker-gene, signature, correlation, dot plot and
composition results are saved with their analysis (in the sidebar panel or a
floating window) and come back without being run again. A result is saved with
what it was computed on: it comes back current, or marked **Stale · Run** when
that has changed on the page it is loaded into. A session is still not the
canonical scientific-results archive: export result tables for publication or
downstream analysis. A saved result that names a gene or category the open
dataset lacks is refused.

## Different screens and browsers

Floating panels are clamped so they remain reachable on the current viewport.
A layout saved on a large monitor can therefore move slightly on a laptop or
mobile-sized window. Reposition the panels and save again if exact placement
matters.

A load does not promise canvas focus, hover state, or pointer-lock state.
Click the canvas once before using keyboard shortcuts.

## Large session files

Session files grow with:

- many highlight groups containing millions of cells;
- many fields you made;
- many kept views; and
- many saved analysis results.

Each step adds only the parts it changed, and an array (a highlight group's
cells, a field's values) is stored once however many steps keep it. The loader
enforces declared
sizes, checks each chunk's gzip stream before inflating it, and stays
responsive so **Cancel** and newer loads still work. You can keep files
lighter by deleting highlight groups and analyses you no longer need, and by
removing trackers that hold work you do not want to send.

## Missing field or changed schema

A matching fingerprint does not make a missing field acceptable. A field or
gene the shown step uses that the dataset lacks refuses the load before
anything changes, naming every missing key. A field that only another step
uses is checked when you show that step: the step is refused, naming the
field, and the step shown stays. Cellucid does not quietly keep the current
field instead.

Load the unchanged dataset file, or set up the changed dataset and use
**Share…** to write a new session. To bring the work across, replay the old
file's steps on the changed dataset (**Load… › Replay a session file here…**,
{doc}`../m_state_history/06_sharing_and_replay`).

## A file from another version of Cellucid

This release reads session format 4 only. A file saved by an earlier Cellucid
(format 3 or older, which held one saved state rather than the history) or by
a newer one is refused before anything is read, naming its format. A format-4
file from another version is restored part by part: a part of the state this
version does not have is left out of every step, and a history (tracker) whose
records this version cannot read is left out, each named in the notice that
ends the load (**Session restored but for …**); everything else is restored.
Open the file in the version of Cellucid that saved it to see those parts.

## Catalog starting views and manual sessions

A catalog's starting view is a format-4 file of one tracker of one root, with
no cells to open and no layout, applied only when its dataset is opened from
that catalog. There is no need to download it and load it; loaded by hand it
is added beside your trackers as one tracker, **Starting view**.

Choose the dataset from the catalog first, then use **Share…** if you need a
shareable session. See {doc}`04_official_sample_states`.

## Next steps

- {doc}`10_troubleshooting_sessions` (symptom → diagnosis → fix)
- {doc}`12_reference` (the exact file format)
