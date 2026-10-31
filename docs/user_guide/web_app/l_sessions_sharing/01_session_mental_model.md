# Session mental model

**Audience:** everyone (wet lab + computational + power users)  
**Time:** 15–30 minutes  
**What you’ll learn:**
- What a Cellucid “session” is (and what it is not)
- The three persistence layers: dataset vs session bundle vs browser preferences
- Why restore is all-or-nothing under the current format
- How to choose between sharing a link, sharing a session, and sharing the dataset file

**Prerequisites:**
- A dataset loaded in the web app (demo datasets are fine)

---

## Mental model (one sentence)

A Cellucid session is a **portable snapshot of the app’s state** (camera, fields, filters, highlights, views, and some analysis UI) saved as a single `.cellucid-session` file.

---

## What a session is (plain language)

If you’re not technical, think:

- **Dataset** = the “slide under the microscope” (the actual data).
- **Session** = “how you set up the microscope” (where you looked, what you colored by, what you highlighted).

When you click **Save** under **Session file**, Cellucid downloads a file like:
`cellucid-session-2025-12-30T02-00-14.cellucid-session`

When you click **Load**, Cellucid uses that file to restore the UI back to the saved state.

:::{important}
A session bundle does **not** contain the dataset itself.

It can restore state only after you load the dataset it was created from (or a dataset Cellucid considers the same).
:::

---

## The three persistence layers (this prevents most confusion)

Cellucid uses multiple “memory” layers. They are easy to mix up:

1) **In-memory state (volatile)**
   - Lives only while the page is open.
   - Lost if you refresh, crash, or close the tab.

2) **Session bundle (`.cellucid-session`)**
   - Explicit file you download via **Save**.
   - Intended for:
     - reopening later,
     - reproducibility (“these exact views/filters/highlights”),
     - sharing with collaborators.

3) **Browser preferences (localStorage, small UI conveniences)**
   - Things like theme/background may persist across refreshes.
   - Not a reproducibility mechanism.
   - Not meant to be shared.

Practical implication:
- If you want to **share** or **reproduce**, rely on **session bundles** (and the dataset file).
- If you only want the UI to “feel the same” next time on your machine, preferences are enough.

---

## Fast path (wet lab / non-technical): save your work and send it to someone

1) Load your dataset.
2) Explore until you have what you want (color-by, filters, highlights, snapshots).
3) Open **Session** → **Session file**.
4) Click **Save**.
5) Find the downloaded file (usually in your Downloads folder).
6) Send the file to your collaborator. A session saved on a sample, a data link
   or a server records where the dataset is and opens it; for a dataset opened
   from your own computer, send the dataset file with it.

To reopen later, click **Load** and pick your `.cellucid-session` file.
Its dataset opens first if it is not open already; a session saved on a local
file needs that file opened again before it is loaded.

---

## Practical path (computational users): what “state” actually means

Cellucid has **state at different scopes**:

- **Global-ish state** (affects what you see everywhere):
  - some UI controls,
  - some registries (renames/deletes, user-defined field definitions).
- **Per-view state** (live view and each snapshot has its own):
  - camera and navigation mode,
  - which embedding it shows (its dimension follows),
  - which field is active (obs vs var),
  - which filters are active,
  - view layout settings.

Sessions exist to capture and restore those scopes in a predictable way.

### Ordered restore with one terminal outcome

Sessions use two scheduling classes:

- **Eager chunks** contain the state needed by later chunks, such as field
  definitions, active views, and highlight metadata.
- **Lazy chunks** contain heavier arrays, such as highlight membership and
  user-defined field codes.

This ordering keeps the browser responsive, but it is **not partial success**.
The public load operation waits for both classes, prepares and commits every
feature owner, refreshes the final UI, and only then reports **Session fully
restored**. If any later chunk fails, Cellucid rolls back the earlier changes.

Wait for the terminal success notification before treating the restored view as
ready. Pressing **Cancel**, closing its progress card, or starting a newer
restore aborts the older operation without publishing false success.

### Dataset identity (the first restore check)

Each session carries a **dataset fingerprint**: the dataset id (and its
display name, used only in messages), cell count, variable count, a
**cell-order record**, and the **data mapping** — which stored embeddings and
arrays the dataset was shown with. Beside it the session lists every
observation field and gene it uses, so a dataset missing one is refused before
anything changes, naming every missing key.

The first three say *which* dataset. The cell-order record says *which ordering
of its cells*, which the first three cannot: a session stores your selections by
each cell's position in the dataset's own order, and a dataset republished under
the same name with the same number of cells and genes but its rows in a
different order matches all three scalars while every position now points at a
different cell. The record is a short digest of the stored coordinates the
session was saved against, each cell named by its position, plus the embedding
(and its dimension) the digest was taken over — a digest only means anything
within one embedding. Cellucid checks the digest on that embedding whichever
one is on screen when you load; how the cells were opened (all of them, a
share, grown step by step) does not change it.

If any part of the fingerprint does not match the currently loaded dataset:
- Cellucid rejects the complete restore;
- no layout-only subset is accepted; and
- the existing app state remains intact or is recovered by rollback.

You will see one of these messages, each naming what actually happened:

| What you see | What it means | What to do |
|---|---|---|
| The session was saved on a different dataset (with both cell and gene counts named) | it is not this dataset | open the dataset the session was saved on |
| The session was saved on a different copy of this dataset (with both identities named) | the sizes agree exactly and the publication identity does not — one dataset reached two ways, such as a locally served copy of a published export | open the dataset the session was saved on |
| The session uses an embedding (or neighbour graph, modality, gene names, expression matrix, expression transform, highly variable genes, latent space, or vector field) this dataset does not have or cannot use, naming it | it was saved on another version of the dataset (a saved choice this dataset does have is restored, not refused) | open the dataset version the session was saved on |
| The dataset has the same name and counts but its cell coordinates differ | either the cells are listed in a different order, or the same cells were exported again from a re-computed embedding — a one-way digest cannot say which | load the version the session was saved on, or re-create the selections here |
| The session uses observation fields or genes the dataset does not have, naming every one | it was saved on another version of the dataset | open the dataset version the session was saved on |
| The session was saved by an earlier version of Cellucid (session format 2) | its check that it marks the same cells was taken over coordinates as they were drawn, which this version cannot recompute | re-create the session in this version and save it again |
| The file is not a Cellucid session file of format 3 | it does not start with the format-3 signature `CELLUCID_SESSION/3` | check that you picked a session file saved by this release |

The format version is part of every session file's signature, and this
release reads format 3 only.

See {doc}`07_versioning_compatibility_and_dataset_identity` for the exact
matching rules and safe collaboration recipes.

---

## Deep path (power users / developers): what’s in the file (high level)

A `.cellucid-session` is a single file that contains:

- a JSON **manifest** (metadata + a list of chunks), and
- many length-prefixed **chunks** (JSON or binary, often gzip-compressed).

Chunks are owned by “contributors” (features). This keeps sessions modular:
- core state, highlights, analysis windows, user-defined fields, etc.

Sessions are treated as **untrusted input**:
- the loader exact-validates the manifest and every current chunk profile;
- enforces stored, decoded, and gzip structure limits before native
  decompression;
- rejects unknown, missing, duplicate, aliased, or reordered chunks; and
- commits the complete transaction or rolls it back.

If you want the code-level reference, see {doc}`12_reference`.

---

## When to use sessions vs other artifacts

Use this decision rule:

- **Session bundle**: when you want to preserve the interactive state (views, filters, highlights).
- **Exported figure**: when you want a publication artifact that does not depend on the app to interpret.
- **Community annotation**: when you want many people to collaboratively label/vote in a conflict-free way (GitHub-backed).

Sessions are great for:
- “Please look at this exact selection and camera angle.”
- “Here are the highlight groups I think are interesting.”
- “Here’s the state right before I ran DE; can you verify?”

Sessions are not ideal for:
- multi-user, long-running labeling workflows (use community annotation),
- distribution of the dataset itself (share the dataset file or hosting link).

---

## What changes in app state when you load a session

Loading a session is intentionally “strong”:

- It **overwrites** most of your current view state (camera, active fields, filters, etc.).
- It may **cancel** an in-progress session restore if you load another session.
- It explicitly replaces empty state too: an empty Camera Path clears the
  destination's path.

If you want to preserve your current work, save a session *first*, then load the other session.

---

## Interface reference

```{figure} ../../../_static/screenshots/data_loading/data-loading-session-panel.png
:alt: The Session panel with the Pancreatic endocrinogenesis sample loaded: the Dataset box with its name, cells and genes, format and host, description and a collapsed Details; Open data with the Sample dataset dropdown, the Open… button and the Link or server bar; Save and Load under Session file; and the collapsed Data mapping group.
:width: 532px

The Session panel gathers every loading path under **Open data** and keeps the **Save** and **Load** buttons under **Session file**, directly below it.
```

---

## Next steps

- {doc}`02_what_gets_saved_and_restored` (explicit inclusion/exclusion list)
- {doc}`03_save_restore_ux` (**Save** / **Load** step-by-step)
- {doc}`05_share_workflows_links_sessions_data_files` (how to share with collaborators)
- {doc}`10_troubleshooting_sessions` (if anything goes wrong)
