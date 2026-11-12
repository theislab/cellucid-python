# Changelog

All notable changes to Cellucid will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.9.8]
<!-- CELLUCID_VERSION -->

Version 0.9.8 is a patch release on 0.9.7. Files prepared, sessions saved
(format 4) and dataset catalogs written by 0.9.7 are read as before. The
viewer protocol is still schema 5: use the Python and R packages of this
release with the web app of this release. The main changes are:

### Web app

#### Opening and reading data

- **Every request to a data host goes through the one transport; a network
  failure is one only when no answer came.** A person's "Try again" of a
  dataset catalog, the requests to a Cellucid server (health, info,
  catalog) and GEO's check of whether NCBI is failing were sent with a plain
  `fetch`: outside the host's request budget and `Retry-After`, left
  hanging forever by a host that stopped answering, and, when no answer
  came, failed with the browser's bare `TypeError`, said as a problem inside
  Cellucid instead of who is at fault (a "Try again" of a catalog on a host
  without CORS said a different sentence from its first opening). They now
  go through `fetchOnce` (asked once where they were: a server's failures
  stay terminal). A link's opening called any `TypeError` whose message
  named "fetch" or "network" a network failure; it now asks the transport
  (`isUnanswered`). The raw `fetch` calls left (probes that must not count
  in the book of hosts, OAuth token requests, keepalive revocations, the
  model APIs, the page's own files) are each named with their reason, and a
  new one fails the test suite.
- **A Python or R viewer's connection lets go of what ran when its channel
  reopened.** Disconnecting a viewer host kept the callbacks run when its
  dropped channel opened again (it let go of those run at its end), so a
  viewer connected again would also run the old connection's.
- **A sample's start says its name and no download that is none.** Every
  sample's start showed "↓ Loading Starting view 0 B/s Cancel": a starting
  view is a part of its dataset's opening, already in the page, not a
  transfer. It now has no line of its own (its opening card says the
  opening), and a load's line names a speed only once bytes were measured
  arriving, never "0 B/s" before the first bytes or during work that reads
  none. The opening card said the id ("Opening garcia") until it closed; it
  now says the name the catalog gives the dataset as soon as the catalog has
  answered, and the dataset's own name once it is read. A dataset catalog
  (`datasets.json`) entry may carry its dataset's `name`; the samples' catalog
  does, and the picker lists it before the files are read.
- **Cellucid starts where the browser keeps nothing for the site.** With
  site data or cookies blocked (and in some private windows) Cellucid did not
  start at all ("Cellucid could not start"): the theme read `localStorage`,
  whose very read throws there. Every part now reaches the browser's storage
  through one helper (`utils/web-storage.js`); without it, a setting you
  change (theme, background, antialiasing, history settings) holds until the
  page closes, and Cellucid says once, after the welcome, **Nothing is kept**.
  A full store no longer refuses a background or antialiasing change, and a
  stored preference this version does not know starts the page on its
  defaults instead of stopping it.
- **A link Cellucid cannot read is said as the link's.** A parameter given
  twice or empty (`?dataset=a&dataset=b`, `?dataset=`), spaces around a
  value, `anndata=1`, an unknown `source`, a server and data at once, a viewer
  address Python or R did not write, or `source=remote` on a page no server
  serves stopped Cellucid on "Cellucid could not start — Something went wrong
  while starting", **Try again** first (which met the same refusal). Cellucid
  now starts without the link, and the card over the empty viewer names the
  parameter and what is wrong ("This link names "dataset" 2 times; a link
  names it once."), **Open without this link** first; `?source=foo` alone,
  which opened nothing and said nothing, is refused the same way.
- A `dataset` the link names that its source does not hold is named, with
  what the source holds ("This link names the dataset "x", which the server
  this link names does not hold; it holds "a" and "b"."), and a link to a
  source of several datasets says which ones to choose from; both lead with
  **Open without this link**.
- `?remote=http://localhost:8765/` (a trailing slash, an upper-case scheme)
  opens the server, as the address bar does; a `remote` that is no web
  address is refused as the link's.
- A link's `sampleCatalog` that is no folder address, or given twice, is said
  as such; it was "what is published there is not in a format Cellucid can
  read. Write it again with anndata". **Open without this link** drops it.
- A browser with WebGL2 that would not give the viewer a context (a reset
  graphics process, too many pages drawing) is told **The browser's graphics
  stopped working**, close other tabs and reload; it was "Something went
  wrong while starting".
- A link's `annotations=owner/repo` is said in **Community annotation** until
  you are signed in ("This link names the repository …: Connect repo connects
  it, signing in with GitHub first."); it was silent. A value that is no
  repository, and the errors of applying one, name the parameter as the link
  writes it (`annotations`), not `annotation_repo`, which no link has.
- A request asked again after a failure says so on the line that shows its
  progress, in one sentence for every transport: "GitHub is answering with
  errors (HTTP 502, a gateway error); trying again in 1.4 s (attempt 2)…". A
  community annotation Pull against repeated 502s read "Checking the
  repository… 5%" while it retried in silence.
- **A stopped community annotation Pull or Publish is no failure.** Its
  **Cancel** was reported as a failed Pull ("Pull canceled; …" in red), and a
  Pull or Publish stopped because the dataset, repository or GitHub user
  changed, the repository was disconnected or annotation moved to another tab
  failed its card beside the message that said why. The stop is now a
  cancellation wherever it is caught (`community-annotation/sync-stop.js`):
  your Cancel ends the card with its words, any other stop withdraws the card
  (unless GitHub had already accepted part of a Publish: the card then names
  the files it accepted and asks for a Pull before publishing again), and only
  a connection lost mid-sync is reported as failed. An agent's
  **Read the community annotation** is told the Pull was stopped.
- **Checking your role on an unreachable GitHub is asked once, not forever.**
  Each request was already asked again while that could help; on top of that,
  the whole check was asked again on a timer of its own without end ("Cellucid
  asks again in 37 s"). It now says GitHub could not be reached and checks
  again with the next Pull (yours or an automatic one) and as soon as this
  computer is back online. A fault in Cellucid's own code during the check (a
  `TypeError`) was read as "GitHub unreachable" and retried; it is now
  reported as what it is.
- While a large dataset opens progressively (Suo, the default sample: a
  first view, then the rest of its 561,947 cells in a few steps), the
  collapsed Subsampling header says the share it grows to and the cells read
  so far ("100 % cells · reading 197K"), as the Subsampling control does; it
  said the share open so far ("35 % cells"), which read as a subsample
  someone had chosen.
- **Opening another share of the cells with connectivity edges shown works.**
  With **Show connectivity edges** on, **Open data › a share of the cells**
  (or every cell again) failed half done: "could not be read again with these
  cells", the highlights out of step with the cells on every frame, and every
  later opening failing the same way. The neighbor graph shown belongs to the
  earlier cells; it is now let go in the same step the cells change, and read
  again on the new ones, with the edges still shown and your **Max edges**
  kept (it also went back to its default on every re-read of the graph).
- The newest dataset you ask for is the one that opens. A link typed while
  the default sample was still being chosen lost to the sample, and one typed
  while a large sample (Suo) was still growing in place failed with "Could
  not open the dataset" ("Dataset selection stage was superseded before
  publication") while the sample stayed: the startup's own opening counted as
  a choice of yours, and each growth step of the dataset shown counted as
  another dataset chosen.
- A category with no name is shown as *(empty)* (an empty name) or *(blank)*
  (spaces alone), quietly, and a name with spaces at an end in quotes, by one
  rule everywhere a category is named: legend, hover card, filters, menus,
  Group categories, figure legends, analysis groups and agent reads (R9D-8).
  Two such rows were indistinguishable blanks. The stored value is unchanged,
  and a step can name the category as shown.
- The hover card says *Missing* for a cell with no value, as the legend and
  the filters do, for a number too (it said `NaN`, and *missing* in lower
  case for a category) (R9D-9).
- Every CSV the app writes leaves a missing value as an empty cell: the
  composition and dot plot tables wrote `NA` where every other table writes
  nothing, and one quoting rule now serves every CSV (a carriage return in a
  label is quoted too).
- **A notice about what was open goes when the next dataset opens.** A file
  refused before it could open ("Cannot open this file"), a sample view that
  could not be applied and a control's failure (a field that could not be
  read) were kept until dismissed, so they stayed beside the next datasets as
  if about them (two samples in a row left two "Sample view unavailable"
  cards). They are now withdrawn when the next opening begins, as an
  opening's own failure already was: one list of such notices
  (`untilNextOpening`) for the whole page.
- A gene read waiting (offline, a host's request budget) is said once, by the
  gene's own notification with its Cancel: a second "Reading genes · waiting:
  This computer is offline" card no longer said it again beside it.
- An action still running when another dataset opens ends as its cancellation
  ("… was canceled before it finished; nothing was kept."): a markers run
  stopped by opening Pancreas ended with "The history has no value
  c7be2ab0cc2b", because its change was put back through the new dataset's
  history instead of the one that recorded it; a change ending after the
  switch records nothing in the new dataset's history.
- An opening says what opened: a file "Opened <name>: 1.2K cells" (it said
  "User data ready: 1.2K cells", naming nothing, so two files opened in a row
  left two identical cards), a server's one dataset "Opened <name>", and
  "Connected: 3 datasets found…" (no hyphen as a dash). The Dataset box's
  Kept row of a file read whole names the control that removes it, **Session ›
  Stored › Manage…** (it named a Clear button that no longer exists).
- An obs column of 64-bit integers holding whole numbers past 2^53 (database
  ids, barcodes hashed to integers) is read as the exact text of each value,
  a text column, in every format (AnnData `.h5ad` and Zarr, plain or
  nullable, MuData modalities, Loom); it was rounded to the nearest double
  without a word (F116). The column is text or numbers by every cell it
  stores, whatever cells are open (a subsample scans the stored column once
  when the field is first shown). Listed with the numbers until its first
  load finds such a value, it then moves to the categorical fields, and its legend says
  so in one line ("IDs past 2^53, read as text."), the reason on the line's
  tooltip. Columns whose integers are all within 2^53 stay exact numbers. A
  text table beside a matrix holding such a number is a text column too, and
  an R `bit64::integer64` column is kept as 64-bit integers instead of
  doubles.
- The portable element-name rule (obs/var columns, obsm keys) takes the
  longest file name every platform takes: at most 255 bytes in UTF-8. A
  longer name was written into a Zarr folder that failed with
  `OSError: File name too long` (or into a `.zarr.zip` that cannot be
  unpacked on Linux); it is refused with the rule's sentence.
- A name a writer of the web app did not choose is written in its portable
  form and said, never dropped or nested: a converted MatrixMarket folder's
  metadata column `CD4/CD8` (and a column its table repeats) was left out
  without a word, and is now `CD4_CD8` (`dup_2`), a column an earlier table
  already gave is noted; a copy of an `.h5ad` whose obs column is `a/b`
  ("Save AnnData with this state") wrote a nested group anndata could not
  read, and now writes `a_b`; the highlighted cells' `.h5ad` follows the same
  rule (letter case and length included). One rule, `portableElementName`.

- The five built-in samples are rebuilt by their pinned builds with the 0.9.8
  writer: their stores now hold every chunk, so a web host no longer answers
  404 for the priority seed and single-category columns at every opening.
  Each records its cells' total counts (`obs['total_counts']`, summed over
  the whole source before any gene was selected: Pancreas's 27,998 genes, the
  atlases' 8,192 integration genes; `uns['cellucid']['totalCounts']`), so
  normalizing a sample's genes reads one column instead of every gene.
- The opening's transfer rate leaves out the time it waits on you (an `.rda`
  asking which object to open, a sign-in): beside the question it fell for as
  long as you read it, and it stayed low after you chose. No rate is shown
  while you are asked.
- A host that stops answering is said for what it is, without its address: "…
  stopped answering this site, which it answered before: it may be refusing
  requests from this network for now, or failing, or this computer lost its
  connection" (a sample's listing said only "Try again").
- **A web page in place of data is said for what it is.** A Zarr folder on a
  host that answers a path it has nothing at with its site's page (`200
  text/html`: Cloudflare Pages without a `404.html`, Netlify or Vercel
  rewrites) opens: each chunk zarr-python left out reads as its fill value, as
  a 404 says, where every field, gene or embedding with one failed as "may be
  incomplete or damaged". A link answered with a sign-in or preview page, or a
  page saved under a data file's name, is said to be a web page, with "open
  the file's direct download link" (it was "not a valid HDF5 file … save it
  again with anndata"); a catalog's `datasets.json` answered with a page is
  not found, not a network failure to try again.
- **A link is read in the form it is pasted.** An address without `https://`
  (`example.org/atlas.h5ad`, `zenodo.org/records/1`; this computer's
  `127.0.0.1:8765/…` and `lab-server:8080/…` at `http://`; the address bar
  took `127.0.0.1:8000/atlas.h5ad` for a server), one between quotes, angle brackets or
  backticks, a GitHub repository's own page or clone address
  (`github.com/owner/repo`, `….git`, `git@github.com:owner/repo.git`), a file
  page's `?raw=true`, a Zenodo file's preview page, `hf.co`, and an address
  inside a Zarr store folder (its `.zmetadata`, a group) open what they name.
  `example.org/atlas.h5ad` was read as the GitHub repository `example.org`, a
  repository's page was refused as "not a file Cellucid opens", and a preview
  page as private. A path on this computer (`/Users/…`, `C:\…`, `file:`) is
  refused saying to use Open… › File…, and an `ftp://` address saying no page
  can read one, both where "does not recognize" was said.
- Each growth step of a Zarr folder on a web host reads at once: the check
  that the folder is unchanged (a `HEAD` of its `.zmetadata`) is asked beside
  the step's reads, and a step read from a folder found changed is discarded
  and the folder opened anew. Each step waited one round trip for the check
  first (about 4 % of an opening on a 3G link).
- A viewer stopped from Python or R (`viewer.stop()`, `cellucid_stop()`, a
  notebook cell run again) says **Viewer stopped** on its frozen last view, as
  a viewer whose server went away now does too; it froze without a word (and
  the notice of a server gone was drawn under the frozen picture). Offers the
  frozen page can no longer take (Continue where you left off) go. A failed
  notebook command no longer ends its reason twice (`"sample".. The viewer is
  unchanged`).
- A Zarr folder or a dataset catalog on a host that does not let other
  websites read it is said so, as a file on it is ("… does not let other
  websites read its files … Download the folder and open it (Open… ›
  Folder…)"); it said "could not be opened. Try again; check the address",
  and a catalog that the host "stopped answering". A host that lets this site
  read one file and not another (CORS set per folder) is said to refuse that
  file, not to refuse this network for now.
- My datasets forgets a kept file or folder only when it is gone (moved,
  renamed, deleted) or changed: one in use by another program, or a cloud
  folder's file not downloaded yet, is said to be unreadable for now and kept
  for the next try (any read failure said "moved, renamed or deleted" and
  forgot it).
- A Zenodo record, or the CELLxGENE Census release list, refused with 401 or
  403 is said as the host refusing this network for now (both are public to
  everyone); a Zenodo record said "could not be read (HTTP 403). Check the
  record number".
- A link whose host sends part of a file and stops is said so ("… started
  sending atlas.h5ad and stopped (no more bytes for 3.0 s), also when asked
  again: the host, or the connection to it, is failing or very slow"): it was
  taken for a host that content-codes its files and downloaded whole, which
  stalled again, or said not to let other websites read it.
- A host that asks for a sign-in (HTTP 401) for one file is said to refuse
  that file (private, or an expired link), also after it answered this site
  for another; it was said to refuse "this network for now".
- **A file whose AnnData parts break the format is explained**, at opening and
  at the read of the field or gene concerned: "This .h5ad file breaks the
  AnnData format: Categorical code 5 for 'cell_type' is outside
  missing-or-category bounds. Save it again with anndata …" (obs without its
  index, codes past their categories, sparse pointers past the stored values,
  a matrix whose shape disagrees with obs and var, …). It said only "Could not
  open …", with the rule under Details.
- **A file is refused for what it is, not for an embedding.** A Loom file
  without a two-dimensional `/matrix`, or a 10x matrix missing `data`,
  `indices` or `indptr`, was read as AnnData and refused as "has no cell
  embedding"; it is now refused as a Loom (10x) file not laid out as its
  format requires, naming the missing part, whatever its name. An HDF5 file
  or Zarr store of none of the formats Cellucid reads is said to be that,
  listing them, and an AnnData store without obs or var is refused for them
  before its embedding is looked for.
- A damaged or truncated R file (`.rds`, `.rda`; gzip, bzip2, xz) is explained
  with what to do, as an `.h5ad` or Zarr store is; it said only "Could not
  open …". A truncated HDF5 file names both sizes in bytes.
- **A dataset with no starting view** (your own file, a link, a served
  object) opens colored by its best categorical field, in a ranked order: a
  cell type, another annotation, a lineage, clusters, the tissue, a
  condition, the cell-cycle phase; never a sample, batch, id or quality flag;
  else the first field with 2 to 72 categories. The same ranks pick the
  default field of Composition, Dot plot, Marker genes (**Group by**) and
  Differential expression (**Field**), the field a subsample is judged by,
  and the agent chat's starters (on Suo they offered "Annotate the clusters
  of sample_ID"; now cell_type) (`app/utils/field-roles.js`).
- Opening a sample while **Volumetric smoke** is on applies its starting
  view: it failed with "Sample view unavailable" and the dataset opened
  uncolored. A dataset whose starting view fails still opens colored by its
  best field.
- The Suo sample opens on a view of its whole cloud: its starting view had
  the camera inside the cloud, a close-up of giant points.
- A Zarr store on a web host is asked for an object it answered as absent (a
  chunk holding only the fill value, which zarr-python leaves out) once per
  version of the store in a tab, not at every read: the default sample asked
  for its priority seed 6 to 11 times at start, once per growth step.
- The sample picker reads the samples' files for their names and counts only
  once the startup dataset has opened: the listing began a moment before the
  opening and read every sample's file, the one being opened too, against its
  first view. Closed, the picker names the open dataset alone (the Dataset box
  says its cells); opened, it shows the counts.
- When a sample cannot be listed, the notice names it and says why
  ("Sample datasets could not be loaded: “garcia” could not be listed. Try
  again; …"), not "could not be opened" with the store's address; the
  address stays under **Details**.
- Choosing another dataset while a Subsampling **Apply** reads stops that
  reading at once (on a slow link it read the whole first-view embedding
  first), and the new opening is shown on its own.
- A server on this machine (`cellucid serve`, a local web server) is never
  paced as a host that blocks networks: a restart of it no longer slowed every
  later visit for a day.
- Nor is a host paced when your own connection drops for a moment: requests
  lost on this side (this computer offline, or failing sooner than the host
  ever answers) teach it nothing. A 3-second drop while Suo opened from GitHub
  Pages slowed the rest of the opening to 70 KB/s over a fast link ("limits how
  many requests it answers: the next go in 2.4 s") and kept that pace for a
  day; each step now reads in under a second once the link is back.
- While this computer is offline, an opening or read in progress says "This
  computer is offline: reading goes on as soon as it is back; what was read is
  kept." (it said the host "may be refusing requests from this network") and
  goes on the moment the connection returns: after an 18-second outage Suo
  waited up to 20 s more for a back-off meant for a refusing host.
- A refused file's notice goes when the next opening begins; it stayed
  through every later opening, reading as if the dataset shown had failed.
- **A file with a sparse `obsm` array can be subsampled** (R9M-8): every
  share of such a file (an scvi-tools or converter output with a CSR array
  in `obsm`) failed with "could not be read again with these cells", because
  listing the array as one that cannot be drawn read row pointers of the
  share that were not read yet.
- **A reading of the cells that fails is a failure** (R9M-9): Subsampling's
  Apply, an agent's `cells_open` or a notebook's `open_cells()` whose reading
  failed was recorded in the history as done ("Opened 50% of the cells") and
  answered as such; it now fails, said once, and records no step. The agent
  or notebook gets the sentence the person read, with its reason (R9M-12),
  not a generic one.
- **A subsample the dataset cannot be read with is a refusal, not a
  connection problem** (R9M-12): a stratified field or spatial embedding the
  dataset lacks (or one that is not numeric, or a column past its width) said
  "Try again; if it keeps failing, check the address and your connection".
  Each is now a refusal written for the person ("A spatial subsample needs
  the embedding obsm['X_sparse'], which this dataset does not have as a dense
  array."), settled (asking again meets it again), by the same classification
  of read failures every notice uses; whether a failure may pass now lives in
  that one table.
- **Spatial's Embedding menu offers what Data mapping shows** (R9M-13): it
  listed every `obsm` array, vector fields (`velocity_umap`) and arrays the
  mapping leaves off included. It now lists the embeddings Data mapping shows
  in 2D or 3D, its 2D default first, and follows the mapping when it changes;
  an agent's `cells_open` naming another is refused with the ones it shows,
  and one naming none thins the 2D default. Opened before any mapping (a link
  naming `sample=spatial`), the reader thins the best-ranked array of two
  columns or more, never a vector field.

#### Sessions and sharing

- **Opening a section is no work to lose.** In a browser that keeps nothing
  for the site, opening the Agents section (a step) made the next dataset
  ask *Your work is not kept*. Steps that only arrange the page (the layout)
  are no longer counted as unkept work.
- **A session file dropped on the page loads.** Dropping a `.cellucid-session`
  file (or choosing it in **Open… › File…**) was refused with "Open its
  dataset, then load it with State › Session › Load…"; it now loads as
  **Load… › Open a session file…** does, its dataset opened first when
  another is shown. The sentence Python and R say for a session file given to
  `cellucid serve` names the web app's way in (`input_formats.json`).
- **Every saved state of an earlier version reads, whatever key was added
  since.** A new key in one record of the state (the Camera Path's captions)
  refused every older session and every sample's starting view. Every state
  domain now reads its records through one declaration of their keys
  (`app/state/domains/records.js`): a key added since is filled with the value
  its absence means, a key dropped since (the analyses' remembered sizes,
  moved to the layout) is left out, and only a truly wrong record is refused,
  naming where. Each sample's starting view and saved states written by
  0.9.7 and during this release are read by every domain in the tests.
- **A saved state written before a part of it gained a key still opens.**
  Adding the camera path's **Captions** made every sample's starting view and
  every older session fail to restore the camera path ("Camera Path state
  must contain exactly …"). A key added to a part of the state after files
  were written is now read as the value its absence means and written only
  when it differs, so a record of before reads, and is saved again, byte for
  byte; only a truly wrong record is refused (as F181 did for the settings).
- **A setting saved by another version never fails the others.** Saved
  settings are read one at a time: a value this version cannot use (a control
  removed since, a choice no longer offered, a Figure Export value no version
  has) keeps the setting as shown and one notice names it, with its reason
  under Details; everything else is restored. Before, one such value failed
  every setting: the Figure Export record's earlier version (an SVG strategy
  left unchosen) lost an hour-old session's settings and a dataset's starting
  view. Figure Export records of the earlier version are read in today's terms
  (record version 2: unchosen is Full vector).
- **A saved state is restored part by part** (F155): **Continue** on a history
  kept minutes before the Expression values were renamed (`normalise` →
  `normalize`) failed with "a problem inside Cellucid stopped it" and dropped
  everything, the views, the steps and the conversation with the one part it
  could not read. Each part of the state (views, fields, filters, highlights,
  the Expression values, …), each tracker's steps and each conversation is
  now restored on its own: a part that cannot be put back is left as it was
  shown, the rest is restored, and **Session restored but for N parts** names
  each in a sentence of its own, its words under **Details**. A part saved by
  an earlier or a newer version of Cellucid that this one cannot read says
  so ("Expression values could not be restored: it was saved by an earlier
  version of Cellucid (0.9.7), and this version cannot read this part"),
  never "a problem inside Cellucid", and so does a load that fails whole on
  such a file. The history records the parts left as shown as a step of its
  own ("Not restored from the saved state: … (kept as shown)"); a tracker
  this version cannot read is left out and named, and when it was the one
  shown, what it showed starts a new tracker of its name with its
  conversations. The saved state continued is then not written into, and the
  work goes on as a new saved state beside it. A kept conversation the chat
  cannot show is said instead of dropped, and an agent's
  `history_saved_open` returns the same sentences (`not_restored`).
- **The welcome comes first.** On a start that shows the welcome, "Continue
  where you left off?", a sample view's notices and an opening's question
  (which object of an R file to open) were shown on top of it at once; they
  now come once it is closed, in that order (R9U-14).
- **A session file or saved history from another version is restored part by
  part too** (F164): a file holding a part of the state this version does not
  have, or one tracker whose records it cannot read, was refused whole ("This
  history holds state this version of Cellucid does not have … Open it in the
  version of Cellucid it was saved with"). That part is now left out of every
  step and named ("A part of the state this version of Cellucid does not have
  ("…") is left out: it was saved by a newer version of Cellucid (…)"), such a
  tracker is left out and named, and everything else is restored: Load…,
  **Continue**, a saved state, an agent's `history_saved_open` and a dataset's
  starting view alike.
- **A load says when it sums the cells' totals** (F167): restoring normalized
  values (or a mapping, or a replayed step, that needs them) reads every gene
  of the matrix once for each cell's total, with nothing on the load's
  progress line but the file it had read. The line now says "Summing each
  cell's counts for normalized values…" while it does, with **Cancel**, on
  every load: a session file, a saved state, a file of steps replayed, an
  agent's load (asked first, as before).
- Loading a session file names its dataset as the Dataset box does ("the
  session's dataset "Developing human immune system (Suo)" from the
  samples"), not by its id. It ends in one notice ("Loaded "Tracker 1
  (shared)" beside yours — Session fully restored.", with **Undo**), not two
  saying the same, and the dataset it opens no longer offers "Continue where
  you left off?" in the middle of the load.
- A dataset opened whole is shared whole: **Share** (link, badge, embed),
  Subsampling's **Copy link** and the page address no longer name the first
  cells a progressive opening shows while it grows
  (`cells=131072&sample=uniform`), which reopened a quarter of Suo. A session
  or autosave saved while such an opening grew names the cells it grows to,
  so it no longer reopens only its first cells.

#### Data mapping

- **Gene names and highly variable genes are checked on their own matrix's
  var axis** (R9M-1, R9M-2): a session saved on `X` or a layer with Highly
  variable genes `highly_variable` was refused on the same file ("… which this
  dataset does not have"), because the file opens on `raw.X` (it holds more
  genes) and the column was looked up in `raw.var`; and a change naming only
  the matrix (an agent's `mapping_set`, a notebook's
  `set_data_mapping(expression="raw")`) applied a mapping that kept a column
  `raw.var` does not have. Each value is now checked in the context of the
  mapping it belongs to: its matrix's var axis (`var` for `X` and every
  layer, `raw.var` for `raw.X`) and its modality.
- **The genes groups follow the matrix by one rule, in the panel and in
  `mapping.set`** (R9M-2, R9M-3): a change of Expression keeps the gene-name
  and highly-variable column when the new matrix's var axis has one of that
  name, else takes that axis's default (its var index; `highly_variable`
  when it is true/false there, as opening on that matrix does). Moving from
  `raw.X` to `X` now starts Highly variable genes on `highly_variable`
  instead of leaving it None. A change of Modality takes the modality's
  defaults. A change naming those groups itself keeps what it names.
- **A MuData opens with its modality's neighbor graph** (R9M-7): a MuData
  whose graph lives in a modality (`mdata.mod['rna'].obsp['connectivities']`,
  where scanpy writes it and where `prepare()` writes each modality's graph)
  opened with Neighbor graph None, so `prepare()`'s own MuData drew no edges
  until the graph was chosen by hand. Without a global
  `obsp['connectivities']` the dataset now opens with the graph of the
  modality it opens on (`rna:connectivities`). The rule (`defaultGraph`,
  `default_graph()`, R's `.default_graph`) is pinned by the shared conformance
  fixture, version 8, in all three repositories.
- **Every width gets its default by one ranked order** (R9M-11): a file's
  UMAP named only its own width, while a file without one got a guessed
  default for every width, so `X_umap` beside `X_umap_3d` opened with no 2D
  default and `X_umap` beside a 3-column `X_phate` with no 3D one. Now every
  slot the UMAP (or a launch) leaves empty takes the best-ranked array of
  exactly its width, marked guessed; the best-guess notice is kept for a file
  that names no default at all. The rule is pinned in the shared conformance
  fixture (version 9) and now applied by `cellucid.data_mapping.resolve_embedding_slots()`
  and the R package too, each checked against every `embeddingSlots` case.
- **Default marks what keys 1–3 show first** (R9M-5): a dimension the mapping
  names no default for (a file with `X_umap` and a 3-column `X_phate` named
  only the 2D one, before R9M-11; a session or an agent's mapping can still)
  showed no **Default** although key 3 showed its first embedding, and any
  later edit of the checklist wrote a default for it that nobody chose. The checklist now marks that embedding, and Apply names a
  default only when one is pressed or a dimension gets its first embedding.
- **A MuData modality's embeddings have readable names** (R9M-10): they were
  listed by their keys (`rna:X_phate3`, `prot:X_umap`) beside `UMAP` and
  `t-SNE`; they are named as the dataset's own, their modality after them
  (`UMAP (prot)`, `PCA 3,4 (rna)`), and a Seurat key with dots is split at
  them (`X_wnn.umap` is `WNN UMAP`).
- **Two embeddings are never named alike** (R9U-6): Garcia, He and Kanemaru
  listed `X_umap` and `X_umap_2d` both as "UMAP · 2D". Arrays a menu would
  list the same way (the same short name, namespace and width) are named by
  their keys, "UMAP (X_umap)" and "UMAP (X_umap_2d)", by the one naming
  function every menu, the Data mapping checklist, the Dataset box, the
  agents' reads and a figure's title use (`representationNameOf`); arrays the
  dimension tells apart (`X_umap_1d`, `X_umap_3d`) keep the short name.
- **A data-mapping change stops on Stop** (F147): an agent's `mapping_set`
  (or an undo or a jump putting a mapping back) ran on until every embedding,
  gene list and cells' total it read was done, minutes on an atlas, holding
  the history. The action's signal now reaches every wait and read of the
  change: stopped, it ends at once, publishes nothing (or puts back what it
  had), and the next change runs.
- **Every way into the sidebar opens it before its sections** (R9K-2's
  siblings): Data mapping's **Show** (the best-guess notice) and the "Genes
  read for" link opened their Session section in the moment the sidebar's
  own step started, which left a "Changed outside the history" step when the
  sidebar was hidden. They and `?` now take one path (`revealInSidebar`): the
  sidebar first, then the sections and the focus.
- The Expression values line for a matrix stored as `log1p` names its place
  as every other message does: "raw.X holds log1p values already (Data mapping
  › Genes › Values)." (it said "raw holds … (Data mapping, Genes: Values)").
- **`raw.X` is named `raw.X` everywhere** (R9M-14): messages, analyses (the
  result line, a CSV's `expression_matrix`, Details), the Expression values
  line, the Scanpy reason under Genes and the session messages said `raw`
  where the Expression menu says `raw.X`. One name, the menu's, is now given
  by one function (`expressionMatrixName`).

#### Legend, filters, highlights and the view

- **A double-click names a view in Safari too.** A click on a view's row in
  **Compare views** selects it, and the selection rebuilt every row: the
  second press of a double-click on a row not yet selected could land on a
  row already replaced, and WebKit then sent no double-click, so the name
  field never opened. A selection now sets the rows' state in place; a row
  is rebuilt only when what it says changes.
- **Point size ends at what can be drawn.** The slider ran to 200 while the
  largest cell drawn is 160 pixels, so its last 3.4 % drew the same sprite
  while the readout climbed. Position 100 is now 160 (the curve keeps 0.25
  at position 0, so a stored position means a size at most 3.4 % smaller);
  the reachable sizes run from 0.053 to 160.
- **A kept view no longer crops a wide map.** Every pane of a grid took the
  field of view of the whole area beside the sidebar, which is landscape, so
  two views side by side (upright panes) cut Pancreas's 2-D UMAP off at both
  sides. Each pane now takes its own shape's (45° across an upright pane, as
  on a phone), and is drawn, picked, zoomed at the cursor, framed for an
  agent and exported with it.
- **Grid compare keeps its panes a usable size.** Two or three views always
  stood side by side, so on a 390 px phone held upright three panes were about
  120 px wide each. The grid now takes the arrangement whose panes are largest
  along their shorter side, in the area the sidebar leaves: side by side on a
  wide window, two rows of two for three views beside an open sidebar on a
  1440 × 900 screen, one above the other on an upright phone. A figure of
  every view (**Export all views**) keeps the arrangement on screen.
- **A phone turned sideways with the sidebar closed uses the whole width.**
  The closed sidebar slides to its new place when the window widens past the
  phone layout; the area beside it was measured during that slide (counted as
  shown) and never again, so the panes stood in 574 of 844 px with an empty
  band at the left. Any measure taken while the sidebar slides is taken again
  once it stops.
- **A category label you give is one line.** A rename by an agent or a
  notebook made "Alpha\nnewline" (two lines in the legend, a category the
  community annotation cannot key), and Create categorical's actions took a
  page label with spaces at its ends that its panel refused. Every label given
  now (a rename, a merged category's label, a group, a page's, an overlap's
  or the uncovered cells' label in Create categorical) follows one rule, the
  field names': no leading or trailing whitespace and no control or
  zero-width character ("… cannot be a category label: it holds the control
  character U+000A …"). Labels a dataset brings are read as they are.
- **Smoke while views are kept is dimmed, and says why.** Choosing Volumetric
  smoke with a view kept was taken back afterwards with "Volumetric smoke
  requires a single view. Clear snapshots first." (the app calls them kept
  views), and an agent was refused "cannot be changed now: it keeps its value
  until what it needs is there". The choice is now off while views are kept,
  its tooltip and an agent's refusal say "volumetric smoke draws one view;
  remove the kept views first (Clear, beside Layout)", as Keep view is off in
  smoke. A select's choice turned off says its reason like any control.
- **A highlight reads against any coloring.** Garcia colored by `LVL1`, its
  10,993 germ cells highlighted, looked exactly as before: a highlighted cell
  was a cyan dot in a cyan ring, and the germ cells are sky blue. A
  highlighted cell is now a cyan dot inside a black ring (the ring is black or
  white, whichever contrasts with the dot, `highlightRingColorFor`), so any
  cell or background color is far from one of the two. One rule for every
  shader quality, 2-D and 3-D, every view of a grid, and every figure format:
  the vector SVGs, which drew the cell's own color in a cyan ring wider than
  the screen's, now draw the dot and ring the screen draws, at its radii.
- **Colors the app picks are told apart.** A grouped field's group took the
  color of its first category even when another group's looked the same
  (Kanemaru's atrial cardiomyocytes and the rest: amber and bright orange); a
  new `Unassigned` took the palette color of its place, which a kept category
  could already have; a new highlight page took the default of the page count,
  so after a deletion it repeated a remaining page's color (and the ninth,
  rose, read as the second, red). One rule now gives every such color: kept
  when it is at least 0.08 from the others in OKLab, else the next palette
  color that is (`distinct-colors.js`). Colors a person or a file chose are
  never changed, and a file's categories keep the color of their place.
- **The cell card shows a long label whole.** A category label, a field name
  or another field's row longer than the card's column was cut with "…", and
  the card (which follows the pointer and takes no clicks) offered no way to
  read the rest; each now wraps inside the column, and a card with many rows
  keeps inside the window, naming the rows it leaves out.
- **A clustering of numbers colors the cells.** Coloring by a field whose
  categories are numbers or booleans (a Loom `ClusterID`, as velocyto and
  loompy write it) failed in every dataset: "Could not load the field:
  [CommunityAnnotationSession] category at index 0 must be nonblank …". The
  legend handed community annotation the numbers themselves; it now keys
  categories by their text, as the legend shows them ("0", "1", …), and a
  field whose categories a repository cannot key (a blank one, one with
  spaces at its ends) is shown and colored as any other, just not open to
  voting (an agent's read of it says why).
- **Views can be named.** Double-click a view's row in Compare views (or
  press `F2` on it) and type a name: `Early`, `Day 7, treated`. The name
  shows beside the view's number in its row, its pane's title chip, the
  Navigation panel's Applies to row and the embedding picker, and in every
  history step on that view (`Moved the camera of View 2 "Early"`). Naming is
  one step (`view.rename`, undoable), kept in sessions (the `views` value's
  view entry carries `name` only when the view has one, so sessions and
  starting views without names are unchanged); `view.keep` takes a `name` to
  keep and name a view in one step. The number stays the view's address in
  actions, so a name never has to be unique; a name is one line and cannot
  be of the form `View 3`. An agent reads each view's name in `views.state`
  (`label`) and `view.picture` (`label`). An agent building a multi-panel
  figure of a time course named no panel and put the ages in the caption
  instead (bug hunt F130).
- **A redraw of the legend keeps the list where you scrolled it.** Every
  redraw of the field shown (a rename, a vote, a count arriving) built the
  list anew from its top: in a long legend a rename far down threw you back
  to the first row. The list keeps its place while the field stays the same.
- **The outlier filter's slider keeps the value you gave while it is
  applied.** A redraw of the control between your change and the filter
  taking effect (its action waits behind a coloring still under way) put the
  thumb back at the old threshold for a moment, and a second commit then
  filtered by that one. The thumb stays where you put it until the change is
  in effect, and a value given but not committed yet is left alone too.
- **What stands above a dialog goes back to its own place.** The
  notification center and the agents' dock are lifted into an open dialog
  and put back when it closes; side by side, the first to go back found the
  other away and landed at the end of the page, so the two swapped places
  (and which is drawn over the other). Each goes back where it stood,
  whichever goes first.
- **The sidebar legend lists a field's categories in the field's own order**
  (the file's, Scanpy's), as the Group dialog, the figure export's legend and
  the analysis plots do: Pancreas `clusters` reads Ductal → Ngn3 low EP → …,
  its trajectory, where the legend sorted alphabetically and hid it; a
  renamed category keeps its place, groups stay in the order they were typed.
- **A page's category takes its page's color in Create categorical**: one
  color for one page in the builder's page row, its preview, the legend and
  every view (the preview showed the palette's color beside the page's own).
  The field keeps these colors as its own (its `palette`), in a session too
  (`field-overlays` user fields gain `palette`: one `#rrggbb` a category, or
  null), and a category kept through a merge or a delete keeps its color.
- **"Unassigned" is the one spelling of cells without a label**: deleting a
  category moves its cells into `Unassigned` (it was `unassigned`; the
  action `category.delete` is version 2, and a version 1 step that named no
  label still moves them into `unassigned`), as Create categorical names the
  cells on no page.
- No floating window covers the sidebar's toggle: windows are placed, moved,
  resized and maximized below it (on a phone, a tall floated panel covered the
  open panel's close button).
- **Dock all**: while any sidebar section floats, a button at the top of the
  sidebar ("Dock all 3 panels"; "Dock Navigation" for one) docks every
  floating section back into its place, as one step
  (`layout.panel.dock-all`; `layout_panel_dock_all()` in Python,
  `cellucid_layout_panel_dock_all()` in R). With every section undocked the
  sidebar was a column of logo and footer, and the way back was each window's
  own button.
- Floating windows fit the visible screen: their height cap and their place
  are taken from one measure, the visual viewport (a phone's on-screen
  keyboard no longer hides a window behind it), and an analysis's expanded
  view sized on a large window is cut to a smaller one and moved inside it
  (it ran past the edges, its close button out of reach); widening gives the
  size back.
- A category name cut short in its legend row shows whole on a click or a tap
  (and leads its hover text).
- The **Filters** list names the smaller side of a legend ("cell_type:
  showing only CD4 T, CD8 T" after **Hide all** and two ticked, not the 67
  hidden), says "hiding every category" when none is shown, and counts
  Missing only when some open cell has no category, as the legend lists it.
  Figure captions and view labels read the same.
- A highlight group is named by what it holds, its count once beside it: an
  Annotation selection by its categories, most cells first ("NK + ILCP",
  "A, B, C +2 more"), or by its field and range ("n_genes 120 – 340"); a
  Lasso, Proximity or KNN group by its tool ("Lasso", "Lasso (2)"); one made
  without a label from Python, R or an agent "Highlight". It read
  "Annotation (2 steps) (36,214 cells) (36,214)", cut short. In the
  Highlights and Filters lists the count stands in a column of its own,
  never pushed out by a long name, and **Create categorical** names its
  categories the same way.
- **Group categories**, a new button of the Categorical obs row (after
  **Duplicate**), makes a coarser annotation in one step: type the group each
  category joins (on Pancreas, *Ductal* and *Ngn3 low EP* as *Progenitor*,
  *Alpha*, *Beta*, *Delta* and *Epsilon* as *Endocrine*), and a new field of
  those groups is made and shown, the field it comes from untouched. A
  category left empty keeps its name, or joins **Others into** (*Other*). One
  **Undo** removes it (it took nine steps and nine undos: duplicate, merge one
  pair at a time, rename). It is the action `field.group`. The Categorical
  obs menu names the grouped field at once (it kept naming the field grouped
  from while the legend showed the groups).
- **Create categorical's Each intersection takes any number of pages.** It
  refused a 33rd page ("Each intersection supports at most 32 pages"), as did
  `field.create-from-pages` (an agent's 47 pages): a cell's pages were one bit
  each of a 32-bit number. Each cell's set of pages is now interned however
  many pages there are (`utils/page-membership.js`), and an intersection is
  named by its pages' positions (`0+2`) in the session
  (`userFields[].intersectionLabels`, no other key read: the field's step
  names its pages by name and replays as before). The intersection categories
  come in order of how many pages they combine, then by their pages (a field
  of three or more pages may list them in another order than 0.9.7's), and
  one combination labeled twice in `intersection_labels` is refused instead
  of keeping the last label.
- **Create categorical** shows the field it makes at once, as **Group** and
  **Duplicate** do (`field.create-from-pages` now changes the coloring, and
  one Undo restores it). Its preview is the field made: the uncovered cells'
  category in the color it takes (it was drawn gray, then colored from the
  palette), and no `Unassigned` category of 0 cells when the pages hold every
  cell (the preview said one category, two were made). Its column name follows
  the one rule of every field name: `:` is allowed (it was refused, though a
  rename took it) and a zero-width character is refused in the builder (it
  was accepted, then refused by a notification after **Create**). What it
  still waits for ("Name the new column.") is a gray hint, no longer a red
  "Column name must be a valid non-empty trimmed field name" before anything
  was typed; its other refusals are plain sentences. A page row shows the
  page it comes from only when its label was changed, and its label box has
  the room.
- Highlight page tabs are as wide as their names, up to the row: three pages
  read "Pa…", "Pa…", "P…" beside empty space.
- The legend's alphabetical order is the one a person reads: case aside and
  numbers as numbers (`merged Beta + Alpha` sorted after every capitalised
  label, `Cluster 10` before `Cluster 2`).
- The merge and delete-category confirmations say **Undo** takes the change
  back; they told the person to restore "the original column from Deleted
  fields", which an edit of a field of one's own never puts there.
- The gene box names the gene shown after Escape or a click elsewhere, never
  a name typed that chose no gene.
- **A gene colors a kept view.** After **Keep view**, choosing a gene on the
  kept view (the gene search, or `color.by` from an agent, Python or R) did
  nothing and the action answered null: every view holds its own copies of
  the fields, and the gene search chose the gene as the field object of the
  view it was listed in, which the kept view's list does not hold ("no longer
  on the gene list"). The gene is now read from the list of the view in
  focus, so every kind of field (a column, a gene, a copy of either) colors
  any view. A coloring that ends without being shown (something else holds
  the fields, a newer choice, Cancel) is refused with the sentence the person
  is told, never answered with null as if done.
- **A view's × removes that view.** With exactly one kept view, the × on the
  live view's row removed the kept view instead (and with it the kept view's
  name), and `view.remove` of View 1 did the same. Removing the live view now
  hides it, however many kept views there are; the kept view left shows alone,
  as the last view it has no ×, and **Clear** brings the live view back.
- A viewer background you never picked follows the theme (Studio (light) ↔
  Studio (dark), Atelier ↔ Nebula, …); one you pick stays. A sample's
  starting view no longer sets the background (it had stored Studio (light)
  as your pick on the first visit).
- A message naming a long link breaks inside the sidebar. The **Stored** line
  names the datasets kept and what they take; the other kinds are in its
  title and in **Manage…**.
- **The palette menu marks the palette shown.** After an Undo (or an agent)
  changed a continuous field's palette, the color bar's menu still marked the
  one chosen before, and choosing that one did nothing. The menu is drawn as
  it opens, and choosing a palette changes only the palette: it no longer
  rebuilt the legend, which turned **Live filtering** back on and dropped
  bounds moved but not yet filtered (R9F-1).
- **Live filtering stays as you set it.** Turned off, it stays off for every
  continuous field, gene, view and palette of the page; it turned itself back
  on at every field or view change, which the docs warned about (R9F-2).
- **Renaming in place keeps what you typed.** A field's, a category's and a
  saved state's name are taken when you leave the box, as Enter does (as a
  view's and a highlight page's already were); leaving it dropped the name.
  A name refused ("A field with this name already exists") is said on the box
  until you change it; it was only a hover title on a red border (R9F-6).
- **Two deleted fields of one name can each be restored.** Delete a field,
  rename another to its name and delete that too: both rows of **Deleted
  fields** refused Restore and Delete for good ("Two deleted fields are
  named …; restore or purge one from Deleted fields"). The rows read
  `score (1 of 2)` and `score (2 of 2)` and act on their own field; the
  actions `field.restore` and `field.purge` take `which` for it (R9F-4).
- **Delete for good says what it does.** Its confirmation said the field
  "can no longer be restored, in this session or from a saved state", but
  Undo brings it back; it now says so, the history names the step `Deleted
  <field> for good` (it read `Purged …`), and the row's button is a bin, so
  a deleted field's name is no longer cut to `S_s…`.
- History steps name a palette as the legend does (`Drew Ins1 in Plasma`,
  not `plasma`) and write a range's bounds as the Filters list does
  (`-0.448703`, not `-0.44870302081108093`) (R9F-5, R9F-11). A view name of
  the form "View 3" is refused with "…, or none to number it again" (it said
  "no name (null)").

#### The history

- **Work that is not kept is never dropped without a question.** With
  Autosave in this browser off (or a browser that keeps nothing, or no room
  left), opening another dataset threw away every step on the one shown
  without a word, and so did a reload. Now, while the history shown holds
  steps made since it was last saved and they are not kept, every way of
  replacing the dataset (the sample list, My datasets, a file, a link, a
  server, a saved state or session file of another dataset, None, the
  benchmark) asks first: **Save…** (Save file), **Open anyway** or
  **Cancel**; an agent's `dataset.open` / `history.saved.open` asks through
  its permission card, under a new kind **Replacing work that is not kept**
  (`replaces-unkept`, asked by default); a reload or a closed tab gets the
  browser's leave-page question — in that case only. Canceling a dataset
  change yourself (the opening's Cancel, or keeping the dataset here) no
  longer says "A newer dataset selection is now active.", and the Dataset box
  no longer keeps saying "Loading …" after a change begun outside the list
  (an agent's opening you denied) was stopped.
- **Every slider names its values in its own words.** The velocity sliders,
  the navigation speeds, the forced level of detail, the outlier filter, the
  camera path's Default speed and every slider of a panel built from the
  shared slider now define their readout (`ui/core/control-readout.js`), so
  no comparison of steps names a slider by its position; a page test checks
  every slider with a readout.
- **Comparing two steps names each setting as the sidebar reads it.**
  `history.diff` (and every comparison built on it) said "Grid density: 37",
  "Point size: 16.5" and "Render mode: smoke": the slider's position and the
  option's value, while the sidebar and the step's own label read "48³",
  "0.75" and "Volumetric smoke cloud (beta)". A slider's owner now defines
  the readout of any value (`ui/core/control-readout.js`) and writes its own
  readout with it, so the readout and the history say one thing: Point size,
  the depth cues, every smoke slider and the connectivity sliders. A choice
  is named by its option; a slider without a readout of its own by its value
  when its readout is the value, else as a position ("slider at 37 of
  0–100").
- **Stale results name each input with its own cause, and follow the view
  they were made in** (F156): a test went "stale since step 4 ("Showed UMAP ·
  2D in View 1"): Filters changed … Hidden cells changed … Views changed"
  when step 4 only kept a view. Every reason cited the first change of any
  input, and keeping a view changes each per-view value. A step's inputs held
  per view (filters, hidden cells, the field colored by, the views) are now
  compared at the view it acted in, so keeping another view, focusing one or
  changing another view's filters marks nothing stale, and each reason names
  the step that changed that input ("Filters changed at step 7 ("Hid Alpha of
  cell_type in View 1")."), in `history_stale`, `history_step` and the step's
  details.
- **With several views, a step says which view it changed** (F157): "Colored
  by gene Sox9" went to the view in focus and named none, in the history and
  in the methods text. A step that changed some views, not all, now names
  them as the Compare views list does, with a view's name ("Colored by gene
  Sox9 in View 2 "Early"").
- **Save AnnData with this state…** keeps a field's colors: a categorical
  field you or an agent made is written with the colors it is shown in, as
  Scanpy keeps them (`uns['<key>_colors']`) -- in the view in focus, the
  colors given to its categories, the rest as it opens (a field made from
  highlight pages in its pages' colors, any other in Cellucid's own) -- so
  the copy opens, and plots in Scanpy, as the view showed it; the Colors row
  of its plan says the fields' colors went with their columns. Every column
  the copy writes, renames or relabels takes
  its colors along by the packages' rule: none for a column written without
  colors of its own (stale ones removed), carried by a rename, a merged
  category keeping its first category's color. Before, no colors were
  written, and a rename or a merge left colors that described another
  column.
- **Save AnnData with this state…** also writes the colors of a column of the
  dataset you or an agent recolored in the view in focus (`uns['<column>_colors']`,
  under its name and labels after the save's renames: the categories recolored
  in their new colors, the rest as the column opens), one `column-colors` row
  of the plan each (`write`, `unchanged`, or why not), as Python and R apply
  it; and a live viewer's categorical results (viewer protocol `result`,
  `categorical` meta) carry `colors`, one lowercase `#rrggbb` a category in
  its order, the colors it is shown in, in the view in focus, sent again on a
  recolor or a change of the view in focus, so a pull writes them (protocol
  fixture `c117822a…`). Before, recoloring a dataset column reached no file,
  and a pull wrote no colors.
- Reopening a saved analysis whose live view was removed (three kept views
  only) no longer throws "A hidden live view cannot be the active view": a
  restore applies each view's coloring and filters with that view focused,
  the removed live view too, and awaits field reads between them; the view
  list and the sidebar drawn in a frame between them read the restore's
  cursor as the person's focus. They now wait for the visit to end.
- The methods paragraph, the history report, a figure's provenance and the
  agents' dataset summary name the dataset's citation (its `uns['cellucid']`
  `citation`): they looked for it in the dataset's source, where it never is,
  so none of them cited any dataset, and an agent writing methods said the
  sample "has no citation recorded". A figure's provenance JSON now carries it
  as `dataset.citation`, beside `dataset.source`.
- **Continue** opens the history kept in this browser after a group was made,
  extended, removed or restored: the grouped steps were not written again, so
  Continue refused the whole history ("a group cuts the run of turn …"). An
  edit restored from Recently deleted is no longer listed there again after a
  reload.
- **Verify** of a branch that makes a field (a grouping, a new annotation) and
  colors or filters it says it reproduced: the field made again under a new
  id is the same field (it said "not reproduced: Colors").
- An agent's turn goes on through its own **Verify**, a tracker it starts,
  uses or replays: the chat stopped the turn and said its call was lost ("The
  page closed before this call finished"). An agent (or a notebook) starting
  another tracker no longer fails. A squash or a branch delete the history
  would refuse is refused before the person is asked, never after they
  allowed it.
- A named save or a session file made during an agent's turn (its own
  `history_save` or `history_share`, "with our conversations") keeps that turn
  as it was when saved, said so: opened, it said "The page closed before this
  call finished; nothing was applied" of the very save.
- The card for an agent's delete, squash or tracker removal no longer says
  "the history cannot bring back" what goes to Recently deleted.
- A step whose name wraps differently at a new width (the sidebar resized,
  the tree's scrollbar shown or gone) is measured again: such a row overlapped
  the step under it, and the last step was cut off by the tree's own height.
- Closing a progress notice while its work runs no longer fails that work (a
  field read failed, and a dataset could open uncolored).
- A theme change that moves a background you never picked is one step (it left a
  "Changed outside the history" step after it), and undoing it no longer keeps
  that background as your pick: only a background you pick yourself is kept.
- **Verify** says why a branch got no verdict in a sentence: a problem inside
  Cellucid keeps its words under Details, and a first state that cannot be
  shown no longer reports the raw read error.
- A drag (a slider, a camera move) brought back to where it began leaves no
  step: it left one that changed nothing, and Undo then seemed to do nothing.
- Steps name what they changed as the app names it: an analysis option as
  its form does ("Dot plot: Color set to Mean, scaled 0–1 per gene", "Showed
  Values, Layout options of Dot plot"), an option of the expanded view as its
  options column does ("Composition: Plot set to Per sample"), and what a
  view shows as the picker does ("Showed UMAP · 2D in View 2"). Several
  options set at once are named by the values they changed to, the first
  five ("Marker genes: Group by clusters; top genes per group 5"), not by
  their labels alone ("Marker genes: Group by, Compare against, Groups to
  test, Direction, Top genes per group set"), and only those that changed (an
  agent sending five options of which one differs: "Opened Marker genes: Top
  genes per group set to 5"); **Score as signature** and the other sends
  likewise ("Opened Gene signature: Gene sets set to T cell (2 genes), …").
- A camera path keyframe is never marked stale when the camera moves on. An
  action may now declare what it captures (`captures`, in the shared action
  catalog).
- **Load…** names an autosave by its latest step ("Colored by gene CD33",
  with its **Auto** badge), not "Autosave" on every row. An autosave written
  before this release reads "Autosave" until it is written again. Opening one
  names it the same way ("Opening the autosave at “Colored by gene CD33”…").
- **Stop reaches the history's moves.** An agent's undo, redo, jump, branch
  switch or history edit that waits for work under way (an analysis running)
  never runs once you press Stop; it used to run when that work ended. A move
  stopped while it puts a state back (re-reading cells) puts back what was
  shown.
- **Stop reaches a replay's read and a share's wait.** An agent's
  `history_replay` stopped while it read the saved state, and a
  `history_share` or `history_anndata_save` stopped while it waited for the
  actions running or for the session operation before it, ended only when
  that wait did; each now ends at once (a read and a wait change nothing),
  and the next move runs. The agent's replay now reads the saved state as
  Replay… does, every record checked and with the step parameters only its
  later part holds (it read the first part alone).
- **The guided tour ends where you act.** A change, an undo or a jump made
  while the tour plays (by you or an agent) ends the tour there; its next step
  used to move away from what you just did, and its end jumped back over it. A
  note or bookmark added meanwhile lets it play on.
- **A field you made has no original name** (R9F-8): renamed, a copy or a
  group of yours showed `S score *` with **Put name back**, and lost both
  once the session was saved and opened again. It now shows neither, in the
  page and after a round trip; a dataset's own field keeps both.
- **An agent's question is answered over an open dialog.** An agent that
  asked something (opening a link not seen in this tab, a large download)
  while a confirmation or settings dialog was open put its question under
  the dialog, where it could not be answered, and waited for good. While a
  question waits, the agent's dock now stays over any dialog.
- **One way of saying how long ago.** Comments of the annotation voting
  ("5m ago") and Stored ("used today") now say it as every other list does
  ("5 minutes ago", "used 1 minute ago"); a step's details say how long it
  took as every other duration does ("took 1.3 s", "about 1 min 30 s").
- **A folded run says what it did.** Consecutive steps that only present
  the analysis folded into a row named by their kind ("Colors ×2" for two
  legend switches); the row now reads as its steps do ("Hid T cells of
  cell_type · Showed T cells of cell_type", "Moved the camera of View 1
  ×6"), naming the first two changes and counting the rest.
- **What you do while an action runs stands when it is stopped** (R9W-5):
  closing an analysis window with **×** while its run ran stopped the run,
  whose rollback then put the window back, and no "Closed" step was kept.
  A stopped action now puts back only what it changed itself: a change of
  your own made meanwhile (a window closed or moved, a panel docked, a
  setting, the sidebar's width) stays, as its own step.
- **An agent's replay is priced before it starts** (F166): `history_replay`
  had no estimate, so a replay of a step that summed the cells' totals for
  normalized values ran that pass without asking, though applying a branch or
  verifying one asked. It is now priced by the time its steps took when
  recorded, so it is heavy work when they were, and heavy work a replayed step
  (or a step a delete applies again) finds is its replay's to ask about.
- A notebook command naming a language other than Python or R is refused
  before it runs (the history would have refused its step after it changed
  the view).

#### Analyses

- **Expression values: normalized and log1p, one setting for the whole app.**
  **Session › Data mapping › Expression values** (collapsed, after Genes) has
  two switches, both off by default: **Normalize to library size** (each
  cell's values scaled to 10,000 counts, `x · 10,000 / cell total`, scanpy's
  `normalize_total(target_sum=1e4)`) and **Log transform (log1p)** (`ln(1 + x)`;
  both: Seurat's `LogNormalize`). Coloring by a gene, the legend (its units:
  `Values: log1p(norm.)`; a figure's legend titles `Ins1 · log1p(norm.)`), a
  gene's color range and range filter, highlights by value, the cell card,
  an agent's `cells.values`, `field.summary` and `dataset.summary` (each names
  the values), and Marker genes, Differential expression with cells as
  observations, Dot plot, Gene signature, Detailed and Correlation all
  use these values (Quick reads obs fields only); obs fields never change,
  and pseudobulk always sums raw counts. A cell's total is the obs column the file names
  (`uns['cellucid']['totalCounts']`, `total_counts` as `prepare()` writes it),
  else one named like a total (`n_counts`, `nCount_RNA`, `nUMI`), else every
  gene summed once — the section's line says which and what that reads at the
  link, and an agent is asked first (heavy work); the totals are kept and
  shared, and cells added later are summed alone. Offered for a mapped matrix
  of counts whose values are as stored (otherwise disabled, saying why). A
  view colored by a gene keeps it, read again; its color range and filter
  are dropped (they were in the old units). It is one history step, saved in
  sessions (a session or saved state without it means as stored, also when
  restored onto an open dataset), the action
  `expression.transform.set {normalize, log1p}`, and every result says the
  values it computed from (`about.values`) and goes **Stale** when the values
  it read change. Before, means, fold changes, tests and scores were of the values as
  stored only — on a matrix of raw counts, confounded by each cell's
  sequencing depth — and no view could show normalized values. Held to scanpy 1.11.5's `normalize_total`, `log1p` and both, and
  `rank_genes_groups` (Wilcoxon, t-test) and per-group means on them, to
  10⁻⁶.
- **Summing the cells' totals is asked for wherever it happens.** With
  normalized values on, a change of the data mapping onto another counts
  matrix (or modality), genes read for more cells, the cells opened again,
  and an undo, redo, jump or saved state restoring normalized values each
  sum that matrix's totals when no obs column holds them: the change reads
  them first, priced at the link, before anything shown changes. A person
  sees the price before acting (Data mapping's status line: "Apply sums all
  27,998 genes of layers['counts'] once for each cell's total (normalized
  values): ~74 MB, ~12 s."; Genes read for counts them in its price); an
  agent's action is asked as heavy work (`mapping.set` and `cells.open`
  estimate them; work an action finds only while it runs is asked then,
  through the dispatcher's new `ctx.heavy`). Before, a mapping change, a
  growth or a restore summed the whole matrix unasked, at whatever read came
  first. Each matrix keeps its totals for the open dataset, so mapping back
  onto it reads nothing again; rows read for added cells by every loaded gene
  at once sum the matrix in those cells once (they were summed once per gene);
  a mapping restored with its Expression values shows its matrix in them at
  once, never summing it for values the restore does not end on. The source
  line's time no longer reads "~~12 s".
- **Gene values are named with their units everywhere.** With Expression
  values on, every label of a gene's values carries its units
  (`Ins1 · log1p(norm.)`, `· norm.`, `· log1p`): Detailed's axes, heading,
  summary statistics and CSV column; Correlation's axes, hover and CSV (new
  `x_label`, `y_label` columns); the dot plot's color bar and CSV `mean`
  column; the DE CSV's `mean_a`/`mean_b`; the marker heatmap's color bar; a
  gene's range filter; a range highlight or annotation group drawn on a gene;
  the legend, figure legend and cell card. An analysis names the values its
  run read, so a stale result keeps its own; obs fields and genes made in the
  app stay bare. One helper (`geneValueLabel`) builds every label, and a
  contract test finds any built without it.
- **A result records the values it read, not the setting.** Each analysis
  declares the expression values its genes are read in (the matrix's shown
  transform; raw counts for pseudobulk and UCell, whatever the setting; none
  for Composition), and its provenance (`computed_on.values`, now `null` for
  an analysis that reads no gene values), values line and staleness follow
  it. A pseudobulk result said "normalized · log1p" and went stale when the
  person switched normalization on, though it sums raw counts (an agent re-ran
  it to the same numbers).
- **A saved state restores its Expression values.** Loading a session file or
  a saved state (Load…, an agent's `history.saved.open`) onto the dataset
  open puts back the values it holds, and a state saved without them (written
  before them, or by hand) means the values as stored: they kept whatever was
  shown before, so one file showed other values depending on what was open.
  A state domain says what its absence means (`absent()`); the others keep
  what is shown.
- An analysis's **Expression matrix** takes the matrix as the dataset names
  it: `analysis_differential_set {matrix: "X"}` on a dataset mapped to `X` is
  the mapped choice (it was refused, "Expression matrix has no choice X", and
  an agent retried).
- **Gene signature → Group by.** A signature is scored per category of a
  field in one run (with **Groups to score** and **Cells**), compared across
  the categories as across pages; scoring a program per cell type needed a
  highlight page per type (8 pages in one agent turn on He). Its CSV names
  the rows' category column `group` (`group,set,cell_index,score`), as its
  result table does, not `page`.
- **Pseudobulk with edgeR's TMM and limma's moderation.** Differential
  expression with Replicates sums every gene of the matrix per sample from the
  raw stored values (whatever the expression transform), and fits the samples
  as the limma-trend pipeline does: for a matrix of counts, TMM factors as
  `edgeR::calcNormFactors(method = "TMM")` on the samples' library sizes (their
  cells' stored `total_counts` summed, else their summed counts over every
  gene) and log-CPM as `cpm(log = TRUE, prior.count = 2)`; then `lmFit`
  (`~ group`, or `~ donor + group` paired) and `eBayes(trend = TRUE)`, the
  prior df and variance estimated from every expressed gene; the score is the
  moderated t, `log2FC` the group coefficient. A matrix that is not counts is
  fit on each sample's mean. Genes no compared cell expresses are left out of
  the fit and reported untested. The result names the library sizes' source,
  the TMM factors, the design and the prior (`about.values`,
  `about.model`). Every number is held to R edgeR 4.10.7 / limma 3.68.5 at
  10⁻⁶. Before, the samples were compared by Welch's t without TMM or
  variance moderation: with a few samples a side, a gene whose values
  happened to agree had a tiny variance and was called significant.
- **Paired samples.** Differential expression with Replicates (**Pair samples
  across groups**) and Composition (**Pair samples across conditions**) take
  the donors that gave cells to both groups paired: Differential expression
  blocks by donor in the linear model (`~ donor + group`, above); Composition
  takes Student's t of each donor's difference (`scipy.stats.ttest_rel`), and
  for its ranks the Wilcoxon signed-rank test (R `wilcox.test(paired =
  TRUE)`). Donors on one side only are left out, and a cell type whose cells
  are all in such donors is not tested (one Langerhans cell from a spleen-only
  donor read as "significant" from seven pairs without it). Each cell type is
  tested over the pairs holding any of its cells, and not below 3 of them
  (`too few pairs`, a chip, out of the multiple-testing correction): in a pair
  without the type on either side its two values are the samples' sizes, so
  myelocytes with cells in one donor's spleen read as adjusted p 0.03. A
  paired run needs 3 donors in both;
  every result says how many donors are in both, paired or not, so an agent
  learns when a paired test is possible (Suo: 6 donors gave liver and bone
  marrow).
- Composition's results CSV is the table an agent reads whole, from one
  builder: each type's cells and pooled proportion and its proportion in every
  sample beside the tests, with `n_samples_<condition>` (the samples each mean
  is over) and, paired, `n_pairs` (the pairs holding each type's cells, which
  its test took). Before, a paired
  run's file gave `n_samples_Liver 14` beside a paired t over 7 donors and left
  out the per-sample proportions; the p columns are `p_value` and `p_adjusted`
  as in the result.
- The dot plot names every gene and group without one on another in every
  browser: a name's room is the box the browser draws it in (in Firefox the
  names overlapped).
- The volcano plot labels, and its CSV lists, genes in the order its table
  ranks them when p-values tie (on hundreds of thousands of cells many
  underflow to 0): the stronger test statistic first. Its y axis reads round
  values beside the +∞ band.
- A dot plot names every gene and group, with no name drawn over another, at
  any count: the panel's preview grows to the rows' needs and the expanded
  view scrolls (the generic tick thinner overrode the dot plot's own naming).
- The expanded view of every analysis lists its display options beside the
  plot's (Dot plot, Composition and Gene signature had none there), and they
  stay through a change of the plot's options.
- An analysis's **Expand** opens once its result is drawn; a click as the
  preview appeared failed with "Expanded analysis view requires a completed
  analysis result".
- Detailed, Correlation and Gene signature draw each page in its swatch's
  color, whatever order the pages were picked in (they followed the pick
  order).
- Quick insights summarizes integer fields kept exactly (such as `n_genes`);
  it failed on them.
- A dot plot draws a gene that is in two gene sets in each set's column, as
  scanpy does: the second column was folded into the first (the plot showed
  23 columns for 24 genes and the second set's bracket and dots sat on the
  wrong genes).

#### Vector field

- **A field chosen while the overlay turned on no longer leaves it "loading"
  for good.** A session or an undo that put back an overlay showing another
  field than the page's, or an agent that switched the overlay on and then
  chose a field, left `Show overlay` and `Vector field` dimmed ("the vector
  field is loading"), no particle drawn and the switch refused. The newest
  choice now takes the load over and is the one shown; a restore puts the
  field back before the switch that shows it, so the switch loads it once.
- The only field bound to the embedding a view shows is selected; with
  several, Cellucid still asks (the dataset's default aside) instead of
  offering a one-item list to pick from.
- The block names embeddings as **View N shows** does ("Vector fields are
  bound to UMAP · 3D, UMAP · 2D. Show one of them in this view to enable."),
  not by their `obsm` keys.
- A field with no motion at all says so when it loads ("No particle moves in
  velocity on UMAP · 2D: every vector is zero."); it turned on silently and
  drew nothing.
- A field that cannot be read says so once; the loading notice failed with the
  same words beside the error.
- Cells left out for an infinite value read "(NaN or infinite in the file)",
  for vectors and positions alike; it said NaN.
- Particle size reads every half step apart (11.5, 12); a rejected Trail
  persistence value names it as the control does (it said "Trail fade").
- **A figure of a view with the vector field shown draws the field** (PNG,
  SVG, the export preview and an agent's export picture): a still picture of
  the flow, arrows on a grid as scVelo's `velocity_embedding_grid` (the
  Gaussian-weighted mean of nearby cells' vectors at 40 nodes across the
  plot), from the cells, vectors, camera and filters the panel shows, colored
  by speed with the overlay's color scheme, with a key naming the field and
  the speed range in the embedding's units. The export no longer refuses
  with "Velocity overlay not exported"; the figure's description and JSON
  record name the field each panel draws.
- The vector field's HDR composite keeps each color's hue: it tone-mapped
  channel by channel, so a bright purple flow (viridis's slow end) passed
  through magenta, and a teal one through cyan, on its way to white. The
  luminance is mapped now, the color scaled with it and desaturated toward
  white along its own hue only past the display's range; the contrast curve
  (which assumes display values) runs after the tone mapping, and the bloom's
  screen blend no longer flips sign on values past 1. Every pass reading the
  trail or the bloom also clamps it to a finite ceiling (a half-float sum can
  overflow at the dials' far ends), and both bloom blurs take their weights
  from one normalized Gaussian: the bloom's hand-typed kernel summed to
  1.0008, brightening the bloom 0.08 % each pass, peaked off its center and
  read two texels at weight zero. At every dial's maximum the flow is mostly
  white, ringed by its own colors at full brightness.
- `vector-field.summary` with no field named reads the one the overlay shows,
  else the dataset's default, else the only one bound to the view's
  embedding; with several and none of those it asks for one (it took the
  first by name). A field the overlay cannot draw is refused as unavailable.

#### Figure export

- **A figure is sized as it is printed.** You set its **width** — a journal
  column (**Single column · 89 mm**, **1.5 column · 120 mm**, **Double column
  · 183 mm**, the default) or millimetres of your own — and the figure fills
  it exactly, at the height its plots make it, with its legend, axes and title
  inside the width. Text is set in **points** (Auto text: 7 pt, ticks and
  centroid labels 6 pt, the title 8 pt, as journals ask) and keeps its printed
  size at any width and DPI. The line under the width says the printed size
  and a PNG's pixels (`183 × 130 mm · 2163 × 1541 px` (its pointer says the DPI)); an
  SVG states its size in millimetres, so editors open it at that size; the
  export picture (`view.picture` with source "export") says it too, and when
  a legend makes the figure wider than asked. Before, the size was a box in
  CSS pixels the figure was fitted in, and text in pixels: an agent set 1800
  for "180 mm" and got a figure 476 mm wide at 300 dpi, then did
  pixel-to-point arithmetic wrong twice (live-7 L25-10).
  `figure_export_set` takes `size_preset` and `width` in millimetres and
  `font_sizes` in points (version 3; a version-2 step's pixels are converted,
  as are the settings of a session saved by 0.9.7, except a size that followed
  the window — **Screen**, **Screen (half)**, 0.9.7's default — which was
  nobody's choice and opens at the default, double column). **Fit plot** is gone (the
  figure always takes the frame's shape); **Lock aspect** keeps the frame's
  own shape. A legend beside the plot is at least 10 ems of its text wide (it
  was a fixed 240 px, 63 mm, which left a single column's plot 15 mm wide),
  and a title or legend font is drawn at the size asked (it was floored at
  14 px and rounded to whole pixels). The preview's **Fit the cells** line
  now says when a plot is far from the shape of its cells
  (`Plot narrower than the cells:`).
- **The vector field's key is part of the legend.** It stood in the plot's
  corner whatever the plot's size: in a single-column figure its line ran off
  the figure's left edge, and wrapped it covered the cells. With **Legend** on
  it is now two rows of the legend (an arrow and what the arrows are, the
  colors' ramp and their speed range), laid out with it in the PNG, the SVG,
  the preview and the export picture; panels of one key share it, panels
  whose keys differ each carry theirs. Only a figure without a legend draws it
  in the plot, wrapped to a narrow plot's width. A single view's legend asked
  beside the plot that would widen the figure past its width (a single column
  with the key or long names beside a short plot: 101 mm for 89) stands below
  the plot instead, as a grid's legends do, and the export picture says so.
- **Fit the cells.** A figure keeps the screen's framing, so a narrow pane
  (the chat open, a grid of views) gave a narrow plot of a wide embedding,
  with empty bands above and below; an agent once hid the person's sidebar to
  widen the panes. **Fit the cells** (beside **Frame export**) frames every
  cell the figure draws with a small margin, so each plot takes the cells'
  shape; with **Export all views** one frame holds every panel's cells. An
  agent asks for it with `figure_export_set(frame = "cells")`
  (`figure-export.set` version 2). The preview says in one line when a plot
  is far (1.5 times or more) from its share of the size asked for, and offers
  the button.
- **An SVG is Full vector until you choose another strategy.** The first
  Export of an SVG was always refused ("Choose an SVG point strategy") and
  people and agents then chose Full vector every time; the strategy select
  now opens on it. `svg_strategy` is never null: a recorded step that left it
  unchosen replays as Full vector.
- **Each panel of an all-views figure is titled by its view's name**
  (`A. Early`, `B. Mid`, `C. Late`; a view without a name keeps its field and
  filters), and a named view exported alone is titled `<dataset> • <name>`
  until you write a title. The provenance names it too: each `views` record
  of the embedded JSON carries `name` (null when the view has none), the
  description's clause and the file name use it.
- **A figure's legend lists what the figure shows.** A category that a filter
  on another field leaves without cells is no longer listed as "(hidden)": He
  colored by LVL2 with only the epithelial cells kept listed 30 of its 41
  types as hidden when you had hidden none. Categories you hid keep their
  marked entry. A single view's automatic title names the filters on other
  fields as the view's label does ("Developing human lung (He) • LVL2 • LVL0:
  showing only Epithelial"); a grid panel's label already did.
- **Grid panel legends no longer shrink the plots to strips.** Each panel
  legend's side was decided against the plot the grid has without legends;
  kept beside a plot the legends then shrank, it was planned again in more
  columns, and two wide views in an 800 × 600 figure came out as 182 × 114
  strips beside 862 px legends, 1550 px wide. When the legends beside their
  plots do not fit the size, they stand below them (363 × 227 plots, within
  800 × 600).
- An exported view is framed exactly as on screen: the plot takes the view's
  shape inside the size you ask for, the cells reach its edges and the axes end
  where the drawn view ends (cells were missing in bands inside the axes). Axis
  ticks never collide on a narrow plot, and the orientation badge keeps its
  on-screen share of the plot.
- PNG figures record their DPI (`pHYs`): a 300 DPI export opens at its size
  in Word, Illustrator or Inkscape instead of at 72 dpi.
- A Hybrid SVG draws its points at 300 dpi at the figure's size (it was a
  fixed 2×, about 192 dpi).
- A figure's legend names every category whole: it is as large as its
  entries need and the figure grows around it (Suo's 69 `cell_type` names
  were cut with "…" at 720 × 500). In an **Export all views** grid, a panel
  legend that cannot name its entries beside its plot moves below it, rows
  grow to hold it, and the shared legend is sized to its entries (a grid
  panel's 69 entries were drawn as swatches and "…"). A legend below a
  single plot no longer covers the x axis. A PNG larger than the browser can
  draw is refused before rendering, with one sentence naming its size and
  the remedy, instead of a legend capped at the plot's size.
- A figure's default title names the categories only when some are hidden,
  as the Filters list says it ("cell_type: hiding B1 cells"); with every
  category shown it is the dataset and the field. A grid figure whose panels
  show different fields is titled by the dataset alone, and each panel names
  its own field.
- A view's label and a figure's filename name the colored field once
  ("cell_type: hiding B1 cells", not "cell_type • cell_type: hiding B1
  cells").
- The preview's status line is one short count (how the preview sample was
  drawn is in its title); the background choice reads "White".
- **An Export all views figure is a composed figure.** Every panel's plot has
  one height and its view's shape, the plots of a row and of a column line
  up, and the figure is exactly its panels, labels and legends at the largest
  size that fits the width × height asked for: one side reaches it, nothing
  pads the other. A figure of two tall views was the full 1400 × 800 with the
  panels as narrow strips in its corner, of different heights when their
  legends differed. Legends take room from the plots down to half the size
  they have without legends; past that the figure grows. A panel legend that
  would be wider than its plot stands below it. A legend's frame hugs its
  entries (it was stretched to the plot's height).
- **Axes only where they are rulers**: on a 1-D or 2-D view seen flat
  (planar). A 3-D view, or a flat view tilted by orbiting, has no axes and no
  axis margins: its ticks were camera units on one depth plane, read as UMAP
  coordinates. Its directions are in the **3D orientation** badge.
- Full Vector and Optimized Vector SVGs draw each 3-D cell at the size the
  screen draws it, near cells larger than far ones (size attenuation): they
  were a third to a quarter of the screen's spheres, all the same size.
  Centroid discs (every format) and the preview follow the same rule.
- A legend's entries stand at the same place in SVG and PNG, each label
  centered on its swatch; in a larger legend font the first entry rose into
  the legend's title.
- **A figure's legend lists the highlighted cells it draws.** The cells of
  the highlight page you are on are drawn over the coloring in the
  highlight's cyan, and nothing in the figure said so: the germ cells'
  highlight over NANOS3's colormap was read as "the cells with high NANOS3".
  The legend now carries a cyan swatch and `Highlighted: <page> (<N> cells)`
  wherever highlighted cells are drawn (each panel's own legend, or the
  shared one once), in PNG and SVG alike, and the layout makes room for it;
  a view colored by no field gets a legend of the highlight alone. Legend
  off, or no highlighted cell drawn, no entry.
- **The axis labels never run into the tick labels.** The y-axis label stood
  a fixed distance from the axis, so at 17 px ticks and a 20 px label
  "UMAP 2" touched "-0.5", and larger text ran over the ticks; the axis
  margins were 62 px whatever the text size, so large text was cut at the
  figure's edge. The axes are now laid out from their measured text: the
  y label past the widest tick label, the margins as large as the tick and
  axis-label sizes need (1 to 500 px), a corner tick label or an axis label
  longer than its side kept clear of the title, the legend and the other
  axis. The PNG put the y tick labels 2 px from their ticks and the SVG 8 px:
  both formats now draw one layout of the axes, and the selection badge's
  corners (rounded in SVG only) and the legend swatches' outline match too.

#### Reduced motion

- With reduced motion asked for, nothing transitions. Every element was given
  a 0.01 ms transition of every property, so each change held its old value
  for a frame: an open menu drew a frame at the window's corner before its
  place under its button, and the history tree a frame at its old height.
- With reduced motion asked for, nothing in the viewer moves by itself. The
  vector field overlay draws a still picture of the flow: every particle
  starts from a seeded cell and is traced 48 fixed steps at once, each faded
  as the moving trail would have left it, so each streak shows its cell's
  direction and speed; it is traced again only when the camera, the visible
  cells, the field or a setting changes, the same picture every time, and a
  still scene then asks for no frames (the particles animated without end).
  Highlighted cells keep a steady ring instead of pulsing, volumetric smoke
  stops drifting at any `Animation speed`, and the camera moves only while it
  is moved: a drag let go fast no longer glides on, a zoom arrives at once,
  and free-fly stops with its keys. Changing the system setting
  while the app is open switches at once, both ways. Every mover drawn by
  script reads the one preference (`utils/reduced-motion.js`; the camera
  path, the tour, the viewer and the backdrop looks each read it on their
  own).
- An agent is told when reduced motion is on. A camera path went to its last
  keyframe at once, but `camera_path_play` answered `{seconds: 18.3, ended:
  "end"}` as for a path that played, and an agent told the person "the film
  ran all 18.3 s". It now answers `ended: "reduced-motion"` with the seconds
  the camera moved (0; the seconds it played when the setting is turned on
  mid-path), and its step reads "Went to the end of the camera path at once
  (reduced motion)"; `views_state` says `reduced_motion` with what it changes.
- The documentation's velocity images are captured as that still picture and
  so repeat byte for byte (they were pinned to a number of frames and never
  did), and "Canceled after 1.5 s" in the canceled analysis image is the
  screenshot tool's held page clock, not the capture host's speed. The
  tool's own uniform pin for the highlight pulse is gone: the app draws the
  steady ring itself.

#### Camera path

- **Captions: a camera path shows its keyframes' names as titles** (owner,
  F110). While a path plays, the label of the keyframe the film is at
  appears as a compact caption low over the plot, from halfway into the move
  to that keyframe until halfway out of it, fading in and out (shown and
  hidden at once under reduced motion). Only names someone gave are captions
  (not the automatic "KF 3", not the return to start). **Captions** is a new
  playback setting, on by default and saved with the path in sessions
  (`camera-path.playback.set` `captions`; `camera-path.get` reads it); a
  picture of the screen draws the caption, a figure says which was up. The
  agent knows a label becomes a caption, the director's use. A path saved
  before captions existed plays with them.
- **The playback bar shows where the film stopped, and Play goes on from
  there.** After a film played to its end, an agent's **Stop the camera
  path** or **Stop** on the agent's turn, the camera stayed where the film
  stopped but the bar read 0:00. The bar now shows that place, and **Play**
  goes on from it (from the end it starts again), as it already did after a
  pause or a drag of the timeline. The bar's own **Stop** still returns to
  the first keyframe; an agent's **Play the camera path** still starts there.
- **Stop on an agent's turn stops a loop it started.** A path the agent set
  looping answered at once and played on after **Stop** on its turn. Work an
  action leaves going now belongs to the turn that started it (the action
  layer's `ongoing`), so the turn's Stop stops the loop where it is; a loop
  the turn left playing when it ended plays on until you pause or stop it.
- **A parameter refused for its form says the form.** A keyframe label with a
  space at its end, a view not named "View N" or a malformed color was
  refused "… does not have the expected form"; it now says "must be a name
  with no space at its start or end", "must be a view name such as "View 2"",
  "must be a color as #rrggbb", or the form by an example ("must be like
  "s3"") — for every parameter of every action.
- **Stopping an agent leaves its film where it stopped.** **Stop** on an
  agent's turn while it played a camera path put the camera back where the
  path had started (a stopped change is rolled back). The camera now stays
  where the film stopped, as a pause and **Stop the camera path** leave it,
  and the play's one step holds that camera, so **Undo** or **Undo turn**
  still puts it back. An action declares this exception in its definition
  (`keepsWhenStopped`); every other change stopped is still rolled back.
- **A cell near the camera is drawn up to 160 pixels across** (it was 128),
  or the GPU's own largest point where that is smaller. Near cells are not
  faded or shrunk: however near the camera comes (a camera path's close-up,
  **Frame cells with the camera**), every cell is drawn whole and opaque and
  can be picked, on screen and in every figure format. Every shader quality,
  the highlight rings, the category markers, picking and the figure export
  size a cell by the one rule (`rendering/point-sprite.js`), where four shader
  copies of the clamp and the highlight ring's own could drift apart.
- **Stop the camera path** stops a path that **Play the camera path** is
  still playing, at once: it waited behind the play (both changed the
  camera, so the stop queued until the path had ended) and was then refused
  "The camera path is not playing". A stop is now a control of work under way
  (a new effect, `control`, also of **Stop an analysis**): it records no
  step, the play answers `ended: "stopped"` with one step holding the camera
  where the path stopped, and a stop with nothing playing answers
  `{"stopped": false}` instead of a refusal.
- The camera path's bar says the path's length as soon as it is ready (it
  read 0:00 / 0:00 until played).
- **A path plays as it is timed.** Each keyframe is reached at its own time:
  easing is now how the path starts and ends (from rest, to rest), never a
  warp of the whole timeline. An eased path with 5 s transitions played its
  middle ones in 2.6 s and its first in twice its time, and the bar's marks
  were not where the keyframes were reached.
- The Catmull-Rom spline is timed by the keyframes and clamped where a
  coordinate turns: the speed no longer jolts at a keyframe between a short
  and a long transition, the camera never comes closer to the cells than the
  keyframes it moves between nor tilts past one over the top (a flip), and a
  position saved twice is a hold (the camera eases in, rests, eases out); it
  drifted away and back.
- A looping path that returns to its start at a constant pace (Linear) loops
  seamlessly: its last transition flows into its first (the speed jumped at
  the seam).
- **A keyframe holds.** Under each keyframe of the Camera path list a hold
  field sets the seconds the camera rests there before moving on (empty for
  none; the last keyframe's holds the end); the playback bar draws each rest
  after its keyframe's mark, and every transition keeps its time. Actions:
  `camera-path.keyframe.add` takes `hold` (version 2), **Set a camera path
  hold** (`camera-path.hold.set`) sets one keyframe's or every keyframe's,
  and **Read the camera path** reports each `hold`. A hold took a second
  keyframe of the same camera and a timing call (an agent's tour of 7 holds
  made 25 calls). Sessions keep holds; a keyframe without one is saved as
  before, so format 4 sessions of 0.9.7 load unchanged.

#### Agents

- **Every call of a reply shows how it ended while the turn runs.** The
  results of a reply's calls are kept in one message as each comes; the chat
  redrew only when a message was added, so a result after the first (a
  refused call beside a step) stayed "running" in the row and in the summary
  ("5 actions so far · 1 step", without "· 1 refused") until something else
  changed.
- **The agent is told what you did between its turns, also after a reload
  and when another history is shown.** The chat told the model only of steps
  on the same tracker, forgot it on a reload, and said nothing when you
  loaded a session file beside yours (an agent then took "same questions
  please" for a first look at a new dataset). Each message now carries a
  compact account from the history: your steps and its undone ones by label,
  or which history is shown now, where it came from, and its latest steps;
  what the agent last saw is kept with the conversation (a reload, a saved
  state or a shared session keeps it). Agents connected through Claude,
  ChatGPT or `cellucid mcp` are told the same with their next result.
- **Splitting one label is one call.** `field.create-from-pages` (version 2)
  takes `rest`, a categorical column whose category the cells on no page
  keep, and **Create categorical** offers **Keep their `<column>` category**
  for the cells on no page: an agent had made 45 temporary pages, a field of
  48 pages and deleted 56 pages (101 calls) to split one cell type.
  `highlight.page.delete` (version 2) deletes several pages in one step
  (`pages`).
- **The agent tool list no longer changes with your settings** (bug hunt
  live-8, L26-9): a connector's and a local `cellucid mcp` agent's tool list
  said per tool what your permissions made of it ("asks", "runs without
  asking", "refused") and what History settings leave out, so every change
  of a setting changed the list: the page sent `agent_tools_changed`, and
  `cellucid mcp` (Python and R) sent MCP `notifications/tools/list_changed`,
  also whenever a page took the viewer again. For Claude 5.5 models that
  rewrote the client's cached prompt (about 55k tokens) and dropped the
  conversation's earlier reasoning: measured on Opus 5.5, the connector's
  next turn after a permission change cost $0.41 against about $0.10
  cached, about 4×. The list is now the same whatever the settings, for
  every way in; the agent is told your settings ("[Cellucid: the person's
  settings for you now. …]") after its connection's first tool result and
  after any result once they changed, and the local instructions say so.
  The viewer channel's `agent_tools_changed` event is gone (viewer protocol
  fixture and endpoint events in all three repositories), and both packages'
  `cellucid mcp` advertise `tools.listChanged: false` and never send a list
  change.
- **Removing a connection stops its agent's call** (bug hunt live-8, L26-1):
  **Manage… › Remove** (and **Forget this browser**) while Claude or ChatGPT
  waited on a call left the call running in the tab: the removed agent's
  change was still recorded as its step and what the call read was still
  sent to it, up to the relay's 225 s wait. The relay now ends a
  connection's calls in flight with the connection (removed, forgotten,
  closed after a refresh token was used twice, expired): the tab is told to
  stop each, and the agent reads at once that the connection was removed
  during the call.
- **A tab reloaded during an agent's turn no longer ends the turn** (bug hunt
  live-8, L26-5): a call in flight when the tab reloads is answered at once
  (whether it was done is not known; read the history), but the agent's next
  call was told "No Cellucid tab is connected" 4 ms later, and Claude ended
  its turn asking the person to reopen Cellucid while the page was 2 s from
  being back. The relay now waits for a reloaded tab (up to 20 s after its
  last tab closed) and runs the call there.
- **An agent's change names its turn** (bug hunt live-8, L26-3): an agent
  told to end its turn with a note on the turn could not name the turn
  without reading the history first (Claude through the connector guessed
  `{"turn": "current"}`, was refused, read the summary, tried again). Each
  change an agent makes now answers with its step's turn (`step.turn`), the
  id a note on the turn or Undo turn takes.
- **A ratio of log-scale means is not an "X-fold" change** (F159): two live
  passes saw an agent call a ratio of means of log1p values "16-fold". The
  agent instructions now say that a difference of log-scale means is a log
  fold change and a ratio of them is none, and to compute a fold change on
  the original scale before saying "X-fold". A differential expression or
  marker result whose fold change is "of the means as stored" says that, on
  log-scale values, the ratio is no fold change of expression.
- **Claude's notes between its calls show in the chat** (bug hunt
  agent-loop lane): Claude Opus 5.5 and Sonnet 5.5 (and the Fable models)
  write their short notes between tool calls ("the clusters are read;
  coloring by cell_type next") as `thinking` blocks that come back empty
  unless asked for, so a long turn showed nothing but "Thinking…". The chat
  now asks those models for them (`thinking.display: "updates"`, Anthropic's
  beta `thinking-display-updates-2026-08-18`) and shows each note as a quiet
  line before the call it introduces, folded with that call's row when the
  turn ends and kept in the saved and exported chat. The notes are written
  and billed the same whether shown or not, so showing them costs nothing
  more (measured 2026-11-11 on a live 10-call annotation turn, the setting
  off then on: Sonnet 5.5 728k/1,333 and 725k/1,187 prompt/output tokens,
  Opus 5.5 803k/1,555 and 802k/1,555; the request is accepted from the
  browser, beta header included). Other
  models are asked nothing (a setting that comes with adaptive thinking would
  turn thinking on for a model that runs without it); a request refused for
  the setting is sent again without it.
- **An agent's question about heavy work no longer says "about 88 h" for a
  one-second analysis** (bug hunt live-7 L25-1): after the small live-agents
  fixture (6 genes, 120 cells), Marker genes on Kanemaru (131,636 cells,
  3,691 genes) was asked about as "about 88 h 36 min on this device" and
  Differential expression would have said about 102 h; both took about a
  second. What a run measured was scaled per gene and, from a single run, all
  per cell, so the fixed cost of a tiny run grew 600 × 1,100-fold. A run's
  measured cost now counts only for runs over the same genes, and one run
  measured on fewer cells than are open gives no estimate (two sizes give
  the line). The measurements kept in this browser start afresh.
- **A field made, renamed or deleted is at once a choice of every
  analysis** (bug hunt live-7 L25-6): the analysis forms never heard of
  field changes, so a field an agent had just made was "no choice" of
  Differential expression (and missing from a person's Field select) until
  another change redrew the form. Every analysis, in the panel and in each
  window, now reads the fields and their categories again on any field
  change.
- **Verify reproduces a branch with marker genes** (bug hunt live-7 L25-16):
  a marker genes result records when it was computed, how long it took and
  how this device's reads were tuned (`computedAt`, `computeDuration`,
  `tuning`), so every replay of it differed and Verify said "not reproduced".
  These are now run records, left out of the comparison as durations and the
  cache flag already were.
- **A parameter of the wrong kind is shown as it was sent** (bug hunt live-7
  L25-15): an object where text was asked for read "not [object Object];
  field is optional: leave it out", so a model that sent a field as
  `{source, key}` where an analysis option takes its name was told to drop
  it. The refusal now shows the object and, for a field, the name to give.
- **A conversation shared by the agent's own call no longer says that call
  did not finish** (bug hunt live-7 L25-14): a session file written by the
  agent's `history_share` (or a state saved by `history_save`) holds the
  conversation as it stood while that call ran, and answered the call "still
  running … what it did is not in this save"; the colleague's agent read it
  as "the session-file export didn't finish" of the file it was reading. The
  saving call is now answered as the save the conversation is kept in.
- **Summarize a field per category over a subset, the view left as it is**
  (bug hunt live-7 L25-13): `field.summary` (Python `field_summary()`, R
  `cellucid_field_summary()`) takes `within`: the cells of some categories of
  another categorical field, or in a range of a continuous field or gene.
  Agents hid the person's cells with filters to read donors within one cell
  type or cell types per age range, and put the view back after.
- **Paired Composition: the mean proportions and log2 ratio are over the
  pairs** (bug hunt live-7 L25-12): they were over every sample, so on Suo
  (fetal liver vs bone marrow, 6 donors with both) one donor with 2 liver
  cells turned LYVE1+ macrophages' ratio toward the liver while the paired
  test found more in the marrow in all six. They are now over the samples in
  both conditions, as the test reads them; each sample's own column still
  shows the rest.
- **A list of groups keeps one form** (L25-2): groups an agent named in its
  own order were rewritten into the field's order by the form's next
  refresh, which the history recorded as "Changed outside the history".
  Every list of groups (Marker genes, Differential expression, Gene
  signature, Dot plot, Composition) is kept in the picker's order.
- **Read the history's summary gives the notes in full** (L25-9): the
  tracker's, each branch's, group's and turn's notes, not their counts
  (`history_summary` / Python `history_summary()` / R
  `cellucid_history_summary()`; a turn's `summary` is now the first of its
  `notes`). A colleague's agent opening a shared session could not read the
  branch notes its sender had left.
- **Exclude genes says its values** (L25-5): gene families (mito, ribo,
  hemo); a gene by name goes in Exclude pattern, a regular expression.
- **The chat's agent says what it is about to do, and keeps you posted**
  (bug hunt F169): the chat's instructions now carry Anthropic's suggested
  line for a person watching: a line before the work starts, brief updates
  while it runs, a recap that stands on its own at the end. Measured on the
  four showcase tasks with Claude Opus 5.5 (two passes each way): every turn
  opened with a line (6 of 8 before), the longest stretch of calls with
  nothing to read on the annotation task fell from 4–6 calls to 2, and the
  cost did not change ($0.99 against $1.02 for the four tasks). Connector
  and notebook connections keep their instructions (their clients show
  their own progress).
- **An agent can name every embedding** (bug hunt agent-loop lane, L21-6):
  **Read the dataset summary** (`dataset.summary`, Python
  `dataset_summary`, R `cellucid_dataset_summary()`) lists every embedding the
  data mapping offers, on or off, with the obsm `key`, the `columns`, the
  array's `width` and whether it is `enabled` (and usable), the names
  **Show an embedding** and **Change the data mapping** take. Asked for a 1D
  embedding, an agent saw only "UMAP, 2" and guessed `X_umap`.
- **A refused embedding of the data mapping says what is missing** (L21-7):
  "This mapping enables obsm['X_gone'] columns [0, 1], which this dataset does
  not have."; it said "X_gone is enabled, but this dataset does not have."
- **Add steps to a group** (`history.group.extend`) to a step the group
  holds already says so (L21-8): a turn's later steps join its group by
  themselves, and extending the group to one was refused "Steps 2 to 4 are
  not consecutive steps of one branch", which they were.
- **Pause recording** no longer tells an agent to ask the person itself
  (L21-3): Cellucid asks, as the person's Permissions say; a chat model asked
  in its own words and never called it. The action reference's examples no
  longer merge a category into itself or rename to the same name (L21-5:
  `into`, the new `name` and a page's new name have examples of their own).
- **Where cells are heading, as numbers** (bug hunt live pass 5):
  **Summarize where cells are heading** (`vector-field.summary`, new; Python
  `vector_field_summary`, R `cellucid_vector_field_summary()`) reads a vector
  field (RNA velocity) per category of a field in the embedding a view
  shows: the speed against the average cell, how aligned the arrows are,
  their mean direction (in words on a 2-D view), and the categories the
  arrows lead into (each cell's arrow shared among its 30 nearest cells by
  how directly it points at them, as scVelo's velocity graph, here in the
  embedding). Asked "where are the cells heading, evidence from the data",
  an agent had answered "the tools don't give me the velocity numbers".
- **An agent sees the modalities of a MuData.** **Read the dataset
  summary** (`dataset.summary`) names the modality the genes are read from
  and every modality as Data mapping › Modality lists them ("prot · 5
  features · 36 of 44 cells"); opening a CITE-seq `.h5mu`, an agent had said
  it had "no tool that lists the modalities" and missed the proteins.
- **An agent continuing a shared conversation knows whose it was** (bug hunt
  live pass 6). A conversation received with a session file or a saved state
  is sent to the model with each received message marked as another
  person's, answered by their agent in their browser (and so summarized).
  Asked "what did my colleague find?", an agent had answered "all 37 steps
  were taken by me at your request" and told the colleague it was their own
  session. A received turn's **Sent** says nothing of it left this browser
  and names what the sender's browser sent, where it had said "your message
  … left this browser".
- **An agent sees the figure it is about to export** (bug hunt live pass 6):
  **Take a picture of a view** (`view.picture`) takes `source: "export"`, the
  figure exactly as **Export** would download it now with the Figure export
  settings (size, frame, background, fonts, title, panel titles, legends,
  every view while Export all views is on), scaled into the size asked for;
  nothing is written. Preparing a paper's Figure 2, an agent had told the
  person three times "my previews don't show the export settings, please
  check the file", and never saw a cut-off panel title.
- **A panel title is never cut off.** A multi-panel figure's panel title
  (its letter and its view's name) wider than its column widens the column,
  and the plots shrink to keep the figure in its box; "C. Neurog3 (endocrine
  pro" had been clipped at the panel's edge. A view's name is described to
  agents as printed after the panel's letter ("B. Early"), so names carry no
  letters of their own (a figure had read "A. A · Lineage"), and each Figure
  export text size says what it sizes (`base`: the panel titles; `title`: the
  figure's title).
- **An agent sees the flow it tunes** (bug hunt live pass 5): **Take a
  picture of a view** (`view.picture`) takes `source: "screen"`, the view as
  the screen shows it now, the vector field overlay and the connectivity
  edges included (a figure picture leaves them out and now says how to see
  them). Asked where the cells are heading, an agent tuned six overlay
  settings blind and left the flow a white glare; with screen pictures it
  checked each change.
- **A slider set by an agent, a session or a jump reads back the value it
  was set to, and applies it.** The vector field's sliders answered
  `setting.set` with the value before ("Intensity: set to 0.25" for 0.9, in
  the step and in what the agent was told), and Particle density never
  reached the viewer (the readout said 60K, 15K particles were drawn): a
  restore sent `input` alone, and those sliders apply on release. It now
  ends as a person's change ends, `input` then `change`.
- **Reading the community annotation reads the repository.** Connected but
  never pulled, **Read the community annotation**
  (`annotation.suggestions`) answered "no suggestions" for a repository full
  of them, and the agent had no way to pull. It now pulls first when this
  browser has not pulled the repository yet, and whenever it is asked
  (`pull`); `pull: false` reads what this browser has. Nothing is
  published. It also gives the repository's rule for the field
  (`consensus_rule`: the voters needed and the net share), which an agent
  had guessed for a methods sentence.
- **An agent's votes and suggestions of a turn are one commit.**
  `annotation.publish` (new) takes the votes and the new labels together and
  writes the person's file once, as published now plus all of them (one
  pull request update when they cannot write to the repository), asked once
  ("It publishes 3 votes and 2 suggestions (each with your upvote) to … as
  one commit"); the tool asks an agent to put everything of a turn in one
  call. Eight votes of one turn were eight commits ("Update annotations for
  @…" ×8) and eight questions. Every change is checked before the person is
  asked: one refused (a label already suggested, a suggestion that is not
  there, a vote given twice differently), none is published and the refusal
  names each. `annotation.vote` and `annotation.suggest` are removed.
- **The leader is listed first.** The voting dialog listed suggestions by net
  votes alone and the legend named every label of the top net vote, while the
  consensus line ranks by net votes, then upvotes: "Alpha (immature)" (net 1,
  ▲1) was listed above "Alpha cell" (net 1, ▲2), and the legend called them a
  tie the banner did not. One rank everywhere (the session's order): the
  dialog, the legend and an agent's read follow the consensus.
- **Publish after publishing from another device just works.** When your own
  file changed on GitHub since your last Pull, Publish stopped with "Remote
  user file changed since your last Pull. Pull first…"; it now pulls (your
  file there and what is here are merged) and publishes again by itself, as
  an agent's publication did. An author with repository settings or merges
  not yet published is still told to Pull first, since a Pull brings the
  repository's settings over theirs.
- **A failure of Pull or Publish is said once, in words.** A refused sign-in
  showed "GitHub PUT /repos/…/ghid_….json failed with HTTP 401" beside
  "Your GitHub sign-in expired or was revoked…", and a Pull's failure showed
  the same sentence twice; the card of a failure said in its own sentence is
  now dismissed. A write whose outcome is not known says "The connection
  dropped while GitHub was writing, so whether the write was made is not
  known." (not "GitHub mutation outcome is unknown after a transport
  failure").
- The voting dialog's text boxes no longer show "0 characters" under every
  suggestion; the count shows while a box holds text.
- The consensus rule is stated whole: no voter yet is Pending too, and a
  tie is equal net votes and equal upvotes (docs, the reference repository's
  README and the agent's `consensus_rule`).
- **An agent continues your saved work** (`history.saved.open`, new): it
  opens a saved state of this browser by its id as Load…'s **Open** does,
  its dataset first when another is shown, then its history, which goes on
  where it left off. Asked to "continue my saved work", an agent could only
  replay the autosave into a new tracker (He, live). A whole-file download of
  its dataset is asked about as `dataset_open`'s is; a saved state it cannot
  open is refused with why (a file of your computer the page offers you to
  **Reopen**, a state another version kept, what failed), and the chat turn
  goes on with the history it opened instead of stopping as for a load by
  someone else. Load…'s **Open** no longer says "Opened …" when nothing
  opened (a newer choice, a failure, a file to reopen).
- **Highlight categories** (`highlight.category`, version 2) takes
  `categories`, one or more labels of one field, each made a group of the page
  in one step: one act, one undo, one step named for them ("Highlighted 10
  categories of donor (48,210 cells) on page …"). A page of ten donors took
  ten calls, ten steps and ten undos (He, live). Nothing is added when one of
  them has no cell visible; one already on the page is brought up to date in
  its place. A history of version 1 replays (`category` becomes
  `categories: [category]`); the Python and R wrappers
  (`highlight_category(field=..., categories=[...])`) follow.
- A permission card names what would be done, never as done: "Claude asks
  to export an analysis table: markers." (the action and the names it is
  given), where it quoted the step's past-tense sentence ("asks to export
  …: “Exported …”") for an act not yet allowed. Every kind of asking reads
  the same way, and so does the question a local app (Claude Desktop) asks.
- The chat's **Token budget per chat** is a cap: the chat makes no call that
  would pass it. Each call reads the whole prompt again, and the budget was
  checked against what was spent alone, so a chat with a 150,000-token budget
  ended at 209,751 (the third call, of 70,000, was made at 139,000). Before a
  call the chat now adds what it will read (the last call's prompt and what
  was added since) and stops before it, saying how much it would have read.
- The legend's highlight button names a category by its label, as every
  other legend action does: on a field of numbers or of true/false it sent the
  raw value, and the highlight was refused ("There is no category "1" in …").
  `highlight.category` takes labels only (text: "3", "true").
- **Read the figure export settings** (`figure-export.get`, `figure_export_get`):
  every Figure export setting as `figure-export.set` names it. An agent that
  set up a three-panel figure and came back to it after a dataset switch said
  "I can't read export settings back" and asked the person to check the panel.
- An analysis result says what `top` counted: "the top 40 genes by score up in
  A, then the top 40 up in B: top counts each direction" (Differential
  expression), "each group's top 10 shown markers" (Marker genes). An agent
  asked for 40 genes each way first read all 3,753.
- A table in the chat breaks its cells between words, never inside one: a
  six-column result table in the dock read "Adjuste|d p" and "0.04|%"; the
  columns of numbers now keep their words whole and the long first column
  wraps.
- When a long answer ends, the chat shows it from its first line, not its
  last; a person who scrolled during the turn keeps their place, and a new
  question is followed to the end.
- `dataset_open` takes the cells to open as `cells_open` does (`count`, or
  `fraction`, new in both, with the method, seed and the method's settings):
  "open Suo at 10%" reads Suo once at 10%, where it opened the default view
  and read again. Its `total_cells` is the dataset's, not the cells open.
- A column made from one whose groups an agent drew by naming cells (grouped,
  or a file's column edited by a merge or a deleted category into a new one)
  counts as such a column, as a copy and a rename did, so "names and
  summaries only" covers it too.
- The open chat is a side panel on a wide window: full height on the right,
  and the plot moves over to make room for it, so the agent never covers the
  cells it is changing (it opened over the foot of the plot). Its left edge
  sets its width (drag, or ←/→; a double click resets it), remembered in this
  browser; a narrow window keeps the chat over the plot's foot. The agent's
  calls one after another read as one message, its name said once.
- **A research answer is no longer buried under its calls** (bug hunt
  L5-13: 30–90 rows above the answer, every read saying only its action's
  title, "Read the categories of a field" four times). Each row now says what
  it was about: a read by its action's own sentence from its arguments
  ("Read the categories of cell_type", "Looked up 12 gene names: SOX9,
  SFTPC, ETV5 +9"), a running, refused or stopped call by its title and what
  it named ("Run an analysis: differential · refused"), a change asked for
  what was so already "· already so". A turn with more than three calls
  folds under one summary row ("Worked through 42 actions · 7 steps · 3
  analyses", a disclosure); once it ends its answer is what shows, with its
  refused or stopped calls (but not one the agent then made again and got
  through) and its questions; while it runs, its latest three calls. Show, Sent, Undo turn, Copy and Export are unchanged, and the
  exported chat names every call (`app/agent/chat/work.js`) by what became
  of it: a step, a read, done (a file written, a note, a stop: no step, but
  not a read), no change (so already), refused, or stopped by you (was
  "refused").
- A playing camera path is no longer recorded as your camera moves: an undo,
  a jump or an agent's turn while a path played (a loop for a talk) recorded
  the moving camera as **Moved the camera of View 1** by you, so an undo
  undid that step and the next one branched the history. Nothing is recorded
  while a path plays; it is one step where it comes to rest. An agent's
  looping play is a step when it starts (it was answered "it was so
  already").
- On a narrow window the chat folds to its **Ask an agent** button while a
  camera path plays (Play, or the agent's **Play the camera path**) and comes
  back when the path stops or pauses, its message field focused again if it
  was: it covered the lower part of every frame of the film. It stays open
  while you are writing a message or a question of the agent waits, and
  opening it during the play keeps it open until the path stops. The side
  panel of a wide window stays.
- Chat settings › **Prompt cache** (Claude): Kept 5 minutes (the default) or
  Kept 1 hour, remembered in this browser. An hour costs more when a chat
  starts and spares the full price of the instructions and tools after a
  pause of more than five minutes.
- The empty chat's first prompts name the right fields: labeling the
  dataset's clustering with cell types (`Label each cluster in "leiden" with a
  cell type`), or, with no clustering, the marker genes of its cell types —
  never "annotate" a field that already names the cell types (Pancreas offered
  "Annotate the clusters of cell_type" beside `clusters`).
- The picture an agent asks of a view (`view.picture`) reports the plot's real
  size: width × height is the box the plot fits in, keeping the view's shape
  as framed on screen (it reported the size asked for).
- A change an agent asks for that is so already (coloring by the field the
  dataset opened colored by) is answered with why no step was made: "this
  changed nothing: it was so already", or that the change is not recorded
  (recording paused, a kind your History settings leave out, a merged change
  put back). The agent said "I colored the cells by cell_type" for a change
  that made no step and showed no **Undo turn**.
- A turn that ends with no words for you (a model may end after its
  reasoning alone) says so in a note with **Try again**, instead of leaving
  the turn without an answer.
- Undoing a turn that made no step (it only read, or changed nothing) says
  so: "The tracker … has no record of the turn …: a turn is recorded with its
  first step, so one that only read or changed nothing has nothing to undo";
  it said "There is no turn" and listed every turn's id. Adding a note to
  such a turn says the same.
- Your next message tells the agent what changed in the history since its
  last turn: its steps you undid (**Undo turn**, kept for redo) or removed
  (only its steps, to Recently deleted), and the steps you, a notebook or
  another agent took, by their labels. After **Undo turn** it said "The Ductal
  cells are still hidden from the earlier step". A branch switch is said as one
  ("the branch shown is now …"; its steps stay on their branch), not as an
  undo. A step you take while its turn runs (coloring mid-turn) is told too,
  as taken while that turn ran (bug hunt L18-4: it counted as already seen,
  so neither that turn nor the next heard of it).
- The **Undo turn** question counts the same chat's later turns as such ("9
  from 3 later turns of this chat"); it said "9 from another Claude turn".
- A reload between two tool calls of one reply keeps the result of the call
  that finished; it was answered "The page closed before this call finished".
- An agent's **Stop** of its looping camera path records where the camera
  rests as a step of its turn (Undo turn takes it back), not as your move.
- Each model call carries ~2,000 fewer prompt tokens (63,848 → 61,898 before
  the conversation, measured with Anthropic's token count): a change that is
  one undoable step no longer says so in its tool's description, since the
  agent instructions state the rule once and every exception still says
  itself. The same holds for agents connected through the connector and
  `cellucid mcp`.
- An agent can export an analysis's result table (`analysis.export`; Python
  `analysis_export`, R `cellucid_analysis_export()`): the file the expanded
  view's CSV button writes (differential expression every gene; marker genes
  the shown markers or, with `all`, every gene of every group), after you
  allow writing a file. Asked "export that table", it had to say it could
  not. A result that is out of date is refused, saying why.
- An agent merging a category into another can name the merged category
  (`category.merge` `label`, version 2): "merge A into B" keeps B; it was
  named "merged A + B" (the legend's drag naming), so grouping cell types
  into a lineage took a rename after every merge.
- An agent's run of an analysis the Analysis panel does not show (it shows
  one at a time) shows it first, as you open an analysis to press its Run:
  setting up Differential expression and then Composition and running both
  was refused "not shown, open it first". Changing the options of an
  analysis with a result does the same: after Marker genes, setting
  Differential expression's groups was refused, `show: true` or not.
- **An agent reads how an analysis was tested.** **Read an analysis result**
  (`analysis.result`) gives `about` with the rows: for differential
  expression the two groups and their cells (with **Replicates**, their
  samples and those left out), the test, the correction, the expression
  matrix (`· counts` when it is) and the values tested, the fold change and how many
  genes pass the shown cutoffs each way; for marker genes each tested
  group's cells and what it is compared with; for gene signature each set's
  genes used and not found and the tests between the pages (t-test and
  Mann–Whitney for two, ANOVA and Kruskal–Wallis for more). An agent
  guessed the values were raw counts from their being whole numbers, and
  could not say whether a score differed between pages.
- The composition result an agent reads gives each cell type's cells, its
  pooled proportion and its proportion in each sample (or condition), and,
  tested, its test and p values, with the samples per condition: without a
  condition it gave the cell types' names alone.
- The count of significant genes in Differential expression's Statistics
  applies **Min % expressing in either group**, as the status line's count
  does: the two differed.
- An agent checks a list of genes in one call (`genes.search` `names`: the
  ones the dataset has, as it names them, and those it lacks): checking a
  marker panel against a sample of 5,000 genes took one search per gene,
  31 in one turn.
- **Summarize a field per category** (`field.summary`) leaves out the
  categories with no cells (counted in `empty`) and takes `categories` to
  summarize only those: on a view filtered to the epithelium of the He
  sample each summary by cell type listed about 100 empty groups.
- An analysis option's parameter says what its values are where the label
  does not: the cell scopes (`all`, `visible`, `page:<page name>`, now by the
  page's name; Within also `category`), **Replicates** (pseudobulk),
  Composition's **Samples**, **Condition** and **Test**, Gene signature's
  **Method**, **Reference group**, **Expression matrix**. An agent guessed
  Within's values twice ("Within has no choice categories").
- **An agent can direct a camera path.** **Frame cells with the camera**
  (`camera.framing`; Python `camera_framing`, R `cellucid_camera_framing()`)
  says where the whole view, each category of a field and the categories
  named together sit (center, extent) and gives a camera framing each from
  an azimuth and elevation, filling a share of the view's frame; **Play the
  camera path** / **Stop the camera path** (`camera-path.play`,
  `camera-path.stop`) play it as Play does, the playback one step; **Read
  the camera path** (`camera-path.get`) gives its keyframes, their times and
  its length. The agent instructions teach the craft: the story from the
  data, an establishing shot, close shots from changing directions, holds, a
  two-shot between far close-ups so the camera rises over the cloud, a
  pull-back. Asked for a documentary flythrough of Pancreas, Opus 5.5 took
  15 pictures to triangulate where the cell types were, still put close-ups
  inside the cloud, and could play only by turning Autoplay on (saved with
  the session).
- An optional parameter an agent sends as an empty stand-in (`camera:
  "null"` or `{}` for "the current camera") is refused with how to fix it:
  "camera is optional: leave it out (send no value at all, not null, "null"
  or {}) for the selected view's camera now", for every optional parameter
  of every action.
- An agent names Differential expression's Group A and Group B by the
  highlight page's name, as every other action names pages: the options read
  gave the page's internal id ("restof__page_1") and a page's name was
  refused ("Group A has no choice …"), so the agent had to guess "page_2".
- An agent's dataset summary gives the dataset's own description (the
  card's paragraph) and its vector fields (RNA velocity, drift): on the
  scVelo pancreas sample an agent concluded "there is no RNA velocity in this
  file" and gave the tissue from its own knowledge.
- An analysis run an agent was waiting for that is stopped (Stop, or Cancel
  in the panel) is told to the agent as stopped: it was told "this changed
  nothing: it was so already", which stayed in the conversation.
- A turn the page was reloaded in says "The page closed during this turn, so
  it stopped here." with **Try again**: reloaded between two calls, it ended
  with no word and no mark, its work half done.
- The chat's selection chip attaches highlighted cells once: a selection
  already sent with a message, or highlighted by the agent in its turn, is
  not attached again (the agent was told with every message that "the person
  attached their selection", its own highlight among them).
- **Stop undoes what it stopped.** A change an agent's step was making when
  you pressed Stop is put back and leaves no step, also when the change
  finished regardless (an analysis's display option waiting for its redraw):
  it was kept and recorded. A read stopped ends at once.
- **Stop ends promptly, and the work stops with it** (bug hunt F147). A
  stopped change is put back only once it has ended, so a step that waited
  on work not tied to Stop held every later Undo, jump and action until that
  work finished: an analysis's options, pages, view or new window waiting
  for the update they set off (a re-run reading thousands of genes), a run
  stopped before it had started, the cells read again for **Cells to open**
  or for an Undo that puts other cells back, genes read for more cells, a
  view switched to another embedding waiting for its coordinates, a field
  copied while it was read, a view kept with a branch's state, and the
  history report showing each step for its pictures. Each now stops its work
  when you press Stop — the analysis's update or run is canceled, the
  reading is canceled as its Cancel would and the cells open before stay,
  the report comes back to the step shown — so Stop answers at once and the
  history is free again.
- The Cellucid relay (`mcp.cellucid.com`) no longer redirects a malformed
  connection request of an app other than Claude or ChatGPT straight to that
  app's address: it shows the error with a link back to the app. Anyone can
  register an app with an address of their own, so the relay's address
  could be used to bounce a browser to any site (RFC 9700 §4.11.2).
- The relay keeps a digest of a connected app's id rather than the id
  itself (a registered app's id is up to 8 KiB of its own choosing, kept a
  month per connection); a tab names an approval only by a well-formed
  connect code; and its `subscriptions/listen` answer carries the same
  headers (CORS for a Cellucid page, `no-store`, `nosniff`) as its others.
- **Connect Claude?** warns when you came by a link too: anyone can start a
  connection from their own Claude or ChatGPT and send you its link, and
  **Allow** would connect their agent to your browser. A Cellucid page inside
  another site's frame drops a link's connect code and never asks.
- **A chat's Markdown export runs nothing in any Markdown viewer.** Text the
  export kept as code where a viewer does not (`` `<img …>`` `` closed by two
  backticks, a code span opened on the line before, fences out of step, a
  fence ended by its list) reached the viewer as live HTML, and link
  definitions (`[a]: javascript:…`) and images were not looked at. The export
  now writes the reply as the chat draws it — its code, bold, italic and web
  links — with every other character that starts Markdown escaped, and each
  fenced block closed by a fence no line of its code can match.
- **Double underscores inside a word stay as written, in the chat and its
  export.** The chat read `__` as bold anywhere (`snake__case__name`,
  `window.__xss`), so a reply showed part of a name in bold and the export
  wrote it back as `**`, changing what you and the agent had written; like a
  single underscore, a double one marks bold only at a word's edge now, as
  every CommonMark viewer reads it.

#### Python and R viewers

- **The viewer names what opened it.** Ending `cellucid mcp --open` froze the
  page with "The notebook or R session that opened this viewer stopped", the
  words every viewer said whatever served it (`cellucid serve` in a terminal,
  a script's `AnnDataViewer`, an agent app's MCP server). `/_cellucid/info`'s
  `runtime` now tells them apart -- `jupyter` (a notebook kernel), `python`
  (a script or session), `cellucid-serve`, `cellucid-mcp` (`cellucid mcp
  --open`), `r` (`cellucid_view()`, `cellucid_serve()`) and `r-mcp`
  (`cellucid_mcp(open = …)`) -- and one table of the web app
  (`VIEWER_HOST_RUNTIMES`) gives each one's words, label and language: every
  sentence of the viewer channel ("The cellucid serve command, which opened
  this viewer, stopped (its server no longer answers). … Run cellucid serve
  again."), a failed host command ("Python command failed"), a selection
  the host was not told of and the Dataset box ("cellucid serve") name the
  host by it, and "the program serving this viewer" before it said. History
  steps stay labeled `python` or `r`. The shared viewer protocol file lists
  each runtime with its label and language (SHA `28dbb572…` in all three
  repositories).

#### Security

- **GitHub sign-in takes only the answer to a sign-in started in this tab.**
  A link carrying someone else's GitHub token in its address
  (`#cellucid_github_token=…`) signed the tab in as them, so your votes and
  annotations would have been written as their account. **Sign in** now keeps
  a random state in the tab (`cellucid:github-app-auth:pending`) that the
  answer must name, once.
- A group, category or field named `__proto__` is kept as itself where results
  and saved state are rebuilt (what an agent is sent with **Agents see names
  and summaries only**, a saved analysis, a worker's result, a figure's
  legend): it vanished, its fields read as the record's own.
- A Details row links only to an http(s) address. A test lists every place the
  app's code turns a string into HTML, code or an opened address, with why it
  is safe; a new one fails until it is listed.
- **Only the browser that shows a connection's waiting page can approve it.**
  Anyone could start connecting their own Claude or ChatGPT, send you the
  page's **Open Cellucid to approve** link (or its code), and one **Allow** in
  your tab connected *their* agent to your Cellucid, which could then read
  your dataset (R9Z-1, the remote phishing of a device code). The relay's
  waiting page now leaves a cookie in its browser
  (`__Host-cellucid-approver`: `HttpOnly`, `Secure`, `SameSite=Strict`, ten
  minutes), and the tab claims a code with a request that carries it
  (`POST /oauth/claim`, no longer on its socket; relay protocol 2): a link or
  code from another browser is refused in your tab (*This connect code was
  started in another browser…*). **Code from another device** is now
  **Connect code**, for a tab of the same browser; a local relay (`wrangler
  dev`) sets the same cookie, kept by Chromium and Firefox on loopback
  addresses when the page and the relay share their host name.
- **An app the relay cannot vouch for is named by its own site.** An app
  that called itself "Claude" from `evil.example` asked as *Connect Claude?*
  and its steps were Claude's in the history (R9Z-2). Only Claude and ChatGPT
  are named so, when `claude.ai` or `chatgpt.com` vouches for them; any other
  app is its site (or *an app on this computer*), with the name it gives
  itself in quotes (*evil.example, which calls itself "Claude"*), on the
  waiting page, in the approval, and by one name (**Manage…**'s, *evil.example
  (calls itself "Claude")*) in its activity line, confirm cards and history
  steps, so two such apps connected at once are told apart (L26-4).
- **An agent asks before opening a link you have not seen in this tab.**
  `dataset_open` with a link, an accession or a server ran as an ordinary
  change, so an agent misled by text it read could make your tab request any
  address, carrying what it read to any host (R9Z-4). A new permission kind,
  **Opening links not seen in this tab** (Ask in both presets that let agents
  change anything, refused in Read only), asks before anything of the address
  is requested, unless you typed or pasted it into the address bar or the
  chat, it is open now, or it is in **My datasets**. The card names the
  address; the agent waits and is told why. A connected agent is told what
  your setting makes of it with your other settings.

- **Plain HTTP only on a loopback host a browser trusts.** Whether a host is
  this computer was written twice with different answers: the Worker's
  `ALLOWED_ORIGINS` rule (and the relay's copy) did not lower-case a host or
  know a bare `::1`, the app's own rule also counts `0.0.0.0`. Both now come
  from `utils/loopback.js`: `isSecureContextLoopbackHost` (`localhost`,
  `*.localhost`, 127.0.0.0/8, `::1`, the hosts browsers make a secure context)
  for where a secret may go over plain HTTP, `isLoopbackHostname` (also
  `0.0.0.0`) for "this computer". A development relay address at
  `http://0.0.0.0` was accepted for the tab's secret; it is refused, as the
  Worker always refused it in `ALLOWED_ORIGINS`.

#### Accessibility and keyboard

- **A dialog under a confirmation is out of reach on every browser.** A
  modal dialog opened over another makes the one below inert; Linux's WebKit
  still let focus into a dialog under a confirmation. The app now marks a
  dialog under the topmost `inert` itself and gives it back when it is the
  topmost again.
- **Focus moved after a dialog closes stays where it went.** A dialog closed
  with Escape handed focus back to its opener a task later, taking it from a
  control focused meanwhile: in the history, F2's rename closed by Escape and
  then Delete on a step row acted on the branch row instead (an intermittent
  test failure under load).
- **Keys work after a slider, a checkbox or a button.** After moving
  Navigation › Keyboard speed or ticking a box, W A S D, H, R, O, P, G and the
  other shortcuts did nothing until the plot was clicked: six key listeners
  each kept their own "is the person typing?" rule, five of them treating
  every input as a text field. One rule now decides who a key belongs to
  (`utils/key-ownership.js`): a text field, a dropdown (its letters choose its
  options), a menu, a list or a dialog keeps its keys; every other control
  leaves them to the page. The workarounds that took focus away after a
  checkbox or the Mode dropdown changed (a keyboard user lost their place) are
  gone; a mode chosen from the dropdown with the pointer still hands the keys
  to the plot. A key a focused control already took (the History list's `[`)
  no longer also zooms the camera.
- **One table of every key.** The Keyboard shortcuts panel, the page's key
  listener and the new docs page *Keyboard shortcuts* (accessibility section)
  are made from one table (`ui/onboarding/shortcut-table.js`); the camera keys
  come from the ones the viewer moves by (`rendering/navigation-keys.js`).
  The panel gains what it left out: Ctrl+Y (redo, off a Mac), the History
  list's keys (Home/End/Page keys, Backspace, the Menu key and Esc too) and
  the keys that move a panel, a page or a category (pick up, aim, drop).
- **A camera key moves the camera for as long as it was held.** The viewer
  read "is the key down" at each frame, so on a slow frame (a large dataset
  on a modest GPU) a quick tap could move nothing and equal presses moved the
  camera by different amounts; each key is now credited with the time it was
  held in each frame.
- **`?` with the sidebar hidden records no "Changed outside the history"
  step.** The panel opened in the same moment as the sidebar's action; it
  now opens once the sidebar is shown.
- **The graphics-reset notice is an app dialog**, as every other: the renderer
  reports the lost context and the app shows it (`#graphics-reset-dialog`);
  it keeps Escape and the backdrop while the graphics are restored. The
  analysis panel's field picker is a non-modal popover (no `aria-modal`, no
  Tab trap).

- **The colorblind preview shows what a colorblind reader sees, and no page
  claims more than that.** Figure export's Colorblind preview now applies the
  full-severity model of Machado et al. (2009) in linear light; the old matrices
  on encoded sRGB showed vivid cerulean and iris purple (the default palette's
  categories 2 and 10) as plainly different to a deuteranope, who sees one
  color, so an unreadable figure passed the preview. The accessibility page,
  the accessibility guide and the palette said the first 24 category colors
  stay distinct for red–green color blindness; measured with the same model,
  several pairs look alike (blues and violets, orange and amber). They now say
  the category colors are not chosen for color blindness, and how to check and
  change them; a test holds the claim to the measurement.
- **Tab leaves an open menu where you were.** In every menu (**Open…**, **My
  datasets**, the History panel's ⋯, a history row's or branch's menu, the
  history search's **Filter**, a highlight page's **Combine**), Tab sent focus
  to the top of the page, the plot: an open menu sits above the page, at its
  end. It now closes and goes on to the control after the menu's button
  (Shift+Tab, the one before it). The highlight pages' **Combine** menu is the
  shared menu now, with its look, its keys and its touch size.
- A menu button says whether its menu is open: the history search's **Filter**
  always said closed, and a branch's ⋯ and the tracker said nothing.
- The color swatches in the legend, on highlight page chips, in **Create
  categorical** and on an analysis's page chips show the focus ring: the
  color input under them took focus invisibly.
- `?` with the sidebar hidden shows it; `?` moves focus to the **Keyboard
  shortcuts** panel, and its scroll no longer animates with reduced motion
  asked for. A way into the sidebar from outside it focuses its control once
  the sidebar is shown: with the sidebar hidden, a gene sample's **Genes read
  for** link (from an analysis window) and Data mapping's way to its choices
  opened the sidebar but focus reached nothing (it was still inert).
- **Create categorical**'s rows name their controls by their page ("Category
  label for Page 1", "Highlight color for Page 1", "Remove Page 1"): the label
  field had no name, and every row's remove button was "Remove this page".
- The welcome card's quote and its attribution are readable in the light
  theme (2.87:1, now above 4.5:1); the welcome card, the dialogs and the
  license are measured for contrast as the sidebar is.
- The dataset's **Share** is **Share…**: it opens a dialog, as every other
  button that does.
- The × that removes a kept view, a highlight page or a gene from a set, and
  a page's ⋈, are quiet by a readable gray instead of being faded (the view's
  × read at 3.22:1).
- A history row's buttons keep one line, the row's first, at its end, in the
  same place on every kind of row: a turn's ↗ and details went onto a line of
  their own under a long label.
- **Every dialog is one kind of dialog.** The welcome card, the license, every
  confirmation, the figure export's warnings and the community annotation
  dialogs are now the browser's own modal dialog, as Cite, Share and the
  history's dialogs already were: the page behind is out of reach for the
  mouse, the keyboard and a screen reader alike, a dialog opened over another
  (a confirmation over a community annotation dialog or over the branch map)
  is the one in use and gives focus back into the one under it, Escape and a
  click beside the card close it, and focus goes back where it was. Five
  hand-built versions are gone, with their differences: Tab past a dialog's
  last control now goes to the browser's own controls and back in, as on any
  page (it went round inside before in some dialogs and not others); in
  Safari, Tab reaches buttons and links with Option held or full keyboard
  access on, as everywhere in Safari; a text selection dragged out of a card
  no longer closes it; and a history question asked from the branch map or
  Recently deleted opens over it instead of closing it and opening it again.
- **Keys pressed in a dialog stay in it.** With a dialog open, W, A, S, D
  moved the camera behind it, X cleared the highlight, H hid the sidebar, R
  reset the camera, and Escape left a selection or a history tour instead of
  closing the dialog. The notification cards stay above every dialog and
  usable there.

#### Controls

- **One wording while no dataset is open**, from one module
  (`app/ui/no-dataset.js`): "No dataset is open." and "Open a dataset to
  keep views." everywhere, where panels said "No dataset selected", "No
  dataset loaded", "Load a dataset to …" (Load is the session file's button),
  "Open a dataset before …" and two wordings of the history's line. Every
  analysis (and its copy into a window) says "No dataset is open. Open one in
  Session › Open data to analyze it." where gene analyses said "Requires gene
  expression data.", Composition "Requires a categorical field of cell
  types.", and Quick, Detailed and Correlation opened empty forms; KNN drag
  says there is no dataset instead of "No neighbor graph available for this
  dataset", and a dataset without a graph now says what to use instead
  (Proximity drag or Lasso).
- Camera path **Clear all** is dimmed with no keyframe to clear; a click did
  nothing.
- **One read of the neighbor graph is one notice.** Choosing KNN drag (or
  showing connectivity edges) showed "Loading neighbor graph for KNN mode"
  (or "Loading connectivity data") beside the read's own "Reading the
  neighbor graph · 27 MB read". The read now has one card, whoever asked for
  it, with its progress and its end ("Neighbor graph read: N edges"); a
  failure is said there once.
- With **Select mode** on, each selection tool's step line says that
  Intersect, Add or Subtract under the tools combines each gesture, instead
  of the `Alt` keys, and follows the switch at once.
- **Copying always says how it went.** Every Copy (the share link, badge and
  embed, the citation, a chat answer, the connector address, the history's
  text and code, a subsample link, the situation report, Marker genes' **Copy
  genes**, the annotation share link) goes through one helper: copied, or why
  not ("the browser did not allow it (the page must be in front and focused)")
  and, where it helps, the text to select by hand. **Copy genes** said nothing
  either way.
- Whether this computer is online is read in one place (the book of hosts):
  the chat bar, community annotation, the readers' faults, the renderer's
  index build and the situation report all agree.
- **The analyses' expanded view is a real modal dialog.** Its backdrop covers
  and closes it, yet the page behind stayed reachable by keyboard and screen
  reader (it said `aria-modal` over a hand-made Tab trap). It is now a
  `<dialog>` opened as every dialog of the app: the page behind is inert, Tab
  and Escape are the browser's, focus goes back to what opened it; a
  multi-select list inside it is no longer cut off from it. Closed, it lets
  the page go at once and fades out inert (it stayed a modal for its fade,
  the page unusable and its own buttons in reach after you closed it).
- A link's malformed subsample is said once: "The link's cell subsample was
  left out: "seed" needs "cells" …" (the reason repeated its own "The link's
  subsample:").
- **The empty viewer says what it waits for.** With no dataset open (None
  chosen, a link opened without, an opening canceled) the viewer was a blank
  plane; it now reads "No dataset open — Choose a sample, open a file or
  folder, or paste a link in Session › Open data; a file can also be dropped
  here." with **Show Open data**, which opens the sidebar there. The chat's
  empty state offers "Open a sample dataset and give me a first look at it"
  instead of a first look at a dataset that is not there.
- **A double click is one press.** A double click (or a held Enter) on a
  button that turns into its own Stop or Cancel undid what its first click
  started, silently: **Export** started the figure export and its second
  click, on **Cancel export**, stopped it with nothing downloaded; the same
  for **Read** / **Stop** (Genes read for), **Cancel** / **Run again**, the
  chat's **Send** / **Stop** and an agent turn's **Stop** / **Undo turn**.
  The further clicks of a double click and a held key's repeats no longer
  press a primary button or one of these (`ui/core/press-guard.js`, from the
  first paint), wherever the second click lands; other controls keep their
  double clicks.
- Every dimmed input and select says why it is off, on its tooltip and in the
  reason a refused `setting.set` gives, with what to do: `Show overlay` (the
  dataset has no vector field, none is bound to the embedding the view shows,
  naming those that are, none is selected, or it is loading), the field and
  gene pickers, Genes read for, Force LOD level, the embedding and layout
  selects, a legend category with no cells left, the figure export's grid,
  frame and text controls, the data mapping, the dataset list, the
  outlier filter, the chat's provider and model, and the dialogs' fields. A
  refused setting said "switched off by another setting" for any disabled
  control — on a dataset without a vector field no setting had.
- The **Keyboard shortcuts** panel lists `/` (**Open the chat**): the key
  opened the chat, and the accessibility page listed it, but the panel that
  page calls the authoritative list did not. A contract now fails for any key
  the page answers everywhere that the panel leaves out.
- **Actions inside a row are quiet.** Controls get quieter the deeper they
  sit: an action inside a list item's row or details is plain text with no
  border or fill, on one line, its main one in the accent color, never a
  filled button. This covers a history step's details (**Show**, **Branch
  here**, **Compare…**, **More…**, which was `⋯`; **Add note**, which was
  "+ Note") and notes, the rows of **Load…** (**Open**, **Replay here**),
  **Recently deleted** (**Restore**), the history's memory list, **My
  datasets** (**Open**, **Continue**), **Stored data** (**Remove**,
  **Show…**), **Deleted fields** (**Restore**, and **Delete for good…**,
  which was **Confirm**, as **Load…** asks "Delete for good?"), the agent
  connections and permissions, the chat's message and turn actions (**Edit**,
  **Copy**, **Try again**, **Show**, **Sent**, **Undo turn**), an agent's
  activity line (**Stop**, **Undo turn**), **Add columns** in Data mapping,
  and a Community annotation suggestion's **Edit**, **Delete**, **Merge…**,
  **Detach** and their forms. A view's × in Views is a muted glyph like a
  highlight page's, no longer a bordered pill that filled red. Typed
  characters became drawn icons with names: the branch menus' and Quick
  insights' field pickers' ⋯, the ✕ that dismisses a history banner or an
  agent's line, the ✕ that empties an entry of Recently deleted (now the
  delete icon), the × that leaves out the chat's attached selection.
- **No action looks like a link.** The analysis windows' inner actions
  (**Copy genes**, **All genes CSV**, **GMT**, **Recompute**, **Score as
  signature**, **Dot plot**, **Relax filters**, **Swap A and B**, **Counts
  CSV**, **Reset to defaults**, **All** / **None** / **Invert**, a gene set's
  **Clear**, **Remove not found**, **Add set**, **Import**, **Export**) and
  Community annotation's comment **Edit** / **Delete**, **View merged**,
  **View full evidence** and **View all comments** were underlined text, a
  second quiet look; the comment actions were faded below readable contrast
  and had no touch size. They are the quiet buttons of every other nested
  area. An underline now marks only a link.
- **A notification you are reading is no longer pushed aside by work that is
  not on screen.** Work under way shows its card only after a second, but it
  took a place among the cards at once: on a phone (two cards), a field read of
  a few milliseconds sent the best-guess notice and its **Show** to "Show 1
  earlier message" while no card was left on screen. A card now makes room only
  when it appears. A card retired to make room keeps its follow-up (**Show**,
  **Retry**) in the earlier messages, until it is taken or its owner withdraws
  it.
- A dataset picked while the sample list was still reading its files, which
  then failed to open, kept its explanation and **Try again** beside the
  selector only until the list finished: the list's arrival took them away.
- **A failure is said in a sentence, its technical words under Details.**
  Notifications and status lines quoted an error's own words: a history that
  did not hold together failed as "Loading failed · Saved state: tracker
  4f3c…: a group cuts the run of turn …", a renderer setting as "Culling and
  the tree's draw order are off …(GPU permute could not allocate its source
  texture (WebGL error 0x505))", a vector field as the overlay's own error.
  Every failure (exports, analyses, the vector field, a session or saved
  state, agents' connections, the benchmark, a control's change) now says
  what failed in plain words — a broken invariant inside Cellucid as "a
  problem inside Cellucid stopped it" — with the error itself under the
  card's **Details**; a session file that breaks the format reads "the file
  is damaged, so it was not opened; the history shown is unchanged". A
  figure export that failed showed two cards for it; it shows one. Advice to
  "try refreshing" (which loses the work open) is gone. The wording of
  failures and the ellipsis "…" are held by a contract.
- **One spelling: color.** The labels said "Color by" while step names,
  notices, tips and the actions' titles said "Colored by …" in British
  spelling (`Coloured`, `Colour cells by`). Every word a person or an agent
  reads — step names (kept in saved histories and conversations), notices,
  tips, the action catalog, the documentation, the Python and R packages'
  messages and help — now says color, colored, coloring.
- **One spelling: US English throughout.** The same mix stayed with
  `neighbour`/"neighbor", `grey`/"gray", `normalised`/"normalized",
  `centred`, `catalogue`, `labelled`, `summarise`, `recognise`,
  `organisation` and the like: the page, notices, tips, step names, the action
  catalog, the documentation and the packages' messages and help now use US
  spelling ("neighbor graph", "Gray for zero", "Normalized · log1p").
  `tests/failure-wording-contract.test.mjs` refuses the British spellings in
  the page's prose and the action catalog. Names are US too (F149): the
  parameter `normalize` of `expression.transform.set` (version 2 of its
  parameters; a step recorded with `normalise` replays as `normalize`), the
  expression values `{normalize, log1p}` saved in sessions, histories and
  `uns['cellucid']['history_applied']` (and a result's values key
  `normalized`, `normalized-log1p`), `neighbors` of `vector-field.summary`
  (parameter and result), `center` of `camera.framing`'s groups, and the
  setting kind `color` of `settings.list`. A contract holds every parameter,
  result key and value of the catalog and every field of the history model to
  US spelling.
- **"Canceled", one l.** Notices, step names and job lines said both
  (`Figure export was cancelled.`, `Pull cancelled`, "Canceled …"): every
  one now says "canceled"/"canceling" ("Move canceled.", "Differential
  expression canceled", "Canceling the new branch"). The spelling contract
  now also reads a template's words and a one-word sentence ("Canceled."),
  which it skipped, and the other British forms it missed (`travelled`,
  `totalled`, `licence`, `defence`, `judgement`, `metre`, …). The answer key
  `cancelled` of `history.branch.cancel` stays as it is (a name agents and
  scripts hold).
- **One spelling in the packages and the documentation too.** The US
  spelling held only the web app's strings: this package's documentation
  quoted "Compute the `neighbour` graph", its messages and docstrings said
  `catalogue`, `recognise` and `cancelled` (a request's
  `This session load was cancelled.`), and its documentation said
  `neighbour`, `greyed`, `centred`, `labelled` in prose and in the
  screenshots' alt texts; cellucid-r's code,
  help and vignettes did the same. One British word list
  (`tests/fixtures/spelling.json`, shared byte for byte by the web app, this
  package and cellucid-r at one pinned digest) is now read by
  `tests/test_spelling_contract.py` over every message, docstring, `#:`
  attribute doc, the documentation, the notebooks, the screenshot
  descriptions and this change log's unreleased section; by cellucid-r's
  `test-spelling.R` over its messages, help pages, vignettes, README and
  NEWS; and by the web app's contract over the app's strings and multi-line
  templates, the relay's words, every README and guide of the repository,
  and the expected text of the node tests as well as the browser specs (a
  test still expected `your organisation's administrator`). The repr of a
  waiting request says `canceled`. Words added to the list: `learnt`,
  `honour`, `tumour`, `acknowledgement`, `lustre`, `foetal`, `spelt` and the
  `un-`/`re-`/`de-` forms (`unrecognised`, `reorganised`, `deserialise`).
- **One compact count, one p-value, one first-letter case.** Two compact
  formatters wrote one count "4.0M" and "4M" (Benchmark notices, the
  performance panel, plot labels and the plugin helpers beside the dataset
  card, legends and the highlight pages); every one is now
  `formatCompactCount` of `utils/format-units.js` ("4M", "2.5M", and "1M"
  where 999,960 read "1000K"). Five rules wrote a p-value ("< 0.0001",
  `4.12e-8`, `4.123e-8`, `4.1e-8`, `0.0000412` across the DE table, the genes
  panel, correlation, composition, the hover card and the test summaries);
  every shown p-value is now `formatPValue`: two significant figures, "< 1e-300"
  where the double underflowed. A CSV keeps the number (the scatter plot's
  statistics export wrote `p_value` rounded to four figures; it writes the
  full number now). About thirty copies upper- or lower-cased a first letter,
  some lowering an acronym ("pCA computed"); `upperFirst` and `lowerFirst`
  (`utils/text-case.js`, a name in capitals kept) replace them, and the
  filter-change lines and step details write numbers by `formatBound`.
  `tests/one-implementation-contract.test.mjs` refuses a new copy of each.
- **Names said with one spelling too.** The refusal code an agent, a
  notebook or an R session reads when the person stops a dataset's or a
  saved state's loading is `canceled` (it was `cancelled`); the read log's
  line says "(canceled)". `history.branch.cancel` still answers `{cancelled}`
  (a key agents and scripts hold). Files renamed: the screenshot
  `analysis/de-canceled.png` (capture `analysis-de-canceled`), the tutorials
  `46_export_an_h5ad_subset_and_reanalyze.ipynb` and
  `12_subset_and_reanalyze.ipynb`, and the web app's tests
  `history-catalog-contract` and `action-catalog-coverage`.
- **A percent is written one way: "50%".** The app wrote both "50%" and
  `50 %` ("100 % cells · genes 6 %" in Session › Subsampling, "Hid outliers
  of x above 95 %" as a step name, "Reads about 12 % of the file" beside
  "Loading fields… 40%"), from about thirty hand-built strings that each
  rounded a share its own way ("5.0%", "5%" and "5.5%" for one share). Every
  percent is now `formatPercent` or `formatShare` (`utils/format-units.js`,
  beside the count nouns): the sign attached; a share whole from 10% (never
  100% short of the whole), two significant figures below.
  `tests/percent-format-contract.test.mjs` refuses `${n}%`, `n + '%'` and a
  number followed by " %" in the app's text. The methods paragraph of this
  package (`cellucid.history`) and cellucid-r writes a share of the cells by
  the same rule ("5.5%", where it said "6%"; "99%", never "100%" short of
  every cell).
- **A full stop ends a sentence, never a fragment.** Status lines that are
  sentences end with one ("Agents cannot connect on this page.", "Codes look
  like ABCD-EFGH.", "Reload the page to reconnect."), status words do not
  ("Repository disconnected"); "Missing dataset context." reads "No dataset
  is open.". The rule is written in the copy guidance (`assets/css/README.md`).
- **A failure is said once.** A refusal wrapped in another sentence no longer
  says its verb twice ("“Beta” could not be opened: The linked data could not
  be opened: …" reads "“Beta” could not be opened: …"): every wrapper of a
  failure's sentence goes through one helper (`utils/errors.js`
  `joinFailure`), held by the failure-wording contract.
- **A problem inside Cellucid points to Details.** The sentence for a broken
  invariant says, on a card that has Details, that its words for a bug report
  are there ("a problem inside Cellucid stopped it; Details has its words for a
  bug report."). The code's own failures (the vector field overlay's
  "VelocityOverlay could not allocate transform-feedback ownership.", the
  renderer's, the plot slots') are said as one, their words under Details.
- **A damaged file says it is damaged.** A session file or a My datasets
  library file whose framing breaks (truncated, trailing bytes, a manifest out
  of the format) read raw "Invalid session file (manifest JSON parse
  failed)." or "This library file is damaged: chunk 1 must be exactly …"; it
  now says the file is damaged, with what breaks it under Details (the
  library's Import and Export status lines on hover).
- The chat's **Effort** menu says **Extra high** and **Maximum** (it showed
  the provider's values, "Xhigh" and "Max"), and so does its status line.
- A conversation saved during a turn ends on a **saved** note, which offers no
  **Try again**: in a session file someone received, it would have sent the
  sender's message to the reader's agent.
- A conversation opened from a session file or a saved state is received: no
  note of its turns offers **Try again** (the sender's failed or stopped turns
  did, and acted again on what the history it came with holds), and a question
  its agent asked offers no choices (it was put to the sender; you can still
  write). Turns you send after opening it are yours. The saved conversation
  names the turns it brought (`received`).
- Counts and marks: "…and 3 more categories" in the Gene signature statistics
  (it said "1 more categories"), counts grouped there and in "Loaded 12
  fields"; the edge limit reads "250K" like every other count; the
  performance panel's details lose their colons and Title Case ("GPU memory",
  "Frame time"); notification titles are sentence case ("Generating data");
  in-progress text ends in "…" everywhere ("Preparing projectiles…", which
  was "Preparing projectile collision system...").

#### Performance tools

- **Analyze performance finishes again.** It stayed on "Measuring
  performance…" for good, the button disabled and the viewer paused, after
  the renderer began to require the sprite cap the analysis did not pass; a frame that fails now ends the analysis with a notice saying why
  and resumes the viewer (it used to be an `alert`).
- Its **Problems found** no longer names "Vertex Processing: LOD scaling
  factor: 9921.55" on every run: the level-of-detail sweep read level 0 (the
  smallest sample) as the most detailed, so a GPU with room to spare scored
  as the most vertex-bound. Both sweeps now report how many times longer a
  frame takes with the largest points, or the most points, and **What to do**
  gives frame-time advice only below 55 frames a second.

#### Development

- **Each package's CI runs the cross-repository tests.** cellucid-python's
  workflow checks out cellucid-r and the web app's code beside it, and
  cellucid-r's fetches the web app and cellucid-python at the commit it
  pins, both with Node and `CELLUCID_SIBLINGS_REQUIRED`: a sibling or Node
  missing fails those tests instead of skipping them (a parity test that
  read cellucid-r failed every job without it). Every sibling read goes
  through one helper (`tests/_web_root.py`, `helper-repository.R`).
- **The web contracts pass on Windows.** Seven read a file URL's `pathname`
  as a path or compared `\`-separated names with `/`; the contracts that
  walk the sources use `tests/helpers/source-files.mjs`, and a contract keeps
  the `pathname` form out of every test and script.
- **One implementation of hex colors and of the decoded-piece cache.** The
  `#rrggbb` parse and format were written out a dozen times with as many
  strictnesses (the legend, figure export's and the analysis' `color-utils`,
  the colors and field-overlays domains, the field actions, the cell exports,
  the palettes, the dataset reads, the dot plot, the Plot theme, the plots'
  `hexToRgba`, the connectivity color input, eleven copies of the pattern);
  they are one module, `assets/js/utils/hex-color.js`, held by
  `tests/hex-color-contract.test.mjs`. The HDF5 file's `ByteCache` and the
  Zarr store's `ByteLru` were the same class twice: one `ByteCache`
  (`formats/io/byte-cache.js`).
- **One part of the page that fails no longer stops the others hearing a
  change.** 44 sets of listeners in 31 modules had five ways of meeting a
  listener that throws: the rest went untold (`DataState`'s events, the camera
  path, the memory budgets, a worker's failure, a sign-in revoked), the error
  was swallowed, logged, or collected and thrown back at the change, so an
  opening whose dataset was shown was said to have failed when one display of
  it threw. They are one listener set (`assets/js/utils/listener-set.js`):
  every listener hears every change, one that throws is reported once on the
  console, naming what it listened to, and a listener removed while a change
  is told is not called after (held by `tests/one-implementation-contract.test.mjs`).
- **One way bytes become text, one SHA-256, one string hash.** Hex digests,
  random hex states, base64 and base64url were written out in eighteen
  modules (sign-in states, the community-annotation Worker's PKCE, the tab
  secret, OneDrive share ids, pictures sent to a chat, kept stores, the cell
  record of a session, GitHub contents, Zarr fill values); they are one
  module, `assets/js/utils/byte-text.js`. Three modules called
  `crypto.subtle.digest` themselves, which a plain-HTTP origin (the Python
  package served on a LAN address) does not have: **Save AnnData with this
  state…** failed there while digesting the state it saves, as did a kept
  folder's identity and a kept store's directory name. Every digest now goes
  through `utils/sha256.js` (WebCrypto, else the same digest in script),
  which moved with `opaque-id.js` from `app/utils/` to `utils/`, where the
  readers may import it. FNV-1a over text was written four times (the string
  table, the retry spread, the figure seed, the pull-request branch name):
  one `utils/fnv1a.js`, same values. Sweeps: `tests/byte-text-contract.test.mjs`,
  `tests/one-implementation-contract.test.mjs`.
- **Download names, one rule.** Every download leaves through `downloadBlob`,
  which now makes its name safe for every file system (`safeFileName`): the
  history report, a session file and an analysis table named after a dataset
  such as `a/b` no longer carry a `/`. A figure keeps its dataset's name in
  every script (an accented or CJK name was dropped to `untitled`). The time
  stamp in a file name is one, local time `YYYYMMDD_HHMMSS` (`fileTimestamp`):
  the session saved after a graphics reset, the My datasets library and the
  exported chat wrote UTC `YYYY-MM-DDTHH-MM-SS`.
- **Smaller copies made one.** "This computer" in the address box and the
  Share dialog is `utils/loopback.js`'s (it missed `0.0.0.0`, `*.localhost`
  and `127.0.0.2`); HTML escaping is `escapeHtml` (the report and the embed
  code had their own); "1 repo available", "1 note, 1 bookmark" and "1
  category label" agree (the count-noun contract knows those nouns now); the
  level-of-detail and clustering refusals give sizes as every other message
  does ("3.8 MB", never "4 MiB"); a superseded gene-sample read is said as
  superseded, not as your Stop, and an abort's reason is always an error;
  the next-task scheduler is one (`onNextMacrotask`; the history's own made a
  channel per call, which WebKit can collect before it fires); `Retry-After`
  is read once (`hosts.js`: an empty header no longer means "wait 0 ms").
- **No cache keeps a count or an age limit; every one is bounded by memory.**
  The analysis page-data cache kept at most 100 entries beside its share of
  the memory budget, so a page × variable sweep of small results read again
  while memory was free; its share is now its one bound (no entry cap, no
  age), and the cleanup pass that only pruned expired entries is gone.
- A remote file's access (Google Drive, OneDrive, Figshare: the provider's
  closure holding its item and the token getter) is registered only once the
  file is described (a refused file registered one) and let go with the
  source that reads it: a dataset replaced or closed, an opening not adopted,
  a link whose opening failed before a reader took it. They were kept for the
  page's life, one per opening.
- A DOM id minted on the page holds a reserved separator (`<prefix>--<n>`),
  which no hand-written id holds (a contract test reads `index.html` and the
  sources): it is recognized by its form, so the set of every id minted (kept
  for the page's life) is gone, and an id made from a minted one
  (`<minted>-trigger`) is left out of a session as the minted one is.
- The page no longer deletes the IndexedDB database `cellucid-highlight-drafts`
  that 0.9.6 left behind: a one-time cleanup of an earlier version is a
  migration step, which the project does not keep.

- **One debounce, one test of a signal, one wait in words.** Twenty
  debounces written by hand (`clearTimeout(t); t = setTimeout(…)`: the
  history search, figure export's preview, the wheel's rest, the history
  graph's hint, Quick insights, the camera path's bar, idle workers, kept
  parts' journal, …) are `utils/debounce.js` `debounce` (with `cancel`,
  `flush`, `pending`). Eleven validators of an `AbortSignal` (several took any
  object with an `aborted` boolean for one) are `isAbortSignal`,
  `requireSignal` and `requireSignalOrNull` in `utils/abort.js`, by the
  platform's brand check. The annotation lock's "Lease remaining: ~12 s" and
  a gene sample's "a gene in ~20 s here" are written by `formatRemaining`.
  `tests/one-implementation-contract.test.mjs` finds a new copy of each.
- `npm run serve` runs the app locally with `Cache-Control: no-store`, so a
  reload runs the module just edited; the docs no longer recommend
  `python -m http.server`, under which a reload kept running the old code.
  Its second port serves every file of a dataset (a `.loom`, an `.rds`, a
  MatrixMarket folder) to open by link.
- The documentation builds with sphinx-design 0.8 (`docs` extra:
  `sphinx-design>=0.8,<0.9`): the home page's cards use the `card-header` and
  `card-footer` directives in place of the `^^^`/`+++` separators 0.8
  deprecates, which failed the `-W` build.

### Python

- **The viewer token stays out of errors and logs.** A browser that could not
  be opened raised `Could not open the browser for <the address with its
  viewerToken>`, which lands in tracebacks, CI logs and bug reports; and the
  server's request log (at debug level) wrote the viewer page's request line,
  token and all. Both now say the address without the token
  (`public_viewer_url()`, one helper beside the banner, which the debug report
  uses too); the address with it is said only where you are told to open it
  (the banner, `display()`, `cellucid mcp --open` on stderr). A contract test
  classifies every sentence that writes a viewer's address.
- **A store whose declared id or name cannot be read is refused with that
  sentence, unless you name that field.** `serve_anndata()`, `cellucid serve`
  and the viewers read `uns['cellucid']['id']` and `['name']` even when
  `dataset_id` / `dataset_name` were given (the declared id applies a session
  back to the object read from the store), and a damaged one failed the view
  with the raw decoder error. It now says `The store f.zarr holds
  uns/cellucid/id and it cannot be read (…). Write the store again, or give
  dataset_id= yourself.`, and a field you name is not needed from the store,
  as in cellucid-r.
- **`cellucid serve` says its Viewer URL at once into a file or a pipe.**
  `cellucid serve data.h5ad --no-browser > serve.out` (a process manager,
  `nohup`, an IDE's run panel) wrote its steps, but the Local URL / Viewer
  URL box stayed in Python's output buffer while the server ran, so the
  viewer could not be opened (`PYTHONUNBUFFERED=1` showed it). Every line the
  command, `serve_anndata()`, the viewer and `prepare()` print now reaches
  its reader at once (`console_print` flushes each one), and `cellucid mcp
  --open` writes its address to stderr the same way.
- Each way of opening a viewer names itself to it (`/_cellucid/info`'s
  `runtime`): a notebook kernel's `jupyter`, a script's or session's
  `python` (`serve_anndata()`, `AnnDataViewer`), `cellucid serve`'s
  `cellucid-serve` (it was `python`) and `cellucid mcp --open`'s
  `cellucid-mcp` (it was `python` or `jupyter`), which the viewer's
  sentences name; the agents folder's records carry the same.
- The viewer's `field_restore()` and `field_purge()` actions take `which`:
  when two deleted fields have one name, which of them, in the order the web
  app's **Deleted fields** lists them (named by the name alone, the request
  was refused with no way to name either).
- The viewer channel never allocates a body at the length a frame names
  (R9Z-3): a selection, a result, an action's outcome or an agent's result
  grows as its bytes arrive (a frame header naming 256 MiB, of which a few
  bytes came, zeroed 256 MiB at once). A selection longer than the dataset's
  cells, a cell beyond them (a selection's, a hover's, a click's) or a result
  for a dataset of another size closes the channel (`1009`, `1008`), as
  cellucid-r refuses them; before, the cells were delivered and a
  selection's obs names failed with an `IndexError` in the notebook.
- A file written again or replaced on disk while a server shows it (an
  `.h5ad`, `.h5mu`, `.loom`, 10x `.h5` or `.zip` served in priority order) is
  refused with HTTP 409 and why ("… changed on disk since the server opened
  it: serve it again …"), which the viewer says as it is and then reads
  nothing more: every chunk was answered 500 with an HDF5 error ("wrong
  B-tree signature"), said as the server failing.
- Offline, a viewer start serves the cached web app, as the R package does:
  when the web source cannot be reached, a cached build whose every file
  matches its checksums is that exact build, and the start says "Could not
  reach …; serving the cached web app (build …)." (also when quiet). With no
  cache, or one that does not check out, it stops saying which and how to get
  the app; it stopped with the bare network error and never served a cache.
  Both packages print the same words (`web_cache.WEB_APP_*`, pinned against R's
  by the parity tests); `WebCachePrefetchSummary.unreachable` says why the
  source was not read.
- A viewer start no longer downloads the web app it already has: the source's
  inventory is read at every start, and when the cache holds that very build,
  every file verified, nothing else is fetched (each `cellucid serve`,
  `serve_anndata()` and `AnnDataViewer` start fetched all 1,209 files again,
  over six minutes from www.cellucid.com). A new build copies the files it
  did not change from the cache and downloads the others six at once (the
  first download took 6 min 49 s one file at a time, now under a minute).
  The R package does the same.
- In a notebook, a session or an action's answer collected in a later cell
  before the viewer sent it says so, and why (the viewer is not shown, has
  not connected, or is still working on it): it said "ask in one cell and
  collect it in the next", which is what had been done.
- `show_anndata()`, `serve_anndata()` and `AnnDataViewer` show an AnnData
  holding a record array in `uns`, as every Scanpy object does after
  `sc.tl.rank_genes_groups`: they failed with "Field elements must be 2- or
  3-tuples" (the served store read a structured Zarr dtype as a type string).
- An obs column of integers holding identifiers past 2^53 is read as the
  viewer reads it (F116): `apply_cellucid_session_to_anndata()` takes its
  categories as their exact decimal texts and writes a renamed label back as
  text, `push` sends it as a categorical of those texts, and a MuData
  modality's such column is placed on the cells as text, never as rounded
  doubles. A boolean column's categories are `false` and `true` when a
  session applies, as the viewer shows them: a label renamed in the viewer is
  written (the column becomes text) instead of "holds numbers, not
  categories".
- `prepare()` refuses an element name longer than 255 bytes in UTF-8 (the
  longest file name every platform takes) with the rule's sentence instead
  of zarr's `OSError: File name too long`, and refuses two names that differ
  only in Unicode form (`é` as one code point or as `e` and an accent), as
  the web app's writers do.
- `serve_anndata()`, `show_anndata()`, `AnnDataServer`, `AnnDataViewer`, `Widget`
  and `cellucid serve` name the gene-name column and the opening embeddings as
  `prepare()` and cellucid-r do: `var_gene_id_column=` (CLI
  `--var-gene-id-column`) and `embedding_defaults=`. The former
  `gene_id_column=`, `embedding_keys=` and `--gene-id-column` are gone, with
  no second name.
- `prepare(total_counts_over=)` (new) says what the column `total_counts`
  names covers: `"all genes"` (the default) or `"source panel"`, every gene of
  the source's gene panel -- a targeted or integration panel, more genes than
  are stored but not the whole transcriptome. The Cellucid-ready AnnData
  specification's `totalCounts` `over` gains `"source panel"`, which the web
  app reads and says ("Cell totals from obs['total_counts'] (source panel)"),
  and a file prepared again keeps it. With a MuData, one value or a mapping
  by modality. cellucid-r's `cellucid_prepare(total_counts_over =)` is the
  same. The four atlas samples' totals, over their sources' 8,192-gene
  integration panels, are recorded so from their next build (they said
  `"all genes"`).
- `prepare(data=adata)` writes the expression matrix the viewer opens the
  AnnData with (`cellucid.data_mapping.default_expression`, the web app's
  `defaultExpression`, pinned by the shared conformance fixture): `raw.X` with
  `raw.var` when `raw` holds more genes than `X` or `X` is scaled -- a Scanpy
  object after `adata.raw = adata`, the highly variable genes and
  `sc.pp.scale` -- else `X`; without `X`, `raw`, else a layer. It wrote `X`
  whatever it held, so such an object's prepared file colored genes by the
  scaled values of the highly variable genes where the object itself (opened,
  served or shown) showed the log-normalized values of every gene; without
  `X` it wrote no genes. A message says which matrix was written and why, and
  how to write another; a `gene_expression=` given in its place goes with
  `var`. The note for a `raw` left out no longer suggests
  `adata.X = adata.raw.X`, which fails when `raw` holds more genes.
  A MuData's modalities follow the same rule: a Scanpy-processed `rna`
  modality is written from its `raw` (the root `var` and that modality's
  `varmap` then hold `raw`'s genes), and a modality with no `X` but a `raw`
  or a layer is written from it instead of refused.
- `prepare()` writes a column's colors (`obs_colors=`, or Scanpy's
  `uns['<column>_colors']` of `data`) only when each category's is hex
  (`#rrggbb`, `#rgb`, `#rrggbbaa`), the forms the viewer reads: a palette
  holding a named color (`"red"`) was written and the viewer opened the
  column in its own colors without a word. Given as `obs_colors=`, such a
  palette is refused; held by `data`, it is left out with a message.
  cellucid-r's `cellucid_prepare()` checks the same.
- `cellucid serve`, `serve_anndata()` and `show_anndata()` refuse a file
  Cellucid does not open in the web app's words, as cellucid-r does: "<path>
  is not a file Cellucid opens (an .h5ad file, a .loom file, …, an .rds file,
  … or a .zarr.zip archive)." It listed one hand-picked ending per format
  (`.rds.gz`, `.rdata.gz`, and the refused `.h5seurat`), and `cellucid serve`
  wrapped it in "Unable to detect format for: <path>.", naming the path twice.
  A missing path is "Nothing exists at <path>." (it was "Path not found:").
- A server that cannot bind its address says why in one sentence, the one
  `cellucid serve` says, with each setting spelled as the caller gives it
  (`host='127.0.0.1'`, `port=9000` in Python; `--host`, `--port` on the
  command line): a `host=` that is no address of this machine (it was the
  resolver's `gaierror` or "Can't assign requested address"), an explicit
  `port=` already in use (it was "[Errno 48] Address already in use"), a port
  the system refuses. It is an `OSError` (`ServerBindError`, its `errno` the
  system's). cellucid-r says the same.
- Descriptive text in `uns['cellucid']` -- `name`, `description`, the
  citation's text, the source's name -- is read and written with white space
  at either end removed: a file whose description ended in a space was
  refused whole by the web app ("must not start or end with white space").
  `prepare()`, `serve_anndata()`/`show_anndata()` (`dataset_name`) and the web
  app's browser copy write it trimmed; identifiers (`id`, keys, fields) keep
  the strict rule. One rule in the specification, the web reader and both
  packages; cellucid-r's `cellucid_prepare()` and `cellucid_view()` are the
  same.
- A refused `total_counts` column quotes the offending value as a number
  (`-2.0`), not as NumPy's `np.float64(-2.0)`.
- A change the viewer applied without a step of its history now says why.
  The viewer's `action_result` carries `unrecorded` (`"unchanged"`: it was so
  already; `"retracted"`: it put back the step just before it, which was
  taken out; `"off"`: the History settings leave this kind of change out;
  `"paused"`; `"elsewhere"`: another tracker was shown meanwhile;
  `"no-history"`: the dataset's history has not started), kept as
  `ActionRequest.unrecorded`; `request.result()` returns the outcome and
  warns `cellucid.jupyter.ActionNotRecordedWarning` (a `UserWarning`, new),
  as in "viewer.color_by(...): it changed nothing (it was so already), so no
  step was recorded." A read records no step by nature and warns nothing.
  Before, `request.step` was `None` with no reason, whether the change
  happened or not. The viewer protocol fixture (`action_result` header,
  `action.unrecorded` sentences) is re-pinned in the three repositories;
  cellucid-r's `request$unrecorded` and `cellucid_action_not_recorded`
  warning are the same.
- `viewer.view_rename(view="View 2", name="Early")` (new, generated from the
  viewer's `view.rename`) names a view or, with `name=None`, takes its name away;
  `viewer.view_keep(name=...)` keeps and names a view in one step;
  `viewer.views_state()` and `viewer.view_picture()` give each view's name
  (`label`).
- `viewer.annotation_publish(votes=[...], suggestions=[...])` (new, generated
  from the viewer's `annotation.publish`) votes and suggests community
  annotation labels together, published as one commit of the person's file
  and asked of them once; `viewer.annotation_vote()` and
  `viewer.annotation_suggest()` are removed with the viewer's actions.
- `apply_cellucid_session_to_anndata()` writes every categorical user field
  with the colors it is shown in, where Scanpy keeps a column's colors,
  `uns['<key>_colors']` (one `#rrggbb` a category, every category), so
  `sc.pl.*` and the viewer show the field alike: in the view in focus at the
  applied step, the color the person or an agent gave a category (the
  session's `colors` value; a field grouped from clusters and recolored gray,
  gold and teal is written so), else the color it opens with -- its `palette`
  (read and checked: a field made from highlight pages opens in its pages'
  colors), else Cellucid's own category colors by position (the web app's
  `getCategoryColor`, ported and pinned by the shared fixture
  `category_default_colors.json` in the three repositories). The plan's
  Colors row stays `not-applied` -- every other color is how the data was
  shown -- and says so when a field's colors went with its column. A
  column's colors follow it in every write of an apply and a pull: a column
  without categories (a highlight group, a continuous field) has none (stale
  colors under its name, which Scanpy and the viewer would have used for the
  wrong categories, are removed); a renamed column carries its colors; merged
  categories keep their first category's color (colors that did not cover
  the categories are removed); a text column relabeled loses them; a column
  a viewer pull removes takes them along. `pull()` of a session file writes a
  field's colors the same way, under the column name it writes. Before, no
  colors were written but a highlight-page field's palette (a grouped or
  merged field's chosen colors never reached the object), and a rename or a
  merge left colors that described another column. The web app's **Save
  AnnData with this state…** and cellucid-r write the same.
- A column that came with the dataset (`clusters`, `cell_type`) recolored in
  the view in focus is written with the colors it is shown in too: one
  `column-colors` row of the apply plan per column (a new row kind of the
  history model), the categories recolored, the rest in the colors the
  column opens with (its own `uns['<column>_colors']` as the viewer draws
  them, else Cellucid's), under the column's name once the apply's renames
  are made; `unchanged` when the column holds them already; `not-applied`,
  saying why, for a column that is not categorical, a MuData modality's
  column, or one with none of the categories recolored. Before, the colors a
  person gave a dataset column never reached the object. And `pull()` from a
  live viewer writes a categorical field's colors: each categorical result
  of the viewer protocol carries the colors the view in focus shows it in
  (`colors`, one `#rrggbb` a category), checked when it arrives. Before, a
  viewer pull wrote no colors. cellucid-r and the web app's copy do the same.
- A viewer given another `dataset_id=` than the one its data declares
  (`AnnDataViewer(adata, dataset_id=...)`, `show_anndata`, a widget) takes
  its own session back: `viewer.apply_session_to_anndata(adata, request)`
  refused it ("This session was saved on the dataset 'pancreas-notebook',
  and this object is 'pancreas'"), since the session is saved under the id
  served and the object declares its own. The viewer now names the id it
  served the data under for an object declaring the id its data declared;
  `apply_cellucid_session_to_anndata()`, `CellucidSession.apply()` and
  `pull()` of a session file take `dataset_id=`, the id the object was shown
  under, and the refusal says to name it. A session of another dataset is
  refused as before, and the cell count, data mapping and cell order are
  checked as always.
- `viewer.category_delete()` moves the cells into `'Unassigned'` by default
  (the one spelling of cells without a label, as Create categorical's; it was
  `'unassigned'`).
- `pull()` of a field's obs-column CSV keeps the field's categories and their
  order: the web app writes a third row, `categories`, holding each
  categorical column's categories in the field's order as a JSON array of
  strings, and `pull()` builds the column from it. They came back in the
  order the rows first showed them, so Pancreas `clusters` lost its
  trajectory order and a category no cell held was lost. A label the row does
  not list, a categorical column without a list and a CSV without the row (the
  two-row layout of earlier exports) are refused; export the field again.
- `viewer.history_saved_open(saved=...)` (new, generated from the viewer's
  `history.saved.open`): the viewer refuses it, as `dataset_open`, since a
  notebook's viewer shows the dataset its notebook serves.
- `viewer.figure_export_get()` reads the Figure export settings, by the
  names `viewer.figure_export_set()` takes.
- `viewer.highlight_category()` takes `categories`, one or more
  labels of one field (was `category`, one): each becomes a group of the page
  in one step of the viewer's history, as the viewer's `highlight.category`
  version 2 does.
- `prepare()` records each cell's total count (its library size) when `X`
  holds counts: `obs['total_counts']`, the sum of each row of `X` as written,
  with `uns['cellucid']['totalCounts'] = {"field": "total_counts", "over":
  "stored genes"}`, so the web app normalizes a gene from one column instead
  of reading every gene once (Suo: the whole 74 MB store). `total_counts=`
  names the caller's own column, computed over every gene before any was
  selected (`sc.pp.calculate_qc_metrics` before subsetting), recorded `"over":
  "all genes"` and refused when it is below a cell's stored genes' sum (a log
  of the totals) or `X` holds no counts; with a MuData, one column per
  modality (`{"rna": "total_counts"}`), and each counts modality's own sums in
  `mod/<name>/obs/total_counts`. A prepared file prepared again keeps a named
  column named. The Cellucid-ready AnnData specification's new `totalCounts`
  section is the contract, which cellucid-r and the web app follow.
- `prepare()` keeps Scanpy's colors of the columns it writes: each categorical
  `obs` column's `uns['<column>_colors']` is written with it, so a
  Scanpy-colored AnnData (or one a session was applied to) opens in its own
  colors instead of Cellucid's. A MuData's are kept where they were: the
  root's columns' in `mdata.uns`, a modality's own in `mod/<name>/uns` (the
  viewer reads `<name>:<column>`'s there). `obs_colors=` (new) gives or
  replaces them (`{column: [color, ...]}`, a modality's column as
  `<name>:<column>`; `{}` writes none), so the same components passed by hand
  write the same file again. Colors that cannot color the written column (not
  a list of text, fewer than its categories, a column not written as a
  categorical) are left out with a message, or refused when given. Before,
  every column's colors were dropped. cellucid-r's `cellucid_prepare()` and
  `cellucid_view()` of an object keep them by the same rule; the
  Cellucid-ready AnnData specification states it.
- `apply_cellucid_session_to_anndata()` records the expression values the
  genes were shown in at the step applied: each record of
  `uns['cellucid']['history_applied']` (`cellucid.history_applied()`, the
  history model's `Provenance`) has `expression`, `{"normalize": bool,
  "log1p": bool}` as the step's **Expression values** held it (None for a
  session without them), so a notebook knows what the person looked at.
  Recorded, never applied: `X` and the layers keep their values. cellucid-r's
  `cellucid_apply_session()` and the web app's **Save AnnData with this
  state…** write the same field (history model, shared golden provenance).
- `prepare()` stores every chunk of a Zarr store, also one that holds only
  the fill value (the priority seed of 0, an all-zero column), as
  `cellucid_prepare()` and the served store do: a chunk left out was answered
  404 by a web host at every opening, reported by the browser as an error.
  The Cellucid-ready AnnData specification states the rule.
- The action catalog shipped with the package (`cellucid/action_catalog.json`)
  is the web app's of this release: every action declares what it captures
  (`captures`; `camera-path.keyframe.add` captures `camera`), and the
  generated action reference says so.
- `viewer.camera_path_keyframe_add(hold=...)` and
  `viewer.camera_path_hold_set(seconds=..., keyframe=...)` give a camera path
  keyframe a hold (seconds the camera rests there), generated from the
  catalog as every action is; `cellucid mcp` lists `camera_path_hold_set`.
- `viewer.highlight_cells(cells=...)` without `label` makes a group named
  "Highlight" (the viewer's list shows its cell count), not "Highlight (N
  cells)"; the method's docstring and the `cellucid mcp` tool say so.
- `viewer.field_group(field=..., groups=[{"label": ..., "categories": [...]}],
  rest=..., name=...)` makes a new categorical field of groups of a field's
  categories in one step; `viewer.dataset_open(...)` takes `count` or
  `fraction` (and `method`, `seed`, `field`, `allocation`, `embedding`) to
  open a subsample at once, and `viewer.cells_open(...)` takes `fraction`
  beside `count`; generated from the catalog, listed by `cellucid mcp`.
- `viewer.camera_path_stop()` stops a path `viewer.camera_path_play()` is
  still playing, at once, and returns `{"stopped": true}` (`false` when
  nothing played; it was refused). It and `viewer.analysis_cancel(...)` are
  documented as stopping work under way (the catalog's new effect
  `control`); `cellucid mcp` describes them so.
- A hook for an event the viewer never sends is refused:
  `viewer.register_hook("selections", ...)` (and `@hooks.on`,
  `unregister_hook`, `clear_hooks`) was taken and never fired, and
  `viewer.wait_for_event("selections")` waited out its timeout. The events a
  hook takes are `cellucid.jupyter.HOOK_EVENTS` (every event of the viewer
  channel but an agent app's, and `message`), pinned to the shared protocol
  fixture, the names R's `cellucid_on()` takes.
- `expression_source="layer:<name>"` (`serve_anndata()`, `show_anndata()`,
  `cellucid serve --expression`) refuses a layer name holding `/`, which the
  web app's `defaults.expression` (`layers/<name>`) cannot name; R's
  `cellucid_serve()` refuses it alike.
- **The local servers hold nothing for a connection that never says who it
  is.** A viewer channel socket that says no `hello` within 60 s is closed
  (1008, "No hello frame was received"), as cellucid-r closes it
  (`HELLO_SECONDS`); a binary first frame is refused from its header, unread
  (it was read to its end, however long it said it was); and a connection to
  the local agents endpoint that sends nothing for 120 s is closed. Each held
  a server thread of the notebook's process for as long as any program on
  the computer kept it open, with no token needed. A binary event's header
  is read as its bytes arrive (over a WebSocket or a notebook widget), never
  allocated at the length its first four bytes name: up to 4 GiB, zeroed,
  before any of it came.
- **JSON nested past what Python follows is refused in a sentence**: a
  `hello` or event of the viewer channel (1008, "… nested more deeply than
  this package reads"), a command to the agents endpoint (400), and a session
  file (`CellucidSession`, `apply_cellucid_session_to_anndata()`, a session
  bundle): "This session file holds a value nested more deeply than this
  package reads …", which cellucid-r says too. Each raised the interpreter's
  `RecursionError`: a traceback on the kernel's stderr per such socket, and a
  session bundle's temporary file (the person's session) left behind with
  the wait running to its timeout. Whatever else ends a bundle's read now
  removes its file too.
- `debug_connection()` reports `viewer_url` without the viewer token (the
  report is what you paste when asking for help; the page's own debug
  snapshot left it out already). `viewer.viewer_url` is the address to open.
- `cellucid mcp` with a relative `CELLUCID_AGENTS_DIR` says so in
  `cellucid_status` and to each call ("CELLUCID_AGENTS_DIR must be an absolute
  path"), as cellucid-r's `cellucid_mcp()` does; it ended at start with a
  `ValueError`.

## [0.9.7] - 2026-11-07

Version 0.9.7 is a patch release on 0.9.6. Files prepared and dataset
catalogs written by 0.9.6 are read as before. **Session files are format 4**:
a session file saved by 0.9.6 or earlier (format 3 or before) is refused,
saying which format it is. The viewer protocol is schema 5: use the Python and
R packages of this release with the web app of this release. **Community
annotation repositories take the version-3 schemas** of
`cellucid-annotation`, and its Worker is contract 3 (a self-hosted Worker must
be redeployed). The main changes are:

### Highlights

- **History v2.** Every change a person, an agent or a notebook makes is one
  step of one history, shown as a branch graph in the new **State** section:
  branches and trackers, recording on and off, groups and squashed steps,
  notes and bookmarks, search, Recently deleted, comparing and replaying a
  history on another dataset.
- **Share and load whole analyses.** **State › Session › Share…** writes the
  history, the current state only or the steps only, with exactly the parts
  chosen; **Load…** adds a file's trackers beside yours and can be undone.
  Python and R read every file the web app writes.
- **Save AnnData with this state.** The History **⋯** menu writes a new
  Cellucid-ready AnnData `.zarr.zip` of the open dataset with your labels,
  fields, highlights and gene sets written in, by the same rules as applying
  a session in Python and R.
- **Agents.** A chat bar with your own Anthropic or OpenAI key, a connector
  for Claude and ChatGPT, and a local MCP server for agent apps on your
  computer (`cellucid mcp` in Python, `cellucid_mcp()` in R). All three use the
  same tools and permissions, show confirm cards, record their steps as
  turns, and let you take a turn back with one Undo turn.
- **What agents are sent stays yours.** A one-time notice says what an agent
  reads. The per-dataset switch **Agents see names and summaries only** keeps
  single-cell values, extremes and pictures out of what any agent is sent.
  Text written by others reaches an agent marked as data, never as
  instructions.
- **Every action, from Python and R.** The web app's action catalogue (177
  actions) is shipped with both packages as one generated function or method
  per action, each one step of the viewer's history.

### Web app

#### The history

- **One history of every change.** Every change is a step: colouring,
  filters, labels, highlights, views, the camera (one step per gesture),
  sidebar settings, panels, windows and theme, analyses (their view, plot,
  size and window), figure exports, camera-path playback and **Cells to
  open**. Going back shows labels, fields and analysis results as they were,
  without computing anything again. **Undo** / **Redo**, `Ctrl/⌘+Z` and
  `Ctrl/⌘+Shift+Z` (`Ctrl+Y` off a Mac) also work with a slider or checkbox
  focused. Text fields keep their own undo. Each undo or redo is announced to
  screen readers.
- **The branch graph.** Lanes, fork stubs, agent and notebook turns as one
  row, groups as one node, folded runs ("Camera ×6"), a paused marker, and ★
  bookmark, ✎ note and ⚠ stale marks. A click selects a row, a double-click
  (or a double tap) shows its state, and **›** opens its details (action,
  parameters, the parts of the state it changed, who took it). Right-click
  opens the row menu at the pointer. Shift+click and Shift+↑/↓ select a range,
  and `Ctrl/⌘+G` groups it. Steps are numbered in one way everywhere (the
  opening is not counted: "4 steps · 1 undone"). The graph draws only the rows
  in view: with 10,000 steps it opens in about 30 ms, and a frame takes about
  8 ms while scrolling.
- **Branches and trackers.** Acting after an undo starts a branch; the undone
  steps stay on theirs. The branch picker shows the branches, and **Branch
  map…** shows every branch as a lane graph and a table: find, archive or
  delete several at once. A branch can be renamed, archived, locked or made
  main, compared with another side by side as kept views, or checked with
  **Verify** (replay it and compare every part of the state). **Start another
  tracker…** keeps an independent history of the same dataset. A step's
  **Compare…** lists what differs between it and another step, part by part.
- **Recording on and off.** The **Recording** switch pauses the history: your
  changes still apply, and resuming adds them as one squashed step. Going back
  while paused asks **Keep / Discard / Cancel** (History settings can turn the
  question off). An agent asks before pausing ("Pausing the history"). What
  changed while paused is kept in this browser, so **Continue** after a closed
  tab or a crash brings it back.
- **Groups and squash.** **Group** puts a named folder over consecutive steps
  of one branch ("QC", "Clustering"), with notes of its own; the row and
  selection menus offer **Add to group "X"**. **Merge into one step…**
  (squash) replaces consecutive steps with one step holding their net change
  (refused across a fork), and its question offers **Group instead**.
  **Clear history, keep the state…** starts afresh from what is shown. A
  squashed step says why it was made and who acted; the notes and bookmarks
  of the steps merged move to it.
- **Notes, bookmarks and search.** Steps, branches, groups and turns take
  notes (a small Markdown subset, shown as text; an agent's turn gets a
  one-line summary), and a step takes a bookmark. **Search** (`/`) finds
  steps, notes and bookmarks; **Filter** narrows them by who acted (an agent,
  Python, R, a replay), by kind, by the part of the state changed, by day, or
  to steps with a note, bookmarked, or stale or flagged. A group's row says
  where the step shown is among its steps ("at 4 of 7"), and a step's details
  offer **Show this step's view** (its cameras, even while camera moves are
  not recorded; the tour flies to them too).
- **Recently deleted.** Deleting a step re-applies the later steps without it;
  a step that needed it is kept and flagged with the reason. Deleted steps,
  branches, trackers, squashes and loads go to **Recently deleted…**, where
  each can be restored. Every history edit shows a notice with **Undo**, and
  each new notice replaces the last.
- **Replay on another dataset.** **Replay a saved history…** applies a saved
  history's steps to the dataset shown, finding fields, categories and genes
  by name. The **Replay report** says what became of each step (applied,
  unchanged, unmatched, skipped or stopped) and what was missing, and opens by
  itself when something did not apply. A field, category or gene not found
  offers **Use instead ▾**, and **Replay again with these names** replays the
  same steps from the same source in a new tracker that keeps the names
  chosen. Genes and a category highlight's category are mapped too.
- **History settings** (⚙): **Camera moves**, **Panels, windows and theme**,
  **Record history for new datasets**, **Ask before going back while paused**,
  **Your name in shared histories**, **Autosave in this browser** and
  **Pictures of the steps chosen**. Each is kept in this browser.
- **Kept in this browser as you work.** Steps are written to IndexedDB
  `cellucid-history` as they are made, and each write holds only what changed.
  Reopening a dataset offers **Continue where you left off?** (Continue /
  Start fresh). Two tabs never mix their steps: a saved state continued in
  another tab opens as a copy there. When storage fills, autosave pauses and
  says so; removing saved states anywhere resumes it. The highlight drafts
  0.9.6 kept in this browser are gone: saved states replace them, and their
  IndexedDB database `cellucid-highlight-drafts` is deleted when the page
  starts.
- **Memory warnings, not walls.** Once the history holds more than the device
  has room for, the status line shows `⚠ 1.2 GB` and says once which parts
  take the most; the **History memory** dialog lists every part, largest
  first, with **Show** and **Recently deleted…**.
- **Tour, pseudocode, methods and report.** The **⋯** menu plays the history
  as a guided tour, copies it as pseudocode or a methods paragraph, and prints
  a report. Each of these reads any branch and names groups, squashed steps,
  bookmarks, notes and the tracker's abstract. The report draws its pictures
  from the live view when it is printed; no step stores a picture. **Show the
  state at…** shows the last step of the branch shown at or before a time.
  Under reduced motion the tour flies nowhere: each step is shown at once.
- **Run again only where it brings a result up to date.** A stale step's
  details offer **Run again** only for an analysis (`analysis.run`) and a
  highlight of a category or a range (`highlight.category`,
  `highlight.range`), the actions the catalogue marks `rerun`. A category or
  range group run again, or asked for again on its page, is updated in place
  with the cells the view shows now ("Updated the cells of …"), not added a
  second time; the same cells make no step, and no visible cell is refused
  with the group kept as it was. Other stale steps (a copied field, a kept
  view, a keyframe) no longer offer it: doing them again added a second one.
- **The selection tool travels with sessions** but is never a step: undo and
  jumps leave it alone. A tool this page lacks is left out, saying so.
- **Keep N cells of a category.** Click a category's count in the legend and
  type a number to keep that many of its cells, a seeded sample (`all` shows
  every cell again). This is a filter like any other: undoable and saved.

#### Sessions and sharing

- **The session file is the history (format 4).** **State › Session ›
  Share…** writes **History** (the trackers and branches chosen; a branch
  brings the steps it starts from), **Current state** (one step, "Shared
  state"; notes, conversations, who acted and times left out by default) or
  **Steps only** (a replay plan with no cell-level values, which opens as a
  replay). Include notes and bookmarks, conversations and Recently deleted
  only when chosen. Times can be exact, rounded to the day or left out, and
  the dataset as its source or its name only. Analysis results are kept or
  left out (left out, their settings are kept, and each analysis shown is
  run again as a step once the file opens, saying why when its inputs differ
  from the shared result's).
  **What's inside** lists exactly what the file holds. The file records what
  wrote it (the Cellucid version, the catalogue digest and each action's
  parameters version). A history recorded with older parameters keeps
  replaying.
- **Load… adds beside.** **State › Session › Load…** lists the saved states
  of this browser and opens a session file. A loaded file's trackers are
  added beside yours (renamed apart on a clash) and its active tracker is
  shown; nothing of yours is replaced. The notice after a load has **Undo**.
  A file is checked whole before anything changes. A session saved with
  Autoplay on plays its camera path once loaded.
- **What a share leaves out is left out of its values too.** Without exact
  times, a user field's creation time, an analysis result's computing time
  and a conversation's own times are rounded or left out as the steps' are;
  without conversations, no turn names a conversation or message.
- **Session files scale.** A 100,000-step history saves and opens. Shared
  sub-values are stored once, each chunk is compressed only where that pays,
  and arrays that only step parameters name are read after the step shown is
  on screen.
- **Starting views** of the five sample datasets are format-4 files (one
  tracker of one root).

#### My datasets

- **My datasets** (Session › Open data): the datasets opened in this browser,
  each reopened with one click with the cells it was last open with. Samples,
  links, repository datasets, Cellucid servers and files and folders are
  listed; where the browser keeps no handle, you choose the file again. The
  menu shows the five most recent. **All…** shows where each comes from, its
  storage and saved states, with **Continue** and **Forget**. Kept only in
  this browser (IndexedDB `cellucid-library`): a link in its public form,
  never its credentials, a file's path or its content. The startup's default
  sample and a notebook's dataset are not listed. **Forget** removes a
  dataset's saved states only when ticked. The Stored dialog shows My
  datasets as one row and lists the saved states.
- **Export… and Import…** write and read a `.cellucid-library` file with the
  datasets and saved states chosen, for Cellucid in another browser.
  Importing adds beside what this browser keeps, everything or nothing. A
  library file dropped on the page or chosen in **Open… › File…** is
  imported; a session file there is refused, saying it goes to **State ›
  Session › Load…**.

#### Agents

- **The chat bar.** Bring your own Anthropic or OpenAI key, kept for the tab,
  or in this browser with **Remember**. Requests go straight from the browser
  to the provider. Choose the model and effort from what the provider offers.
  The chat bar streams answers, shows expandable tool rows, and has Edit, Try
  again, Stop, Copy and Export. It counts tokens per turn and per chat
  against a budget, offers custom instructions and a `/` shortcut, and has a
  **Sent** view of what went to the provider. Conversations are kept with the
  history and, when chosen, in shared files. Opening another dataset or
  loading a session mid-turn stops the turn and keeps it with its dataset,
  unless the agent opened the dataset itself. The chat bar fails on its own,
  without stopping anything else.
- **Selection chip and choices.** A chip ("Selection: N cells ×") attaches
  the current highlight's page, group names and counts to each message, never
  the cells. An agent can offer 2 to 5 choices as buttons; a click sends the
  choice and keeps your draft.
- **Connect Claude or ChatGPT.** The **Agents** section gives the connector
  address, exactly `https://mcp.cellucid.com` (the Cellucid relay answers MCP
  at its root), and a box for a code typed on another device. An approval
  dialog pairs the tab, and **Manage…** lists, removes and forgets
  connections. Calls go to the tab focused last and survive a dropped link or
  a relay restart. A picture card (`ui://cellucid/view`) shows the view and
  the step it was taken after, with Refresh and Open in Cellucid. Pictures are
  sent as images. An agent can list datasets and open one (`datasets.list`,
  `dataset.open`): it is asked before a whole-file download, which names its
  size; `datasets.list` gives each of My datasets' size and what this browser
  keeps of it. An agent, Python or R can delete a saved state of this browser
  by its id (`history.saved.delete`, asked as deleting what was saved); the
  Load dialog's Delete goes through the same action.
- **Agents can do what a person can outside the history.** Play and stop the
  guided tour (`history.tour.start`, `history.tour.stop`); make the printable
  report as an `.html` file or as text (`history.report`); export the figure
  (`figure.export`), the highlighted cells as CSV or `.h5ad`
  (`cells.export`) and a field as CSV (`field.export`); list and remove what
  this browser keeps (`stored.list`, `stored.remove`) and forget a dataset
  of My datasets (`my-datasets.forget`); give Share's public link
  (`dataset.share-link`); read community annotation suggestions, vote and
  suggest a label (`annotation.suggestions`, `annotation.vote`,
  `annotation.suggest`, once the person has connected the repository and
  signed in). These are of a new effect, `outside`, with
  `history.share`, `history.save`, `history.saved.delete` and
  `history.anndata.save`: they record no step, Undo does not take them back,
  and each is asked first as writing a file, deleting what was saved or
  publishing. With **names and summaries only** on, an export still runs:
  the file goes to the person, and the agent is told its name, size and
  counts.
- **Permissions.** For each kind (changes in this tab, heavy work of more
  than about 10 s, writing files, deleting what was saved, large downloads,
  publishing, pausing the history) choose **Ask**, **Allow** or **Deny**, or
  one of three modes (**Ask before heavy or irreversible**, the default;
  **Allow everything undoable without asking**; **Read only**). An answer on
  a card can hold for "This chat" or "This connection". A confirm card names who asks and what it
  costs; Deny, a closed card or a closed tab refuse with one sentence. Each
  connection's tool list states what Cellucid will do for each tool. An
  activity line shows each working agent with **Stop** (also while a card
  waits) and **Undo turn**, and starts with `⏸` while recording is paused.
  The tool list says what History settings leave out of the history
  ("Camera moves are not recorded in this history"). A turn alone at the end of the branch is undone
  to the step before it (kept for redo); when someone else acted after it
  began, one dialog asks whether to take out only the turn's own steps (they
  go to Recently deleted) or everything since the turn began.
- **Agents see names and summaries only.** A notice says once what agents
  are sent. The switch, kept per dataset in this browser (`cellucid:agent-data`),
  keeps the values of single cells, extremes (minimum, maximum, median) and
  pictures out of every agent's results: the chat bar's, connected agents'
  and local agents'. Summaries of groups an agent shaped by naming cells are
  refused. Earlier results that could hold cell values are withheld when the
  chat re-sends the conversation. It is not anonymisation: a summary of a
  very small group can still reveal its cells' values.
- **Prompt-injection guard.** Every tool result, every refusal that quotes a
  name and the attached selection reach the model inside
  `<cellucid-data>…</cellucid-data>` with every `<` escaped. The instructions
  say early that text from others (labels, notes, dataset text) is data,
  never instructions.
- **Local agent apps.** An agent app on your computer (Claude Code, Claude
  Desktop, Cursor) drives a viewer page through the Python or R package's MCP
  server, with the connector's tools, permissions, confirm cards (also
  answerable in the app), turns, Stop and Undo turn, shown as
  "Claude Code (local)".

#### Notebooks

- **A notebook's commands are the viewer's actions.** Every action a Python
  or R notebook sends is dispatched as the notebook's: one step of the
  history, labelled as the notebook's, answered with its outcome and step
  (`action_result`) or its refusal (`command_error` naming the request). An
  action sent inside a turn carries the turn's id and label. The viewer
  reports every change of its history (`history_changed`), and a notebook can
  push a session file into the viewer (`load_session` / `session_loaded`).
  None is carried out while the person shows another dataset. A notebook's
  own data mapping change is reported back like any other (`data_mapping`,
  origin `notebook`).
- **Undo, redo and jumps are actions** (`history.undo`, `history.redo`,
  `history.jump`, effect `history`). They move the history and record no step.
  The history's own edits (branches, trackers, groups, squash, notes,
  recording, share, save and more) are actions too, for agents, Python and R.
- **An AnnData shown again after a pull or an apply opens.** `uns['cellucid']`
  may hold the packages' records (`history_applied`, `gene_sets`, `pulled`),
  which the viewer accepts and never reads. The `session` key is no longer
  accepted: nothing writes it.

#### Reading and exports

- **Save AnnData with this state…** (the History **⋯** menu; the action
  `history.anndata.save`, asked of an agent as writing a file) writes a new
  Cellucid-ready AnnData `.zarr.zip` (a MuData for several modalities) of
  the open dataset, every cell of it, with the state shown applied by the
  apply plan of cellucid-python and cellucid-r: renamed obs columns and
  categories, user-defined fields and highlight groups as obs columns, gene
  sets in `uns['cellucid']['gene_sets']` and as boolean var columns
  `cellucid_gene_set__<name>`, and a record in
  `uns['cellucid']['history_applied']` (package `web`). Cells a subsample
  did not open are missing or `<NA>`. It never writes over the dataset's own
  file: a save dialog where the browser has one, else a download. Cancel
  aborts the file chosen (no partial file) and removes the spool. What was
  left out is listed in the notification's details. A name a file cannot
  hold, a matrix still streaming in and the dataset's own file are refused
  in one sentence. The writer (`data/anndata/write/`, behind
  `BaseAnnDataAdapter.writeAnnDataCopy`) keeps the data as the source stores
  it, writes every expression matrix gene-major, and checks every name
  before a byte is written; it never writes, renames or replaces
  `obs['cellucid_original_position']`.
- **`uns` tables** of older AnnData files read correctly, and names in either
  Unicode form are found. Figure-export ZIP files take ZIP64 where they need
  it. Data kept in this browser by an earlier version is removed once, since
  its layout changed.
- **PNG export** works under the site's content security policy (`blob:`
  allowed).

#### Community annotation

- **One broken file no longer stops a Pull.** A user file or merges file
  that is not valid is left out of the consensus on its own, kept with its
  reason by blob id (not downloaded again until it changes), and reported in
  one warning with each file and its reason under its details; your own file
  is named first, and Publish will not overwrite it. For the open dataset only
  the profile and that dataset's section of a user file are read, so a broken
  section of another dataset no longer fails the file. Paths under
  `annotations/users/` that are not `ghid_<id>.json`, and case variants such
  as `annotations/Users/`, are ignored and listed.
- **Pull reads one commit and only what changed.** The branch's commit and
  the listings of `annotations/` come in one request; everything after is read
  by Git object id from that commit, so a Pull never mixes two commits. A Pull
  where nothing under `annotations/` changed takes one or two requests (was
  about seven); changed files are read many per request, in batches packed by
  their sizes (a first Pull of 5,000 files took 5,001 requests). Each batch is
  kept as it arrives, so an interrupted Pull resumes where it stopped (before,
  a repository under 32 MB cached nothing until the end). The file cache keeps
  each file's dataset view by blob id in IndexedDB; the index of file SHAs in
  `localStorage` is gone, and the session keeps only the SHAs of the files it
  publishes onto.
- **No fixed timeouts.** Every request to the Worker, GitHub's sign-in and
  repository lists, CAP and ORCID goes through the app's one transport: it is
  given up only when nothing arrives for the link's stall timeout, retried
  when that is transient, and paced by GitHub's `X-RateLimit-*` and
  `Retry-After` (a wait GitHub names is shown in the progress notification).
  A 1 MB publish on a slow uplink, or a steady 21-second download, no longer
  fails at 15 or 20 seconds.
- **Equal labels are one suggestion.** Suggestions with the same label in a
  category (compared after Unicode NFC; nothing else folded) count as one: the
  earliest proposal is shown, and each person's votes across them count once
  (their majority; a tie is no vote). Two people proposing "T cell" were
  reported as disputed, and a 4-to-3 majority split over two "T cell" cards
  could elect the minority label. Adding or renaming a suggestion to a label
  already shown is refused. A bucket nobody voted on is always pending, also
  with `minAnnotators: 0`. Consensus with moderation merges is computed once
  per change (it re-validated every merge on every call: 2.9 s for 5,000 users
  and 2,000 merges), and what Cellucid writes or exports is ordered by code
  units, the same in every browser.
- **Fork publishing no longer loses or conflicts annotations.** A
  contributor without write access has one pull request per base branch
  (head branch `cellucid-annotations/ghid_<id>/…` in their fork), carrying
  their whole file: every dataset they publish goes into it. Each Publish
  starts from the open pull request's file, resets the head branch to the
  base branch's current tip and writes the file on it, so the diff is that
  one file and an earlier squash-merged pull request never makes it conflict
  (before, each dataset had its own branch and pull request, and a reused
  branch stayed at an old tip). A direct publish that races another device
  creating the same file is reported as a conflict.
- **Cellucid no longer creates forks; the GitHub App asks for no
  Administration permission.** Creating a fork through the API needs
  Administration write, so the Worker's fork-creation route and the client's
  creation and five-minute readiness wait are removed. A contributor without
  write access forks the repository on GitHub; Cellucid finds the fork (its
  own name, else the repository's forks). Without one, Publish stops with
  `You have no fork of <owner>/<repo> yet: fork it on GitHub, then Publish
  again.` and a **Fork on GitHub** link to
  `https://github.com/<owner>/<repo>/fork`; a fork GitHub is still copying says `Publish again in a
  minute` (no fixed wait); a fork the app cannot write to (not installed on
  the contributor's account with access to it) says so with an **Install the
  app** link instead of disconnecting the repository as a lost access. App
  permissions: Contents read/write, Pull requests read/write, Metadata read.
- **Repository contract version 3** (`cellucid-annotation` main). User files
  and the merges file are version 3 (`user-v3`, `merges-v3`; config stays
  `config-v1`). Every suggestion id names its owner, `ghid_<id>_<rest>`
  (Cellucid writes `ghid_<id>_<random UUID>`), and must start with its
  `proposedBy`, so no contributor can write an id that collides with or
  claims another's; deletion markers name only the owner's own suggestions.
  No text has a length limit any more (labels were 120 characters, evidence
  2,000, comments 500, …): the only bound is the file, 1,000,000 UTF-8 bytes,
  the most GitHub's Contents API serves as JSON. Integers are at most
  9,007,199,254,740,991 (Min annotators was clamped to 50 in the panel). One
  whitespace class (JavaScript's `\s` together with Python's
  `str.isspace()`) is used by the schemas, the validator and the app, and
  JSON nesting is bounded at 512 levels in both. Cellucid accepts a repository
  only when its three schema files are byte for byte the template's, checked
  by Git blob id without downloading them; a schema that differs is named
  with how to copy it. Version-2 repositories are refused with that message;
  a version-2 user file is left out of the consensus with its reason (there
  is no migration). Local annotation state saved by an earlier Cellucid is
  reported, not read: clear it and Pull (published annotations stay on
  GitHub).
- **The template repository** (`cellucid-annotation`): authors copy
  `annotations/`, `scripts/validate_user_files.py`,
  `.github/workflows/validate.yml` and `.gitattributes`. The shipped
  `validate.yml` runs only the validator, with read-only permissions (it ran
  the template's own tests first, which stopped it in every repository set up
  from the instructions); the template's tests run from a separate
  `tests.yml`. The validator also checks the schemas byte for byte, rejects
  stray files under `annotations/users/` and case variants, and reports huge
  integers and deep nesting as errors instead of crashing. One shared
  385-case accept/reject corpus holds the validator and the app to the same
  verdict.
- **Worker contract 3.** New routes `/graphql/snapshot` and `/graphql/blobs`
  (GraphQL queries the Worker writes itself), the users folder read by tree
  id (the recursive tree route, cut at 7.5 MB, is gone), and a forced move of
  the contributor's own `cellucid-annotations/ghid_<id>/` branch in their own
  fork. Every write is bound to the signed-in user (their own user file, their
  own branches, pull requests from them); a refused write is `403`
  `not-applied`. Installations and repositories are listed one GitHub page
  per request with no 10,000 cap. GitHub's `ETag`, `Retry-After` and
  `X-RateLimit-*` headers are relayed, reads accept `If-None-Match`, and
  preflights are cached for 7200 seconds. A request has no total deadline: it
  is stopped only after 13 seconds without progress. A failed token exchange
  shows GitHub's reason (GitHub answers it under HTTP 200), and the redirect
  after installing the app ends on a short "installed" page instead of
  `400 Missing state`. Writes are checked with native JSON parsing, so a
  1 MB publish stays within the Workers Free plan's CPU. An older Worker is
  reported as one its owner must redeploy, a newer one asks for a reload, and
  an unreachable one is no longer called incompatible.
- **The Worker moved to `https://api.cellucid.com`** (a Cloudflare custom
  domain), so a WAF rate-limiting rule can cover it. The GitHub App's callback
  is `https://api.cellucid.com/auth/callback`.
- ORCID: one untidy record no longer empties a search, and the search text
  has no length limit. External links (installing the GitHub App, its
  settings) open reliably.
- **An agent's vote publishes only itself.** `annotation.vote` writes the
  person's file as published now plus that one vote (in a pull request when
  they cannot write to the repository, added to the file that pull request
  already proposes); their unpublished drafts and an author's staged settings
  and merges stay in the browser. A vote that changes nothing publishes
  nothing, the question names the vote, where it goes and what stays
  unpublished, and an agent's call never opens a dialog. New
  `annotation.suggest` proposes a label the same way (with the person's
  upvote; a label already suggested is refused with its id). Python's and R's
  `annotation_vote` and `annotation_suggest` ask the person first like an
  agent. Every confirmed action states its parameters in the card and in the
  connector's question. Stop or a dropped connection during the write says
  the vote may already be published, and Cellucid pulls at once.
- **"N changes not published · Publish"** in the panel counts the person's
  unpublished votes, suggestions, comments, deletions, profile changes and an
  author's settings and merges, item by item against the file as last pulled
  or published.
- **A local rename is a display name.** Votes, suggestions, the legend, the
  voting dialog, the agent's read and vote and the consensus column key on the
  dataset's own field and category names, so a category renamed in this tab no
  longer hides others' suggestions or publishes under the new name, and the
  annotated field stays locked through history jumps and session loads.
- **A second tab on the same dataset** pauses annotation there ("Open in
  another tab — annotation paused here", **Use here**), and the dataset
  controls keep working. A sign-in GitHub refuses (expired) signs out and
  keeps the repository and the unpublished work.
- **The consensus column** (`field.consensus`) records the source field's
  category → label mapping, the rule (threshold, minimum annotators) and
  where the votes came from (repository, commit, pull time, voters), and is
  rebuilt from them on replay; it defaults to the field's repository rule and
  `community_<field>`. Categories agreeing on one label (the NFC rule
  suggestions are bundled by) are one category of it. History steps that
  carry one code per cell replay only on the dataset they were recorded on.
- Names other people wrote (labels, field and category names) stay marked as
  data in every refusal an agent reads. The voting dialog's case-folding
  "Duplicate label" marker is gone: equal labels are one suggestion, and
  labels that differ are merged by an author. The cross-tab lock has no text
  caps, and the per-scope SHA indexes earlier versions kept in
  `localStorage` are deleted once when community annotation starts.

#### Fixes

- Enter in the agent chat takes the message at once, with Stop on offer; it used to wait for the model list first, with the message still in the input.
- A guided tour stopped before the panel draws again returns to the step that was shown.
- `analysis.run` right after setting an analysis's options waits for the groups or categories still being read, and runs. It was refused with "Reading the groups…", the Run row's progress line, when the read had not ended. It is refused only with a reason a read cannot change.
- On a slow GPU, a camera settling between pointer moves no longer redraws every pane.
- A number key pressed while a view is still switching embedding is done next, not dropped; with no embedding of that dimension it says so at once.
- In Safari, on a link where the stylesheet arrives after the page, Cellucid
  no longer starts before its styles are in (the analysis tools could fail to
  load there, "The analysis tools could not load"): the page is read past its
  stylesheet only once the stylesheet is in, in every browser.
- A stylesheet that does not load is said on the startup failure screen
  ("Cellucid could not load", Try again) instead of starting the app without
  its styles.
- When an agent's turn ends, its last reply is in the chat the moment Send
  comes back; it was drawn a frame later.
- A tab closed while the history was being written could, in Firefox, leave
  a saved state with a part of that write: it then could not be continued
  ("This saved state is missing messages of a conversation; it cannot be
  opened"). A write is now kept whole or not at all in every browser, and one
  handed to the browser is finished even when the tab closes right after.
- The Autosave tip and the documentation say what a closed tab can lose: the
  change of its last moment, a fraction of a second, and nothing written
  before it.
- Downloading the web app for a local viewer accepts a JavaScript file a
  host serves as `text/javascript` (RFC 9239) where the asset list says
  `application/javascript`, in Python and R alike; a file of another kind (an
  HTML error page) is still refused.
- Undo after deleting categories in place, and undo after colouring, keeping
  a view and renaming a category, work again.
- Reading an analysis's options, result or estimate no longer changes the
  analysis, so undo returns exactly to where it started. An analysis action on
  a hidden panel is refused with a sentence saying it is not shown.
- Analysis windows' places, sizes and folding are layout, as every other
  panel's; opening or copying a window waits for its first update.
- Deleting the last step of a branch no longer breaks the history.
- Redo pressed right after Undo (`Ctrl/⌘+Z`, then `Ctrl/⌘+Shift+Z`) is no
  longer dropped while the undo is still being shown: whether there is
  anything to undo or redo is asked of the history, not of the panel.
- **Replay again** no longer starts on top of the first replay, and a
  replay's mapping now maps genes and a category highlight's category.
- The State panel no longer rebuilds a step's details on a redraw that
  changed nothing, so a note being written is never lost.
- A session's cell sets are read in the three encodings Cellucid writes;
  the reader of category-relative sets that nothing wrote is removed (web,
  Python and R).
- Switching to another tab no longer breaks what follows. A hidden page lets go
  of the fields no view shows, and the next colouring, undo, jump or Verify
  failed with `Observation category field 4 requires exact dataset-length
  codes` (or the same for a continuous field). A field let go of is read again
  when used; a category filtered by cell count now stays in memory, as every
  other filter does.
- A field you made (a duplicate, a pushed column) is no longer lost when you
  switch tabs. It exists nowhere else, so it stays in memory.
- **Verify** gives no verdict when a step could not run because its field could
  not be read: it says so in one sentence and puts back what you were looking
  at, where it could stop on the scratch tracker and report the wrong thing.
- A history operation that fails on an internal check says that a problem inside
  Cellucid stopped it, with the check's text under Details, and no longer shows
  a stack trace.
- A hosted dataset redeployed while open keeps working. A redeployed site gives
  every file new validators although its bytes are the same, and every open tab
  said the dataset had changed; what was read is now compared, and only other
  bytes count as a change. A real change is refused as before, and no longer
  offers **Try again** or **Read the rest**, which could never succeed: the
  message says to open the dataset again.
- A dataset whose published sample view is of another session format, as happens
  for minutes during a release, opens on its default view with one short line
  instead of an error.
- A gene chosen while a MatrixMarket folder's matrix is still being read is
  read with its values once they arrive. It was refused as having "no value in
  any cell", because the page counted genes as read for every cell while the
  open cells grew; the gene sample now stays "none read" until the matrix says
  otherwise, and a gene read meanwhile waits for it.
- A field chosen from a menu stays chosen while it loads. A refresh of the
  menus meanwhile (the gene list arriving late, an annotation connection
  changing) showed **None** until the read ended; a stopped choice no longer
  leaves the menus disabled.
- When the browser's graphics stop working as Cellucid starts (a WebGL context
  already lost, as Safari gives after a GPU reset), startup says so and asks
  for a reload instead of failing on a shader. A field or gene with nothing to
  draw is refused before the active field changes.
- After the browser restores a lost graphics context (a GPU reset, a driver
  update), a black, white or grid background comes back as it was instead of
  transparent, and so does the chosen backdrop look with its motion.
- On a slow GPU, smoke, flowing velocity particles and a backdrop's motion no
  longer slow the page: their frames give way when the page's own frames run
  late, so clicks and changes are shown promptly. A GPU that keeps up draws
  them at full rate as before.
- **Search and link previews.** The description, canonical, Open Graph and
  Twitter tags are the first thing in the page head: a link preview reads only
  a prefix of the page (Slack the first 32 KB), and they used to start after
  the 116 KB import map, so a link to www.cellucid.com in Slack showed its
  title only. The structured data adds a WebSite node (the site name a search
  result shows), and the crawlable text, the JSON-LD and `llms.txt` name every
  format the app opens (they named AnnData only). `llms.txt` gains the key
  documentation pages, the MCP connector and how to cite. `relay/server.json`
  is the relay's entry for the official MCP Registry.

### Python

- **Every action of the viewer (viewer protocol schema 3 to 5).** The action
  catalogue ships with the package (`cellucid/action_catalog.json`, byte for
  byte the web app's, SHA-256 pinned in the three repositories). The viewer
  has one generated method per action (`viewer.color_by(field=...)`,
  `viewer.highlight_cells(cells=...)`, `viewer.history_get()`,
  `viewer.history_undo()`, `viewer.dataset_open(...)`, …) and
  `viewer.act(type, **params)` for any action named in data. Each returns an
  `ActionRequest` at once: `result(timeout=60.0)` is the viewer's answer and
  raises `ActionError` with its reason; `step` is the history step the action
  made. Cells are given as rows, a boolean mask or obs names. A notebook's
  viewer refuses `dataset.open`. Optional parameters are left out (`OMIT`)
  unless given, and `cellucid.jupyter.action_catalog()` returns the
  catalogue.
- **Wire.** The viewer's commands are `action` (text, or binary with its
  typed arrays), `requestSessionBundle` (answered with session format 4),
  `ping`, `debug_snapshot`, `freeze` and the agent commands; schema 2's
  `setColorBy`, `setVisibility`, `highlight`, `clearHighlights`,
  `resetCamera`, `setDataMapping`, `setGeneSet` and `setField` are gone. The
  new binary event `action_result` carries an action's outcome and step;
  `command_error` names the failed command's `requestId`, and `data_mapping`
  and `result` no longer carry one. `GET /_cellucid/protocol` lists the
  binary commands and events.
- **Removed, no aliases:** `highlight_cells(cell_indices, color, label)` (now
  the generated `highlight_cells(cells=..., color=..., label=..., page=...)`),
  `clear_highlights()` (→ `highlight_clear_all()`), `set_color_by(field)`
  (→ `color_by(field={"source": "obs", "key": ...})`), `set_visibility(...)`
  (→ `cells_visibility_set(...)`), `reset_view()` (→ `camera_reset(view="live")`),
  `set_data_mapping(**fields)` (→ `mapping_set(...)`), `push_gene_set(...)`
  (→ `gene_set_save(name=..., genes=...)`), `request_session_bundle()` and
  `SessionBundleRequest` (→ `request_session()` and `SessionRequest`),
  `CellucidSessionBundle` and `cellucid.session_bundle`. `viewer.push(values,
  key)` stays and returns its request.
- **Reading web sessions.** `cellucid.read_session(path)` (a
  `CellucidSession`) reads every format-4 file the web app writes (history,
  current state only, steps only; bounded chunks, parts, lazy parameter
  arrays) and checks all of it first against the history model
  (`cellucid/history_model.schema.json`, pinned in the three repositories). A
  damaged file says what is wrong; a file from a newer Cellucid asks you to
  update the package; a format-3 file is refused with the web app's sentence.
  The history is read with `cellucid.history`: `session.trackers`,
  `session.tracker(name_or_id)`, its `steps` table (with `kind`, `v` and
  `group`), `branches`, `groups`, `turns`, `notes`, `journal`, `stale()`,
  `methods(branch)` and `pseudocode(branch)` in the web app's wording, and
  `tracker.step(id)` with `state()` and `diff(other)`. `Session.chrome` holds
  the selection tool, read only. `notes_table()`, a squashed step's
  `squashed` (`cause`, `count`, `by`, `carried`) and each tracker's `paused`
  are read too.
- **Applying a history, a plan first.**
  `apply_cellucid_session_to_anndata(adata, session, *, tracker=None,
  branch=None, step=None, dry_run=True, overwrite=False, keep_session=False)`
  (and `session.apply(...)`) returns an `ApplyPlan`, one row per change with
  its status. `plan.apply()` (or `dry_run=False`) writes it and returns an
  `ApplyReport`. Applying the same step again changes nothing; a column you
  made is never replaced, and `overwrite=True` replaces only an earlier
  apply's columns. Each apply is recorded in
  `uns['cellucid']['history_applied']` (`cellucid.history_applied(adata)`).
  Works on AnnData and MuData. A steps-only file is read, not applied. Only
  what the chosen step needs is checked: the cells (count, subsample and
  order), the dataset id when `uns` names one, the data mapping and the
  fields the step uses. `pull()` of a session file writes the state it shows,
  checked the same way.
  Each gene set written or replaced is also a `bool` column of `var`,
  `cellucid_gene_set__<name>` (the highlight columns' naming rule), on the
  gene axis the data mapping reads (`X` or a layer; a MuData: the mapping's
  modality), `True` on the set's genes; missing genes, a mapping reading
  `raw` and a `var` column of that name the dataset owns are warned. The
  columns are pinned by the shared golden history
  (`apply.sequence[].var_columns`), as the web app's copy writes them. A plan
  that would write, rename or replace `obs['cellucid_original_position']` is
  refused, dry runs included.
  Removed, no aliases: `state_at()`, `steps()`, `apply_to_anndata()`,
  `dataset_source`, `required_fields`, `expected_dataset_id` (also of
  `pull()`), `apply_labels`, `add_highlights`, `inplace`,
  `ApplySummary` and the other options of the old apply.
- **The viewer's history, live.** `with viewer.turn("Normalise and
  cluster"):` makes every action sent inside one turn of the history (turns
  do not nest). `viewer.history` (`LiveHistory`) follows `history_changed`,
  and `on_history_change` hooks it. `viewer.load_session(path)` pushes a
  session file into the viewer. `viewer.apply_session_to_anndata(adata,
  request, *, tracker, branch, step, dry_run=True, overwrite, keep_session,
  timeout, cleanup_file)` plans or applies the session a `request_session()`
  brought, as above.
- **`cellucid mcp`** (`pip install 'cellucid[mcp]'`; `cellucid.mcp.run_mcp()`):
  a local MCP server over stdio for Claude Code, Claude Desktop and Cursor.
  It drives the newest viewer open on this computer, `--viewer` names one, and
  `--open PATH` serves a file in a viewer of its own. Every viewer with a
  channel, `cellucid serve` / `serve_anndata()` pages included, is published
  to agent apps on this computer through a record readable by your user only
  (`CELLUCID_AGENTS_DIR` overrides the folder; `CELLUCID_AGENTS=off` publishes
  nothing) and a token-checked loopback endpoint. The tools are generated
  from the catalogue (`cellucid/mcp_tools.json`). The `cellucid_status` tool
  lists the viewers found, and one sentence says what is missing when none
  is reachable. A call the app cancels is stopped in the viewer; a page
  reloaded or closed mid-call fails the call at once, and a page back makes
  the app list the tools again. The Viewer URL `cellucid serve` prints
  carries `viewerId` and `viewerToken` (the bare `/?anndata=true` address
  still opens the dataset). An agents folder that is not your user's alone
  is refused.
- `cellucid serve` and `show_anndata()` refuse a `.cellucid-session` or
  `.cellucid-library` file, saying what it is.
- On Windows, a file `cellucid serve` is serving can be saved over (a new
  file moved over it, as `write_h5ad` to another name and a rename does) or
  removed while the viewer reads it; it failed with "Access is denied". The
  viewer is then told the file is a new version.

### Releases

- **One way to release every package**: publish a GitHub Release named
  `vX.Y.Z` in its repository. `pypi-publish.yml` (PyPI) now runs on a
  published release only, not on a tag push, and `validate_release.py
  --release-event` refuses a draft, a release marked pre-release for an
  `X.Y.Z` version, and a release naming another tag; the tag, `main`, source
  and artifact gates are unchanged. `readthedocs.yml` no longer runs on a tag
  push (it built `latest` only). `cellucid-r`'s `release.yaml` ("Publish to
  CRAN") lost its manual trigger and makes the same refusals. The npm
  name-reservation package is released under the Cellucid version by
  `cellucid`'s new `npm-publish.yml` (npm trusted publishing, provenance),
  and `sweep_versions.py` holds its version to the release (the "independent"
  class is gone). A failed run is re-run with **Re-run jobs**.
- The headings of the released 0.9.3 to 0.9.6 carry their publication dates;
  the bump that opens a version dates the one released before it.
- The PyPI release run no longer requires the source distribution to match
  the conda recipe's `source.sha256` (`validate_release.py` lost
  `--match-recipe-digest`): every source edit made that digest stale and only
  the release noticed. It is recorded from the published PyPI file when the
  recipe goes to conda-forge (`scripts/publishing/PUBLISHING.md`).

### Docs

- **State and History** (web app) rewritten for History v2: the branch graph
  and its row rule, branches and trackers, recording, groups and squash,
  notes, bookmarks and search, Recently deleted, replay on another dataset,
  sharing and loading, History settings and saved states. **Save AnnData
  with this state** (what the copy holds, its refusals, Cancel, reading it in
  Python, showing it from R) and **Run again** are described there.
- The Python and R apply pages and the Cellucid-ready AnnData specification
  describe the gene-set `var` columns and the priority-order column rule.
- **My datasets** and **Keep N cells of a category** are documented.
- New chapter **Agents** (web app): the chat bar, connecting Claude or
  ChatGPT, permissions, and what agents are sent. New pages **Local agents
  (MCP)** for Python and R, the `cellucid mcp` CLI section and API page, and
  developer page **Local agent connection**.
- **Sessions & Sharing**, the privacy model and the site's privacy policy
  describe format 4, Share… and Load…, every browser store
  (`cellucid-history`, `cellucid-library`, …) and every agent storage key.
- **Community annotation** (web app): the version-3 repository contract
  (owner-prefixed suggestion ids, no text-length limits, schemas checked byte
  for byte, the files to copy and the shipped CI), what a Pull leaves out and
  what stops it, how a Pull and a fork publish work, equal labels as one
  suggestion, the Worker's contract 3, the GitHub App's token expiration
  (it must stay on; the author guide said to disable it, which breaks
  sign-in) and its permissions (no Administration: contributors fork on
  GitHub and install the app on their account with access to the fork), the
  GitHub App slug a self-hosted deployment must set, the GitHub
  login a user file publishes, and a new *Limits* section (GitHub's file,
  tree-listing and hourly bounds, Cloudflare's Free plan, direct-mode trust,
  writes not bound to the connected branch). Removed: a 10,000-file / 64 MB
  Pull limit and an "Overwrite remote" choice that do not exist.
- **Releasing Cellucid** (Contributing): the one release habit, the three
  repositories in order, the checklist, re-running a release, how changelog
  dates are set, and the one-time PyPI, CRAN and npm setup. The Python and R
  release pages, `PUBLISHING.md`, `cellucid-r/publishing.md` and
  `cellucid/npm/publishing.md` describe the GitHub Release only; tag pushes
  and uploads by hand are gone from them.
- Developer docs: the action layer, state domains and history (the journal,
  reads and staleness, the drift guards, the split control walk), the
  generated action reference, the agent data guard, the browser AnnData copy
  in the Cellucid-ready AnnData specification, and the Python and R session
  readers.
- **Search results, link previews and AI assistants.** Every page has a meta
  description and Open Graph and Twitter card tags (`sphinxext-opengraph`,
  pinned in the `docs` extra), so a shared docs link shows a card with the
  app's preview image; the home page is titled *Cellucid documentation*, and
  every guide page's title names its guide, so the Python, R and web app pages
  of one name are no longer duplicate titles. `sphinx-llms-txt` writes
  https://docs.cellucid.com/en/latest/llms.txt and `llms-full.txt`.

## [0.9.6] - 2026-11-01

Version 0.9.6 is a patch release on 0.9.5. Files prepared, sessions saved and
dataset catalogs written by 0.9.5 are read as before. The viewer protocol is
schema 2: use the Python and R packages of this release with the web app of
this release. The main changes are:

### Highlights

- **Subsampled reading.** Open a seeded share of a dataset's cells; only the
  parts of the file that hold them are read, and every part of Cellucid
  (view, filters, analyses, genes, exports, notebooks) works on those cells.
  The same cells in the web app, Python and R, in a link and in a session.
- **Files Cellucid writes store their cells in priority order**, so a share
  of the cells is the first rows of every array.
- **Datasets too large for the device open fitted, and large ones open
  progressively**, growing in place; a larger or smaller share is applied in
  place, reading only the cells added.
- **One notebook widget for Python and R** (`show_anndata()`, `cellucid.Widget`,
  R `cellucid_view()` in Jupyter): no port, no proxy.
- **Round trip in one call**: `viewer.pull(adata)` / `cellucid_pull()` write
  every highlight, label, score and gene set back onto the object;
  `push()` / `cellucid_push()` send yours in.
- **Up to 100 million cells in 3D** on a 16 GB laptop, with levels of detail
  that are uniform samples of the cells.

### Web app

#### Subsampled reading

- **Cells to open** (Session › Open data): one slider (or a typed share) with
  the exact count it stands for ("228,245 of 2,282,447 cells"). Below 100 %
  the method, its seed and what the method needs appear, with **Apply**.
  Datasets opened next in this tab take the same share; a method a dataset
  cannot use (a stratifying field or embedding it lacks) opens the same share
  by the file's default method, and says so.
- **Four methods, one specification.** `uniform` (every cell equally
  likely), `blocks` (runs of neighbouring cells aligned to how the file
  stores them, so parts holding only left-out cells are never read),
  `stratified` (proportional or the same count per category of a categorical
  or text field, optional minimum per category) and `spatial` (balanced over
  a grid on an embedding, keeping sparse regions). Every method is seeded;
  with `uniform` and `blocks` a smaller share is inside a larger one. The
  specification is published, with a conformance fixture shared by the web
  app, Python and R.
- **The file's default.** A file Cellucid wrote opens a share as `uniform`
  under the file's own seed (the first rows of every array); any other file
  as `blocks`, whose block size is set by the file's layout alone (about four
  chunks of its dominant per-cell array, or 64 KiB of it when it is read byte
  by byte), so the same file and share always open the same cells. The panel
  says what `uniform` would cost on such a file ("fetches 45 % of the matrix
  for these 10 % of the cells") and how far a `blocks` share's categories are
  from every cell's.
- **At least two cells.** One cell has no coordinate variation to draw (the
  rule a one-cell file is refused by): the control does not offer a share of
  one cell, a link asking for `cells=1` is refused with that reason, and
  fitting to the device never goes below two.
- **In a link.** `?data=…&cells=10%&sample=blocks&seed=0` (also `by`,
  `alloc`, `min`, `embedding`, `grid`, `block`) opens exactly that share,
  overriding fitting to the device; `cells=100%` opens every cell. Parsed
  strictly: a malformed value is named and the dataset opens as if none were
  given. A share chosen in the panel is written to the address as its exact
  count, method and seed (never one fitted to the device, never the cells),
  **Copy link** copies one, and the panel shows the cells' digest so two
  people can confirm they see the same cells.
- **Sessions keep their cells.** A session saved on a share records the cells
  (compactly), how they were chosen and their digest; loading it opens
  exactly those cells, reading the dataset again first when needed. Per-cell
  session data names cells by original position, so it applies whatever
  order a file stores its cells in.
- **Counts say "k of n".** The Dataset box ("228K of 2.3M cells"; Details:
  the method and seed), the filter line (*Showing 120 of 240 of 1,200
  cells*), categorical legends, the category builder, the screen-reader
  status line and the dataset picker all say how many of every cell are open.
- **Neighbour graphs on a share** keep the edges between open cells; the
  Dataset box says how much of the open cells' graph that is, and Data
  mapping offers every graph on the open cells.
- **A figure of a share says so**: SVG and PNG provenance name the open
  cells, how they were chosen and their digest.

#### Too large for the device, progressive opening, Apply in place

- **Fitted to the device.** When every cell would not fit this device's
  memory (77 bytes a cell, as measured with levels of detail on, against half
  the device's memory; or, for a matrix stored cell by cell, the index that
  makes its genes readable), a share sized to fit opens by the file's default
  method, with a notice and **Open all cells**. About 104 million cells fit a
  16 GB machine.
- **Progressive opening.** When reading the cells asked for is measured to
  take more than about a second (the link's speed and latency, the stored
  bytes of the embedding's parts the first view needs and their decoding, and
  how long this device took to show cells before), and a first view would
  arrive at least twice as soon, a first view of fewer cells shows at once and
  grows in place, doubling, to the cells asked for. Every step is a sample of
  them. A file stored in parts much coarser than its blocks can grow part by
  part, each part read once, ending exactly as opening the share directly.
  The same rule prices a file picked from this computer.
- **Apply in place.** A larger or smaller `uniform` or `blocks` share under
  the same seed (and every cell of a file Cellucid wrote) changes the open
  cells in place: raising reads only the added cells (their positions, their
  values of every loaded field and gene, the gene index for them alone),
  lowering reads nothing. Camera, every view's field, filters and settings,
  highlights, labels, pushed columns, hand-hidden cells and analysis windows
  stay; highlights and labels made from a column take the new cells as
  opening them would. Each per-cell structure ends bit for bit what opening
  at that share directly holds. On the graphics card only the added cells are
  sent, and lowering cuts the end. Another method or seed reads the dataset
  again through the parts already open, carrying the same state by original
  cell position.
- **Nothing read twice.** Another share of the same file, in the same or a
  new reader, takes what was already fetched and decoded (within the tab's
  memory budget); Zarr folders on a web host, picked folders and file sets
  are reused while their version is unchanged.

#### Reading at the floor

- **A share reads only its cells' parts**, element by element (obs of every
  encoding, embeddings, vector fields, graphs, `X` as CSR, CSC or dense,
  layers, `raw`, stored outlier quantiles, MuData modalities): a permanent
  test holds every element, at every share, to the stored chunks that hold a
  kept cell, each fetched once. Opening 1 % of CELLxGENE's 2.3M-cell Human
  Lung Cell Atlas: first view 3.8 MB in 143 requests at 8.2 s (was 20.3 MB in
  2,645 requests), at the floor computed from the file's layout.
- **The first view never waits for gene names.** The gene list (the var
  index, or the column Data mapping names) is read right after it, and no
  other var column; meanwhile the gene search, gene-based analyses and their
  pickers say *Reading gene names…*, and sessions, exports with genes and
  data-mapping changes wait.
- **Faster reads over a network**: HDF5 chunk-index nodes no longer read
  before a gene (the lung atlas read 12,382 of them first); index nodes of
  one task are read as one burst (all HDF5 formats; 100 Loom genes 11.1 →
  7.8 s over 150 ms / 5 MB/s); every S3 address of a public bucket is used;
  neighbouring reads are merged when that saves time; HTTP/1.1 hosts are
  recognised; sharded Zarr v3 reads all its shard indexes at once, then all
  inner chunks; uncompressed Zarr chunks are read in part; two columns of a
  latent space read only the open cells' rows.
- **Stricter, steadier reads.** A host answering another byte range, a file
  whose size, `ETag` or modification time changed while open, is reported,
  never read into what was read before. Busy hosts (429, 502, 503, 504) are
  waited for with spaced retries and fewer requests at once. A failed gene is
  read again on the next attempt. A shared read is cancelled only when every
  reader waiting on it has withdrawn.
- **Decoders.** Safari reads Zstandard-compressed Zarr about seven times
  faster (55 → ~395 MB/s); a Blosc frame with the delta flag is now decoded
  as c-blosc2 writes it (it was misread).

#### Genes on every path

- **Genes leave the viewer by identifier** (the var index) everywhere: gene
  sets, the gene coloured by, GMT files, DE and marker CSVs, `.h5ad`
  exports. When Data mapping names the genes by a var column (a CELLxGENE
  file's `feature_name`), the name shown goes alongside; named by the var
  index, nothing does.
- **The first gene of a cell-major matrix says what it costs** before it is
  read (bytes and time at the measured speed), with **Cancel**; on a share,
  its colours fill in as cells are indexed, each step a sample of the open
  cells. No silent rescans; errors say what failed.
- **The gene index is kept across visits.** For a matrix stored cell by cell
  (CSR, or dense row-chunked) the gene index is kept in this browser under
  the file, the matrix and how the cells were chosen (never the share):
  opening the same cells again, another day, reads genes at once, and a
  larger `uniform` or `blocks` share indexes only its new cells. Kept only for
  a file whose version is known; checked before use against its record, the
  file's row layout and, unless the version names the content, a fingerprint
  of the file; each gene's piece carries a CRC-32.
- **Fixed:** the first gene of a share of a Cell Ranger `.h5`; a text field
  for stratifying; SingleCellExperiment graph weights (were all 1); a
  MatrixMarket first gene or cell that looked like a header (`IDH1`,
  `cell0`).

#### Files that cannot be read in parts

- **Said up front, kept for next time.** An R object, a MatrixMarket folder
  or any file on a host that serves no parts of files: the opening says it is
  read whole, how much, and that a Cellucid-ready AnnData is read in parts
  instead. What it is read into is kept in this browser's storage across
  visits (for a file with a known version) and opens again in milliseconds,
  reading nothing; converted stores are written in priority order with a
  gene-major `X` (a gene of a 50,000-cell Seurat object: 7 ms). **Session ›
  Open data › Stored** counts what is kept, with **Clear**.
- **Checked as it is read.** Everything kept carries a CRC-32 of every
  16 KiB; damaged bytes fail the read as damaged, the copy is removed and the
  next opening reads the file again.
- **Safari keeps converted files and gene indexes across visits** like
  Chromium and Firefox; only private windows keep them in memory for the tab.
- **Safari: the Stored list reads a kept dataset without a failed load in
  the console.** The list read each kept dataset's last-used time while a
  reader opening that dataset was writing it; Safari refuses such a read
  now and then and logs it as a failed load. The read and the write now
  wait for each other.

#### Rendering at scale

- **Levels of detail are uniform samples.** A level is the first cells of
  the open cells' priority order (1,000 cells, then a quarter more at each
  level), so a reduced level shows the cloud's density and mix of categories
  as all cells do. Cells already in that order (a file Cellucid wrote, a
  share) draw every level at once (100M 3D: 0.33 s after levels of detail
  were turned on, was about 29 s); the spatial order culling uses is taken a
  few million cells a frame, without stopping a frame.
- **100 million cells in 3D.** Levels of detail at 100M 3D were silently
  unused (a GPU texture Chromium refused); fixed, with a notice if a GPU ever
  refuses. Positions and colours live in segments of 16,777,216 points, so
  no WebGL buffer limits the cell count; an array the browser will not
  allocate in one piece is made in WebAssembly memory (3D up to 357.9M
  cells); the velocity overlay's textures are layered. Memory per cell
  measured at 100M 3D: 54 bytes drawn, about 50 after levels of detail are
  built (was 65); 150M cells in 2D fit 16 GB.
- **Keep view and recolouring above 16.7M cells on Chromium** work (a
  vertex attribute released without a buffer kept stale bindings); released
  kept views free their GPU memory; GPU errors name the call that failed.

#### Analyses on a share

- **First answers on a sample.** Differential expression, Marker genes, Gene
  signature, Correlation, Composition, Dot plot, Detailed and Quick Insights
  show a first pass on a uniform sample of the open cells when a run is
  predicted to take more than about a second, then refine in place to every
  cell; the final result is bit for bit a single run, and only it is saved
  or exported. A first pass reads only its sample's part of each gene where
  the file allows.
- **Analyses follow the open cells.** When the open cells change, results are
  marked **Stale · Run** (its tooltip gives both counts); a stale result is
  not saved in a session, not exported and withdrawn from notebooks until it
  is run again. Differential expression and Marker genes no longer kept the
  earlier cells' groups after an Apply.
- **Stored outlier quantiles hold on a share** (they describe the whole
  dataset); quantiles computed on the open cells say they are an estimate.
- **Sparse statistics.** Differential expression and Marker genes read each
  gene's stored values, never a dense column (100M cells at 10 %: peak
  2.4 GB → 0.6 GB).
- Every analysis is checked against an independent reference on every way
  cells open (each method, progressive, Apply, sessions, served data) and
  every format and source.

#### Round trip and exports

- **Viewer protocol schema 2** (shared fixture re-pinned in all three
  repositories): every cell either side names is an original position; new
  events `result`, `result_removed`, `results`, `color_by`; new commands
  `setField`, `setGeneSet`; `ready` carries `n_obs`. Results carry an epoch,
  so a reloaded page replaces an older page's results.
- **Exports without a notebook**: the highlighted cells as an obs-name CSV
  or an `.h5ad` subset (written in the browser; anndata, h5py and Cellucid
  read it), any obs field as a CSV (**Export** beside **Duplicate**), gene
  sets as GMT. Exports list cells in the dataset's own order; genes by
  identifier, with the names shown beside when a var column names them; DE
  and marker CSVs name genes by identifier.
- **From a notebook**: pushed columns and gene sets show as fields and under
  **From notebook** in Gene signature and Dot plot, and are saved with the
  session.

#### Genes, legend and Session panel

- **Gray for zero**, off by default: a switch in the legend (for genes with
  no negative value) draws a gene's zeros in the theme's gray beneath the
  expressing cells, everywhere the colouring is drawn (2D/3D, smoke, kept
  and compared views, the cell card, figures). Per field and saved in
  sessions; a session saved with 0.9.5 restores with it off.
- **Genes read for**: under Cells to open, a slider of the share of the open
  cells a gene is read for (now, 1, 5, 10, 25, 50, 100 %), priced in cells,
  download and time at the connection's measured speed. **Read** grows the
  sample, reading only the added rows; **Stop** keeps what is read.
- **A compact Session panel**: Dataset, Open data, **Stored** (one line beside
  Manage…), Session file, then **Subsampling** and **Data mapping**,
  collapsed with a one-line summary. Minimum and Maximum of the visible range
  are one row each.
- **A gene chosen is shown, or the app says why.** Enter on a gene while more
  cells are opening colours by it; a gene chosen while the gene list is
  replaced follows to the new list; a gene or field that cannot be shown yet
  says so ("Gene names are still being read: CD3E can be chosen once they are
  listed."). Every gene search uses one lookup: case-insensitive, dotted ids
  as written, Ensembl accessions with or without a version. In the analysis
  gene picker, Enter chooses the gene typed (it chose None).
- **A dataset is called by its name**, never by its address, in openings,
  failures, the opening progress and sessions.

#### When a host or the connection fails

- **Whose problem it is, for every source.** An opening, a gene or an
  analysis names the provider, its host and the status: not answering (5xx),
  limiting this network (429, with the wait it asks for), refusing this
  network, private or expired (401/403), nothing at the address (404/410),
  or this computer offline. The requests made are under **Details**, with
  **Try again** where asking again may help; every place that reports a
  host's answer uses the same words. When unrelated hosts all stop answering
  at once, the message says it is most likely this computer's connection.
- **GEO no longer calls a link wrong while NCBI is down**: before a listing
  is taken to mean a file is missing, NCBI itself is asked.
- **A worker that cannot load its code says why** ("this browser is
  offline"), and a build that failed offline is tried again once online.

#### Steadier reading and drawing

- **Closing a dataset stops every read it started**, background reads
  included, and an opening that fails cancels the other requests it made.
  Closing the page ends what the dataset on screen was waiting for. Work
  cancelled before it began no longer runs (a MatrixMarket pass, the Google
  file picker, a catalog open).
- **A MatrixMarket folder reads its matrix once**, however its cells open
  (1M cells: 90.6 MB read once, genes in 6.5–6.8 s, was 24.8–25.8 s); a pass
  serves every read waiting for it. A folder whose matrix is still arriving
  opens in every engine.
- **Level of detail and culling are kept only where they are proven
  faster.** A frame's GPU time is read through the GPU timer where the
  browser has one, else bounded by the fence polls; the page turns them on
  when a frame of every cell takes over 33 ms on the GPU and keeps them only
  for a gain larger than it can measure, checking again after a restored
  graphics context, cells added or cut, or another representation. A
  caller's own choice (a notebook, the benchmark) is kept. The Performance
  tip says what was measured, and nothing is measured where nothing could
  be drawn fewer.
- **Nothing under the pointer moves when the page changes its own drawing**:
  the highlight count, the selection bar, the `Force LOD level` row and the
  Dataset box hold their place.
- **A failed Apply of Cells to open says so** and keeps the cells open
  before; a picked file is read again through the copy already held.
- **Exported figures keep their tick step in every browser** (a power of ten
  rounded one bit low in some engines).
- The welcome card fits phones and tablets, upright or on its side; the cell
  info card waits out a lost graphics context; no Firefox warnings about
  unflushed GPU fences.

### Python

- **`cellucid.subsample_cells()`** gives exactly the cells the viewer opens
  for a method, size and seed (any count from 1; the viewer opens at least 2).
- **`prepare()` stores cells in priority order** (`priority_seed=0`;
  `None` keeps the order given) in `.h5ad`, Zarr v2/v3, `.zarr.zip`, AnnData
  and MuData: `uns['cellucid']['priority'] = {seed, specVersion, nObs,
  levels}`, `obs['cellucid_original_position']` and, beside a sparse `X`,
  `varm['cellucid_prefix_counts']`. Chunks are smaller for shares: per-cell
  arrays 2^16 values, `X` data and indices 2^14. Preparing an object read
  from such a file puts it back in its original order first.
- **`cellucid.restore_original_order()`** puts an AnnData or MuData read from
  such a file (also a backed one) back in the order its cells were given. A
  record that no longer describes the cells (a subset, a reordering) is
  never trusted.
- **`show_anndata()` returns a notebook widget** (`cellucid.Widget`, an
  anywidget; `cellucid.Widget(adata)` makes one without showing it). The page
  and the viewer channel travel over the notebook's own connection: no port,
  proxy, `client_server_url` or `allowed_hosts`. Tested end to end in
  JupyterLab. Calls that wait for the viewer raise at once inside the cell
  that would wait, saying what to do. `anywidget` is a dependency;
  `jupyter-server-proxy` is no longer one. `AnnDataViewer` is the
  browser-tab viewer (its unused `height` is gone).
- **Served in priority order.** `show_anndata()`, `serve_anndata()`,
  `AnnDataViewer` and `cellucid serve` serve an in-memory AnnData, and any
  AnnData, MuData, Loom or 10x file `prepare()` did not write, as a Zarr store
  with the cells in priority order, made chunk by chunk as the viewer asks: a
  share transfers only its cells (1 % of a 1M-cell object: 0.18 MB, was
  2.1–10.6 MB). The temporary `.h5ad` is gone. An in-memory object is
  snapshotted before the call returns (memory-mapped, about 2 s per 100M
  values), its matrices stored as `prepare()` stores them. Files are served
  at disk speed (CSR `X` 7 → 109 MB/s). Each served file is listed with a
  strong `ETag`, so the browser keeps what it read across visits and sees a
  rewrite.
- **Round trip.** `viewer.pull(adata)` / `cellucid.pull(adata, viewer or
  paths)` writes every highlight, label, score and gene set onto `obs` /
  `uns` in one call (`<NA>` for cells a subsampled viewer did not open), from
  a live viewer or from the web app's CSV, GMT and session files;
  `viewer.push()`, `push_gene_set()`, `subset()`, `gene_sets()`, `results`,
  `selection` (obs names), `write_gmt()` / `read_gmt()`. Pull and subset
  refuse an object whose cells are in another order or that lacks the
  viewer's cells.
- **One session API.** `request = viewer.request_session_bundle()`, then
  `request.result()` or `viewer.apply_session_to_anndata(adata, request)`.
  `get_session_bundle()` is removed. Sessions saved on a share apply (cells
  not open are `<NA>`); `set_data_mapping()` no longer waits or returns.
- Viewers of a `prepare()` file name cells correctly in Python (obs names in
  original order); count messages read "1 cell".

### Docs

- **Fourteen runnable round-trip notebooks**, seven Python–R twins with
  "Open in Colab" badges: show a dataset in any notebook, pull results back,
  push clusters and gene sets, read the web app's exports, re-analyse an
  exported `.h5ad`, subsampled views; a gallery page lists them side by side.
- New pages: **Subsampled loading** (web app) and the **Subsampling
  specification**; the Cellucid-ready AnnData specification gains **Priority
  order**; the notebook pages describe the widget; the published-data
  tutorials open CELLxGENE's 2.3M-cell Human Lung Cell Atlas.
- Every page was checked against this release's behaviour (served stores,
  kept stores, sessions, analyses on a share, rendering limits).

## [0.9.5] - 2026-10-26

Version 0.9.5 is a patch release on 0.9.4 that breaks with every earlier
release: no export, session file or dataset catalog written for 0.9.4 is read.
**AnnData is the one data model.** The prepared export format is removed
from the web app, the Python package and the R package, with the export
server, the manifests, quantization and the `exports-v3` layout of the sample
datasets. `prepare()` and `cellucid_prepare()` now write a
**Cellucid-ready AnnData**, a standard AnnData that anndata reads and the web
app opens directly, and take a required `format=` argument. Session files are
**format 2**; a format-1 file is refused. Dataset catalogs (`datasets.json`)
are **schema 2** and list AnnData files. The contract between the writers, the
web app and the dataset repositories is the
[Cellucid-ready AnnData specification](https://docs.cellucid.com/en/latest/user_guide/python_package/c_data_preparation_api/09_cellucid_ready_anndata_specification.html)
(`docs/user_guide/python_package/c_data_preparation_api/09_cellucid_ready_anndata_specification.md`).
The main changes are:

- **One data model**: `.h5ad`, AnnData Zarr v2 and v3 folders and
  `.zarr.zip`, read by the web app's own HDF5 and Zarr readers, and MuData,
  Loom, 10x Cell Ranger, MatrixMarket folders and Seurat or
  SingleCellExperiment `.rds`/`.rda` files, each opened as AnnData; one input
  format registry names them in the web app, Python and R. The Python and R
  writers produce identical Zarr stores from the same input.
- **Open from a link**: CELLxGENE, GEO, GitHub, Zenodo, Figshare, Hugging
  Face, Dropbox, Google Drive, OneDrive and any web address, read in place with
  range requests that follow the measured link, from one address bar that also
  connects to a server; a private Google Drive or OneDrive file opens for a
  collaborator after one sign-in, with no Cellucid server involved.
- **`cellucid serve` and Jupyter hand the file to the browser**, which reads it
  itself; the server computes nothing per request. They take every format the
  web app reads (AnnData, MuData, Loom, 10x, MatrixMarket, R objects) and
  in-memory AnnData and MuData.
- **One viewer channel**: a notebook viewer and Python talk over one
  authenticated WebSocket; the notebook output is a plain iframe.
- **MuData**: `prepare(data=mdata)` writes a Cellucid-ready MuData, and the
  `modality` launch default names the modality a dataset opens on.
- **Stored outlier quantiles** and an outlier filter that reads or computes
  them on first use, cancellably.
- **Analysis, selection, figure export and edges off the main thread**, every
  analysis cancellable.
- **Sessions** save analysis results, Figure Export settings and the theme,
  in smaller files that restore in short tasks.

### Added

- **Python**: `prepare(modalities=[...])` chooses which modalities of a MuData are written, in that order, as cellucid-r's `cellucid_prepare(modalities = )` does: several are a Cellucid-ready MuData of those modalities (the root `var`/`varmap` hold theirs alone); one is a plain Cellucid-ready AnnData with the root's cells, obs, embeddings, graph and identity and that modality's `X`, `var` and `obsm` (it must have measured every cell; `gene_expression`, `var` and `gene_identifiers` then apply to it).
- **R**: `cellucid_apply_session()` writes a saved session's highlight groups and user-defined fields onto a Seurat object or an experiment, with the rules of `apply_cellucid_session_to_anndata()`; both are held to one web-app-written session file.

- **Web app, Python, R**: The Cellucid-ready AnnData. `uns['cellucid']` names the dataset: `schema` (1), `id` (the id sessions and community annotations key on), `name`, `description`, `citation` (`text`, `url`), `source` (`name`, `url`), the launch `defaults` (`embedding1d`/`2d`/`3d`, `knn`, `expression`, `geneIdColumn`, `latentKey`, `vectorFieldDefault`, `obsKeys`, `normalizeEmbeddings`, `centroidMinPoints`) and `outlierQuantiles`. Every key is optional; an unknown key refuses the file with a message. A server's explicit launch choices override the file's defaults, and the dataset panel shows the citation.
- **Python, R**: The writers store the outlier quantiles. With `latent_key`, `prepare()` and `cellucid_prepare()` compute every categorical obs column's latent-space outlier quantiles (categorical, text, boolean and nullable boolean columns; a column with no category of `centroid_min_points` cells, 10 by default, gets no entry) and store them as `uns['cellucid']['outlierQuantiles'] = {latentKey, minPoints, fields: {<column>: float32[n_obs]}}`, NaN where a cell has none. Python (`cellucid.prepare_data._outliers`), R (`src/latent_outliers.c`) and the web app (`data/latent-outliers.js`) are bit-identical, checked against a numpy conformance fixture. On Suo (562k cells, 14 categorical fields, a 24-column scANVI latent space) they add 25.5 MB and take 3 s in `prepare()`.
- **Web app**: The Session panel's **Link or server** address bar, and the `?data=` page parameter, open data where it lives, with no sign-in: a CELLxGENE dataset (its id or Explore link, read from the CZ CELLxGENE Census copy of the latest release), a GEO series (`GSE…`, its `.h5ad` supplementary file), GitHub (`gh:owner/repo/path`, `github.com` tree, blob and Raw links, `raw.githubusercontent.com`, Git LFS files read from `media.githubusercontent.com`), Zenodo (`zenodo:<record>`, record pages, file links and DOIs; a concept DOI opens the latest version), Hugging Face (`hf:<user>/<repo>[@rev][/path]`, pinned to the commit) and Dropbox file links, or any web address of an `.h5ad`, a `.zarr.zip`, a Zarr folder or a catalog folder. A GitHub link without `@branch` reads the repository's default branch, whatever it is called. Single files are read in place with HTTP range requests in the reader worker, so only what a view needs is downloaded. Each provider is one module behind one contract (`assets/js/data/sources/providers/`).
- **Web app**: The **Folder…** item of the Open… menu, and a folder dropped on the page, open an AnnData `.zarr` store folder, as `prepare()` and `cellucid_prepare()` write it, a Cell Ranger `outs/` folder or a MatrixMarket folder.
- **Web app**: AnnData files written before anndata 0.8 open: anndata 0.7 and 0.6 files and SeuratDisk `Convert()` output (0.7 dataframes with categories under `__categories`, Zarr stores without a root `.zattrs`, 0.6 compound `obs`/`var`/`obsm` with categories in `uns`, `raw.X`/`raw.var`, h5sparse-marked matrices, and every column a SeuratDisk file leaves out of its `column-order`). Each is read as a view of the current layout, so data mapping, annotation and sessions see the same dataset; files written by the real anndata 0.7.8 and 0.6.22 read value-for-value like their current-layout twins.
- **Web app**: A session opens its own dataset. The session manifest records where the dataset came from (`source`: `{kind: 'sample'|'link'|'remote'|'file', …}`, or null for a notebook); **Load State** reads it first and, when another dataset or none is open, opens that one, then restores. Only the fingerprint decides whether the session applies.
- **Web app**: Every analysis can be cancelled. Differential expression, Marker genes (and the heatmap genes it reads afterwards), Correlation, Gene signature, Quick insights, Detailed and the outlier quantiles run as one kind of cancellable job (`app/jobs/cancellable-job.js`): a box shows what the run covers, its progress, bytes read, rate and time left, with **Cancel**. Closing the window or switching dataset cancels the run. A cancelled run keeps nothing, and its reads stop in the reader worker unless another caller still waits for them. Cancelling one analysis worker task no longer ends the whole worker pool.
- **Web app**: The Legend's outlier control is shown for every categorical field and says what it can do: no latent space (with a *Choose latent space* button that opens Data mapping), no category large enough (the fewest cells it needs), not stored (what the first use reads, in MB), reading (progress with *Cancel*) or failed (why; moving the slider tries again). Data mapping marks the stored latent space "outliers stored".
- **Python**: `CellucidSessionBundle` reads session format 2: `dataset_source`, `required_fields`, and `cellucid.session_codecs` (`decode_cell_set` for all five highlight encodings, `digest_cell_set`, `decode_user_defined_codes`, `decode_typed_json`), vectorized with numpy. `apply_cellucid_session_to_anndata` refuses an AnnData missing any obs column or gene the session requires, naming all of them, before any write.
- **Python**: `cellucid serve` gains `--normalize-embeddings` / `--no-normalize-embeddings` and `--centroid-min-points N`, the launch choices `serve_anndata` and `show_anndata` already took.
- **R**: `cellucid_transition_drift(transition_matrix, embedding, normalize_rows)` computes each cell's drift in an embedding from a cell-cell transition matrix (CellRank's, for example), operation for operation as Python's `compute_transition_drift()`. A shared conformance fixture (`transition_drift_conformance.json`) holds both to the same drift: bit for bit where the arithmetic is exact, otherwise within the forward-error bound.
- **Datasets**: The five samples in `cellucid-datasets` are Cellucid-ready AnnData Zarr v2 folders, `anndata/<id>.zarr/`, written by `prepare(format="zarr")` and listed by `anndata/datasets.json` (schema 2, default `suo`, each with a format-2 starting view). They hold the same cells, genes, fields and coordinates as before (every 3D cell-order digest is unchanged), raw counts as float32 instead of 8-bit quantized values, the latent space with its stored outlier quantiles, and Pancreas's velocity under scVelo's names. Sizes: Suo 215.5 MB, Garcia 107.2 MB, He 32.3 MB, Kanemaru 50.2 MB, Pancreas 3.8 MB; no file over 1.6 MB. The app's sample catalog is `https://theislab.github.io/cellucid-datasets/anndata/`.
- **Demo datasets**: `cellucid-demo-custom-datasets` publishes its three datasets as a `catalog/` of Zarr v2 folders with a schema-2 `datasets.json`, and the 2-D dataset as an `.h5ad` and a `.zarr.zip` in `files/`, all written by `prepare()` and rebuilt byte for byte by `generate_datasets.py --check`. The README describes publishing with `.nojekyll` and the link forms.
- **Demo datasets**: `cellucid-demo-custom-datasets` adds a multimodal sample in two forms: `catalog/synthetic-citeseq.zarr/`, a Cellucid-ready MuData Zarr (60 cells, modalities `rna` and `prot`) written by `prepare(data=mdata, format="zarr")`, and `catalog/synthetic-citeseq-seurat.rds`, the same cells as a Seurat v5 object (assays `RNA` and `ADT`, reductions `pca` and `umap`, graph `RNA_snn`, `misc$cellucid`). The generator writes the `.rds` without R, as R serialization format 3 in gzip; R 4.5.3's `serialize()` of the object SeuratObject 5.4.0 builds from the same values is byte for byte the same. Both are listed in the catalog, so its **Sample dataset** list shows an R object.
- **Annotation**: Community annotation keeps every dataset separate. User files are version 2 (the profile once, then one section per dataset id) and every merge in `merges.json` (version 2) names its dataset; schemas `user-v2` and `merges-v2`. Repositories on the version-1 schemas are refused with a message.
- **Docs**: "Host a Dataset on GitHub" (why a Zarr folder is the way to publish on GitHub, `.nojekyll`, measured speed of GitHub Pages and `raw.githubusercontent.com`, GitHub's limits and Git LFS), "Which format?" in `prepare()`, `?cellucid_prepare`, the specification and the preparing guides, and a Memory budgets section in the developer docs.
- **Docs**: Tutorials, each run end to end before publishing: "Tutorials: open published data as it is" (a GEO MatrixMarket sample by link, a Seurat `.rds` on Zenodo by link, a Cell Ranger `outs/` folder from disk); Python "CITE-seq as a MuData: prepare and show" and "Serve Loom, 10x and MatrixMarket files"; R "Explore a Seurat object from R" (`cellucid_view()`, commands, a selection back in R), "Share a SingleCellExperiment on GitHub" and "CITE-seq: a Seurat object with RNA and ADT". The example-link tables of "Open from a link" tag each file Python or R, and list an R object and a GEO MatrixMarket sample.
- **Web app, Python**: Two data-mapping kinds: **Values**, whether each expression matrix (`X`, `raw.X`, each layer) holds `log1p` values and in which base, and **Highly variable genes**, the true/false `var` column that marks them (`raw.var` for `raw.X`). They are chosen under Data mapping › Genes and saved in sessions. The defaults are Scanpy's: with `uns['log1p']`, `X`, every layer and a non-integer `raw.X` are `log1p` in its `base` (none is e), else every matrix is as stored; `var['highly_variable']` when it is true/false, else none. A default read from Scanpy's records is marked guessed. The data mapping (sessions, the Jupyter wire, `cellucid.data_mapping.validate_data_mapping`) gains the keys `log1p` (`{"<matrix>": {"base": null|number}}`) and `highlyVariable`, pinned in the shared conformance fixture. Analyses read them through the gene loader: `log1p({matrix})` and `highlyVariableGenes()`, which reads the column on first use. The sample presets in `cellucid-datasets` carry the two keys.
- **Web app**: A categorical obs column opens in the colours its file stores the way Scanpy keeps them, `uns['<column>_colors']` (one hex colour per category, as `sc.pl.umap(adata, color="leiden")` leaves it), read as Scanpy reads it: a shorter list is not used, a longer one gives its first colours. The colours are read with the column's categories when it is first shown. Without them, or once categories are merged, the colours are Cellucid's palette. A session records a colour only where it differs from the one the field opens with.

- **Web app**: Every opening shows its progress from the first moment: a link (`?data=` or the address bar, every provider), a file, a `.zarr.zip` or a folder, a sample, a catalog dataset, a `cellucid serve` server and a notebook. While no dataset is shown, a panel over the empty viewer names what is opened and the step under way (*Finding the file on NCBI GEO…*, *Reading the file's structure…*, *Reading embedding X_umap…*), with a progress bar, how much has arrived, how fast and how long is left when that is known, and **Cancel**. While another dataset stays on screen, the same appears as a notification with Cancel. A request made again after a failure shows as *Retrying…* with why. Cancel stops the opening where it is (its reader worker and every request end) and keeps the dataset on screen; how an opening ended is kept in the notification center. Every layer that opens data reports into one opening (`data/opening-progress.js`, shown by `app/opening-progress-view.js`), so a link, the file it names and the dataset change that draws it show as one.
- **Web app**: **Loom, 10x Cell Ranger and MuData files open directly**, read in place like an `.h5ad` (only what a view needs), through layout views that present each file as AnnData (`data/anndata/views/`): no conversion step. **Loom** (`.loom`, Loom 2.0.1 and 3.0.0, loompy 2 and 3, velocyto, SCope/pySCENIC, loomR's orientation): `/matrix` and the layers as `X` and `layers` (a transposed view), cell attributes as `obs`, two-dimensional ones (`X_umap`) and coordinate columns (`_X`/`_Y`, `_tSNE1`/`_tSNE2`, `umap_1`/`umap_2`) and SCope's compound embeddings as `obsm`, `col_graphs` as `obsp`, Loom 2's XML character references decoded. **10x Cell Ranger HDF5** (`.h5`, Cell Ranger 1.2 to 9, one group per genome or `matrix/features`): `X` is the stored matrix as it is (CSC over barcodes is CSR over cells), var columns named as `scanpy.read_10x_h5` names them (`gene_ids`, `feature_types`, `genome`). A matrix has no embedding, so a lone `.h5` is refused saying to open its **Cell Ranger `outs/` folder** (Open… › Folder…, or a dropped folder), which opens with Cell Ranger's own `analysis/` projections (`X_umap`, `X_tsne`, `X_pca`) and clusterings (`graphclust`, `kmeans_<k>_clusters`) matched by barcode; `molecule_info.h5` is refused naming the matrix to open. **MuData** (`.h5mu` with its 512-byte user block, and MuData Zarr v2 and v3; mudata 0.2, 0.3 and 0.4): the global `obs`, `obsm` and `obsp`, each modality's obs columns and embeddings as `<modality>:<key>` (mudata 0.4 no longer copies them to the global `obs`), a MuData sharing features instead of cells refused. A Loom file with no embedding is refused with the two lines that compute one in Python and in R. The file picker, drops, links and the format registry list `.loom`, `.h5` and `.h5mu`. A corpus written by loompy 2.0.17 and 3.0.8, as Cell Ranger writes (with scanpy's Cell Ranger 1.2 test files), and by mudata 0.2.4, 0.3.2 and 0.4.1 (`scripts/formats-corpus/`, `tests/fixtures/{loom,tenx,mudata}/`) is read value for value against what loompy, scanpy and mudata read (`tests/formats-corpus-conformance.test.mjs`).
- **Web app, Python, R**: **Modalities**: a dataset of several feature spaces over the same cells (a MuData's modalities, the feature types of a 10x matrix such as Gene Expression and Antibody Capture, the genomes of a Cell Ranger 2 file) opens with its gene-expression modality, and Data mapping › Genes shows one **Modality** select above Expression for it (nothing changes for a dataset of one feature space). Gene names, expression, its `log1p` and the highly variable genes are then the modality's; gene search, colour by gene, Marker genes, Differential expression and Dot plot read it. A modality that measured only some cells reads NaN for the others (not 0), through MuData's `obsmap`. The 10x feature types are column windows of the one matrix and share its gene index. The data mapping (sessions, the Jupyter wire, `cellucid.data_mapping.validate_data_mapping`) gains the key `modality`, present only for a multimodal dataset, so every existing session, preset and fixture is unchanged; the launch defaults gain `modality`. The shared conformance fixture holds the rules (`defaultModality`, the modality a dataset opens with: the first named `rna`, `gex`, `gene expression`, `gene_expression` or `sct`, else the first; `mappingModality`, a mapping's modality checked against the dataset's), pinned in all three repositories; Python adds `default_modality()` and `require_mapping_modality()`. Every multimodal store is presented in MuData's layout (`/mod/<name>`, `mod-order`, `obsmap`), which is the one contract a writer of several modalities follows.
- **Web app**: **Seurat and SingleCellExperiment objects open directly** from `.rds`, `.rda` and `.RData` files: no conversion step and no R. The web app reads R's serialization itself (`assets/js/formats/rds/`: versions 2 and 3; gzip natively, bzip2 and xz with its own decoders, no WebAssembly; every item type R writes, including closures, byte code, environments, references, S4 objects and the ALTREP compact sequences, deferred strings and wrappers), in one pass over the file that keeps only what is shown (`scale.data`, `commands`, `tools`, neighbours, loadings and models are passed over; a Seurat v5 layer, whose name comes after its values, is passed over while the file streams and read in a second pass only if it turns out to be needed), into a store kept in the browser's origin-private file system for the tab's life (memory where there is none), opened as AnnData (`data/anndata/transcode/`). An uncompressed `.rds` (`compress = FALSE`, `readr::write_rds()`) is not copied: its matrices are read in place, from disk or a link. **Seurat** (SeuratObject 3.0 and later, `Assay`, `Assay5` with split layers joined, `SCTAssay`, Signac's `ChromatinAssay`): `meta.data` and `active.ident` as obs, every reduction as `X_<name>`, `graphs` as `obsp`, the default assay's `data` (else `counts`) as `X` with the other as a layer, variable features as `highly_variable`, other assays as modalities (the default assay first), Visium and FOV positions as `spatial`, `misc$cellucid` as `uns['cellucid']`. **SingleCellExperiment** (1.8 and later) and **SummarizedExperiment**: `colData`, `reducedDims` (`X_<name>`, case kept), `logcounts` (else `counts`) as `X`, `rowData`, `altExps` as modalities, `colPairs` as `obsp`, `metadata$cellucid`. A `save()` file of several objects keeps all of them in the one pass and the opening asks which to open (one button per object); choosing another later in the tab reads nothing again. Refused with the R line that fixes it: matrices kept outside the file (BPCells, HDF5Array/DelayedArray, `SaveSeuratRds(move = TRUE)`), Seurat v2 objects, h5Seurat files, objects with no embedding, zstd, ASCII and native-binary saves. On a 100,000-cell Seurat object (2.45 GB gzip, 5.3 GB decompressed, counts and data with 1,500 genes per cell, `scale.data`, PCA, UMAP, an SNN graph) the pass takes 1.3 times plain decompression (16.6 s in Node on an M1 Pro), keeps 3.5 GB on disk and peaks at 277 MB of memory; reopening it in the tab takes milliseconds. The opening shows the pass as one step with bytes, rate, time left and Cancel; a cancelled or failed pass leaves nothing on disk, and a closed tab's stores are removed at the next start.
- **Web app**: **MatrixMarket folders open** with **Open… › Folder…** or a drop: a count matrix (`.mtx`, `.mtx.gz`) with its feature and barcode lists as Cell Ranger 2 and later, STARsolo, kb-python, alevin-fry and GEO uploads lay them out, its orientation read from its size, written once into a cells × features matrix as the entries stream. A matrix that arrives in an order Cellucid cannot keep (a MatrixMarket file in no order, Seurat v5 split layers, an R `dgTMatrix`) is put in order in blocks of 4 Mi entries in browser storage, never whole in memory; one written by feature (kb-python, alevin-fry) is kept as it arrived, as a `csc_matrix`. Several feature types (`features.tsv` column 3: Gene Expression with Peaks or Antibody Capture) open as one modality per type, Gene Expression first, as a 10x `.h5` does. The folder's own cell tables give the embedding: Cell Ranger's `analysis/` projections and clusterings, and any `.csv`/`.tsv` table of the cells (`UMAP_1`/`UMAP_2`-style columns become embeddings, the rest obs fields; a clustering column such as `Cluster`, `seurat_clusters`, `leiden`, `louvain` or `kmeans_10` is categorical even when its labels are numbers); without one the folder is refused with how to compute one. A Cell Ranger `outs/` folder that also holds the `.h5` opens through the `.h5`, in place.
- **Web app**: The Dataset box says what became of an R object or a MatrixMarket folder read once into the browser: one line under its format, *Read once into this browser's storage for this tab*, *Opened from this browser's storage, read earlier in this tab* or *Read once into memory for this tab* (no origin-private file system, as in a Firefox private window), with *· N parts not read* when something was left out, and in Details the rows **Kept** (where, until the tab closes) and **Not read** (one line per part, named by its R path and its size in the file: `assays$RNA@scale.data (15 KB): Cellucid does not use it.`). The R parser records the bytes each dropped part spans and whether it held anything, so the empty `scale.data`, `commands` and `neighbors` every Seurat object carries are not listed; a Seurat v5 `scale.data` layer is listed too. The **Read** row of such a dataset says *Once, from start to end*, and its origin line names the format alone (*R object (.rds) · this computer*).
- **Web app**: The R and MatrixMarket corpora (`tests/fixtures/rds/corpus.zip`): objects written by SeuratObject 4.1.4 and 5.4.0 and SingleCellExperiment 1.32 in every compression and both serialization versions, SeuratObject's own `pbmc_small`, refusal cases, and MatrixMarket folders in each producer's layout (`scripts/rds-corpus/`, `scripts/mtx-corpus/`), read value for value against what R's accessors and scanpy read (`tests/rds-corpus-conformance.test.mjs`); every transcoded store also reads in anndata and mudata (`scripts/rds-corpus/verify_with_python.py`). The parser is checked item type by item type, at every split of its input (`tests/rds-parser-contract.test.mjs`), and the kept stores' lifetime with an origin-private file system on disk (`tests/transcode-space-contract.test.mjs`).
- **Web app**: **Figshare, Google Drive and OneDrive** open from a link (the address bar or `?data=`), read in place with range requests ([Google Drive, OneDrive and private sharing](https://docs.cellucid.com/en/latest/user_guide/web_app/b_data_loading/16_private_drive_onedrive.html)). **Figshare**: an item page (figshare.com or an institution's portal), a Figshare DOI or any DOI that resolves to one, or `figshare:<id>[v<n>][/<file>]`; the page asks `ndownloader.figshare.com` for the file's signed S3 link (its bot check refuses requests from a reader worker in some browsers) and the reader reads S3 with it for the link's stated ten seconds, then asks the page for the next with a one-byte request; private links (`figshare.com/s/…`) and bare download links are refused with the reason. **Google Drive**: a file shared with anyone opens with no sign-in through the Drive API and Cellucid's public, referrer-restricted API key, its details and first bytes asked for together (one round trip); a private file shows "Private file on Google Drive. [Sign in to Google]", signs in with the `drive.file` permission only and opens Google's file picker already showing that file, so one click on **Select** grants it. **OneDrive and SharePoint** (personal, work and school accounts): every file opens after **Sign in with Microsoft** (OneDrive serves nothing anonymously; PKCE in the page, no client secret, no library), with `Files.Read`, and `Files.Read.All` asked for only for a colleague's file (which many organisations let only an administrator approve); a share link is redeemed durably, so sessions and `?data=` name the file as `onedrive:<drive>/<item>`; the short-lived download address is renewed through Microsoft Graph once for all reads waiting on it, and a file changed meanwhile is refused. Folder links are refused (share one `.h5ad` or `.zarr.zip`). Sign-ins live in the tab's `sessionStorage` through one store shared with the GitHub sign-in of community annotations, never in a link, a session, the console, analytics or an error; an expired one shows a calm **Sign in again** under the address bar and as a notification, and the reads waiting on it resume. A read behind a sign-in or an expiring link asks the page for its URL and headers once for all parallel reads (`httpByteSource`'s new `access`), and any redirect that states its lifetime (`X-Amz-Expires`) is reused for that lifetime and no longer. Google answers only the origins Cellucid registers (www.cellucid.com and a local `cellucid serve` on port 8765, `provider-config.js` `registeredOrigins`), for public files and sign-ins alike; on any other page Cellucid says where Drive works before asking Google (a local server on another port: start it without `--port`; a notebook viewer on another port: stop the server holding 8765 and show the data again; another host: www.cellucid.com). Microsoft sign-in is refused up front, with the reason, on a plain `http://` page of another computer and inside another site's page (a notebook at `localhost` showing a viewer at `127.0.0.1`, whose storage the browser keeps apart from the sign-in window); a sign-in window found closed offers the step again, and a sign-in finished in it after all still counts until the step ends. The site's privacy policy is `www.cellucid.com/privacy.html`, linked as **Privacy** in the sidebar footer; the public client IDs, API key and project number are in `assets/js/data/sources/provider-config.js` (registered for www.cellucid.com and `http://127.0.0.1:8765` / `http://localhost:8765`), and a copy of the app without them says "not configured on this site yet". Checked against the real services: public Drive files open with no sign-in and are read in place, and Microsoft's sign-in page accepts all three return addresses. Google holds back downloads without a sign-in from a network that made many of them (an HTML "unusual traffic" page, HTTP 403, without CORS headers, while the file's details still answer and signed-in requests go through; seen three times on 2026-10-25, each after some 200–400 MB from one address within minutes, whether in 30 requests or 400, for about 50 minutes): its head read is asked once, and a failure confirmed with one byte (Google sends its page after about seven seconds, so this is told in about eight seconds); a file whose first bytes are held back opens after "Google Drive sends this file only with a Google sign-in for now. [Sign in to Google]" and the picker, and a read held back later says so, in place of the generic "no connection, or no CORS headers" (`registerSourceAccess`'s `unreadable`). A Drive file granted to Cellucid is read with the sign-in while one is held, public or not. `drive.google.com/u/<n>/uc?id=…` links (as gdown writes them) are Drive links, and the address of a OneDrive or SharePoint page that is not a share link is refused with how to copy the file's link. The Dataset box's read line ("5.5 MB read in 43 requests") is shown for a file read at an address with a query string (Google Drive's `…/files/<id>?alt=media`, a Dropbox `?dl=1` link): the box names a file by its address without the query, and the Transfers log now finds its entry by that address. The documentation lists public Google Drive files (Python and R) with what each first view read.
- **Web app**: Google Drive and OneDrive in every browser. Time spent in a sign-in window, a consent screen, Google's file picker or the choice of an R file's object never counts as a stalled download: every wait on the person registers itself (`utils/person-wait.js`), the startup deadline stands still meanwhile and starts over once the person has answered, and the opening goes on by itself after the sign-in; a data link's stall says "while opening the data link". Signed in to Google, a Drive link to a file not yet chosen for Cellucid (the `drive.file` sign-in reads only those) opens Google's picker on it at once and is read with the sign-in; closing the picker reads a public file with the key unless Google holds it back, and a read held back later says to open the file again with the Google account. A Drive download held back after the file opened is told in about 9 s, as at opening: the first read that gets no answer asks Drive, from the file's access, for the file's details and one byte, and the details answering while the byte fails ends the read in Drive's words (`unreadable.confirm()` in `sources/access.js`, asked by the reader worker with `source-unreadable`); a check that does not confirm the hold leaves the read's failure as it was. Signed Drive reads go to one stable address per request in flight (`lane=<n>`, a parameter Drive ignores) under the HTTP cache mode, so Chromium and Firefox preflight each address once: on a Drive-shaped host (0.65 s a request, 0.2 s a preflight) a signed first gene takes 5.3 s, as with no sign-in. The head a provider reads beside a file's details answers every read that starts inside it, so an R object from Drive or OneDrive is not read twice (the test kit's Seurat `.rds`: 9.4 MB in 3 requests). A OneDrive item is asked for whole, since Microsoft Graph leaves the download address out of an answer whose fields are selected (`$select`), even when the selection names it; an answer without what reading needs is said in Graph's terms ("Microsoft Graph described … without a download address"). **Sign in with Microsoft** opens its window empty, on Cellucid's own address, in the click, and sends it to Microsoft once the attempt's PKCE verifier is made, so Safari on a Mac with Microsoft's single sign-on extension (Company Portal) shows Microsoft's ordinary sign-in page for any account instead of handing the window to the extension; a second click brings back a sign-in window still open, and a pop-up blocker's closed window is said as a blocked window. Sign-in failures are said for what they are, with the next step: an organisation that lets only an administrator approve Cellucid (`AADSTS90094`, `90095`, `65001`, naming the organisation), an account it did not assign (`AADSTS50105`), its sign-in policy (`AADSTS53003`), a decline (`AADSTS65004`, which is also how Microsoft reports leaving its "Need admin approval" page, so both are named), a site registration Microsoft or Google does not accept, a Google Workspace administrator's block (`admin_policy_enforced`), and Google's window closed (Google shows "Access blocked" there). Where the same account would fail again the button offers **Use another account**, and the line links the administrator steps in the documentation.
- **Web app**: **Composition** (Analysis → Cells) asks whether cell-type proportions differ between conditions, accounting for samples. It counts the cells of each type per sample in one pass in a worker (a sample spanning conditions counts once per condition) and, with a condition, tests every type across the samples: logit proportions `(count + 0.5) / (total + 1)` with Welch's t-test or one-way ANOVA (after propeller, without its empirical-Bayes moderation), ranks of proportions (Wilcoxon or Kruskal-Wallis) or centred log-ratios, with Benjamini-Hochberg or Bonferroni over the types. Without samples it runs χ² tests (Fisher's exact test for a sparse 2 × 2 table) and says the cells were the replicates. It shows stacked bars per sample or condition, per-sample points per type, a table, the results as CSV and the counts as CSV; the statistics match SciPy on the shared fixture.
- **Web app**: **Dot plot** (Analysis → Genes) draws, for each gene of one or more named gene sets and each group of a categorical field, the mean expression as color and the fraction of cells expressing as dot size, as scanpy's `pl.dotplot` (or as a matrix plot, or as a stacked violin of each gene's distribution per group): scaled 0–1 or to z per gene, a free expression cutoff, the mean of expressing cells only, swapped axes, groups ordered by size or a dendrogram, Reds, Viridis or Blues. Only the listed genes are read, each summed per group in one pass over its non-zero entries in a worker; display changes never read again. Genes the dataset lacks are left out on a chip. It exports the long table (`group, gene, gene_group, mean, fraction_expressing, n_cells`) and opens Gene signature with its sets (**Score as signature**). Means and fractions match scanpy's `dot_color_df` and `dot_size_df` on the shared fixture.
- **Web app**: Marker genes, Differential expression and Dot plot have an **Expression matrix** option: the matrix the data mapping shows (default), or `X`, a layer or `raw` of the same file for that run alone, without changing what is shown. Each matrix's genes are named on its own axis (`raw` has `raw.var`), read through the same `var/{axis}/{names}` routes as a mapped matrix on every way of opening a dataset (a file, a link, `cellucid serve`, a notebook); a gene the matrix lacks is named as not found; Automatic fold change follows that matrix's own `log1p`. The choice is part of the marker cache key, the result line names it, and sessions save it with the options and the result.
- **Web app**: Differential expression has a **Replicates** option: with a sample field the run is a pseudobulk test, each sample's cells of a group one observation (its mean expression, from one pass over each gene's non-zero entries) and Welch's t between the two groups' samples, as `scipy.stats.ttest_ind(equal_var=False)` on the per-sample means; no variance moderation. Each group needs cells in at least 2 samples.
- **Web app**: Clicking a gene in a Marker genes or Differential expression result colours the cells by it through the sidebar's gene selector (**Coloring & filtering → Gene expression**), whose search box then names it.
- **Web app**: Analysis forms share one layout: the essentials, then one collapsed **More options** holding a collapsed group per topic, each counting the options changed from their defaults, with **Reset to defaults**; then **Run**, with the reason above it when it cannot run. Every option has an (i) tip and a default stated in the documentation; the options and which groups were open are saved in sessions. After a run, a change that would change the numbers dims the result with **Stale · Run**; display options (filters, ranking, labels, plot, rescaling) redraw it without computing again. Correlation and Gene signature run on **Run** instead of on every change.
- **Web app**: The site has **terms of use** (`www.cellucid.com/terms.html`: BSD 3-Clause license, no warranty, your data stays yours, other services keep their own terms) and an **accessibility statement** (`accessibility.html`: what works, the WebGL plot's limits, how to report a barrier), reading pages in the privacy policy's design, linked as **Terms** beside **Privacy** and **Accessibility** beside the licence in the sidebar footer; the privacy policy now names Google Analytics' two cookies. `/.well-known/security.txt` (RFC 9116) names the security policy's contact, its `Expires` renewed a year after each build by `scripts/generate-web-asset-inventory.mjs`.

- **Python, web app**: The viewer channel. A notebook viewer opens one WebSocket to its server, `<server>/_cellucid/ws`, named by `GET /_cellucid/protocol` (`{"channel": "_cellucid/ws", "commands": [...], "binaryCommands": [...], "events": [...]}`), and says `hello` with its viewer id and token; the server checks the `Host` and `Origin` headers (the server's own host, a loopback name or an `allowed_hosts` name), compares the token in constant time and answers `welcome`. Commands are JSON text frames, and a command that names cells (`highlight`, `setVisibility`) one binary frame: a `uint32` header length, the JSON header padded to a multiple of 4 bytes, then the cells as little-endian `uint32`, strictly ascending, which the viewer reads in place and checks as strictly as the hosts check its binary events (a host whose protocol document lists no `binaryCommands` is told to upgrade). Commands are queued until the viewer connects and sent in order, never twice; events are JSON text frames, and a selection (uint32 cell indices) and a session bundle (streamed to a temporary file) binary frames. Close codes: 1001 stopped, 4401 unauthorized, 4409 replaced by a newer page of the same viewer, 1008 protocol violation. The server pings every 20 s; the viewer reopens a dropped socket once, at once, and one health request tells a stopped host from a proxy that does not pass WebSockets. The rules are one document, `viewer-protocol.json` (web app `data/viewer-protocol.js`), kept byte for byte in cellucid-python and cellucid-r. Two million selected cells reach a hook in 89 ms.
- **Python**: `cellucid serve`, `serve_anndata()`, `show_anndata()` and `AnnDataViewer` take every format of the input format registry by path: `.h5ad`, AnnData Zarr (folder or `.zip`), `.h5mu`, MuData Zarr, `.loom`, 10x `.h5` files and Cell Ranger `outs/` folders, MatrixMarket folders (with GEO's `GSM…_` prefixes) and `.rds`/`.rda`/`.RData` (gzipped or not); `anndata.file.kind` and the health `format` carry the registry id. A folder of files is served file by file and listed in `/_cellucid/datasets` (`file.files: [{path, size}]`); the health and info `n_cells`/`n_genes` are `null` for a format Python does not read. An in-memory `MuData` is written once to a temporary Cellucid-ready `.h5mu`.
- **Python**: `prepare(data=...)` takes an AnnData or a MuData whole. Every component it holds (obs, X with var, obsm embeddings and vector fields, `obsp['connectivities']`, `uns['cellucid']`) is the default, an explicit argument replaces one component, and the output is byte-identical to passing the components by hand. `embeddings`, `dataset_id` and `dataset_name` are required only without `data`. The caller's arrays are never changed: a matrix that needs canonicalising (sorted indices, summed duplicates, dropped stored zeros) is copied first, and value ranges are checked without whole-matrix temporaries, so peak memory beyond the input is about one converted matrix. A backed AnnData (`backed="r"`) is prepared in two passes over row blocks of its file, never held whole.
- **Python**: `prepare()` writes a Cellucid-ready MuData from a MuData: `format="h5mu"`, or MuData Zarr with `"zarr"` / `"zarr.zip"`. Each modality is a Cellucid-ready AnnData, the shared cells, embeddings and graph are at the root, and modalities are written one at a time to bound memory (30,000 cells, 20,000 genes and 200 proteins: 4.8 s and 316 MB peak to a Zarr folder). New optional extra: `pip install "cellucid[mudata]"`.
- **Python**: New launch default `modality` (`prepare(modality=...)`, `uns['cellucid']['defaults']['modality']`, `serve_anndata`/`show_anndata(modality=...)`, `cellucid serve --modality`): the modality a multimodal dataset opens on.
- **Python**: `apply_cellucid_session_to_anndata` applies a session to a MuData; the added columns go to the global `mdata.obs`, and each modality's arrays are placed among the cells only when read.
- **Web app**: GEO links open MatrixMarket samples: a series' `GSM…_matrix.mtx.gz` with its `_features`/`_genes` and `_barcodes` files is one dataset, opened by its matrix (`geo:GSE…/GSM…_matrix.mtx.gz`), read from the sample's own folder with the sizes `filelist.txt` states. `GSM…` sample accessions are accepted. An incomplete or ambiguous sample is refused naming the file it lacks or the two that compete.
- **Web app**: `.rds.gz`, `.rda.gz` and `.RData.gz` (GEO's gzipped R objects) open from links, providers, drop and the File… picker; the content decides whether a file is gzipped once or twice, and it is decompressed as it streams. A host that sends an R file as `Content-Encoding: gzip` even for a part of it (NCBI's GEO file host does, for `saveRDS()`'s default gzip) cannot be read in parts by any browser: its file is downloaded whole once, decoded by the browser, and the notification says why (GSE327442's 25 MB Seurat object opens in 7–10 s in Chromium, Firefox and WebKit).
- **Web app**: A folder of MatrixMarket files sharing a sample prefix (an extracted GEO `_RAW.tar`) opens with Open… › Folder…; with several samples, Cellucid asks which to open.

### Changed

- **Python, R, web app**: A host program opens its viewer at `?viewer=true&viewerId=…&viewerToken=…` (it was `?jupyter=true`), and `/_cellucid/info` names the program serving it, `runtime`: `"python"` (`cellucid serve`, a viewer outside a notebook), `"jupyter"` (a viewer in a notebook kernel) or `"r"`. The viewer shows it in the Dataset box ("h5ad file · Jupyter", "Seurat object (.rds) · R") and in the dataset list ("R data"); its startup steps name "the viewer's host". The page keys and the runtime names are part of the shared viewer protocol file (`page`, `runtimes`; SHA `642c04e1…` in all three repositories). The old key opens nothing.
- **Python**: A notebook viewer (`show_anndata`, `AnnDataViewer` with `port=None`) serves on 8765, the one local port Google sign-in accepts, or on the next free port when another server holds it, saying so; it used a port the operating system picked, where Google sign-in never worked. cellucid-r's `cellucid_view()` already did this.
- **Python**: `prepare(format="h5ad")` (and `.h5mu`) compresses chunks on every core, written with HDF5's direct chunk write in grid order, so the file is laid out as HDF5's own writer lays it out: 1M cells (150 values each, 1.15 GB of X) in 13.3 s at 4.5 GB peak (63.9 s when HDF5 compresses one chunk after the other on one thread), 2M cells in 22.1 s, 200,000 cells in 2.8 s.
- **Python**: Applying a session copies a continuous field's missing values (NaN) as the viewer copied them, and accepts field names and highlight fields with `:` (ATAC peaks, MuData `<modality>:<column>` columns), as the viewer does; both were refused.
- **Python, development**: The `dev` extra brings the `mudata` extra, so a development install, CI and the release gate run the MuData tests instead of skipping them. On Python 3.11, where pip installs mudata 0.3, the suite (`tests/conftest.py`) sets mudata's `pull_on_update=False`, so every MuData it builds or reads behaves as under mudata 0.4. `tests/test_minimum_python_grammar.py` parses every Python file with the grammar of the oldest supported Python and refuses the f-strings only Python 3.12 accepts.

- **Web app, Python, R**: One input format registry decides what a file or folder is, everywhere (`assets/js/data/input-formats.js`). Each format has a stable id (`h5ad`, `zarr`), the name endings of each layout it arrives in (`file`: `.h5ad`; `zip`: `.zarr.zip`; `folder`: `.zarr`), container signatures and folder markers. The File… picker's `accept` list, a drop, the address bar's routing, the link providers (GitHub, Zenodo, GEO, Hugging Face, any web address), the dataset catalog, the served-store validator, the server health checks and the reader worker all ask it; the rules each of them spelled out before are gone. A name decides whether a file is taken; its content decides what it is, so a misnamed file (an `.h5ad` saved as `.zip`, a Zarr ZIP saved as `.h5ad`) opens as what it is, and an HDF5 file with a user block of up to 2 KiB is recognised by its signature. `tests/fixtures/input-formats.json` is the registry's document; cellucid-python ships the same bytes (`src/cellucid/input_formats.json`), from which `cellucid serve` classifies a path (`classify_anndata_path`) and names its `kind`, and cellucid-r ships them as `inst/input_formats.json`, which `cellucid_view()` and `cellucid_serve()` read; the three repositories pin one SHA-256. Every message that names what Cellucid opens (a link that is not a dataset, a catalog path, a Drive, OneDrive or Dropbox folder, a listing with nothing to open) lists the registry's single-file formats (`openableNamesText`).
- **Web app**: One data source reads every dataset file, `AnnDataFileSource` (`data/anndata-file-source.js`), addressed by one scheme, `anndata://<dataset id>/`; it replaces `H5adDataSource` and `ZarrDataSource` and their `h5ad://` and `zarr://` schemes; the data layer addresses an open file by its AnnData elements (`obs`, `var`, `obsp`, `obsm/{i}/{cols}`, `outliers/{l}/{f}`, …) instead of export file names. Its loaders are `loadFile` (a local or remote file), `loadFolder` and `loadRemoteFolder`; the format is a fact of the open dataset (`format`, `describeStorage()`: `{format, layout, zarrVersion, remote, url, name, bytes}`). A dataset that states no id is still `<format>_<name>` (`h5ad_<file>`, `zarr_<store>`), so every saved session keeps matching. A local file's source line names its format (`AnnData (.h5ad)`, `AnnData Zarr`).
- **Web app**: Opening has one entry, `openDatasetStore({source, expect})` in `data/anndata/open.js`, over one file `{file}`, a picked folder `{files, rootName}` or a web folder `{url, name}`; it opens the container and presents the root through a layout view when one recognises it (`presentStore`, `data/anndata/views/`), or transcodes a single-stream format into an AnnData Zarr store (`data/anndata/transcode/`). The view nodes `older-layouts.js` is built from are shared in `data/anndata/views/view-nodes.js`, with transposed, row-window, column-window and computed arrays for the layout views.
- **Web app**: The usage analytics' local load methods are named by how a dataset arrived, whatever format it holds: `local-user-file`, `local-user-zip` and `local-user-folder` (they were `local-user-h5ad`, `local-user-zarr-zip` and `local-user-zarr-folder`); `data_load` reports the input format it was read as (`data_format`: `h5ad`, `zarr`, `loom`, …) on success. The Dataset box names the format (`Loom (.loom)`, `MuData Zarr v3 folder`, `Cell Ranger outs/ folder`) and, for a multimodal dataset, its modality.
- **Python, breaking**: `prepare(out_path, *, format, embeddings, dataset_id, dataset_name, …)` writes a Cellucid-ready AnnData and returns its path. `format` is required, with no default: `"h5ad"` (one file, for single-file hosts such as Zenodo, Figshare, Hugging Face, S3 or a lab server, local files, `cellucid serve` and Jupyter), `"zarr"` (a `.zarr` folder of small files, the way to publish on GitHub: no file over 100 MB, no LFS) or `"zarr.zip"` (the Zarr store in one file). `out_path` must end in the format's suffix. `zarr_format` (2 or 3; 2 when None) applies to the Zarr formats only. New arguments: `citation_text`, `citation_url`, `latent_key` (an `embeddings` key), `normalize_embeddings`. Every obs and var column is stored; `obs_keys` sets which columns the viewer lists.
- **R, breaking**: `cellucid_prepare(out_path, format, embeddings, dataset_id, dataset_name, …)` writes the same Cellucid-ready AnnData as a Zarr store: `format` is `"zarr"` or `"zarr.zip"`, required, matched as `match.arg()` matches; `zarr_format = 3` writes Zarr v3. There is no `"h5ad"`: the package writes no HDF5 and depends on no HDF5 library.
- **Python, R**: The layout both writers share: `X` gene-major (CSC from sparse input; dense stays dense) in float32, or float64 when the genes kept hold integer counts above 2^24; `obs`/`var` with anndata's own encodings and pandas' category code widths (int8 below 127 categories, int16 below 32,767); `X/indices` as `<i4` unless the axis it indexes has 2^31 entries, `indptr` as `<i8` from 2^31 stored values; values that float32 rounds to zero dropped rather than stored as zeros; every embedding whole in `obsm` (float32, unnormalized); `obsp['connectivities']` as symmetric float32 CSR; chunks of 2^16 values for `X` and 2^18 for arrays read whole; deflate level 6, with byte shuffle on every multi-byte array except the expression values: in Zarr v2 numcodecs' `shuffle` filter and zlib, in Zarr v3 only codecs of the v3 specification, the `blosc` codec (byte shuffle, zlib level 6) on shuffled arrays and gzip on the others (`numcodecs.shuffle` is outside the specification; on He the v3 store is 29.00 MB, against 37.62 MB unshuffled); consolidated metadata. A `.zarr.zip` stores chunks first in grid order and the consolidated metadata last, with fixed timestamps and permissions, so it opens with the read of its end and the same input gives the same bytes. Python's `.h5ad` declares every element before writing values, so a reader gets the whole tree and the values opening needs from its first 512 KiB (431 KB on Suo); the file reads the same in anndata 0.8 to 0.12, scanpy and h5py. R and Python write identical Zarr v2 and v3 stores (every metadata document and every value), checked by cellucid-r's `tests/testthat/test-cross-writer.R`; R writes every NaN and NA as numpy's canonical quiet NaN.
- **Web app, Python, R, breaking**: A vector field binds by scVelo's naming on any embedding: `obsm['<field>_<basis>']` is the field `<field>` over `obsm['X_<basis>']` (`velocity_umap` over `X_umap`, `velocity_branches` over `X_branches`), the longest existing basis winning, and its id is `<field>` (`velocity`), so one field spans every basis. `prepare(vector_fields=...)` and `cellucid_prepare(vector_fields = ...)` write each array under that key, with as many columns as the embedding; arrays of another width, colliding keys and keys the viewer would read otherwise are refused before anything is written. `vector_field_default` names the field and must be drawn on an embedding the viewer shows at open. The rule is pinned in the shared data-mapping conformance fixture (web `parseVectorFieldKey`, Python `cellucid.data_mapping.parse_vector_field_key`, R `.parse_vector_field_key`).
- **Python, breaking**: `add_transition_drift_to_obsm(adata, T, basis=..., field_prefix=...)` computes the drift over `obsm['X_<basis>']` (every column) and stores `obsm['<field_prefix>_<basis>']`; `dim` is removed (use `basis="umap_3d"`). `compute_transition_drift` takes an embedding of any width, and a zero row sum is refused.
- **Python**: `cellucid serve <file>`, `serve_anndata()` and `show_anndata()` hand the AnnData store to the viewer, which reads it in place with range requests, as it reads a file or a link. The server serves the store's bytes at `/_cellucid/anndata/<dataset id>.h5ad` (or `.zarr.zip`, or `.zarr/<key>`) with exact `Range`, `HEAD` and strong `ETag` answers, over HTTP/1.1 with one thread per connection and `sendfile`, and names the dataset, its store and the launch choices in `/_cellucid/datasets`. Name and id come from the store's `uns['cellucid']` (`--dataset-name`/`--dataset-id` are needed only when it has none), and every launch choice is checked against the store before the server binds. An in-memory AnnData, or a store with a codec the browser cannot decode, is written once to a temporary `.h5ad`, one matrix at a time, and removed when the server stops. On a 562k-cell × 8,192-gene `.h5ad`: startup 0.16 s instead of 3.0 s, a UMAP about 40 ms instead of 170–250 ms, 100 genes 0.1–0.6 s instead of 115 s.
- **Python, breaking**: `AnnDataViewer` is the one notebook viewer class (`cellucid.jupyter.BaseViewer` is removed; same methods and signature).
- **Web app, breaking**: `?github=owner/repo/path` is replaced by `?data=gh:owner/repo/path`, and the **GitHub** control by the **Link or server** address bar; the source type `github-repo` is `link`. `?exportsBaseUrl=` is `?sampleCatalog=`, and the `cellucid-exports-base-url` meta element is `cellucid-sample-catalog`.
- **Web app, breaking**: `datasets.json` schema 2: `{"schema": 2, "default": "<id>", "datasets": [{"id", "path", "session"?}]}`, `path` a Zarr folder, an `.h5ad` or a `.zarr.zip` below the catalog folder, `session` (`{path, bytes, sha256}`) the dataset's starting view. The sample list and every folder link read this one format through one loader; listing reads each file's metadata only, and a catalog link without `dataset=` opens its default. A catalog in another schema is refused.
- **Web app, breaking**: Session format 2. Files start with `CELLUCID_SESSION/2\n`; any other file is refused before anything is read ("This file is not a Cellucid session file of format 2 (it does not start with CELLUCID_SESSION/2)."), in the web app, the Jupyter upload and `CellucidSessionBundle` alike. Sessions now save Figure Export settings, the docked analysis panel with analysis results (differential expression, marker genes, gene signature and correlation, as typed columns; floating windows carry theirs), per-view centroid toggles and the app theme, and no longer save derived labels, duplicated cameras, antialiasing or per-cell gene caches. Highlights are stored in the smallest of sorted gaps, a bitset, runs or relative to their category (checked by count and digest on restore); user-defined codes raw, as runs or bit-packed. The fingerprint records the dataset's identity id instead of how it was opened, and the manifest lists every obs field and gene the session uses (`requires`), so a dataset lacking any is refused naming every one. Restores request every field together, rebuild each kept view once and finish in background tabs and off-screen notebook frames. Suo (562k cells), a heavy session of 7 views: file 317 → 193 KB, save 189 → 101 ms, restore 1.89 → about 0.83 s, with no task over 50 ms.
- **Web app, breaking**: A data mapping is always `kind: "anndata"`; the `export` kind is deleted and a mapping or session naming it is refused.
- **Web app**: Own HDF5 and Zarr readers (`assets/js/formats/`) under one AnnData reader (`assets/js/data/anndata/`) replace the vendored h5wasm and the two separate loaders. They read Zarr v2 and v3 (sharded included), from a folder, a `.zarr.zip` or a web host; Zarr codecs gzip, zlib, Blosc, Zstandard, LZ4, shuffle and CRC-32C; HDF5 deflate, shuffle, Fletcher-32 and LZF. A codec or HDF5 filter outside these (Blosc or Zstandard inside an `.h5ad`, bz2 in Zarr) is refused naming it and the dataset, with how to save the file again. A Zarr store on a web host that hides paths starting with `_` or `.` (GitHub Pages without `.nojekyll`) is refused naming `.nojekyll` instead of read with empty cell and gene names. Suo (562k cells × 8,192 genes, CSR `.h5ad`) read locally in Node, against h5wasm: open 210 → 27 ms, first gene 1.9 → 0.85 s, neighbour graph 7.1 → 0.54 s.
- **Web app**: Reading follows the measured link (`formats/io/link.js`): every byte reader measures round trip and bandwidth as it reads, and request size, gap merging, read-ahead and the number of requests in flight follow the measurement and the host's protocol (six to an HTTP/1.1 host; a host that hides its protocol, as GitHub does, starts at HTTP/2's 100 streams until a burst shows otherwise). An `.h5ad`'s size probe carries its first 512 KiB in the same round trip; range reads bypass the HTTP cache so Chromium sends a file's ranges together; a link delivering at least 64 MB/s, a local file and a picked Zarr folder are read ahead like a disk. Suo `.h5ad` over a simulated 80 ms / 50 MB/s link: probe and open 98 ms (1.2× the round-trip bound), dataset ready in the app 1.3 s.
- **Web app**: An `.h5ad` on a web host is opened with 16 KiB metadata pages, since anndata writes each element's headers within a few KiB beside its data: a 1M-cell `.h5ad` reads 1.8 MB before its first view (1.4 s at 80 ms and 2 MB/s).
- **Web app**: A part of an unfiltered HDF5 chunk (anndata's default sparse `.h5ad`) is fetched as the bytes it needs, not the whole chunk: 1,000 cells of a 1M-cell CSR `.h5ad` at 30 ms and 5 MB/s read 4.4 MB (the whole chunks are 121 MB).
- **Web app**: Over a slow link, the first gene of a CSR `.h5ad` keeps the matrix values it reads (its entries spread over the whole matrix), so later genes read nothing (a 1M-cell file at 5 MB/s: the second gene in 6 ms, not 133 MB in 27 s). A cell subset weighs the bytes it would fetch by the link: an analysis of two large groups reads from the gene index in memory (5.5 s) instead of 400 MB of rows (117 s).
- **Web app**: Opening a file reads no category dictionary. Every obs column is described from its metadata, with its category count, and a categorical's categories (and the colours Scanpy stored for them) are read when the column is first shown. A GEO file whose barcodes are stored as a categorical (21,751 categories, 1.1 MB) no longer reads them before its first view.
- **Web app**: Every request into a Zarr store bypasses the HTTP cache (`cache: 'no-store'`), as range reads of a file do: WebKit never finished delivering to the reader worker the bodies of store chunks served from its cache, and a store's metadata is never paired with chunks cached from an earlier version.
- **Web app**: A refusal written for the person (a record with several files, a file with nothing to plot, a folder link with no `datasets.json`, an AnnData too old to read) is shown as it is instead of a generic network sentence. A folder link without `datasets.json` says so for every host: "The linked data could not be opened: https://example.org/atlas/ holds no datasets.json, so it is not a dataset catalog. Link to the folder that holds datasets.json …".
- **Web app**: Analysis computes each gene from its non-zero entries, read as `{length, indices, values}` and never installed in field state: Wilcoxon and Welch statistics treat the zeros as one rank block, bit-identical to the dense kernels for differential expression and Wilcoxon markers (Welch marker p-values within 2.2e-9 relative). Differential expression sends the groups once and each gene as its column; a cell subset is read by whichever of rows or columns costs fewer bytes. The worker pool runs one worker per core less the page and reader (no cap of 8). In Chromium on Suo: marker discovery by `cell_type` (69 groups × 8,192 genes, Wilcoxon) 12.3 → 0.86 s with no long task; differential expression of one 266-cell type against the rest 14.1 → 2.1 s.
- **Web app**: Gene clustering in the Genes panel is exact agglomerative clustering on the analysis workers (Müllner's nearest-neighbour-candidate algorithm, double precision; the same merges as before, trees equal to scipy's to 5e-15): 3,000 genes × 69 groups 5.0 s on the main thread → 0.15 s. Memory is 8 bytes per gene pair (was about 104), so the default budget fits about 11,580 genes instead of 2,270. Gene signatures read each gene once as its sparse column (100 genes 450 → 32 ms), CSV export of page data and signatures is written by the workers, categorical page data is returned as codes, and Correlation computes pages concurrently and refuses a categorical X or Y.
- **Web app**: Box and violin points, and the volcano plot's genes, are one image per trace drawn by a worker at exactly the positions Plotly lays out, instead of one SVG element per point; hover and SVG/PNG export are unchanged. Suo, 544,130 cells in two pages: box with all points 5.0 → 0.13 s, violin 5.2 → 0.23 s, hover no longer blocks for 0.7–0.8 s.
- **Web app**: Figure export computes bounds, draws the point passes, writes vector circles and encodes PNGs in a worker; the page only composes. Where a browser's workers have no WebGL2 -- WebKitGTK and WPE ([WebKit bug 209809](https://bugs.webkit.org/show_bug.cgi?id=209809)), Safari 16.4–16.6 and 17.0–17.1 on macOS Ventura ([WebKit bug 183720](https://bugs.webkit.org/show_bug.cgi?id=183720)) -- the worker says so once and the same code draws the point passes on the page. Output is byte-identical apart from the timestamp, except that Full Vector now honours a depth sort at every point count (it skipped one above 200,000 points). Suo: main-thread long tasks of up to 804 ms → none; 600 dpi PNG about 1.1 → 0.6 s; Full Vector SVG peak memory +239…+581 MB → +15…+114 MB.
- **Web app**: The lasso projects the pane's visible cells once per press and tests only the cells a new vertex can flip, bit for bit the full even-odd test: 0.06–0.11 ms per pointer move on Suo with a 900-vertex lasso (was 0.93–1.3 s). Proximity re-tests only the radius band that changed; KNN expands on the neighbour lists the reader worker built. Selections, highlights, undo history and the Python selection wire share one representation, an ascending `Uint32Array`; selected cell sets are identical before and after over 510 exact comparisons.
- **Web app**: The reader worker builds a neighbour graph's page payload as it checks the matrix: edges arrive already in the renderer's fixed random order as 12-byte records, with one ascending neighbour list per cell, and a transform-feedback pass gathers edge positions on the GPU, so a filter change gathers nothing on the CPU. With edges on (Suo, 6.3M edges, K = 250,000), a legend click takes 6–8 ms instead of 19–57 ms; frames are pixel-identical. Embeddings and vector fields are checked, normalized and bounded in the worker (2D embedding main-thread time 71.6 → 7.3 ms), and no spatial tree is built on the main thread.
- **Web app**: Colouring and filtering touch each cell once, with byte-identical output over 45 interactions. Main-thread time on a 562,000-cell dataset: picking a loaded gene 18 → 5 ms, a legend click 12 → 6 ms, an outlier-slider step 22–24 → 9–12 ms, keeping a view 56 → 40 ms. A view whose filter inputs match an array another view holds reuses it, and GL error reads happen only after new storage is allocated, so restoring a 7-view session runs with no task over 50 ms. Smoke rebuilds only when its inputs change.
- **Web app**: Outlier quantiles are read from the file, one field at a time, on the filter's first use; for any other file or latent space they are computed in the reader worker, every field asked for together in one pass over the latent space, which is then kept within its budget. Quantiles that no longer fit the file's categories are computed instead. The latent space is no longer downloaded when a dataset opens. Sessions never drop an outlier threshold silently: one that cannot be restored is left at 100% and named.
- **Web app**: Memory budgets are disjoint shares of one allowance, half the estimated device memory (`navigator.deviceMemory`, raised past an 8 GiB report by four times the JavaScript heap limit; 8 GiB where neither is reported). Decode workers are sized by cores and that budget; switching expression matrices keeps each built gene index while they fit, so switching back reads nothing; the default analysis budget is a sixteenth of device memory.
- **Web app**: A server or notebook dataset (`cellucid serve`, `serve_anndata`, `show_anndata`) is opened with the browser's AnnData reader, as a linked file is, under the server's dataset id, name and launch defaults; its dataset list is described from metadata and a dataset is opened when chosen. The `remote://` and `jupyter://` schemes and the server-side mapping routes are gone.
- **Web app**: Sessions and community annotations follow the dataset, not the route: a dataset opened from a sample, a link, a server or a local file is the same dataset. A local file is named by `uns['cellucid']['name']`, else `uns['title']`, else its file name.
- **Web app**: The sidebar gene search and the analysis gene selector share one windowed list (`app/ui/components/windowed-listbox.js`) with the same keys: arrows wrap at either end, Home/End jump to the ends, and every option carries `aria-posinset`/`aria-setsize`. Field renames, deletions and restores touch only what changed (Suo: gene delete 20 → 2.3 ms).
- **Datasets**: The datasets are built by `scripts/build_atlas.py` and `scripts/build_pancreas.py` through `scripts/pinned_build.py`, which pins the environment, inputs, analysis and the producer modules `prepare()` executes, reads every staged store back with anndata, and compares its digest with `sources/<id>.json` (schema 4). Starting views were re-captured on each store as format-2 sessions; new `scripts/capture_session_preset.mjs`, `scripts/verify_catalog.mjs` and `scripts/local_site.mjs`. CI opens every store with anndata 0.12.19.
- **Reproducibility**: `cellucid-reproducibility` states that every recorded run measured the export-based Cellucid (web `0ab545aed`, cellucid-python 0.9.1) and what a re-run on the current version needs; its setup checks out the measured commits. Nothing was re-run.
- **Web app**: The Session panel is redesigned. The **Dataset** box describes only the dataset that is open, rewritten on every dataset change: its name, cells and genes, format and host (`Zarr v2 folder · theislab.github.io`, `h5ad file · this computer`), its own description (two lines, **More**), how much has been read in place, and a collapsed **Details** with source and citation links, file, size, obs fields, embeddings, layers, raw, neighbour graph, latent space and vector fields; nothing is cut mid-word. One **Open…** button opens a File… / Folder… menu (a browser dialog picks files or a folder, never both) in place of the **H5AD**, **Zarr ZIP** and **Zarr folder** buttons. One **Link or server** address bar with **Open** replaces **Server** and **Open from link**; what is typed decides (a port, `host:port` or `localhost` is a server, a provider form or a file address is a link, any other web address is a server when one answers `/_cellucid/health`), and a connection shows as one line with **Disconnect**. The **Transfers** block is removed: how much has been read is one line in the Dataset box, and a retried read is a notification. The sample list shows `—` while a file, a single-dataset link or a server is open, and lists a connected catalog's datasets, named "Dataset from <host>", instead of mixing sources. **Data mapping** no longer opens on its own: the best-guess notification's **Show** opens it, its rows are one flat list, and arrays that are off and wider than three columns, or cannot be drawn, are listed under a collapsed **Other arrays**. Without a latent space the outlier filter is one **Choose latent space** button; with one, the slider and at most one short line. A Zarr v2 folder chosen in Safari, whose folder dialog leaves out hidden files, is refused with how to open it instead (drag the folder onto the page, or its `.zarr.zip`).
- **Web app**: While the first dataset is still opening, the sample list already offers every entry of `datasets.json` under the name its path gives it, the one being opened marked *· opening…*; names and cell counts follow once the files are read.
- **Web app**: A dataset that does not open is said once, in one notification of one sentence (*Could not open "Garcia".*), its raw error folded under the notification's **Details** and in the console; the explanation and **Try again** stay next to the selector. Opening and connection failures keep their raw error under **Details** too, never in the sentence.
- **Web app**: The Highlighting tool's gesture line (and *Color by a field first…*) and the Antialiasing line (*Chosen automatically: …*) are the first line of their (i) tips; under the Tool buttons the sidebar shows only what a gesture did.
- **Web app**: The *Signed in to Google* line with **Sign out** goes when the sign-in expires: the token store says when a sign-in is stored, removed or expires.

- **Web app**: A link to a host that ignores range requests opens: its file is downloaded once into memory, with its progress shown, and opened from there. A host that does not let other sites read its files says so and what to do: "`<host>` did not let this site read `<file>`: it does not allow other websites to read its files (CORS), or it could not be reached. Download the file and open it from your computer instead."
- **Web app**: Opening a link starts the AnnData reader worker while the link is resolved, so the worker and its modules are loaded by the time the file is found. A CELLxGENE link probes the file in the Census release the browser read last beside the release list, so it usually opens one round trip sooner; a newer release is followed. A Zenodo link that names its file asks for the file's first bytes beside the record, so it opens one round trip sooner too. Opening a link again in the same page (Disconnect and open, a session reopening its dataset) asks its provider nothing and starts from the first bytes read the first time.
- **Web app**: Analysis windows and floating sidebar panels resize from all four edges and corners, by mouse, pen, finger and keyboard (arrow keys on the focused edge, 16 px or 64 px with Shift, Home/End to the least and largest size; `Alt`+arrows move the window). They are held only to their content's least size and the viewport: the 500 px height cap and the twice-the-opening-width cap are gone. Edges snap to other windows, the sidebar, the middle of the free area, the screen and the content height, with a guide line, and show the size while they move. A double-click on an edge fits it; **Maximize** (or a double-click on the header) fills the free area beside the sidebar. A new window of an analysis opens at the size the last one was given. On a screen 720 px wide or less, analysis windows are bottom sheets at 40%, 75% or 100% of the screen, with a grip, a collapsing header and a tab strip for the others.
- **Web app**: The expanded analysis view moves and resizes by mouse, pen and finger (it was mouse-only) from every edge and corner, its splitters and bottom-right corner answer arrow keys, a double-click on its header maximizes it, its least size is its content's (not 600 × 500), and its size and splits are remembered per analysis. On a screen 768 px wide or less it fills the screen with a touch grip on its split.
- **Web app**: The Analysis panel groups its rows: **Genes** (Marker genes, Differential expression, Dot plot, Gene signature), **Cells** (Composition) and **Explore** (Quick, Detailed, Correlation).
- **Web app, breaking**: The `analysis/windows` session chunk is `{windows, sizes}`: each window records its own box (not the header-high box of a collapsed window, which it used to save), whether it is open and its stacking layer, and `sizes` holds the remembered window and expanded-view sizes per analysis. A chunk without them is refused.
- **Web app, development**: The web contract suite needs Node.js 22.15 or newer (`engines.node`); its Zstandard contracts compress their inputs with `node:zlib`'s Zstandard encoder, which Node.js 20 (end of life) lacks. CI runs Node 22, 24 and 26.
- **Web app**: Marker genes compares every group of a field against the rest or a chosen reference group, for every group or a chosen subset, inside all cells, the visible cells or a page. A group under **Min cells per group** is never an error: it is skipped and listed (default), pooled with the other small groups, or tested anyway and flagged low n, and chips on the result name them. The test follows scanpy's `rank_genes_groups`: Wilcoxon (score = the z of the rank sum, with AUROC), Welch's t-test or the t-test with overestimated variance; p-values exact or continuity-corrected (SciPy, R) or scanpy's normal approximation with or without tie correction; Benjamini–Hochberg or Bonferroni over the genes tested; fold changes as stored or with `log1p` undone (automatic when the data mapping marks the matrix `log1p`). The result keeps every gene's statistics for every group, so direction, filters (adjusted p, |log2FC|, % in group, % in reference), ranking (score, adjusted p, log2FC, AUROC, % difference) and the number shown re-apply without testing again; it shows as a list, a heatmap or a clustered heatmap. Genes can be excluded by family (mitochondrial, ribosomal, hemoglobin) or pattern, or limited to the highly variable genes or a list. It exports the shown markers and every gene as CSV with scanpy's column names plus `auroc`, `small_group` and `tested`, and the markers as GMT; **Score as signature** and **Dot plot** hand them on. The statistics match scanpy and SciPy on the shared fixture.
- **Web app**: Differential expression compares two pages or two sets of categories of a field (disease against healthy), optionally within the visible cells, a page or categories of another field, with **Swap A and B**. A group under **Min cells per group** is tested and flagged low n instead of refused; only overlapping groups and a group under 2 cells disable **Run**, with the reason. It has Marker genes' test options, free adjusted-p, |log2FC| and % expressing filters that redraw the volcano without testing again, and exports scanpy's column names with `pct_nz_group`, `pct_nz_reference`, `auroc`, the group means and counts.
- **Web app**: Gene signature scores one or more named gene sets — typed, pasted with any separators, imported (`.txt`, `.csv`/`.tsv` in long or wide shape, `.gmt`) or taken from Marker genes — by mean, median, mean z-score, scanpy's module score (`tl.score_genes`: control genes drawn from expression bins with a seed) or UCell (ranks, with a max rank). A name is matched exactly, or else to the one gene matching it ignoring case; names the dataset lacks are left out and listed, never an error. Rescaling (z-score, 0–1), the plot and the layout of several sets redraw without scoring again. The CSV is `page,set,cell_index,score`. Module score matches scanpy and UCell a numpy reference on the shared fixture.
- **Web app**: Gene-streaming analyses (Marker genes, Differential expression, Gene signature) size themselves: a run reads ahead as many genes as the workers consume while one gene is fetched, from the measured link and each gene's measured size and test time, within the analysis memory budget, and adapts while it runs (more read-ahead when the workers wait on reads, kept only when it makes the run faster; less on a throughput drop or memory pressure). Analysis workers start one at a time, so building the pool never holds the page. How many requests go out at once is the reader's scheduler's to decide. Nothing it decides changes a statistic.

- **Web app**: One design for every screen of the site: the startup failure screen, the code-load failure screen, the no-JavaScript screen, the 404 page, the privacy policy and the sign-in return page are the app's brand card (wordmark, heading, one sentence, clear actions) in the theme chosen in the app; the standalone pages load a small `page.css` instead of the whole app stylesheet. The design principles are in `assets/css/README.md`.
- **Web app**: The startup failure screen says what kind of problem it was (the data could not be opened, a server did not answer, this browser cannot run Cellucid, Cellucid could not start), folds the exact message (address, host, status) under **Details**, and offers **Try again** and **Open Cellucid without this link**. A module that fails to download shows the same screen, not only a red status line.
- **Web app**: Data a link names that does not open no longer stops Cellucid. A `?data=` link, a `?remote=` server, a `?dataset=` or a catalog's default dataset that cannot be read leaves the app running without it: a panel over the empty viewer says the data could not be opened, with the exact message under **Details** (and in the notification log), **Try again**, and **Open without this link**, which removes those settings from the address. The full-screen card is kept for the app itself failing (no WebGL2, its code did not load, its own server). The panel says the reason itself, as its notification does ("This 10x Cell Ranger .h5 file has no cell embedding …"; "… holds no datasets.json …"; "The owner of x.rds has turned off downloading for it"), and when asking again cannot change the answer (a file or folder that is not there, a refusal, a file Cellucid does not read), **Open without this link** comes first, before **Try again**. Google Drive and OneDrive refusals and a Drive read held back are shown as written wherever they arrive.
- **Docs**: The documentation site takes the app's look: Inter, the app's ink and cyan colours in light and dark, square surfaces, the app's favicon, and the wordmark without a duplicate "Cellucid" beside it.
- **Docs**: Every screenshot is recaptured from the current app by its scenario, with a `captures.json` record naming the clean web commit it shows, so none is left without provenance; captions and alt text describe what each image shows. New figures for the R objects, Loom, Cell Ranger `outs/`, MatrixMarket and `.h5mu` Modality pages, the opening progress and failure panels, the Google Drive and OneDrive sign-in lines, Composition and Dot plot, the R viewer, and the phone layout. The JupyterLab images are replaced by the viewer `show_anndata()` serves and a transcript of a real `viewer.debug_connection()` report, and the Firefox and WebKit captures are gone. The screenshot tool takes a per-scenario `viewport` (with `touch` for a phone) and grows a crop to hold the whole drawn pointer.
- **Docs**: A new home page: the tagline and one sentence, **Open Cellucid** and the 60-second tour as calls to action at a documentation size, the formats it opens, a fresh light and dark capture of the Suo sample in 3D, the three ways in (web app, Python, R), four features with their own captures, compact capability and learning cards, and citation and contribution; it has no sidebars. The outdated close-up screenshot is gone. The images are WebP captures from new `home-*` scenarios of the screenshot tool, which now writes WebP (`format: 'webp'`, a fixed `quality`) as well as PNG, with provenance in `captures.json`.
- **Web app**: One `favicon.ico` (the root one) for every page.
- **Web app**: A published sample view that cannot be applied says so in one sentence; the reason is under the notification's **Details**.

- **Python, R, web app**: Highlights and visibility changes carry their cells in a binary frame, not a JSON list. `viewer.highlight_cells()` and `viewer.set_visibility()` take a list of ints or a one-dimensional NumPy integer array (the fast form), in any order, each once; the cells are checked and sorted once and written straight into the frame. `set_visibility([])` is refused (pass `None` for every cell). Two million cells are a 7.6 MB frame instead of 15 MB of JSON: from NumPy, Python builds it in 3 ms (the JSON took 0.35 s) and the viewer has applied the highlight 0.22 s after the call; ten million cells, 13 ms (1.7 s) and 0.86 s. R builds the frame of ten million cells in 16 ms (the JSON took 4.2 s).
- **Python**: `viewer.send_message()` returns `True` when the command was sent and `False` when it waits for the viewer to connect; it raises only after `viewer.stop()`. `get_session_bundle()` raises at once when the viewer reports it could not capture the bundle, or sent an invalid one. Hook callbacks run on the viewer channel's thread, one event at a time and in order; a hook that raises is logged and the channel stays open. The event is recorded in `viewer.state` and handed to `wait_for_event()` before any hook runs, and the other hooks still run: one buggy `on_message` hook no longer drops the event, so `get_session_bundle()` no longer times out behind it (several failures are logged together, as one `ExceptionGroup`). `debug_connection()` reports `channel_connected`; the `debug_snapshot` event no longer has `parentOrigin`, and no event carries `viewerId`.
- **Docs**: "Architecture: message routing (HTTP vs postMessage)" is "Architecture: the viewer channel" (`e_jupyter_hooks/05_architecture_viewer_channel`).
- **Web app**: In Full quality, the highlight pulse redraws only the layer over the points. While highlights pulse, each pane's base layer (backdrop, grid, every point, the vector-field overlay) is kept and every pulse frame draws the highlights, edges, centroids and orbit anchor over a copy of it, pixel for pixel the frame drawn whole (`rendering/pulse-layer-cache.js`; `tests/browser/highlight-pulse-layer-parity.spec.mjs`, one view and a grid, three engines). GPU time of a pulse frame (Chromium, Metal, antialiasing off): 10M cells with 20,000 highlighted 9.97 → 0.48 ms, with 5M highlighted 38.9 → 24.1 ms; 1.5M cells with 15,000 highlighted 4.91 → 0.69 ms. With antialiasing on, a multisampled pane keeps its own samples instead (WebGL2 cannot copy one multisampled buffer into another): the ring pass marks the samples it covers in the stencil, and a pulse frame resets and draws again only those, every point still submitted but only those under a ring shaded (`rendering/pulse-sample-mask.js`; the same spec, the pulse moved between two phases, three engines). 10M cells, antialiasing on: 20,000 highlighted cells in one region 40.3 → 7.5 ms, 300,000 in one region 70.0 → 15.7 ms, 20,000 scattered over the embedding 37.2 → 25.3 ms, 1M scattered 44.3 → 45.0 ms (their rings cover the embedding).
- **Web app**: Applying a highlight costs the highlighted cells, not the dataset: the page counts no longer copy a page's cells into an array to count them, nor to digest them (5M of 10M cells: 85 ms of the apply); 1M highlighted cells of 10M apply in 24–26 ms (57 ms).
- **Web app**: Startup no longer runs as one ~90 ms task: it yields between its first stages (the interface, the WebGL programs), and the page's first layout (about 25 ms) is done once, where the welcome dialog takes focus, instead of also being forced by checks that measured the page (the dialog's focus check and the notification center's empty card stack). Chromium: no main-thread task over 50 ms while the app starts. A PNG figure's pixels are read in bands of rows with a yield between them, so the read-back of a 300 dpi figure is no longer one long task.
- **Web app**: Publishing a large dataset no longer holds the main thread. Its GPU stores are allocated in a task of their own just before the publication, once the GPU has done the work queued before (a WebGL allocation is checked by a `getError` round trip that waits for the GPU: up to 30 ms of the publishing task), so the publication reads nothing back; its positions and colors go to the GPU with the first frame (the colors were sent twice); and each per-cell array is written once: the first visibility pass adopts the new dataset's all-shown visibility, clearing the active field leaves gray colors gray and republishes no unchanged transparency, the highlight list of a dataset with no highlight is empty instead of being derived from ten million bytes, a fresh mask is proven binary as it is made, and the 0/1 per-cell bytes (validation, the R8 alpha texture of the live and kept views) are read four cells to a word. The first slice of positions for the spatial-tree worker is copied in its own task. Chromium (Metal), 10,000,000 cells in 3D: opening 227–237 ms → no task over 50 ms, switching to it 337–349 ms → none (largest 45 ms).
- **Web app**: Coloring by a category of a large dataset no longer holds the main thread: before the field is shown, its per-category totals and its centroids for the shown representation are computed in slices of about 8 ms (the same additions in the same order, so the same centroids), the codes are proven against their categories once when they load, the colors are written without the transparency the visibility pass writes next, and the visible count of each category is read from the totals when the only filter is the field's own (a category is then all shown or all hidden). Hiding or showing one category reads no mask when every cell is shown. Chromium, 10,000,000 cells: coloring by a 30-category field 187–195 ms → no task over 50 ms; hiding one category in the legend 319–332 → 38 ms; showing it again 232–237 → 42 ms.
- **Web app**: A 10M-cell tree landing no longer holds the main thread: the spatial-tree worker sends the tree as typed arrays (a structured clone of its 38,568 node objects took 35 ms in one task; now 10–11 ms with the install), and the store's adoption of the tree's draw order copies its source words on the GPU, drawn from the store into the source texture, instead of uploading 120 MB of positions and 40 MB of colours from the page (adoption 32–35 → 14–19 ms). Chromium (Metal), 10M cells: no task over 50 ms from opening to the tree in use (was a 77–90 ms frame of two tasks, once one 101 ms task).
- **Web app**: A highlight of one cell list costs that list: a page with one enabled group is painted (runs of consecutive cells as one fill) and published from the group's own list, counted as its length, and a cell list is proven ascending once (`app/utils/cell-set.js`), so the page selectors and analysis no longer walk it on every change; the highlight layer writes each record's color as one word, in a packing loop of its own. Chromium, 10,000,000 cells, highlighting 5,000,000: 114 → 47 ms, next frame 72 → 30 ms; all ten million: 177 → 72 ms, next frame 164 → 30 ms. A highlight of many cells is drawn from the view's point store, already on the GPU in draw order, with a colour word a slot written on the GPU from the highlight and visibility bytes, whenever that writes fewer bytes than records of its own (the CPU gathered ten million positions at random: 72 ms, and a new 160 MB store was cleared by the GPU process before it could be checked: 54 ms); the pulse frame of ten million highlighted cells takes 9 ms of GPU time instead of 54 (`tests/browser/highlight-store-draw-parity.spec.mjs`, three engines).
- **Web app**: Opening an analysis loads Plotly at the page's next idle moment, so its evaluation (one script, about 190 ms of main thread) runs while the analysis is set up instead of inside the first run: a differential-expression run started two seconds after its panel opened no longer has the 234 ms task (Chromium, 1.5M cells).
- **Web app**: A MatrixMarket folder opens in about one and a half times the time plain inflate of the same bytes takes: the entries are parsed in blocks of about 1 MiB on the decode workers into typed arrays and collected a block at a time, a large gzip text or R file is decompressed on a worker of its own (`formats/rds/gunzip-worker.js`; `DecompressionStream` inflates on the calling thread), the transcode space is written 1 MiB at a time, and the cell tables are read while the matrix streams. 100,000 cells (44M counts): Chromium 1.6 s (1.5× plain inflate), Firefox 1.7 s (1.5×), WebKit 1.0 s (1.33×); 500,000 cells in Chromium 7.6 s. `npm run bench:formats` runs the decode workers on `worker_threads` (Node, 100k cells: 1.80×). An entry whose row or column is not a whole number is refused.
- **Web app**: One shared module (`utils/exact-record.js`) holds the record and key checks the app's modules each defined for themselves (about forty copies, with different messages and different ideas of a record): `requirePlainRecord`/`requireKeys`/`requireExactKeys` for data from JSON, a session or a worker, `requireObjectRecord`/`requireObjectKeys` for an owner the app built, `requireClonedRecord` for worker payloads. A key mismatch reads the same everywhere: "… must contain exactly: a, b, c (missing: b; unexpected: d)." The community-annotation wire contract keeps its own path-addressed messages, the ones the annotation repository's validator writes.
- **Web app**: One formatter for sizes and durations (`utils/format-units.js`): the notification center wrote "2 MB" and "850ms" where the dataset card and read progress wrote "2.0 MB" and "850 ms". Every message now writes "512 B", "4.2 KB", "12 MB" and "850 ms", "4.2 s", "3 min 5 s", and no unit rounds into the next one ("1024 KB", "1000 ms", "60.0 s").
- **Web app**: A file whose bytes are not HDF5 names the format it was opened as: "The selected file is not a valid HDF5 file, so it cannot be a Loom (.loom) file. Write it again with loompy …"; a damaged file or Zarr store gives its own format's route to write it again.
- **Web app**: The viewer channel's modules and names are neutral: `data/viewer-host-source.js` (`ViewerHostDataSource`), `app/viewer-host-command-handler.js` and `app/viewer-host-pointer-delivery.js` (they were `jupyter-*`, `JupyterBridgeDataSource`, for a channel that serves Python, Jupyter and R hosts), and their messages say "Viewer host …". The `jupyter` source type of sessions and analytics, and the `jupyter` runtime, are unchanged.
- **Web app**: The "Link or server" ⓘ tip names Figshare, Google Drive and OneDrive, and every other link provider (a contract holds it to the provider registry). The cell-info tip and page state the spatial-index threshold the viewer uses (views of 250,000 cells or more, built in a worker); they said "above two million cells" and "above 256 MB of positions no background build".
- **Web app, development**: `npm run bench:formats` measures every non-AnnData format (Loom, a Cell Ranger `outs/` folder, a MatrixMarket folder, Seurat and SingleCellExperiment `.rds`, a Seurat `.rda`, `.h5mu`, MuData Zarr) against the same data as a Cellucid-ready `.h5ad`: first view, first gene, peak memory, and a stream format's pass against plain inflate of its own bytes (`scripts/benchmark-formats/`).
- **Web app, development**: A monthly CI job (`.github/workflows/r-corpus-drift.yml`, also on changes to the R readers; locally `sh scripts/rds-corpus/drift.sh <Rscript>`) writes the R object corpus with the current CRAN SeuratObject and Matrix and Bioconductor SingleCellExperiment, reads every file with the transcoders against what R read back, and fails on a slot those packages add or remove (`scripts/rds-corpus/known-slots.json`).
- **Python, R, web app**: An h5Seurat file is refused the same way everywhere, in the words of the input format registry (`refusal`, in all three repositories' `input_formats.json`): "This is an h5Seurat file (SeuratDisk), which Cellucid does not open. In R, read it with seu <- SeuratDisk::LoadH5Seurat(path) and save it with saveRDS(seu, "object.rds"), then open the .rds, or show it with cellucid::cellucid_view(seu)." `cellucid serve x.h5Seurat` and `cellucid_serve()` refuse it before serving.
- **Python**: Stopping a server closes the connections it still holds, so a viewer started next on the same port never answers from the previous dataset.
- **Web app**: A gzip or zlib Zarr chunk, and a ZIP entry, are decompressed with their input fed in 32 KiB slices, so the size the format declares bounds the output while the chunk decodes, in every engine (WebKit's `DecompressionStream` inflates everything written to it at once). The session codec and the format readers share this one input feeder (`formats/io/inflate.js`).
- **Web app**: A read's progress counts the bytes of every request as they arrive, so the progress, rate and time left of a long request (an 850 KB gene-name heap over a slow link, an R object or a MatrixMarket folder read once) stay live; the reader worker posts a running read's progress at most every 250 ms and always posts the last change of an interval, and the bytes of an attempt that fails are taken back. A GEO sample file its host sends compressed is downloaded with its bytes counted.
- **Web app**: Picking or dropping a file, or making a connection, is always a new opening, so the file already open opens again (a newer copy with the same dataset id is shown); only choosing a listed dataset again keeps what is open. An `.rds` opened again in the tab reads what was kept the first time.
- **Web app**: A format refusal's advice is what can be done from where the file was opened (the reader worker's `open` is told `fromLink`): a 10x Cell Ranger `.h5` with no embedding opened from a link says to download the whole `outs/` folder and open it from this computer, or to compute an embedding (a link opens one file, and no provider shares a folder Cellucid reads); a linked molecule table says to link `filtered_feature_bc_matrix.h5` instead; a 10x MatrixMarket matrix from a GEO sample, which has no `outs/` folder, says to compute an embedding; and an `outs/` folder whose `analysis/` has no projection says Cell Ranger ran without its secondary analysis.
- **Web app**: An address typed without `http://` on this computer (`127.0.0.1:8000/atlas/`), or a server's address, that answers but is no Cellucid server (a static web host serving a folder: no `/_cellucid/health`, or a web page there) is said to be no Cellucid server, with how to open such a folder: the address of its `datasets.json`.
- **Web app**: The startup phase that opens the dataset is published as `dataset` (`data-startup-phase`, and the `startup_phase` analytics field), not `identity`.
- **Web app, development**: `npm run bench:anndata` benchmarks the AnnData reader, and `scripts/benchmark-real-links.mjs` opens files from this computer through the File… and Folder… buttons (`{local}/path`), serves a link over HTTP/2 (`--link-protocol h2`), and reports each step's peak browser memory and main-thread long tasks.
- **Web app**: Opening a link (`?data=`, the address bar) or a server reads no sample file: the dataset list names the samples from `datasets.json` (one request) and reads their files (title and cell count, about 85 requests) only once the list is reached for. A page opened on a sample, or with nothing connected, lists them as before.
- **Web app**: The address bar names a link the page was opened with in its provider's form, as one opened from the bar is named: a pasted Google Drive link becomes `?data=gdrive:<id>`, a Figshare page `figshare:<id>`, a OneDrive share link `onedrive:<drive>/<item>`, a DOI of another publisher the Figshare item it names, and a Dropbox link its share link (`dropbox:https://www.dropbox.com/…`).

### Removed

- **Python**: The prepared export format and everything that wrote or served it: `prepare()`'s `latent_space`, `out_dir`, `obs_categorical_dtype`, `var_quantization`, `obs_continuous_quantization`, `compression`, `var_gene_name_columns`, `source_citation` (now `citation_text`/`citation_url`), `created_at` and `centroid_outlier_quantile` arguments; `Embedding` (store the whole array and choose columns with `embedding_defaults={"3d": ("X_pca", [0, 1, 2])}`); `generate_datasets_manifest()`; `serve()`, `CellucidServer`, `show()` and `CellucidViewer`; the prepared-export branch of `cellucid serve`; the export lock, transaction journal, manifests, quantization and payload helpers.
- **Python**: `AnnDataAdapter` and the server-side virtual export it computed (manifests, payload and mapping routes, server-side gzip, centroids at startup), the backed read mode, `serve_vector_fields` / `cellucid serve --vector-fields` (the viewer reads a vector field when its overlay is shown), `/_cellucid/health`'s `is_backed`, and `centroid_outlier_quantile` from `serve_anndata`, `show_anndata` and `cellucid serve`.
- **R**: `cellucid_prepare()`'s `latent_space` (name an embedding with `latent_key`), `out_dir`, `obs_categorical_dtype`, `var_quantization`, `obs_continuous_quantization`, `compression`, `var_gene_name_columns`, `source_citation`, `created_at` and `centroid_outlier_quantile` arguments, `cellucid_embedding()`, and the export lock and transaction journal.
- **Web app**: Prepared exports (`dataset_identity.json`, the obs, var and connectivity manifests, quantized and gzip payloads), catalogs of exports, `state-snapshots.json`, the **Prepared** button and the prepared-folder drop; with them gene download-ahead and its payload cache, the quantization worker, the data-mapping "Read from file" section and the analytics `data_version` parameter.
- **Web app**: `defaults['centroidOutlierQuantile']` (a file or server naming it is refused as an unknown default) and the per-centroid outlier array of views and kept-view snapshots; the marker cache's IndexedDB tier (`cellucid_marker_cache`) and its `maxAgeDays`/`cacheVersion` options; the bulk gene-expression cache and its API (`getGeneExpressionSubset`, `fetchBulkGeneExpression`); the data layer's session cache export and import.
- **Datasets**: `exports-v3/` in `cellucid-datasets` (299 MB of prepared exports and their snapshots) and the eight preparation notebooks; `exports/` and the Zarr v3 folder in `cellucid-demo-custom-datasets`.

- **Web app**: The per-file notices "Opened H5AD directly in the browser…" and "Opened a Zarr … directly in the browser…", the loaders' own download notifications and their `showProgress` option, and `ZarrDataSource.loadFromFileMap` (nothing called it): the opening progress shows every file.
- **Web app**: The analysis **Performance settings** (batch size, memory budget, network and compute parallelism, and their presets); Marker genes' mode select (Ranked, Clustered, Custom genes), its **Use cached results** checkbox (the cache is always on, marked **Cached** with **Recompute**) and its ranked and heatmap CSV files; Gene signature's sum; the run-cost line under the run button.

- **Python, web app**: The notebook `postMessage` command path (the script in the notebook output, `AnnDataViewer.viewer_origin`) and the routes `POST /_cellucid/events`, `/_cellucid/selection` and `/_cellucid/session_bundle`: commands, events, selections and session bundles travel on the viewer channel, and the server answers every `POST` with 405. The viewer's health poll every 3 s is gone: the host closing the channel is the stop signal.
- **Web app**: The remote source's unused update socket (`ws_port` in `/_cellucid/info`) and the connection-loss recovery that only it could trigger; a server advertising `ws_port` is refused.

### Fixed

- **Python**: The viewer's web build is checked once per file version: every asset request re-read the inventory and hashed the asset under one lock, 30 ms each (6 s for 200 modules before a first frame); now 0.2-0.75 ms.
- **Python**: An error rebuilt from a `DataMappingError` raised `TypeError` instead of the error.
- **Python**: Stopping a viewer ends every wait still pending on it (`get_session_bundle`, `wait_for_event`, `set_data_mapping`) at once with `RuntimeError` naming the stop; a wait without a timeout used to hang forever.
- **Python**: `get_session_bundle(timeout=None)` waits for the bundle as long as it takes: its request stayed answerable for an hour only, so a bundle that arrived later was refused while the call kept waiting.
- **Python**: `CellucidSessionBundle.save()` writes atomically and saving a bundle onto its own file leaves it as it is.
- **Python**: The viewer channel closes the RFC 6455 way, so the page reads the close code and reason it was sent. The server shut its socket down right after its close frame, with what the page was still sending (the rest of a text frame refused as over 1 MiB, its own close frame) unread, which resets the connection; on Windows the reset discards the close frame before the page reads it, so the page saw an abnormal closure (1006) instead of 1009, 1008 or 4401. Every close (a broken rule, a refused hello, a replaced or stopped viewer) now sends the close frame, shuts the write side, and reads and drops what the page still sends until it closes its side or 5 s pass (`CLOSE_DRAIN_SECONDS`), on the connection's own thread; a close made from another thread reaches a thread waiting to read within 0.25 s. Stopping the server closes an open channel with 1001 ("The server stopped"); it sent no close frame at all. A session bundle cut off by a dropped connection no longer leaves its partial temporary file behind.

- **Web app**: On a portrait screen (a phone held upright) the camera no longer cuts the data off at both sides. The 45° field of view applies to the plot's shorter side (`rendering/viewer/field-of-view.js`), so every camera — the default of 1D, 2D and 3D, a sample's preset, a session's — shows across a portrait screen what it shows down a landscape one; landscape screens draw exactly as before.
- **Web app**: A count and its noun agree everywhere ("1 cell", "1 row", "1 entry"; the legend read "1 cells"): every count with a noun goes through `formatCountOf` or `countNoun` in `utils/format-units.js`, with `formatCount` moved there from the sidebar's formatter so every layer uses one; numbers are grouped in English whatever the browser's locale (a refusal read "40.000.000 points" on a German or Turkish machine). `tests/count-noun-agreement-contract.test.mjs` refuses a `${n} cells` template, a hand-rolled `n === 1 ? '' : 's'` and an unlocalised `toLocaleString()`.
- **Python**: Messages, errors and warnings make a count and its noun agree ("1 row", "1 column", "1 cell"; they read "1 rows"): the noun comes from `cellucid._plural.plural`, and `tests/test_count_noun_agreement.py` refuses an f-string that writes a count before a plural noun and a hand-rolled `"" if n == 1 else "s"`.
- **Web app**: A sign-in line that links the administrator steps (a Google Workspace account its administrator blocks, a Microsoft organisation that lets only an administrator approve Cellucid) shows its reason on a row of its own, then the link, **Use another account** and the (i); beside them the reason was squeezed into a column one letter wide.
- **Web app**: Links in the light theme (dataset details, ⓘ tips, link buttons) are cyan-700, 5.9:1 on white, instead of cyan-600 at 4.1:1, below WCAG AA; the contrast contract now measures the link tokens.
- **Web app, Python, R**: Stored zeros in a sparse neighbour graph are not edges: an `obsp['connectivities']` (or a `prepare()` / `cellucid_prepare()` graph) with explicit zero entries refused the whole dataset; the zeros are dropped. Duplicate coordinates are still refused.
- **Web app, Python**: "Cellucid found nothing to plot" ends with how to fix it: "Add a 2-D embedding to obsm (for example with sc.tl.umap in Scanpy), then open the data again." A file with no cells says so ("This AnnData has no cells (n_obs is 0), so there is nothing to plot. …") instead of asking to try again, and a catalog dataset that cannot be read names its file and catalog entry.
- **Web app**: Opening another dataset with kept views or a gene shown, or after a session restore, works; it failed with "cannot replace … live points … while snapshot views remain published", and the half-done rollback left the old dataset open. A replaced dataset is released: after a switch from Suo the page's array buffers are back at the pre-open 9 MB.
- **Web app**: A deleted WebGL texture or program is never bound again. A smoke rebuild deleted a texture an owner later restored as its saved binding, raising INVALID_OPERATION (0x502), which a dataset switch then reported as its own failure; every texture and program is deleted through one tracked path, and restores drop deleted objects.
- **Web app**: No request to a web host is left hanging. Each request is made again after a transient failure, with growing waits or as long as `Retry-After` asks: HTTP 408, 429 or 5xx, or a failed or dropped connection. A 503 that GitHub Pages sends without CORS headers arrives as a failed connection. A request whose answer or body stops arriving is aborted and made again. It counts as stalled when neither it nor any other request of its link has received anything for TCP's retransmission timeout (RFC 6298), measured on the link. A failure that outlasts its repeats names the request's URL ("…/obsm/X_umap_3d/3.0: the host could not answer (HTTP 503). Asked 7 times."). Each repeat is logged as a `[cellucid:transfer]` warning naming the URL. A failure that will stay is reported at once: the same error status twice (Zenodo's 500 for a missing record), or no answer twice (a host that does not let this site read it). Asking 7 times took about 16 s. Catalog files and the CELLxGENE, GEO, Hugging Face and Zenodo APIs use the same rules. Cancelling still ends a request, and any wait to repeat it, at once.
- **Web app**: The startup error reads "Startup stalled while opening the dataset (its metadata, field lists and the embedding its first view draws): nothing arrived from the host of dataset "suo" for 60 s, while reading …". It names the reads still open and the last failed request with its URL. Opening a dataset at startup now stalls only when nothing at all arrives for 60 s. Before, any open that took 60 s failed, even one still receiving bytes over a slow link.
- **Web app**: A dataset opens on the embedding its first view draws. The other dimensions named in `uns['cellucid']['defaults']` are read when first shown, not one after another before anything is drawn. The opening embedding is read beside the field lists. A consolidated Zarr v2 store on a web host opens with one request (`.zmetadata`; `zarr.json` is asked for only when it is absent), so no request at open answers 404. The sample list reads the other catalog datasets only after the dataset being opened has opened, and lists that one from its opening, without reading its store twice. Suo from the sample catalog over HTTP/2, first render before → now (median of three, WebKit / Chromium / Firefox): at 30 ms and 5 MB/s, 3.3 / 3.1 / 3.5 s → 1.9 / 1.9 / 2.2 s; at 80 ms and 2 MB/s, 6.9 / 6.8 / 7.2 s → 3.8 / 3.7 / 4.1 s. Before first render: 139 → 27 requests, 10.9 → 5.5 MB, six 404s → none. The 0.9.4 app with the prepared export took 2.1 / 2.0 / 2.4 s and 4.3 / 4.2 / 4.5 s (11 requests, 6.3 MB).
- **Web app**: A gene-column index that a closed or crashed tab spilled to the origin-private file system is removed when the app starts, not only when the next index is written.
- **Web app**: A file opened after another dataset opens on its own view. Its camera is the default framing of the dimension it opens in; it kept the pose the previous dataset was left in, so a 28,000-cell 2D file opened after Suo (whose preset puts the camera inside its 3D cloud) filled the canvas with giant points. The point size follows its cell count, as before; a session or published preset still sets both afterwards.
- **Web app**: A file without a default UMAP opens on the default whose name ranks best (UMAP, then t-SNE, force layouts, PHATE, PaCMAP/TriMap, spatial, diffmap, any other name, then PCA and latent spaces), the higher dimension among equals, instead of always on its highest dimension. A segmentation export with 2-column `X_umap_scale*_latents` and a 3-column t-SNE or spatial array opened in 3D on the 3-column array; it opens on the UMAP in 2D, with the 3D default still a key press away. A file whose defaults are all UMAPs opens in its highest dimension as before.
- **Web app**: In Data mapping, an `obsm` array wider than three columns (`X_pca`, a latent space thousands of columns wide, a 4-column `coordinates`) is offered at its first two columns in 2D and listed after the embeddings; it was offered as a 3D embedding at its first three columns beside the file's real 2D UMAP.
- **Web app**: A dataset chosen while another was still opening stops the replaced opening at once (a catalog dataset's read, or a file's read and its reader worker). Its progress stayed on screen until the stopped read answered, the newer opening's outcome (*User data ready: …*) could be lost, and a sample or file chosen while the startup dataset loaded was announced under the replaced dataset's name.
- **Web app**: The `Full`-quality highlight pulse no longer holds the main thread on a slow renderer. A frame drawn only for the pulse now waits until what the frames before it held the main thread is measured (from the first of them to the first task after them), and then until the page has been free for as long (`rendering/ambient-frame-pacing.js`), so pulse frames take at most half of the main thread; where a frame costs less than a display interval nothing changes. A frame anything else asks for is drawn at once. While the measurement is owed the viewer's loop stays unarmed, so the measuring task always runs, also where the browser runs the next animation frame before queued tasks (Firefox, and WebKit at times). `tests/ambient-frame-pacing-contract.test.mjs` proves the rule on that event loop with 450 ms frames, in both task orders and under a stream of worker replies. Headless Chromium with 2 analysis workers: a 400-gene marker run with highlighted pages went from 7.0 s to 0.5 s (9.7 s to 2.0 s at 4× CPU throttling); on the CI runner the same run had not finished in 30 s.
- **Web app**: A category no cell holds has no centroid. The centroid pass returned one for every category in a field's list, at the origin with `n_points: 0`, and the viewer drew a point and built a label for each: a field listing 70,000 categories of which three hold cells built 70,000 label elements per view (140,000 after restoring a session with a kept view), stacked at (0, 0, 0), and the stats line read "Centroids: 70000". The pass now returns the categories with cells, in category order; its consumers find a centroid's category by identity, never by its place in the list (the visibility pass keyed categories by their text, so a `1` and a `"1"` shared one count). A consensus column rebuilt while it is shown is drawn again from its new codes, colours and centroids, so a category that gains its first cell gets its centroid at once. Colouring by the 70,000-category browser fixture took 0.74 / 1.26 / 1.21 s (Chromium / WebKit / Firefox) and now 0.31 / 0.69 / 0.59 s; restoring its session 2.2 / 3.0 / 3.3 s, now 0.65 / 1.3 / 1.0 s; at 6× CPU throttling in Chromium, 4.3 → 1.8 s and 13.0 → 3.5 s.
- **Web app**: Every view's copy of a field you made is the same field. Deleting or merging one of its categories, or rebuilding a consensus column, replaced the categories, codes and centroid tables of the edited view's copy and of the stored field only. A kept view then switched to an embedding not yet loaded failed with "Normalized user-defined centroids of layout "X_umap_2d" were not published", went on drawing the old categories, and an edit made in a kept view gave the live view the new categories under colours meant for the old ones. Each copy now takes the new definition at once and keeps its own colours and hidden categories, by category; a view that shows the field, or hides some of its categories, is drawn again. A field you made that is deleted, restored or confirmed deleted in one view is so in every view once that view is selected, as a dataset field already was, and Confirm releases its codes from every copy (from the Deleted fields list it only set flags and kept them).

- **Web app**: A data mapping applied while another dataset is opened belongs to the dataset it was asked of. Data mapping's **Apply** reads its neighbour graph, gene list or embedding over the network before it publishes anything; a dataset chosen meanwhile received the old dataset's graph, gene list and mapping label, its refusals and, on failure, the old dataset's state restored into it. The transaction now stops at the first await after another dataset is published (an `AbortError`), publishes and restores nothing, and the panel shows nothing of it.
- **Web app**: An outlier filter being computed stops when another dataset is published. Its read of the old latent space ran on (hundreds of MB at 2M cells), then announced "Outlier filter … ready" for the old dataset, and a latent space it could not use was marked unusable in the new one.
- **Web app**: An analysis releases what it left in the analysis workers however it ends (finished, cancelled, superseded or failed): a marker run's groups no longer stay in every worker or get replayed into replaced workers, a run cancelled during its setup leaves no groups in the workers that took them, and a UCell signature's per-cell ranks (about 18 bytes a cell per set gene) go with it.
- **Web app**: An analysis scoped to **Visible cells** reads the filters: a filter change marks its result **Stale · Run** and reads its groups again; it changed neither. A highlight change also marks stale the result of an analysis mode that is not shown (only the shown mode heard it), which catches up when shown.
- **Web app**: Gene search lists the gene whose name is exactly what was typed first, then names starting with it, then the rest (in the sidebar and in the gene-set editor), so `TNF` + Enter colours by TNF, not by whichever of `C1QTNF1` or `TNFAIP3` came first in the file.
- **Web app**: The notebook's `clear_highlights()` (Python) and `cellucid_clear_highlights()` (R) clear every highlight page, as documented; they cleared only the active one. The Clear button still clears the active page (`clearActivePageHighlights`).
- **Web app**: Main-thread work per interaction at scale. Marker genes' display options (direction, filters, ranking, top N) pick each group's top genes with a bounded heap instead of sorting every gene of every group: 100 groups of 30,000 genes, 716 → 36 ms per change. Differential expression resolves its two groups (a pass over every cell, 8 ms at 2M cells) once per change of the options that choose them, not several times on every option, highlight and page change. The highlight page tabs count a page again only when its enabled cell lists change (each highlight or page event, and every input event of a colour drag, recounted every page); a colour drag publishes once per frame. Mean z-score signatures reuse gene moments a marker run already computed and yield between genes. Marker genes no longer reread its groups on a highlight change when **Cells** is All cells.
- **Web app**: No GPU-pipeline drain on a frame path. In Firefox and WebKit, `isEnabled` and `getParameter` of the active texture unit, viewport, scissor box, colour mask, blend functions and clear colour each wait for every queued GPU command (6–8 ms in Firefox, 16–24 ms in WebKit after one 64 MB upload; bindings and pixel-store parameters answer at once in every engine). The backdrop look asked three of them per pane per frame, the GPU gather one per draw-order change, the smoke density, noise and edge uploads several per build. Every renderer now reads that state from a page-side shadow of the context (`rendering/gl-state-shadow.js`), kept by the context's own setters.
- **Web app**: The viewer redraws less. A smoke volume waiting for its noise textures no longer redraws the scene every frame until they land (the generation wakes the viewer once); reads, listener and handler registrations, `withBorrowedViewData` (every figure-export preview build), `stopInertia` and the navigation preferences no longer mark a frame, which redrew every pane of a grid; a live filter publication that moved no cell keeps the vector-field particles and their trails.
- **Web app**: A background spatial-tree build whose positions were replaced (another dataset, representation or a removed kept view) is dropped with its worker: it ran to the end first, so the new data drew without culling or LOD meanwhile, and kept the old positions alive. A replaced dataset also releases its smoke volume (its 3D textures and the copies of its colours and visibility) at once.
- **Web app**: Keyboard focus survives a rebuild. Toggling or removing a filter, or a highlight group, keeps the focus on that row's control (on the row now in its place after a removal, on the list after the last one); each row's remove button is named for its item ("Remove filter: cell_type: 2 of 3"). A long categorical legend keeps the rows still in view when it scrolls (it rebuilt every row, dropping the focus and an open rename). The phone sheet tabs take arrow keys, Home and End, and keep the focus on the tab chosen.
- **Web app**: The orbit anchor is shaded for the pane's real camera: its eye used another spherical convention, which at the horizon put the camera overhead, so its specular highlights and rim light faced nowhere.
- **Web app**: Smaller fixes: a failed group read in Marker genes or Differential expression says why on Run instead of reading again and posting its error on every redraw; the retired-notification log appends a row instead of rebuilding the list; the drop overlay names no format list that can fall behind the registry; analysis plot names and options are in sentence case.

- **Web app**: Messages name the sidebar's panels and blocks as their headings do ("Turn "Show overlay" off under Vector field", not "under Vector Field Overlay"; "in the Highlighting section", not "Highlighted Cells"): they take the names from one module (`app/ui/place-names.js`), and a contract test holds each name to its heading in `index.html` and fails on a message that names a place in words of its own.
- **Web app**: Analysis plots lay their labels out for the size they are drawn at, and again when resized, so none lands on another at any panel width. The correlation's per-page R² boxes were stacked 15% of the plot apart, less than a box is tall: they stand in a column, else a row, else one compact box beside the legend (which gives way in a narrow preview), and each colour category is listed once across pages, in one colour. The dot plot's colour bar and size legend stand in a column right of the plot, placed in pixels: the size legend's title ran into the colour bar in the sidebar preview and its dots into each other in the expanded view. Its gene names shrink and turn to fit their columns, and a panel too narrow to name every gene scrolls the preview across. A legend above a plot gets the rows it needs, or is left out of a preview too short for it.
- **Web app, Python, R**: An embedding under the R world's conventional UMAP name is the default UMAP, as Scanpy's `X_umap` is: after `X_umap_1d`/`X_umap_2d`/`X_umap_3d`, the first 1- to 3-column one of `X_umap` (also a Seurat object's `umap` reduction), `X_UMAP` (a SingleCellExperiment's `UMAP` reducedDim, as scater's `runUMAP()` names it), `X_wnn.umap` and `X_ref.umap` (a Seurat object's weighted-nearest-neighbour and Azimuth reference UMAPs) fills the slot of its width. The rule is pinned in the shared data-mapping conformance fixture (version 6, SHA `cfc84c8c…`, in all three repositories; web `defaultEmbeddingSlots`/`UMAP_KEYS`, Python `cellucid.data_mapping.default_embedding_slots`/`UMAP_KEYS`, R `.shown_embedding_keys`/`.umap_keys`).
- **Web app**: A sample whose file cannot be read while the list is made (a busy host, a file taken down) no longer empties the dataset picker with "The dataset list could not be loaded": the open dataset stays shown, the notice names the sample that failed with Try again, and a dataset that opens does not clear that notice.
- **Web app, development**: The two design-token validators (`scripts/validate-tokens.js`, `scripts/validate-token-types.js`) run with the contract suite; `types/design-tokens.d.ts` had drifted by sixteen names. Dead stylesheet rules for removed controls, two unused sidebar-bound tokens that disagreed with the real clamp, a duplicate `.gene-chip` rule and raw opacities are gone.

- **Web app**: On a phone the open full-screen menu no longer shows a band of the plot above it: the strip its ✕ sits in is painted as part of the menu (it was transparent). The ☰/✕ toggle no longer stays grey and pressed after a tap on a touch screen; its hover look is kept for a pointer that hovers.

- **Web app**: After a deploy, a browser cache no longer hands the page modules of two builds. A page could run a module its browser had cached from the build before a deploy beside one fetched from the build after it (GitHub Pages caches every file for ten minutes, a static server lets the browser guess a freshness of its own), and the app stopped on "The requested module '../../../data/data-mapping.js' does not provide an export named 'UMAP_KEYS'" with "Starting Cellucid…" shown forever. Every script is now requested at a URL that names its bytes (`<path>?v=<first 12 hex digits of its SHA-256>`): `index.html` (and the sign-in callback page) carries an import map of every module under `assets/`, the `src` of each script tag says the same, and the Content-Security-Policy allows the map by hash. `scripts/generate-web-asset-inventory.mjs` writes all of it, and the JSON-LD hash, beside the inventory; nothing needs a build step to run. A file that did not change keeps its URL and stays cached. A startup module that fails to parse, link or run now shows the "Cellucid could not load" card with **Reload** instead of waiting forever, and when the server already serves a newer build than the page's the card says "Cellucid has been updated". A page left open across a deploy that then loads a part on demand (the benchmark panel, figure export renderers, Plotly) and fails says so in a notification with **Reload**. Workers, which get no import map in any browser, are covered by the next entry.
- **Web app**: A deployed page starts every worker from one bundle named by its bytes, so a worker can no longer pair its own modules from two builds. `scripts/generate-web-asset-inventory.mjs` now builds each worker (AnnData reader, decode, gunzip, spatial tree, point raster, figure export, analysis, benchmark, and the two modules decode workers import by URL) with esbuild 0.28.2, pinned as a dev dependency, into `assets/workers/<name>.<first 12 hex digits of its SHA-256>.js` with a source map to the sources, and writes `assets/workers/scripts.js` naming them; they are committed, listed in the inventory (the Python and R web caches download and verify them) and `npm test` refuses a stale one. On `localhost`, `127.0.0.1` and `[::1]` (development, the browser tests, the Python and R servers) workers still start from their sources, with nothing to build. A worker that fails to start on a page that outlived its build shows **Cellucid has been updated** with **Reload**.

- **Web app**: Opening a file no longer hangs in WebKit. The inflation of a compressed chunk on the reader thread (`formats/io/inflate.js`) released its stream reader's lock, and WebKit (314030@main, also landed on a Safari release branch) runs `ReadableStreamDefaultReader.releaseLock()`'s JavaScript while it holds a lock its garbage collector takes when it visits the reader: a collection that started inside the release never finished, and every later collection in the tab, the page's included, waited behind it. The file stayed "Opening …" for good and closing the page hung (sampled native stacks: the reader thread parked in that lock inside `releaseLock`, the main thread in `JITWorklist::suspendAllThreads`); under CPU load it hung within 1 to 5 opens of the 70,000-category fixture, and the browser suite's WebKit shard failed on it. No browser code releases a stream reader any more (`tests/browser-stream-reader-release-contract.test.mjs`): a read-to-the-end or cancelled reader is dropped, and 100 opens in a row under the same load all opened. Closing the app also stops the dataset change still opening and the file or connection it reads (`retire()` of the dataset-change coordinator) instead of waiting for them to answer, and an opening that makes no progress for 60 s, the startup stall interval, says which step it waits on and the reads still open with what has arrived of each, in its notification and in the console (`[cellucid:opening] Opening a.h5ad stalled while reading the file's structure: nothing has changed for 1 min 0 s, while reading the structure of a.h5ad (0 B of 512 KB arrived).`), without stopping it; while a startup step's own deadline runs, that deadline reports the stall instead, once.
- **Web app**: The genes an analysis computes in its workers (marker genes, differential expression, the gene-moments pass of a module score, the dot plot, UCell) no longer pass through the page. Each gene went from the AnnData reader worker to the page and from the page to an analysis worker, and WebKit counts every array a worker sends the page as the page's own allocation: a run over every gene drove the page's garbage collector, whole-heap collections included, while the analysis ran, each collection holding the page's thread. A run now opens a gene route, one `MessageChannel` per analysis worker, and the reader delivers each task's genes straight to the worker that computes them (`data/gene-delivery.js`, `app/analysis/data/gene-route.js`); the page sends which genes each task takes, hears their bytes and takes in each task's statistics, packed. A task takes several genes (sixteen tasks per worker over the run), one read serves a task on each worker, and a worker the pool replaces during a run has its genes read again. Measured on the 20,000-cell, 200-gene markers run: 13 MB of gene values in and 13 MB out of the page per run became none, the page's allocation per run halved, and 40 reader requests became 13 (8 workers).

## [0.9.4] - 2026-10-21

Version 0.9.4 is a patch release on 0.9.3. The export format stays 3, so
exports written by 0.9.3 open unchanged; `prepare(obs_categorical_dtype="uint32")`
is the one addition to what an export may hold. The web app, the Python
package and the R package now share one version and are released together.
The main changes are:

- **Six new backgrounds** (Atelier, Nebula, Abyss light and dark, Aurora light
  and dark), prepared in the background and sharpened while the GPU keeps up.
- **No fabricated limits**: any number of kept views, selections of any size,
  direct H5AD and Zarr files of any size, categorical fields with more than
  65,535 categories, and session files bounded only by their format.
- **Analysis on the workers**: summaries, tests, histograms, violins and
  density plans leave the main thread; gene downloads go through one scheduler
  with a payload cache; marker-gene clustering is bounded by memory, not a
  count.
- **Jupyter selections** of any size reach Python as binary indices.
- **Sessions** restore the vector-field controls after the embeddings they
  depend on, and keep every view's embedding.

### Added

- **Web app, version**: the web app declares the release version it shares with the Python and R packages: `version:` in its `CITATION.cff`, a `cellucid-version` meta element and `softwareVersion` in the page's structured data, and the footer prints it before the build id. The three are released together under one version.
- **Web app**: Six new backgrounds, in light and dark pairs: **Atelier (light)** and **Nebula (dark)**, **Abyss (light)** and **Abyss (dark)**, **Aurora (light)** and **Aurora (dark)**, listed under the Light and Dark groups of the **Background** menu; **Studio (light)** stays the default. Each is a painted stage that turns with the camera, with a nearer layer that slides against it for depth: Atelier's paper studio, where the cells cast a soft contact shadow that follows filters (and glows through the paper seen from below); Nebula's deep field, tinted by the colours of the field the cells are coloured by; Abyss's sunlit shallows and deep ocean with glowing marine snow; Aurora's dawn and night skies. A look is prepared in the background in idle time the first time it is chosen (its sky coarse first, then sharpened level by level while the GPU is measured to keep up, so it appears within a fraction of the full work on every machine and reaches 1024² per face on a capable one), costs nothing while the scene is still, keeps its large areas outside the lightness range of the category palettes, and appears without its fade under the system's reduce-motion setting. Only the look on screen keeps GPU memory (at most about 27 MB). Sessions and dataset presets can name the looks, and figure export redraws the look under the points in PNG and SVG (as an image) when the figure's background is Match viewer. (backgrounds)
- **Web app**: Keep any number of views: the four-view ceiling is gone. GPU memory decides instead: when an allocation is refused, the colour and visibility stores of kept views that are not on screen are released first and rebuilt, unchanged, when shown again. The benchmark harness's view-count axis takes any number too.
- **Web app**: Selections of any size: the annotation, lasso, KNN and proximity tools, category highlights and Jupyter commands hold selected cells as one bit a cell instead of a JavaScript `Set`, which Chromium refuses above 16,777,216 entries, so selecting a large cluster in a 30M+ cell dataset failed. Combining selections (union, intersection, difference) works a word at a time.
- **Web app**: Categorical fields with more than 65,535 categories (barcodes, free text) load: they use 32-bit codes (export format `uint32`, missing marker 4,294,967,295) in prepared exports, `cellucid serve`, Jupyter, direct H5AD/Zarr and saved sessions. A direct H5AD/Zarr or served AnnData with such a column used to be refused as a whole. Existing exports are unchanged. (limits)
- **Web app**: Direct H5AD files have no size ceiling: files above 512 MiB used to be refused before they were read. The file is read from disk as each view needs it instead of being copied into memory, in a background worker, so opening it no longer blocks the interface (a 250 MB, 1-million-cell CSR file: about 0.8 s of main-thread loader work and 1.1 GB peak before, about 0 s and 0.8 GB now; it opens in about 0.25 s instead of 0.7 s). A row-major (CSR) matrix is reorganised by gene in the background right after opening, with progress shown, and a gene asked for meanwhile is answered from the same pass; when the reorganised matrix does not fit the memory budget it is kept on disk (the browser's private file storage). The fixed working-set ceilings on sparse matrices, dense matrices and observation columns are removed. (limits)
- **Web app**: Zarr directories and Zarr ZIP archives load through the same background reader, with several chunks decoded at once, and without the fixed per-array and per-archive size ceilings; every chunk is still checked against the exact size its shape and dtype declare. (limits)
- **Web app**: Memory budgets for kept layouts, inactive fields, the AnnData cache, decoded blocks and the gene index scale with the memory the browser reports (Chromium) instead of fixed values, and they decide only what is kept for reuse: no dataset is refused for them any more. Browsers that report nothing keep the previous values. (limits)
- **Web app**: Latent-space outlier quantiles of a direct H5AD/Zarr file work for a latent space of any size: the 512 MiB latent-space budget is removed. The latent space is streamed from the file in row slabs in the reader worker, with the same bit-identical rule, and kept in memory only while it fits the reader's budget. (limits)
- **Web app**: Categorical fields made from highlight pages work at any cell count (they were refused above 16,777,216 cells, a JavaScript `Set` limit), an intersection field combines up to 32 pages (was 12), and the limits of 20 user-defined fields and of 256 characters for a field name or category label are removed. (limits)
- **Web app**: Session files have no 16 MiB manifest or 512 MiB chunk limit; the only limit is the format's 4 GiB − 1 per frame (its length is a u32), and every inflate is still bounded by the chunk's declared size. (limits)
- **Web app**: Prepared exports have no 512 MiB size ceiling per file: a 3D dataset above 44.7 million cells was refused. Every payload is checked against the exact size its manifest declares, and one the browser cannot allocate is reported as too large, naming the bytes. Compressed payloads above 4 GiB load too (their gzip trailer is compared modulo 2^32).
- **Web app**: A sample dataset's published starting view has no size limit (it was refused above 32 KiB, or 64 KiB decompressed): the catalog now declares the state file's exact length in `state_bytes` beside `state_sha256` (`exports-v3/datasets.json`, and `cellucid serve`'s dataset listing), and the viewer reads exactly that many bytes and checks the digest. A catalog entry advertising a state must carry all three fields. (limits)
- **Web app**: Camera Path: no limit of 100 keyframes, 60 seconds per transition or 40 characters per label; a transition still lasts at least 0.1 s. (limits)
- **Web app**: Highlight selection undo reaches every earlier step: the oldest steps go only when every tool's recorded selections together outgrow a share of the device's memory (was 5 steps per tool). (limits)
- **Web app**: The cell info card lists every loaded categorical field (was 4), adding columns when they do not fit the window's height. (limits)
- **Web app**: The gene search list lists every match (built 100 at a time as you scroll or arrow; Home and End jump to the ends) instead of the first 100 with "…and N more". (limits)
- **Web app**: Notifications retired early to keep five cards on screen are kept: a "Show N earlier messages" control under the stack lists every one of them. (limits)
- **Web app**: Community annotation: a Pull has no limit of 10,000 user files or 64 MB in total. Changed files are fetched and cached about 32 MiB at a time, and the merged view is compiled by streaming the cache one file at a time; an interrupted Pull keeps what it cached. CAP lookup results are kept whole (results with more than 100 synonyms, 200 markers or 512-character text were dropped), every fetched result is shown (was 20 of 25), and synonyms, markers and author handles are shown in full. (limits)
- **Web app**: Analysis: the gene list lists every matching gene (built a page at a time while scrolling or arrowing) instead of the first 50; analysis windows and saved plot options have no count limit. (analysis, limits)
- **Web app**: Marker genes: clustering is no longer limited to 500 genes or groups. It is admitted by the memory it needs (about 100 bytes per pair of genes: about 100 MB for 1,000 genes, 1.6 GB for 4,000) against the analysis memory budget. A request past the budget is refused with that size and what to change, instead of "maximum is 500". (analysis, limits)
- **Web app**: Jupyter: selections of any size reach Python (a million cells arrive in about 120–200 ms). A confirmed selection is posted to the notebook's server as binary cell indices (`POST /_cellucid/selection`, little-endian `uint32`, ascending) instead of a JSON event sent with `keepalive`, which the browser refused above 64 KiB (about 8,000 cells). Validating and encoding one million cells takes a few milliseconds instead of 150–300 ms, and no step builds a per-cell `Set`, so selections above 16,777,216 cells are no longer refused. No event is sent with `keepalive` any more. Works with the matching `cellucid` Python version, released next. (DX-17)
- **Web app**: Gene downloads (analysis and the sidebar) go through one download scheduler: a gene the user chose starts at once and no new download of a long run (DE, markers, heatmap) or of a prefetch starts until it has arrived; the browser gets the same order through the fetch `priority` hint, and a prefetched gene the user then chooses is raised, not downloaded twice. Observation-field loads count as foreground too. A finished or aborted run cancels its queued and in-flight downloads, and replacing the dataset ends every download of the old one, so a stalled read can never hold back the next dataset. (analysis, network)
- **Web app**: A dataset load keeps its downloaded gene payloads (compressed, up to 64 MiB) and decodes an evicted gene again from them: a second DE or marker run over the same genes makes no requests (pancreas over 4G: 3,753 revalidation round trips before). Long runs download the payloads of the genes after their decode window ahead, as many at once as the connection needs (the window starts at 8 and doubles while throughput grows, never more than half the cache), so a run is bound by bandwidth instead of round trips. (analysis, network)
- **Web app**: A gene the user waits for shows its download progress in its loading notification. (analysis, network)
- **Web app**: The documentation moved to https://docs.cellucid.com (the old Read the Docs addresses redirect there); the footer, Cite, info-tip and welcome links point at it. Search and sharing metadata: a description that names UMAP, 3D, Scanpy/AnnData and Seurat; `og:image:alt`, `twitter:image:alt` and `og:site_name`; corrected JSON-LD (`WebApplication`, operating system "Web browser", licence, documentation, `sameAs` for GitHub, PyPI and the R package), still allowed by its CSP hash; a crawlable "About Cellucid" text with links to the docs, GitHub and citing; a `sitemap.xml` whose `<lastmod>` is the build date, written by `scripts/generate-web-asset-inventory.mjs`; a custom `404.html`; an `llms.txt`; the manifest name "Cellucid"; and a preload of the interface font.
- **Web app**: Data mapping, Embeddings: one compact row per array -- checkbox, name and dimension; an embedding that is on adds a second line with **Default** and **guessed**, one that is off stays one line. A long name wraps to a second line instead of losing its columns to an ellipsis. The dimension of an AnnData array is a chip with a chevron (`3D ⌄`) that opens the row's editor, closed until opened: **Dimension**, one labelled column menu per axis on one line, and, for an array wider than three columns, **Add columns** (it replaces the "+ Columns of …" pill under every wide array). Nothing in the panel overflows the sidebar on a desktop or a phone. (data mapping)
- **Web app**: Data mapping: every `obsm` array of two or more columns can be shown in any dimension its columns allow, at any of its columns -- a 2-column UMAP in 1D or 2D, a 3-column array in 1D, 2D or 3D, as a wide array already could. The guessed default is unchanged; views, sessions, highlights and figure export take the new selections like any other (`UMAP 1`, `Spatial 1,3`). (data mapping)
- **Python**: Jupyter: `POST /_cellucid/selection` receives selections as binary cell indices (little-endian `uint32`, strictly ascending, with no size limit; a body the server cannot hold is refused with HTTP 413), authenticating before it reads the body; `@viewer.on_selection` still receives `cells` as a `list[int]`, now always ascending. `POST /_cellucid/events` refuses `selection` and `session_bundle` records (HTTP 400): only their own routes produce them. Selections were previously capped by the event endpoint's 1 MiB limit. (DX-17)
- **Python**: `prepare(obs_categorical_dtype="uint32")` writes 32-bit categorical codes (`obs/{index}.codes.u32`, missing marker 4,294,967,295) for fields with more than 65,535 categories, inside export format 3; existing exports are unchanged and an older web app refuses such an export by name. The width error now names the width that fits ("Pass obs_categorical_dtype='uint32'."). `serve()` and `show_anndata()` serve a string or boolean obs column with more than 65,535 distinct values with uint32 codes instead of refusing the whole AnnData, and session bundles carry uint32 user-defined codes (encodings 4 and 5). (limits)
- **Python**: `cellucid serve` lists a published default session with its exact length (`state_bytes`) beside `state_manifest` and `state_sha256`, and reads a bundle of any size (the 32 KiB limit is removed), hashing it as it streams. `debug_connection()` accepts the three fields together. (limits)
- **Python**: Documentation: served from https://docs.cellucid.com, the Read the Docs canonical domain; every link in the package, README, CITATION.cff and docs uses it. Every page declares `<link rel="canonical">` on that host (`html_baseurl`), `sphinx-sitemap` (pinned in the `docs` extra) writes https://docs.cellucid.com/en/latest/sitemap.xml, and a `robots.txt` at the root of the default version advertises it.

### Changed

- **Web app**: When several views' spatial trees land at once (four views of four representations), each view switches to its tree in its own frame instead of all in one long frame.
- **Web app**: Slider, legend and camera-path changes are drawn in the frame they happen in, not the next one: they used to wait for an animation frame that ran after the viewer's, which added one frame of latency to every slider step (about 50–90 ms at 10 million cells in four views).
- **Web app**: Starting a drag or a wheel zoom in another view no longer waits for the sidebar: the camera follows the view at once, and the sidebar switches to it (legend, counts, controls) when the gesture ends (on release, or once the wheel has rested for 150 ms). A wheel turn that cannot zoom further, or a wheel in free-fly mode, no longer redraws the scene.
- **Web app**: The full-window selection overlay (lasso, proximity and neighbour rings, the projectile charge ring) is hidden while nothing is drawn on it, so the browser no longer composites an empty window-sized layer over the plot every frame.
- **Web app**: Category counts in the legend are recounted only while the legend can be seen: with its panel collapsed or the sidebar hidden, a filter change no longer recounts every category; the counts update when the legend is shown.
- **Web app**: A camera drag captures the mouse until it is released, so dragging across the sidebar or an analysis plot no longer makes the page work out hover styles and plot tooltips under the pointer at every move.
- **Web app**: Shader programs are compiled in batches: the point programs compile 9 shaders instead of 14 (shared vertex shaders) and the velocity overlay 8 instead of 14, with every link issued before any status is read, so the driver overlaps the work and the page waits once per batch instead of once per compile and link. Turning velocity on and starting the app pause less.
- **Web app**: Multiview: a view whose camera and contents did not change keeps its pixels instead of being redrawn, so dragging, zooming or flinging the camera of one view of four redraws only that view, and the settle frame after a change in a single view draws nothing.
- **Web app**: Zooming through levels of detail no longer rebuilds the figure-export preview at every step; it rebuilds once the camera settles.
- **Web app**: An enabled velocity overlay with nothing to draw (no field for the layout, nothing visible, particles off when zoomed out, smoke mode) no longer keeps the scene redrawing every frame.
- **Web app**: Filtering no longer blocks on the GPU: writing alpha, colours and draw orders into storage the GPU already holds no longer waits for the previous frame to finish (a filter tick at 10M cells in four views spent about 40 ms waiting).
- **Web app**: Framed (cropped) figure exports are drawn at their full resolution: a frame's points are rasterized at the pixels they land on instead of being enlarged from the whole view, so a framed PNG, or the point layer of a Hybrid SVG, is as sharp as an unframed one.
- **Web app**: When a Remote server connection drops, the viewer says so at once, even while another dataset you chose is still loading; the notice is completed with "Reconnect to continue" once recovery has run.
- **Web app**: Filtering a kept view uploads only the rows of its visibility texture that changed, reuses its own transparency array instead of allocating a new one per update, and skips its highlight, edge and overlay work when nothing moved.
- **Web app**: Recolouring a kept view rewrites its colour store in place instead of allocating a new store and vertex array each time (about 70 MB transient a view at 10M cells).
- **Web app**: A kept view's colour store is RGBA8 on a 4-byte stride, like the live one: Metal (Chrome and Safari on macOS) needs that alignment, and a 3-byte RGB store made the browser keep and re-convert a hidden 4-byte copy.
- **Web app**: Figure export of the Studio background no longer has a ceiling on the dots it considers: a flat embedding's dot surface seen at a grazing angle is enumerated only where it shows, and the exported marks are unchanged.
- **Web app**: An idle page no longer wakes 60–120 times a second: the frame loop stops once the scene has settled and starts again when something asks for a frame (an interaction, a data or setting change, a resize, a held navigation key).
- **Web app**: The KNN selection tool builds its neighbour graph as compact index ranges instead of a JavaScript map per cell, so it works on graphs with more than 16.7 million cells and builds much faster.
- **Web app**: Spatial trees (level of detail, frustum culling, hover) are built in the background worker at every dataset size; above 512 MiB of positions (44.7 million cells in 3D) they were built on the main thread, freezing the page for about 15 s at 50 million cells. Positions go to the worker in 64 MiB slices, so the page never holds a second full copy.
- **Web app**: Analysis keeps the interface responsive on large pages: the Detailed and Gene signature summaries and tests (from 200,000 values), the Quick Insights median and quartiles and the signature-score median (from 500,000 finite values), histogram bins (from 200,000 values) and density-plot plans are computed on the analysis workers, with the same functions and the same numbers as before. Categorical counts are taken once while a page is built and shared by every consumer; gathering a page and aligning two variables for a scatter plot yield to the interface in 12 ms slices. (analysis)
- **Web app**: Gene search, in the sidebar and in analysis: a gene kept under the pointer or highlighted with the arrow keys for 150 ms is loaded in the background, so choosing it usually draws at once; moving on cancels the load and choosing the gene joins it. (analysis)
- **Web app**: Violin plots of large pages no longer freeze the interface: the box statistics Plotly's violin calc needs are computed on the analysis workers with the density, and the vendored Plotly takes that calcdata entry instead of building a point record per value. A ten-million-cell violin held the main thread for 1.5 s and 800 MB; it now holds it for under 50 ms and 90 MB, drawing the identical violin. Violin outliers, like box-plot outliers, draw at most 16 copies of an identical value. (analysis)
- **Web app**: Analysis worker dispatches copy at most 96 MiB per main-thread task; a ten-million-cell density plot used to copy about 1 GB in one 114 ms task. A density over more than 4,096 distinct values sends the page's own values with unit weights instead of a float64 copy and an array of ones (a quarter of the bytes). Quick Insights' category composition yields to the interface like its continuous summary. (analysis)
- **Web app**: Latent-space outlier quantiles are 4 to 6 times faster (10 million cells by 10 columns: 4.6 s to 0.8 s), bit-identical to `cellucid-python` and `cellucid-r`: the latent space is read row by row and large categories are ranked by radix ordering. (analysis)
- **Web app**: Removed the analysis bulk loaders that nothing called (`DataLayer.fetchAnalysisData`, `fetchComprehensiveAnalysisData`, `loadAnalysisBulkData`, `loadLatentEmbeddings`); every analysis gene load goes through the field loader and its cache. (analysis)
- **Web app**: The outlier filter control stays hidden for a field whose every outlier quantile is missing (no category reaches `centroid_min_points`), as for a field without quantiles. (limits)
- **Python**: `prepare()` exports a categorical field where no category reaches `centroid_min_points` with quantized outlier quantiles, publishing every quantile as missing, instead of refusing the export. The mapped latent-space route of `serve()`/`show_anndata()` has no 512 MiB latent-space budget any more and uses `obsm` without caching a copy. Session bundles accept field names and category labels of any length (the 256-character limit is removed). The release build's sdist limit is PyPI's 100 MB upload limit instead of 1 MiB, and the 4,096-member limit is removed. (limits)
- **Python**: `.cellucid-session` reading has no 16 MiB manifest or 512 MiB chunk limit, matching the web app: a frame holds up to 4 GiB − 1 and a length past the end of the file is refused. (limits)
- **Python**: Jupyter: `get_session_bundle()` uploads have no 512 MB limit; the upload streams to disk and is bounded by the free space of the temporary directory (HTTP 413 beyond it). `wait_for_event()` and `set_data_mapping()` are handed their event as it is recorded, so a burst of other events (thousands of hovers) can no longer push it out of the 512-event log they used to scan. (limits)

### Fixed

- **Web app**: Sessions: the vector-field controls (overlay on or off, the field, its colour map and sliders) are restored after the embeddings, the kept views and the active view, because the field select lists the fields bound to the embedding the active view shows. A session saved while that view showed an embedding without vector fields was refused by a page whose view showed one (`Session UI select option "" is not available in current control "velocity-field"`), and the other way round the saved field was overwritten when the embedding came back. The block is marked `data-state-serializer-after-embedding` in `index.html`; any control whose choices follow the embedding can use the same marker. (K3)
- **Web app**: The page preconnects to the sample-catalog origin only when it will fetch from it: a page told another catalog (`?exportsBaseUrl=`, as `cellucid serve`, Jupyter and the browser tests do) no longer opens a connection to `theislab.github.io` it never uses. WebKit reported that unused connection's failed TLS handshake as a console error when the network was poor. (`assets/js/app/ui/core/preconnect-init.js` replaces the static `<link rel="preconnect">`.)
- **Web app**: Benchmark harness: a configuration waits until every presented view draws from its final geometry (its spatial tree, with levels of detail when LOD is on, built and adopted by its store) before a measurement window starts. On a slow host the tree built off the main thread landed inside the window's first frames, so the window measured the landing's uploads and error reads and reported sync stalls that no frame of the configuration makes. The viewer exposes this as `isPresentedGeometrySettled()`.
- **Web app**: Direct H5AD/Zarr: the embedding payload is handed to the viewer instead of being shared with the reader's cache, so evicting a layout now frees its memory and a reload re-reads the file rather than renormalizing an array it already rewrote. Embeddings with a positive range below 1e-8 now get the export format's scale, the same one `cellucid prepare` and `cellucid serve` use, instead of a floor of 1. (DX-14)
- **Web app**: Firefox draws datasets above 30 million points. Every point draw (the point passes, the GPU gather, highlights and figure export) is now split into calls of at most 16,777,216 points; Firefox refused a single draw of more than 30,000,000 and showed nothing. The image is unchanged.
- **Web app**: Per-point textures (visibility, draw-order permutation, the GPU gather's source) are stacked in layers instead of one 2D texture, so no point count is refused for texture size. On GPUs that report a 8,192 texture limit (Firefox on macOS) a 3D dataset above 22.4 million points could not build its draw order, and a dataset above 67.1 million points (268 million elsewhere) was refused at load; the only device limit left is the texture layer count, at least 2^34 points.
- **Web app**: Connectivity edges load and draw for any graph the browser has memory for: a k = 15 neighbour graph was refused from about 12 million edges (1.1 million cells) and aborted the whole dataset. Only the edges a panel shows reach the GPU (the slider's K, 28 bytes each), so GPU memory no longer grows with the graph, and the edge shader reads no texture. The shown edges are the first K visible edges of one fixed random order computed on the fly (a keyed bijection), not a shuffled copy of the graph, so a different, equally random sample of edges is shown than before. Counting all visible edges for the slider no longer blocks a filter change. A graph the browser cannot allocate leaves the dataset loaded without edges.
- **Web app**: Switching back to a representation stays free whenever memory allows: geometry and colour stores kept for switching back are released only when the GPU refuses an allocation (then everything kept is released and the allocation runs again), never ahead of one. A spatial tree whose build failed for memory is built again on the next request once memory has been released, instead of staying off until level of detail or culling was toggled.
- **Web app**: Figure export draws its point passes on WebGL contexts kept for the page instead of creating and losing a context for every pass: exports start faster, and Linux WebKit (WebKitGTK/WPE) no longer crashes the page during an export (it aborted on a GL call made without a current context while a context was torn down).
- **Web app**: Figure export works on synthetic benchmark clouds (it refused them: "normalization for view live is incomplete").
- **Web app**: PNG figure export draws each plot's point raster at exactly the device pixels the plot covers: at 300 dpi a plot's edges fell between device pixels (for example at x = 256.25), so the raster was resampled with a quarter-pixel blur and its bytes depended on the browser's resampling. The raster is now rendered at the snapped size and copied one to one. (figure export)
- **Web app**: The Python documentation names `POST /_cellucid/selection` wherever it describes how the viewer reports back to Python; it had described `/_cellucid/events` alone, which carries every event except selections.
- **Web app**: Analysis loops that hand the event loop back (statistics, scoring, clustering) now yield through the app's one shared macrotask channel. Each yield used to create its own MessageChannel held only by its handler, which WebKit can garbage-collect with the message still pending, leaving the loop waiting forever. (analysis)
- **Web app**: Analysis reads a categorical field at the width it is stored in. A field with more than 65,535 categories (stored as uint32 codes) was copied into 16-bit codes, so every category above 65,535 was counted as another one, in Quick Insights, Detailed counts and marker discovery, without an error; marker discovery also refused more than 32,767 groups. Codes are no longer copied to another width, and a uint8 field is no longer widened for analysis. (analysis, uint32 categories)
- **Web app**: Sessions: restoring a session with several views keeps the live view on the embedding it was turned to. Restoring the views marked the live view as following its dimension's default, so when the session's mapping differed from the one the file opened with (for example a 3-column array saved in 2D), the live view was moved to the default. (data mapping)
- **Web app**: Data mapping: choosing something in the panel (an embedding, a column, Genes, Neighbour graph, Latent space, Vector fields) no longer scrolls the sidebar or drops the focus. The checklist was rebuilt on every change, and Safari's scroll anchoring then scrolled the sidebar by about a screen; rows are now updated in place, and unchanged text, options and lists are not rewritten. Pressing **Default** keeps the keyboard focus on it. (data mapping)
- **Web app**: The antialiasing tooltip and documentation now quote current measurements: in a dense 2D view antialiasing costs about twice the GPU time of a frame (it said about ten times).

## [0.9.3] - 2026-10-20

Version 0.9.3 is a patch release on 0.9.2. It introduces export format 3, so
exports written by 0.9.2 or earlier must be prepared again: a viewer reads
exactly one format. The main changes are:

- **Export format 3**: any number of embeddings per dimension, each with a
  label and its source columns; several gene-name sets; vector fields bound to
  an embedding id; centroids keyed by embedding.
- **Data mapping**: choose, per dataset, which embeddings are on and which one
  each dimension opens on, the neighbour graph, the gene names and expression
  matrix, the latent space for outlier filters, and the vector-field bindings,
  in the browser, from `cellucid serve` and from a notebook.
- **Several representations at once**: every view shows one of the enabled
  embeddings of any dimension, with its own camera and navigation mode.
- **Rendering**: points are stored in the draw order of their spatial tree,
  built in a worker; levels of detail draw a prefix, culling draws tree-node
  ranges, and colours are gathered on the GPU.

### Added

- **Web app, data loading**: drop one `.h5ad` file, one `.zarr.zip` archive or
  one prepared export folder onto the viewer to open it. An outline and a
  *Drop to open* card cover the viewer while you drag; the drop is read like
  the matching button and replaces any load in progress. Anything else is
  refused with a message saying what can be dropped, and a file dropped
  anywhere on the page, even while Cellucid is starting, never makes the
  browser leave it.
- **Web app, cell info**: **Show cell info** (Coloring & filtering, on by
  default) shows a card for the cell under the pointer: its value on the
  field the view is coloured by, where that sits on the colour scale, and for
  a category how many of its cells are visible in that view (what selecting it
  would pick); **Show recent fields** adds the cell's labels in other loaded
  fields. It works in every view of a grid; press and hold on a touch screen.
  A lookup never loads a field or redraws the plot (1M cells: 1.4 ms mean).
- **Web app, citing and sharing**: a **Cite** dialog (APA and BibTeX, read
  from `CITATION.cff`, with a placeholder note until the preprint has a DOI),
  and a **Share** dialog for datasets that can be opened from a link: the
  link, an "Open in Cellucid" README badge and an embed snippet.
- **Docs**: new pages on cell info, share links, badges and embeds, how to
  visualize single-cell data interactively, and how Cellucid relates to other
  viewers.
- **Citation**: a placeholder citation for the forthcoming Cellucid bioRxiv
  preprint, the same in every repository. It is written once, as the
  `preferred-citation` of `theislab/cellucid`'s `CITATION.cff`, and copied from
  there to this repository's `CITATION.cff`, the README, and the new
  [Citing Cellucid](https://docs.cellucid.com/en/latest/citing.html) page;
  the test suite fails if a copy drifts. Until the preprint has a DOI, every copy
  says it is a placeholder.
- **Web app**: cells can be highlighted on a touch screen. Every selection tool
  started only from an `Alt` press, which a phone or tablet does not have, so touch
  users could not select at all. A **Select mode** switch in the Highlighting panel
  now stands in for a held `Alt`: while it is on, a one-finger tap or drag runs the
  active tool through the same path an `Alt` press takes, and its Intersect / Add /
  Subtract buttons stand in for `Shift` and `Ctrl`. Two fingers always pinch to zoom
  and, newly, drag to pan; a second finger landing mid-selection abandons that
  gesture rather than committing it. `Esc` with nothing in progress leaves the mode.
  Mouse and keyboard behaviour is unchanged while it is off.

- **Web app**: a **selection bar** at the bottom of the plot carries Select mode's
  controls, so a selection can be made, combined, confirmed, cancelled and cleared
  without opening the menu, which covers the whole screen on a phone. It shows the
  same switch and Intersect / Add / Subtract choice as the panel, the active
  tool's cell count with its Undo, Cancel and Confirm while a selection is in
  progress, the active page's highlighted count with Clear, and a one-line note
  when a tap on empty plot changed nothing. Every control on it is the panel's
  control shown again, so the two always agree. It is always shown on touch
  screens (just the switch while the mode is off), shown to a mouse while the
  mode is on, centred in the part of the plot the sidebar leaves uncovered, and
  notifications stack above it.
- **Web app, backgrounds**: two new viewer backgrounds, **Studio (light)** and
  **Studio (dark)** (neutral grey and charcoal): a 3D view stands on a floor
  of constant-width hairlines that fades out before its edge, the same ruled
  floor from every angle above it (below it, the stage is plain); 1D, 2D and
  Planar views sit in front of dotted paper. Which one a view shows follows
  its dimension and navigation mode, never the camera angle, so rotating a 3D
  view never swaps it to dots. The pattern rescales with zoom without popping, and
  figure export reproduces it in SVG and PNG. **Studio (light) is the new
  default** (the sample datasets open on it too), and the Background menu is
  grouped into Light and Dark.
- **Web app, shader quality**: a fourth quality, **Sculpted (soft 3D
  light)**. Every cell is a soft-lit, glossy ball — a wrapped key light, a
  sky-and-ground ambient, a broad gloss instead of a hard white highlight, and
  darkening toward the edge that keeps touching cells apart — computed in
  linear light so the centre of each ball is exactly its legend colour. With
  antialiasing on, its edges are smooth in 3D views too, without the
  background showing through nearer cells; highlighted cells get the same 3D
  ring as `Full`, held still. It costs about what `Full` does (3–5% more in
  3D views on an Apple M1 Pro, the same in 2D). `Full` stays the default.
- `var_gene_name_columns=` on `prepare()`: further `var` columns offered as
  alternative gene names (for example symbols beside an Ensembl index).
- `embedding_keys=` (`--embedding-1d/2d/3d KEY[:c0,c1]`), `expression_source=`
  (`--expression X|raw|layer:NAME`) and `--gene-id-column` for serving AnnData.
  Without `--expression`, a served object reads `X`, or, without `X`, `raw`,
  its only layer, or its best-named layer (`norm`/`log` first, `counts` last),
  and reports that guess -- the rule a file opened in the browser follows.
- `cellucid.data_mapping`: the data-mapping model shared with the web app and
  the R package, with a conformance fixture
  (`tests/fixtures/data_mapping_conformance.json`) and a cross-writer byte
  parity fixture (`tests/fixtures/writer_parity_v3.json`).
- A sixth reconciled surface, **Embedding**: the writer proves `embeddings/`
  holds exactly the files `dataset_identity.json` declares.
- `cellucid serve`, `serve_anndata()` and `show_anndata()` with AnnData publish
  the data mapping they serve at `/_cellucid/mapping/candidates.json` (as
  `defaults`, beside every candidate `obsm`, `obsp`, layer and `var` column,
  read from metadata only), and `AnnDataAdapter` gains `get_data_mapping()` and
  `get_mapping_candidates()`. The web app reads it to record the mapping a
  session was saved with. Every other narrow `obsm` array is turned on beside
  the served defaults and served on request at
  `/_cellucid/mapping/embeddings/{obsm}/{columns}.bin` (see the data-mapping
  routes below), so a served file shows
  any of its embeddings and its checklist works as for a file read in the
  browser.
- **Web app, data mapping**: a **Data mapping** group at the bottom of the
  Session panel lists every embedding the dataset can be shown with as a
  checklist: turn any number on, of any mix of 1D, 2D and 3D, and mark one per
  dimension as its **Default**. For an `.h5ad` or Zarr ZIP opened in the
  browser every numeric `obsm` array is listed, a wider one at any selection
  of 1 to 3 of its columns (`X_pca` columns 2 and 3 as `PCA 3,4`), as many
  selections as you like; a prepared dataset lists the embeddings it was
  exported with, all on. Turning an embedding on loads nothing until a view
  shows it. Choices are staged and applied as one step that loads and checks
  every embedding a view must switch to first; the neighbour graph, genes,
  latent space and vector fields have disclosures of their own (below). The
  data mapping's embeddings are now `{enabled, defaults}`, in
  the web app, `cellucid.data_mapping` and the R package's conformance fixture.
- **Web app, several embeddings at once**: every view shows one of the
  embeddings that are on, of any dimension — two 3D embeddings side by side, a
  3D one beside a 2D one, a t-SNE beside a UMAP — chosen in **View N shows**
  under **Compare views**, grouped by dimension with the default first. Keys
  `1`, `2` and `3` show that dimension's default in the selected view and,
  pressed again, its next embedding ("2D · t-SNE (2 of 3)"). Each view row
  names what its view shows (`UMAP · 3D`), its field, its navigation mode, and
  whether its camera is linked. Views of one embedding share its coordinates
  and GPU buffers; coordinates no view shows are kept in a cache bounded by a
  quarter of the device memory (512 MiB to 2 GiB) and released least recently
  shown first. Up to four views (the live view and three kept views).
- **Web app, per-view cameras**: each view owns its camera and navigation
  mode. **Navigation** is its own accordion and edits the selected view; its
  **Applies to** row names that view and the views linked to it. **Link
  cameras** (on by default) links only views that show the same embedding, and
  is off and disabled while no two views do.
- **Web app, rendering**: points are stored on the GPU in the draw order of
  their spatial tree (by LOD band, then tree leaf or coarse node, then cell),
  in every dimension, whenever the tree exists. A reduced level of detail
  draws a prefix of the stored points, frustum culling draws whole
  tree-node ranges (one multi-draw call of at most 1,024 ranges, 256 one call
  a range without `WEBGL_multi_draw`, bridging the gaps that hide the fewest
  points, with no per-point index upload), and full detail draws faster
  because neighbouring points are neighbours in memory: at 10M points about
  three times faster in 3D and in 2D. LOD no longer keeps a second copy of the
  positions: about 200 MB of GPU memory at 10M points instead of 320. The
  tree, its levels of detail and its draw order are built in a worker, so
  turning LOD or culling on, loading a dataset or showing another
  representation never blocks the page (it took about 2.2 s at 10M points);
  the view draws full detail until the tree lands, then adopts its order in
  one frame, gathered on the GPU. Colours are gathered on the GPU too, so a
  recolour costs one upload of 4 bytes a point and one GPU pass whether LOD is
  on or off: 5-8 ms at 10M points (Chromium), where a recolour with LOD on
  used to take about 68 ms. A view switched back to a representation it
  showed keeps its colours without an upload, and a representation switch
  back to a layout already shown no longer rescans its positions. Translucent
  points in 2D overlap in tree order, the same in every view, every level of
  detail and figure export.
- **Web app**: the in-browser position budget (512 MiB where the browser
  reports no device memory) never refuses the embeddings the views show plus
  the one being shown; it only evicts cached embeddings no view shows.
- **Web app**: a notification can carry one action button; the best-guess
  notification uses it (**Show**) to open the Data mapping section.
- **Web app, usage analytics**: `data_load` reports how the embeddings were
  chosen (`mapping_origin`: `default` or `best_guess`) and how many the
  dataset could show (`embedding_candidates`: `1`, `2-3`, `4-9`, `10+`);
  how many it opened with turned on (`enabled_representations`); `feature_used`
  counts `data_mapping` (`embedding`, `embedding_columns`, `default`) and
  `dimension` when the picker switches a view to another dimension. No key or
  embedding name is ever sent.
- **Further data-mapping routes** beside that catalog: on request, the server
  serves any catalog entry -- an embedding at chosen columns,
  a vector field bound to an embedding, a neighbour graph, the gene list and
  values of `X`, `raw` or a layer named by any `var` column, and outlier
  quantiles in any latent space. Entries are addressed by catalog position,
  never by key. An unknown index is 404; an entry whose values fail their check
  is 422 with `{"error", "key", "reason"}`, so a graph that is not symmetric is
  refused, not repaired. `reason` is the short clause (`NaN at cell 17`) a file
  opened in the browser gets in the same words, and the catalog lists an array
  the best guess already refused as unusable with that reason, as the browser
  does. HEAD answers lengths without reading payloads, and one
  byte-bounded, lock-guarded cache serves the server thread and the notebook
  together. `AnnDataAdapter` gains the `get_mapped_*` methods behind them.
- **Jupyter:** the `setDataMapping` command and the `data_mapping` event, both
  listed by `/_cellucid/protocol`. `viewer.data_mapping` (and
  `viewer.data_mapping_origin`) hold the mapping the viewer last reported,
  `@viewer.on_data_mapping` runs on each report, and
  `viewer.set_data_mapping(**fields)` changes the embeddings, graph, gene names,
  expression, latent space or vector bindings from the notebook, validated before
  anything is sent.
- **Web app, served data mapping**: a file served by `cellucid serve` or shown
  with `show_anndata()` can show any numeric `obsm` array, at any columns, from
  **Session › Data mapping**, as a file opened in the browser can; the server
  reads the array on request, and its centroid labels are computed from its
  positions in the viewer, as for a file read in the browser. A choice
  the server refuses is listed unusable with the server's reason and the
  previous mapping stays.
- **Web app, neighbour graph**: for an AnnData file (in the browser, served or
  in a notebook), **Session › Data mapping › Neighbour graph** chooses the
  graph connectivity edges and the KNN drag tool use: None or any square
  `obsp` matrix with one row per cell. The graph is read and checked when it is
  applied; one that is not symmetric (`obsp['distances']`) is refused in the
  server's words and stays unavailable, never repaired. Shown edges are redrawn
  from the new graph and **Dataset › Connectivity** counts its edges. A session
  restores the graph it was saved with when the dataset has it. In the browser,
  a file keeps its own graph and at most one other in memory, and their edges
  are byte-identical to what `cellucid serve` sends for the same file.
  `feature_used` counts `data_mapping` / `knn`.
- **Web app, gene names and expression**: for an AnnData file (in the browser,
  served or in a notebook), **Session › Data mapping › Genes** chooses the
  expression matrix (`X`, any layer, `raw.X` with its own gene list, or None)
  and the gene names (the var index or any text, categorical or integer var
  column; a floating-point or true/false column is listed with the reason it
  cannot name genes). Empty names take the var-index name and repeated names
  are made unique (`MATR3`, `MATR3-1`) by the rule `prepare()` and the server
  use, and the panel says what that changed and when Ensembl ids have a symbol
  column beside them. Apply reads the new list first; a view coloured by a
  gene keeps it with its colour range, log scale and filter (by position while
  the var axis stays, as between `X` and its layers; by name to or from
  `raw.X`), a gene not on the new
  list loses its colouring and is named, cells and highlights are untouched,
  and analysis windows close. A prepared dataset switches between the
  gene-name sets it was exported with. A session restores the gene names and
  matrix it was saved with and checks its gene count against that matrix.
  The browser reads one matrix at a time (`X`, `layers/<name>` or `raw/X`),
  and its gene lists and values are byte-identical to what `cellucid serve`
  sends for the same file. `feature_used` counts `data_mapping` /
  `gene_names` and `expression`.
- **Web app, latent space**: for an AnnData file, **Session › Data mapping ›
  Latent space** chooses the `obsm` array (or None) a categorical field's
  outlier filter measures distance in. Apply reads and checks the space and
  computes the quantiles of every loaded categorical field before anything
  changes, off the main thread (`data/latent-outliers.js`), bit for bit what
  `prepare()` and `cellucid serve` compute; a space with a bad value is refused
  in the server's words and stays unavailable. Every view keeps its outlier
  threshold; None removes the outlier filters. A space larger than 512 MiB as
  float32 is listed with that reason, as the server lists it.
- **Web app, vector fields**: **Session › Data mapping › Vector fields** binds
  vector arrays to embeddings with one checkbox per array that fits each
  embedding: as wide as its dimension (read whole) or as wide as its array
  (read at the same columns, `velocity_pca` on `X_pca[2,3]`). Arrays named for
  the embedding are marked suggested. A prepared dataset lists the bindings it
  was exported with, any of which can be turned off. Apply rebinds and turns
  the overlay off; vectors are read, checked and scaled only when the overlay
  shows them, through `_cellucid/mapping/vectors/...` on every AnnData path,
  and the overlay holds GPU data only for what it shows. `feature_used` counts
  `data_mapping` / `latent` and `vector_fields`.
- **Web app, Jupyter**: the viewer reads `GET /_cellucid/protocol` when it
  connects. It applies `setDataMapping` as one transaction and answers it with
  a `data_mapping` event naming the request; it reports every other mapping
  change, and reports a refused notebook command as `command_error`, to a
  `cellucid` that lists those events.
- `/_cellucid/info`'s `version` is the documented signal a newer viewer uses to
  tell a user their server predates the viewer's export format ("This Cellucid
  server (0.9.2) is older than the viewer; update it with pip install -U
  cellucid.").

### Changed

- **Latent-space outlier quantiles** follow one stated rule in `prepare()`,
  `cellucid_prepare()` in R (now in C), `cellucid serve` and the web app:
  float32 values widened to float64, centroid sums accumulated in cell order,
  squared distances in column order, and the share at that distance or closer
  as float32 (numpy's pairwise summation and R's long-double sums are no longer
  part of the result). `tests/fixtures/latent_outlier_conformance.json`
  (written by `scripts/generate_latent_outlier_conformance.py`) pins the bytes
  in all three repositories. Quantiles can differ from 0.9.1's in the last rank
  only where two distances tie within rounding.
- **`prepare()` gene names from an integer column**: a `var` column of integers
  names genes by their digits (`7157` is `"7157"`), the rule the viewer,
  `cellucid serve` and `cellucid_prepare()` in R apply to the same column,
  instead of being refused. A floating-point or boolean column is refused with
  the catalog's reason (`floating-point values cannot name genes`).
- **`cellucid serve` / `show_anndata()` vector fields**: without
  `--vector-fields` (`serve_vector_fields=True`) the served data mapping binds
  no vector field, matching what the server declares; with it, the default
  rule binds `<field>_umap_<d>d` as before. The startup line now reads
  `Vector fields: not bound (... to bind them, or bind one under Session > Data
  mapping)`.

- **Export format version 3.** One export can now carry several embeddings per
  dimension and several sets of gene names. This is a clean break: exports in
  the earlier format are not read, and format-3 exports need a web app that
  reads format 3. Re-export with `prepare()` (or `cellucid_prepare()` from the
  matching R release, which writes the same bytes). `dataset_identity.json`
  has `version: 3`, and `embeddings` is now `{default_dimension, defaults, items}`: every exported
  embedding is an item with a `label`, a `dimension`, the exact payload `file`
  `embeddings/{index}_{dimension}d.bin(.gz)` and its `source` (`key` and the
  0-based `columns` it was read from). Items are ordered by dimension, then id;
  `defaults` names the item each dimension shows first. The `points_*d.bin`
  files at the export root are gone.
- **Vector fields bind to embedding ids.** A field's `files` are keyed by
  embedding id, `vectors/{fieldIndex}_{embeddingIndex}.bin(.gz)`, scaled by that
  embedding's normalization factor. `basis`, `available_dimensions` and the
  per-field `default_dimension` are removed; `label` is the field id.
- **Manifests are `compact_v2`.** A categorical field's element `[5]` is
  `centroidsByEmbedding`, keyed by embedding id. `var_manifest.json` replaces
  `var_gene_id_column` with `gene_name_sets`: the default set (from
  `var_gene_id_column` or the var index) and any alternative sets, stored as
  metadata only, so a reader can switch gene names without reading expression
  again.
- **`prepare()` API:** `X_umap_1d`, `X_umap_2d` and `X_umap_3d` are replaced by
  `embeddings={id: array or Embedding}` (required) and `embedding_defaults`
  (required for every dimension with more than one embedding). The new
  `cellucid.Embedding(values, columns=None, label=None, source_key=None)` reads
  chosen 0-based columns of a wide array, such as `X_pca` columns `[0, 1]`.
  `vector_fields` is now `{field: {embedding id: array}}`.
- **Gene names are made unique instead of refused.** Repeated names become
  `MATR3`, `MATR3-1`, ... by the rule anndata's `var_names_make_unique()` uses,
  and an empty name takes its row's var index name. The same rule runs in the
  Python writer, the R writer and the AnnData server, and what changed is
  printed and recorded in `var_manifest.json`. `gene_identifiers` matches the
  unique names.
- **`cellucid serve` / `serve_anndata()` / `show_anndata()` for AnnData** serve
  identity v3. `serve_connectivity` and `--connectivity` are replaced by
  `knn_key` and `--knn KEY`: the named `obsp` matrix is read and checked when it
  is named, and one that is not symmetric (such as `distances`) is refused with
  that reason instead of being symmetrized.
- When an object has no `X_umap_1d`/`X_umap_2d`/`X_umap_3d` and no 1- to
  3-column `X_umap`, the server no longer stops with rename hints. It makes a
  deterministic best guess from the other `obsm` arrays (t-SNE, force-directed
  layouts, PHATE, spatial coordinates and diffusion maps are preferred over
  latent spaces such as `X_pca`) and reports it; when nothing can be plotted it
  says what `obsm` holds and why each entry was not used.
- Sessions: `datasetFingerprint` must carry `cellOrder` (now naming the
  `embedding` its digest was taken on) and the new `dataMapping`. User-defined
  categorical fields no longer store `centroidsByDim` or `normalizedDims`.
  `apply_cellucid_session_to_anndata()` checks the mapping against the AnnData,
  compares `varCount` on the mapped gene axis, and resolves var-sourced fields
  through the mapped gene names and expression matrix.

- **Web app, data**: the viewer reads export format 3 and only format 3 — from a
  prepared folder, the sample datasets, GitHub, a server, and the notebook —
  and opens `.h5ad` and Zarr files as the same format. Every embedding an
  export carries is loaded only when a view shows it; each dimension opens on
  the embedding `defaults` names. An `.h5ad` or Zarr file's vector field now
  reports its id as its label (`velocity_umap`, not `Velocity (UMAP)`), as the
  format requires, and the reader now refuses any export whose vector-field
  label is not its id, from every source.
- **Web app, data**: an `.h5ad` or Zarr ZIP without a UMAP is no longer
  refused. It opens on the same deterministic best guess `cellucid serve`
  makes (t-SNE, force-directed layouts, PHATE, spatial coordinates, …, and last
  the first two columns of a wide array such as `X_pca`), marked *guessed*,
  with a notification naming what is shown. Only a file with nothing drawable
  is refused, with the message `cellucid serve` prints: "Cellucid found
  nothing to plot…", listing what `obsm` holds and why each entry was not used
  (the error code is `NO_PLOTTABLE_EMBEDDING`). Every embedding is checked
  before it is drawn — NaN, infinity, a value outside 32-bit float range, or
  no coordinate variation is refused naming the first cell — so a one-cell
  file can no longer be opened. A one-dimensional `obsm` array is read as one
  column, and a wide array is read only around the columns a view needs.
- **Web app, vector fields**: the overlay is offered by the embedding a view
  shows, not by its dimension: a field is drawn only in views showing the
  embedding it is bound to, and the info line names those embeddings
  ("Vector fields are bound to … Show one of them in this view to enable.").
- **Web app, sessions**: a session records its data mapping (the embeddings
  that are on and each dimension's default), the embedding its cell order was
  taken on, the embedding every view shows, every view's camera and navigation
  mode, and whether cameras are linked. Loading checks the dataset, its data
  mapping and the cell order of the embedding the session names before
  anything changes on screen, loading that embedding if another is shown — so
  a session saved in 2D now opens while 3D is shown, instead of asking you to
  switch back. A session whose embeddings this dataset can show turns them on
  before its views are restored, and one naming an embedding the dataset does
  not have is refused naming it ("This session shows obsm['X_tsne'] columns
  [0, 1], which this dataset does not have."). A
  session saved by an earlier release is refused with a message saying so.
  User-defined categorical fields no longer store their centroids; they are
  recomputed from the codes.

### Removed

- `vector_fields.missing_embedding_message()`, `resolve_embedding_keys()`,
  `embedding_declaration_hints()` and the other helpers of the one-UMAP-per-
  dimension rule, together with `AnnDataAdapter.get_points_binary()`.

### Fixed

- **Web app, analysis performance**: analysis no longer freezes the page on
  large selections, and every result is the same number it was. Medians,
  quartiles, ranks and sorted copies (Quick, Detailed summaries and tests, CSV
  export, Spearman correlation, Kruskal–Wallis, high-cardinality Wilcoxon in
  differential expression and Marker Genes, gene signatures) come from
  linear-time radix ordering instead of comparison sorts; correlation pairs its
  colour field by merging sorted cell indices instead of a per-cell map; box,
  violin and density plots of large pages compute their statistics on the
  analysis workers; a histogram hands the plot a few values per bin instead of
  every value; and a box plot draws at most 16 copies of an identical outlier
  (all on one spot) instead of one SVG element per cell. At one million cells:
  box plot of a zero-inflated gene 1.1 s → 0.12 s (longest main-thread task
  70 ms), histogram 0.17 s → 0.06 s, Detailed summary 0.10 s → 0.03 s,
  three-page Kruskal–Wallis 0.67 s → 0.09 s, Spearman 0.86 s → 0.14 s,
  high-cardinality marker statistics 0.33 s → 0.06 s per gene. The Plotly
  notification watcher no longer runs a subtree query for every element any
  plot inserts.
- **Web app, data loading**: when the connection to a server is lost, Cellucid
  clears only the dataset that came from it, and not while another dataset you
  chose is still loading; it used to clear over that choice. Recovering from a
  lost connection while its dataset was on screen also failed outright
  ("Remote dataset URL requires an active connection"), because the open
  dataset's address was re-read from the closed connection; it is now the
  address the dataset was loaded from. A server with one dataset that answers
  **Connect** after you applied a session no longer replaces that dataset and
  session; its dataset is listed in the picker instead.
- **Web app, data loading**: a dataset chosen while another is still loading
  now replaces it. At startup the data controls stayed disabled until the
  startup dataset had downloaded, and **Load State** did nothing at all in
  that window; afterwards, a file still being read published over a sample or
  **None** chosen after it. Now the controls open once the dataset catalog has
  been read, and every dataset change (the startup dataset, a sample, **None**,
  a file or folder, a server or GitHub repository with one dataset, the
  benchmark) supersedes the load in flight the moment you act: its downloads
  are cancelled, it never appears late, and the dataset picker shows what is
  loading. A connection that answers after a newer choice keeps its dataset in
  the picker instead of replacing that choice, and a session chosen during a
  load waits for that dataset, says so, and is then applied.
- **Web app, rendering**: switching from **Grid (light)** to **Grid (dark)**
  no longer leaves the scene redrawing on every frame (the grid fade only
  moved upwards, so it never reached the lower target). With reduced motion
  requested by the system, background fades and crossfades are instant.
- **Web app, rendering**: connectivity edges declare their per-vertex
  attribute at location 0. Left to the linker, Firefox on macOS could draw
  them without an enabled attribute 0 and emulate one per vertex.
- **Web app, no dataset**: after choosing **None**, Keep view, figure Export,
  Save and Connect repo are unavailable and say why instead of throwing;
  Differential expression, Gene signature and Marker genes show `Requires gene
  expression data.` (as on a dataset without genes) instead of a page error;
  Analyze performance and Copy situation report ask for data; and Load says no
  dataset is open instead of quoting `0 cells and 0 genes`.
- **Web app, analysis**: marker-gene results are cached per gene list, so
  results computed with one set of gene names or one expression matrix are
  never shown for another.
- **Web app, data**: a dataset re-exported at the same address can no longer
  show values from its previous version. Gene and observation payloads were
  reused from the browser cache without asking the server, and payload files
  are numbered, not named: a re-export that added one gene left
  `var/3.values.u8.gz` holding another gene of exactly the same size, so the
  cached file drew the old gene's values under the new gene's name. Every
  dataset file — catalogs, manifests, payloads and saved views — is now
  revalidated with the server before use; an unchanged file costs a short
  "not modified" reply.
- **Web app, interface**: one calmer design. One font (Inter) in sentence case
  with a small type scale; soft borders; compact spacing. Every panel has
  short headings, and every explanation moved from the panel into an (i) tip
  that links to the exact docs page. Navigation is its own panel and shows
  which view it drives when several are open; the benchmark panel is called
  Benchmark; Visualization, Navigation and Compare views start collapsed; the
  Outlier filter sits under Legend; the logo follows the theme. Slider thumbs
  no longer overhang the rows around them, and the legend tick no longer
  touches its swatch.
- **Web app, rendering**: a hover lookup and any click on the page no longer
  redraw the scene for a second.
- **Web app, picking**: hovering and clicking now pick the dot drawn under the
  pointer, at the size it is drawn, instead of the first cell near the view
  ray, which could be a nearer cell beside the pointer. Where dots overlap the
  front-most wins; filtered-out cells and cells the current level of detail
  does not draw are never picked. Cell info, Alt+click selection, proximity
  and KNN all use it, each in the pane under the pointer.
- **Web app, selection**: in a grid of views, Alt+click and Alt+drag in a
  pane other than the focused one now select by that pane's field. They used
  the focused view's field, so a click in a pane colored by `cell_type`
  beside a focused view colored by a continuous field selected a value range
  of the wrong field. The category, the value range, the drag preview and its
  range label now all come from the pressed pane, among the cells it shows;
  view focus still does not move. A pane colored by no field now says so.
- **Web app, selection**: whether an Alt press can start, and the armed
  cursor, now follow the pane under the pointer instead of the focused view,
  so Alt+click works in a colored pane while the focused view has no field,
  and the panel asks for a field only when no pane on screen has one.
- **Web app, cell info**: the card no longer disappears for good. A key whose
  release the page never saw (a Cmd shortcut such as the macOS screenshot
  keys, or a key held while the window lost focus) read as camera motion
  forever, and so did a drag released outside the page; the card waits for
  the camera to settle, so it never came back. Held keys are now only the
  navigation keys, are dropped on blur, and Shift alone no longer counts as
  motion. The card no longer repeats the pane's name in a grid of views.
- **Web app, benchmark**: the presets are now 500K, 1M, 5M, 10M, 20M and 40M.
- **Web app, welcome**: the welcome dialog opens with no choice pre-selected;
  focus rests on the dialog and the first Tab reaches **Choose a dataset**.
  The sidebar footer links GitHub next to Docs and Cite.
- **Web app, onboarding**: closing the welcome dialog, by any path, no longer
  focuses the **Sample dataset** dropdown; **Choose a dataset** opens the
  **Session** panel, and from the keyboard moves focus to its heading.
- **Web app, data loading**: dataset catalogs, identities and manifests are
  revalidated on every load, so for a few minutes after a dataset is updated
  a stale cached copy can no longer clash with the new files (for example a
  sample view failing with a missing control).
- **Web app, pointer**: hover, click and selection read the camera as it is,
  not as the last frame drew it, so a pointer read right after the camera
  moved no longer lands on the previous picture.

- **Web app**: automatic antialiasing now accounts for the view's dimension.
  A 1D or 2D view has no depth test, so every point is shaded and antialiasing
  costs about ten times more there (1M cells: 33.8 ms with it, 3.4 ms without).
  It is now on below 5,000,000 cells in 3D and below 200,000 cells in 1D/2D,
  and the automatic choice is re-decided when the active view switches
  dimension. The status line names both (`Chosen automatically: 3,696 cells in
  3D.`); an explicit choice is never changed.
- **Web app**: less memory per cell, so larger datasets fit. Per-cell visibility
  is one byte (exactly 0 or 1) instead of four, the outlier quantiles are no
  longer copied into the state, categorical codes keep their stored width, and
  a kept view shares the visibility arrays copy-on-write instead of cloning
  them. Loaded fields are released by bytes rather than by count, so browsing
  genes no longer grows memory without bound. A loaded 3D dataset went from
  36 to 22 bytes per cell, and each kept view from 26 to about 6; 150M cells in
  2D and 120M in 3D load and render on a 16 GB machine. Pixels, filters and
  exports are byte-identical. The benchmark's pre-publication GPU estimate now
  follows the dimension (13 B/point in 2D, 17 in 3D) instead of a fixed 28.

- **Web app**: usage analytics rebuilt around one declared schema
  (`assets/js/analytics/schema.js`). It no longer sends dataset names, ids built
  from local file names, exact counts, error text or full page URLs (the page
  view carried `remote=` and `github=` addresses and their tokens); a dataset is
  described by count ranges, dimensions and writer version, with an id only for
  the public samples. Each load is reported once (it was reported twice), timed
  from the click or file pick (local loads excluded the file read), and reported
  when it fails before activation (such failures sent nothing). New events cover
  startup success and failure, feature adoption, and device classes. Every
  payload is validated against GA4's limits, and analytics is never awaited.

- **Web app**: faster, quieter full-detail rendering with identical pixels.
  Hovering the plot or typing into a sidebar field redrew the whole scene for a
  one-second settle window (237 full frames over two seconds of hover at 1M
  cells); only input that changes something now draws. In 2D views — drawn
  without depth testing — the point shaders use a single `discard`, and `Light`
  quality uses a no-depth program that blends a zero-alpha fragment instead of
  discarding: at 1M cells with antialiasing, `Light` 30.7 → 17.0 ms and `Full`
  33.8 → 30.6 ms per frame, byte-identical in 45 compared scenes. Orbit and pan
  inertia, smooth zoom and held-key acceleration now scale with elapsed time, so
  a 120 Hz display no longer moves the camera at double speed, and the render
  loop allocates about half as much per frame.

- **Web app**: with an ad or tracking blocker (uBlock Origin, Brave, AdGuard, iOS
  content blockers — anything shipping EasyPrivacy, which lists
  `/analytics/tracker.js`), Cellucid did not start and showed “Cellucid could not
  load its code…”. The analytics module was a static import, so its refused
  request aborted the whole module graph. Analytics is now loaded with a dynamic
  import and stays off when refused; the dataset-load vocabulary it shared with
  the controls moved to `app/data-load-methods.js`. No other application file
  matches the EasyList, EasyPrivacy, uBlock Origin or AdGuard lists.

- **Web app**: choosing volumetric smoke on a device without float blending
  (Safari on iPhone and iPad, many Android browsers) failed mid-build with
  “Smoke density unavailable: GPU smoke density requires EXT_float_blend for
  exact R32F accumulation.”, and restoring a smoke session there threw. Smoke's
  requirements are now checked up front; such a device gets a plain message
  (“Volumetric smoke isn't available on this device: …”) and stays in Points,
  sessions restore into Points, and other smoke errors no longer show internal
  detail.

- **Web app**: the volumetric smoke cloud is now **beta** rather than alpha. It
  is coloured, bounded and idle-aware after this release's rewrite; the tag beside
  `Render mode:` reads `Beta` in the calmer info tone, and its hover text no
  longer says the mode is still changing.

- **Web app**: analysis plots of large pages draw several times faster, with
  identical results. Plotly's violin density evaluated its kernel once per value
  at every curve point; a local, test-pinned patch to the vendored Plotly
  evaluates it once per run of equal values (expression data is mostly ties) and
  skips terms that are exactly zero, adding the rest in the same order. Box plots
  that draw points hand Plotly exact precomputed statistics; density curves are
  computed once (not twice), reused across cosmetic changes, and split across the
  worker pool when large. At 2,000,000 cells: a gene violin 3.7 → 0.67 s, a
  continuous box 1.15 → 0.26 s, a continuous density 1.39 → 0.56 s. Page values
  keep their own numeric type, the page-data cache is bounded at 256 MiB as well
  as 100 entries, and decodes go to the least-loaded worker.

- **Web app**: on browsers that composite in software (blocklisted GPUs, virtual
  machines, remote desktops), part of the plot vanished beside a floating panel
  whose content scrolled: Chromium culled the canvas under the panel's opaque
  scroller and drew the hole at its mirror image. Each floating window is now its
  own compositing surface.

- **Web app**: multiview and the data path are faster. Kept views drew their
  reduced Level-of-Detail levels through a gather from the full-detail streams;
  they now draw the live view's compact positions with their own compact colors,
  and share the live positions instead of copying them. Measured on an Apple M1
  Pro (Chromium, ANGLE Metal), 10,000,000 cells in 2D with LOD: two views
  15.9 → 6.9 ms a frame, four views 43.6 → 13.5 ms, and each kept view now adds
  38 MB of GPU memory instead of 114 MB. Hiding a category took 225 ms instead
  of 413 ms, a cached recolor 225 ms instead of 321 ms, and loading a prepared
  export 1.30 s instead of 1.43 s (30,000,000 cells: 2.89 s instead of 3.45 s).
  The documentation's claim that a four-view grid costs about as much as one
  view held only for large, overlapping points and is corrected.

- **Web app**: a panel's undock, dock, copy and close buttons no longer sit inside
  its header's open/close toggle. A `<summary>` is itself a button, so each was a
  control nested in a control — announced poorly by screen readers and flagged by
  accessibility checkers. They now stand beside the header, look and sit exactly
  where they did, are named after their panel ("Undock Visualization panel into a
  floating window"), and come just before the header in `Tab` order.

- **Web app**: multiview panes are laid out in the part of the window the sidebar
  leaves uncovered. The grid was cut from the whole canvas, so with the sidebar
  open the left pane was partly hidden under it (at 1440 px with two panes, 430 px
  of the left pane was visible against 720 px of the right one). Every pane now
  has the same size beside the open sidebar, and the layout follows at once when
  the sidebar is closed, opened, resized, or the window or its pixel density
  changes. One layout drives drawing, picking, hover, focus, every selection
  tool, projectiles, the pane title chips, centroid labels, figure export, and
  the selection bar, so what you click is always the pane you see. Wheel zoom
  now focuses the pane under the cursor, free-fly is centred in the uncovered
  area like Orbit and Planar, and the hidden sidebar is `inert`, so `Tab` no
  longer walks into it.

- **Web app**: the sidebar's ✕ works from the first frame. It was wired only once
  the first dataset had loaded, so during a slow first load a phone user, who
  sees the menu full-screen, tapped it and nothing happened.

- **Web app**: smaller layout fixes.
  - On a short screen (a phone on its side), startup notifications covered the
    welcome dialog's buttons; while it is open they now wait beneath it.
  - A long view title in a narrow multiview pane ran into the next pane's
    title; it now ends in an ellipsis inside its own pane, with the full name
    on hover.
  - The benchmark's Frame Time reading wrapped onto three lines. The second
    number now sits on its own line, labelled `CPU` for what it measures: the
    main-thread time spent issuing the frame, not GPU time.
  - Figure export's Background and Font selects now have visible labels, with
    the same trailing colon as every other label.

- **Web app**: in multiview, a finger dragging a pane now focuses that pane first,
  as a mouse press does, instead of moving whichever pane last had focus, and a
  touch drag saves the camera to the focused pane when cameras are unlocked.

- **Web app**: the volumetric smoke cloud is rebuilt and fixed.
  - It drew near-black soot whatever the cells were coloured by. The cloud now
    takes the active field's colours: each region shows the mean display colour
    of the cells there, lit by a key light and sky occlusion pre-computed once
    per build rather than by a shadow ray from every sample.
  - Switching dimension or recolouring (field, colormap, category colour, log
    scale) left the previous volume on screen. Both now rebuild it, and a 2D or
    1D embedding is spread into a soft slab instead of a sheet rays stepped over.
  - `Ray quality` changed almost nothing. It now sets samples per voxel, 0.5 to
    2 (`High`, 1.41, by default), and empty space is skipped at no cost.
  - The cloud shimmered and the GPU never idled: the ray-start jitter moved
    every frame and every frame was redrawn. The jitter is now fixed, and at
    `Animation speed` 0, which reads `Still`, a still view is not redrawn.
  - `Render resolution` reached 2.00x and `Noise detail` 256³, either of which
    could stall one frame long enough to lose the WebGL context. Render
    resolution now runs 0.25x–1.00x (default 0.50x) and a frame marches at most
    2560×1440 rays; noise detail is 32³–128³ (default 128³).
  - `Cloud density`, `Fine detail` and `Light absorption` printed readouts that
    changed with grid and noise resolution, so a recorded value did not
    reproduce a look. Every smoke slider now has one fixed scale.

  The documentation's control tables, recovery order and pipeline description
  are corrected (grid size was never the cause of context loss: the grid is
  capped at 128³, an 8 MiB volume), and the sample datasets' saved views now
  open with the new smoke defaults.

- **Web app, user-defined fields**: a categorical made with
  `Create Categorical`, a category merge or delete, or a field duplicate now
  shows its category labels, keeps them when a view switches dimension (that
  switch threw `Raw 2D positions are unavailable for centroid computation`),
  and gets them back from a saved session exactly where they were; restoring a
  session whose active field is one of these no longer fails with
  `Active categorical field is absent from its rendered selector`.
- **Web app, sessions**: a dataset without genes can be saved to a session;
  capture used to fail with `Dataset fingerprint requires the current var
  field inventory`. Such a dataset now has an empty gene inventory rather than
  none.
- **Web app, camera**: at a low frame rate (a large dataset with level of
  detail off) inertia, zoom and the grid fade now run in real time instead of
  slowing down below 20 fps; only a stall over a quarter second is cut short.
- **Web app, rendering**: switching to another window or tab no longer redraws
  the whole scene for a second; focus changes redraw only when they end a held
  key or a selection in progress.
- **Web app, touch**: on an iPhone or iPad, pressing and holding the plot no
  longer selects the whole page as text or opens the callout menu; the press
  reaches the viewer (Cell info's hold, Select mode, the camera). Safari reads
  only `-webkit-user-select`, so every non-selectable control now declares it
  too.
- **Web app, highlights**: a highlighted selection no longer keeps the whole
  scene redrawing every frame in `Light` and `Ultra-light` shader quality,
  where the ring does not animate; the view goes idle as it does without a
  selection. `Full` still pulses its ring.
- **Web app, fog**: points, highlight rings, connectivity edges, the orbit
  anchor and the velocity overlay share one fog model, evaluated once per
  point instead of per pixel, so every layer of a scene fades at the same rate.
  Fog now thins a point toward whatever is behind it instead of tinting it
  with a fixed fog colour, so fogged points fade correctly over a grid or a
  vignette; over a plain background the picture is unchanged. With fog off a
  dense 10M-cell view draws about 8% faster (Apple M1 Pro).
- **Web app, 2D views**: point edges in 1D and 2D views are smooth (a
  one-pixel antialiased rim) in `Full` and `Light`, instead of a jagged
  staircase that the antialiasing setting could not smooth; `Full` draws 2D
  views as cheaply as `Light`, and figure export draws them the same way.
- **Web app, rendering**: connectivity edges take about a third less GPU time
  (all 6.3M edges of the Suo sample: about 190 to 115–128 ms a frame) and the
  volumetric smoke cloud about 15% less (17.2 to 14.5 ms), on an Apple M1 Pro.
  Every shader now pins a vertex attribute to location 0, which Firefox
  otherwise emulates per vertex. Projectiles no longer allocate new vertex
  arrays every frame.
- **Web app, projectiles**: projectiles bounce off the **Studio** floor under a
  3D view, as they bounce off the grid's walls, and land on it instead of
  falling through. A projectile that had come to rest against a cell no longer
  falls through the grid's floor either.
- **Web app, rendering**: after a change the viewer draws two more frames
  and goes idle, instead of at least twelve frames and a full second; an
  adaptive LOD still stepping toward its level keeps the scene drawing on its
  own account, so nothing is cut short. A fast display no longer draws 60–120
  identical frames after every click, and where a frame takes seconds (a
  software renderer) every change costs seconds less.
- **Web app, smoke**: generating the cloud's noise volumes no longer queues
  GPU work faster than the GPU finishes it. Each batch of slices waits for the
  previous one and shrinks while the GPU lags, so on a slow or software
  renderer the page stays responsive while smoke mode prepares (before, a
  click there could wait half a minute behind the queued noise work).
- **Web app, selection tools**: while a lasso, proximity or KNN selection is
  drawn, a pointer move that leaves the previewed cells unchanged no longer
  re-uploads the highlight or redraws the scene; only a change does.
- **Web app, sessions**: a session saved while a view was coloured by a
  continuous field could not be loaded ("Session UI control … does not
  exist"): it recorded the legend's range sliders, whose ids are made afresh
  on every page. Controls with page-made ids are no longer recorded; the
  colour range and filter they show are saved with the field, as before.
- **Web app, shader quality**: `Full` shading is a little brighter on the
  side of each cell away from the light (ambient 0.5, was 0.4); the lit side
  still shows the colour itself. Highlight rings and centroid markers use the
  same lighting.
- **Web app, Studio backgrounds**: the floor's rules and the dots are lighter
  and finer (major rules no thicker than minor ones), so the stage reads as a
  calm, clean surface behind the data.
- **Web app, rendering**: the plot is drawn at exactly the screen pixels it
  covers. Its drawing buffer used to be rounded down from a fractional layout
  width (a plot beside the sidebar is often half a pixel wider), and the
  browser then rescaled every frame by a fraction of a pixel, which softened
  one-pixel lines and dots.

## [0.9.2] - 2026-10-17

Version 0.9.2 is a patch release on 0.9.1.

### Fixed

- A refused upload or event now reaches the client as the status the server sent.
  Every early refusal of a POST -- an expired or unknown session-bundle request, bad
  viewer credentials, a wrong `Content-Type`, an oversized body, an unknown path, a
  foreign `Host` -- was answered before the request body was read, and the socket
  was then closed with that body unread. A socket closed with unread input is reset
  rather than shut down, so a client still sending or not yet reading lost the
  response: Windows reported `WinError 10053`, macOS and Linux a reset or a broken
  pipe. The server now drains a refused request's declared body, within the
  largest accepted upload and for at most five seconds, before it closes. Both
  servers share one `do_POST` that does this.

## [0.9.1] - 2026-10-17

Version 0.9.1 is the PyPI submission release.

### Added

- Scanpy's own `X_umap` is served as written: when no `X_umap_1d`/`X_umap_2d`/`X_umap_3d`
  key exists, a bare `X_umap` of 1, 2, or 3 columns is read at the dimension its column
  count states, without renaming or writing to the caller's object. A dimensional key
  still decides its own dimension and still takes precedence; any other width is refused
  with its shape named. The same rule governs `add_transition_drift_to_obsm` and the
  browser's h5ad and Zarr readers, so one file opens the same way everywhere.
- A long build reports its phases. The direct-AnnData server runs five reported steps, the
  adapter reports obs classification, embedding resolution and which key resolved it,
  vector-field discovery and the graph read, and the manifest and centroid pass is a step
  of its own instead of silence after a success line. `AnnDataAdapter` takes `quiet`.
- `GET /_cellucid/protocol` answers with the wire capabilities the installed package
  accepts, so a web build can ask instead of assuming.
- A notebook learns when a command it sent was refused, rather than reading `None` as
  success.
- Release-contract validation for package, citation, documentation, and downstream recipe
  metadata; installed wheel and sdist gates on Linux, macOS, and Windows before
  publication; reproducible sdist normalization.
- Antialiasing is a live setting, applied on the next frame instead of at context
  creation, and defaults from the cell count until the checkbox is touched.
- A dataset opens at a point size derived from its cell count; the slider now
  reaches 0.050 and every stored session keeps the size it recorded.
- The volumetric smoke render mode is marked `alpha` in the viewer.

### Changed

- **Web app**: the synthetic benchmark generator's custom point-count ceiling is
  raised from 50,000,000 to 250,000,000. The previous value was a declared product
  limit chosen well below the hardware's; the new one is the hardware's, sitting
  just inside the `MAX_TEXTURE_SIZE²` bound on the renderer's interaction-state
  texture (268,435,456 on a context reporting 16,384). A request that the machine
  cannot allocate now fails as an allocation error rather than being refused in
  advance, because the limit that actually binds is the machine's memory — 16 bytes
  per point, so 4 GB of source arrays at the new ceiling — and the control cannot
  know what that is. The point-count input, both validators and the worker boundary
  all read the one constant, so they cannot drift apart.
- The neighbor graph and the vector-field overlay are read only when asked for:
  `serve_connectivity=True` / `--connectivity` and `serve_vector_fields=True` /
  `--vector-fields` on the direct-AnnData path. Each is otherwise never touched,
  and a dataset serving neither publishes exactly what an object holding neither
  publishes. Both are read and validated in full before the socket binds, because
  the manifests the viewer fetches first declare what they contain — and on a
  large object that is the longest part of startup, for a capability most sessions
  never turn on. A 50-neighbor graph over 18 million cells is roughly 900 million
  stored neighbors. A graph asked for and absent is reported; vector fields asked
  for and undeclared serve none, which is ordinary; `vector_field_default` is
  refused unless the fields it selects among are served. The startup report names
  each state separately. `prepare()` is unchanged in Python and R — both already
  took `connectivities=` and `vector_fields=` as data — and now announces all
  three optional capabilities alike when they are not requested.
- Connectivity validation is measurably cheaper: degrees are counted with `np.bincount`
  rather than the unbuffered `np.add.at` (about twenty times faster on that step), one
  coordinate expansion feeds every check, and a graph already in compressed form is
  compared against its own transpose instead of being rebuilt first. Measured on a
  9M-stored-neighbor graph: 1.2x faster end to end and 36% lower peak memory. A graph
  above 50,000,000 stored neighbors logs what it is about to cost.
- The tested Python range is 3.11 through 3.14, and direct AnnData operation targets
  AnnData 0.12, Zarr 3, and numcodecs 0.16.
- PyPI publication uses GitHub OIDC trusted publishing and accepts only a pushed tag
  matching package version metadata exactly.
- Export format: every scientific payload is named by its integer index or a fixed
  neutral name, so no filename carries an observation key, a gene name, or a field id.
- The notebook viewers, the session applier, the data server, and the export writer are
  packages rather than single files. Every public name is imported from where it was, and
  the servers no longer import numpy at startup.
- Adaptive level of detail is bounded by a point budget rather than a fixed ratio.
  A cloud already inside the budget is never reduced, and one above it stops at the
  budget instead of falling to a 44x reduction on any pull-back: at 18.1M cells that
  is 2.4M points rather than 412,000.
- Served gene columns are cached under a byte budget rather than a count of a
  hundred columns, which on an 18.1M-cell object was 7.3 GB of resident memory.

### Fixed

- `cellucid serve --help` no longer fails with `UnicodeEncodeError` on a Windows console
  using a single-byte code page such as CP1252. argparse wrote help, usage, version and
  error text straight to the stream, bypassing the escaping every other console line
  already gets; the CLI's parsers now route that text through the same writer, so any
  character the console cannot show is printed as a backslash escape instead of
  aborting the command. The `serve` help text itself is now plain ASCII.
- A bind address is no longer printed as a URL. `--host 0.0.0.0` used to print
  `http://0.0.0.0:<port>` and try to open it; the banner now prints the loopback origin,
  which is what a browser on the serving machine and the near end of an SSH tunnel both
  use, plus the machine's own name for clients elsewhere. An IPv6 literal bind is
  bracketed as a URL requires.
- An import failure names the condition it is: a distribution absent from the environment
  (`pip install`), one present but too old to hold the name asked of it
  (`pip install --upgrade`), or a standard-library module the running interpreter cannot
  provide -- reported with that interpreter's version, its executable, and the supported
  range, because no package supplies it.
- `cellucid serve` answers an operator mistake with one actionable line instead of a
  traceback, and keeps the traceback for genuine defects.
- Both writers refuse the same numeric input everywhere it can enter. Float32
  underflow (a nonzero value below `2^-149` becoming `0.0`) is now refused for vector
  fields here, and for quantized fields, embeddings, and `latent_space` in the R
  writer, where normalization hid the caller's magnitude from the write-time check.
- A label the browser cannot draw, including an unpaired surrogate, is refused by name.
- `prepare()` reconciles the export root, so it cannot publish coordinates the viewer is
  never told to fetch, and every manifest is checked as written bytes rather than only as
  an in-memory payload.
- Every published payload is little-endian whatever the host is.
- `prepare()` and the web cache both refuse a directory nobody set aside for them,
  including the working directory, the home directory, and paths reached through a link.
- A constant continuous field -- a gene at one level in every cell, or a column a subset
  flattened -- has an encoding instead of failing the export.
- `debug_connection()` reports the dataset probes for a published sample state.
- Publishing a sample's starting view restates its identity. A view is captured
  against a locally served copy of the generation, and the catalog serves the
  same bytes under its own id, so a capture published as-is was refused on every
  open and the sample showed none of its view. The publish step now renames the
  capture -- only after re-deriving the generation's counts and cell-order digest
  from what it ships -- and repoints the catalog digest itself, so the artifact
  and its advertisement cannot part company.
- A refusal that names two identical cell and gene counts as proof of a
  different dataset now names the two dataset identities instead.
- `prepare()` no longer rejects an export because of a gene it was never asked to export.
- A payload that is not entirely finite is refused with counts -- how many NaN,
  infinite, or outside float32, and the first positions -- instead of a bare `500`.
  The server answers `422` with that diagnosis as JSON; the exporter and the R writer
  raise the same report, and the viewer shows it.
- An official sample's published view applies over plain HTTP. Its integrity check used
  `crypto.subtle`, which is secure-context only, so a LAN bind opened the dataset
  uncoloured; the digest is computed the same way with or without WebCrypto.
- `Auto` level of detail follows the camera at every dataset size. Its floor was an
  absolute point budget, unsatisfiable below two million points, so on most datasets
  `Auto` answered full detail at every distance and did nothing; the floor is now the
  smaller of that budget and an 8x reduction cap.
- Level of detail no longer draws square patches. Sorting on a bit-reversed Morton
  code made each level a set of whole lattice cells; sorting on the plain code and
  bit-reversing the rank takes the same fraction from every neighbourhood. At 44x
  reduction: density error 138% to 5%, empty populated regions 42% to none.
- The velocity overlay's `Sync with LOD` read the level scale backwards, running at
  full particle count on the coarsest level and disposing itself at full detail.

### Security

- Jupyter events, session uploads, web asset inventories, and server artifacts use
  authenticated schemas, bounded payloads, and path confinement.
- Both servers validate the `Host` header ahead of routing, file reads, and event
  delivery. A non-loopback bind is an explicit request for network exposure and has no
  authentication of any kind: every user who can reach the port reads the dataset.

### Documentation

- The Python, R, and web guides are organized around the current executable contracts,
  with real browser screenshots for primary workflows.
- Serving from an HPC cluster is documented end to end, including why a compute node
  commonly refuses SSH, the login-node port forward that works instead, why that requires
  `--host 0.0.0.0`, and why the browser address stays `127.0.0.1`.
- Corrected embedding guidance: an all-identical embedding has no finite normalization
  range and is rejected.
- Added a complete standard scVelo recipe for vector fields.

## [0.0.9] - 2026-01-01

This release graduates Cellucid Python out of alpha (still pre-1.0).
The GitHub release was tagged `v0.9.0`, while its Python distribution metadata
declared `0.0.9`. Release 0.9.1 restores one exact version across both systems.

### Added

- AnnData-first viewing and serving with `.h5ad` and `.zarr` support (lazy/backed loading where possible).
- Unified CLI serving with auto-detection for AnnData files, Zarr stores, and
  pre-exported dataset directories.
- Jupyter notebook integration (`CellucidViewer`, `AnnDataViewer`, `show()`, `show_anndata()`) with event hooks and session export.
- Session bundle support (`.cellucid-session`) including `CellucidSessionBundle` and `apply_cellucid_session_to_anndata()` for round-tripping highlights and user-defined fields back into AnnData.
- Multi-dimensional embedding exports (1D/2D/3D) and vector-field overlays (RNA velocity / drift) via `prepare()` and `vector_fields` helpers.
- Hosted web UI proxy mode with on-disk caching helpers (`get_web_cache_dir()`, `clear_web_cache()`).

### Changed

- Export format now includes explicit dataset identity metadata (`dataset_identity.json`) for reproducible sharing and session compatibility checks.
- Reduced export size and improved load performance with optimized manifests, improved connectivity edge export, and optional quantization + gzip compression knobs.

### Security

- Session bundles are treated as untrusted input with bounds checks and dataset mismatch policies when applying to AnnData.

### Documentation

- Major Read the Docs expansion and restructuring (Python package + web app guides), plus new publishing and contributing documentation.

## [0.0.1a0] - 2025

### Added

- Initial alpha release
- AnnData visualization with UMAP embeddings (1D, 2D, 3D)
- Gene expression overlays with sparse matrix support
- Cell metadata coloring (categorical and continuous)
- Interactive filtering and cell selection
- KNN connectivity visualization
- Multiple deployment modes:
  - Local demo (browser-only)
  - Browser file picker (h5ad and exported formats)
  - Server CLI (`cellucid serve`)
  - Python API (`serve()`, `serve_anndata()`)
  - Jupyter integration (`show()`, `show_anndata()`)
- Export functionality for web deployment

[0.9.8]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.8
[0.9.7]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.7
[0.9.6]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.6
[0.9.5]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.5
[0.9.4]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.4
[0.9.3]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.3
[0.9.2]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.2
[0.9.1]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.1
[0.0.9]: https://github.com/theislab/cellucid-python/releases/tag/v0.9.0
[0.0.1a0]: https://github.com/theislab/cellucid-python/releases/tag/v0.0.1a0
