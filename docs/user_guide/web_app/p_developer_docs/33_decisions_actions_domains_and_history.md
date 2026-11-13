# Decision record: actions, state domains and the history

This page answers: **"Why is every change of the open dataset one cataloged
action, why is state captured in domains, and why is the history a tree of
immutable steps over a pool of values, with soft deletes, declared inputs and
an autosave that keeps all of a write or none?"** It records the decisions
with their dates, the measurements that justify them, what was tried and
taken out, and what is still open. Every field, rule and how-to of the layer
is in {doc}`29_action_layer_domains_and_history`, and this page does not
repeat it. Update this page in the same change as the code it describes.

## At a glance

**Audience**
- Developers changing the action catalog, the dispatcher, a state domain, the
  history's core, its autosave or the session file, and anyone asking why a
  rule of that layer exists.

**Time**
- 20 minutes.

**Prerequisites**
- {doc}`29_action_layer_domains_and_history` (the reference) and
  {doc}`10_sessions_persistence_and_serialization` (the file a history is
  saved in).

:::{note}
Paths are under `cellucid/assets/js/` unless they name another place. A date
written "by 2026-11-12" is the day the entry was first written into these
docs, where the day of the decision itself was not recorded. The history
described here is its second design, of 2026-11-04; test files name it
(`history-v2-core.test.mjs`).
:::

---

## What it does now

1. **Every change is one action, defined once** (`app/actions/catalog.js`,
   `app/actions/definitions/`). The definitions import nothing from the page,
   so the catalog loads in Node and every other surface (the Python methods,
   the R functions, the agents' tools, the action reference) is generated
   from it.
2. **One dispatcher applies it** (`app/actions/dispatcher.js`): it checks the
   parameters, orders actions by the domains they touch, records one step
   when something changed, and puts the before-values back when a handler
   fails or is stopped.
3. **Some changes are recorded after the fact.** A settings control and a
   camera gesture change their owner first and are recorded by observation
   (`recordObserved`), synchronously from the owner's change event, so a
   change the history itself makes while showing another step is never a
   step; a gesture is recorded once it settles.
4. **State is sixteen domains** (`app/actions/domains.js`,
   `app/state/domains/`), each with one owner that captures its value,
   restores it and describes a difference. Values name things (a field by its
   key, a cell by its original position), so they apply to whatever cells are
   open and replay on another dataset.
5. **The history is trackers of step trees over a value pool**
   (`app/history/`). A step is an immutable record with the hashes of the
   values before and after it; a value is kept once, by the SHA-256 of its
   canonical form. Showing any step puts values back and computes nothing
   again.
6. **Branches, groups, notes and bookmarks are records beside the steps**,
   and every structural edit is journaled, so a delete is restorable from
   Recently deleted until the journal is emptied.
7. **Staleness is derived.** An action declares its result inputs; a step
   records their hashes; a step whose input differs later is stale, with the
   step that caused it.
8. **A history replays by name** on another dataset, as a new tracker, with a
   report of what applied.
9. **Every change is autosaved as it is made** (IndexedDB
   `cellucid-history`): only what changed, in one transaction that the
   browser runs whole or not at all. A session file (format 4) holds the same
   records.
10. **A person chooses what is recorded and shared**: the recording switch,
    the two detail switches of History settings, and a share of the history,
    the state only or the steps only.

## Decisions and why

### The design of the history

The owner asked on 2026-11-03 for a more rigorous action layer and history:
branches easy to manage, delete and see fork; notes; agents tied into the
history. The design was written and decided on 2026-11-04.

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-11-03 | Out of scope: merging branches (a three-way merge of values), exporting a history as a Python or R script, and searching for datasets inside the app | Decided on 2026-11-03 and kept out of the design | A merge of two branches' values |
| 2026-11-04 | A branch is a record with identity (`b<n>`: a unique name, its fork step, its tip, main, locked and archived), and every step records the branch it was made on | In the first design a branch was the name of its current last step: it could not be deleted except step by step, two branches could share a name (the agent and a view kept from a branch then refused both), and the panel showed only the branch shown, the others as a list that said where each left and nothing else | A branch as a label on a tip |
| 2026-11-04 (owner) | Deletes are soft: what is deleted goes to Recently deleted and is restorable until the person empties it; nothing expires by time | Lost analysis is the worst outcome, and git and Galaxy both keep a recovery path. The cost is memory for the kept values | A hard delete behind a clearer confirmation |
| 2026-11-04 | Steps are immutable records; branches, trackers and bookmarks are labels with identity; undo never destroys (selecting a parent shows it, acting from there branches) | Deleting a branch then removes only the steps no other branch needs, and those stay recoverable | Rewriting a step in place |
| 2026-11-04 | Notes, bookmarks and groups are metadata outside a step's identity, and editing them is journaled | Annotating a shared workflow never changes a value hash or invalidates a replay | Notes inside what a step hashes |
| 2026-11-04 | Loading a session file adds its trackers beside the person's own: a tracker whose id is taken gets another, a name taken is suffixed "(shared)" | Loading a colleague's workflow must never hide your own work | A load that replaces every tracker |
| 2026-11-04 (owner) | "Save AnnData with this state" writes a `.zarr.zip` only | Zarr is in the Cellucid-ready specification and a Zarr ZIP writer is small; an `.h5ad` would need an HDF5 writer in the browser, a large new module. Python and R convert in one line | An `.h5ad` written by the browser |
| 2026-11-04 (owner) | "Apply these steps here" is in scope (`history.branch.apply`): a branch's own steps replayed by name on top of the step shown, as a new branch | It is the natural way to try one branch's labels on top of another and costs little on top of replay. No values are merged; a step that does not apply is a flagged step | A three-way merge |
| 2026-11-04 (owner) | A recording switch: while paused, changes apply and make no step, and resuming keeps what differs as one squashed step | A person can work freely without a step for every change, and what was done meanwhile is still handled correctly | Dropping what was done while paused |
| 2026-11-04 (owner) | Consecutive steps can be grouped, and squashed into one step | A long run reads as one node; the paused period is shown the same way | — |
| 2026-11-04 (owner) | A state can be shared without its history (the state only), and a history without its values (the steps only, a replay plan) | Sharing is chosen, not total | A share that always carries every tracker |
| 2026-11-04 (owner) | History settings can leave two details out of the history: camera moves, and panels, windows and theme. Nothing else can be switched off | A person can stop recording trivial details; it must never reach the analysis state, and never clash with histories or branches: when a detail is off, its changes are simply not taken into the history | A switch for any domain that is part of the analysis |
| 2026-11-04 (owner) | No picture of a step is stored: not in the history, the autosave or a file. A picture is drawn from the live view when asked for (the printable report, `view.picture`) | Stored step pictures would bloat the history | A `picture` field on a step; a picture chunk in a session file |
| 2026-11-04 | Staleness is derived from the input hashes a step records, never stored, and values and results use one digest (SHA-256) | Analysis results had flagged themselves stale with a second digest, and nothing in the timeline said which step made a result stale. The same declarations make deleting a step apply again only the steps that depend on it | A stored stale flag; a second hashing path |
| 2026-11-04 | A step records the version of its action's parameters, and upgraders bring earlier parameters forward | A history recorded with older parameters keeps replaying | A file that stops replaying when an action's parameters change |
| 2026-11-04 | Session format 4 was changed in place for this design: no format 5 and no migration | Format 4 had not shipped in a release (its first is 0.9.7, of 2026-11-07) | Two readers for one format number |
| 2026-11-04 | The value hash rule stays as it is; a big sub-value two values share is stored once on the storage side only (`$part`) | Hashing a result once by its identity would change the rule the session file, the autosave and the Python and R readers verify | A hash that depends on how a value is stored |
| 2026-11-04 | Budgets are measured regressions, never caps: nothing limits a history's size | The stress suite measures what a history of a size costs; a number is there to be compared with the next run, not to refuse a person's work. The browser half of the suite prints its figures and asserts none; the Node half asserts what the suite returns (a few hundred steps in `npm test`), never a timing (`tests/history-stress-harness.test.mjs`) | A cap on steps, branches or bytes |

### The catalog, the dispatcher and recording

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-11-12 | The page selector of an analysis (`app/analysis/ui/shared/page-selector.js`) tells its listener whether a person's click made the change (`byPerson`), and only those become `analysis.pages.set` steps | The selector's own refresh while a reload restored an analysis was recorded as a person's step, so undoing `cells.open` closed the gene signature analysis and dropped its page | Recording every selection event as a person's step |
| 2026-11-10 (owner) | An input whose value a result holds itself is declared in `captures`: still recorded, and a later change of it never makes the step stale | A keyframe is the camera it was saved at: a camera move after it does not make the keyframe out of date | Marking a step stale for an input it carries a copy of |
| by 2026-11-12 | An action that may fetch an address its parameters name (`opensLinks`: `dataset.open`) asks, before anything is requested, about an address the person has not seen in this tab (the permission kind `unseen-link`) | An agent misled by text it read could otherwise make the tab send a request to any host | Opening an address only an agent has seen without asking |
| by 2026-11-12 | An action can say it keeps what it did when stopped (`keepsWhenStopped`: `camera-path.play`) | A stopped path's camera stays where the film stopped, as a pause leaves it; a rollback would jump it back | Rolling back every stopped action |
| by 2026-11-12 | The recorder that captured an action's before-values is the one told how the action ended | Another dataset opened meanwhile has a history of its own, which must never put back or record the replaced dataset's change; the replaced history, retired, does neither | Telling whichever recorder is installed when the action ends |
| by 2026-11-12 | Every wait of a handler ends or hears its signal, and hearing it means stopping the work it waits for | The dispatcher waits for a stopped change's handler to end, so a handler that never ends keeps `dispatcher.whenApplied()` (every history move, a save) waiting with it; work that ran on after the rollback would land on the state put back | A wait that only stops waiting |
| by 2026-11-12 | A change observed while a dispatched action on the same domains is applying is its own step when its caller says it is the person's own act (`{person: true}`), and a rollback of the running action leaves those domains as the person made them | A window closed with × while its run ran came back with the stopped run's rollback | Folding a person's act into another action's step |
| by 2026-11-12 | An input held per view is compared at the view the step acted in (`readPart`) | Keeping a view or changing another view's filters marked steps stale | Comparing a per-view input across every view |
| by 2026-11-12 | With several views, a step that changed a per-view domain in some views, not all, names them in its label | The label, the methods text and the pseudocode then say which views a step changed | A label that is the same for any view |
| by 2026-11-12 | A share of the cells is written by one rule in the web app (`formatShare`, `utils/format-units.js`) and in the history words of Python and R: whole from 10%, never 100% short of every cell, two significant figures below | Python and R said "6%" for 5.5% where the web app said "5.5%", and below 1% the web app's own two copies differed | A rounding rule of its own in each implementation |
| by 2026-11-12 | `annotation.publish` is one commit per call | Eight votes of one agent turn were eight commits of the person's file and eight questions to the person | One commit per vote |
| by 2026-11-12 | A parameter refused for its form says the form in words (`PATTERN_WORDS`, `app/actions/schema.js`); an optional parameter sent as null, "null" or {} (for "no value") is refused with how to fix it; an object value is written out; a field given as a name where `{source, key}` is asked, and the reverse, is told how to give it | "Does not have the expected form" left a person and an agent to guess what a keyframe label, a view name or a color must look like. Agents sent camera "null" (Haiku 5.5), camera {} (Opus 5.5), a guessed eye with yaw and pitch 0 and a field object, and the refusals said "[object Object]" (seen live, 2026-11-11); models sent the one form of a field for the other | A refusal that names no form |
| 2026-11-04 (owner) | Examples of every version of the parameters of every action whose steps are recorded are kept (`tests/fixtures/action-versions/<type>@<version>.json`, made by `scripts/generate-action-versions.mjs`) | A history recorded with older parameters keeps replaying | A version with no example to replay |
| by 2026-11-12 | The parameter versions in force: `camera-path.keyframe.add` and `camera-path.playback.set` v2 add `captions`; `field.create-from-pages` v2 adds `rest`; `highlight.page.delete` v2 and `highlight.category` v2 take several pages and several categories; `expression.transform.set` v2 spells `normalize` in US English like every catalog name; `figure-export.set` v2 (2026-11-11) never has a null `svg_strategy` and its `frame` takes "cells", v3 (2026-11-11) makes the figure physical (width in millimetres, no height, font sizes in points); the options of `analysis.differential.set`, `analysis.composition.set` and `analysis.signature.set` v2 (2026-11-11) add Pair samples across groups, conditions, Group by, Groups to score and Cells | Each older step keeps replaying through its upgrader (the row above) | Changing what a parameter means under the same version |
| by 2026-11-12 | `highlight.page.delete` takes several pages at once and `highlight.category` several categories | Clearing the pages an agent had made took 56 calls; naming several categories is one act, one step and one undo for a page of ten donors. A v1 step named a category by its raw value (a number or a boolean, which the legend sent and a boolean or numeric field refused) | One call per page or category |
| by 2026-11-12 | `field.create-from-pages` takes `rest`: the cells on no page keep their category of a categorical field | Splitting one label is then a page per new label and one call. An agent had made 45 temporary pages, a field of 48 pages, then deleted 56 pages: 101 calls, $3 | A page for every label of the field |
| by 2026-11-12 | A failure of a control's dispatch is about what was open when the control was used, so the next opening withdraws it (`app/actions/control-dispatch.js`) | A field's failure stayed through three datasets opened after it | A failure shown for a dataset that is gone |
| by 2026-11-12 | The golden history's apply plans and apply rules are written out in the packages' session documentation ({doc}`../../python_package/h_developer_docs/11_hooks_events_protocol_and_schema`), and Python, R and the browser's AnnData copy compare with the golden's rows field by field (`scripts/generate-history-golden.mjs`) | The rules are written once, where the packages that read them are documented | A rule kept only in the generator |
| by 2026-11-12 | A parameter's description is left to the schema's `const` where it says it | A description repeated the field's own in every tool, about 230 tokens a call | A restated constant |

### Persistence

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-11-04 | The autosave writes only what a change names (steps, turns, journal entries, messages, the tracker's own record), each value once by its hash and each cell set in the cell-set encoding | In the first design every step rewrote every chat transcript and walked the whole value pool, and IndexedDB held a brush of about 1M cells as 4.0 MB where the session file held it in 2.7 KB (Measurements, below) | A write whose cost grows with the whole history |
| 2026-11-04 | One page writes a saved state (a Web Lock, and the store's writer check): a second tab continuing it opens a copy | A saved state continued in two tabs would otherwise be written by both | Two writers of one saved state |
| 2026-11-10 (owner) | A starting view leaves the viewer background as shown (`STARTING_VIEW_SHOWN_CONTROL_IDS`); a background never chosen follows the theme (`backgroundForTheme`) | The restore had dispatched the background's `change`, which kept Studio (light) as the person's pick on every first visit, so the background never followed the theme again | A starting view that sets the background |
| by 2026-11-06 | A transaction of `utils/indexed-db-database.js` reads first, then makes every write at once, none awaited, and is told to commit (`IDBTransaction.commit`); a read after a write throws | The browser commits a transaction once it has no request left to run, and a page that goes away while its code waits on one request never makes the next. Writes made one at a time, each awaited, were therefore committed as far as they had got: Firefox kept 98 and 252 of 400 such writes of one transaction when its tab closed (3 closes in 30), and a saved state came back with a conversation's count and not its messages ("This saved state is missing messages of a conversation; it cannot be opened", hosted macOS Firefox, on a hosted runner) | An awaited write inside a transaction that must be whole |
| by 2026-11-12 | A group edit reaches the history kept in this browser: the autosave writes the grouped steps again | A group made, extended, removed or restored changes its steps' group in place, and the journal records the runs, never the member steps; the history's `edited` event named only the steps an edit kept, so a step made later rewrote its parent with the group and left the steps before it without one, and Continue after a reload refused the whole history ("a group cuts the run of turn ... between steps s5 and s6"): an agent's two branch steps grouped as "Hormone-only comparison" on Pancreas | An edit event that names only the kept steps |
| 2026-11-04 | Each chunk of a session file is gzip'd where gzip saves at least a sixteenth, decided per chunk by measured gain, never for the whole file | It pays for skeletons (they inflate about 6.8x for the lazy values chunk) and little for float arrays (1.27x) and cell sets (1.37x) (Measurements, "Persistence by size") | One codec for every chunk |
| 2026-11-04 | Analysis result floats stay float64: the differential expression statistics are not stored as float32; the gene names are stored once through `$part` instead | Float32 would change the values the hash names (a lossy write), which is the analysis code's to decide; float64 statistics gzip little (1.6 MB raw to 1.12 MB per run). Not done | A lossy write of a result |
| by 2026-11-12 | A part's reference is 76 characters and a part costs its name once more where it is stored, so a sub-value below `PART_MIN_CHARS` (512) is not stored as a part | Below it a part gains little when it repeats and costs most of its size when it does not (`app/history/value-parts.js`) | A part for every sub-value |
| by 2026-11-12 | The cell-set codec (`session/codecs/cell-set.js`) picks, for each set, the smallest of three encodings (sorted deltas, bitset, runs) and never stores a set relative to another | The Python and R readers decode these three; on Suo (562k cells, 22 highlight groups of 5-170k cells) the per-set choice stores 188 KB where sorted deltas alone stored 197 KB, both gzipped | A set stored relative to another |

### Stopping work and the history's moves

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-11-11 | A control of work under way (the `control` effect: `camera-path.stop`) records no step and touches no domain (`isSteplessEffect`, `app/actions/catalog.js`) | Declared on the camera domain it waited in order behind the `camera-path.play` it was sent to stop and was refused "not playing"; a stop of a play that waits for the path's end now ends it well before the path's length | A control ordered behind the work it stops |
| by 2026-11-12 | A looping camera path stopped in the frame of its wrap records no rest step | The camera is exactly at the first keyframe, where the loop's step already holds it | A step for a rest that changed nothing |
| by 2026-11-12 | A dispatch the person stopped (Stop, the page closing) while a handler ran on to its end ends as a cancellation all the same, and what it changed is put back; one that a newer choice of its type replaced and that still finished keeps its step | "Stopped" must never leave a change or a step, and a handler that never watched its signal had left both. The newer choice runs right after the replaced one, and putting back a domain it is about to replace (the open cells) is work for nothing | A stopped action that leaves a change |
| by 2026-11-12 | Work an action starts that goes on after it answers (`ctx.ongoing`: a looping camera path) belongs to the turn that started it and is stopped with it | A film an agent set looping outlived the turn the person stopped | Ongoing work with no owner |
| by 2026-11-12 | An exclusive history task (a move, a share, an edit) or a session operation stopped while it waits for its turn never runs, and its caller hears so at once; stopped while its restore runs, what was shown comes back and the head stays | An agent's `history.undo` that the person stopped while it waited ran anyway once its turn came, a move stopped while it put a state back went on to the end, and a stopped agent's share would have started after Stop | A move or an edit that Stop does not reach |
| by 2026-11-12 | `history.replay` hears the Stop at every wait (reading the saved state, waiting for its turn, between steps) and the exclusive lock is free however it ends | In the every-action sweep with live agents `history.replay` never answered, the Stop did not end it and the undo after it waited for ever | A replay that holds the lock |
| by 2026-11-12 | The signal goes down to the reading an action waits for (reopening the open dataset, keeping its views, the cells domain's reads of genes), and that reading is cancelled as the Dataset box's Cancel cancels it | The dispatcher waits for a stopped action to end before it puts back what it changed, so a reading not tied to the signal held every later history move for as long as the dataset was read (a minute on an atlas); a stopped markers run on Garcia ended with "The history has no value ..." when Pancreas opened | A reading that outlives its Stop |
| 2026-10-28 | One way to stop work (`utils/abort.js`): a cancellation is an AbortError with a kind (`user`, `superseded`, `closed`, `released`, `page-left`), made by `cancelled()` and `cancel()` and told apart by `isCancellation()` and `cancellationKind()` alone, and one brand check (`isAbortSignal`) says what a signal is. A grep guard (`tests/cancellation-guard-contract.test.mjs`) forbids the old forms outside that module | About 45 private copies of `throwIfAborted`, `abortError` and `isAbort` disagreed on what they threw, 13 subsystem codes each had a predicate that checked only its own, and an abort with a plain Error or no reason was shown as a failure by whichever check it passed through. Eleven validators said "is this an AbortSignal" in different ways, and the duck-typed ones (gzip, cell sets, bundles, community annotation) took any object with an `aborted` boolean for a signal | Private abort helpers; a plain-Error abort; a duck-typed signal |
| by 2026-11-12 | The guided tour (`history.tour.start`) ends at once when someone else changes the history, and does not return to the step shown before it began | Anyone else's step, move or edit while it played was overwritten by the tour's next jump (the person's coloring vanished from view, kept on a branch nobody was told of) and the tour's return jumped away from it again; an agent's `history_tour_start` ended the tour it had just started until the history told of the chat before the jump moved the head | A tour that overwrites the person's work |
| by 2026-11-12 | The selection-commit counter of the data source does not move when the shown dataset grows or is cut in place (`data/data-source-manager.js`) | That growth runs until another dataset is published and must never invalidate the change that replaces it: a link opened while Suo grew failed "superseded before publication" and Suo stayed | Counting resizes with selections |

### Steps, staleness and verification

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-11-04 | Ctrl/Cmd+Z with the focus on a slider, a checkbox or a select undoes the step they made; only text fields keep their own undo (the History tip: "except while you type in a field") | A slider, a checkbox and a select have no undo of their own to keep | Ctrl+Z reserved for the control that has the focus |
| by 2026-11-12 | Staleness follows the inputs that actually changed, and each reason names its own input's cause | A test read "stale since step 4 (Showed UMAP · 2D in View 1): Filters changed ... Hidden cells changed ... Views changed" when step 4 only kept a view: every reason cited the first change of any input, and keeping a view changes the per-view domains' values though no cell the test used changed | One cause for every reason |
| 2026-11-10 (owner) | A step whose name only resolves in a domain is never stale when that domain changes | Building a path (save, move the camera, save) must mark no keyframe stale; the keyframe's own camera is the `captures` rule above | Staleness from a domain the step only names things in |
| 2026-11-10 | A step that changes nothing is not recorded | The outlier filter moved 100 → 99 → 100 in one drag left "Showed every outlier of cell_type" (touching nothing) in the history, and Undo then changed nothing on screen (seen on Suo) | A step for a drag that ended where it began |
| by 2026-11-12 | A turn's step after its group joins the group by itself, and extending the group to it is not refused | It was refused "Steps 2 to 4 are not consecutive steps of one branch", which they were (seen in the live sweep) | A refusal of consecutive steps |
| 2026-11-11 | Verify compares the analyses in a canonical form that leaves out what a result records about its run (durations, `computedAt`, the reads' tuning, whether it came from the cache) | Verify said a colleague's Figure 2 branch was "not reproduced: Analysis differential: result changed" at the step that ran the test, though running it again gave every number the same: the duration was 627 ms in the sender's browser and 82 ms in the replay, and the analysis domain had no canonical form. It then said it of every branch with marker genes | Comparing a result's record of its run |
| by 2026-11-12 | Verify renames a user field's generated id in every domain, not only in `field-overlays` | A step that makes a user field gives it a new generated id each time it runs, and other domains name the field by it, so every branch that made a field and colored or filtered it was reported "not reproduced", and a delete applied such steps again for an input that had not changed | A rename in one domain only |
| by 2026-11-12 | `history.summary` carries each note in full, not only a count | A colleague's agent opening a shared session read "notes: 1" on each branch, had no way to read the text, and said the donor check the notes reported "is not recorded" | A count in place of the text |
| by 2026-11-12 | `history.diff` writes a value as its owner's readout does (`ui/core/control-readout.js`): each owner defines the readout of any value and uses it for its own readout | `history.diff` said "Grid density: 37" and "Render mode: smoke" where the sidebar and the step's label read "48³" and "Volumetric smoke cloud (beta)" | A second wording of a value |
| by 2026-11-12 | The phone sheet is page chrome, not state: on the phone layout the toggle, the peek bar and the H key open and close the sheet and record nothing, and every way to a control in the sidebar goes through `revealInSidebar` (sidebar first, then sections and focus) | A panel opened in the same moment as the sidebar's action reached the history as a change made outside it | A page-chrome change recorded as a step |
| by 2026-11-12 | A camera path's playback the person started, or a loop, is recorded once where it comes to rest; while it plays, nothing it does is recorded as the person's camera move | Each undo or turn during a loop recorded the moving camera as "Moved the camera", and an undo then undid that step and branched the history on the next (seen in use) | A step for every frame of a film |
| by 2026-11-12 | A saved state is named by its name and an autosave by its latest step in every list (Load..., Stored, My datasets); the panel's Opening message says nothing when the page shows the outcome | Autosaves had been "History kept in this browser" or "Autosave" alike in two of the three lists, so eight of them read the same; the message said "Opened" whatever happened (2026-11-11) | One name for every autosave |

### Restoring, sessions and saved states

| Date | Decision | Why | Rules out |
|---|---|---|---|
| by 2026-11-12 | A history another version of Cellucid wrote is restored part by part: a domain whose owner refuses its value keeps what is shown (and a part refused half way is put back to what it showed), a tracker this version cannot read is left out (the shown one kept as a new tracker), a state domain this version lacks is left out of every step. The person is told it was saved by an earlier or newer version, never that the file is damaged or that something broke inside Cellucid. A history of a saved state that could not be restored whole is never written into | A history kept minutes before the rename of the expression values (`expression` held `{normalise, log1p}`) was continued, its owner refused it and the whole restore was rolled back: the views, the steps and the conversation were lost with the one part this version could not read. A part refused half way: the slider was set to the saved 100, its renderer refused, and the notice said the settings were left as shown while the slider showed 100 | One unreadable part that fails the whole load |
| by 2026-11-12 | A history that cannot be loaded is said in one sentence, with its ids under Details | After a reload Continue failed as "Loading failed. Saved state: tracker 4f3c...: a group cuts the run of turn ... between steps "s5" and "s6"; a group holds whole runs of a turn": the model's own check, with the tracker's uuid, as the sentence | The model's check as the message |
| 2026-11-10 | Loading a session file of another dataset offers no "Continue where you left off?" in the middle of the load, and the welcome comes first; Load... tells autosaves apart by their latest step | The load opened that dataset, which offered it, and the offer went once the session landed; nine rows all read "Autosave [Auto]" | An offer made while a load is under way |
| by 2026-11-12 | A guard asks (Save..., Open ...) before the dataset shown is replaced or the page is left while steps made since the last save are not being kept (Autosave off, no storage, no room) (`app/unkept-work-guard.js`) | With Autosave off, opening another dataset dropped every step of the dataset shown without a word (a reload or a closed tab did the same) | A silent drop of unkept steps |
| 2026-11-10 | One notice announces a restore, says what it added beside the trackers shown and offers Undo; the progress notice gives way to it. A sample's own starting view restores without a word and with no transfer line. A session's notice names the dataset as the Dataset box does, never by its id ("from the sample dataset "suo"") | A successful restore published "Session fully restored." twice (its progress entry and a second toast) and Load State then added "Session loaded successfully", then a third announcement and the offer to undo as a second notice of it; a failed starting view was announced twice (its progress entry and "Sample view unavailable"); the starting view's own line read "Loading Starting view 0 B/s" at every sample's start, for a transfer that was none (`restoreLine`, `session/session-serializer.js`) | Several notices for one load |
| by 2026-11-12 | A session restores across engines: the camera's two halves (orbit and free-fly) are derived from each other through `Math.sin`, `Math.atan2` and `Math.asin`, whose last bits differ by engine, so only the half the navigation mode drives is compared | The multiview restore compared the derived half bit for bit (a WebKit pitch of -0.16539816339745403 against -0.16539816339745406 in Chromium and Firefox), refused the camera and re-asserted the same unreachable camera from its rollback: the He and Kanemaru sample views (saved in WebKit) opened with "Multiview restoration failed and graph rollback was incomplete." in every engine but WebKit, and every engine refused a file saved by another | A bit-for-bit comparison of a derived value |
| by 2026-11-12 | The checked-in published sample presets are compared with the settings control inventory by a contract test (`tests/session-preset-control-inventory-contract.test.mjs`) | A preset records the complete inventory, so adding or removing one control of the sidebar invalidated every preset saved before: the settings domain refused the value ("Session UI state is missing current control ...") and the starting view failed as a whole, the advertised camera never applying and no field colored. The advertised-default browser spec builds its bundle at runtime and could not see a stale checked-in file; the failure surfaced only when someone opened a sample | A stale preset found by a visitor |
| 2026-11-11 | A state domain's record is read through one declaration of required, optional and added keys; a key that is no longer written is not read (no dropped keys) (`state/domains/records.js`, `record-keys.js`) | A reader that demanded the exact key set of today refused every sample's starting view and every older session the moment a key was added ("Camera Path state must contain exactly autoplay, captions, ..." on every opening) | An exact key set |
| 2026-11-06 | The sizes remembered per analysis moved from the `analysis` domain to the `layout` domain: what an analysis window shows is the analysis domain's; where it stands, its size, whether it is folded and the sizes a person gave each analysis's windows and expanded views are `layout.analysis`. A saved state of before that holds the sizes in the old place is no longer read: the state domains read no dropped `sizes` key of analysis state (there are no users yet, so no read-time compatibility) | History settings > "Panels, windows and theme" then covers them as every other panel's geometry | Two homes for a window's geometry |
| 2026-11-03 | The theme belongs to the history (layout, panels, docking and theme are steps); page chrome (the phone sheet, the sidebar's scroll and focus, which panel is in front) does not | — | A step for page chrome |
| 2026-11-04 | The active selection tool (the highlight toolbelt: Annotation based, KNN drag, Proximity drag, Lasso; the pressed button of a `data-state-serializer-pressed-group`) is page chrome: kept in a session, a saved state and the autosave beside the history, written at a save and the autosave (`capturePageChrome`) and put back after the head is shown (`applyPageChrome`), never a step, never moved by undo or a jump, and held by no domain | A session saved in Lasso and reopened in Annotation would make the next drag do what the person did not ask for | A selection tool in the history |
| by 2026-11-12 | A saved record of controls is read one control at a time: a value one version of Cellucid cannot read keeps the control as shown and is said, and never fails the others (`readSavedControls`, `state-serializer/ui-controls.js`) | A session's every setting used to fail on one | A whole record that fails on one control |
| by 2026-11-12 | A control whose choices follow the embedding (the vector-field select) is checked and applied only after the embeddings and the active view are back, and its choice is applied before the switch that shows it | Checked earlier, a session saved on one embedding was refused on a page showing another; switched on first, the overlay showed the choice the page had and the restored one had to supersede it | A control restored before what its choices follow |
| by 2026-11-12 | `#hp-antialias` is a device preference, never a session control | A session that carried it was the last writer of it, and a sample's published view turned antialiasing back on after the person turned it off | A device preference in a session |
| by 2026-11-12 | The four renderer controls (shader quality, adaptive LOD, forced LOD level, frustum culling) have an owner, `ui/modules/render-controls.js`, initialized before anything is restored; a renderer that rejects a restored setting fails the restore loudly | They had been wired inline at the far end of the bootstrap, so a restored session wrote their saved values into the panel and never published them to the viewer | A control wired after the restore |
| by 2026-11-12 | The shared restore path ends a slider as a person's change ends: `input`, then `change` (an agent's `setting.set`, a session restore, a history jump) | An agent's `setting.set` of the vector field's sliders read back the value before ("Intensity: set to 0.25" for 0.9), and Particle density never reached the viewer (readout 60K, 15K particles drawn): the owner applies the density only on `change` (seen live, 2026-11-11) | A restore that dispatches `input` alone |
| by 2026-11-12 | A restore visits views one by one, and what redraws the view list and the sidebar from the focused view waits for the visit's end (`visitingViews`, `state/domains/view-list.js`) | A frame inside the visit of a hidden live view threw "A hidden live view cannot be the active view" when a saved lung analysis with its live view removed was reopened | A redraw in the middle of a visit |
| by 2026-11-12 | A restore raises nothing: `ui/modules/view-controls.js` reads `activeViewId` and `liveViewHidden` from their owners on every use and keeps no copy | It kept private copies of state its owners had already moved: a restore of the views emptied the viewer's snapshot inventory at the start of the rebuild, and the badge render the rebuild scheduled threw an uncaught `RangeError` on a kept view that no longer existed, once per restore. The error mattered because `state-serializer/ui-controls.js` installs a capture-phase window `error` listener while it dispatches each control event, so a badge render reached synchronously from a restored control's listener became a restore failure naming an unrelated control | A copy of what an owner holds |
| 2026-11-11 | A field of any kind is shown on any view, and a choice that ends unshown is refused with a sentence (`color.by` on a kept view) | Every view holds its own copies of the field objects (`_cloneVarData`) and the gene search's entries held the objects of the view the list was built in, so coloring a kept view by a gene did nothing, `color.by` answered null, and a gene chosen on another view was "no longer on the gene list" | A field list tied to the view it was built in |
| by 2026-11-12 | A user-defined categorical's definition (categories, codes, centroids) belongs to the field, and every view's copy is the same field: an edit in place (deleting or merging a category, a consensus column rebuilt) replaces it in every view's copy, and a copy follows the field's deletion when its view is shown | The edit replaced the edited view's copy and the stored field only: a kept view switched to an embedding not yet loaded computed centroids into the stored field's table, which its copy no longer shared ("Normalized user-defined centroids of layout ... were not published"); an edit in a kept view gave the live view's copy the new categories under its old per-category colors ("Category colors must contain one RGB triplet per category"); a field deleted, restored or purged in one view stayed as it was in the others | Per-view copies of a definition |
| by 2026-11-12 | User-defined centroids are computed from the normalized coordinates; the restore's UI settlement rebuilds the field selectors before syncing them; saving a session on a dataset without genes needs no var inventory | Positions are normalized in place as each layout loads, so no raw coordinates are resident: centroids computed only from raw ones left a new field without centroids (no labels) and threw "Raw 2D positions are unavailable for centroid computation" when the view's dimension changed; the settlement synced before rebuilding, so restoring a session whose active field was user-defined failed with "Active categorical field is absent from its rendered selector"; a save on a dataset without genes failed for want of a var inventory | Centroids from raw coordinates |
| by 2026-11-12 | A kept view's fields share the category list and centroid tables of their source, so a rename and every re-read of the renames never write them in place | Restoring the renames in the focused view relabeled the other view's centroids under its unchanged categories, and focusing that view (the `views` restore) failed with "Category centroids of layout ... reference an unknown category identity", leaving seven domains changed and the rollback meeting the same error | Writing the category lists and centroid tables in place |
| by 2026-11-12 | The statistics line (`sr-only`) names the active field of the view being shown and is published from a signal, not from refresh lists | Two kept views carry two different active fields, so switching between them is an active-field change, and the switch path published the legend without the line (it read "score" over a plot colored by cell type), which only a screen-reader user could notice | A line republished from lists that drift |
| by 2026-11-12 | The `layout` domain's value reads a panel's size from its owner (`--pos-width`, `--pos-height`), not off the rendered box, and holds no stacking order; a restored panel is floated through the same entry point as a person's tear-off | `_accordion.css` forces `height: auto` on a collapsed floating panel, so collapsing one before saving threw away the size it reopens at; the stacking order came from a counter private to `dockable-accordions.js`, which a restore wrote past, so a panel torn off after a restore opened underneath the restored ones; a panel floated through the entry point analysis windows use leaves no placeholder, so docking it dropped it at the end of the sidebar | A size measured off a collapsed box |
| by 2026-11-12 | A key that comes from data stays the data's own key through the walkers that rebuild results and saved state (`setOwnDataProperty`) | A category, field or group named `__proto__` changed the output's prototype instead, so the entry vanished from what an agent was sent, from a saved session and from a figure's legend | A walker that assigns a data key with `=` |
| by 2026-11-12 | A minted DOM id is recognized by its form (`ui/components/dom-id.js`), so the page keeps nothing per id; ids are minted without `crypto.randomUUID` | A set of every minted id had grown for the page's life; `randomUUID` exists only in a secure context, and the Python package supports `--host 0.0.0.0` over plain-HTTP LAN addresses | A registry of minted ids |
| by 2026-11-12 | A sign-in window left open for minutes is the person's time, not a request's silence: waiting on the person stops the startup deadline, and it starts over, whole, once the person answers (`app/startup-phase.js`) | An open failed as "nothing arrived for 60 s" while its owner was still in Google's sign-in window | A deadline that counts the person's time |

### Filters and the cells shown

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-10-28 (owner) | Every filter applies until it is removed, whatever a view is colored by: a gene's range filter and a categorical field's outlier filter hide cells exactly as an obs range or category filter does; the Active filters list, the visibility pass and "Visible cells" analyses read the one definition (`state/managers/filter-manager.js`) | — | A filter that applies only while its field colors the view |

### What an agent reads of the dataset, the views and the history

| Date | Decision | Why | Rules out |
|---|---|---|---|
| by 2026-11-12 | The summary read (`field.summary`) takes `within`: the cells of some categories of another field, or a range, with the view left as it is | Agents hid the person's cells with filters to summarize a subset and then put the view back: 6 filter steps in an Opus turn for donors within one cell type, 4 in a Haiku turn for cell types per age range | Hiding cells to read a subset |
| by 2026-11-12 | The gene search takes `names` (which of these are genes) | An agent checking a marker panel against a dataset of a few thousand genes made one search per gene, 28 calls in one turn | A search per gene |
| by 2026-11-12 | What an agent is given of the dataset lists every embedding the data mapping offers with the obsm key and columns `mapping_set` and `view_embedding_set` take, the modalities of a MuData, and the per-cell vectors the overlay draws | An agent given labels only guessed keys; an agent said it had no tool that lists the modalities of a CITE-seq file; an agent told of no vectors said a scVelo sample had no velocity | An inventory of labels only |
| by 2026-11-12 | `view.picture` says what a picture cannot show: the highlight page's cells over the coloring, the preview's color-vision simulation, the figure's printed size in millimetres, the camera path's caption on screen and the vector field's flow as a still picture. It takes `source: 'export'` (the figure exactly as Export would download it, scaled into the box asked for; 2026-11-11) and `source: 'screen'` (the pane as the screen shows it, overlay included); a figure picture taken while an overlay is on says so, and how to see the particles | A newcomer was told the bright cyan patches were high NANOS3: they were the germ cells' highlight over a gene's color scale. An agent preparing Figure 2 could not see the figure it was about to export: the picture used the screen's background and grid, no title, 12 px text and no frame, so the agent told the person three times "my previews don't show the export settings, please check the file", and never saw a panel title cut off; an agent once told the person the canvas is 1800 x 800 for a 954 x 841 figure and did the pixel-to-point arithmetic wrong twice. Another turned the vector field overlay on, tuned it blind and left a flow blown to white | A picture that claims to show what it leaves out |
| by 2026-11-12 | `vector-field.summary` reads per category the speed, alignment, direction and where its arrows lead | An agent asked "where are the cells heading ... evidence from the data" answered "the tools don't give me the velocity numbers" | A flow only as a picture |
| by 2026-11-12 | `camera-path.play` answers `ended: 'reduced-motion'` with the seconds moved (0) under reduced motion, and `views.state` names the preference | It had answered as if it had played ({seconds: 18.3, ended: 'end'}) though the camera went to its last keyframe at once, and the agent told the person "it ran all 18.3 s" (seen live) | An answer that says a film ran when it did not |
| 2026-11-13 | The `analysis` and `layout` domains capture without waiting for the analysis code: before it arrives they are the analyses with nothing opened (`docked: null`, `sizes: null`), the value the code gives on arrival, and the docked panel makes no mode until one is opened or used (`comparison-module.js`); only a restore that shows analyses waits for the code (`startup-analysis-code-later.spec.mjs`) | The history's start (`startTracker`) awaited the analysis capture, which awaited the code, so with the code held back startup stayed at `positions` and actions waited at the history's gate: about 5 s past `positions` in every engine with the code held 6 s. The panel made Quick at creation, so its arrival changed `analysis` and would have been kept as a change outside the history | A start or a step that waits for code the viewer does not need |

## Measurements that justify it

These are the figures of the history stress suite, one table for each question
it answered: its conditions, the first design against the second, what a step
costs by history size and by kind, the bytes a step adds and where they go,
what broke first, a real page, the timeline against the graph, the Python and
R readers, persistence, and the fix pass. `cellucid/tests/stress/RESULTS.md`
says how to run the suite and holds no figures, and `npm run stress:history`
measures it again: add the new figures here, beside the ones they are compared
with, in the same change.

### Conditions of every run

- Machine: Apple M1 Pro (10 cores), 16 GB, macOS 26.7.1, Node v26.11.1. The browser is Chromium of Playwright 1.62.0, launched on the machine's GPU (Metal through ANGLE; renderer string `ANGLE (Apple, ANGLE Metal Renderer: Apple M1 Pro)`); frames run at 120 Hz, so 8.3 ms is a full-rate frame. Each heavy run ran alone among heavy jobs, other light work on the machine at the same time; timings within about 20% are noise. Node runs used `--expose-gc` and read memory after a forced collection; each mixed run was one process with `--max-old-space-size=7168` and a memory guard of 7 GiB that neither run reached
- The suite's page is synthetic: every one of the state domains has a value of the shape its owner captures, on a dataset of 1,000,000 cells and 20,000 genes; analysis results are made by the analysis code's own encoder from results of real shape (differential expression 20,000 gene rows; markers 20,000 genes x 12 groups x 8 statistics; correlation 2 pages x 50,000 cells). A real owner's own capture and restore work is not in the figures (`recorder` below means the history and the pool alone)
- The mix, per 1,000 actions: camera 250, brush of 1,000 cells 50, setting 150, color by field 70, color by gene 40, color map 20, category filter 70, range filter 40, label edit 50, user field with 1,000,000 codes 1, brush of 100,000 cells 5, brush of 500,000 scattered cells 0.5, brush of about 1,000,000 cells (all but 1,000) 0.5, highlight group toggled 20, cells hidden (10,000) 5, all cells shown 3, panel moved 50, section opened 50, analysis panel switched 10, differential expression 1, marker genes 0.5, correlation 1. Past 30 highlight groups the oldest is removed; a branch is made every N/101 steps by undoing 1 to 20 steps and acting; nine side trackers get 1% of the steps each; about 30% of camera steps coalesce into the previous one. At 100,000 steps: 100,816 to 100,824 steps, 111 branches, 10 trackers, an active path 92,247 to 92,446 steps deep
- The first design is measured at `main` of 2026-11-04 (the history of branches named by their last step); "second design" is the same suite after the history core, then persistence (both 2026-11-04) and a fix pass (2026-11-05), seed 1, same machine. Two seeds (1 and 2) were run for the first design, shown "seed 1 / seed 2"; the 1,000-step column depends on which rare kinds fell in the first thousand steps (seed 2 had a marker result early) and is noisy; at 10,000 steps sizes that depend on which large values were made differ by up to 55%; from 50,000 up the seeds agree within about 16%
- Noisy: every maximum and most p95 (garbage-collection pauses land in single steps: the recorder maximum of 87 to 91 ms at 100,000 steps is collection), the 1,000-step column, replay milliseconds per step. Not noisy: byte counts, character counts, the sizes at which things crash, and how the accounting, autosave, jump and delete costs grow

### The first design against the second

Measured on 2026-11-04 (Apple M1 Pro, 10 cores, 16 GB, Node v26.11.1; a mixed
history on 1,000,000 cells and 20,000 genes, seed 1; timings within about ±20%
are noise). The rest of the same runs follows in the tables below.

| What | 1k steps | 10k steps | 100k steps |
|---|---|---|---|
| Memory accounting per step, ms (a walk of the pool → a counter) | 3.2 → 0.0 | 46.5 → 0.0 | 796 → 0.0 |
| `history.get` of every branch (recursion → a loop) | ok → ok | stack overflow → ok | stack overflow → ok |
| Delete a step with 1,000 after it (every later step applied again → only its dependants) | 265 → 595 ms | 1,823 → 41 ms | 42,909 → 491 ms |
| Undo, p50 ms (no walk of the path) | 0.18 → 0.21 | 1.6 → 0.16 | 40.2 → 0.08 |
| `history.get`, characters (everything → a page of 200 steps) | 291,503 → 67,308 | 2,887,596 → 67,654 | 27,330,523 → 73,009 |
| Recorder, p50 / p95 ms, every step | 0.13 / 0.41 → 0.22 / 0.80 | 0.13 / 0.60 → 0.20 / 1.48 | 0.15 / 1.45 → 0.15 / 1.24 |
| Autosave write per step, p50 ms | 1.7 → 0.3 | 47 → 0.2 | 848 → 0.3 |
| Autosave in IndexedDB, MiB | 16 → 8.2 | 217 → 96 | 2,880 → 1,474 |
| Session file of the whole history, MiB | 2.5 → 3.3 | 66.8 → 35.5 | write fails → 743 |

### The first design at every size

First-design figures at 1,000, 10,000, 50,000 and 100,000 steps for both
seeds ("seed 1; seed 2"). The comparison above restates seed 1 of some rows
beside the second design.

| What | 1k steps | 10k | 50k | 100k |
|---|---|---|---|---|
| Recorder p50 / p95 ms, every step since the last checkpoint | 0.13 / 0.41; 0.11 / 0.39 | 0.13 / 0.60; 0.13 / 0.62 | 0.14 / 0.89; 0.14 / 0.92 | 0.15 / 1.45; 0.15 / 1.52 |
| Recorder maximum ms | 3; 17 | 20; 32 | 36; 34 | 91; 87 |
| Memory accounting per step, ms (a walk of the pool on the main thread) | 3.2; 4.0 | 46.5; 62.3 | 342; 333 | 796; 927 |
| Autosave write per step, p50 ms | 1.7; 3.0 | 47; 63 | 345; 335 | 848; 897 |
| Autosave in IndexedDB, MiB | 16; 24 | 217; 290 | 1,379; 1,301 | 2,880; 2,966 |
| Whole history written at once ("Save as", a resumed autosave), s | 0.06; 0.12 | 0.82; 1.00 | 5.6; 5.5 | 13.3; 13.5 |
| A saved state read back and loaded, ms | 94; 202 | 1,307; 1,574 | not measured | not measured |
| Jump to a random step, p50 / max ms | 0.6 / 0.7; 6.9 / 7.8 | 10.1 / 11.0; 9.2 / 10.9 | 22.4 / 27.8; 21.8 / 29.6 | 41.4 / 53.0; 42.4 / 84.6 |
| Undo / redo p50 ms | 0.18 / 0.27; 0.17 / 0.24 | 1.6 / 1.9; 1.4 / 1.8 | 14.1 / 20.0; 13.7 / 18.2 | 40.2 / 56.2; 40.3 / 57.7 |
| Delete a step with 10 / 100 / 1,000 after it, ms | 5 / 33 / 265; 36 / 157 / 563 | 87 / 290 / 1,823; 21 / 264 / 1,782 | 193 / 1,561 / 15,090; 177 / 1,443 / 14,991 | 515 / 4,817 / 42,909; 502 / 5,237 / 41,728 |
| Replay, ms per step | 0.18; 1.2 | 1.3; 0.8 | 0.58; 0.23 | 0.46; 0.42 |
| `history.get` of the shown path, characters | 291,503; 291,803 | 2,887,596; 2,902,814 | 13,651,519; 13,617,050 | 27,330,523; 27,237,749 |
| `history.get` with every branch | 294,651; ok | stack overflow; same | stack overflow; same | stack overflow; same |
| Session file, MiB | 2.5; 14.6 | 66.8; 102.7 | 504; 446 | write fails |
| Session write / load, s | 0.7 / 0.2; 1.1 / 0.4 | 7.7 / 2.5; 10.2 / 3.3 | 50.4 / 18.0; 46.2 / 17.4 | not measured |
| JS heap / ArrayBuffers after collection, MiB | 33 / 43; 75 / 80 | 235 / 325; 318 / 418 | 399 / 956; 394 / 863 | 866 / 1,864; 882 / 1,907 |
| Pool: distinct values / typed-array MiB | 1,098 / 13; 1,104 / 17 | 10,599 / 155; 10,509 / 197 | 51,986 / 943; 51,877 / 846 | 101,585 / 1,835; 101,411 / 1,869 |

The memory accounting and the autosave write both walked the whole pool after every step (the
accounting in `dataset-history.js`, the autosave through `pool.byteLength()` and a copy of every
stored step id), so together they passed one 60 Hz frame at about 2,000 steps and cost about 1.7 s
a step at 100,000. Recording itself stayed flat: hashing was over 90% of the recorder in every kind
whose recorder p50 was above 0.2 ms and 60 to 80% below.

### Recording cost by step kind

First design, recorder p50 / p95 ms, mixed history.

| Kind | 10k, seed 1 | 10k, seed 2 | 100k, seed 1 | 100k, seed 2 |
|---|---|---|---|---|
| camera | 0.05 / 0.07 | 0.05 / 0.08 | 0.05 / 0.08 | 0.05 / 0.09 |
| color by field or gene | 0.17 / 0.23 | 0.17 / 0.21 | 0.17 / 0.24 | 0.17 / 0.53 |
| category filter, range filter | 0.09 / 0.12, 0.14 / 0.19 | 0.09 / 0.14, 0.14 / 0.33 | 0.09 / 0.13, 0.14 / 0.20 | 0.09 / 0.16, 0.14 / 0.40 |
| brush of 1,000 cells, of 100,000 | 0.22 / 0.76, 0.46 / 1.13 | 0.21 / 0.82, 0.43 / 1.10 | 0.23 / 1.32, 0.45 / 1.54 | 0.25 / 1.30, 0.47 / 1.54 |
| brush of 500,000, of about 1,000,000 | 1.36 / 1.75, 2.35 / 2.37 | 1.36 / 1.62, 2.29 / 2.34 | 1.39 / 3.28, 2.37 / 2.74 | 1.41 / 2.91, 2.38 / 2.67 |
| label edit | 0.58 / 0.99 | 0.57 / 1.12 | 1.48 / 2.80 | 1.57 / 3.24 |
| user field (1,000,000 codes) | 1.09 / 1.33 | 1.09 / 2.76 | 2.09 / 3.08 | 2.13 / 3.50 |
| setting | 0.26 / 0.58 | 0.27 / 0.64 | 0.27 / 0.75 | 0.27 / 0.85 |
| panel moved, section opened | 0.07 / 0.09 | 0.07 / 0.12 | 0.07 / 0.10 | 0.07 / 0.14 |
| analysis panel switched | 11.97 / 13.65 | 12.37 / 13.58 | 13.36 / 15.06 | 13.36 / 17.20 |
| differential expression | 8.95 / 14.65 | 14.30 / 20.29 | 15.62 / 21.51 | 15.19 / 19.96 |
| marker genes | 18.54 / 20.35 | 18.68 / 32.06 | 19.82 / 38.71 | 19.90 / 36.79 |
| correlation | 13.09 / 14.33 | 13.67 / 14.89 | 14.15 / 16.71 | 14.36 / 18.49 |

A kind costs what the values it touches cost: the analysis kinds hash the whole `analysis` value
(canonical JSON of every result, 20,000 gene names included, about 1.5 MB) twice, before and after,
so switching the analysis panel, which changes one string, costs 12 to 13 ms; label edits and user
fields grow with the history because `field-overlays` holds every user field made so far, and one
label edit touches five domains; a brush costs its cells' bytes (SHA-256 of a 4 MB array is
2.3 ms). The camera, color, filter, setting and layout kinds are flat from 10k to 100k steps.

### Bytes a step of each kind adds

First-design figures are bytes a step (session file; step record / values / arrays of it; autosave
in IndexedDB); the second design is "→ " after (persistence, 2026-11-04). Caution on the second
design's camera, setting and panel rows: the run measured each kind with a person's origin and a
turn that a person's step does not record, so a kind that coalesces merged every change into the
warm-up step; "225 → 0" for camera and the other coalescing rows measured merges, not steps. The
run was corrected on 2026-11-11 (each kind alone as a notebook cell with a turn each, only appended
steps counted, `merged` must be 0, the turn record of about 270 bytes a step at 5,000 cells included
in `fileBytesPerStep`); the second-design column of the coalescing kinds has not been measured again.

| Kind | Session file | Step record / values / arrays (first design) | Autosave in IndexedDB | Recorder p50 ms | Autosave write p50 ms |
|---|---|---|---|---|---|
| camera | 225 → 0 (coalesced, see above) | 125 / 100 / 0 | 999 → 282 | 0.06 | 0.36 → 0.11 |
| color by field, by gene | 27, 116 → 32, 124 | 25 / 2 / 0, 72 / 44 / 0 | 501, 580 → 960, 1,057 | 0.10, 0.11 | 0.29, 0.30 → 0.09, 0.10 |
| category filter, range filter | 138, 182 → 226, 220 | 79 / 57 / 0, 108 / 74 / 0 | 1,422, 1,404 → 2,033, 1,946 | 0.07, 0.12 | 0.38, 0.37 → 0.14, 0.14 |
| brush of 1,000, of 100,000 cells | 2,413, 65,908 → 2,512, 66,009 | 121 / 239 / 1,746, 124 / 201 / 65,268 | 18,228, 407,975 → 16,656, 108,475 | 0.33, 0.46 | 0.90, 0.63 → 0.53, 1.68 |
| brush of 500,000 scattered, of about 1,000,000 | 125,687, 2,661 → 125,745, 2,774 | 125 / 182 / 125,062, 113 / 174 / 2,061 | 2,003,696, 3,999,997 → 129,551, 7,429 | 1.29, 3.32 | 0.79, 1.73 → 6.65, 8.12 |
| cells hidden (10,000 more) | 979 → 941 | 105 / 83 / 175 | 1,556,493 → 1,406 | 1.02 | 0.56 → 3.39 |
| label edit, user field (1,000,000 codes) | 169, 293,496 → 260, 293,570 | 82 / 87 / 0, 156 / 275 / 292,752 | 4,025, 1,009,124 → 4,488, 1,003,955 | 0.24, 0.93 | 0.54, 0.81 → 0.47, 0.82 |
| setting | 167 | 79 / 88 / 0 | 6,011 | 0.28 | 0.88 |
| panel moved, section opened | 139, 90 | 83 / 55 / 0, 56 / 33 / 0 | 1,288, 789 | 0.07, 0.06 | 0.38, 0.32 |
| analysis panel switched | 1,605 → 31 | 23 / 1,582 / 0 | 8,299 → 1,358 | 12.38 | 15.29 → 0.18 |
| differential expression, markers, correlation | 1,228,632, 7,811,814, 825,373 → 1,133,244, 7,712,512, 722,146 | 63 / 103,829 / 1,121,630, 64 / 103,875 / 7,678,163, 63 / 103,810 / 718,779 | 2,115,669, 10,115,645, 1,315,735 → 1,607,066, 9,630,817, 808,421 | 13.88, 20.36, 13.52 | 24.17, 27.40, 23.44 → 8.31, 14.96, 7.64 |

The session file is gzip'd chunks; "arrays" is the typed arrays each value names, raw or cell-set
encoded, gzip'd; the autosave of the first design held arrays raw and skeletons and step records
uncompressed. The analysis-panel row is small only because the panel was switched among three
values that repeat; a real session's analysis values do not repeat, and each new combination stores
the whole skeleton again (about 0.5 MB uncompressed, 0.1 MB gzip'd, mostly gene names). A brush's
step holds its cells in its parameters as well as in the group, as the same array stored once; in
the first design every parameter array was eager, so at 50,000 steps 65 MB of the 504 MB file had to
be read before the head was shown (the head's own value skeletons were 0.1 MB).
With persistence, encoding a cell set costs 6 to 8 ms a step for 500,000 to 1,000,000 cells, off the
recording path (the write runs in the next task).

### Where the bytes go in the first design

| Domain | Distinct values | Skeleton JSON | Arrays (each once) | Why |
|---|---|---|---|---|
| analysis | 719 | 339 MB | 965 MB | every value repeats every result's 20,000 gene names in its skeleton (about 0.47 MB each); the arrays are float64 and float32 statistics that gzip barely compresses |
| field-overlays | 6,816 | 252 MB | 125 MB | every value holds every user field and rename so far, so a label edit re-stores the list: O(edits x fields) |
| highlights | 16,720 | 106 MB | 648 MB | every value lists every live group (up to 30); arrays are shared by digest |
| settings | 16,275 | 86 MB | 0 | every value is the whole 120-control inventory (about 5.3 KB) re-stored by every slider step |
| filters | 11,377 | 18 MB | 0 | |
| layout, camera, colors | 9,711 / 22,374 / 7,233 | 6.7 / 5.2 / 5.2 MB | 0 | |
| visibility | 631 | 0.1 MB | 63 MB | |

The values held about 0.8 GB of skeleton JSON and 1.8 GB of arrays. In the session file the lazy
values chunk was one JSON string: 357 to 375 MB at 50,000 steps (seed 1, uncompressed), which
overflowed V8's string limit (about 512 MiB) between 50,000 and 100,000 steps, extrapolated near
70,000 to 75,000. A marker-genes result adds 7.8 MB to the file, a differential expression 1.2 MB,
a correlation 0.8 MB; a camera step 225 bytes, a setting 167, a label edit 169.

### What broke first in the first design, and at what size

| Size | What |
|---|---|
| about 400 steps (extrapolated from 1,000 and 10,000) | with the State panel open (as it is by default) a recorded step's timeline redraw passes one 60 Hz frame: 35 to 37 ms at 1,000 steps and 296 to 299 ms at 10,000; at 10,000 steps opening the panel held the page 660 to 750 ms, scrolling ran at about 5 fps, a slider drag at about 3 fps with recording and slowed by paint even without |
| about 510 steps on the shown path | `history.get` passes the connector's 150,000 characters (about 293 a step in both runs; with every branch at about 500 steps in all); at 100,000 steps one read was 27 MB |
| about 2,000 steps | the memory accounting plus the autosave write after each step pass one 60 Hz frame (16.7 ms): 4.9 to 7 ms at 1,000, 94 to 125 ms at 10,000, linear in the pool; at 10,000 steps a slider drag that records a step per frame could not keep up |
| about 5,000 steps deep | two recursions threw "Maximum call stack size exceeded" (Node's default stack; the browser's stack was not measured): `history.get` with every branch (4,997 steps on a bare chain, failing at every checkpoint from 10k up) and deleting a step with more than about 4,600 to 6,000 steps after it (5,036 on a bare chain); a failed delete changed nothing |
| 10,000 to 50,000 steps | undo and redo pass a frame (14 to 20 ms at 50,000, 40 to 58 ms at 100,000); delete a step costs about 15 to 50 ms per later step |
| between 50,000 and 100,000 steps | saving the session file threw `RangeError: Invalid string length`; the autosave kept working because it wrote one record per value |
| 100,000 steps | nothing else crashed; every step cost about 1.7 s of main thread, the autosave held 2.9 GB, and in Node the history held 0.87 GB of heap and 1.9 GB of arrays with no collectable garbage; the heap alone does not reach a tab's limit (about 4 GB in Chromium) but the arrays count against the tab's memory. The pool keeps every distinct value for the dataset's life, so memory growth is unbounded by design |

Moves were O(depth): a jump to a random step 0.6 ms at 1,000 steps and 41 to 42 ms at 100,000; undo
40 ms at 100,000 because every move re-walked the path from the root several times; deleting a step
applied each later step again at O(depth), about 50 ms a step at 100,000.

### In a real page: recording on against bypassed

The dataset: a 1,000,000-cell AnnData (Zarr ZIP) with 30 cell types, 20 genes (CSR) and UMAP; the
history is the page's own (real owners, real autosave in IndexedDB). Two runs through the heavy
queue (5.0 and 1.8 minutes). Bypassed means the dispatcher's recorder set to none
(`setRecorder(null)`, no product switch). Each figure is three alternated rounds each way, run 1
then run 2.

| Scenario | Recording on: frame p50 / p95 ms | Bypassed: frame p50 / p95 ms | Recording on: action ms | Bypassed: action ms |
|---|---|---|---|---|
| Camera gesture (120 pointer moves) | 8.3 / 8.7-8.9, 8.3 / 9.3 | 8.0-8.3 / 9.1-10.8, 8.3 / 9.1-9.3 | gesture 1,023-1,047, 1,025-1,034 | 1,013-1,053, 1,021-1,024 |
| Slider drag (point size, 120 frames) | 8.3 / 8.7-9.2, 8.3 / 9.3-9.4 | 8.3 / 8.6-9.1, 8.3 / 9.3 | 1,018-1,029, 1,019-1,021 | 1,015-1,022, 1,015-1,017 |
| Highlight of about 1,000,000 cells | | | dispatch 33-35, 41-48 | 25-37, 36-41 |
| Marker genes (`analysis.run`) | | | 190-193 first, then 27-41 | 24-54 |

No long task (over 50 ms) in any camera or slider window either way; camera and slider show no
difference beyond noise. The highlight is about 3 to 5 ms slower with recording (34.1 against
31.0 ms in run 1, 43.4 against 38.3 in run 2: within one run's spread; it agrees with the Node
figure of 2.4 ms hashing plus capture). The marker run cannot be compared: only its first run
(recording on) computed (190 ms, one 64 ms long task); every later run in both modes reused the
analysis's cached result.

### The timeline of the first design against the virtualised graph

First design, settings steps, two sliders alternating; run 1, run 2. The graph is a tracker of
10,000 steps and 100 branches (100 branches of 10 steps forking from the first 100 steps, a main
line to 10,000, a third of it camera moves), three runs.

| | Timeline, 1k steps | Timeline, 10k steps | Graph, 10k steps and 100 branches (3 runs) |
|---|---|---|---|
| Opening the panel until the step shown is drawn, ms | 67, 72 (long task 64, 69) | 663, 659 (long tasks 708, 705) | 28, 33, 32 |
| Rows in the page | every step | every step (10,000) | 24 |
| One more step until its row is drawn, ms | 37, 35 | 299, 296 (long task 374, 370) | 78, 82, 79, of which the panel 22 to 25 (the rest is the stand-in core validating the whole tracker per step) |
| One more step, the panel's script / layout / style ms (run 2) | 9 / 10 / 14 | 64 / 112 / 127 | script 14.4-14.5 |
| Opening, main thread script / layout / style ms (run 2) | 45 / 13 / 14 | 366 / 158 / 133 | |
| Scrolling 120 frames: frame p50 / p95 / max ms | 8.3 / 9.3 / 15.8, 8.3 / 8.9 / 9.3 | 199 / 675 / 1,108, 200 / 667 / 1,112 | 8.4 / 17.6 / 41.7, 8.4 / 17.6 / 33.3, 8.4 / 16.9 / 41.5 |
| Scrolling, main thread task / script / layout / style ms (run 2) | 263 / 3 / 0.3 / 3.5 | 31,427 / 309 / 72 / 28 | no long task |
| DOM nodes in the page | about 20,000 to 22,000 | about 120,000 to 160,000 | |

The timeline rebuilt the whole list on every render (`replaceChildren` with every row), so opening
and adding a step grew linearly with the rows; scrolling 10,000 rows ran at about 5 fps with under
2% of the time in script, layout or style (the rest inferred, not traced, to be paint and
compositing of about 100,000 more DOM nodes).

The graph keeps only the rows in view, a few above and
below, the focused row and the open details, so its cost is flat in the steps: its layout is 5 ms at
10,000 steps and 100 branches in Node. Its p95 scrolling frame is about two frames because each of
the 120 frames jumps about 2,100 px (the whole graph in two seconds), so every row in view is
replaced each frame; a person's scroll replaces a few. Folded runs (Camera x6) and collapsed stubs
make 9,100 rows for 10,000 steps.

A new step is laid out from the layout before (only the path's
last row and what follows): in Node, 40 steps appended at 10,000 cost 0.17 ms a step (median) of
layout against 6.6 ms laid out whole, equal row for row; in Chromium the panel's frame for one more
step (20 steps, the stand-in core's validation excluded) is 10.3 ms median, 12.3 p90, 13.1 max
against 22 to 25 ms before, the rest the browser's style and layout of the rows written. Formatting
times reuses one `Intl.DateTimeFormat` per format (a new one per row cost 4 ms a step). The
10,000-step run after the change: opening 31 ms, 24 rows, scrolling p50 / p95 15.7 / 17.3 ms on a
machine shared with other runs.

### A slider drag at 10,000 steps, and the page's memory

| | Recording on | Bypassed |
|---|---|---|
| Wall time for 120 frames, s | 35.8, 39.2, 37.9 | 23.1, 27.2, 24.9 |
| Frame p50 / p95 ms | 326-333 / 824-924 | 16.7-16.9 / 907-942 |
| Long tasks: count / total ms | 61-62 / 21,257-21,672 | 2 / 456-505 |
| Main thread task / script / layout / style, s | 23.2-23.6 / 4.5-4.7 / 6.5-6.6 / 7.4-7.5 | 2.5-2.8 / 0.66-0.69 / 0.06-0.08 / 1.27-1.32 |

Recording cost about 21 s of main thread over 120 inputs (each recorded step's timeline redraw,
300 to 370 ms; about 3 fps). Bypassed, the drag still took 23 to 27 s with 2.5 s of main thread: the
frames waited on painting a page that holds the 10,000-row list, so the visible history slowed the
page even when nothing was recorded.

Memory in the page for settings steps only (417 distinct values for 10,000 steps, so the floor of a
history's heap), run 1 / run 2: JS heap after collection 14.2 / 14.3 MB at the start, 16.7 / 16.6 at
1,000 steps, 26.3 / 26.0 at 10,000; `navigator.storage.estimate().usage` 20.4 / 20.4, 21.1 / 21.1,
24.5 / 24.6 MB; about 1.2 KB of heap and 0.4 KB of origin storage per settings step. Building the
steps with the panel closed cost 2.06 to 2.20 ms per settings step at 1,000 and 2.62 to 2.73 ms at
10,000 (`dispatch`, in the page); the same step in the Node suite is about 0.3 ms, the difference
being the real owners' capture and handler (the settings domain captures the sidebar's whole control
inventory).

### Python and R reading a session file

The time is `read_session` / `cellucid_read_session` alone (the R package loaded once). First design
(session file as the web app then wrote it): run 1 ran beside another job's Playwright run, run 2
alone.

| File | Bytes | Python, run 1 | Python, run 2 | R, run 1 | R, run 2 |
|---|---|---|---|---|---|
| golden (second-design shape), 13 steps | 23,710 | 7, 6, 5 | 20, 9, 7 | 382, 381, 41 | 361, 341, 38 |
| golden + 1,000 steps | 43,232 | 63, 63, 80 | 67, 59, 70 | 738, 734, 814 | 684, 692, 763 |
| golden + 10,000 steps | 209,291 | 674, 730, 717 | 625, 670, 677 | 7,350, 7,343, 7,417 | 6,958, 7,051, 7,033 |
| golden + 100,000 steps | 1,870,373 | 7,480, 8,071, 8,697 | 6,957, 7,476, 8,004 | 259,743, 261,665, 254,863 | 217,849, 237,363, 241,057 |
| the web app's own file then (first-design shape, 1,000 steps) | 168 KB | refused | refused | refused | refused |

After persistence the same files (one run, beside another job's Playwright run):

| File | Bytes | Python | R |
|---|---|---|---|
| golden, 13 steps | 27,085 | 14, 7, 8 | 513, 528, 55 |
| golden + 1,000 steps | 46,745 | 76, 89, 76 | 838, 851, 890 |
| golden + 10,000 steps | 213,106 | 834, 866, 864 | 8,933, 8,507, 7,894 |
| golden + 100,000 steps | 1,883,686 | 9,594, 10,289, 10,522 | 110,353, 102,554, 104,897 |
| the web app's own Save, 1,000 mixed steps | 235,662 | 581, 546, 544 | 3,859, 4,200, 4,321 |

Python is linear at about 70 to 80 microseconds a step (about 0.1 ms after persistence). R was
superlinear in the first design: about 0.7 ms a step at 1,000 to 10,000 steps but 2.2 to 2.6 ms at
100,000 (10x the steps took 31 to 37x the time; a list grown one element at a time or a lookup by
scanning would do it, not profiled); after persistence, 10x the steps takes 12x (0.85 s, 8.5 s,
105 s). The files grew by step records only (their values are the golden's 13 steps' values), so
these figures are the step-parsing cost alone.

### Persistence against the first design, by size

Rows of the first table of this section (autosave write, autosave in IndexedDB, session file size) are not repeated.

| What | 1k | 10k | 50k | 100k |
|---|---|---|---|---|
| Session write, s | 0.7 → 0.81 | 7.7 → 7.2 | 50.4 → 51.8 | write fails → 114.9 |
| Head shown (manifest and eager chunks read), s | → 0.08 | → 0.37 | → 1.42 | → 4.0 |
| Whole load (eager, then lazy, then applied), s | 0.2 → 0.29 | 2.5 → 2.49 | 18.0 → 19.8 | → 75.1 |
| Share, state only (head as one squashed step), MiB / s | → 1.9 / 0.62 | → 7.9 / 3.5 | → 28.4 / 12.2 | → 46.8 / 23.5 |
| Share, steps only (no values), MiB / s | → 0.20 / 0.06 | → 1.26 / 0.32 | → 6.0 / 1.9 | → 12.1 / 3.5 |
| Save with results recomputed, MiB / s | → 2.0 / 0.76 | → 15.9 / 5.3 | → 73.5 / 25.7 | → 143 / 55.2 |
| A named save ("Save as" to this browser, the whole history), s: first design → first run → second run | 0.06 → 0.15 → 0.12 | 0.82 → 1.55 → 1.11 | 5.6 → 15.2 → 18.2 | 13.3 → 30.2 → 19.4 |

Where the 100,000-step file's bytes are: steps 12.1 MB (15 eager pieces, 97 MB inflated), the head's
values 1.8 KB with 46.6 MB of the arrays it names, and lazily 625.9 MB of raw arrays (analysis
results and user-field codes), 62.9 MB of cell sets, 18.3 MB of parts (283 MB inflated: each
sub-value two values share is stored once) and 8.1 MB of value skeletons. The load at 100,000 steps
was a second run with a memory guard of 12 GiB so the process had room to hold the history twice
(the first run skipped it: 5.2 GB resident against half a 7 GiB guard); the head was shown after
4.0 s and the whole 100,813-step history loaded in 75.1 s (50,000 steps: 1.5 s and 18.0 s).

Share
forms: from 10,000 steps a state alone is 4.5 to 16x smaller than the whole history and steps only
28 to 61x (at 1,000 steps 1.7x and 16x); recomputing results makes the file 1.7 to 5x smaller,
because analysis results are the largest values.

Compression per chunk by measured gain: gzip is
kept only where it saves a sixteenth; skeletons inflate 6 to 15x (the lazy values chunk about 6.8x),
raw float arrays 1.28x (1.27x measured earlier), cell sets 1.37x.

Parameter arrays a step names are
lazy unless the head's values name them, so the head is shown after 1.4 s at 50,000 steps while the
rest reads behind it; the whole load stays where the first design's was.

The session file is
written at 100,000 steps because JSON is cut into pieces of 16 MiB, so no string nears the engine's
limit; at 10,000 to 50,000 steps the file is a third smaller than before (`$part`: gene names and
growing user-field lists once).

A named save was 2 to 2.7x slower than the first design's (30 s at 100,000 steps in the first run)
because it encoded every value, part and cell set of the whole history again before the store found
most of them already kept by the autosave (the named save adds 98 MB to the 1.5 GB the autosave
holds). Profiled on a history of 1,386 steps and 200,000 cells (named save 875 ms), measuring
cell-set encodings was half of it and an `await` per node of every value most of the rest. Each
array's cell-set bytes are now measured once (`history-file.js` `encodeArray`, cached per array as
its digest is) and `splitValue` waits once per array and per part, not per node: the same save takes
304 ms (2.9x quicker), and a session file written after the autosave reuses the same encodings. In
the second full run (beside other jobs) a named save takes 0.12 / 1.11 / 18.2 / 19.4 s at 1,000 /
10,000 / 50,000 / 100,000 steps: the snapshot and the stand-in's structured clone of every new
record are what is left at scale.

A collect (what removing a saved state anywhere, and the page's start, runs over every record this
browser keeps, `autosave-store.js` `collect()`, which plans its deletions in `planCollect`) takes 8 ms at 1,000 steps and 55 ms at 10,000 in the
suite's IndexedDB stand-in, linear in the records; from 50,000 steps the stand-in holds no records
(the suite lets them go to fit in memory), so it is not measured there. A browser's IndexedDB adds
the structured clone of every record it reads (the parts' 283 MB of text at 100,000 steps), so a
collect there costs seconds at that size, off the recording path; reference counts per record would
make it O(change) (not built).

### The history core against the first design

Rows the first comparison does not hold (second design, 2026-11-04).

| What | 1k | 10k | 100k |
|---|---|---|---|
| Delete a step with 10 / 100 after it (only dependants applied again), ms | 5 / 33 → 7 / 40 | 87 / 290 → 24 / 44 | 515 / 4,817 → 179 / 213 |
| Jump to a random step (the path between the two), p50 / max ms | 0.6 / 0.7 → 1.3 / 2.1 | 10.1 / 11.0 → 8.9 / 13.4 | 41.4 / 53.0 → 43.0 / 85.1 |
| Analysis panel switched, recorder p50 ms | | 11.97 → 9.73 | 13.36 → 9.62 |
| Differential expression / markers / correlation, recorder p50 ms | | 8.95 / 18.54 / 13.09 → 10.59 / 19.17 / 10.79 | 15.62 / 19.82 / 14.15 → 11.48 / 16.78 / 10.87 |
| Label edit / user field, recorder p50 ms | | 0.58 / 1.09 → 0.70 / 1.51 | 1.48 / 2.09 → 1.19 / 1.85 |
| `history.summary`, characters (no provider in the first design) | → 664 | → 1,941 | → 16,103 (102 branches) |

With 1,000 steps after the deleted one at 100,000 steps, 69 steps were applied again (491 ms against
42.9 s for all of them); what remains of a delete at 100,000 is linear and cheap: reading the values
at the deleted step's parent from the root (one walk) and the dependants' own dispatches. A jump
costs the steps between the two: a random step at 100,000 is about 57,000 steps away, walked once
each way, so a random jump is as long as before; a jump to a near step, an undo and a redo cost a
few steps. `history.summary` grows with branches, groups and turns, never with steps (16 KB at 102
branches) and takes 87 ms at 100,000 steps, most of it the stale check over the path to the head.
Hashing big values is cheaper but not gone: the canonical text is written in one synchronous pass
after the arrays' digests (no wait per node), so switching the analysis panel went from 12 to 13 ms
to 9.6 to 9.7 ms; what remains is SHA-256 over the canonical text of every result the analysis
domain holds (1.5 to 10 MB), which a new result needs anyway. The label edit's and user field's
recorder still grows with the user fields held (0.70 to 1.19 ms from 10,000 to 100,000 steps).
The autosave write per step after the history core alone (the walk of the pool gone, the copy of
every stored step id per write not yet) was 2.2 / 24.3 / 179.6 ms at 1,000 / 10,000 / 100,000 steps.

### After the fix pass

2026-11-05, seed 1, same machine; second design → after.

What the pass changed that the suite sees: every edit runs alone, the pool releases what nothing
refers to, a domain that follows another (the analyses follow the highlight pages) is captured with
it, observed changes are read as their gesture's task ends, a delete compares inputs in canonical
form, and steps taken out at once (a squash, a branch delete) have their depths computed again once.

| | 1k | 10k | 100k |
|---|---|---|---|
| Recorder p50 / p95 ms, every step | 0.22 / 0.80 → 0.29 / 0.92 | 0.20 / 1.48 → 0.31 / 1.43 | 0.15 / 1.24 → 0.24 / 2.33 |
| Brush of 1,000 cells, recorder p50 / p95 ms | 0.36 / 1.46 → 0.53 / 1.36 | → 1.35 / 2.99 | → 2.46 / 5.36 |
| Delete a step with 1,000 after it, ms | 595 → 581 | 41 → 25 | 491 → 414 |
| Delete with 10 / 100 after, ms | 7 / 40 → 4 / 34 | 24 / 44 → 24 / 26 | 179 / 213 → 140 / 156 |
| Undo p50, ms | 0.21 → 0.13 | 0.16 → 0.34 | 0.08 → 0.08 |
| Jump p50 / max ms | 1.3 / 2.1 → 0.9 / 1.4 | 8.9 / 13.4 → 10.8 / 21.5 | 43.0 / 85.1 → 58.3 / 92.1 |
| A 500-step squash near the root, 100,000 steps below (fake page) | | | 6,251 → 89 ms |

Capturing the analyses with every highlight change first cost a hash of every result they held: a
brush's recorder p50 went to 11.7 ms at 100,000 steps in an intermediate run. The analysis domain
now gives back the value it captured last while nothing about it changed, and the pool never
hashes a value object twice (`value-pool.js`), so the follow-up costs a comparison: 2.5 ms p50 at
100,000 steps with the mixed history's analyses held (about 1.6 MB of results), 0.5 ms at 1,000. A
squash or branch delete takes all its steps out first and computes the depths below once (6.3 s to
89 ms). A delete compares its inputs in canonical form: at 1,000 steps the 1,000-later-steps delete
is where it was (most of its dependants change what they read); the others are flat or faster.

### Smaller measurements

| What | Measured |
|---|---|
| A step of `history.get` | About 340 characters (293 before `kind`, `v`, `group` and `squash` were added), so a page of 200 steps is about 68,000 characters, well under a connector's tool result limit (`HISTORY_PAGE_STEPS`, `RESULT_TEXT_CHARS`) |
| The cell-set bytes of each `Uint32Array`, measured once and cached by the array's digest (`encodeArray`, `session/history-file.js`) | Measuring the cell-set encodings was half of a named save: on a 1,386-step history of 200k cells a named save took 875 ms, and with the bytes measured once per array and one await per array and part instead of per node it takes 304 ms (2.9x quicker); a session file written after the autosave reuses the same encodings (2026-11-04) |
| Category state as typed arrays (`utils/category-state.js`) | A field of 1M categories (clone or lineage ids, a text column of one label per cell) holds 4 MB shared by every view, where an RGB array per category plus an object keyed by index held about 90 MB per view and took about 50 ms to check on every action |
| The visibility kernels of `computeGlobalVisibility` as separate functions (`state/managers/filter-manager.js`) | Inside one large method a loop that had not run yet was optimized without type feedback and deoptimized on its first real use: a 40 ms first outlier step at 562,000 cells |
| The write of a hidden-cell list drives its loop bound by the buffer length | V8 keeps the elided bounds-check fast path; a narrower "touched range" bound measured 4x slower |
| One range scan per values array in `ensureContinuousMetadata`, kept with the array (`state/managers/color-manager.js`) | The method runs from coloring, filtering, the legend and the color domain several times per field switch; scanning every value on each run was most of a switch's main-thread time at 5,000,000 cells |
| A label element a category only while the view shows labels (`state/managers/color-manager.js`) | A field of a million categories made a million elements on every switch and filter, shown or not: 7-11 s of the page at 10M cells |
| View contexts share the per-cell arrays and the field inventories instead of copying them (`state/managers/view-context-core.js`, `view-context-viewer-sync.js`) | A copy of four per-cell arrays on each sync was hundreds of megabytes of transient memory while one large dataset loaded, and a copy of every field object (8,192 genes) cost milliseconds on every gene pick. The alpha-texture path of a transparency update is about 16x faster than a buffer rebuild |

### The autosave and the browser

**Cells in the autosave** (same run, bytes a step in IndexedDB): a brush of
about 1M cells 3,999,997 → 7,429; 500,000 scattered cells 2,003,696 →
129,551; 10,000 more cells hidden 1,556,493 → 1,406.

**A tab closed during a write** (Firefox, by 2026-11-06): with every write of
the transaction made at once and the transaction told to commit, 60 closes
mid-write kept 54 writes whole and 6 not at all, never a part. Before, 3
closes in 30 kept a part (98 and 252 of 400 writes).

**From a change to its write.** A write is handed to the browser some tens of
milliseconds after the change (a task, the session header, the digests of new
values): measured 20 to 120 ms on a loaded laptop.

## Tried and reverted

| What | When | Result | Would need |
|---|---|---|---|
| A branch as the name of its current last step | until 2026-11-04 | No delete of a branch, names not unique, fork points hard to find (Decisions, above). Replaced by branch records | — |
| Recursive walks of the tree (`history.get` with every branch, the steps below a deleted one) | until 2026-11-04 | `RangeError: Maximum call stack size exceeded` about 5,000 steps deep. Replaced by loops over steps kept in id order | — |
| A walk of the value pool after every step (the memory accounting, the autosave's size) | until 2026-11-04 | 0.8 s a step at 100k steps for each. Replaced by a byte counter kept as values come and go | — |
| Cells stored raw by the autosave | until 2026-11-04 | 4.0 MB for a brush of about 1M cells. Replaced by the cell-set encoding the session file uses | — |
| A stored picture of each step | removed 2026-11-04 | Never shipped in a release: judged to bloat the history | A reason to keep pictures that outweighs their bytes |
| Hashing a result once by its identity (a sub-value's hash in the canonical text) | 2026-11-04, not adopted | It would change the value hash rule that the session file, the autosave and the Python and R readers verify. Sub-values are deduplicated on the storage side instead | One change of the hash rule across the three repositories |
| A timeline that rebuilt its whole list on every render (`replaceChildren` with every row) | until 2026-11-04 | At 10,000 steps opening the panel held the page 660 to 750 ms and one more step took 300 to 370 ms; scrolling ran at about 5 fps over about 100,000 more DOM nodes (Measurements, "The timeline of the first design against the virtualised graph"). Replaced by a graph that keeps only the rows in view | — |
| A named save that measured every cell-set encoding again and awaited once per node of every value | until 2026-11-04 | 875 ms on a history of 1,386 steps and 200k cells, half of it measuring (Measurements, "Smaller measurements"). Replaced by bytes measured once per array and one await per array and part | — |
| A narrower "touched range" bound for the loop that writes a hidden-cell list | not adopted | 4x slower than the buffer length, which keeps V8's elided bounds-check fast path | — |
| Writes of one transaction made one at a time, each awaited | until 2026-11-06 | A closed tab kept a part of a write (Persistence, above). Replaced by all writes at once and an explicit commit | — |

## Known limits and open items

- **A closed tab can lose the change of its last moment.** A page that goes
  before the write is handed to the browser has not handed it over; a page
  that is left retires its autosave, which writes what is left if the
  browser gives it the time. Nothing can make it.
- **Hashing big values is on the input path.** What a new analysis result
  costs the recorder is SHA-256 over the canonical text of every result the
  analysis domain holds (9.6 to 9.7 ms for a panel switch at 10k and 100k
  steps in the run above). Hashing in a worker is not built: it would take the remaining 2.4 ms of a brush of 1,000,000 cells off the input path.
- **A jump costs the steps between the two.** A jump to a random step at 100k
  steps is about 57,000 steps away and takes about as long as in the first
  design (p50 41 → 43 ms); an undo, a redo and a jump to a near step cost a
  few steps.
- **The pool of values is not bounded.** It keeps every distinct value for
  the dataset's life (at 100,000 steps 0.87 GB of heap and 1.9 GB of arrays in
  Node); the suite measures it and no cap exists.
- **A collect reads every record.** What removing a saved state anywhere, and
  the page's start, run over the records this browser keeps is linear in them
  (8 ms at 1,000 steps and 55 ms at 10,000 in the suite's stand-in); in a
  browser the structured clone of what it reads makes it seconds at 100,000
  steps, off the recording path. Reference counts per record would make it
  O(change) and are not built.
- **R reads a long session file superlinearly.** After persistence 10 times the
  steps take 12 times the time (0.85 s, 8.5 s, 105 s); the cause (a list grown
  one element at a time or a lookup by scanning) is not profiled.
- **The page's choice of "all pages" is derived, not stored.** Storing the
  person's choice in the analyses instead of the derived list would make a
  page change change nothing there and take the follow-up cost of a highlight
  change lower still; deferred.
- **The diff dialog's step pickers list every step twice** and were not
  measured; they would follow the timeline's row costs.
- **Not measured, and why.** The steps-only and state-only share forms in the
  first design (neither existed); groups, squashes and paused periods in the
  first design; a real owner's capture and restore cost in Node (the page is
  synthetic); IndexedDB's own time and on-disk size in Node (the stand-in
  counts bytes and clones, the browser spec reads
  `navigator.storage.estimate()`); a saved state read back above 10,000 steps
  (the stand-in stops holding records past 10,000 so one process does not hold
  the history twice); Firefox and WebKit (the browser half runs in Chromium
  only); frame times in a real page at 50,000 and 100,000 steps (the first
  design's timeline already held the page for seconds a step at 10,000, so
  larger figures would have measured the list; the graph is measured at 10,000
  steps and 100 branches only); paint and composite time by name (Chromium's
  `Performance.getMetrics` gives task, script, layout and style, so paint is
  what remains, inferred, not traced); the analysis run with recording on
  against bypassed (only the first run computed).
- **A browser with no IndexedDB keeps no saved state**, and a saved record of
  another layout is listed, never opened.

## Where the code and tests live

Every guard of this layer, with what makes it fail, is listed once:
{ref}`drift-guards`.

| Concern | Module(s) | Tests |
|---|---|---|
| Catalog and definitions | `cellucid/assets/js/app/actions/catalog.js`, `cellucid/assets/js/app/actions/definitions/`, `cellucid/assets/js/app/actions/schema.js`, `cellucid/scripts/generate-action-versions.mjs`, `cellucid/tests/fixtures/action-versions/` | `cellucid/tests/action-catalog-contract.test.mjs`, `cellucid/tests/action-catalog-document-contract.test.mjs`, `cellucid/tests/action-versions-contract.test.mjs`, `cellucid/tests/action-reads-contract.test.mjs` |
| Dispatcher, stops and rollbacks | `cellucid/assets/js/app/actions/dispatcher.js`, `cellucid/assets/js/app/actions/control-dispatch.js` | `cellucid/tests/action-history-core.test.mjs`, `cellucid/tests/action-dispatch-stop-contract.test.mjs`, `cellucid/tests/action-stop-reaches-work-contract.test.mjs` |
| State domains | `cellucid/assets/js/app/actions/domains.js`, `cellucid/assets/js/app/state/domains/` | `cellucid/tests/state-domains-contract.test.mjs`, `cellucid/tests/browser/state-domains-roundtrip.spec.mjs` |
| Who may change state | `cellucid/assets/js/app/actions/classification.js` | `cellucid/tests/action-method-classification-contract.test.mjs`, `cellucid/tests/action-mutation-owners-contract.test.mjs`, the control walk |
| The history's core | `cellucid/assets/js/app/history/history.js`, `cellucid/assets/js/app/history/model.js`, `cellucid/assets/js/app/history/tree.js`, `cellucid/assets/js/app/history/value-pool.js`, `cellucid/assets/js/app/history/journal.js`, `cellucid/assets/js/app/history/provenance.js`, `cellucid/assets/js/app/history/replay.js` | `cellucid/tests/history-v2-core.test.mjs`, `cellucid/tests/history-model-contract.test.mjs`, `cellucid/tests/history-actions-contract.test.mjs`, `cellucid/tests/history-move-stop-contract.test.mjs` |
| Autosave and saved states | `cellucid/assets/js/app/history/autosave.js`, `cellucid/assets/js/app/history/autosave-store.js`, `cellucid/assets/js/utils/indexed-db-database.js` | `cellucid/tests/history-autosave.test.mjs`, `cellucid/tests/indexed-db-database.test.mjs`, `cellucid/tests/browser/history-autosave.spec.mjs`, `cellucid/tests/browser/saved-state-close.spec.mjs` |
| The history in words, the AnnData copy | `cellucid/assets/js/app/history/exports/`, `cellucid/assets/js/app/history/anndata-copy-save.js` | `cellucid/tests/history-exports.test.mjs`, `cellucid/tests/history-anndata-copy.test.mjs`, `cellucid/tests/browser/history-anndata-copy.spec.mjs` |
| Scale | `cellucid/scripts/history-stress.mjs`, `cellucid/tests/stress/` | `cellucid/tests/history-stress-harness.test.mjs`, `cellucid/tests/browser/history-stress.spec.mjs` (opt-in) |
| The shared documents (catalog, history model, golden history) | `cellucid/scripts/sync-shared-files.mjs`, `cellucid/tests/fixtures/` | `cellucid/tests/shared-files-contract.test.mjs`, `cellucid-python/tests/test_action_catalog.py`, `cellucid-r/tests/testthat/test-actions.R` |

## Related pages

- {doc}`29_action_layer_domains_and_history`: the reference.
- {doc}`30_action_reference`: every action, generated from the catalog.
- {doc}`10_sessions_persistence_and_serialization`: session format 4.
- {doc}`25_decisions_exports`: the writer behind "Save AnnData with this
  state".
- {doc}`34_decisions_agents`: turns, permissions and what an agent may do to
  the history.
- {doc}`38_decisions_testing_and_hosted_runners`: the control walk and the
  specs that close a page.
- {doc}`../m_state_history/index`: what a user sees.
