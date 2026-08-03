# Screenshot coverage audit — `l_sessions_sharing` and `n_benchmarking_performance`

Audited read-only against `cellucid/assets/js` (session subsystem, notification
center, `index.html` markup), `cellucid/tests/browser`, and
`cellucid-datasets/exports`. The app was **not** launched and no benchmark was
run (a sibling unit holds the rendering-measurement lock). Every DOM id, button
label, accordion title, and user-visible string quoted below was read out of the
source, not inferred.

---

## Summary

| Metric | `l_sessions_sharing` | `n_benchmarking_performance` | Total |
|---|---:|---:|---:|
| Pages (`.md` files) | 13 | 10 | 23 |
| `{figure}` directives | 4 | 5 | 9 |
| Unique image files referenced | 1 | 3 | 4 |
| Images the section *owns* (live under its own `_static/screenshots/<topic>/`) | **0** | 1 | 1 |
| Pages with at least one figure | 4 | 3 | 7 |
| Pages with **zero** figures | 9 | 7 | 16 |
| Borrowed-image directives | 4 (100 %) | 3 (60 %) | 7 (78 %) |

:::{important}
**The table above records the state this audit was taken against, and section L
has since been captured.** `_static/screenshots/sessions_sharing/` now holds the
section's own images, embedded on `02`, `03`, `04`, `10` and `11`. The counts
above and the per-page **has** / **verdict** entries below therefore predate
them — read {doc}`/user_guide/web_app/l_sessions_sharing/11_screenshots` for the current gallery, and
the directory listing for the current file set. Do not refresh the counts here:
they are derivable from the docs tree and will go stale again the moment a
capture lands.
:::

At audit time the `sessions_sharing/` directory was empty; all four figures in
section L were the same borrowed file,
`data_loading/data-loading-session-panel.png`; and section N owned exactly one
image, `benchmarking_performance/benchmark-panel.png` (1440 × 1000), used on two
pages, its other three directives borrowing `web_app/multiview-two-panels.png`
and `web_app/dark-theme-multiview.png`.

### Verdicts

| Verdict | Count | Pages |
|---|---:|---|
| REPLACE | 5 | L:``index``, L:`01`, L:`11`, N:`02`, N:`08` |
| ADD | 11 | L:`02`, L:`04`, L:`05`, L:`07`, L:`09`, L:`10`, N:``index``, N:`03`, N:`06`, N:`07`, N:`09` |
| SEQUENCE | 2 | L:`03`, N:`05` |
| NONE (justified) | 3 | L:`06`, L:`08`, N:`04` |
| DIAGRAM | 2 | L:`12`, N:`01` |
| OK (keep as-is) | 0 | — |
| *additionally flagged* STALE-RISK | 7 | L:`01`, L:`03`, L:`04`, L:`07`, L:`10`, N:`05`, N:`08` |

### Resolution note applying to every capture below

~~Existing screenshots are 1×; every new capture in this audit should be taken
at **`deviceScaleFactor: 2`**, and that change of convention must be applied in
one pass or the gallery will visibly mix sharp and soft images.~~ **Adopted.**
The rule now lives in ``00_capture_tooling_and_conventions`` — fixed viewport,
device scale factor 2, and an emitted `:width:` of the captured region's CSS
width × 2 capped at 1440 px, printed by the capture tool. Take the viewport and
the width rule from that page rather than from the per-step specifications
below, which were written before it existed.

---

# LIFECYCLE: session save and restore

This was the headline deliverable: at audit time every page in section L shared
one borrowed sidebar crop and nothing in the docs showed a user what a session
round-trip looks like. **Section L has since been captured** — see the note
under *Summary* — so read the sequence below as the specification it was written
as, and {doc}`/user_guide/web_app/l_sessions_sharing/11_screenshots` for what was actually shot. The
sequence is a single scripted Playwright run producing 14 captures. Steps 6 and
13 are the **proof pair**.

## Fixed run parameters (identical for every step)

| Parameter | Value | Why it is fixed |
|---|---|---|
| Viewport | 1440 × 900 CSS px, `deviceScaleFactor: 2` | Floating panels are clamped to the viewport on restore ({doc}`/user_guide/web_app/l_sessions_sharing/09_edge_cases`), so a BEFORE/AFTER pair captured at different sizes is not a valid comparison. |
| Dataset | `#dataset-select` → option value `dataset:local-demo:pancreas` | 3,696 cells / 3,753 genes — loads in seconds; publishes 1D + 2D + 3D embeddings with `default_dimension: 3`; recognisable endocrinogenesis trajectory for a wet-lab reader. Verified in `cellucid-datasets/exports/pancreas/dataset_identity.json`. |
| Dimension | 3D throughout | The session fingerprint records `cellOrder.dimension`; saving in 3D and restoring in 3D is the success path. Deviating provokes the recoverable error captured separately in `07`. |
| Theme / background | app defaults (light / grid) | Matches the reviewed official starting state, so the reset in step 10 is visually unambiguous. |
| Welcome modal | dismissed with `Escape` (`#welcome-modal`) — see `tests/browser/helpers/welcome.mjs` | Otherwise it covers the sidebar in every frame. |
| Build stamp | `#web-build-version` inside `#copyright` must read the same value in steps 6 and 13 | The docs tell readers to compare footer builds; the proof pair must be self-consistent. |

Output directory for every file: `_static/screenshots/sessions_sharing/`.

## The ordered sequence

### Step 1 — a clean app with the sample picker open

- **Interaction:** load the app, `Escape` to dismiss `#welcome-modal`. Do not
  select anything yet.
- **Crop:** `#session-section` (the accordion whose `<summary>` reads
  **Session**) — full block including `#dataset-info` and `#dataset-select`.
- **Cursor:** yes, hovering `#dataset-select`.
- **File:** `session-lifecycle-01-picker.png` · `:width: 246px`
- **Why:** this is the *only* accordion the whole section refers to and its
  title is **Session**, not "User data" (see the stale-prose finding on `01`,
  `03`, `05`). One correct capture fixes the navigation instruction for nine
  pages.

### Step 2 — the reviewed starting view (the state we will deliberately depart from)

- **Interaction:** `#dataset-select` → `dataset:local-demo:pancreas`; wait for
  the official static state to apply (`cell_type` active, 3D Orbit close-up).
- **Crop:** full page, 1440 × 900.
- **Cursor:** no.
- **File:** `session-lifecycle-02-sample-default-state.png` · `:width: 1440px`
- **Why:** establishes "this is what everyone gets by default", so the restored
  state in step 13 is provably *not* the default.

### Step 3 — change the colouring

- **Interaction:** `#categorical-field` → `clusters_coarse` (5 categories:
  Ductal, Ngn3 low EP, Ngn3 high EP, Pre-endocrine, Endocrine). Deliberately
  **not** `clusters`, whose 8 category labels are identical to `cell_type`'s and
  would make the before/after legend look unchanged.
- **Crop:** `#coloring-filtering-section` (summary **Coloring & Filtering**),
  including `#legend`.
- **Cursor:** yes, on `#categorical-field`.
- **File:** `session-lifecycle-03-colouring.png` · `:width: 246px`

### Step 4 — add a filter (and show the "apply once" pattern in passing)

- **Interaction:** `#continuous-field` → `G2M_score`; in `#legend`, toggle
  **Live filtering** off, drag the Min slider up, click **FILTER** once.
- **Crop:** `#coloring-filtering-section` from the continuous-field row down
  through `#active-filters-container` — must include `#filter-count` (its text
  changes from `Showing all points` to `Showing N of 3,696 points`) and the
  `#active-filters` entry.
- **Cursor:** yes, on the **FILTER** button.
- **File:** `session-lifecycle-04-filter.png` · `:width: 246px`
- **Bonus:** this single frame is also the missing "Live filtering off → FILTER
  once" illustration that section N repeats on five pages and never shows.
  Cross-reference it from `n_benchmarking_performance/03`.

### Step 5 — make a selection (this is what makes the round-trip non-trivial)

- **Interaction:** `#highlighted-cells-section` (summary **Highlighting**) →
  click the highlight-mode button `[data-mode="lasso"]` (**Lasso**) → drag a
  closed loop on `#glcanvas` around the Endocrine lobe → name the resulting
  group `Endocrine lobe` in `#highlighted-groups-list`.
- **Crop:** two captures.
  - 5a: `#glcanvas` with the lasso path still drawn — cursor **on the canvas at
    the closing point of the loop**; `session-lifecycle-05a-lasso.png` ·
    `:width: 1440px`.
  - 5b: `#highlighted-cells-section` showing `#highlight-pages-tabs`,
    `#highlight-count` (now `N cells highlighted`) and the named group row;
    `session-lifecycle-05b-highlight-group.png` · `:width: 246px`.
- **Why it matters:** a lasso group is stored as a `highlights/cells/<groupId>`
  **lazy** chunk. Without a selection the round-trip only exercises eager
  chunks, and the pages' central claim — that eager/lazy is one atomic
  operation, not partial success — is never demonstrated.
- **Privacy:** the group name must be a biology term, never a donor or patient
  label (`06` best practice #4, `08` rule 2).

### Step 6 — BEFORE (half of the proof pair)

- **Interaction:** orbit `#glcanvas` to a distinctive 3D angle (keyboard
  `W`/`A`/`S`/`D` per `#glcanvas-help`, so the pose is scriptable and
  repeatable); click **Keep view** (`#split-keep-view-btn`); switch
  `#categorical-field` to `cell_type`; click **Keep view** again; set
  `#view-layout-mode` to `Grid compare`.
- **Crop:** full page, 1440 × 900, sidebar visible.
- **Cursor:** no (a cursor in the proof pair invites "the pointer moved, so
  something else moved too").
- **File:** `session-lifecycle-06-before.png` · `:width: 1440px`

### Step 7 — open the Session panel and press Save State

- **Interaction:** scroll the sidebar to `#session-state-controls`.
- **Crop:** `#session-state-controls` only — the `Session state:` label, the
  `i` info button, and the two buttons `#save-state-btn` / `#load-state-btn`.
- **Cursor:** yes, **on `#save-state-btn`**.
- **File:** `session-lifecycle-07-save-state-button.png` · `:width: 246px`
- **Note:** this tight crop is what pages `01`, `03`, `11` and ``index`` actually
  need; the borrowed `data-loading-session-panel.png` shows the whole loading
  stack and buries the two buttons.

### Step 8 — the in-app confirmation

- **Interaction:** click `#save-state-btn`.
- **Crop:** `#notification-center`, the single success toast.
- **Cursor:** no.
- **File:** `session-lifecycle-08-saved-toast.png` · `:width: 400px`
- **Exact string** (from `ui/modules/session-controls.js`): `Session saved
  successfully`. On failure the toast reads `Failed to save session` while the
  console logs `Failed to save state:` — {doc}`/user_guide/web_app/l_sessions_sharing/10_troubleshooting_sessions`
  currently tells readers to look for the console string in the panel.

### Step 9 — the artifact the user actually receives

The Save State path is `downloadBlob()` → an anchor `click()` → a browser
download (`downloadBlob()` in `session-serializer.js`, called from
`downloadSession()`). There is **no** in-app artifact, no
link, and no server upload. The browser's download shelf and the OS file dialog
are browser chrome — Playwright cannot capture them, and both would expose
`/Users/<username>/`.

- **Capture instead:** a terminal screenshot, prepared as follows.
  - `cd ~/Downloads` and set `PS1='$ '` so no path or username is rendered.
  - Run: `ls -lh cellucid-session-*.cellucid-session` then
    `head -c 96 cellucid-session-*.cellucid-session | xxd | head -6`.
  - The hexdump shows the ASCII magic `CELLUCID_SESSION\n`, the 32-bit
    manifest length, and the opening of the JSON manifest — exactly the framing
    documented in {doc}`/user_guide/web_app/l_sessions_sharing/12_reference`.
- **Crop:** the terminal window, ~900 px wide, dark theme.
- **Cursor:** no.
- **File:** `session-lifecycle-09-artifact-terminal.png` · `:width: 900px`
- **Redaction checklist before publishing:** shortened prompt (no host, no
  user, no path), no other files in the `ls` output, timestamp in the filename
  is fine (it is generated, not personal).

### Step 10 — reset to a clean app

- **Interaction:** navigate to the app URL again (full reload, not a soft
  reset), `Escape` the welcome modal, re-select
  `dataset:local-demo:pancreas`, and wait for the official starting state.
- **Crop:** full page, 1440 × 900.
- **Cursor:** no.
- **File:** `session-lifecycle-10-reset.png` · `:width: 1440px`
- **Why re-select the same sample:** it re-applies the reviewed static state, so
  the app is provably back at the default of step 2 — single view, `cell_type`,
  no filters, no highlights — and the reader can see there is nothing left over
  to fake the restore.

### Step 11 — press Load State

- **Crop:** `#session-state-controls`.
- **Cursor:** yes, **on `#load-state-btn`**.
- **File:** `session-lifecycle-11-load-state-button.png` · `:width: 246px`
- **Do not capture the OS file picker.** It is out of Playwright's reach and it
  renders the user's home directory. Drive it with `setInputFiles` and state in
  the caption that the picker is the operating system's own dialog.

### Step 12 — the restore in flight

- **Interaction:** immediately after the file is supplied, capture the
  NotificationCenter download card while it is still progressing.
- **Crop:** `#notification-center` — the card titled `Loading session`, its
  `.notification-progress-bar`, `.notification-progress-text`,
  `.notification-speed`, and the `×` button
  (`.notification-dismiss[data-role="cancel"]`, `aria-label="Cancel"`).
- **Cursor:** yes, hovering the `×` — this is the **Cancel** the docs mention
  seven times and never show. It is a `×` glyph, not a button labelled
  "Cancel"; every page that writes "press **Cancel**" is describing this.
- **File:** `session-lifecycle-12-restore-progress.png` · `:width: 400px`
- **Timing:** pancreas restores fast. Throttle by capturing on the
  `startDownload` notification's first paint, or use a larger highlight group.

### Step 13 — AFTER (the other half of the proof pair)

- **Interaction:** wait for the terminal notification.
- **Crop:** full page, 1440 × 900 — identical framing to step 6.
- **Cursor:** no.
- **File:** `session-lifecycle-13-after.png` · `:width: 1440px`

**What must be identical between `06-before` and `13-after`:**

1. `#glcanvas` layout — three panels in `Grid compare`, in the same order.
2. The 3D camera pose in every panel (same azimuth, elevation, and zoom).
3. `#categorical-field`'s selected value and the `#legend` category list, in the
   same order with the same colours.
4. `#filter-count` text, verbatim, including the number.
5. Every row of `#active-filters`, verbatim.
6. `#highlight-count` text, the tab in `#highlight-pages-tabs`, and the
   `Endocrine lobe` row and its count in `#highlighted-groups-list`.
7. The `Dataset`, `Source`, `Cells`, `Genes` rows of `#dataset-info`.
8. `#web-build-version` in the footer.

**What is allowed to differ, and must be called out in the caption so a reader
does not mistake it for a failure:**

- The `#notification-center` stack: `13-after` carries the success toast,
  `06-before` does not.
- Nothing else. If any of the eight items above differ, the capture run is
  invalid and must be re-shot — publishing a mismatched "proof" pair would be
  worse than publishing no pair at all.

### Step 14 — the terminal success message

- **Crop:** `#notification-center` with both the completed download card and the
  info toast.
- **Cursor:** no.
- **File:** `session-lifecycle-14-fully-restored.png` · `:width: 400px`
- **Exact strings, from source:**
  - `session-serializer.js`, `notifications.info(...)` at the end of the restore
    — info toast: **`Session fully restored.`** (with a trailing period; the
    docs write it without one).
  - `session-serializer.js`, `notifications.completeDownload(...)` on the next
    statement — the download card completes with the same text.
  - `session-controls.js`, the `#load-state-btn` handler — a separate success
    toast: **`Session loaded successfully`**.
- **Correction this capture forces:** {doc}`/user_guide/web_app/l_sessions_sharing/03_save_restore_ux` (line 44) and
  {doc}`/user_guide/web_app/l_sessions_sharing/10_troubleshooting_sessions` (line 36) say the **panel status** reports
  `Session loaded successfully`. There is no status element inside
  `#session-state-controls` — the markup contains only a label, an info button,
  and the two buttons. Both strings are NotificationCenter toasts. The prose
  must be corrected in the same change as the capture.

### Optional step 15 — the same file, refused

Run once more without reloading: switch `#dimension-select` to 2D, then
`#load-state-btn` → the same file. This yields the recoverable rejection and is
the single most useful image for a wet-lab reader in the whole section. Details
under {doc}`/user_guide/web_app/l_sessions_sharing/07_versioning_compatibility_and_dataset_identity` below.

---

# Section L — `user_guide/web_app/l_sessions_sharing/`

## {doc}`/user_guide/web_app/l_sessions_sharing/index`

- **explains** — Hub page: what a session is, a four-row "choose your goal"
  table, and the reading order for the section.
- **has** — 1 figure, line 56, `data_loading/data-loading-session-panel.png`,
  `:width: 246px` (line 58). **Borrowed** from `b_data_loading`; the caption
  describes the loading paths, which this page does not teach.
- **verdict** — **REPLACE**
- **needs** — Drop the borrowed panel. Use the proof pair as the hub's hero: a
  side-by-side of `session-lifecycle-06-before.png` and
  `session-lifecycle-13-after.png` (or the two figures stacked), captioned "the
  same three-panel view, saved on Monday and restored on Friday". Fixture:
  pancreas, run as above. No cursor. If only one image is wanted, use
  `session-lifecycle-07-save-state-button.png` (`:width: 246px`) — the tight
  two-button crop is what the "Fast path" table actually points at.
- **notes** — The `Fast path` table row "Send someone an exact view" says to
  send the export folder + session; correct. No stale prose found.

## {doc}`/user_guide/web_app/l_sessions_sharing/01_session_mental_model`

- **explains** — What a session is, the three persistence layers, eager/lazy
  ordering, and the dataset-fingerprint rejection rules.
- **has** — 1 figure, line 215, `data_loading/data-loading-session-panel.png`,
  `:width: 246px` (line 217). **Borrowed**, identical to {doc}`/user_guide/web_app/l_sessions_sharing/index`'s.
- **verdict** — **REPLACE** (+ **DIAGRAM**, + **STALE-RISK**)
- **needs** —
  1. Replace the borrowed panel with `session-lifecycle-07-save-state-button.png`
     (`:width: 246px`, cursor on `#save-state-btn`).
  2. **DIAGRAM** for "the three persistence layers": in-memory (volatile) /
     `.cellucid-session` file / browser localStorage. This is a concept with no
     on-screen representation; a screenshot cannot show it. Mermaid block,
     rendered by the docs build — no capture needed.
  3. **ADD** a fingerprint capture for the "Dataset identity" subsection: crop
     `#dataset-info` (the eight-row block: Dataset, Description, Source, URL,
     Cells, Genes, Obs fields, Connectivity) on pancreas, which reads
     `3,696` cells / `3,753` genes. Cursor on the `Cells` row.
     `session-identity-dataset-info.png` · `:width: 246px`. This makes the
     abstract "five-part fingerprint" concrete: four of the five parts are
     literally on screen.
- **notes** —
  - ~~"Open **User data** → **Session state**". The accordion's `<summary>` is
    **Session** (the `<details id="session-section">` block in `index.html`);
    "User data" is a label on a *different* block inside it (`#user-data-block`,
    label `Local data:`).~~ **Done** — the step now reads "Open **Session** →
    **Session state**", and the string `User data` no longer appears anywhere in
    section L.
  - The third mismatch row asserts "The dataset has the same name and
    counts but stores its cells in a different order". The shipped message
    (the final branch of `describeDatasetFingerprintMismatch()` in
    `session-context.js`) deliberately does **not** assert re-ordering:
    it says the coordinates differ and offers *both* causes ("Either the cells
    are stored in a different order, or the same cells were exported again from
    a re-computed embedding"). The source carries a long comment explaining that
    asserting the alarming cause would be its own integrity failure. The docs
    contradict that intent. Same defect on `07` (line 64) and `10` (line 111).
  - Line 111 and elsewhere: the success string is `Session fully restored.`
    with a trailing period.

## {doc}`/user_guide/web_app/l_sessions_sharing/02_what_gets_saved_and_restored`

- **explains** — The explicit saved-vs-not-saved contract, chunk by chunk, plus
  a five-step "how to verify a restore" checklist.
- **has** — `none`
- **verdict** — **ADD**
- **needs** — Two images, both derived from the lifecycle run so they cost
  nothing extra.
  1. The verification checklist (lines 239–261) is a list of five places to
     look. Capture one **annotated composite** at `:width: 1440px` built from
     `session-lifecycle-13-after.png` with five numbered callouts pinned to:
     `#glcanvas` (camera/layout), `#categorical-field` (active field),
     `#active-filters` (filters), `#split-view-badges-list` (snapshot layout),
     `#highlighted-groups-list` (highlight counts). File
     `session-verify-restore-annotated.png`. Callouts drawn in post, not by the
     app.
  2. For the "Not saved" list, capture the *excluded* modules side by side:
     crop the three accordions carrying
     `data-state-serializer-skip="true"` — `#community-annotation-section`,
     `#figure-export-section`, `#benchmark-section` — as one stacked 246 px
     column, `session-excluded-modules.png` · `:width: 246px`. The exclusion is
     a real, inspectable attribute in the markup, so the image is verifiable
     rather than decorative.
- **notes** — The saved/not-saved lists match `session-serializer.js` and
  `state-serializer/README.md`. The chunk inventory on this page agrees with the
  table in {doc}`/user_guide/web_app/l_sessions_sharing/12_reference`. No stale claims found.

## {doc}`/user_guide/web_app/l_sessions_sharing/03_save_restore_ux`

- **explains** — The click-by-click Save State / Load State procedure, what
  happens under the hood, and naming/milestone best practice.
- **has** — 1 figure, line 162, `data_loading/data-loading-session-panel.png`,
  `:width: 246px` (line 164). **Borrowed**, third use of the same file.
- **verdict** — **SEQUENCE** (+ **STALE-RISK**)
- **needs** — This page is the natural home for the full lifecycle. Embed steps
  1, 7, 8, 9, 11, 12, 13, 14 in order, each under its own sub-heading matching
  the existing "Save State" / "Load State" structure. The proof pair (6 and 13)
  belongs here too, immediately under "How to confirm a restore worked".
  All parameters as specified in the LIFECYCLE section.
- **notes** —
  - ~~"Open **User data** (the accordion that also contains dataset loading
    controls)" — the accordion is **Session**.~~ **Done** — the step names
    **Session**.
  - Still open: "…and for the panel status to say **Session loaded
    successfully**." There is no panel status; both strings are toasts in
    `#notification-center`.
  - The page correctly warns that the catalog's `default.cellucid-session` is
    rejected by the manual loader. Re-verified against the shipped presets: all
    five carry exactly five chunks (`core/field-overlays`, `core/state`,
    `ui/dockable-layout`, `analysis/windows`, `highlights/meta`) in that order
    and no `cinematic/camera` or `analysis/cache-inventory`, which the manual
    profile requires. A Save State writes seven.

## {doc}`/user_guide/web_app/l_sessions_sharing/04_official_sample_states`

- **explains** — How a built-in sample's SHA-256-pinned starting state is
  fetched, verified, and applied automatically, and why it is not a Load State
  artifact.
- **has** — `none`
- **verdict** — **ADD** (+ **STALE-RISK**)
- **Do the presets have visible identities worth showing?** Yes, but not
  individually. All five states are structurally identical — each activates
  `cell_type`, sets a static 3D Orbit close-up, light theme, grid background, no
  filters, no snapshots (confirmed: all five manifests carry the same five
  chunks and `cellOrder.dimension: 3`). Photographing them one at a time would
  produce five near-identical sidebar crops. What *is* worth showing is that one
  click yields a **reviewed biological view rather than a grey blob**, and that
  the five look different from each other because the biology differs.
- **needs** —
  1. **A five-up contact sheet.** One scripted pass: for each of `suo`,
     `garcia`, `he`, `kanemaru`, `pancreas`, select
     `dataset:local-demo:<id>` in `#dataset-select`, wait for the state to
     apply, crop **`#glcanvas` only** (no sidebar), 600 × 600 each. Composite
     into a 5-panel row with the sample id and cell count under each
     (`suo` 561,947 · `garcia` 219,731 · `he` 71,650 · `kanemaru` 131,636 ·
     `pancreas` 3,696 — read from the shipped manifests). No cursor.
     `official-sample-states-contact-sheet.png` · `:width: 1440px`.
  2. **The picker itself**, cursor on `#dataset-select` with the option list
     open showing all five entries: `official-sample-picker-open.png` ·
     `:width: 246px`. Reuse `session-lifecycle-01-picker.png` if the option list
     can be forced open there instead.
  3. **Do not** capture the rejection of a downloaded `default.cellucid-session`
     here — that image belongs on {doc}`/user_guide/web_app/l_sessions_sharing/09_edge_cases`, which is where the
     behaviour is catalogued.
- **notes** — ~~**Stale value.** The page quotes a Suo `state_sha256` literal
  that no longer matches the catalog.~~ **Done** — the page now shows the
  field's *shape*, `"state_sha256": "<64 lowercase hexadecimal characters>"`,
  and tells the reader to read the current value out of
  `cellucid-datasets/exports/datasets.json`. **This record quotes no digest
  either, deliberately:** while the finding stood open the catalog was
  regenerated again, so the "current" value written down here went stale on its
  own. A published state is rebuilt whenever the reviewed view or the session
  format changes; the catalog is the only place the digest can be read.
  Everything else on the page checks out:
  the five-chunk profile, the `{ "states": [...] }` manifest root, the
  `local-demo`-only opt-in, and the no-probe source list all match the source.

## {doc}`/user_guide/web_app/l_sessions_sharing/05_share_workflows_links_bundles_exports`

- **explains** — The three shareable artifacts (link to the dataset, session
  bundle, export folder), the "portable pair" rule, and three sender/recipient
  workflows.
- **has** — `none`
- **verdict** — **ADD**
- **Does the URL bar need to be visible, and is that safe?** **No, and it is not
  needed.** Cellucid's "share a link" means *a link to the dataset* — a public
  URL or a GitHub `owner/repo/path` typed into `#github-repo-url` or the remote
  server field — **not** a shareable app URL that encodes state. There is no
  session-in-URL feature: `restoreFromUrl()` exists in the serializer but is not
  wired to any user-facing share control, and Save State is a plain file
  download. Photographing the browser's URL bar would therefore illustrate
  nothing on this page while risking a `localhost:PORT` or file path. Show the
  **in-app input fields** instead, which is where the link actually lives.
- **needs** — Three images, each answering one of the page's three artifact
  types:
  1. **Link:** crop the GitHub row — `#github-repo-url` pre-filled with a public
     path (use `theislab/cellucid-datasets/exports`) plus `#github-connect-btn`.
     Cursor on `#github-connect-btn`.
     `share-link-github-input.png` · `:width: 246px`.
  2. **Bundle:** reuse `session-lifecycle-09-artifact-terminal.png` — the `ls`
     of the `.cellucid-session` file.
  3. **The portable pair:** one terminal screenshot showing an export folder and
     the session file side by side, e.g. `tree -L 2 project/` rendering the
     exact layout the page recommends. Prompt reduced to `$ `, no home path.
     `share-portable-pair-terminal.png` · `:width: 900px`. This is the page's
     central rule and it is currently pure prose.
- **notes** — ~~"**User data → Session state → Save State**" — again, the
  accordion is **Session**.~~ **Done.** The fingerprint "cannot detect every
  content change" claim is accurate and correctly hedged.

## {doc}`/user_guide/web_app/l_sessions_sharing/06_collaboration_best_practices`

- **explains** — Eight team conventions: the context triple, directory layout,
  milestones, safe highlight naming, review packets, and clean-tab validation.
- **has** — `none`
- **verdict** — **NONE (justified)**
- **needs** — Nothing. Every item on this page is a *convention*, not an
  interface: a folder layout (already an ASCII tree), a filename template, a
  naming policy, and a validation procedure that is deliberately performed in a
  *clean browser profile with extensions disabled* — a state that cannot be
  distinguished in a screenshot from the normal app. The one thing worth
  showing, "what a review packet looks like", is a `README.md` and a question,
  both already rendered as text. Adding a screenshot here would be decoration.
  If the section wants one visual anchor, cross-link to
  `session-lifecycle-06-before.png` under best practice #5 ("send one screenshot
  that shows the expected view") rather than capturing anything new.
- **notes** — Best practice #7's clean-tab check is the correct procedure and
  matches the strict reader's behaviour. No stale prose found.

## {doc}`/user_guide/web_app/l_sessions_sharing/07_versioning_compatibility_and_dataset_identity`

- **explains** — The five-field dataset fingerprint, the outcome matrix, why the
  reader rejects rather than partially restores, and cross-machine recipes.
- **has** — `none`
- **verdict** — **ADD** (+ **STALE-RISK**)
- **Which error states are worth capturing, and how to provoke each.** All four
  fingerprint rejections are produced by
  `describeDatasetFingerprintMismatch()` and surfaced verbatim as a red toast
  in `#notification-center` (the `#load-state-btn` catch in
  `session-controls.js` passes `err.message` straight through to
  `showSessionStatus`). For a wet-lab reader these four toasts are the most
  valuable images in the section — they turn "it didn't work" into "the app told
  me exactly what to do".

  | # | Error | How to provoke it (scriptable) | Worth capturing? |
  |---|---|---|---|
  | 1 | **Different dataset** — names both cell and gene counts | Save a session on `pancreas` (3,696 / 3,753). Reload, load `dataset:local-demo:he` (71,650 / 5,152). Load State the pancreas file. | **Yes — highest value.** The commonest real-world mistake, and the message names four concrete numbers. |
  | 2 | **Different dimension** — the recoverable one | Save on `pancreas` in 3D. Without reloading, set `#dimension-select` to 2D. Load State the same file. | **Yes — second highest.** This is the only failure that is *not* a data problem, and the page says so. Showing it prevents a user from concluding their dataset is corrupt. |
  | 3 | **Coordinates differ** (re-ordered rows *or* re-computed embedding) | Needs a second export with the same id and sizes but permuted rows. Clone `cellucid/tests/browser/fixtures/exports/current-ui-prepared` (120 cells, 6 genes, 2D only, 52 KB) into a scratch fixture, apply the same row permutation to `points_2d.bin` and every `obs/` column, keep `dataset_identity.json` byte-identical. Save a session against the original, load the permuted clone, Load State. | **Yes**, but as the *third* image — it is rarer and its message is long. Note the fixture is 2D-only, so error #2 cannot be provoked on it. |
  | 4 | **Pre-`cellOrder` file** (`SESSION_WITHOUT_CELL_IDENTITY_MESSAGE` in `session-context.js`) | Cannot be produced by the current app. Hand-write a bundle whose `datasetFingerprint` has exactly the four keys `sourceType`, `datasetId`, `cellCount`, `varCount`. Framing is documented in {doc}`/user_guide/web_app/l_sessions_sharing/12_reference` and implemented in `session/bundle/format.js`. | **Optional.** Capture only if a fixture author is willing to maintain a hand-built bundle. The message is self-explanatory in prose. |

- **needs** —
  - Crop for all four: `#notification-center`, the single red toast, full text
    visible and not truncated — the messages are long, so allow the crop to run
    to ~520 px wide, `:width: 520px`.
  - Cursor: **no** on the toast itself. For image #2 add a second frame with the
    cursor on `#dimension-select` showing the value the message asks you to
    switch back to — that pairs the diagnosis with the fix in one glance.
  - Files: `session-error-different-dataset.png`,
    `session-error-different-dimension.png`,
    `session-error-cell-order.png`, `session-error-no-cell-order.png`,
    `session-fix-dimension-select.png`.
  - Also **ADD** the identity crop from `01` (`session-identity-dataset-info.png`)
    here under "Exact fingerprint".
- **notes** —
  - Line 64 of the outcome matrix asserts the cause is re-ordered storage. The
    shipped message names both possible causes and refuses to choose. Fix the
    row to match.
  - The "six compared values" arithmetic (five fields, `cellOrder` contributing
    two) is correct and matches `datasetFingerprintMatches()`.
  - The claim that the Python reader accepts a fingerprint with or without
    `cellOrder` was not verified in this pass — it lives in `cellucid-python`
    and is outside these two sections.

## {doc}`/user_guide/web_app/l_sessions_sharing/08_security_privacy_and_trust`

- **explains** — What a session bundle can and cannot contain, the untrusted-input
  trust model, and four safe-sharing rules.
- **has** — `none`
- **verdict** — **NONE (justified)** — but this page owns the redaction policy
  for every other image in the section.
- **needs** — No screenshot. The page's subject is *what is inside a file* and
  *what not to put in a label* — neither has a visual form, and a screenshot of
  a session bundle's contents would either be a hexdump (already assigned to
  {doc}`/user_guide/web_app/l_sessions_sharing/12_reference`) or a fabricated example of the bad practice the page warns
  against. Adding one would be actively counterproductive.
- **Leak review of every image proposed in this audit** (this is the deliverable
  for this page):
  - **`#dataset-url`** (inside `#dataset-info`) renders `metadata.source.url`
    verbatim (the `datasetUrlEl` assignment in `dataset-controls.js`). For a
    remote or GitHub load this is
    a live URL and may carry a host, port, private org, or query string. Every
    capture that includes `#dataset-info` — `session-identity-dataset-info.png`,
    `session-lifecycle-01-picker.png` — **must** be taken on a `local-demo`
    sample, where the row renders a public catalog value, and the row must be
    checked frame by frame before publishing.
  - **`#github-repo-url`** in `share-link-github-input.png` must contain only the
    public `theislab/cellucid-datasets` path. Never a private org.
  - **The remote-server URL field** must be empty or public in any frame that
    includes it.
  - **Terminal captures** (`session-lifecycle-09`, `share-portable-pair-terminal`)
    must run with a reduced prompt (`PS1='$ '`) from inside the directory, so no
    `/Users/<username>/` string is rendered. Do not use `pwd`, do not show the
    window title bar (macOS Terminal puts the path there), and do not use a
    shell with a path-bearing prompt theme.
  - **The OS file picker and the browser download shelf must never be captured.**
    Both render the home directory. This is why step 11 captures the button and
    step 9 captures a prepared terminal instead.
  - **Highlight and page names** in any frame must be biology terms
    (`Endocrine lobe`), never donor, patient, barcode, or project codes — the
    rule this page states, applied to our own screenshots.
  - **Tokens:** none of the proposed images can contain one. GitHub OAuth tokens
    for Community Annotation live in browser session state and never appear in
    the Session accordion; the Community Annotation panel is not captured here.
- **notes** — The trust-model bullet list matches the reader implementation
  (exact manifest keys, chunk profile validation, gzip preflight with
  `maxOutputBytes`, transactional rollback). No stale claims found.

## {doc}`/user_guide/web_app/l_sessions_sharing/09_edge_cases`

- **explains** — Expected-but-surprising outcomes: differing identity on similar
  files, eager/lazy scheduling, empty-state replacement, cancellation, viewport
  clamping, large bundles, and why the official default is not a manual session.
- **has** — `none`
- **verdict** — **ADD**
- **needs** — Three captures, each pinned to a section that currently has no
  visual:
  1. **Cancellation** (section "Eager and lazy scheduling"): reuse
     `session-lifecycle-12-restore-progress.png` with the cursor on the `×`, and
     add one frame taken immediately after clicking it showing
     `#notification-center` **empty** — the docs' claim that cancellation
     produces "neither success nor a red product failure" is only convincing
     when you can see that nothing was published.
     `session-cancel-dismissed.png` · `:width: 400px`.
  2. **Official default rejected** (final section): download
     `exports/pancreas/default.cellucid-session`, choose Load State, supply it.
     Crop the red toast in `#notification-center`. Cursor: no. The property that
     makes it get refused is the five-chunk profile the manifest names
     `published-default` (`PUBLISHED_DEFAULT_STATIC_CHUNK_PROFILE` in the
     session serializer), not its size — **do not pin a byte count here.** The
     preset is regenerated
     whenever the session inventory changes, so the only current value is the
     file itself in `cellucid-datasets/exports/pancreas/`.
     `session-error-official-default-rejected.png` · `:width: 520px`. This is a
     mistake a curious user will genuinely make — the file is publicly
     downloadable from the datasets repo.
  3. **Viewport clamping** (section "Different screens and browsers"): save a
     session with a floated panel positioned near the right edge at
     1440 × 900, then restore it at 1024 × 700 and capture the clamped result.
     Two frames side by side, `:width: 1024px`,
     `session-panel-clamping-pair.png`. Caption must say this is expected, not a
     restore failure.
- **notes** — Section "Same lightweight fingerprint, changed content" (lines
  24–31) says the fingerprint "checks source type, dataset id, cell count, and
  variable count" and that "Reordering cells … can therefore pass the identity
  guard". **This is now wrong** — `cellOrder` was added precisely to catch
  re-ordering, and pages `01`, `07`, `10`, and `12` all describe it. This
  paragraph is a survivor of the pre-`cellOrder` text and directly contradicts
  the rest of the section. Highest-priority prose fix in section L.

## {doc}`/user_guide/web_app/l_sessions_sharing/10_troubleshooting_sessions`

- **explains** — Symptom → diagnosis → fix for save failures, picker problems,
  identity rejection, format rejection, cancellation, "it looks different", and
  slow restores; plus a bug-report checklist.
- **has** — `none`
- **verdict** — **ADD** (+ **STALE-RISK**)
- **needs** — This page should be an image *index*, not a new capture session.
  Reuse, with tight per-symptom placement:
  - Fast-diagnosis table row "Dataset identity error" → link the four
    `session-error-*.png` files from `07`.
  - "Know the success boundary" → `session-lifecycle-14-fully-restored.png`.
  - "Cancel and replacement" → `session-lifecycle-12-restore-progress.png` +
    `session-cancel-dismissed.png`.
  - "Did you download an official `default.cellucid-session`?" →
    `session-error-official-default-rejected.png`.
  - "It restored, but looks different" → the proof pair, plus the footer crop:
    `#copyright` including `#web-build-version`, cursor on the build string,
    `session-footer-build.png` · `:width: 246px`. Step 1 of that checklist is
    "app build shown in the footer" and no reader currently knows where that is.
  - **Not capturable, and the page should say so:** "File picker never opens"
    and "You canceled the picker" are OS-dialog states. Keep them as prose.
- **notes** —
  - "The Session panel then reports **Session loaded successfully**" — it is a
    toast, not a panel status.
  - "Check the panel status and browser Console for `Failed to save state`" —
    the console logs `Failed to save state:`; the visible toast reads
    `Failed to save session`. Give the reader both strings and say which appears
    where.
  - The "cells are stored in a different order" assertion the shipped message
    avoids is repeated here. Same fix as `01` and `07`.

## {doc}`/user_guide/web_app/l_sessions_sharing/11_screenshots`

- **explains** — Nothing; it is a 10-line stub whose entire body is one borrowed
  figure. Titled "Verified session-state capture" (singular).
- **has** — 1 figure, line 3, `data_loading/data-loading-session-panel.png`,
  `:width: 246px` (line 5). **Borrowed**, fourth and final use.
- **verdict** — **REPLACE**
- **needs** — This page must become the section's canonical gallery, exactly as
  {doc}`/user_guide/web_app/n_benchmarking_performance/08_screenshots` is for section N. Structure it
  as three headed groups and embed every file produced by the LIFECYCLE run:
  - *Save and restore round trip* — steps 1, 7, 8, 9, 11, 12, 13, 14 and the
    proof pair.
  - *Rejections* — the four `session-error-*.png` plus the official-default
    rejection.
  - *Identity and build* — `session-identity-dataset-info.png`,
    `session-footer-build.png`.
  Each caption must record the build stamp read from `#web-build-version`, the
  dataset id, and the viewport, in the style already used by section N's
  screenshots page. Title should become plural.
- **notes** — Nothing on this page is currently verified or session-specific.

## {doc}`/user_guide/web_app/l_sessions_sharing/12_reference`

- **explains** — The exact binary framing, the manifest schema with its nine
  required chunk keys, the singleton and dynamic chunk inventories, the atomic
  restore phases, and how the official profile differs.
- **has** — `none`
- **verdict** — **DIAGRAM**
- **needs** —
  1. **A framing diagram**, not a screenshot: magic bytes → `manifestByteLength`
     (u32 LE) → manifest JSON → repeated (`storedBytes` u32 LE + chunk payload).
     A byte-layout figure communicates this in one glance and cannot go stale
     against a UI change. Mermaid or an inline SVG committed beside the docs.
  2. **Optionally** the hexdump half of `session-lifecycle-09-artifact-terminal.png`,
     cropped to just the `xxd` output, so a reader can match the diagram to real
     bytes: `session-bundle-hexdump.png` · `:width: 900px`.
  - No app capture is appropriate; this page has no UI.
- **notes** — Verified against the shipped artifacts: all five official
  `default.cellucid-session` files begin with `CELLUCID_SESSION\n`, carry a
  five-key `datasetFingerprint` including `cellOrder` with `dimension: 3` and a
  16-hex digest, and list exactly the five chunks the page names, in the stated
  order. The `state_sha256` values in the catalog are present for all five
  samples. The implementation pointers at the end of the page all resolve to
  existing files.

---

# Section N — `user_guide/web_app/n_benchmarking_performance/`

**Constraint honoured:** no app run, no benchmark executed. Every capture below
is *specified*, none is *taken*, and the whole section's capture work must be
queued behind the sibling unit currently holding rendering measurements.

**Honest triage — which of these topics are pictures and which are numbers.**
Of the ten pages, four are fundamentally tabular or procedural and will not be
improved by screenshots: `01` (a conceptual model), `04` (a methodology and a
worksheet), `06` (a symptom catalogue, mostly), and `09` (a report template).
The in-app benchmark panel is the only genuinely photogenic surface in the
section, and it already has one image. The realistic upside here is **four or
five new captures**, not thirty — and one of them (the Live filtering / FILTER
control) is worth more than the rest combined because five separate pages give
that instruction and none of them shows it.

## {doc}`/user_guide/web_app/n_benchmarking_performance/index`

- **explains** — Hub: the three-bottleneck framing, a "if it feels slow right
  now" table, reading order, and a card grid.
- **has** — `none`
- **verdict** — **ADD**
- **needs** — One hero, and it should be the *diagnosis* rather than the raw
  panel: crop `#bottleneck-results` after **Analyze Performance** has run —
  `#bn-verdict-box` (title + detail), the `Current FPS` row, `Problems Found:`
  (`#bn-problem-list`), and `What to do:` (`#bn-fix-list`), with the collapsed
  `Show detailed stats` summary visible at the bottom. Cursor: no. Fixture: the
  synthetic benchmark at the `1M` preset with the default `Model Surface`
  pattern, so the verdict is reproducible without a real dataset.
  `benchmark-verdict-panel.png` · `:width: 246px`. This single image says "the
  app will tell you which of the three bottlenecks you hit", which is the
  page's whole thesis.
- **notes** — The fast-path table's advice matches the controls that exist. No
  stale claims found.

## {doc}`/user_guide/web_app/n_benchmarking_performance/01_performance_mental_model`

- **explains** — GPU vs CPU vs I/O, how each feels, the multipliers, a four-step
  triage workflow, and normal-vs-suspicious.
- **has** — `none`
- **verdict** — **DIAGRAM**
- **needs** — The subject is a mental model; a screenshot cannot depict "which
  resource is saturated". Produce one diagram: three bottleneck lanes (GPU /
  CPU / I/O) with the symptom that identifies each and the one knob that
  relieves it — essentially the "multipliers" table rendered as a figure.
  Mermaid, no capture. If one photograph is wanted, reuse
  `benchmark-verdict-panel.png` from the index under "Step 1 — identify the
  likely bottleneck", since the Analyze Performance verdict is the app's own
  implementation of that step.
- **notes** — "Turn Live filtering off and use **FILTER**" (lines 113, 116, 237)
  matches the control: `legend-renderer.js` renders a toggle labelled `Live
  filtering`, enabled by default (`liveFilteringEnabled = true`), beside a
  `FILTER` button, with the hint `Adjust limits; click FILTER or enable Live
  filtering.` Accurate — but see `03` for the missing image.

## {doc}`/user_guide/web_app/n_benchmarking_performance/02_performance_considerations_what_gets_slow_and_why`

- **explains** — The per-interaction cost model (loading, points, smoke,
  filtering, highlighting, analysis, export) and a symptom → knob cheat sheet.
- **has** — 1 figure, line 241, `benchmarking_performance/benchmark-panel.png`,
  `:width: 1440px` (line 243), under an "Interface reference" heading.
- **verdict** — **REPLACE**
- **needs** — The current figure is a 10-million-point synthetic GLB render. The
  page is about the cost of *real* interactions — filtering, highlighting,
  smoke, export — and never mentions the synthetic benchmark. The image is
  off-topic filler and its caption duplicates `08_screenshots.md` verbatim.
  Replace with the two controls the page actually tells you to reach for:
  1. The Live filtering / FILTER pair — reuse
     `sessions_sharing/session-lifecycle-04-filter.png`, or capture a dedicated
     crop of `#legend` for a continuous field showing the `Live filtering`
     toggle **off**, the Min/Max sliders, and the `FILTER` button, cursor on
     `FILTER`. `perf-live-filtering-off.png` · `:width: 246px`.
  2. The view-count control — crop `#split-view-controls` showing
     `#split-keep-view-btn`, `#view-layout-mode`, `#split-view-badges-list`
     populated with two badges, and `#split-clear-btn`; cursor on
     `#split-clear-btn` (the "clear snapshots" action the page recommends
     seven times). `perf-view-controls.png` · `:width: 246px`.
  Keep the benchmark panel image on `08_screenshots.md`, where its acceptance
  caption belongs.
- **notes** — The smoke-mode grid choices (32³, 48³, 64³, 96³, 128³, starting at
  128³) were not re-verified in this pass — they belong to
  `c_core_interactions` and the render-mode source, outside this section's
  scope. Flagging as unverified, not as wrong.

## {doc}`/user_guide/web_app/n_benchmarking_performance/03_large_dataset_best_practices`

- **explains** — An eight-item large-dataset survival checklist, a wet-lab
  workflow, a computational export workflow, three patterns, and a cliff list.
- **has** — 1 figure, line 252, `web_app/multiview-two-panels.png`,
  `:width: 1440px` (line 254). **Borrowed** from `web_app`, but genuinely
  on-topic — it illustrates the view multiplier the page warns about.
- **verdict** — **ADD** (keep the existing figure)
- **needs** — The checklist's most-repeated instruction, item 4 ("turn **Live
  filtering** off, move the sliders, click **FILTER** once"), has no image
  anywhere in the docs despite appearing on `01`, `02`, `03`, `06`, `07`, and
  `09` of this section. Add `perf-live-filtering-off.png` (spec under `02`)
  directly beside item 4. This is the single highest-value new capture in
  section N: it is the one instruction a wet-lab reader is told to follow and
  cannot currently find.
  - Optionally add item 3's control: `perf-view-controls.png` beside "Keep views
    (snapshots) lean".
- **notes** — The "< 100k / 100k–500k / 500k–2M+" bands are explicitly labelled
  "very rough intuition"; that hedging is appropriate and should not be
  hardened into a screenshot caption.

## {doc}`/user_guide/web_app/n_benchmarking_performance/04_benchmarking_methodology_and_metrics`

- **explains** — What to measure (TTFR, TTI, FPS, latency, memory, network), how
  to control confounders, cold vs warm, trial counts, browser tooling, and a
  copy/paste worksheet.
- **has** — `none`
- **verdict** — **NONE (justified)**
- **needs** — Nothing, and this should be a deliberate decision rather than an
  omission. Everything measurable on this page is **numbers in a table**: the
  worksheet is already a table, and the medians and ranges it asks for are text.
  The only screenshot candidates are *browser* UI — Chrome's Task Manager
  (`Shift+Esc`), the DevTools Performance Monitor, the Network waterfall,
  Firefox's `about:performance`. All four are browser chrome, outside the page
  and outside the Playwright capture tool's reach, and all four are re-skinned
  by browser vendors on a release cadence far faster than this documentation is
  revised. Screenshotting them would create four images guaranteed to go stale
  and impossible to re-shoot with the docs' own tooling. Keep the textual
  instructions.
- **notes** — Line 183 says "Chrome and Edge expose a tab-level task manager
  (`Shift+Esc`)". Correct at time of writing but exactly the kind of claim that
  drifts; it is already hedged by naming the alternative paths.

## {doc}`/user_guide/web_app/n_benchmarking_performance/05_benchmark_tools`

- **explains** — Where the in-app synthetic benchmark lives, its quick start,
  the reproducible workflow, the data patterns, how to read the stats, Copy
  Situation Report, and Analyze Performance.
- **has** — `none` — the page that describes the benchmark panel in detail is
  the one page in the section without a picture of it.
- **verdict** — **SEQUENCE** (+ **STALE-RISK**) — **deferred until the sibling
  rendering-measurement unit releases the app.**
- **needs** — A five-step scripted sequence, all crops confined to
  `#benchmark-section` (the accordion whose `<summary>` reads **Performance
  Benchmark**), viewport 1440 × 900, `deviceScaleFactor: 2`, files under
  `_static/screenshots/benchmarking_performance/`:
  1. **The controls, before any run.** Crop `#benchmark-section` from the
     `Quick test:` label through `#benchmark-run`, including the six
     `.benchmark-preset` buttons (`100K`, `500K`, `1M`, `5M`, `10M`, `20M`),
     `#benchmark-count`, and `#benchmark-pattern`. Cursor on the `1M` preset
     (`button[data-count="1000000"]`).
     `benchmark-controls.png` · `:width: 246px`.
  2. **Running it.** Cursor on `#benchmark-run` (`Load Synthetic Data`).
     Same crop. `benchmark-run-button.png` · `:width: 246px`.
  3. **The live stats.** Crop `#benchmark-stats` — the six-cell grid (`Points`,
     `FPS`, `Frame Time`, `GPU Memory`, `LOD Level`, `Visible Points`), the
     `#bench-timing-details` block (`min:` / `p95:` / `max:`), and
     `#bench-gen-info` (`Generated in N ms`). No cursor.
     `benchmark-live-stats.png` · `:width: 246px`. The page's "How to interpret
     the stats" section walks these fields one by one and shows none of them.
  4. **Copy Situation Report.** Click `#benchmark-report-btn`, then crop
     `#benchmark-report-output` (the read-only `<textarea>`) with its text
     populated. Cursor on `#benchmark-report-btn`.
     `benchmark-situation-report.png` · `:width: 246px`. **Redaction:** the
     report may embed a GPU/driver string and the page URL — read the textarea
     content before publishing and confirm no `localhost:PORT` or local path is
     present.
  5. **Analyze Performance.** Click `#bottleneck-analyze-btn`, wait for
     `#bottleneck-results`, crop as specified for `benchmark-verdict-panel.png`
     under {doc}`/user_guide/web_app/n_benchmarking_performance/index` — one file serves both pages. Optionally a second frame
     with the `Show detailed stats` `<details>` expanded, exposing
     `#bn-lod-overhead`, `#bn-frustum-overhead`, `#bn-point-size-response`,
     `#bn-frame-stability`, `#bn-jank-percent`, `#bn-cpu-health`, since the page
     tells readers to act on the "LOD overhead" verdict and on point-size
     response. There is no `#bn-shader-overhead`; the shader figure and its
     bottleneck were removed as unmeasurable by that instrument.
     `benchmark-verdict-detailed.png` · `:width: 246px`.
- **notes** — **Gap in the prose.** The "Data patterns (how to choose)" section
  describes `Uniform Random`, `Gaussian Clusters`, `Batch Effects`,
  `Atlas-like`, and `Flat UMAP` — but the `#benchmark-pattern` select ships with
  **`Model Surface` (`glb`) selected by default**, and also offers `Octopus` and
  `3D Spirals`. The default option a user will actually run is the one the page
  never mentions. Fix the prose in the same change as the captures, or the
  screenshot of `#benchmark-pattern` will visibly contradict the text beside it.
  Also: "live stats (FPS, frame time, GPU memory, LOD)" understates the panel —
  it also reports `Points`, `Visible Points`, min/p95/max frame times, and
  generation time.

## {doc}`/user_guide/web_app/n_benchmarking_performance/06_edge_cases_performance`

- **explains** — Eleven performance cliffs in what-you-see / why / confirm / fix
  / prevention form.
- **has** — `none`
- **verdict** — **ADD** (two images only — most of this page is not picturable)
- **needs** — Nine of the eleven cliffs are *sensations* (choppy, hot, throttled)
  or *configurations* already described in text. Exactly two produce a distinct
  on-screen artifact worth photographing:
  1. **WebGL context lost.** The viewer builds a modal overlay with the title
     `WebGL Context Lost`, the body `WebGL context lost. Reload required.`, and
     a `Reload` button (the context-lost overlay in
     `assets/js/rendering/viewer.js`). Provoke it in a
     scripted run with the WebGL debug extension
     `WEBGL_lose_context.loseContext()` on `#glcanvas` — `render-failure-state.spec.mjs`
     already reaches into `document.querySelector('#glcanvas')` for a comparable
     purpose, so the technique is established. Crop the overlay panel only.
     Cursor on the `Reload` button.
     `perf-webgl-context-lost.png` · `:width: 520px`. Serves this page **and**
     {doc}`/user_guide/web_app/n_benchmarking_performance/07_troubleshooting_performance`, where the same symptom has its own
     section. Do **not** provoke context loss by loading 20M points — that
     collides with the rendering-measurement unit's work.
  2. **Category explosion.** Colour by a categorical field with a very large
     category count and crop `#legend` showing the unusable legend, with
     `#coloring-filtering-section` framing it. No shipped sample has such a
     field; a scratch prepared fixture derived from
     `tests/browser/fixtures/exports/current-ui-prepared` with a synthetic
     high-cardinality column would be required.
     `perf-category-explosion-legend.png` · `:width: 246px`. Lower priority than
     the context-lost image — it needs a new fixture.
- **notes** — "Grid changes are applied after the approximately 300 ms slider
  debounce" (line 129) is a precise claim about the smoke-mode controls that was
  not verified in this pass; it lives in the render-mode module outside this
  section. Flagging as unverified.

## {doc}`/user_guide/web_app/n_benchmarking_performance/07_troubleshooting_performance`

- **explains** — Seven symptom → cause → confirm → fix → prevention entries plus
  a template for adding more.
- **has** — `none`
- **verdict** — **ADD**
- **needs** — Reuse, do not re-shoot:
  - Symptom "WebGL context lost / canvas goes blank" → `perf-webgl-context-lost.png`.
    This page quotes the message in prose; showing the actual dialog with its
    `Reload` button is the fix, not just the diagnosis.
  - Symptom "Filtering is slow / sliders lag" → `perf-live-filtering-off.png`.
  - Symptom "Fine in single view but slow in grid view" → `perf-view-controls.png`
    with the cursor on `#split-clear-btn`.
  - Symptom "Slow on one machine but fine on another" → `benchmark-controls.png`
    plus `benchmark-live-stats.png`, since the fix is "run the synthetic
    benchmark with the same preset and compare FPS".
  - Not picturable, keep as prose: thermal throttling, extension interference,
    driver differences.
- **notes** — The console one-liner
  `document.createElement("canvas").getContext("webgl2") !== null` appears here
  and on `06`; it is correct and browser-neutral. A DevTools console screenshot
  would be browser chrome — do not capture it.

## {doc}`/user_guide/web_app/n_benchmarking_performance/08_screenshots`

- **explains** — The section's verified-capture gallery: the benchmark panel
  plus a four-row cross-engine acceptance table, and two multiview frames.
- **has** — 3 figures.
  - Line 5, `benchmarking_performance/benchmark-panel.png`, `:width: 1440px`
    (line 7) — the section's only owned image.
  - Line 32, `web_app/multiview-two-panels.png`, `:width: 1440px` (line 34) —
    **borrowed**.
  - Line 40, `web_app/dark-theme-multiview.png`, `:width: 1440px` (line 42) —
    **borrowed**, and its caption is about *theme consistency*, which is not a
    performance topic at all.
- **verdict** — **REPLACE** (+ **STALE-RISK**)
- **needs** —
  - Drop `dark-theme-multiview.png`. Its stated subject — "Theme and background
    changes apply consistently across the active comparison layout" — belongs to
    `c_core_interactions`, not to a performance gallery. It is padding.
  - Keep `multiview-two-panels.png` but re-caption it to the view-count cost the
    surrounding heading claims, or replace it with a genuine 1-view vs 4-view
    pair captured at the same point count with the FPS readout visible in both —
    that would make the "each view multiplies work" claim measurable rather than
    asserted. Files `perf-viewcount-1.png` / `perf-viewcount-4.png` ·
    `:width: 1440px` each. **Deferred** with the rest of section N.
  - Promote the five new benchmark-panel captures from `05` into this gallery so
    the section owns its images.
- **notes** — **Stale-risk by construction.** The caption for
  `benchmark-panel.png` states "Build 2026-07-26.4 … 17 FPS, a 57.67 ms frame
  sample, full LOD, all 10.0M points visible, generation in 1,592 ms", and the
  table beneath adds four more engine rows. Those numbers are baked into both
  the caption and the pixels. Any re-capture of that panel invalidates the
  caption, and any re-run on different hardware invalidates the table. The page
  does hedge correctly ("acceptance evidence, not universal performance
  claims"), which is the right posture — but whoever re-shoots this image must
  re-derive every number in the caption and the table in the same change. The
  same caption is duplicated on `02`, so both copies drift together.

## {doc}`/user_guide/web_app/n_benchmarking_performance/09_reporting_performance_bugs`

- **explains** — A pre-filing checklist, three report kinds, a copy/paste
  template, how to gather machine/GPU info, minimal reproductions, and a
  redaction checklist.
- **has** — `none`
- **verdict** — **ADD** (one image)
- **needs** — The template asks the reporter to "Paste **Copy Situation Report**
  from the **Performance Benchmark** panel", and nothing in the docs shows where
  that button is or what it produces. Reuse
  `benchmark-situation-report.png` (spec under `05`) — the populated
  `#benchmark-report-output` textarea with the cursor on
  `#benchmark-report-btn`. One image, placed in the "Measurements" subsection.
  - The "Analyze Performance result" bullet is served by
    `benchmark-verdict-panel.png` if a second image is wanted.
  - Nothing else here is picturable: the template is a form and the redaction
    checklist is a policy.
- **notes** — This page's redaction checklist ("remove local file paths (often
  contain usernames), remove private org/repo names, remove tokens/URLs with
  secrets") is the same policy this audit applies to its own captures under
  `l_sessions_sharing/08`. The situation-report image is the one capture in
  section N that could itself violate that checklist — its textarea content must
  be read before publishing.

---

## Cross-cutting findings

1. ~~**The accordion is named `Session`, not `User data`.** Sections `01`, `03`
   and `05` all instruct readers to open "User data".~~ **Done** — all three
   name **Session**, and the string `User data` no longer occurs anywhere in
   section L. The element is the `<details id="session-section">` block in
   `index.html`; "Local data:" is a label on `#user-data-block` *inside* it.

2. **There is no Session panel status line.** `#session-state-controls` contains
   a label, an info button, and two buttons — nothing else. `03` ("…for the
   panel status to say **Session loaded successfully**") and `10` ("The Session
   panel then reports…", "Check the panel status…") describe a panel status that
   does not exist; every one of those strings is a NotificationCenter toast.
   **Still open.**

3. **{doc}`/user_guide/web_app/l_sessions_sharing/09_edge_cases` still describes the pre-`cellOrder` world.** Its "Same
   lightweight fingerprint, changed content" section says re-ordering cells
   passes the identity guard. Every other page in the section says the opposite,
   and the source agrees with the other pages. This is the most serious prose
   defect found.

4. **Three pages assert a cause the shipped error message deliberately refuses
   to assert.** `01` (line 142), `07` (line 64), `10` (line 111) all say the
   third rejection means the cells are stored in a different order. The message
   names both possible causes because a one-way digest cannot distinguish a
   permutation from a re-computed embedding, and the source carries an explicit
   comment about why claiming otherwise would be an integrity failure of its own.

5. ~~**A stale SHA-256 ships in {doc}`/user_guide/web_app/l_sessions_sharing/04_official_sample_states`.**~~ **Done** —
   the page shows the field's shape and points at `exports/datasets.json`. The
   general rule this settled: **do not copy a regenerated value into prose.**
   Digests, preset byte sizes, and source line numbers all live somewhere that
   is rebuilt without the docs being touched; cite the file, the symbol, or the
   element id instead. Three literals in this record had gone stale by the time
   the findings were actioned.

6. **The default synthetic benchmark pattern is undocumented.**
   `#benchmark-pattern` defaults to `Model Surface` (`glb`); {doc}`/user_guide/web_app/n_benchmarking_performance/05_benchmark_tools`
   documents five other patterns and omits the default, plus `Octopus` and
   `3D Spirals`.

7. ~~**Resolution convention.** Existing screenshots are 1×; moving to 2×
   capture with 1× `:width:` values must be done as one sweep.~~ **Adopted** and
   owned by ``00_capture_tooling_and_conventions``, which states the viewport,
   the device scale factor, and the emitted-width rule. No count of existing
   screenshots is recorded here — the number changes with every capture run.
