# Screenshot checklist for the web-app guide

Working material, not published. ``docs/conf.py`` excludes this directory from
the Sphinx build.

**Start with ``00_capture_tooling_and_conventions.md``.** It has the tool, the
size rule, the contract the Python suite enforces, and every trap found while
building it. Nothing here is useful until you have read that.

## What is in this directory

| File | Sections | Pages | Form |
| --- | --- | ---: | --- |
| ``00_capture_tooling_and_conventions.md`` | how to capture, size rules, traps, validation | — | current |
| ``01_pages_orientation_data_loading.md`` | A orientation · B data loading · G · Q troubleshooting | 20 | **assignment** |
| ``02_pages_interactions_fields_filtering.md`` | C core interactions · D fields/colouring/legends · E filtering | 25 | status record |
| ``03_pages_selection_analysis_accessibility.md`` | F highlighting/selection · H analysis · O accessibility/privacy | 23 | **assignment** |
| ``04_pages_velocity_annotation_export.md`` | I velocity · J community annotation · K figure export | 23 | **assignment** |
| ``05_pages_sessions_performance.md`` | L sessions/sharing · N benchmarking/performance | 23 | **assignment** |
| ``06_pages_developer.md`` | P developer docs | 20 | **assignment** |

Each per-page record in an **assignment** file carries: the page path, what it
explains, what figure it had *when the audit was written*, a verdict, and — for
anything that needed work — the dataset, the interaction sequence, the crop
target as a DOM concept, whether a cursor helps and on which element, and a
suggested filename.

:::{warning}
**The five files still marked "assignment" describe the tree as it was at
``f304688``, before the capture sweep that same commit delivered.** Their
``has`` records and their counts are pre-work and are now largely wrong — treat
their *specifications* as useful and their *inventories* as historical. Measure
the current state before acting on one. ``02`` has been brought to a measured
status record; the other five have not.
:::

## Where coverage stands

Measured 2026-08-02 across the 136 pages under ``user_guide/web_app/``:
**311 figure directives** referencing **155 distinct images**. At audit time the
same tree carried 86 directives over far fewer images, because one file was
pasted onto page after page.

| Section group | Pages | Figure directives | Pages with none |
| --- | ---: | ---: | ---: |
| A · B · G · Q | 20 | 46 | 8 |
| C · D · E | 25 | 103 | 2 |
| F · H · O | 23 | 71 | 4 |
| I · J · K | 24 | 57 | 8 |
| L · N | 23 | 34 | 12 |
| P (developer) | 20 | 0 | 20 |

Where the gaps now are:

- **``p_developer_docs`` — 20 pages, still zero images.** Correctly so for most
  of them: nine want a diagram, four a terminal or devtools capture, five nothing
  at all, and exactly one wants a picture of the app.
- **``l_sessions_sharing`` and ``n_benchmarking_performance``** carry the thinnest
  coverage of the illustrated sections — 34 directives over 23 pages, 12 of them
  with none.
- **Provenance, not count, is the live risk.** 86 of the 155 published images
  have no record in ``captures.json`` and therefore cannot be re-verified by
  ``node capture.mjs check``. That absence is what let an empty screenshot ship.
  ``tests/test_docs_current_api_contract.py::test_screenshot_provenance_records_describe_the_published_images``
  holds that number where it is; it may only decrease.

## The four lifecycles

The user named community annotation and analysis. Two more are the same shape.
Each is specified step by step in its section file: state before, action, state
after, crop target, cursor position, filename.

| Lifecycle | File | Steps | Notes |
| --- | --- | ---: | --- |
| Selection | ``03_pages_selection_analysis_accessibility.md`` | 17 | idle → lasso drag → combine → undo → KNN → proximity → annotation → a named page reused downstream |
| Analysis | ``03_pages_selection_analysis_accessibility.md`` | 19 | two named pages → six modes → configure → run → result → expanded modal → export |
| Community annotation | ``04_pages_velocity_annotation_export.md`` | 18 | disconnected → connect → authenticated → browse → propose → vote → consensus → publish |
| Session save and restore | ``05_pages_sessions_performance.md`` | 14 (+1 failure) | configure a non-trivial state → save → the artifact → reset → restore → a comparable pair |

Three things to carry from them:

**Selection step 17 is deliberately the same frame as analysis step 1**, so the
two capture in one session with agreeing page names and counts.

**All 18 community-annotation steps are reachable without a real GitHub
repository.** Thirteen need no network at all; five need route interception
against a synthetic worker. The section file names the mechanism per step, and
names the two sub-states that are genuinely unreachable rather than inventing
them. Never point this at a real repository or account.

**The session round-trip is proved by a pair**, not a single image: the section
file lists exactly which items must be identical between the before and after
frames, and the one difference that is permitted.

## What must not be screenshotted

Recorded here so nobody spends a run discovering it:

- **``g_cross_highlighting``** is an empty directory and is not in any toctree.
  The pages were deleted deliberately because they documented an unbuilt feature;
  the module directory in the app holds a placeholder file and nothing else.
  Commission no images for it. The two empty directories on disk are the only
  residue.
- **Flow rate in the velocity overlay** cannot be shown in a still. Do not fake a
  before/after pair for it.
- **Screen-reader behaviour** is not screenshot-shaped. Focus rings and the
  shortcuts overlay are; announcements are not.
- **The browser download shelf and the OS file picker** render real local paths
  and cannot be sanitised. Where an artifact must be shown, use a prepared
  terminal instead.
- Several performance topics are numbers in a table, not pictures. The section
  file says which.

## Defects found while auditing

The audit was supposed to count figures. It also turned up documentation and app
defects, which are recorded in the section files with page and line numbers.
They are listed there rather than fixed, because fixing them is not this
directory's job — but several of them **change what a screenshot must show**, so
read the relevant section file before capturing, not after.

The ones that block a capture:

- A widely reused overview image is captioned as showing a legend that is not in
  frame. The cause is a scrolled sidebar (see the ``scrollIntoView`` trap in the
  conventions file). Recapturing it without fixing the framing reproduces the
  defect at higher resolution.
- Two analysis images show gene identifiers that no published sample dataset
  produces any more. They cannot be recaptured as they are; the scenario has to
  change.
- Two figure-export images predate controls that three pages describe in prose.
- Several captions hard-code an exact build string that has already moved.
  Remove build strings from captions rather than chasing them. **Done for
  ``a_orientation/03`` and ``04``**, whose captions pinned ``Build 2026-07-27.1``
  while the app reported ``2026-07-27.24`` — the footer identity bumps on every
  web release, so a pinned literal documents a build nobody is running. The
  contract test now forbids the pattern rather than pinning a value.
- One documented keyboard modifier for adding to a selection contradicts both the
  rest of the guide and the renderer. Establish which is right before
  illustrating it.
- A page promises a diagram and has an empty gap where it should be. No diagram
  extension is enabled in the build today; the developer-section file names what
  would have to be added and argues for a generated-from-source approach, since a
  hand-drawn architecture diagram rots faster than a screenshot.

## Suggested order of work

1. **The conventions file, in full.** Then capture one image and open it.
2. **Community annotation and the analysis lifecycle**, because the user named
   them and because both need a mock-state harness that everything else can then
   reuse.
3. **The session lifecycle**, the largest bare section.
4. **The re-shoot sweep.** Every existing image is a 1× capture. Re-shooting them
   is mechanical, needs no editorial decision, and doubles the on-page size of
   every panel crop — do it as one sweep so the site does not carry two visual
   generations at once.
5. **Velocity, selection, filtering** — the per-control gaps.
6. **The developer section last**, and mostly not as screenshots.

Retire the per-section gallery pages as you go. Ten sections have one. They
strand images away from the prose that explains them, and several of them consist
entirely of copies of figures already shown on the topic pages.
