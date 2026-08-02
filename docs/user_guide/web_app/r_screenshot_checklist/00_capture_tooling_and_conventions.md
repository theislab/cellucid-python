# Capture tooling and figure conventions

Read this before capturing a single image. It assumes you have not seen this
codebase.

This directory is **not published**. It is excluded from the Sphinx build in
``docs/conf.py``, because it is working material for whoever is adding figures,
not documentation for a reader of the app.

---

## 1. Where the tooling is and why

``cellucid-python/docs/_tooling/screenshots/``

| File | What it is |
| --- | --- |
| ``capture.mjs`` | The engine and the command line. |
| ``scenarios.mjs`` | The registry of named, reproducible scenarios. |
| ``overlay.mjs`` | The pointer and annotation layer, injected into the page. |
| ``serve-exports.mjs`` | A loopback host for a prepared dataset export tree. |
| ``captures.json`` | Generated provenance: one record per captured image. |

It lives in the **Python repository** rather than the web repository for three
reasons. The artifact it produces — ``docs/_static/screenshots/**`` — lives
here. The workspace guide puts web-app documentation in this repository. And the
contract that governs the images (see section 5) is a Python test in this
repository, so the tool and the rule that binds it stay together.

It has no dependencies of its own. Playwright is resolved out of the web
repository's ``node_modules`` at run time. **Never run ``npm install``** — the
browsers and the Playwright build are already provisioned.

---

## 2. Invoking it

```sh
cd cellucid-python/docs/_tooling/screenshots

node capture.mjs list                                   # the registry
node capture.mjs shoot <id> [<id> ...] --out /some/dir  # write PNGs to a scratch dir
node capture.mjs shoot --all --out /some/dir
node capture.mjs shoot <id> --emit                      # write into _static/screenshots/
node capture.mjs check <id> [<id> ...]                  # recapture and diff against captures.json
```

``shoot`` prints, for every image, the exact ``{figure}`` block to paste,
already carrying the correct ``:width:``.

### Ports

The tool starts the web repository's own browser-test server and reads the same
two variables that server does, plus one of its own:

| Variable | Default | Serves |
| --- | --- | --- |
| ``CELLUCID_BROWSER_TEST_PORT`` | 4173 | the application |
| ``CELLUCID_BROWSER_TEST_SAMPLE_PORT`` | main + 1 | the permissive-CORS mirror |
| ``CELLUCID_DOCS_EXPORTS_PORT`` | main + 2 | ``cellucid-datasets/exports/`` |

**Check the defaults are free before you run.** A browser-test run or another
capture run may already hold them, and the server refuses to bind rather than
adopting someone else's. Move all three together:

```sh
CELLUCID_BROWSER_TEST_PORT=4213 \
CELLUCID_BROWSER_TEST_SAMPLE_PORT=4214 \
CELLUCID_DOCS_EXPORTS_PORT=4215 \
node capture.mjs shoot --all --out /tmp/shots
```

### Data

Scenarios load the **real published exports** from
``cellucid-datasets/exports/``, served over loopback. That is deliberate: a
reader has to recognise the picture as the app they are using, and a 120-cell
synthetic fixture does not look like single-cell data. The tree is pinned on
disk, so it is reproducible in the way that matters.

Available ids: ``suo``, ``garcia``, ``he``, ``kanemaru``, ``pancreas``.

**Prefer ``pancreas`` unless the page is about scale.** 3,696 cells, 3,753
genes, 1D/2D/3D embeddings, connectivity, vector fields, two genuinely
continuous observation fields, and eight named cell types. It loads in seconds
and every analysis mode completes on it. ``suo`` is 562k cells and will make
every capture slow for no editorial gain.

---

## 3. Writing a scenario

A scenario is a plain object in ``scenarios.mjs``. Everything a reviewer needs
in order to judge whether the image shows the right thing is in that one object,
and the same object is what regenerates the image a year later.

```js
{
  id: 'filtering-active-stack-two-filters',
  topic: 'filtering',
  name: 'active-filters-two-rows',
  describes: 'The Active Filters list with two filter rows, ...',   // becomes :alt:
  dataset: 'pancreas',
  steps: [ /* ordered actions */ ],
  shot: { mode: 'element', selector: '#active-filters', padding: 6 },
  cursor: { on: '#active-filter-scope-btn', at: 'center', state: 'hover' },
  rings: [{ on: '#active-filters', label: '1' }],
}
```

### Step vocabulary

| Step | Effect |
| --- | --- |
| ``{ waitFor: sel }`` | wait until visible |
| ``{ waitForHidden: sel }`` | wait until hidden |
| ``{ waitForText: sel, text }`` | wait until the element's text contains ``text`` |
| ``{ waitForStable: sel, forMs }`` | wait until the element's markup stops changing |
| ``{ click: sel }`` / ``{ hover: sel }`` | pointer actions |
| ``{ press: 'Escape' }`` | keyboard |
| ``{ select: sel, label }`` or ``{ select: sel, value }`` | choose an option, then **assert it stuck** |
| ``{ fill: sel, value }`` | type into a field |
| ``{ openSection: sel }`` / ``{ closeSection: sel }`` | set a sidebar accordion, idempotently |
| ``{ scrollIntoView: sel }`` | scroll the containing scroller |
| ``{ wait: ms }`` | last resort; recorded in ``captures.json`` as a determinism weakness |
| ``{ run: async page => {} }`` | raw Playwright escape hatch |

### Shot modes

- ``{ mode: 'window' }`` — the whole 1440×1000 viewport.
- ``{ mode: 'element', selector, padding }`` — crop to one element's bounding
  rectangle.
- ``{ mode: 'region', selectors: [...], padding }`` — crop to the union of
  several, for a control group that spans siblings.

**A crop is always a selector, never a pixel box.** The rectangle is read at
capture time, so a control that moves takes its crop with it. The failure mode
of a hand-tuned box — silently framing the wrong thing after a layout change —
cannot happen.

### The pointer

``cursor: { on, at, dx, dy, state }``. ``at`` is one of ``center``, ``left``,
``right``, ``top``, ``bottom``; ``state`` is ``hover`` or ``press``.

The tool performs a **real** ``page.hover()`` at that point and *then* draws the
arrow there. Both, or neither: an arrow floating over an idle button is worse
than no arrow, because it shows an interaction state the app never entered.
``press`` additionally holds the mouse down, so ``:active`` styling is real, and
adds a ring at the hotspot — a click is an event, and a still cannot show an
event without a mark for it.

**Use a cursor when the prose says "click", "drag", "hover", or "open".** Do not
use one on a pure state or result image; it implies an action that is not being
described.

### Annotation

``rings: [{ on: selector, label: '1' }]`` draws an orange rounded outline around
an element, with an optional numbered badge.

Off by default, on purpose. On a cropped panel the crop *is* the annotation and
a ring adds noise. On a full-window shot that discusses one control, a ring is
the difference between a picture and an instruction. Numbered badges are for the
lifecycle sequences, where the reader is following steps.

---

## 4. What determinism actually rests on

Fixed and recorded for every capture: viewport 1440×1000, device scale factor 2,
light colour scheme, ``en-US``, UTC, reduced motion, Playwright's
``animations: 'disabled'`` and ``caret: 'hide'``, and a **seeded
``Math.random``** installed before any application code runs (the particle
overlay and parts of the analysis layer draw from it).

Two things do the real work:

**A fresh browser context per scenario.** The app persists theme, sidebar
geometry, and onboarding state. A shared context makes the result depend on what
ran before it.

**The settle loop.** The tool recaptures the clip region until two consecutive
frames are byte-identical before writing anything. A WebGL canvas is not settled
because a promise resolved; it is settled when it stops changing. Every capture
reports ``settled`` or ``NOT SETTLED``, and ``captures.json`` records it.

For a scene that never settles — an animated particle overlay — declare
``settle: { frames: N }``. The capture is then pinned to N animation frames from
a seeded start, and is recorded as **not** byte-stable. Do not pretend otherwise.

### Provenance: capture from a clean web repository

``captures.json`` records ``web_repository_revision`` — the application the
image was served from. The server reads the web repository **from disk**, so
that value describes the working tree, not ``HEAD``:

| Recorded value | Meaning |
| --- | --- |
| ``<sha>`` | the tree was clean; this commit rebuilds the pictured app |
| ``<sha>-dirty`` | uncommitted changes were served; **no commit contains what the image shows** |
| ``<sha>-unknown-worktree`` | the revision resolved but its cleanliness did not |
| ``unknown`` | no revision at all |

Only the first form is provenance. ``shoot`` and ``check`` both print a warning
when the tree is dirty, and ``check`` prints a note when a *stored* record
carries a suffix — bytes agreeing proves an image is reproducible, but it cannot
recover which application produced it.

**Commit or stash the web repository before capturing**, and be aware that
several units may be editing it at once. A capture taken over somebody else's
half-finished change is stamped ``-dirty`` and has to be taken again.

---

## 5. Size and resolution — the one rule, and the contract behind it

> **Emitted width = the captured region's CSS width × 2, capped at 1440 px.**

The tool computes this. You do not choose it, and you should not override it
with ``emitWidth`` to make an image "look right" on one page.

Why those numbers:

- The documentation theme lays the article column out at ``max-width: 60em``
  with ``1rem`` padding — **928 CSS px** — and ``docs/_static/custom.css`` gives
  every figure ``max-width: 100%; height: auto``. So an image *wider* than 928 px
  is scaled down by the browser and its extra pixels become sharpness on a
  high-DPI display. An image *narrower* than 928 px renders 1:1 and gains
  nothing.
- Capture happens at device scale factor 2, so a region is captured at twice its
  CSS size. Emitting that natively is what makes panel crops sharp.
- The 1440 cap keeps a full-window shot at 1.55× oversampling in the column
  instead of 3.1×, which would cost megabytes per file for no visible gain — and
  1440 is the width the eighteen existing full-window captures already use, so
  full-window shots stay drop-in compatible.

**The visible consequence.** A 246 CSS px sidebar panel used to be emitted at
246 px and rendered at 246 px. That is why the existing panel figures are soft
and too small to read. The same panel now emits at 492 px and renders at 492 px:
twice the size on the page, and sharp. Panel figures will look noticeably
different after a sweep. That is the fix, not a regression.

### The contract you cannot violate

``cellucid-python/tests/test_docs_current_api_contract.py``,
``test_every_documented_screenshot_uses_its_intrinsic_responsive_width``:

1. Every figure line must start with exactly ```` ```{figure} ```` at column
   zero.
2. The target must resolve inside ``docs/_static/screenshots/`` and must exist.
3. Exactly one ``:alt:`` (non-empty) and exactly one ``:width:``.
4. **``:width:`` must equal the PNG's intrinsic pixel width**, written as
   ``<N>px``. There is no freedom here — the emitted width *is* the declared
   width. ``shoot`` prints the correct block; paste it.
5. **``referenced == screenshots``.** Every PNG under
   ``_static/screenshots/`` must be referenced by at least one figure, and every
   figure must reference an existing PNG. This is an **equality**, so:

> **An orphan PNG fails ``python -m pytest``.** You cannot land an image in one
> change and its ``{figure}`` block in the next. They go together, always.

That is why ``shoot`` requires ``--out DIR`` unless you pass ``--emit``. Iterate
into a scratch directory; ``--emit`` only when the prose is ready. ``--emit``
prints a reminder to that effect for every file it writes.

The scenarios with topic ``_proof`` exist only to demonstrate that the tool
works. **Never ``--emit`` them** — they would land as orphan PNGs in a directory
no page references.

``captures.json`` and ``scenarios.mjs`` are **shared state**. More than one unit
edits them at once. Append your scenarios; do not rewrite either file wholesale,
and re-read before you edit.

Only ``--emit`` writes to ``captures.json``. A ``--out`` iteration and a
``check`` run leave it alone, because provenance describes the *published*
bytes: a record written from a scratch image points at a file nobody can see,
and ``check`` would then compare the published image against a digest it was
never supposed to have.

---

## 6. Adding a figure to a page

Paste the printed block and replace the caption line:

    ```{figure} ../../../_static/screenshots/<topic>/<name>.png
    :alt: <one sentence describing what is visible>
    :width: <N>px

    One sentence saying what the reader should take from the image.
    ```

(The block above is indented by four spaces so that the docs contract test does
not parse this page's example as a real figure. **Remove the indentation when
you paste it.**)

- ``../../../`` is correct from a page at
  ``docs/user_guide/web_app/<section>/<page>.md``. Count the levels if you are
  somewhere else.
- ``:alt:`` describes what is **visible**. ``:width:`` comes from the tool.
- The caption says what to conclude. Do not repeat the alt text.

Then run both checks (see section 8). A figure edit can break the Python suite
and a docs edit can break the build.

---

## 7. Traps, every one of which cost time to find

**The Python contract suite reads every ``.md`` under ``docs/``.** Not just
published pages — ``DOCS_ROOT.rglob("*.md")``, with only ``_build``,
``jupyter_execute``, ``.ipynb_checkpoints`` and ``__pycache__`` excluded.
``exclude_patterns`` in ``conf.py`` does **not** exempt a file from pytest. Two
consequences:

- ``test_doc_cross_reference_contract.py`` fails on a single-backtick span
  shaped like a page name — a bare `NN_page_slug` in single backticks, or a
  `<section>/index` path in single backticks. Write it as a ``{doc}`` role, or as a
  double-backtick span, or append ``.md``. The files in this directory were
  mechanically converted to double backticks for exactly this reason.
- Several tests assert that a phrase is **absent** from the concatenation of all
  documentation markdown. Quoting an obsolete UI string in a note can fail a
  test in a file you have never opened. Run the suite.

**The field selectors are rebuilt after the point count settles.** The
observation manifest arrives later than the geometry, and the rebuild silently
discards a selection made in that window — ``selectOption`` reports success and
the control comes back holding a different value, so you get an image captioned
"coloured by clusters" that shows ``Field: None``. Always:

```js
{ waitForText: '#filter-count', text: 'Showing all 3,696 points' },
{ waitForStable: '#categorical-field' },
{ select: '#categorical-field', label: 'clusters' },
{ waitForText: '#stats', text: 'Field: clusters (category)' },
```

The tool now asserts the selection stuck and fails loudly if it did not.

**``#dataset-cells`` is not a ready signal.** It abbreviates to ``4K``. Use
``#filter-count`` (``Showing all 3,696 points``) or ``#stats``
(``Points: 3,696 • Field: clusters (category) • Centroids: 8``).

**``scrollIntoView`` inside the sidebar changes what a full-window shot shows.**
It scrolls the sidebar's own scroller, so the panel you wanted at the top is
gone and the reader sees the app mid-scroll. This is precisely why the existing
``web_app/app-overview-cell-type.png`` is captioned as showing a legend that is
not in frame. For a window shot, use ``closeSection`` to collapse what you do
not need rather than scrolling to what you do.

**The welcome overlay must be waited for before Escape is sent.** It paints
before its key handler is wired, and an early Escape is swallowed, leaving the
overlay covering every later capture. The tool does this for you; a raw
Playwright script must too.

**The ✕ floating over the top-left of the canvas is real UI**, not an artifact.
It is ``#sidebar-toggle`` (``title="Toggle sidebar"``), which renders as ✕ while
the sidebar is open. Do not "fix" it, and do not crop it out as though it were a
capture defect. If a figure would confuse a reader with it, say what it is in the
caption.

**Sidebar accordions are ``<details class="accordion-section">``.** Clicking the
``<summary>`` *toggles*, so a scenario that ran twice would leave the panel
wrong. ``openSection`` / ``closeSection`` set ``.open`` directly and are
idempotent. Open on first load: session, visualization, compare-views,
coloring-filtering, highlighted-cells, page-analysis. Closed:
community-annotation, cinematic-camera, figure-export, shortcuts, benchmark.

**Community-annotation dev hooks are host-gated.** The simulation flags the
browser suite uses are behind ``isLocalDevHost()``. Pointed at anything other
than ``127.0.0.1`` or ``localhost`` they do nothing and you get disconnected
panels **silently**. The capture servers bind ``127.0.0.1``, so this is already
correct — do not change the host.

**Never leak.** No local paths, usernames, tokens, or private repository names
in any committed image. The browser's download shelf and the OS file picker
render real paths and cannot be sanitised — do not capture them.

**Never photograph a terminal. Transcribe it.** Terminal output ships as a
verbatim ```` ```text ```` block, never as an image. This is a rule, not a
preference, and it is not satisfied by a neutralised prompt:

- **The leak is in the output, not the prompt.** ``cellucid serve`` over a
  prepared export prints the **absolute** directory it resolved —
  ``server/_server.py`` calls ``print_detail("Path", str(self.data_dir))`` — and
  AnnData's own ``INFO`` logging prints the absolute file path beside it. Both
  are correct behaviour and neither can be configured away, so *every*
  screenshot of these commands contains a host path by construction. A
  neutralised ``PS1`` hides the prompt and changes nothing about the transcript.
- **An SSH example leaks twice over.** ``ssh -L 8765:localhost:8765
  user@remote-host`` carries a username and a hostname in the command itself.
- **Rendered text cannot be audited.** ``strings`` over a PNG finds nothing,
  because the characters are pixels. Two images shipped with a home directory
  and a pyenv path in them and no scan could see it; they were found only by a
  person opening every file.
- **Text is better anyway** — greppable, diffable, copy-pasteable, kept honest
  by the contract tests, and a few hundred bytes instead of a few hundred
  kilobytes.

Elide the environment-specific lines using the placeholder convention the guides
already teach (``/path/to/…``) and say in the prose that the real run prints an
absolute path there. ``tests/test_data_loading_docs_contract.py`` asserts that
no page under ``docs/`` names a home, scratch, or pyenv path, and that the
server transcripts are ``text`` blocks rather than figures.

A screenshot is still the right medium for a **graphical** application window —
JupyterLab, the app itself — where layout and colour carry the meaning.

---

## 8. Validation

From ``cellucid-python`` with ``source ../.environments/python/venv/bin/activate``.
**Never install anything.**

```sh
sphinx-build -W --keep-going -b html docs docs/_build/html
python -m pytest
```

``-W`` means a broken image reference fails the build. ``pytest`` enforces the
width and orphan rules in section 5. Run both after any figure change; neither
alone is sufficient.

To confirm previously captured images still match the app:

```sh
node capture.mjs check --all
```

A difference is a signal to look, not automatically a failure — the app may have
legitimately changed. Recapture, **open the image**, and confirm it still shows
what the caption claims. An unopened screenshot is unverified.
