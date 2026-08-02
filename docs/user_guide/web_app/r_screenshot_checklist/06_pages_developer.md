# Figure-coverage audit — `user_guide/web_app/p_developer_docs/`

Section audited: 20 `.md` files (19 numbered pages + `index.md`). Every page read in full.
Web app source cross-checked against `cellucid/` (1575 tracked files enumerated).

## Summary

| Metric | Value |
|---|---|
| Pages in section | 20 |
| Pages with an existing figure | **0** (verified: zero `{figure}`, zero `{image}`, zero `![]` in the section) |
| Fenced code blocks in section | 12 language-tagged (6 `bash`, 4 `js`, 1 `html`, 1 `{toctree}`) + 2 untagged ASCII art blocks |

### Verdict counts

| Verdict | Count | Pages |
|---|---|---|
| DIAGRAM | 9 | 01, 05, 06, 07, 09, 10, 11, 12, 16 |
| NONE | 5 | `index.md`, 00, 15, 17, 18 |
| TERMINAL | 2 | 02, 14 |
| DEVTOOLS | 2 | 03, 13 |
| SCREENSHOT | 1 | 08 |
| STALE-RISK | 1 | 04 |

The distribution confirms the premise: **one** page in twenty wants an actual picture of the running app.

### Diagram mechanism: NOT AVAILABLE TODAY

Verified directly, not inferred:

- `cellucid-python/docs/conf.py` → `extensions = [sphinx.ext.autodoc, sphinx.ext.autosummary, sphinx.ext.viewcode, sphinx.ext.intersphinx, sphinx.ext.napoleon, sphinx_autodoc_typehints, sphinx_design, myst_nb, sphinx_copybutton]`. **No mermaid, no graphviz, no plantuml.**
- `cellucid-python/pyproject.toml` → `[project.optional-dependencies] docs = [sphinx>=7.2,<8.1, myst-nb>=1.1.0, sphinx-autodoc-typehints>=1.25, pydata-sphinx-theme>=0.15, sphinx-design>=0.5, sphinx-copybutton>=0.5]`. **No diagram package.**
- `grep -rni "mermaid|graphviz|sphinxcontrib|plantuml" docs/ --exclude-dir=_build` → **zero hits**. The docs tree has never contained a diagram directive.
- `.readthedocs.yaml` installs `.[docs]` only; no `build.apt_packages`.

So today a diagram can only ship as a committed image file. Three ways forward, in the order I would rank them:

1. **`sphinxcontrib-mermaid`** — add `"sphinxcontrib-mermaid>=0.9"` to the `docs` extra in `pyproject.toml` **and** `"sphinxcontrib.mermaid"` to `extensions` in `conf.py`. Renders ```` ```{mermaid} ```` blocks; works with MyST out of the box; diagram source is plain text living next to the prose, so it diffs and reviews like code. This is the only option where a wrong diagram shows up in a PR diff. Note the `-W` build: a malformed mermaid block becomes a build failure, which is the desired behaviour.
2. **`sphinx.ext.graphviz`** — bundled with Sphinx, so **no new pip dependency**, but it shells out to `dot`. That needs `build.apt_packages: [graphviz]` added to `.readthedocs.yaml` and `graphviz` installed on every contributor's machine. More moving parts than option 1 for no benefit here.
3. **Pre-rendered SVG committed to `_static/`** — works today, zero config, but it is exactly the hand-drawn artifact the brief warns about. Acceptable only for the two figures below that describe byte layouts and repo boundaries (things that change on a schema break, not on a refactor).

**The rot problem, and the fix that actually holds.** Nine pages want a diagram, and most of those diagrams are module/import graphs — precisely the artifact that is wrong one refactor later. The repo already has the pattern for this: `cellucid/scripts/validate-token-types.js` and `cellucid/scripts/validate-tokens.js` are CI-run generators/validators for the CSS design system, and `cellucid/tests/` is full of `*-exact-contract.test.mjs` files that assert source facts. The same shape applies here: a `cellucid/scripts/generate-architecture-diagrams.mjs` that walks the real ES-module import graph and emits mermaid source, plus a node contract test asserting the committed output matches a fresh run. That makes a stale diagram a **failing test**, not a reader's problem. Every DIAGRAM verdict below names the source file(s) the generator should read.

### DEAD PATHS AND DEAD IDENTIFIERS

Two categories, both verified against source. Dead identifiers matter as much as dead paths on a config-reference page — a reader who types them gets silence, not an error.

**Dead filesystem paths (1)**

| Page:line | Claim | Reality |
|---|---|---|
| `01_codebase_map_and_entry_points.md:63` | repo-map tree lists `assets/exports/  # (Optional) demo datasets, session snapshots, etc.` | No `exports/` directory exists anywhere under `cellucid/assets/` or at any depth ≤3 in `cellucid/`. `assets/` contains exactly `css external fonts img js`. |

Every other backticked path in all 20 pages resolves. I extracted 255 path-shaped tokens and checked each against the real tree under all plausible base directories (`cellucid/`, `assets/js/`, `assets/js/app/`, `assets/js/app/ui/`, …); 234 resolved, and the 21 non-resolving tokens are all non-paths (MIME types `text/javascript`/`text/plain`, session chunk IDs `core/state`/`analysis/windows`, the URL template `owner/repo/path`, HTML element names `input/select`, and the deliberately hypothetical `my-feature-controls.js`). I also spot-checked all 8 Sphinx `{doc}` cross-reference targets that leave the section — all 8 resolve.

**Dead identifiers (9)** — named in prose, absent from source

| Page:line | Claim | Reality |
|---|---|---|
| `04:65` | `?ga_debug=1` enables GA debug mode, "see `.../ui/core/ga-init.js`" | The file exists; `ga_debug` appears **nowhere** in `cellucid/` (repo-wide grep, excluding `node_modules`). `ga-init.js` contains no `debug` handling at all. |
| `04:86` | localStorage key `cellucid:community-annotations:last-repo-map` | Zero occurrences repo-wide. |
| `04:87` | localStorage key `cellucid:community-annotations:repo-meta` | Zero occurrences in `assets/js`. The only hit repo-wide is `tests/community-annotation-wire-contract.test.mjs:3293`, which asserts the key is **not** written: `assert.equal(values.has('cellucid:community-annotations:repo-meta'), false)`. The docs document a key a contract test proves is gone. |
| `04:90` | localStorage key `cellucid:community-annotations:last-github-user-key` | Zero occurrences repo-wide. |
| `04:103` | sessionStorage key `cellucid:github-app-auth:token:v1` | Zero occurrences. The real key is a single `cellucid:github-app-auth:session` (`assets/js/app/community-annotations/github-auth.js:22`). |
| `04:104` | sessionStorage key `cellucid:github-app-auth:user:v1` | Zero occurrences. There is no separate user key; identity rides in the one `:session` key above. |
| `04:118` | global `window.__CELLUCID_DEV__` ("Dev helper namespace (self-tests, etc.)") | Zero occurrences in `assets/js`. The only hit repo-wide is `tests/dataset-controls-exact-contract.test.mjs:535`, in a test literally named *"dataset controls contain no dev global, reload, or guessed-source path"*: `assert.doesNotMatch(moduleSource, /__CELLUCID_DEV__/)`. |
| `05:54` | `createViewer({ canvas, labelLayer, viewTitleLayer, sidebar })` | The `sidebar` parameter does not exist. Real signature is `createViewer({ canvas, labelLayer, viewTitleLayer, onViewFocus })` (`assets/js/rendering/viewer.js:536`), and the sole call site passes only three (`assets/js/app/main.js:250`). Page 07:28 documents the same function **correctly** — the two pages contradict each other. |
| `14:160` | `window.__CELLUCID_DEV__.runLocalUserInPlaceSwitchSelfTest(...)` | Zero occurrences in `assets/js`. `tests/dataset-controls-exact-contract.test.mjs:536` asserts its absence: `assert.doesNotMatch(moduleSource, /runLocalUserInPlaceSwitchSelfTest/)`. The entire section `14:154–167` ("Dev-only dataset reload self-test") documents removed behaviour. |

---

## `index.md`

**explains** — Section landing page: a "Fast Path (Pick Your Goal)" routing table mapping seven developer goals to a start page and a follow-up page, plus a scope note and a hidden globbed toctree.

**has** — `none`. Verified.

**verdict** — **NONE**

**needs** — Nothing. The page is a routing table of `{doc}` links; the table *is* the navigational artifact and a diagram of it would restate seven rows as seven boxes. The one figure that could belong here — a section-wide reading-order graph — is already encoded by the table's two columns and would need re-drawing every time a page is added.

**notes** — All internal `{doc}` targets resolve, including the two that leave the section (``../index``, ``../b_data_loading/index``, ``../o_accessibility_privacy_security/index``). Line 31's description of the codebase map is consistent with page 01.

---

## `00_how_to_use_this_section.md`

**explains** — Reading-depth guidance (three audience layers with prescribed page sequences), the repo separation between `cellucid/`, `cellucid-python/`, and `cellucid-annotation/`, a "definition of done" checklist for web app changes, and the documentation conventions the rest of the section uses.

**has** — `none`. Verified.

**verdict** — **NONE**

**needs** — Nothing. Three of the four sections are ordered lists of `{doc}` links and policy prose; neither renders as a figure. The fourth ("Key repo separation") is the only candidate — a three-box diagram of the repositories — but it is three bullet points that a reader absorbs faster as text, and the authoritative version of that relationship already exists as the repository table in the workspace `CLAUDE.md`. Adding a picture here duplicates a table without adding information.

**notes** — Line 26 says `cellucid/` holds "community annotation UI, session bundles" — accurate. Line 61's example code pointer `cellucid/assets/js/app/main.js` resolves. Line 63 ("Format and runtime contracts describe only the current exact schema") is consistent with the clean-break policy and with page 10's actual content. No defects found on this page.

---

## `01_codebase_map_and_entry_points.md`

**explains** — Where the web app code lives (layer breakdown + an ASCII repo tree) and what executes first (`index.html` → `main.js` → `ui-coordinator.js` → `data-state.js`), plus a "where do I look for X" navigation table and key terminology.

**has** — `none`. Verified.

**verdict** — **DIAGRAM** — and both diagrams are generatable from source. This is the single highest-value figure target in the section.

**needs** — Two figures.

*(a) Directory tree — generate it, do not draw it.* The page already carries a hand-maintained ASCII tree at lines 46–64, and it is already wrong (see notes). Replace it with output from a generator walking `cellucid/assets/js/` to a fixed depth, emitting one line per directory plus the four entry-point files, with the one-line purpose annotations kept in a small mapping table inside the generator. Nodes: `index.html`, `assets/css/`, `assets/js/{app,data,rendering,utils,analytics,dev}`, and under `app/`: `main.js`, `state/`, `ui/`, `session/`, `state-serializer/`, `analysis/`, `community-annotations/`, `registries/`, `utils/`, `dockable-accordions/`, `notification-center/`, `url-state.js`. Source of truth: the filesystem itself, gated by a contract test in `cellucid/tests/` asserting the committed tree matches a fresh walk.

*(b) Boot / entry-point call graph.* Directed, top-to-bottom, six nodes and seven edges. Nodes: `index.html` → `app/main.js` → {`rendering/viewer.js::createViewer`, `app/state/index.js::createDataState`, `app/ui/core/ui-coordinator.js::initUI`, `data/data-source-manager.js::DataSourceManager`, `app/session/index.js::createSessionSerializer`, `app/analysis/comparison-module.js::createComparisonModule`}. Edge labels are the factory function names. Derived from: the top-level `import` statements of `cellucid/assets/js/app/main.js` filtered to the six factories the page already names at lines 87–94. A generator reading `main.js`'s import list keeps this honest.

**notes** — Real defects on the highest-value verification target:

- **`:63` — dead path.** `assets/exports/` does not exist. `cellucid/assets/` contains exactly `css external fonts img js`. Delete the line.
- **`:46–64` — tree is materially incomplete.** Missing real directories: `assets/js/analytics/` (`performance.js`, `tracker.js`), `assets/js/dev/` (`benchmark.js`), `assets/external/`, `assets/fonts/`, `assets/img/`. Missing real `app/` children: `registries/`, `utils/`, `dockable-accordions/` + `dockable-accordions.js`, `notification-center/` + `notification-center.js`, `url-state.js` (which page 04 and page 05 both reference), `startup-failure.js`, `startup-url-intent.js`, `jupyter-command-handler.js`. Missing top-level repo dirs a developer needs: `tests/`, `scripts/`, `types/` — all three are referenced by page 14 and page 16.
- **`:96–99` — dev globals list is incomplete.** Documents `window._cellucidViewer`, `window._cellucidState`, `window._comparisonModule`. Source also sets `window._cellucidUi` (`main.js:2314`) and `window._cellucidBenchmarkHarness` (`main.js:3520`).
- **Duplication / drift risk.** `cellucid/assets/js/app/README.md` contains a second hand-maintained ASCII tree of the same `app/` layer. The two already disagree (the in-repo README lists `registries/` and `app/utils/`, this page lists neither). Two hand-maintained maps of one tree is one map too many — generating (a) resolves it, or the page should link the in-repo README instead of restating it.
- All 14 code-pointer paths in the "Where do I look for…?" table (`:125–133`) resolve, including the two `.md` targets (`figure-export/README.md`, `community-annotations/REPO_SETUP.md`). The four entry-point files at `:76–117` all resolve.

---

## `02_local_development_setup.md`

**explains** — How to run the web app locally: the two distinct "local" modes (frontend vs data server), the static-server command, the debug-logging toggle, four dataset-loading options (demo / file picker / Python server / Jupyter), and boot troubleshooting.

**has** — `none`. Verified.

**verdict** — **TERMINAL**

**needs** — Two terminal captures, plus one devtools capture.

*(1) The static server, exact command:*
```
cd cellucid
python -m http.server 8000
```
Must be visible in the capture: the `Serving HTTP on 0.0.0.0 port 8000 (http://0.0.0.0:8000/) ...` banner **and** at least three subsequent request lines proving modules are being served over HTTP — `"GET / HTTP/1.1" 200 -`, `"GET /assets/js/app/main.js HTTP/1.1" 200 -`, and one `assets/js/rendering/viewer.js` line. The request log is the point: it is what distinguishes this from the `file://` failure the page warns about at `:47–51`.

*(2) The Python data server, exact command:*
```
cellucid serve /path/to/dataset.h5ad --dataset-name "My dataset" --dataset-id my-dataset
```
Must be visible: the printed viewer URL containing the chosen port (the page tells the reader to enter `http://127.0.0.1:<port>` at `:114`, so the capture must show where that port number comes from). Crop after the first few log lines.

*(3) Devtools, supporting:* Network tab with **Disable cache** checked, matching the instruction at `:136–140`.

**Leak risk — real, must be handled.** Capture (2) prints an absolute filesystem path containing the operator's home directory and username, and a dataset name. Use a synthetic path under a neutral root and a demo dataset. Capture (1) is safe (only `cd cellucid` and localhost). Do not capture a shell prompt that embeds hostname/username; set `PS1='$ '` before recording.

**notes** — `:95` says the file picker can select "a prepared `exports/` folder" — this is the *user's* prepared export directory, not the repo path from page 01, so it is correct and is **not** a dead path. `:72` (`assets/js/utils/debug.js`), `:84` (`index.html`), and the `{doc}` target at `:128` all resolve. Options A–D are consistent with page 09's source list. The `document.createElement("canvas").getContext("webgl2") !== null` check at `:169` is a valid one-liner.

---

## `03_build_run_and_deployment.md`

**explains** — That Cellucid is a static ES-module app with no required bundler; how to serve it; a four-item deployment checklist (MIME types, module paths, CSP, no catch-all rewrites); and four dataset-hosting hazards (CORS, `.gz` double-decompression, HTML-for-404, session-bundle proxies).

**has** — `none`. Verified.

**verdict** — **DEVTOOLS**

**needs** — One capture, and it is the highest-signal figure on the page. **Panel: Network.** Click a failing `assets/js/app/main.js` request and show the **Response Headers** pane with `Content-Type: text/plain` visible, alongside the Console tab's resulting error text (`Failed to load module script: Expected a JavaScript module script but the server responded with a MIME type of "text/plain"`). Ideally a two-panel side-by-side: the broken host on the left, a correct `Content-Type: text/javascript` on the right. This converts the page's most common failure mode (`:40–42`, restated at `:111–123`) from a sentence into something a reader can pattern-match against their own devtools.

A second, optional capture: the Network row for a `.bin.gz` asset showing `Content-Encoding: gzip` present — the misconfiguration described at `:76–84`. Lower priority; the header name alone is enough.

**Leak risk — low.** HTTP response headers on a static asset carry no personal data. Redact any `Server:` header naming an internal hostname, and any `Set-Cookie`.

**notes** — Stability is good: HTTP response headers are the most durable thing that could possibly be screenshotted in this section. `:49` (`index.html` CSP meta) resolves and `index.html` does carry a `Content-Security-Policy` meta. `:51`'s claim about the JSON-LD CSP hash needing an update is a real constraint worth keeping. No dead paths. No stale prose found.

---

## `04_configuration_env_vars_and_feature_flags.md`

**explains** — Every configuration surface the frontend has: URL deep-link parameters (dataset/source/remote/anndata/github/annotations/debug), `localStorage` keys, `sessionStorage` keys, and developer-only `window.__…` global overrides.

**has** — `none`. Verified.

**verdict** — **STALE-RISK**

**needs** — **No figure until the page is corrected.** This is the deliberate call, not a dodge. The obvious figure — devtools **Application → Local Storage / Session Storage** showing the `cellucid*` keys — is exactly the artifact that would rot fastest here: it pins a screenshot to a key inventory that changes whenever a feature is refactored, and it would immediately contradict the tables on the page, because **seven of the documented keys/params do not exist**. A screenshot cannot be added to a reference table that is wrong; fixing the table is the prerequisite, and once the table is right the reader gains nothing from a picture of it.

If a figure is added *after* the corrections land, restrict it to a single Application-panel capture of the Local Storage pane filtered to `cellucid`, showing only the four keys that are stable and locally verifiable (`CELLUCID_DEBUG`, `cellucid_theme`, `cellucid_viewer_background`, `cellucid_last_quote_index` — these are the only four written via a direct `localStorage.setItem('<literal>')` in the whole app). **Leak risk is severe** for the community-annotation and GitHub keys: `cellucid:github-app-auth:session` holds an OAuth token and the user identity payload. Never capture Session Storage on this page. Never capture a signed-in state.

**notes** — Seven dead identifiers, all listed in the section summary above: `?ga_debug=1` (`:65`), `…:last-repo-map` (`:86`), `…:repo-meta` (`:87`, contract test asserts absence), `…:last-github-user-key` (`:90`), `cellucid:github-app-auth:token:v1` (`:103`), `cellucid:github-app-auth:user:v1` (`:104`), `window.__CELLUCID_DEV__` (`:118`, contract test asserts absence).

The keys that **are** correct, verified: `CELLUCID_DEBUG` (`utils/debug.js`), `cellucid_theme` (`utils/theme-manager.js:14`), `cellucid_viewer_background` (`render-controls.js:37`), `cellucid_last_quote_index` (`onboarding/welcome-modal.js`), `cellucid:community-annotations:repo-map` (`community-annotations/*.js:15`), `cellucid:community-annotations:lock:*` (`scope-lock.js:27` `LOCK_PREFIX`), `cellucid:community-annotations:cache:*:files:shas` (`cache-scope.js:98`), `cellucid:community-annotations:cache:*:session` (`cache-scope.js:92`), `cellucid:community-annotations:tab-id:v1` (`:26`), `window.CELLUCID_DEBUG`, `window.__CELLUCID_DEBUG__` (`analysis/shared/debug-utils.js:62`), `?debug=1` and `localStorage.debug` (same file, `:43`/`:56`), `window.__CELLUCID_GITHUB_WORKER_ORIGIN__` (`github-auth.js`).

One omission: the URL-parameter table (`:27–33`) does not list `?exportsBaseUrl=`, even though it is a real query key (`EXPORTS_QUERY_KEY = 'exportsBaseUrl'`) and both page 02 (`:85`) and page 09 (`:64`, `:94`) tell readers to use it. It belongs in this table — this is the page that claims to be the complete parameter inventory.

---

## `05_app_architecture_overview.md`

**explains** — The big-picture mental model: an ASCII boot→state→UI→rendering map, the seven-step startup sequence, the five data sources, what `DataState` owns and emits, how the UI is wired, the three distinct persistence mechanisms, four "hard" performance invariants, and how to write a fixable bug report.

**has** — `none`. Verified.

**verdict** — **DIAGRAM** — generatable.

**needs** — One diagram, replacing the ASCII block at `:24–43` (which is already doing the job badly and is already wrong at `:54`).

Directed graph, three horizontal bands, top to bottom.

- **Band 1 (boot):** `index.html` → `app/main.js`. Edge label: `<script type="module">`.
- **Band 2 (subsystems main.js constructs):** six nodes fanning out from `main.js` — `rendering/viewer.js` (createViewer), `app/state/core/data-state.js` (createDataState), `app/ui/core/ui-coordinator.js` (initUI), `data/data-source-manager.js`, `data/data-loaders.js`, `app/session/session-serializer.js`, `app/analysis/*`.
- **Band 3 (the runtime edges that matter):** `DataState` → `viewer` with three labelled edges `updateColors` / `updateTransparency` / `updatePositions`, and `DataState` → `ui-coordinator` with a dashed edge labelled `emit: visibility:changed, field:changed, highlight:changed, page:changed, dimension:changed, vectorFields:changed`. The dashed-vs-solid distinction (calls vs events) is the whole architectural point of the page and is the thing ASCII cannot express.

Derived from: `cellucid/assets/js/app/main.js` (the import list gives band 2 exactly), `cellucid/assets/js/rendering/viewer.js` (the three `update*` method names), and the event-name list, which is independently stated in `cellucid/assets/js/app/README.md` under "Common state events" — so the generator can cross-check the docs against the in-repo README rather than a hand list.

**notes** —

- **`:54` — dead identifier.** `createViewer({ canvas, labelLayer, viewTitleLayer, sidebar })`. There is no `sidebar` parameter. Real signature: `createViewer({ canvas, labelLayer, viewTitleLayer, onViewFocus })` (`assets/js/rendering/viewer.js:536`); the sole call site passes three arguments (`assets/js/app/main.js:250`). Page 07 `:28` documents the same call correctly — fix 05 to match 07.
- Every other code pointer resolves (`:73`, `:95`, `:115`, `:130`, `:131`, `:132`, `:144`, `:152`), as do both external `{doc}` targets.
- The six event names at `:104–109` all match `cellucid/assets/js/app/README.md`. The four `viewer.*` methods at `:117–120` are real. The boot sequence at `:51–65` is consistent with `main.js`.
- `:78`'s description of `local-demo` (exports base URL via meta tag or `?exportsBaseUrl=`) matches page 09 `:93–95` and the source constant. Good cross-page consistency here.

---

## `06_state_datastate_and_events.md`

**explains** — The state model in depth: the core typed arrays and their stable-index invariant, bounded LRU field caches, multi-view contexts, the multi-page highlight system, how the public API is assembled from mixins and manager getters, the six events and their payloads, batch mode, state→viewer synchronization, and four classes of subtle bug.

**has** — `none`. Verified.

**verdict** — **DIAGRAM** — generatable.

**needs** — Two diagrams.

*(a) Ownership / composition graph.* Center node `DataState`. Incoming edges from five **mixed-in method surfaces** (`DataStateViewMethods`, `DataStateFieldMethods`, `DataStateColorMethods`, `DataStateFilterMethods`, `highlightStateMethods`) drawn as one edge style, and five **manager getters** (`state.views → ViewManager`, `state.fields → FieldManager`, `state.filters → FilterManager`, `state.colors → ColorManager`, `state.highlights → HighlightManager`) drawn as another. Outgoing edges to the five owned typed arrays (`positionsArray`, `colorsArray`, `categoryTransparency`, `outlierQuantilesArray`, `highlightArray`) and the two LRU caches. Derived from: `cellucid/assets/js/app/state/core/data-state.js` — the page already enumerates the exact symbol names at `:104–113`, so a generator can parse the mixin wiring directly from that file and fail CI when a manager is added or renamed.

*(b) Event-flow sequence diagram, one lane per participant.* Lanes: *UI module* → *manager* (`filter-manager.js`) → *typed array* (`categoryTransparency`) → *viewer* (`updateTransparency`) → *emit* (`visibility:changed`) → *ui-coordinator* (debounced) → *three UI re-renders* (filterControls, highlight summary, view badge). This is the exact chain the page describes twice in prose (`:191–195` and page 08 `:137–143`) and it is the mechanism behind the "why did this fire twice?" question at `:147–158`. A sequence figure is strictly better than prose here because the ordering and the debounce are the content.

**notes** — This page verifies clean. All eight code pointers resolve, including the two glob forms (`state/managers/view-context-*` at `:203` matches `view-context-viewer-sync.js`; `state/managers/*-manager.js` at `:204` matches). The six events and the `field:changed` payload shape `{ source, fieldIndex, changeType, detail }` at `:133` match `cellucid/assets/js/app/README.md` exactly. `beginBatch()`/`endBatch()` (`:168–169`) are corroborated by page 05 `:167`. `assets/js/app/utils/event-emitter.js` (`:125`) exists. No dead paths, no stale prose found.

---

## `07_rendering_pipeline_webgl_and_performance_notes.md`

**explains** — How the WebGL2 renderer draws millions of points: the viewer entry point and its WebGL2 hard requirement, six render subsystems (scatter, smoke/volumetric density, connectivity edges, highlight tools, centroids, overlays), the viewer's public update API, multiview rendering, vector-field overlay coupling, three performance footguns, and two GL failure modes.

**has** — `none`. Verified.

**verdict** — **DIAGRAM**

**needs** — One diagram, and it is unusually well-specified because the source already states it. **The smoke/volumetric density GPU transaction**, as a left-to-right pipeline: `bounded position batches` → `R32F splat atlas` → `logarithmic GPU maximum reduction` → `normalized R8 atlas` → `R8 3D texture`. Annotate the first edge with the batch cap (262,144 XYZ points / 3 MiB) and the atlas node with the 128³ worst case (1536×1408, under WebGL2's 2048 `MAX_TEXTURE_SIZE` floor). Add a side branch showing the three preconditions (`clean GL error state`, `EXT_color_buffer_float`, `EXT_float_blend`) gating entry, and a failure edge showing "prior volume preserved" — the atomicity property is the page's central claim and a diagram is where it becomes legible.

Derived from: `cellucid/assets/js/rendering/smoke-cloud/smoke-density.js` — the pipeline is stated verbatim in a source comment at line 658 (`R32F splat atlas -> max reduction -> normalized R8 atlas -> 3D slice copies`), and the two constants are exported (`MAX_SPLAT_POINTS_PER_BATCH = 262_144` at `:16`; `MAX_SMOKE_GRID_SIZE = 128` in `smoke-density-contract.js:1`). A generator or a contract test can assert the figure's numbers against those exports.

A second, lower-priority diagram: the six render subsystems as boxes under `viewer.js` with their backing files. Derived from `viewer.js` imports.

**Explicitly not wanted here: a devtools Performance panel capture.** Frame timings are machine-, GPU-, and dataset-dependent; a flame chart screenshot would be non-reproducible and would invite readers to compare their numbers against an arbitrary one. The page's timing claims are already correctly framed as "not a GPU-completion measurement" (`:66–68`) and should stay prose.

**notes** — This is the most technically precise page in the section and it **verifies clean against source**, which is worth stating explicitly since it is the page most likely to be assumed stale:

- `:56` "1536×1408 at 128³, below WebGL2's required 2048 minimum" — atlas dimensions are computed as `gridSize * slicesPerRow` × `gridSize * numRows` (`smoke-density.js:748–749`) and validated against `maximumTextureSize` (`:766`). Consistent.
- `:55` "accepts exact integer sizes from 8 through 128 and rejects values outside that range before allocation" — exactly matches `smoke-density.js:690–695` (`gridSize < 8 || gridSize > MAX_SMOKE_GRID_SIZE` → `RangeError`), with `MAX_SMOKE_GRID_SIZE = 128`.
- `:55` "The exact UI grid inventory is 32³, 48³, 64³, 96³, and 128³" — correct, and worth guarding: `render-controls.js:826` contains a *different* array `[32, 48, 64, 96, 128, 192, 256]`, but that is the **noise texture resolution** slider (`updateNoiseResolutionSlider`), not the density grid. The density grid comes from `sliderToGridSize` / `SMOKE_GRID_SIZES` (`:769`). The docs are right; a future editor comparing against the wrong array would "correct" them into being wrong.
- `:164` "transient buffer fixed at 262,144 XYZ points (3 MiB)" — matches the exported constant.
- `:27–28` `createViewer({ canvas, labelLayer, viewTitleLayer, onViewFocus })` — **correct**, unlike page 05 `:54`.
- All seven code pointers resolve, including the two globs and `rendering/overlays/velocity/velocity-overlay.js`.

---

## `08_ui_modules_map.md`

**explains** — The UI layer map: the thin `initUI` coordinator, the centralized DOM cache and its domain buckets, a five-table index of sidebar modules grouped by concern, how modules communicate via state events rather than reaching into each other, which UI state sessions exclude and how, and a module-boundary troubleshooting entry.

**has** — `none`. Verified.

**verdict** — **SCREENSHOT** — the one page in this section where a picture of the running app is the right artifact.

**needs** — An **annotated sidebar screenshot**: the app's full sidebar with every accordion section outlined and each outline labelled with its owning module file. This is the "UI module map benefits from an annotated app screenshot with each module region labelled" case exactly — the page's five tables map *module file → what it owns*, and the reader's actual question is the inverse: "I clicked this thing, which file is it?" Only a labelled picture answers that direction.

- **Dataset:** one of the built-in demo datasets from the sample picker (Pancreas — small, fast, and every sidebar section populates). Do not use a local-user dataset; it puts a filesystem path in the dataset info block.
- **Interaction:** load the dataset, then **expand every accordion section** so each labelled region has visible content rather than a collapsed summary row. This will exceed one viewport — capture the sidebar as a tall single image (devtools device-mode with a tall viewport, or two stacked crops), not a scrolled screenshot with cut-off sections.
- **Crop target:** the sidebar column only, full height, plus ~40px of canvas on the right edge for context. Exclude the browser chrome and the URL bar (the URL carries `?dataset=` and possibly `?exportsBaseUrl=`).
- **Cursor:** none — no hover state. Hover would highlight one section and imply it is the subject.
- **Annotation:** outline + callout label per region, labels reading e.g. `ui/modules/filter-controls.js`. Labels in a consistent accent colour that survives both the light and dark docs theme.

Secondary (cheaper, complements it): a **DIAGRAM** of the coordinator's import graph, generated from `cellucid/assets/js/app/ui/core/ui-coordinator.js`'s import list — that graph is what would keep the module tables from drifting again.

**notes** — Real defects on the section's second-highest-value verification target. Every documented file exists — the failure here is **omission and a wrong wiring claim**, not dead paths.

- **`:63–115` — the module index covers 19 of 29 real entries.** The five tables list 19 modules; `cellucid/assets/js/app/ui/modules/` contains 29 entries. The 10 undocumented ones: `benchmark/`, `camera-input-contract.js`, `cinematic-camera/`, `community-annotation/` (a whole directory, distinct from the documented `community-annotation-controls.js`), `community-annotation-modal-owner.js`, `community-annotation-voting-modal.js`, `cross-highlight/`, `field-interaction-owner.js`, `highlight/`, `legend/`. Two of these are user-visible features with no entry anywhere in the section: **Cinematic camera** and **Benchmark**. `cinematic-camera` is not a minor internal — page 10 `:90` lists `cinematic/camera` as a mandatory session chunk, so the section documents its *persistence* while never introducing the feature.
- **`:66–69` — the wiring claim is imprecise.** "Most UI modules live under `ui/modules/` … They are imported/initialized in `ui-coordinator.js`." Three modules in the page's own tables are **not** imported by the coordinator: `dataset-connections.js` is imported by `dataset-controls.js:28`; `field-selector-gene-expression.js` and `field-selector-deleted-fields.js` are imported by `field-selector.js`. Conversely the coordinator imports two things absent from every table: `../modules/cinematic-camera/index.js` and `../components/info-popovers.js` (note: `ui/components/`, a sibling directory the page never mentions, alongside `ui/category-builder/`, `ui/onboarding/`, and `ui/keyboard-move.js`).
- `assets/js/app/ui/modules/cross-highlight/` contains only a `.gitkeep` — an empty directory. Not a docs defect; flagged because a reader who trusts a completed module map would look for it, and because it is a candidate for deletion under the clean-break rule.
- Everything else checks out: `ui-coordinator.js` and `dom-cache.js` resolve; `collectDOMReferences` is the real exported name (`ui-coordinator.js:13`); the `data-state-serializer-skip="true"` mechanism and both state-serializer pointers (`:159`, `:160`) resolve; the exclusion list at `:152–156` is consistent with page 10 `:327–337` and page 12 `:93–107`; the analysis note at `:117–121` is consistent with page 11 `:26–44`.

---

## `09_data_loading_pipeline_and_caching.md`

**explains** — The data-loading pipeline: the two-layer "source selection → base URL → file loading" model, the boot-time source registration and URL-override order, five data sources with their files and failure modes, the loader layer's artifact inventory, four distinct caching layers, and a six-step loading-failure checklist.

**has** — `none`. Verified.

**verdict** — **DIAGRAM**

**needs** — Two figures.

*(a) Source fan-in → baseUrl → loader fan-out.* Left column: five source nodes, each labelled with its file — `local-demo` (`data/local-demo-source.js`), `github-repo` (`data/github-data-source.js`), `local-user` (`data/local-user-source.js`), `remote` (`data/remote-source.js`), `jupyter` (`data/jupyter-source.js`). All five converge on a single node `DataSourceManager` (`data/data-source-manager.js`) emitting `{ sourceType, datasetId, baseUrl, metadata }`. That feeds one node `data/data-loaders.js`, which fans out to the artifact list: `dataset_identity.json`, `obs_manifest.json`, `var_manifest.json` (dashed = optional), embeddings/points buffers, `connectivity_manifest.json` + edges (dashed). The bottleneck shape is the page's entire thesis ("the base URL is the critical bridge", `:40`) and prose cannot draw a bottleneck. Derived from: `data/data-source-manager.js` registration list plus the five source modules' filenames — a generator can list `assets/js/data/*-source.js` directly.

*(b) Caching-layer stack.* Four stacked bands with what lives in each and who invalidates it: HTTP cache (browser, `{cache:'force-cache'}`), in-memory source/manifest caches (`DataSourceManager`), `DataState` bounded LRU field caches (obs ≤50 / var ≤20 entries — numbers from page 06 `:57–58`), IndexedDB + localStorage sha index (`community-annotations/file-cache.js`). Four bands is exactly the kind of thing a reader must hold simultaneously to debug "why is my changed file not reloading", and it is stable — these four layers are architectural, not implementation detail.

*Supporting devtools capture (optional, medium value):* Network waterfall of one real demo-dataset load, filtered to XHR/Fetch, showing the actual request order the page claims at `:69–75`. **Leak risk low** — the exports base URL is a public host — but crop the URL bar and any `Authorization` header.

**notes** — Page verifies clean. All eleven code pointers resolve; both `{doc}` targets resolve. `:93–95`'s two-step `EXPORTS_BASE_URL` resolution order (query field, then meta tag) matches the source constant `EXPORTS_QUERY_KEY = 'exportsBaseUrl'` and page 02 `:83–85`. `:114`'s `owner/repo/path` is a user-input template, not a filesystem path. The four source interface methods at `:81–84` and the two health endpoints at `:150–151` are plausible and consistent with the rest of the section. `DecompressionStream` handling at `:203–205` is consistent with page 03 `:76–84`. No dead paths, no stale prose found.

---

## `10_sessions_persistence_and_serialization.md`

**explains** — The `.cellucid-session` contract in full: design invariants, the code map, the exact container framing and manifest schema, the closed inventory of static and dynamic chunks, capture ordering, the ten-step restore pipeline with rollback, reversible large-state ownership, the bounded-gzip preflight, the dataset fingerprint and `cellOrder` digest with its four distinguished refusals, the official five-chunk sample profile, cancellation outcomes, and the discipline for changing persistence.

**has** — `none`. Verified.

**verdict** — **DIAGRAM**

**needs** — Two diagrams. Both describe byte-level and control-flow contracts, which is the one category of diagram that does **not** rot with refactors — it changes only when the format breaks, and a format break is exactly when the docs must be rewritten anyway.

*(a) Container byte-layout figure.* A horizontal byte-ribbon, left to right: `"CELLUCID_SESSION\n"` (magic) │ `u32 LE manifest length` │ `manifest JSON (UTF-8, exactly that many bytes)` │ then a repeating group `u32 LE stored length` │ `stored payload bytes` × N chunks. Below the ribbon, an expansion of the manifest object showing its exactly-three keys (`createdAt`, `datasetFingerprint`, `chunks`), the fingerprint's exactly-five keys, and a chunk's exactly-nine keys. Derived from: `cellucid/assets/js/app/session/bundle/format.js` (framing) and `session/session-context.js` (`getDatasetFingerprint`). This turns `:51–74` — currently three dense paragraphs a reader must reconstruct mentally — into one glance.

*(b) Restore pipeline with the rollback edge.* Ten numbered steps down the page as stated at `:127–139` (supersede → frame under limits → validate manifest/profile/fingerprint → register reversible snapshots → eager decode/apply → lazy decode/apply → prepare all participants → commit in order → final UI refresh → publish "Session fully restored"), with a **red return edge** from any step back through "rollback in reverse participant order". The all-or-nothing property (`:24`, `:142–143`) is the page's most important claim and the reverse edge is the only way to show it. Derived from `session/session-serializer.js`.

**Explicitly not a screenshot.** There is a tempting figure here — the Session panel showing "Session fully restored" — but it is a toast notification, it belongs in the user-facing sessions chapter, and it would illustrate the least interesting sentence on the page.

**notes** — The most precise page in the section, and it verifies clean. All eight code-map pointers resolve, including `session/bundle/format.js`, `session/codecs/gzip.js`, and the `session/contributors/` directory. The chunk IDs flagged as unresolved by a naive path check (`core/state`, `analysis/windows`, `cinematic/camera`, `highlights/meta`, `ui/dockable-layout`, `core/field-overlays`, `analysis/cache-inventory`) are **session chunk identifiers, not filesystem paths** — correctly rendered in backticks, no defect. `/^[0-9a-f]{16}$/` at `:232` is a regex, likewise not a path. The `cinematic/camera` chunk at `:90` is corroborated by the real `ui/modules/cinematic-camera/` module (see page 08 notes) — though that module is undocumented in the UI map, so this page is currently the only place a reader learns the feature exists.

---

## `11_analysis_architecture.md`

**explains** — The Page Analysis subsystem: why it is initialized in `main.js` rather than by the UI coordinator, the "analysis runs on pages + highlights" model, the analysis directory code map, five core pieces (`ComparisonModule`, `DataLayer`, GPU/Worker/CPU compute backends, Plotly plot infrastructure, memory monitor + cleanup), session interaction, three correctness edge cases, and two troubleshooting entries.

**has** — `none`. Verified.

**verdict** — **DIAGRAM** — generatable.

**needs** — One diagram: the analysis subsystem's ownership graph.

Root `app/main.js` → `createComparisonModule({ state, container })` mounting into DOM id `#page-analysis-section`. From `comparison-module.js`, seven child nodes matching the real subdirectories: `analysis/core/` (plugin contracts, validation, registries), `analysis/data/` (`data-layer.js`, query builder, transforms), `analysis/compute/` (expand into four leaves: `operations.js`, `compute-manager.js`, `gpu-compute.js`, `worker-pool.js`), `analysis/stats/`, `analysis/plots/` (leaves `plot-factory.js`, `plotly-loader.js`, `types/`), `analysis/ui/` (leaf `analysis-ui-manager.js`, `analysis-types/`), `analysis/shared/` (leaves `memory-monitor.js`, `resource-cleanup.js`, `plot-theme.js`). Add one cross-edge `DataState --page:changed, highlight:changed--> ComparisonModule` (dashed, matching the event style used in the page 05 and 06 diagrams) and one `sessionSerializer.setAnalysisRefs()` edge.

Critically, annotate the three compute backends as **mutually exclusive selection, not a substitution chain** — the page states at `:130–131` that a backend is chosen before the operation starts and execution errors never switch backends. A naive reader assumes a substitution chain; the diagram should make the single-selection rule visible (single selection arrow, not a chain).

Derived from: `cellucid/assets/js/app/analysis/index.js` (the barrel export, which the page describes at `:74–75` as documenting the core building blocks — it is the natural machine-readable source) plus `comparison-module.js`'s import list.

**notes** — Page verifies clean. All 17 code pointers resolve, including every subdirectory at `:81–87` and both shared utilities. `#page-analysis-section` exists in `cellucid/index.html` (1 occurrence). `window._comparisonModule` is real (`main.js:2359`). `analysis/ui/analysis-types/` exists (referenced from page 17). The `analysis/shared/debug-utils.js` debug toggles at `:228` (`?debug=1` or localStorage `debug=1`) match source (`debug-utils.js:43`, `:56`). No dead paths.

One cross-page note: `:199` points at "local-user switch self-test wiring in `ui/modules/dataset-controls.js`" — the *file* exists, but the self-test it refers to is the removed `runLocalUserInPlaceSwitchSelfTest` (see the page 14 finding and `tests/dataset-controls-exact-contract.test.mjs:536`). This sentence should go in the same change that fixes page 14.

---

## `12_figure_export_architecture.md`

**explains** — The Figure Export subsystem: its code map, three design goals, the three SVG export modes (full-vector / optimized-vector / hybrid) and their fidelity-vs-scalability tradeoff, why PNG export uses WebGL2 plus its `tEXt` provenance metadata, the shared legend/axes/orientation builders, why export UI is excluded from session serialization, and two troubleshooting entries.

**has** — `none`. Verified.

**verdict** — **DIAGRAM**

**needs** — One diagram: the export engine's control flow.

`figure-export-ui.js` (sidebar panel, mode selection, preview) → `figure-export-engine.js` (**snapshots current view state**: view id, camera, visibility mask, colors, dimension level — the five things listed on page 18 `:77–82`) → branch to two renderers. The SVG branch (`renderers/svg-renderer.js`) splits three ways into `full-vector` / `optimized-vector` / `hybrid`, annotated with what each emits (every point as a vector element / density-reduced vector / rasterized `<image>` + vector annotations). The PNG branch (`renderers/png-renderer.js`) shows the two-stage composition: WebGL2 point rasterization → Canvas2D composition → `tEXt` chunk injection. Shared builders (axes, legend, orientation indicator, centroid overlay) drawn as a bar feeding **both** branches — that shared-ness is the page's stated extensibility goal at `:43` and its correctness constraint at `:89`.

Derived from: `cellucid/assets/js/app/ui/modules/figure-export/index.js` (module wiring), `figure-export-engine.js`, and the module's own `figure-export/README.md` — which the page already designates as the starting point at `:26`, making it the natural source for the generator to reconcile against.

**Not a screenshot.** A picture of the export panel, or a sample exported figure, belongs in the user-facing figure-export chapter (`_static/screenshots/figure_export/` already exists as a topic directory). This page is about the code path.

**notes** — Page verifies clean. All eight code pointers resolve, including `figure-export/README.md`, both renderer files, and `state-serializer/README.md`. The `data-state-serializer-skip="true"` exclusion mechanism at `:102` is consistent with page 08 `:149–156` and page 10 `:334`. The three SVG mode ids (`full-vector`, `optimized-vector`, `hybrid`) are used consistently by page 18. No dead paths, no stale prose found.

---

## `13_debugging_playbook.md`

**explains** — The procedural "if it's broken, do this" playbook: record the environment, enable debug logging, then a symptom→subsystem router into ten diagnostic sections (boot, data loading, fields/legend, filtering, highlights, dimension, rendering/GPU, sessions, analysis, export, community annotation), each with ordered likely causes, a testable confirmation step, and a fix.

**has** — `none`. Verified.

**verdict** — **DEVTOOLS** — the highest-value devtools target in the section, and the page where figures most change what a non-expert can accomplish.

**needs** — Three captures. The page's whole method is "look at this and compare"; readers who have never opened devtools cannot execute it from prose.

*(1) Console, debug logging on.* After `localStorage.setItem('CELLUCID_DEBUG','true'); location.reload();` — show the prefixed startup log stream with `[Main]` and `[Viewer]` lines visible, enough rows to make the shape obvious (~15 lines). This is step 2 of the playbook (`:48–55`) and it is what tells a reader their toggle worked.

*(2) Network panel, two failure signatures side by side.* (a) A request blocked by CORS, with the red status and the CORS error text visible in the row/detail; (b) a request that returned `200` with an HTML body for a `.json` URL — show the Response preview tab with the `<!DOCTYPE html>` visible. These are the two causes behind the page's most-repeated symptom, "Unexpected token `<`" (`:134–140`, and again on pages 03 and 09). The HTML-body-with-200 case is genuinely hard to describe in words and trivial to recognize in a picture.

*(3) Console, the exact renderer error.* `[Viewer] updatePositions: position count mismatch` — the string the page tells readers to search for at `:221` and page 07 `:240`. One line, cropped tight.

**Leak risk — the highest in the section; handle deliberately.**
- The Console prints the dataset base URL. On a local-user dataset this can include a filesystem path with the operator's username. Use a **built-in demo dataset only**.
- The Network panel's request rows show the full origin. Capture from `http://localhost:8000`, and crop the browser URL bar out — it may carry `?remote=`, `?github=`, or `?exportsBaseUrl=`.
- Never capture with community annotation signed in: the UI displays the GitHub identity (`:313`) and Session Storage holds the OAuth token.
- Devtools request headers can contain `Cookie` / `Authorization`. Capture the **Response** headers pane, not Request, or redact.
- Use a clean browser profile (the page itself recommends this at `14:188`) so no unrelated extensions or personal tabs appear.

**notes** — Page verifies clean; both code pointers at `:332–333` resolve (`community-annotations/REPO_SETUP.md`, `community-annotations/scope-lock.js`), as do `utils/debug.js` (`:54`) and `analysis/shared/debug-utils.js` (`:64`). All five `{doc}` targets resolve including ``../n_benchmarking_performance/index``. The debug toggles documented at `:48–65` are the ones that actually exist in source — this page is *free* of the dead-flag problem that afflicts page 04, which is worth noting: the two pages document overlapping surfaces and only one of them is stale.

---

## `14_testing_ci_and_release_process.md`

**explains** — How to validate web app changes: the validation contract (`npm test`, `npm run test:browser`), running two browser suites concurrently via `CELLUCID_BROWSER_TEST_PORT`, a nine-part manual smoke checklist (A–I), automated checks (protected smoke-density contracts, CSS token validation, a dev-only dataset reload self-test), a docs build sanity check, and a six-item release checklist.

**has** — `none`. Verified.

**verdict** — **TERMINAL**

**needs** — Two terminal captures.

*(1) The contract suite, exact command:* `npm test` (run from `cellucid/`; resolves to `node scripts/run-web-tests.mjs`). Must be visible: the tail of the run showing the pass/fail counts and the total — a reader needs to know what "green" looks like before they can recognize red. Crop to the last ~20 lines.

*(2) The port-collision message, exact command:* run `npm run test:browser` twice concurrently, and capture the **second** run's output. Must be visible: the message naming the taken address `127.0.0.1:4173` and the `CELLUCID_BROWSER_TEST_PORT` variable, then the successful move via `CELLUCID_BROWSER_TEST_PORT=4183 npm run test:browser`. This is the page's most specific operational claim (`:26–40`) and it is verifiable — `tests/browser-test-port-contract.test.mjs` asserts `DEFAULT_BROWSER_TEST_PORT === 4173`, sample port `4174`, and that `describePortInUse` names `127.0.0.1:4173` and the variable. A capture here is directly backed by a contract test, so it will not silently rot.

**Leak risk — low but present.** Test output prints repo-absolute paths in stack traces on failure; capture a **passing** run, and set `PS1='$ '` so the prompt carries no username or hostname.

**notes** — Three real defects, one of them a documented-but-removed API.

- **`:154–167` — documents removed behaviour.** The whole "Dev-only dataset reload self-test (analysis cache reset)" section, including `window.__CELLUCID_DEV__.runLocalUserInPlaceSwitchSelfTest(...)` at `:160`. `tests/dataset-controls-exact-contract.test.mjs:534–541` is a test named *"dataset controls contain no dev global, reload, or guessed-source path"* whose body asserts `assert.doesNotMatch(moduleSource, /__CELLUCID_DEV__/)` and `assert.doesNotMatch(moduleSource, /runLocalUserInPlaceSwitchSelfTest/)`. The docs describe a helper the test suite exists to keep out. Delete the section; fix the matching sentence at `11:199`; fix `04:118`.
- **`:19–21` — the validation contract is incomplete.** It lists `npm test` and `npm run test:browser`. `cellucid/package.json` defines a third: `"test:worker-bundle": "npx --yes wrangler@4.115.0 deploy --dry-run --config wrangler.community-annotations.jsonc"`. The workspace `CLAUDE.md` lists all three as the authoritative local checks for this repo. A contributor following this page alone will not run the worker-bundle check — and `:193`'s release step ("confirm worker origin + CORS allowlist are correct") is precisely what that command validates.
- **`:176` — the suggested docs build does not reproduce CI.** The page gives `make -C cellucid-python/docs html`. That target is `$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html $(SPHINXOPTS)` with `SPHINXOPTS ?=` empty — **no `-W`**. CI runs `sphinx-build -W --keep-going -b html docs docs/_build/html` (`cellucid-python/.github/workflows/docs-check.yml:30`). A contributor runs the make target, sees it pass, and CI fails on a warning. Either document the CI command verbatim or pass `SPHINXOPTS="-W --keep-going"`.
- Everything else verifies: `cellucid/assets/css/README.md`, both `scripts/validate-token*.js`, and `dataset-controls.js` all exist; the `4173`/`4174` defaults and `test-results/<port>/` isolation match `tests/browser-test-port-contract.test.mjs`; the smoke-density contract claims at `:120–135` are consistent with the verified constants on page 07.

---

## `15_extension_points_overview.md`

**explains** — The supported extension points: five golden rules (crisp boundaries, no per-frame DOM work, no hot-path allocations, decide persistence on day one, document edge cases), five extension categories each routing to a dedicated page and listing its key files, and a four-branch "should this be state, UI, or viewer?" decision checklist.

**has** — `none`. Verified.

**verdict** — **NONE**

**needs** — Nothing. The page is a router: rules → five categories → three `{doc}` links. The only figure it could carry is a decision tree for "state vs UI vs viewer", but that checklist is already four if-then lines (`:115–118`) and a tree drawing would contain the identical four branches with more visual overhead and no added dimension — there is no ordering, no cycle, and no fan-in to reveal. The file lists per category (`:52–56`, `:68–72`, `:84–86`, `:93–98`, `:103`) are better served by the generated module graphs on pages 08, 11, and 12, which this page should link rather than re-illustrate.

**notes** — Page verifies clean. All 13 code pointers resolve, including the four glob forms (`ui/modules/*`, `analysis/plots/types/*`, `analysis/compute/*`, `figure-export/renderers/*`) and `data/data-source.js`. The `session/contributors/` directory referenced at `:103` exists and matches page 10 `:43`. No dead paths, no stale prose found.

---

## `16_extension_point_add_ui_module.md`

**explains** — A seven-step mechanical recipe for adding a sidebar UI module: decide state-vs-UI-vs-viewer, add the accordion DOM scaffold in `index.html` (with an early session-persistence decision), register selectors in `dom-cache.js`, create the module file with the established initializer shape, wire it into the UI coordinator, style it with the design system, decide session persistence, and add troubleshooting hooks — plus two common-mistake entries.

**has** — `none`. Verified.

**verdict** — **DIAGRAM**

**needs** — One small diagram: the **wiring chain**, showing where a contributor's code goes and in what order. Five nodes, left to right, each labelled with the file and the specific thing added:

`index.html` (`<details class="accordion-section" id="my-feature-section">`, stable id) → `ui/core/dom-cache.js` (new entry in a `dom.*` domain bucket) → `ui/modules/my-feature-controls.js` (`export function initMyFeatureControls({ state, viewer, dom, callbacks })`) → `ui/core/ui-coordinator.js` (import + call inside `initUI`) → `DataState` (method calls out, events back in — draw the return edge dashed to match the event convention used on pages 05, 06, and 11).

Add one branch off node 1: the session-persistence fork (`data-state-serializer-skip="true"` → excluded, or default → captured by DOM id). The page raises this decision twice (`:63–72` and `:142–154`) precisely because contributors miss it, and a visible fork is a better reminder than two separate prose passages.

Derived from: an existing simple module traced end to end — `render-controls.js` is the good exemplar (it appears in `dom-cache.js`, is imported at `ui-coordinator.js:16`, and owns a discrete accordion). A generator could emit this chain for any named module, which also makes the figure a template rather than a one-off drawing.

This is the lowest-priority DIAGRAM in the section, and the honest alternative verdict is NONE — the seven steps are already ordered and each names its file. It earns the figure because a first-time contributor's actual question is "how many files do I have to touch, and in what order", and the answer (five, in a fixed order, with a fork) is a shape, not a list.

**notes** — Page verifies clean. All code pointers resolve: `assets/css/README.md`, `types/design-tokens.d.ts` (`:137` — verified present), `state-serializer/README.md`, `dom-cache.js`, `ui-coordinator.js`, `session/contributors/`, `assets/js/rendering/`. `cellucid/assets/js/app/ui/modules/my-feature-controls.js` (`:91`) is explicitly a hypothetical example, not a dead path. `input/select` (`:67`) is HTML element names, not a path. The `NotificationCenter` reference at `:168` corresponds to real `app/notification-center.js` / `app/notification-center/`. No dead paths, no stale prose found.

---

## `17_extension_point_add_analysis_mode.md`

**explains** — Three distinct analysis extensions and how to pick: Option A adds a UI mode (create a factory under `ui/analysis-types/`, register it in `_registerAnalysisTypes()`, ensure the `data-mode` container exists, route data through `DataLayer`, implement cleanup); Option B adds a plot type; Option C adds a compute operation, transform, or stat test — plus a validation list and two troubleshooting entries.

**has** — `none`. Verified.

**verdict** — **NONE**

**needs** — Nothing. Every figure this page could want is the analysis subsystem graph, which belongs on — and is specified for — `11_analysis_architecture.md`; page 17 is its procedural counterpart and should link that figure rather than carry a second copy. The three options are numbered file-editing recipes whose content is symbol names (`_registerAnalysisTypes()`, `onPageSelectionChange(pageIds)`, `destroy()`, `purgePlot(...)`, `data-mode="..."`) — text that must be typed exactly, which is the one form of content a figure actively degrades.

**notes** — Page verifies clean. All nine code pointers resolve, including `analysis/ui/analysis-types/` (`:25`), `analysis/compute/operation-handlers.js` (`:108`), `analysis/stats/statistical-tests.js` (`:115`), `analysis/data/transform-pipeline.js` (`:118`), and the glob `analysis/plots/types/*.js` (`:87`). The mode ids named at `:45` (`simple`, `detailed`, `correlation`, `differential`) are plausible registration keys and are consistent with page 11 `:97`'s list of registered analysis types. No dead paths.

---

## `18_extension_point_add_export_renderer.md`

**explains** — How to extend Figure Export: decide what you are adding (SVG mode / PNG option / new format) under three design constraints, locate the entry points and renderer/utils directories, implement the renderer (SVG reuse rules; PNG high-DPI and `tEXt` metadata constraints), add UI wiring, confirm WYSIWYG against five view properties, add performance safeguards, plus a testing checklist and two troubleshooting entries.

**has** — `none`. Verified.

**verdict** — **NONE**

**needs** — Nothing. The one useful figure — the export engine control flow with its SVG/PNG branch and shared builders — is specified for `12_figure_export_architecture.md`, which this page already requires as prerequisite reading at `:5–6`. Duplicating it here would create the same two-maps-of-one-tree drift the section already suffers from between page 01 and `app/README.md`. The remainder is constraints and checklists.

A tempting alternative — a side-by-side comparison image of `full-vector` vs `optimized-vector` vs `hybrid` output at the same zoom — is a genuinely good figure, but it is a *user-facing* one about choosing a mode, not a developer-facing one about adding a renderer. It belongs in the figure-export user chapter (`_static/screenshots/figure_export/` already exists), not here.

**notes** — Page verifies clean. All six code pointers resolve, including `figure-export/utils/` (`:35`) and `figure-export/README.md` (`:133`). The five WYSIWYG properties at `:77–82` match the engine responsibilities described on page 12 `:31`. The three SVG mode ids match page 12 `:47–60` exactly. No dead paths, no stale prose found.

---

## Cross-cutting observations

**1. Figure paths, for whoever produces these.** From `docs/user_guide/web_app/p_developer_docs/`, the correct relative reference is `../../../_static/screenshots/<topic>/<name>.png` (three levels: `p_developer_docs` → `web_app` → `user_guide` → `docs`), matching the convention in the rest of the docs. There is currently **no** `_static/screenshots/developer/` topic directory — the 13 existing topics are `analysis`, `benchmarking_performance`, `community_annotation`, `data_loading`, `figure_export`, `filtering`, `highlighting_selection`, `jupyter`, `jupyter_hooks`, `server`, `sessions_sharing`, `vector_field_velocity`, `web_app`. Terminal and devtools captures are PNGs like any other and would live under a new topic (`developer/`); mermaid diagrams, if the extension is added, live inline in the `.md` and need no `_static` entry at all.

**2. The `-W` build applies to figures too.** Every `{figure}` needs `:alt:` — a missing alt is not currently a Sphinx error, but a broken image path is, and `sphinx-build -W --keep-going` turns it into a CI failure. Verify each new path builds before committing.

**3. The defect pattern is one-directional.** Nine of the ten stale claims in this section are on pages 01, 04, 05, 08, and 14 — the pages that *enumerate* things (files, keys, modules, commands). The pages that *explain mechanisms* (06, 07, 09, 10, 11, 12, 13) verify clean, including page 07's exact GPU constants and page 10's byte-level framing. Enumerations rot; explanations do not. That is the strongest argument in this audit for generating the enumerations (page 01's tree, page 08's module index) rather than adding figures to them, and it is why page 04's verdict is STALE-RISK rather than DEVTOOLS.

**4. Two contract tests already encode doc truths.** `tests/dataset-controls-exact-contract.test.mjs:534–541` and `tests/community-annotation-wire-contract.test.mjs:3293` each assert the absence of something these docs still document. The test suite is ahead of the documentation. Any generator added for the diagrams should follow the same pattern — assert, in `cellucid/tests/`, that the committed figure source matches a fresh walk of the real tree — so that the next divergence fails a test instead of surviving to the next audit.
