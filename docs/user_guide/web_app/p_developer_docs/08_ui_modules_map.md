# UI modules map

This page is a **map of the UI layer**: which modules exist, what DOM they own, and how they communicate with state/rendering without creating tight coupling.

It is especially useful when:
- you want to add a new sidebar module,
- you want to find where a UI behavior lives,
- you are debugging “why did this update twice?” issues.

## At a glance

**Audience**
- Computational users: skim “UI coordinator” and “Which module owns what?”
- Developers: read fully; use the module table as a navigation map.

**Time**
- 20–40 minutes

**Prerequisites**
- {doc}`05_app_architecture_overview`

---

## UI coordinator (the orchestrator)

The UI layer is orchestrated by:
- `cellucid/assets/js/app/ui/core/ui-coordinator.js` (`initUI`)

`initUI` is intentionally thin:
- collects DOM references once (`collectDOMReferences`)
- initializes feature modules (`initX` functions)
- wires cross-module callbacks (e.g. “visibility changed → refresh legend + badges”)
- subscribes to a small set of `DataState` events

This separation exists to prevent:
- a single monolithic UI file that becomes unmaintainable
- accidental performance regressions from repeated DOM queries

---

## DOM cache (single source of truth for selectors)

Instead of every module calling `document.getElementById(...)` repeatedly, Cellucid uses:
- `cellucid/assets/js/app/ui/core/dom-cache.js`

`collectDOMReferences(document)` returns a nested object keyed by UI domains:
- `dom.dataset`
- `dom.fieldSelector`
- `dom.filter`
- `dom.highlight`
- `dom.render`
- `dom.camera`
- `dom.view`
- `dom.session`
- `dom.communityAnnotation`
- …

Design rule:
- When adding new UI elements, prefer adding them to the DOM cache rather than scattering selectors across modules.

---

## Module index (what exists today)

Most UI modules live under:
- `cellucid/assets/js/app/ui/modules/`

Most — not all — are imported and initialized by:
- `cellucid/assets/js/app/ui/core/ui-coordinator.js`

:::{important}
“Lives under `ui/modules/`” and “is initialized by the coordinator” are two
different facts, and four of the modules below are initialized by a *sibling
module* rather than by the coordinator:

- `dataset-connections.js` ← initialized by `dataset-controls.js`
- `dataset-drop.js` ← initialized by `dataset-connections.js`
- `field-selector-gene-expression.js` ← initialized by `field-selector.js`
- `field-selector-deleted-fields.js` ← initialized by `field-selector.js`

Two more are not initialized by the coordinator at all:
`sidebar-controls.js` and `select-bar-placement.js` are wired by
`ui/core/sidebar-layout.js`, which `main.js` calls right after the viewer is
created. The sidebar toggle and, on a phone, the panel's peek bar
(`#phone-sheet-bar`: the dataset's name and a handle; a tap or an upward swipe
opens the panel) are painted with the first frame, so they work before the
first dataset has loaded, which is when the coordinator runs. A phone starts
with the bar: `ui/core/viewport-init.js` marks the root before the stylesheet
applies (`data-phone-sheet="peek"`), unless this tab last left the panel open
(`cellucid:phone-sheet` in session storage, `ui/core/phone-sheet-memory.js`).

If you are tracing why a control is inert, check which of these paths owns it
before assuming the coordinator does.

The `ui/` layer also has directories this page does not otherwise cover:
`ui/components/` (shared widgets — the coordinator initializes
`components/info-popovers.js` directly; `components/menu.js` is every menu's
look and keys: `openMenu` builds a menu for the items of the moment, at a
button or a point, and `menuKeyDown` is the keys the markup menus below answer
too; `components/app-dialog.js` is every modal dialog: `wireAppDialog` opens a
native `<dialog>` with `showModal()` and adds the close button, the backdrop
click, focus back to the opener, keys kept from the page and the notification
cards kept above it (`keepAboveModals`), and `components/confirm-dialog.js`
is the one confirmation built on it), `ui/category-builder/`,
`ui/onboarding/` (the welcome dialog and the global keyboard shortcuts), and
`ui/keyboard-move.js`.
:::

### Core “shell” modules

| Module | File | Owns |
|---|---|---|
| Sidebar toggle/resize | `ui/modules/sidebar-controls.js` | sidebar visibility (`inert` while hidden) + `--sidebar-width`; publishes the sidebar's left occlusion of the canvas, from which the viewer resolves the presented layout |
| Selection bar placement | `ui/modules/select-bar-placement.js` | centers `#select-bar` in the viewer's presented area (`viewer.onPresentedLayoutChange`) |
| Stats display | `ui/modules/stats-display.js` | top stats line (dataset/field/cell counts) |
| Cell info | `ui/modules/cell-info.js` | the card describing the cell under the pointer (**Show cell info**) |
| Cite dialog | `ui/modules/cite-dialog.js` | the footer's Cite Cellucid dialog, read from the deployed `CITATION.cff` |
| License dialog | `ui/modules/license-dialog.js` | the footer's license link and the BSD 3-Clause text it opens |
| Place names | `ui/place-names.js` | the names of the panels and blocks a message sends someone to (`PLACES`, `placePath`), as their headings write them; `tests/ui-place-names-contract.test.mjs` holds each to its heading in `index.html` and fails on a message that names a place in words of its own |

### Dataset + session modules

| Module | File | Owns |
|---|---|---|
| Dataset controls | `ui/modules/dataset-controls.js` | dataset dropdown (the samples, or a connected catalog's datasets as "Dataset from <host>"; `__dataset_unlisted__` while a file, single-dataset link or server dataset is open) |
| Dataset box | `ui/modules/dataset-card.js` | the **Dataset** box (`#dataset-info`, `data-state` empty/loading/ready): name, cells · genes, format · host, description, bytes read in place, **Details**; the retry notification |
| Open… menu | `ui/modules/open-data-menu.js` | a menu written in the markup bound to its button (the **Open…** button's **File…** / **Folder…**, My datasets, the History panel's ⋯), with the keys of every menu (`menuKeyDown`, `ui/components/menu.js`) |
| Address routing | `ui/modules/address-routing.js` | what the **Link or server** bar opens: a `cellucid serve` server, linked data, or (asked over `/_cellucid/health`) either |
| Dataset connections | `ui/modules/dataset-connections.js` | local file and folder pickers, the **Link or server** bar (server and link connectors), the **Server** / **Link** connection lines with **Disconnect**, the sign-in step |
| Dataset drop | `ui/modules/dataset-drop.js` | a file or folder dropped anywhere on the page, opened as the matching picker would (initialized by `dataset-connections.js`; `ui/core/file-drop-guard.js` stops the browser navigating to a dropped file; a folder listed concurrently by `data/local-folder.js`; in Chromium the item's handle kept to reopen, `app/library/local-file-handles.js`) |
| Share dialog | `ui/modules/share-dialog.js` | the Dataset box's **Share…**: link, badge and embed snippet for a dataset a URL can reopen |
| Session controls | `ui/modules/session-controls.js` | Loading a `.cellucid-session` file (format 4: the history; State › Session › **Load…** › **Open a session file…**) and opening a saved state of this browser: the session's own dataset and cells are opened first, then its history is loaded. Writing one is the history's action `history.share` (State › Session › **Share…**, `ui/modules/state/share-dialog.js`) or `history.save` (the Load dialog) |
| My datasets | `ui/modules/my-datasets/` (`index.js` `initMyDatasets`, `menu.js`, `dialog.js`, `model.js`) | Session › Open data › **My datasets**: the five most recent datasets opened in this browser (`RECENT_DATASETS`), the one open now marked, then **All…**, a dialog of every dataset (Open, Continue its latest saved state, Forget); drawn from `app/library/my-datasets.js`, `app/library/local-file-handles.js` and the saved states, keeping no list of its own; opening runs inside the person's click. Not a session control and not a history step |
| Session sections | `ui/modules/session-sections.js` | the Session panel's collapsible sections, **Subsampling** (`#reading-section`) and **Data mapping** (`#data-mapping-controls`): collapsed by default, each a one-line summary in its header while collapsed (`setSummary`; `subsamplingSummary`: "100% cells · genes 6%"), the open state kept per viewer in `localStorage` (`cellucid_session_sections`, every access guarded; never a session or the address), and `expandSessionSection`, which whatever leads into one calls first |
| Cells to open | `ui/modules/subsample-control.js` | Session › Subsampling › **Cells to open**: the share slider, method and seed, **Apply**, one status line and its **Details** table |
| Genes read for | `ui/modules/gene-read-control.js` | Session › Subsampling › **Genes read for**: the slider of the share of the open cells a gene is read for (steps: the share now, then 1, 5, 10, 25, 50, 100% up to what the gene index fits, `geneReadSteps`), its one status line (`geneReadStatus`: the chosen cells and their cost, repriced every 5 s while shown; while reading, the cells read and the time left), **Read** / **Stop** (`readGenesForMoreCellsSaid`). The one place the gene sample changes; `createGeneReadNote` / `createGeneReadLink` are the lines in a gene's legend and an analysis's result that link to it (`revealGeneReadControl`) |
| Stored line | `ui/modules/kept-stores-line.js` | Session › **Stored**: what this browser keeps, one line beside **Manage…** (`summarizeKeptStores`: `{text, title}`, the title naming every kind) |
| Data mapping controls | `ui/modules/data-mapping-controls.js` (+ `ui/modules/data-mapping/`) | Session › Data mapping (its summary through `onSummary`): the Embeddings checklist (on/off, defaults, column selections; arrays rarely drawn under **Other arrays (N)**), the Genes (**Modality** for a dataset of several, Expression, Gene names, Values, Highly variable genes), Neighbor graph, Latent space and Vector fields sections, Revert/Apply |

### State panel

`ui/modules/state/` (its `README.md` is the full map). Header
`History (i) 🔍 ↶ ↷ ⚙ ⋯`, a status line, a virtualized branch graph; the
session's **Share…** and **Load…** sit above it.

| Module | File | Owns |
|---|---|---|
| State section | `index.js` (`initStateSection`) | The history of the dataset shown, drawn from the page's current history (`app/history/current.js`): the header (search, undo, redo, History settings, the menu, whose **Save AnnData with this state…** dispatches `history.anndata.save`, which tells its own progress, outcome and refusals), the status line (recording switch, branch picker, steps, stale), the graph's row operations, each a `history.*` action dispatched as the person's; one operation at a time; while recording is paused, Keep / Discard / Cancel before the history moves |
| The seam to the history | `history-source.js` | `trackerView` (the tracker in use as accessors and derived facts), `createHistoryActions` (the only way the panel acts), the plans a confirmation is worded from |
| Branch graph | `graph-layout.js` (pure), `graph.js` | Rows, lanes and segments (path, fork stubs, expanded branches, turns, groups, folded runs, the paused marker); only the rows in view are in the DOM. Click selects, Shift+click extends, double-click (or a double tap) shows, `›` opens details, ▸ expands |
| Details, notes, search | `step-details.js`, `notes.js`, `search.js` | A step's, branch's, turn's, group's or run's details and actions (**Run again** only on a stale step of an action marked `rerun` in the catalog); notes (a Markdown subset from text nodes); the search box and Filter menu |
| Branches | `branch-picker.js`, `branch-map-dialog.js`, `compare.js` | The branch picker (trackers, branches, archived, New branch); the branch map (every branch as a lane graph and a table; archive or delete several); two branches side by side through the views' own actions (`view.keep` with each branch, `view.layout.set`) |
| Share and Load | `share-dialog.js`, `saved-states-dialog.js` | Share this history (History / Current state only / Steps only, and what goes in; `history.share`); Load (saved states of this browser: open, replay here, rename, delete; save as history, state or steps, `history.save`; open or replay a session file) |
| Recently deleted, banners, replays | `recently-deleted.js`, `workflow-banner.js`, `replay-report-dialog.js` | Every undoable edit of the history with Restore and Empty; the line above the graph for a loaded, shared-state or replayed tracker; a replay's report |
| Compare two steps | `diff-dialog.js` | What differs between two steps (`history.diff`) |
| Tour, copy, report | `tour.js`, `share.js`, `dataset-facts.js` | The guided tour (each step shown by `history.jump`); copy pseudocode or methods (`history.pseudocode`, `history.methods`); the printable report, pictures drawn on demand (none stored) |
| Undo keys, words | `keyboard.js`, `model.js`, `shared.js` | Ctrl/⌘+Z and Ctrl/⌘+Shift+Z (left to text fields, dialogs, menus and a selection in progress); the menus at a row or the pointer are the shared `openMenu` (`ui/components/menu.js`); the words of rows and the status line; the list's keys are the History list group of the one shortcut table (`ui/onboarding/shortcut-table.js`); every question in the app's one confirmation dialog, the edit toast with Undo |

The panel keeps no history of its own: it reads the current history and
redraws when it tells of a change, at most once a frame and only while it can
be seen; it calls no history mutator, and its controls are not session
controls (`data-state-serializer-skip` on the section). See
{doc}`29_action_layer_domains_and_history`.

### Agents

The agent dock at the foot of the plot (the activity lines, the confirm
cards, the chat bar) and the Agents panel are drawn by `app/agent/`, not by
`ui/modules/`: `app/agent/index.js` (`initAgents`, from the coordinator),
`app/agent/presence/` (cards, lines, the permission settings dialog),
`app/agent/chat/` (the chat bar), `app/agent/connect/` (Agents › Claude and
ChatGPT: the connector's status line, address, approval dialog and
Manage…). `ui/modules/dataset-open-action.js` provides `datasets.list` and
`dataset.open` for agents and hosts over the Open data controls. See
{doc}`32_agents_in_the_tab`.

### Actions installed by the page

These modules own no DOM of their own; they provide the actions of the
controls they stand behind and record what a person changes by hand:

| Module | File | Provides |
|---|---|---|
| Settings actions | `ui/modules/settings-actions.js` (`installSettingsActions`, from `main.js`) | `setting.set` and `settings.list` over every control of the inventory but the layout's; a person's change of a settings control recorded by observation |
| Layout actions | `ui/modules/layout-actions.js` (`installLayoutActions`, from `main.js`) | The `layout.*` actions (sections, panels, sidebar, theme) and their observed changes |
| Field actions | `app/state/actions/field-actions.js` (`installFieldActions`, from the coordinator) | The state's own actions: fields, categories, colors, filters, cells hidden by hand, gene sets, the data mapping |

### Field/filter/highlight modules

| Module | File | Owns |
|---|---|---|
| Field selector | `ui/modules/field-selector.js` | categorical/continuous/gene selectors + copy/rename/delete |
| Gene expression helper | `ui/modules/field-selector-gene-expression.js` | gene dropdown/search behavior |
| Deleted fields panel | `ui/modules/field-selector-deleted-fields.js` | restore/purge UI |
| Filter controls | `ui/modules/filter-controls.js` | filter stack UI + count |
| Legend renderer | `ui/modules/legend-renderer.js` | legend DOM + outlier UI container; a continuous field's **Gray for zero** (genes with no negative value), **Log color scale** and **Rescale to range** switches, a gene sample's one-line layer note and its link to Genes read for (`createGeneReadLink`), and **Visible range** (label, slider and number box a row each, `legend/range-bounds.js` `boundText`) |
| Outlier filter control | `ui/modules/outlier-filter-control.js` | the legend's **Outlier filter**: **Choose latent space** without a latent space, else the slider with at most one line, progress and **Cancel** |
| Highlight controls | `ui/modules/highlight-controls.js` | highlight pages + groups + selection mode tools |

### View/render/camera modules

| Module | File | Owns |
|---|---|---|
| View controls | `ui/modules/view-controls.js` | multiview (“Keep view”), layout mode, **Link cameras**, one row per view (what it shows, field, navigation mode, linked), the Navigation “Applies to” row |
| Representation controls | `ui/modules/representation-controls.js` | the `View N shows` picker (`#view-representation-select`) + keys 1/2/3 (default, then cycle) |
| Render controls | `ui/modules/render-controls.js` | points vs smoke, smoke sliders (scales in `smoke-control-scales.js`; the unavailable-device message in `smoke-availability.js`), connectivity toggles |
| Velocity overlay controls | `ui/modules/velocity-overlay-controls.js` | vector field overlay settings + field selection |
| Camera controls | `ui/modules/camera-controls.js` | the Navigation accordion: mode select of the focused view, pointer capture, projectiles |
| Navigation settings | `ui/modules/navigation-settings.js` | per-mode settings boxes, speed sliders and option checkboxes |

### “Large feature” modules

These are bigger subsystems but still treated as modules:

| Module | File | Owns |
|---|---|---|
| Figure export | `ui/modules/figure-export/` | export UI + engine + SVG/PNG renderers |
| Community annotation | `ui/modules/community-annotation-controls.js` | auth + pull/publish + voting/moderation UI |
| Community annotation (split parts) | `ui/modules/community-annotation/` | the controls module's own submodules (browser integration, connection flow, consensus column, modal shell, auto-pull storage) |
| Community annotation dialogs | `ui/modules/community-annotation-voting-modal.js`, `ui/modules/community-annotation-modal-owner.js` | the voting dialog and its ownership/lifecycle |
| Cinematic camera | `ui/modules/cinematic-camera/` | camera-path authoring, transport bar, playback, and the `camera-path.*` actions. Initialized by the coordinator; its keyframes and playback settings are the `camera-path` state domain, and a loaded session plays the path when it was saved with Autoplay on (`startAutoplay`) |
| Benchmark | `ui/modules/benchmark/` | synthetic data generation contract + worker, the configuration-matrix harness, GL upload counters, frame recorder, camera path. The panel’s wiring lives in `main.js`; the measurement code is here and in `assets/js/dev/benchmark.js` |
| Highlight subsystem | `ui/modules/highlight/` | the four selection tools, their copy, and the highlight-page UI |
| Legend subsystem | `ui/modules/legend/` | categorical and continuous legend rendering |
| Visualization reset | `ui/modules/visualization-reset.js` | reset actions (clear filters/highlights/etc.), and the global `R` reset-camera key |
| Field interaction owner | `ui/modules/field-interaction-owner.js` | arbitrates which module owns a field interaction |
| Camera input contract | `ui/modules/camera-input-contract.js` | the validated shape of camera control input |

:::{note}
`ui/modules/cross-highlight/` exists but contains only a `.gitkeep`. There is no
cross-highlight UI module today — cross-view highlight synchronization is
handled in the state layer (`app/state/managers/highlight-manager.js`) and in
`ui/modules/highlight/`. A reader who trusts a complete module map would go
looking for a module that is not there.
:::

:::{note}
Analysis (plots / DE / correlation / etc.) is initialized directly in `cellucid/assets/js/app/main.js` (not via `ui-coordinator.js`).
Its UI lives in `cellucid/assets/js/app/analysis/ui/`.
See {doc}`11_analysis_architecture`.
:::

---

## How modules communicate (actions + state events + callbacks)

Most modules follow this pattern:

- Inputs (DOM events) dispatch the action the control stands for
  (`dispatchFromControl`, `app/actions/control-dispatch.js`); the action's
  implementation calls the `DataState` (or viewer) methods. A control calls no
  mutating state method itself, so every change is one history step, undone the
  same way whoever made it ({doc}`29_action_layer_domains_and_history`).
  Settings controls and camera gestures change their owner first and are
  recorded after (`recordObserved`).
- `DataState` emits events (`field:changed`, `visibility:changed`, etc.).
- The UI coordinator reacts to those events and triggers targeted UI updates (legend/stats/badges).

Why this matters:
- It prevents “UI module A reaches into UI module B’s internal state”.
- It makes it easier to reason about what needs to update when.

Example: visibility changes
- Filtering module dispatches a `filter.*` action → its handler updates the
  filter → state recomputes visibility → state emits `visibility:changed`.
- UI coordinator debounces a visibility UI update and triggers:
  - filterControls render
  - highlight summary render
  - view badge render

---

## The control inventory (which UI state is a setting)

Every `input` and `select` with a stable DOM id and every accordion section of
the sidebar and its floating panels is in the **control inventory**
(`app/state-serializer/ui-controls.js`, read from the markup), and
`controlDomainOf` is the one map that says which state domain holds each:
the accordion sections and the theme are the `layout` domain's, the centroid
toggles the `views` domain's, every other control is a setting (the `settings`
domain), set by `setting.set` and recorded by observation. A control added to
the markup is a setting the history records, a session saves and an agent can
list and set, with no code. A pressed-button group
(`data-state-serializer-pressed-group`, the selection tool) is not in it: it
is page chrome, kept by a save and the autosave as the session's `chrome` and
put back after a load, never a step (`capturePageChrome` / `applyPageChrome`).

**Some UI subtrees are intentionally excluded** with
`data-state-serializer-skip="true"`: what they show is held by another domain,
or is not state at all.

A second attribute, `data-state-serializer-after-embedding="true"`, marks a
subtree whose controls stay in the inventory but are restored only after the
embeddings, kept views and active view are back, because their choices follow
the embedding the active view shows (the vector-field block: a field is bound
to one embedding). See `cellucid/assets/js/app/state-serializer/README.md`.

Examples:
- The Figure export section is excluded from the inventory; its settings are
  one record of the `settings` domain (`figureExport`,
  `exportSessionSettings()`/`importSessionSettings()` of the panel).
- The docked analysis panels and the floating analysis windows are excluded;
  they are the `analysis` domain.
- The Camera path is excluded; it is the `camera-path` domain.
- The State section and My datasets are excluded (they show the history and
  the library, which are not state of the dataset).
- `#hp-antialias` is never captured (a device capability); the theme is.
- Benchmark section is excluded.
- Dataset connection controls are excluded (sessions assume the dataset is already loaded).
- Community annotation UI state is excluded (auth/network-driven).

Code pointers:
- `cellucid/assets/js/app/state-serializer/ui-controls.js`
- `cellucid/assets/js/app/state-serializer/README.md`
- `cellucid/assets/js/app/state/domains/settings.js`, `layout.js`

---

## Troubleshooting (UI module boundary issues)

### Symptom: “One click triggers multiple UI updates”

Likely causes:
- Two modules subscribed to the same `DataState` event and both re-render the same DOM region.
- A module is initialized twice (usually due to calling `initUI` twice on reload paths).

How to confirm:
- Set a breakpoint in the module init function and reload.
- Log subscriptions by temporarily instrumenting `EventEmitter.on`.

Fix:
- Keep DOM ownership crisp: one module “owns” a DOM subtree.
- Prefer coordinator wiring over cross-module subscriptions.

---

Next: {doc}`09_data_loading_pipeline_and_caching` (data sources) or {doc}`15_extension_points_overview` (how to add modules safely).
