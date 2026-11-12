# Extension point: add a UI module

This guide is a step-by-step recipe for adding a **new sidebar UI module** to the Cellucid web app.

It is intentionally “mechanical”: follow the steps and you’ll end up with a module that:
- has clear DOM ownership,
- changes state only through actions, so every change is one undoable history
  step that an agent and a Python or R host can make too,
- does not regress performance,
- is held by a state domain (so steps, undo and sessions keep it), or is
  explicitly excluded.

## When to use this extension point

Use this when you want to add:
- a new accordion section (sidebar panel),
- new controls inside an existing section,
- a new UI workflow that primarily manipulates `DataState` and/or viewer settings.

If you are adding:
- a new analysis mode → use {doc}`17_extension_point_add_analysis_mode`
- a new export renderer → use {doc}`18_extension_point_add_export_renderer`

---

## Step 0: Decide “state vs UI vs viewer”

Before writing code, answer:

1) Does the feature affect **which points are visible** or **how they’re colored**?
   - If yes: it likely needs `DataState` changes (typed arrays + events).

2) Does the feature affect **how pixels are drawn** (shaders, buffers, per-frame render)?
   - If yes: it likely belongs in `cellucid/assets/js/rendering/`.

3) Is it a **pure UI preference** (toggle, layout, display option)?
   - It can be a UI module over existing state/viewer APIs; a stable sidebar
     control is a setting with no code (Step 1, "Settings and state").

4) Which **action** does each control stand for, and which **state domain**
   holds what it changes? An existing action (`app/actions/definitions/`) is
   reused; a new one is defined first
   ({doc}`29_action_layer_domains_and_history`, “How to add an action”).

If you are unsure, read:
- {doc}`05_app_architecture_overview`
- {doc}`06_state_datastate_and_events`

---

## Step 1: Add the DOM scaffold in `cellucid/index.html`

Most sidebar panels are `<details class="accordion-section">` elements.

Add a new section with:
- a stable `id` (used by DOM cache + serializers + tests)
- a `<summary>` label
- a `.accordion-content` container where your module will render or attach controls

Example pattern (illustrative):

```html
<details class="accordion-section" id="my-feature-section">
  <summary>My Feature</summary>
  <div class="accordion-content">
    <!-- controls go here -->
  </div>
</details>
```

Keep the `<summary>` to plain text. A summary is itself the button that opens
and closes the panel, so a button, link, input or info `i` inside it is a
control nested in a control: screen readers announce it poorly, and its clicks
and keys reach the toggle as well. The browser suite
(`tests/browser/panel-header-controls.spec.mjs`) fails if any summary in the app
contains an interactive element.

At startup `cellucid/assets/js/app/dockable-accordions.js` stands every sidebar
section in a frame, which is also what docks and floats:

```html
<div class="accordion-panel">
  <div class="accordion-panel-actions"><!-- undock / dock buttons --></div>
  <details class="accordion-section" id="my-feature-section">…</details>
</div>
```

The action bar is laid over the right end of the summary row, and the summary
keeps an empty slot of the same size so its title never runs under the buttons.
The section keeps its `id` and its open state; whether it floats, its position
and its stacking layer belong to the frame. If a panel needs header controls of
its own, pass them to `register(details, { actions: [...] })` rather than
appending them to the summary — that is how an analysis window gets its copy
and close buttons.

### Readiness decision (also early)

Ask when your control’s listener is attached.

If it is wired anywhere **after an `await` in `main.js`’s bootstrap** — which is
true of anything that needs a loaded dataset, the data-source manager or the
renderer — the control must be authored `disabled` in `index.html` and
registered with `cellucid/assets/js/app/ui/core/deferred-control-readiness.js`,
which enables it once the listener exists.

This is not cosmetic. A control that is clickable before anything listens
swallows the click silently, and the user cannot distinguish that from a broken
feature. Add `aria-describedby="data-source-status"` too, so a screen reader
hears why it is inert.

### Settings and state decision (early)

Decide what holds the value of each control:

- **A setting, with no code.** Every `input` and `select` with a stable DOM id
  outside a skipped subtree is in the control inventory
  (`app/state-serializer/ui-controls.js`, read from the markup), and
  `controlDomainOf` makes it a setting (the `settings` domain) unless it is an
  accordion section or the theme (`layout`) or a centroid toggle (`views`).
  A setting is then, with no code of yours: recorded as a `setting.set` step
  when a person changes it (`ui/modules/settings-actions.js` observes its
  `input` or `change`; a drag is one step), undone and
  redone, saved in the history and in a session file, listed by
  `settings.list` and set by `setting.set` for an agent or a host (its label,
  section, kind, range or choices read from the markup). Its owner's listener
  must already be attached when a restore runs: a restore writes the value and
  dispatches the event the owner listens for, and a control whose owner refuses
  the value fails the restore.
- **Another domain's state.** A control whose value is held elsewhere (a
  field picker, a filter, a highlight) dispatches its action, and the domain's
  owner captures the value from the state; add its id to the inventory's
  `DOMAIN_OWNED_IDS` (and to the list in `app/state-serializer/README.md`,
  which `tests/session-docs-current-contract.test.mjs` holds to it) or mark its
  subtree skipped, so it is not also a setting.
- **Page chrome.** A toolbelt of buttons of which exactly one is pressed and
  which decides how later gestures are read (the selection tool) is marked
  `data-state-serializer-pressed-group`: kept by a save and the autosave and
  put back after a load (a session's `chrome`), never a setting and never a
  history step (`capturePageChrome` / `applyPageChrome`).
- **Not state.** Transient state (previews, auth tokens, network state, a
  dialog) is excluded: add `data-state-serializer-skip="true"` on the section
  or a subtree.

### Help decision (also early)

A control group explains itself with an info button beside its title, not with
help text in the panel: a `button.info-btn` whose `aria-controls` names an
`info-tooltip` element carrying `data-info-tip="<key>"` (see any
`control-group-header` in `index.html`). Its text is one entry of
`cellucid/assets/js/app/ui/components/info-tips.js`: a title, a short body and
at least one link into these docs by page and heading anchor.
`tests/info-tips-contract.test.mjs` checks every entry and, when this docs tree
is present, that every linked heading exists, so renaming a docs heading a tip
links to fails it.

### Accessibility decision (also early, and cheap here)

The browser suite asserts that **every control in the analysis panels has an
accessible name** and that icon-only buttons resolve to exactly one name. Give
your controls real labels now; retrofitting them means fighting a test that is
already green elsewhere.

If you add a persistent store or a network call, see the “Definition of done”
list in {doc}`00_how_to_use_this_section` — both have documentation obligations
outside this section.

Reference:
- `cellucid/assets/js/app/state-serializer/README.md` (the control inventory)
- {doc}`29_action_layer_domains_and_history`

---

## Step 2: Add DOM references to `dom-cache.js`

Add new elements to:
- `cellucid/assets/js/app/ui/core/dom-cache.js`

This keeps selectors centralized and avoids repeated `getElementById` calls scattered across modules.

Pattern:
- extend the relevant domain bucket (e.g. `dom.render`, `dom.view`, or create a new bucket).

---

## Step 3: Create the module file under `ui/modules/`

Create a new module file, for example:
- `cellucid/assets/js/app/ui/modules/my-feature-controls.js`

Follow the established shape:

- export a single initializer:
  - `export function initMyFeatureControls({ state, viewer, dom, callbacks }) { ... }`
- do all DOM event wiring inside init
- return an object containing:
  - render functions (optional)
  - cleanup/destroy function (recommended if you add listeners outside stable DOM nodes)

Design rules:
- **One module owns one DOM subtree**.
- **A control dispatches the action it stands for and calls no mutating state
  method**: `dispatchFromControl(type, params, {failure, category})`
  (`app/actions/control-dispatch.js`), which says a refusal or a failure once in
  a notification and resolves `null` for it; for a slider or a color picker,
  `createLatestDispatch(type)`, so a drag dispatches only its latest value and
  the action's `coalesce` makes it one step. Name targets as the action's
  schema does (a field by `{source, key}`, a category by its label, a view
  `"live"` or `"View N"`), never by position.
- If your module implements an action (its orchestration lives here), provide
  it when the module starts (`getActions().provide(type, handler)`) and call
  the returned remover in `destroy`.
- A change your module's owner makes before it can be recorded (a gesture)
  is recorded with `getActions().recordObserved(type, params)`, synchronously
  from its event.
- Avoid direct cross-module DOM writes; communicate via state events or coordinator callbacks.
- If you subscribe to `state.on(...)`, store unsubscribe functions and clean up if you support re-init.

---

## Step 4: Wire the module into the UI coordinator

In:
- `cellucid/assets/js/app/ui/core/ui-coordinator.js`

1) Import your init function.
2) Call it during `initUI(...)` with the dependencies it needs:
   - `state` if it manipulates app state
   - `viewer` if it controls rendering/camera
   - the relevant `dom` bucket from `collectDOMReferences`
   - optional callbacks for cross-module updates

If your module changes visibility, fields, or pages:
- dispatch the action, and let state events trigger existing coordinator updates.

---

## Step 5: Style it using the design system

Cellucid uses a CSS design system:
- `cellucid/assets/css/README.md`

Rules of the road:
- no inline styles in `index.html`
- use utilities + semantic tokens
- avoid raw hex outside token files
- build buttons with `createButton` / `createIconButton`
  (`app/ui/components/controls.js`), and make them quieter the deeper they
  sit: a panel's, a dialog's or a card's own actions are `secondary` (the one
  that commits `primary`); an action inside a list item's row or details, a
  chat turn or a card's row is `quiet` (its one main action `quiet-primary`,
  accent text, never filled), one line; a well-known action there (delete,
  rename, close, more) is an icon button with an `ICONS` glyph, never a typed
  `⋯` or `✕`. A module that draws such rows is listed in
  `tests/nested-actions-quiet-contract.test.mjs`, which fails on a pill or a
  filled button in it ("Nested actions are quiet" in `assets/css/README.md`)

Every screen of the site follows the same design principles (the full list is
"Site Design Principles" in that README), whether it is a panel, a dialog, a
standalone page (`404.html`, `privacy.html`, `terms.html`, `accessibility.html`,
`assets/auth/callback.html`) or a screen shown when the app cannot run:
- tokens only, in light and dark, in the theme the visitor chose in the app;
- the brand card's composition (`components/_welcome.css`): wordmark, one calm
  heading, one short sentence, clear actions, quiet links;
- tidy and never verbose: say what happened and what to do; technical detail
  goes under **Details**, to the console and to the notification log;
  explanations go behind (i) tips and in the documentation;
- never a dead end, and terminal only when the app itself cannot run;
- phone widths, a visible keyboard focus, contrast in both themes, and
  `prefers-reduced-motion`.

Standalone pages load `assets/css/page.css`, a subset of `main.css` (tokens,
themes, base, the brand card and `layouts/_page.css`);
`tests/site-pages-design-contract.test.mjs` holds them to it.

If you add new tokens:
- update `types/design-tokens.d.ts`
- run validation scripts (see {doc}`14_testing_ci_and_release_process`)

---

## Step 6: Decide whether it is undone and saved

Ask:
- If a user undoes, jumps to an earlier step or saves a session, would they
  reasonably expect this feature to come back?

If yes:
- ensure the state is represented declaratively (not as DOM-only ephemeral
  state) and is held by a state domain: a setting (Step 1) needs nothing more;
  anything else is a domain's value or part of one, captured and restored by
  its owner in `cellucid/assets/js/app/state/domains/`
  ({doc}`29_action_layer_domains_and_history`, “How to add a piece of state”);
  the session file stores it with no encoder of its own
- add the action to its area's browser spec, which proves its `touches` by
  apply, undo and redo on the page's own history

If no:
- explicitly exclude the UI subtree with `data-state-serializer-skip="true"`
- document the exclusion in {doc}`10_sessions_persistence_and_serialization`
  (“What a session does not keep”)

---

## Step 7: Add troubleshooting hooks (always)

Before shipping, write down the top failure modes:
- missing DOM ids
- event handler attached twice
- state changes not reflected in viewer
- a change that is not a history step (a control calling a state method), or
  one undo does not put back (`touches` too narrow)
- session restore mismatches

Add at least:
- a console warning when required DOM elements are missing
- a clear user-facing notification when an action fails (use NotificationCenter)

---

## Troubleshooting (common mistakes)

### Symptom: “My new panel shows up but controls do nothing”

Likely causes:
- module init never ran (not wired into `ui-coordinator.js`)
- DOM ids in `index.html` don’t match `dom-cache.js`

Fix:
- verify init function is imported and called
- verify DOM cache returns non-null references

### Symptom: “My control works once, then doubles every reload”

Likely cause:
- event listeners are being registered multiple times (re-init without cleanup).

Fix:
- ensure `initUI` runs only once per page load
- if re-init is possible, return cleanup functions and call them on re-init paths

---

Next: {doc}`17_extension_point_add_analysis_mode` or {doc}`18_extension_point_add_export_renderer`.
