# Extension point: add a UI module

This guide is a step-by-step recipe for adding a **new sidebar UI module** to the Cellucid web app.

It is intentionally “mechanical”: follow the steps and you’ll end up with a module that:
- has clear DOM ownership,
- integrates with `DataState` cleanly,
- does not regress performance,
- can be persisted (or explicitly excluded) in sessions.

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
   - It can be a UI module that calls existing state/viewer APIs.

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

### Session persistence decision (early)

Decide whether your UI controls should be persisted by sessions:

- If you want generic inputs persisted, **do nothing special** (the serializer can capture many `input/select` values by DOM id).
- If your section contains transient state (previews, auth tokens, network state), exclude it:
  - add `data-state-serializer-skip="true"` on the section or a subtree

Reference:
- `cellucid/assets/js/app/state-serializer/README.md`

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
- prefer calling state methods and let state events trigger existing coordinator updates.

---

## Step 5: Style it using the design system

Cellucid uses a CSS design system:
- `cellucid/assets/css/README.md`

Rules of the road:
- no inline styles in `index.html`
- use utilities + semantic tokens
- avoid raw hex outside token files

If you add new tokens:
- update `types/design-tokens.d.ts`
- run validation scripts (see {doc}`14_testing_ci_and_release_process`)

---

## Step 6: Decide whether it should persist in sessions

Ask:
- If a user saves a session, would they reasonably expect this feature to restore?

If yes:
- ensure state is represented declaratively (not as DOM-only ephemeral state)
- add session contributor coverage if needed:
  - `cellucid/assets/js/app/session/contributors/`

If no:
- explicitly exclude the UI subtree with `data-state-serializer-skip="true"`
- document the exclusion in {doc}`10_sessions_persistence_and_serialization`

---

## Step 7: Add troubleshooting hooks (always)

Before shipping, write down the top failure modes:
- missing DOM ids
- event handler attached twice
- state changes not reflected in viewer
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
