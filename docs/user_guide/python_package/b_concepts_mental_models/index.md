# Concepts and Mental Models (Python-side)

This chapter is the “how to think about Cellucid” part of the `cellucid` Python package documentation.

It is written for **mixed audiences**:
- **Wet lab / non-technical**: clear mental models and “what happens when I click?” explanations (no assumption of Python expertise).
- **Computational users**: precise terminology, data-flow diagrams, and reproducibility guidance.
- **Power users / developers**: lifecycle, state, network + security model, and debugging playbook.

If you prefer learning by doing, you can skim this chapter and jump to:
- {doc}`../a_landing_pages/04_quick_start_3_levels`
- {doc}`../d_viewing_apis/index`
- {doc}`../c_data_preparation_api/index`

---

## The 30‑second mental model

Cellucid is a **web app** (the viewer UI). `cellucid-python` is a helper toolkit that:
- prepares data for the viewer (`prepare(...)`),
- serves data to the viewer (`serve(...)`, `show(...)`, `show_anndata(...)`),
- and optionally lets Python and the browser talk to each other (hooks/events, sessions).

In notebooks, a “viewer object” is just a **Python handle** for a running viewer instance.

```text
Python (kernel) ──starts──▶ local HTTP server ──serves──▶ web viewer (iframe/browser)
     ▲                                   │                     │
     │                                   │ fetch()             │ UI interactions
     │                                   ▼                     │
     └────────── events via POST ◀── /_cellucid/events ◀────────┘
```

---

## Reading order (recommended)

1) {doc}`01_what_is_a_viewer_object` (what you get back from `show(...)` / `show_anndata(...)`)
2) {doc}`02_data_flows` (what runs where: Python vs browser vs server)
3) {doc}`03_state_persistence_and_scope` (what persists, what resets, and why)
4) {doc}`04_dataset_identity_and_reproducibility` (stable IDs, exports, and “paper-ready” workflows)
5) {doc}`05_sessions_to_anndata_bridge` (pulling session state into Python and mutating AnnData safely)
   - Optional deep dive: {doc}`05_sessions_to_anndata_design`
6) {doc}`06_privacy_security_and_offline_vs_online` (what leaves your machine, and how to run offline)
7) {doc}`07_performance_mental_model_and_scaling` (what gets slow, why, and what to do)
8) {doc}`08_debugging_mental_model_where_to_look` (a systematic checklist for “it doesn’t work”)

---

## Pages in this chapter

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`device-desktop;1.5em;sd-mr-1` Viewer objects
:link: 01_what_is_a_viewer_object
:link-type: doc

What a viewer is, what it controls, lifecycle/cleanup, and how events map to Python.
:::

:::{grid-item-card} {octicon}`arrow-switch;1.5em;sd-mr-1` Data flows
:link: 02_data_flows
:link-type: doc

Where data lives (disk/server/browser), how the viewer loads it, and how hooks route events.
:::

:::{grid-item-card} {octicon}`history;1.5em;sd-mr-1` State & persistence
:link: 03_state_persistence_and_scope
:link-type: doc

Live viewer state vs durable session bundles vs reproducible export folders (and what resets when).
:::

:::{grid-item-card} {octicon}`shield-check;1.5em;sd-mr-1` Dataset identity
:link: 04_dataset_identity_and_reproducibility
:link-type: doc

How to choose stable dataset IDs, version exports, and avoid “session applied to wrong data”.
:::

:::{grid-item-card} {octicon}`repo;1.5em;sd-mr-1` Sessions → AnnData
:link: 05_sessions_to_anndata_bridge
:link-type: doc

Capture a `.cellucid-session` bundle into Python (no download) and apply highlights/fields to AnnData.
:::

:::{grid-item-card} {octicon}`law;1.5em;sd-mr-1` Sessions → AnnData (design)
:link: 05_sessions_to_anndata_design
:link-type: doc

Protocol and guardrails behind no-download session capture and safe AnnData application.
:::

:::{grid-item-card} {octicon}`lock;1.5em;sd-mr-1` Privacy & offline
:link: 06_privacy_security_and_offline_vs_online
:link-type: doc

What stays local, what can be exposed on a network, and how to make offline usage predictable.
:::

:::{grid-item-card} {octicon}`zap;1.5em;sd-mr-1` Performance
:link: 07_performance_mental_model_and_scaling
:link-type: doc

Export-time vs view-time costs, lazy loading, caching, and practical scaling advice.
:::

:::{grid-item-card} {octicon}`bug;1.5em;sd-mr-1` Debugging
:link: 08_debugging_mental_model_where_to_look
:link-type: doc

The fastest path to diagnosing viewer/server/hook issues, with concrete checks and fixes.
:::

::::

---

## Screenshots and figure placeholders

When this chapter references a UI step in the Cellucid web app, it uses **screenshot placeholders** with explicit capture instructions.
If you want to add images, follow the playbook in:
- `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
