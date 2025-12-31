# Troubleshooting (Web App): Start Here

This page is a **triage map**:
what you see → what to check → where the real fixes live.

If you can load a demo dataset but not yours, it’s usually **data loading / dataset format**.

---

## 2‑minute triage (do this first)

1) Confirm WebGL2 works on this machine: {doc}`../a_orientation/02_system_requirements`  
2) Check you’re not in a “looks broken but isn’t” state:
   - **Visible cells**: do you have `0 visible` due to filters? → {doc}`../e_filtering/07_troubleshooting_filtering`
   - **Active view**: are you editing the intended panel (live vs snapshot)? → {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples`
3) Reload in a clean state (private window; extensions off)  
4) If only one dataset fails, try a known-good demo dataset (if available)  
5) Open DevTools → Console and copy the **first** error (for bug reports)

---

## Choose your problem (symptom picker)

| What you’re seeing | Start here |
|---|---|
| App won’t open / blank page / WebGL2 message | {ref}`Installation & environment <install-env>` |
| Dataset won’t load / spinner forever / missing embeddings | {ref}`Data loading failures <data-loading>` |
| Blank canvas after load / “context lost” / choppy FPS | {ref}`Rendering & GPU issues <rendering>` |
| Selection tools feel wrong / highlights disappear / counts don’t match | {ref}`Selection & highlight issues <selection-highlighting>` |
| Analysis is empty / DE or markers look wrong / windows don’t restore | {ref}`Analysis issues <analysis>` |
| Export fails / downloads blocked / SVG/PNG looks wrong | {ref}`Export issues <export>` |
| Community annotation Pull/Publish/consensus issues | {ref}`Community annotation issues <community-annotation>` |

---

## Deep-dive troubleshooting (where most fixes live)

- Data loading: {doc}`../b_data_loading/08_troubleshooting_data_loading`
- Rendering/performance: {doc}`../n_benchmarking_performance/07_troubleshooting_performance`
- Fields & legends: {doc}`../d_fields_coloring_legends/05_troubleshooting_fields_legends`
- Filtering: {doc}`../e_filtering/07_troubleshooting_filtering`
- Highlighting & selection: {doc}`../f_highlighting_selection/06_troubleshooting_highlighting`
- Cross-highlighting (if enabled in your build): {doc}`../g_cross_highlighting/05_troubleshooting_cross_highlighting`
- Analysis: {doc}`../h_analysis/10_troubleshooting_analysis`
- Vector field / velocity overlay: {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`
- Sessions (save/load/auto-restore): {doc}`../l_sessions_sharing/10_troubleshooting_sessions`
- Figure export: {doc}`../k_figure_export/07_troubleshooting_figure_export`
- Bug-report / developer playbook: {doc}`../p_developer_docs/13_debugging_playbook`

---

(always-check)=
## Before you assume it’s a bug

Most “it’s broken” reports reduce to one of these mismatches:

1) **Filters vs visibility**: you hid the cells you expected to interact with.  
   - {doc}`../e_filtering/07_troubleshooting_filtering`
2) **Active view**: you’re editing a different panel than the one you’re looking at (live vs snapshot).  
   - {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples`
3) **Membership vs visibility**: highlight pages store membership; filters only change visibility.  
   - {doc}`../f_highlighting_selection/01_highlight_mental_model`
4) **Dataset identity mismatch**: sessions/annotation won’t apply cleanly to a different dataset export/version.  
   - {doc}`../b_data_loading/06_dataset_identity_why_it_matters`

---

(install-env)=
## Installation & environment

Use this section when Cellucid doesn’t start, your environment blocks required capabilities (WebGL2, file picking, downloads), or you’re running locally and the app won’t boot.

**Fast checks**

1) Try a modern desktop browser (Chrome/Edge are often the least frustrating).  
2) Confirm WebGL2 works: {doc}`../a_orientation/02_system_requirements`  
3) Try a private window (rules out extensions and stale cached state).  
4) If embedded in an iframe/Jupyter, open Cellucid in a standalone tab (iframes can block file access, pointer lock, fullscreen).  

**Common symptoms → where to go**

| Symptom | What it usually means | Where to go |
|---|---|---|
| `WebGL2 is required but not supported in this browser.` | WebGL2 unavailable (device/policy) | {doc}`../a_orientation/02_system_requirements` |
| Blank page or “module script failed to load” | Boot/deployment or local-server issue | {doc}`../p_developer_docs/13_debugging_playbook` |
| File/folder picker doesn’t open | Browser restrictions / unsupported file API | {doc}`../b_data_loading/08_troubleshooting_data_loading` |
| Downloads blocked (“Save State”, “Export”) | Browser download permissions / extension interference | {doc}`../l_sessions_sharing/10_troubleshooting_sessions`, {doc}`../k_figure_export/07_troubleshooting_figure_export` |

**If you’re running locally (developer / power user)**

- Local setup/build: {doc}`../p_developer_docs/02_local_development_setup`, {doc}`../p_developer_docs/03_build_run_and_deployment`
- Minimal repro + logs: {doc}`../p_developer_docs/13_debugging_playbook`

---

(data-loading)=
## Data loading failures

This is for “my dataset won’t load”, “spinner forever”, “missing embeddings”, or “loaded but missing fields/genes”.

**First: identify your loading path**

- Export folder (recommended) → {doc}`../b_data_loading/03_browser_file_picker_tutorial`
- Server mode (`cellucid serve ...`) → {doc}`../b_data_loading/04_server_tutorial`
- Jupyter notebook integration → {doc}`../b_data_loading/05_jupyter_tutorial`
- GitHub-hosted exports → {doc}`../b_data_loading/02_local_demo_tutorial`

If you’re not sure which applies: {doc}`../b_data_loading/01_loading_options_overview`.

**Fast checks**

1) If demos fail too, this is likely an environment/browser issue → {doc}`../a_orientation/02_system_requirements`  
2) If you loaded a large `.h5ad` directly in the browser, switch to server mode → {doc}`../b_data_loading/04_server_tutorial`  
3) DevTools → Network: look for 404/CORS/timeout errors (don’t guess)  
4) “Loaded but empty” can be filters/visibility → {doc}`../e_filtering/07_troubleshooting_filtering`

**Go to the deep dive (primary)**

- {doc}`../b_data_loading/08_troubleshooting_data_loading`

Related:
- Missing fields/legends → {doc}`../d_fields_coloring_legends/05_troubleshooting_fields_legends`
- Missing velocity/vector overlay → {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`

---

(rendering)=
## Rendering & GPU issues

Use this section when Cellucid is slow, the canvas is blank, or you see WebGL-related errors (including “context lost”).

**Fast checks (safe first)**

1) Confirm WebGL2 is hardware-accelerated: {doc}`../a_orientation/02_system_requirements`  
2) Switch to **Points** render mode while debugging (smoke + heavy overlays can hide root causes).  
3) Reduce GPU load:
   - clear snapshots (single view),
   - disable heavy overlays (vector fields, bloom),
   - make the window smaller.
4) If you see “WebGL context lost”: reload, then re-try with reduced load.

**Go to the deep dives**

- Performance symptoms (FPS, stutter, context lost): {doc}`../n_benchmarking_performance/07_troubleshooting_performance`
- Velocity/vector overlay rendering/perf: {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`

If the issue reproduces on a demo dataset (or across multiple datasets), capture a bug report:
- {doc}`../p_developer_docs/13_debugging_playbook`

---

(selection-highlighting)=
## Selection & highlight issues

Selection/highlighting problems are often mismatches between:
- *membership* (what a highlight page/group contains),
- *visibility* (what filters make visible),
- and *view context* (live vs snapshot).

**Fast checks**

1) Check **Active filters** (most “missing highlights” are just hidden): {doc}`../e_filtering/07_troubleshooting_filtering`  
2) Confirm the active highlight page and that the group checkbox is enabled.  
3) Confirm you are interacting with the intended panel: {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples`  
4) Confirm the tool’s modifier keys (`Alt`, `Shift+Alt`, `Ctrl/Cmd+Alt`).  

**Go to the deep dives**

- Highlighting/selection (primary): {doc}`../f_highlighting_selection/06_troubleshooting_highlighting`
- Cross-highlighting (analysis plots → highlights, if enabled): {doc}`../g_cross_highlighting/05_troubleshooting_cross_highlighting`

---

(analysis)=
## Analysis issues

Analysis issues usually come from one of four root causes:
1) no (or empty) highlight pages/groups to analyze,
2) confusing visibility (filters) with membership (highlight pages),
3) missing required data (gene expression, fields),
4) dataset/session mismatch.

**Fast checks**

1) Do you have a non-empty highlight page to analyze? → {doc}`../f_highlighting_selection/index`  
2) Are you on the view you think you are (live vs snapshot)? → {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples`  
3) Does gene search work (if your analysis is expression-based)? → {doc}`../d_fields_coloring_legends/05_troubleshooting_fields_legends`  

**Go to the deep dives**

- Analysis symptom map (primary): {doc}`../h_analysis/10_troubleshooting_analysis`
- If restore/sharing is involved: {doc}`../l_sessions_sharing/10_troubleshooting_sessions`

---

(export)=
## Export issues

Use this section when exports fail, downloads are blocked, or your PNG/SVG output doesn’t match what you expected.

**Fast checks**

1) Export a tiny PNG (e.g., ~600×450 at 150 DPI) to confirm downloads work.  
2) Confirm you are exporting the intended panel (live vs snapshot vs grid view).  
3) Confirm visibility: exporting “0 visible points” produces an empty figure.  
4) If large exports fail, cut size/DPI in half and try again.  
5) Check DevTools → Console for `[FigureExport]` messages.

**Go to the deep dives**

- Figure export troubleshooting (primary): {doc}`../k_figure_export/07_troubleshooting_figure_export`
- Quality knobs and sizing strategy: {doc}`../k_figure_export/04_quality_knobs_and_best_practices`

---

(community-annotation)=
## Community annotation issues

This section is for problems in **Community Annotation** (GitHub sync, Pull/Publish, consensus, “why don’t I see other people’s votes?”).

**Fast checks (fixes many reports)**

1) Confirm you are online (annotation is offline-first, but Pull/Publish require network).  
2) Confirm you are connected to the intended repo + branch (teams must agree on this).  
3) Click **Pull latest** before debugging anything else.  
4) If others can’t see your work: you must **Publish** (or your PR must be merged).  
5) If you see dataset mismatch behavior, stop and fix dataset identity first:
   - {doc}`../b_data_loading/06_dataset_identity_why_it_matters`
   - {doc}`../l_sessions_sharing/07_versioning_compatibility_and_dataset_identity`

**Deep dives (where the real fixes live)**

- Overview + “fast fix” map: {doc}`../j_community_annotation/index`
- Annotator workflow (vote/suggest/publish): {doc}`../j_community_annotation/01_annotator_guide`
- Author workflow (repo setup + config): {doc}`../j_community_annotation/02_author_guide`
- UI reference + error states: {doc}`../j_community_annotation/03_ui_reference`

Related:
- Session/identity troubleshooting (often confused with annotation state): {doc}`../l_sessions_sharing/10_troubleshooting_sessions`
- Bug report checklist / DevTools capture: {doc}`../p_developer_docs/13_debugging_playbook`

---

(bug-report)=
## Reporting a bug (what to include)

If you’re asking for help (or filing an issue), include:

1) **Where you ran Cellucid** (hosted URL vs local app vs embedded in Jupyter)  
2) **Browser + OS** (include versions)  
3) **Dataset source and loading path** (export folder / server mode / Jupyter / GitHub)  
4) **Dataset identity** (id/fingerprint, if available)  
5) **Exact steps to reproduce** (numbered; include whether snapshots/filters/highlights were involved)  
6) **First Console error** + any failing Network request (status code + URL)

Developer-focused capture checklist:
- {doc}`../p_developer_docs/13_debugging_playbook`
