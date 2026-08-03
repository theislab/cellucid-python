# Troubleshooting (Web App): Start Here

This page is a **triage map**: what you see → what it means → where the real fix lives.

There are two ways in.

- **If Cellucid told you something**, go to {ref}`find-the-message`. It lists the
  things the app actually says, grouped by what they mean, and points each one
  at the page that explains it. This is the fastest route.
- **If nothing is on screen and it just feels wrong**, use the
  {ref}`symptom picker <symptom-picker>` below.

If you can load a demo dataset but not yours, it’s usually **data loading / dataset format**.

---

## 2-minute triage (do this first)

1) Confirm WebGL2 works on this machine: {doc}`../a_orientation/02_system_requirements`  
2) Check you’re not in a “looks broken but isn’t” state:
   - **Visible cells**: do you have `0 visible` due to filters? → {doc}`../e_filtering/07_troubleshooting_filtering`
   - **Active view**: are you editing the intended panel (live vs snapshot)? → {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples`
   - **Still starting**: if controls look greyed out, read the line under the
     dataset picker — see {ref}`starting-up`.
3) Reload in a clean state (private window; extensions off)  
4) If only one dataset fails, try a known-good demo dataset  
5) Open DevTools → Console and copy the **first** error (for bug reports)

---

(find-the-message)=
## Find the message you saw

Cellucid’s failure messages are deliberately specific: several messages that
used to say the same unhelpful thing now say different things, because they
mean different things. Matching the wording is worth the effort — it usually
identifies the cause on its own.

### The app refused to start at all

A full-screen card headed **“Cellucid could not start”** means the application
never reached a usable state. The line beneath it is the reason and the footer
asks you to correct the launch configuration or server response and reload.

| The reason says | What it means | Go to |
|---|---|---|
| `WebGL2 is required but not supported in this browser.` | The browser gave the page no WebGL2 context: old browser, disabled hardware acceleration, a managed-device policy, or a remote-desktop session | {doc}`../a_orientation/02_system_requirements` |
| Anything about a manifest, a dataset, or a server response | The launch parameters pointed at something that could not be read | {doc}`../b_data_loading/08_troubleshooting_data_loading` |
| A blank page instead of a card, or “module script failed to load” | The app’s own files were not served correctly (wrong MIME type, opened via `file://`) | {doc}`../p_developer_docs/02_local_development_setup` |

If JavaScript is disabled you get a plain sentence instead: *“Cellucid requires
JavaScript and WebGL2 to run.”*

---

### A dataset would not load

These four messages look similar and mean four different things. This
distinction is the single most useful thing on this page.

| The message contains | Cause | What to do |
|---|---|---|
| **“was cancelled”** | The request was aborted — you switched datasets, disconnected, or navigated away mid-load. **Nothing is wrong.** | Retry the load. If it cancels on its own, something else is switching the source underneath you. |
| **“is larger than the metadata size ceiling”** | A size limit was hit, not a corrupt file | See {ref}`size-limits` below |
| **“response body transfer failed”** | The connection dropped partway through | Network, VPN, proxy or a server that closed early. Retry, then check the server. |
| **“must contain valid JSON”** | Genuinely malformed content — this now means what it says | The file is not what it claims to be. {doc}`../b_data_loading/08_troubleshooting_data_loading` |

:::{note}
Until recently the first, second and third of those all reported themselves as
an invalid format. If you are working from older notes or an older bug report
that says “invalid format”, re-run it: the message you get now is likely to be
different and more specific.
:::

When a connection or catalogue fails, the interface writes a full sentence of
the form *“&lt;what&gt; &lt;failed how&gt;: &lt;cause&gt; &lt;what to do&gt;”*.
The cause phrase tells you which class of failure it was:

| The sentence says | Meaning | Go to |
|---|---|---|
| “nothing is published at that address” | 404 — the path is wrong | {doc}`../b_data_loading/08_troubleshooting_data_loading` |
| “the server refused access” | 401 or 403 — it exists but you cannot read it | Check that the data is shared publicly, or that you are signed in |
| “the server rejected the request” | Some other 4xx | Check the address |
| **422**, with a body beginning `{"error": "non_finite_continuous_payload"` | The address was right. The server examined a gene or continuous column, found values it cannot publish, and refused to send them | {doc}`../d_fields_coloring_legends/05_troubleshooting_fields_legends` |
| “the server reported a problem of its own” | 5xx — not your fault | Wait and retry |
| “not in a format Cellucid can read” | The bytes arrived but do not match the export format | Re-export with `cellucid prepare`; {doc}`../b_data_loading/07_folder_file_format_expectations_high_level_link_to_spec` |
| “nothing readable came back” | Offline, DNS failure, **CORS**, or an address that simply is not a Cellucid export | {doc}`../b_data_loading/08_troubleshooting_data_loading` |

:::{tip}
**CORS is never named on the dataset path**, because the browser does not tell
the page whether a request was blocked by CORS or simply answered with something
unreadable. That is what “nothing readable came back” covers. If you host your
own exports and the address is definitely right, check the CORS headers first —
open DevTools → Network and look for a request with no response body.
:::

:::{note}
**422 is the one 4xx that is not about the address.** When a gene or continuous
obs column is requested from the Python AnnData server, the server validates it
before sending anything, and refuses a column containing `NaN`, an infinity, or
a value outside the `float32` range. The refusal is HTTP 422 with a JSON body —
`error`, `kind`, `name`, `counts`, `examples` and a `message` — and the viewer
quotes that body verbatim above its own explanation. Repair the values at the
source; nothing in the viewer can place an infinity on a colour scale.
{doc}`../d_fields_coloring_legends/05_troubleshooting_fields_legends` walks
through it.
:::

(size-limits)=
### It says something exceeded a limit

Browser tabs have finite memory, so Cellucid refuses work it cannot finish
rather than crashing halfway. Every one of these is a deliberate refusal with a
fixed ceiling.

| Limit | Applies to |
|---|---|
| **512 MiB** | one direct-in-browser `.h5ad` file; any single materialised array; a declared payload; a Zarr array; the connectivity working set |
| **64 MiB** | one metadata JSON document; one decoded Zarr chunk |
| **4 MiB** | Zarr metadata |
| **64 GiB** | a Zarr ZIP archive, and its total uncompressed size |

The fix is almost always the same: **stop loading it in the browser**. Use
{doc}`../b_data_loading/04_server_tutorial` for large `.h5ad` or `.zarr`, or
prepare an export with the Python or R package.

:::{note}
Several of these messages print the ceiling as a raw byte count —
`67108864` is 64 MiB, `536870912` is 512 MiB, `68719476736` is 64 GiB.
:::

---

### Something on the plot did not respond

| The message says | What it means | Go to |
|---|---|---|
| “That click did not land on a cell, so there is nothing to select. Zoom in to make cells easier to hit: only the cells the viewer is currently drawing can be clicked.” | Picking searches a fixed radius around the ray in data space. It has nothing to do with how big the points look, so **changing Point size will not help** — zooming will, and so will turning off level of detail | {doc}`../f_highlighting_selection/06_troubleshooting_highlighting` |
| “That cell has no value on the active field, so there is nothing to select.” | The cell is real but missing on the field you are colouring by | {doc}`../d_fields_coloring_legends/03_color_by_behavior` |
| “That gesture ended outside the view, so nothing was selected.” | The drag was released off the panel | {doc}`../f_highlighting_selection/02_selection_tools_document_each_tool` |
| “That gesture did not change the selection.” | The gesture worked; it just selected what was already selected | — |
| “Annotation selection needs an active field.” | Choose a categorical or continuous obs field first | {doc}`../d_fields_coloring_legends/02_field_selector_ux` |

:::{warning}
If you find advice anywhere telling you to **raise Point size so clicks land**,
it is wrong and out of date. Point size changes what is *drawn*; the pick radius
is fixed. Zoom in instead.
:::

---

### WebGL stopped

| The message | Meaning | Go to |
|---|---|---|
| **“WebGL Context Lost”** dialog, “WebGL context lost. Reload required to continue.”, with a **Reload** button | The GPU dropped the context, usually from memory pressure. A reload is genuinely required — the viewer does not resume from a restored context | {doc}`../n_benchmarking_performance/07_troubleshooting_performance` |
| “WebGL context restored, but Cellucid requires a reload…” | Your browser handed the context back; Cellucid still needs a reload | Reload |
| “Smoke rendering failed: …” / “Smoke density could not be built: …” | Volumetric rendering could not run; the view keeps drawing points | {doc}`../c_core_interactions/03_render_modes_points_vs_volumetric_smoke` |
| “Volumetric smoke requires a single view. Clear snapshots first.” | Smoke mode is single-view only | {doc}`../c_core_interactions/04_view_layout_live_snapshots_small_multiples` |
| “Velocity overlay unavailable” | The overlay could not be built for this dataset or dimension | {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay` |

---

### Fields, genes and legends

| The message | Meaning | Go to |
|---|---|---|
| **“No gene matches”**, followed by your query and a count of how many gene names the dataset publishes | The gene is not in *this export*. The count is the point: an export often carries a chosen subset of the source data, so a real gene can be legitimately absent | {doc}`../d_fields_coloring_legends/05_troubleshooting_fields_legends` |
| “This dataset publishes no gene expression fields, so no gene panel can be built.” | There is no expression data at all here, so no marker or gene analysis is possible | {doc}`../h_analysis/08_analysis_mode_genes_panel` |
| `Gene "<name>" not found in dataset. Check spelling or try another gene.` | That specific gene is absent | {doc}`../d_fields_coloring_legends/05_troubleshooting_fields_legends` |
| **“Gene cannot be shown”**, a multi-line notification naming the gene, counting its `+Infinity` / `-Infinity` / `NaN` values and listing the first affected cells | The values arrived (or the server refused to send them) and cannot be placed on a colour scale. The active field never changed. This is a values problem, not a transport one — repair the data at the source | {doc}`../d_fields_coloring_legends/05_troubleshooting_fields_legends` |
| `Field "X" has no value in any of N cells, so it has no colour scale.` | Every value in that column is `NaN`, so there is nothing to scale — usually a column that was never computed | {doc}`../d_fields_coloring_legends/05_troubleshooting_fields_legends` |
| “Failed to load field: …” | A field could not be loaded. Read the rest of the line: if it counts non-finite values it is the case above, otherwise the data could not be fetched | {doc}`../d_fields_coloring_legends/05_troubleshooting_fields_legends`, then {doc}`../b_data_loading/08_troubleshooting_data_loading` |
| “Failed to load gene: …” / “Failed to load field: …” on the progress toast | The load never finished — a dropped connection, a stopped server, a refused payload. It is the terminal state of the `Loading …` toast; the reason arrives as a separate notification beside it | {doc}`../b_data_loading/08_troubleshooting_data_loading` |
| “No neighbor graph available for this dataset” | kNN selection needs connectivity, which this export does not include | {doc}`../f_highlighting_selection/02_selection_tools_document_each_tool` |

:::{note}
If you ever see a message containing a literal `{symbol}`, `{name}` or similar
in curly braces, that is a bug — please report it with
{doc}`../n_benchmarking_performance/09_reporting_performance_bugs`. One such
case existed and was fixed; the message that used to do it now reads *“This
dataset publishes no gene expression fields…”*.
:::

---

### Sessions, exports and analysis

| The message | Meaning | Go to |
|---|---|---|
| “Failed to load session”, or a sentence about a session manifest | The bundle did not validate. These messages are written for developers; the useful question is whether the bundle matches this dataset | {doc}`../l_sessions_sharing/10_troubleshooting_sessions` |
| “The dataset is loaded, but its published sample view could not be applied” | The dataset opened; only its saved default view did not | {doc}`../l_sessions_sharing/07_versioning_compatibility_and_dataset_identity` |
| “Export failed: …” | Figure export could not complete | {doc}`../k_figure_export/07_troubleshooting_figure_export` |
| “Exact point export requires WebGL2…” / “…requires camera matrices…” | A capability the exact-point renderer needs is missing right now | {doc}`../k_figure_export/03_export_formats_and_renderers` |
| “Analysis failed: …”, “Invalid form values” | An analysis could not run. “Invalid form values” does not name the field — check each input | {doc}`../h_analysis/10_troubleshooting_analysis` |
| “Low memory detected…” / “High memory usage detected…” | The analysis layer is shedding cached work | {doc}`../h_analysis/10_troubleshooting_analysis` |

---

### Community annotation and GitHub

| The message | Meaning | Go to |
|---|---|---|
| “Couldn’t reach the GitHub sign-in server from &lt;origin&gt;…” | Your deployment’s origin is not allowed by the authentication proxy. This one **does** name CORS | {doc}`../j_community_annotation/02_author_guide` |
| “GitHub session expired or was revoked.” | Sign in again | {doc}`../j_community_annotation/01_annotator_guide` |
| “Lost access to &lt;repo&gt;. GitHub returned 403” | Permissions changed on the repository | {doc}`../j_community_annotation/02_author_guide` |
| “Annotation repo not accessible (deleted/renamed or access removed)” | The repository moved or you lost access | {doc}`../j_community_annotation/02_author_guide` |
| “Repo &lt;repo&gt; is missing …” | The repository does not have the expected annotation layout | {doc}`../j_community_annotation/02_author_guide` |
| “Sign in required.”, “No annotation repo connected.”, “Missing dataset context.” | A precondition is not met yet | {doc}`../j_community_annotation/03_ui_reference` |
| “Connected with a dataset mismatch…” | The repository was set up against a different dataset identity | {doc}`../b_data_loading/06_dataset_identity_why_it_matters` |

---

### Notebook / Jupyter

| The message | Meaning | Go to |
|---|---|---|
| **“Notebook command failed”** — “The notebook &lt;type&gt; command did not run: … The viewer is unchanged; correct the notebook cell and run it again.” | A command sent from Python was rejected by the viewer. Nothing changed, so re-running a corrected cell is safe | {doc}`../b_data_loading/05_jupyter_tutorial` |
| **“Jupyter server disconnected”** — “The Python server stopped or returned an invalid health response. The last complete view has been frozen.” | The kernel or server went away. What you see is the last good state, not live | {doc}`../b_data_loading/05_jupyter_tutorial` |
| “Selection was saved locally, but Python notification failed…” | The selection is real in the browser; Python was not told | {doc}`../b_data_loading/05_jupyter_tutorial` |

---

(starting-up)=
### “The buttons are greyed out” / “I clicked and nothing happened”

The data-source, renderer and benchmark controls now ship **disabled** and are
enabled only once the code that responds to them is wired up. While the app is
starting, a line under the dataset picker reads:

> Starting Cellucid — the data controls below open as soon as it is ready.

When it changes to *“Data sources are ready…”* the controls are live. Screen
readers announce this, and each control is described by that line.

So:

- **Greyed out for a moment on a slow start** — expected. Wait.
- **A click that did nothing at all** — this should no longer be reproducible.
  If it happens, it is a bug worth reporting.
- **Never becomes ready** — that is a startup failure. Look for the “Cellucid
  could not start” card, or check the console.

---

(symptom-picker)=
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
| Keyboard shortcuts do nothing | {doc}`../o_accessibility_privacy_security/01_accessibility` |
| “Is my data being sent anywhere?” | {doc}`../o_accessibility_privacy_security/02_privacy_model` |

---

## Deep-dive troubleshooting (where most fixes live)

- Data loading: {doc}`../b_data_loading/08_troubleshooting_data_loading`
- Core interactions (camera, views, dimensions): {doc}`../c_core_interactions/06_troubleshooting_core_interactions`
- Rendering/performance: {doc}`../n_benchmarking_performance/07_troubleshooting_performance`
- Fields & legends: {doc}`../d_fields_coloring_legends/05_troubleshooting_fields_legends`
- Filtering: {doc}`../e_filtering/07_troubleshooting_filtering`
- Highlighting & selection: {doc}`../f_highlighting_selection/06_troubleshooting_highlighting`
- Analysis: {doc}`../h_analysis/10_troubleshooting_analysis`
- Vector field / velocity overlay: {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`
- Sessions (save/load): {doc}`../l_sessions_sharing/10_troubleshooting_sessions`
- Figure export: {doc}`../k_figure_export/07_troubleshooting_figure_export`
- Accessibility: {doc}`../o_accessibility_privacy_security/01_accessibility`
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
5) **What the export contains**: a gene or field that exists in your `.h5ad` need
   not exist in the export made from it.  
   - {doc}`../d_fields_coloring_legends/07_genes_in_the_built_in_samples`

---

(install-env)=
## Installation & environment

Use this section when Cellucid doesn’t start, your environment blocks required capabilities (WebGL2, file picking, downloads), or you’re running locally and the app won’t boot.

**Fast checks**

1) Use a current stable Chrome, Edge, Firefox, or Safari release.
2) Confirm WebGL2 works: {doc}`../a_orientation/02_system_requirements`  
3) Try a private window (rules out extensions and stale cached state).  
4) If embedded in an iframe/Jupyter, open Cellucid in a standalone tab (iframes can block file access, pointer lock, fullscreen).  

**Capabilities Cellucid requires, and what happens without them**

| Capability | If missing |
|---|---|
| JavaScript | A plain sentence replaces the app |
| **WebGL2** | Hard stop: “WebGL2 is required but not supported in this browser.” There is no software renderer |
| **Gzip decompression** (`DecompressionStream`) | Compressed payloads are refused before the request is made, with a message naming the missing support |
| Web Workers | Only the benchmark’s off-thread data generation needs this |
| `localStorage` | Community Annotation refuses to run; preferences stop persisting |
| IndexedDB | Caches degrade **silently** — repeat analyses just get slower |
| Clipboard write | “Copy” actions report that you must select the text manually |

Cellucid does **not** use the File System Access API, `SharedArrayBuffer`, or
`OffscreenCanvas`, so a browser or policy that blocks those is not the cause of
your problem.

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

1) Read the message and match it in {ref}`find-the-message` — the four
   “would not load” causes look alike and are not.  
2) If demos fail too, this is likely an environment/browser issue → {doc}`../a_orientation/02_system_requirements`  
3) If you loaded a large `.h5ad` directly in the browser, switch to server mode → {doc}`../b_data_loading/04_server_tutorial`  
4) DevTools → Network: look for 404/CORS/timeout errors (don’t guess)  
5) “Loaded but empty” can be filters/visibility → {doc}`../e_filtering/07_troubleshooting_filtering`

**A load that half-succeeded**

Messages of the form *“The dataset is loaded, but …”* mean exactly that: the
data is usable and one dependent step did not complete. Do not reload
reflexively — read which step failed.

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
2) Make the window smaller. This is the highest-yield single test.
3) Switch to **Points** render mode while debugging (smoke + heavy overlays can hide root causes).  
4) Clear snapshots. Then, if you need a comparison, note that the four-view
   grid is often cheaper than a two-view row —
   {doc}`../n_benchmarking_performance/06_edge_cases_performance` explains why.
5) If you see “WebGL context lost”: reload, then re-try with reduced load.

:::{note}
**Filtering cells away makes drawing faster, not slower.** Hidden cells are
rejected before anything is drawn for them. If hiding most of your data changes
nothing, the points were not the bottleneck — look at window size, overlays and
smoke settings instead.
:::

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
5) If clicks miss cells, **zoom in** — the pick radius is fixed and only drawn
   cells are hittable. Point size is not the knob.

**Go to the deep dives**

- Highlighting/selection (primary): {doc}`../f_highlighting_selection/06_troubleshooting_highlighting`

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
3) Does this export contain expression at all? If the panel says it publishes no
   gene expression fields, no amount of configuration will help → {doc}`../h_analysis/08_analysis_mode_genes_panel`

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

**Before you share the file**

Exported figures embed provenance metadata that can include the dataset name and
sometimes a local path or source URL. See
{doc}`../k_figure_export/05_metadata_and_provenance` and
{doc}`../o_accessibility_privacy_security/02_privacy_model`.

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
5) If you are in a private window, annotation will refuse to run — it needs
   browser storage.  
6) If you see dataset mismatch behavior, stop and fix dataset identity first:
   - {doc}`../b_data_loading/06_dataset_identity_why_it_matters`
   - {doc}`../l_sessions_sharing/07_versioning_compatibility_and_dataset_identity`

**Deep dives (where the real fixes live)**

- Overview + “fast fix” map: {doc}`../j_community_annotation/index`
- Annotator workflow (vote/suggest/publish): {doc}`../j_community_annotation/01_annotator_guide`
- Author workflow (repo setup + config): {doc}`../j_community_annotation/02_author_guide`
- UI reference + error states: {doc}`../j_community_annotation/03_ui_reference`
- End-to-end walkthrough: {doc}`../j_community_annotation/04_lifecycle_walkthrough`

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
6) **The exact wording of the message**, not a paraphrase — the wording is what
   distinguishes the causes above  
7) **First Console error** + any failing Network request (status code + URL)

For performance specifically, use the template in
{doc}`../n_benchmarking_performance/09_reporting_performance_bugs`.

Developer-focused diagnostic checklist:
- {doc}`../p_developer_docs/13_debugging_playbook`
