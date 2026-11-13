# Developer Documentation (Web App Architecture)

This section is a **developer-facing map of how the Cellucid web app works**: where the important code lives, how state flows through the system, and how to debug/extend features without breaking performance.

It lives in the `cellucid-python` docs because most developers arrive from the Python side, but **the frontend code is in the `cellucid/` repo**.

:::{important}
If you are a wet-lab user or a general “how do I use the UI?” reader, this section is *not* the right starting point.
Start with the web app user guide landing page: {doc}`../index`.
:::

---

## Fast Path (Pick Your Goal)

| You want to… | Start here | Then read |
|---|---|---|
| Run the app locally and reproduce a bug | {doc}`02_local_development_setup` | {doc}`13_debugging_playbook` |
| Understand the big-picture architecture | {doc}`05_app_architecture_overview` | {doc}`06_state_datastate_and_events` → {doc}`07_rendering_pipeline_webgl_and_performance_notes` |
| Learn where the code lives (file map + entry points) | {doc}`01_codebase_map_and_entry_points` | {doc}`08_ui_modules_map` |
| Add a new sidebar feature/module | {doc}`16_extension_point_add_ui_module` | {doc}`15_extension_points_overview` |
| Add/extend analysis (plots, compute, UI modes) | {doc}`11_analysis_architecture` | {doc}`17_extension_point_add_analysis_mode` |
| Add or change an action, a piece of state, undo/redo or the history | {doc}`29_action_layer_domains_and_history` | {doc}`10_sessions_persistence_and_serialization` → {doc}`16_extension_point_add_ui_module` |
| Work on sessions (share/load) | {doc}`10_sessions_persistence_and_serialization` | {doc}`29_action_layer_domains_and_history` → {doc}`13_debugging_playbook` |
| Work on agents (chat bar, permissions, data guard, the relay, `cellucid mcp`) | {doc}`32_agents_in_the_tab` | {doc}`31_local_agent_connection` → {doc}`30_action_reference` |
| Work on figure export | {doc}`12_figure_export_architecture` | {doc}`18_extension_point_add_export_renderer` |
| Understand loading sources (local/remote/GitHub/Jupyter) | {doc}`09_data_loading_pipeline_and_caching` | {doc}`../b_data_loading/index` (user-facing view) |
| Change how data is read, priced, kept or drawn (bytes, link, share, gene sample) | {doc}`20_read_path_architecture` | {doc}`19_subsampling_specification` → {doc}`09_data_loading_pipeline_and_caching` |
| Find out why a subsystem is built the way it is (decisions, measurements, what was reverted) | [the decision records](#decision-records), {doc}`21_decisions_decoders_and_codecs` to {doc}`28_decisions_rendering_layers` and {doc}`33_decisions_actions_domains_and_history` to {doc}`38_decisions_testing_and_hosted_runners` | the reference page of the subsystem |
| Run the tests, or write one that holds on every runner | {doc}`14_testing_ci_and_release_process` | {doc}`38_decisions_testing_and_hosted_runners` |

---

## What This Section Covers

- **Codebase map**: how `cellucid/index.html` boots the app and which folders correspond to state, UI, rendering, and data loading.
- **State model**: `DataState`, events, per-view contexts, multi-page highlights, caches, and batch updates.
- **Rendering pipeline**: the WebGL viewer, buffers, LOD, overlays (smoke + velocity), and performance footguns.
- **Actions, state domains and the history**: the one action catalog
  every control, agent and Python or R host dispatches, the sixteen state
  domains, the history's trackers, undo, branches and autosave, and how the
  catalog reaches Python, R and agents.
- **Persistence**: the one current `.cellucid-session` format (session format
  4: the history in one file, or the state only, or the steps only), its
  all-or-nothing load beside the person's trackers, starting views,
  dataset-identity rejection, and My datasets' `.cellucid-library` file.
- **Agents**: one tool call for the chat bar, Claude and ChatGPT through the
  relay, and local apps through `cellucid mcp`; the person's permissions and
  confirm cards; what agents are sent of the data; the generated action
  reference.
- **Testing**: the commands, what CI runs, the rules of a contract test and
  of a browser spec, the fixed device every spec runs on, the drift guards.
- **Decision records**: for each subsystem, what it does now, the dated
  decisions behind it, the numbers measured, what was tried and taken out,
  and what is still open.
- **Debugging**: how to capture minimal repro steps + the exact logs that help fix issues quickly.
- **Extension points**: the safe, repeatable steps for adding a UI module, an action, a piece of state, an analysis mode, or an export renderer.

:::{note}
Privacy/security and telemetry live in a separate user-facing section:
{doc}`../o_accessibility_privacy_security/index`.
This developer section links to it, but does not replace it.

Two things there that contributors get wrong often enough to flag here:
the app **does** ship Google Analytics, gated on three production hostnames and
inert everywhere else; and two of the browser caches persist **dataset-derived**
content, not preferences. If you add persistence or egress, that page has to
change in the same commit.
:::

:::{important}
**This section has no screenshots, on purpose.** Every page here documents code
structure, data flow and invariants — things a picture of the interface cannot
show and that a picture would misdate on the next CSS change. Where structure is
genuinely hard to hold in the head, these pages use ASCII diagrams instead: they
live in the same file as the prose they explain, they diff, and they cannot drift
from the text around them. If you find a page here where a real screenshot would
help, it probably belongs in the user guide rather than in this section.
:::

(decision-records)=
## Decision records (why the code is the way it is)

A reference page says what the code does and the rule it keeps. A decision
record says why: every record has the same sections (what it does now,
decisions and why, measurements that justify it, tried and reverted, known
limits and open items, where the code and tests live). History that would
otherwise sit in a code comment (a CI run, a date, an incident, a number
measured on a machine) is written there, and the comment names the page.
Update a record in the same change as the code it describes.

| Subsystem | Reference | Record |
|---|---|---|
| Decoders and codecs (WebAssembly, fflate, the measured probe, the decode pool) | {doc}`20_read_path_architecture` | {doc}`21_decisions_decoders_and_codecs` |
| Link model and request planning (round trips, rate limits, kept parts) | {doc}`20_read_path_architecture`, {doc}`09_data_loading_pipeline_and_caching` | {doc}`22_decisions_link_model_and_request_planning` |
| Partial reading (the share-cost rule, the block ladder, the gene sample, kept indexes, progressive opening) | {doc}`20_read_path_architecture`, {doc}`19_subsampling_specification` | {doc}`23_decisions_partial_reading` |
| Graph edges at scale (forward adjacency, edges at a share, k-NN) | {doc}`20_read_path_architecture`, {doc}`07_rendering_pipeline_webgl_and_performance_notes` | {doc}`24_decisions_graph_edges_at_scale` |
| Exports (cells, fields, the copy of the open dataset) | {doc}`20_read_path_architecture` | {doc}`25_decisions_exports` |
| Served stores and notebook widgets (Python, R) | {doc}`20_read_path_architecture`, {doc}`19_subsampling_specification` | {doc}`26_decisions_served_store_and_widgets` |
| Memory model | {ref}`memory-budgets` | {doc}`27_decisions_memory_model` |
| Rendering layers (unread cells, zero gray, alpha layers, levels of detail) | {doc}`07_rendering_pipeline_webgl_and_performance_notes` | {doc}`28_decisions_rendering_layers` |
| Actions, state domains and the history | {doc}`29_action_layer_domains_and_history`, {doc}`10_sessions_persistence_and_serialization` | {doc}`33_decisions_actions_domains_and_history` |
| Agents (the chat bar, permissions, the data guard, the relay, local agents) | {doc}`32_agents_in_the_tab`, {doc}`31_local_agent_connection` | {doc}`34_decisions_agents` |
| The interface (sidebar and controls, dialogs, legends, phone layout, keyboard, notices, community annotation) | {doc}`08_ui_modules_map` | {doc}`35_decisions_interface` |
| Figure export | {doc}`12_figure_export_architecture` | {doc}`36_decisions_figure_export` |
| Analysis | {doc}`11_analysis_architecture` | {doc}`37_decisions_analysis` |
| Testing, hosted runners, the build and the code cache | {doc}`14_testing_ci_and_release_process`, {doc}`03_build_run_and_deployment` | {doc}`38_decisions_testing_and_hosted_runners` |

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
