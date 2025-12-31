# What is Cellucid?

**Audience:** everyone (wet lab → computational → developer)  
**Time:** 5–10 minutes  
**What you’ll learn:**
- What Cellucid is (and is not)
- The “mental model” that makes the UI predictable
- How the web app relates to `cellucid-python`
- What state exists (and what gets saved)

---

## Cellucid in one sentence

Cellucid is a **WebGL2-powered, browser-based viewer** for exploring large single‑cell (and other point‑cloud) datasets through **embeddings, fields, filtering, highlighting, analysis, and publication‑grade export**.

---

## What it does (three explanations, pick your depth)

::::{tab-set}

:::{tab-item} Wet‑Lab / Non‑Technical

You can think of Cellucid as a “microscope for embeddings”:

- you load a dataset (someone prepared it for you, or you use a demo),
- you color points by a biological label or a gene,
- you zoom/rotate to understand structure,
- you mark groups (highlights),
- and you export a figure or save a session to share.

Cellucid is designed so you can get results without writing code, but you can still collaborate with computational colleagues who prepare the data.
:::

:::{tab-item} Computational / Power User

Cellucid is an interactive viewer for an **exported dataset bundle** (or a server that serves one). It supports:

- multiple embeddings (1D/2D/3D) with fast switching,
- categorical + continuous metadata fields (`obs`) and gene expression (`var`),
- per‑view filters and per‑view highlights (to compare hypotheses),
- analysis modules that operate on the current highlights/selection,
- deterministic exports (figures, sessions) intended for reproducibility.

The core principle is that *the UI is a thin layer over a state machine*: if you know what state is “global” vs “per view”, you can predict the behavior.
:::

:::{tab-item} Developer / Deep

Cellucid (web) is a stateful UI + WebGL viewer:

- **Data model:** dataset → embeddings → fields → per‑view contexts → artifacts (figures/sessions).
- **Performance model:** big arrays (positions/colors/visibility) are GPU‑resident; UI changes try to avoid hot‑path allocations.
- **Persistence:** “Save State” writes a single-file `.cellucid-session` bundle containing versioned chunks.
- **Integration:** `cellucid-python` can export the required folder layout and can embed/drive the viewer from notebooks (hooks/events).

If you’re trying to reason about edge cases, start from the UI glossary (`a_orientation/04_ui_glossary_terminology`) and the “Core interactions” section (`c_core_interactions/index`).
:::

::::

---

## The Cellucid mental model (learn this once)

This diagram is the quickest way to understand “what is a thing” in Cellucid and what changes when you click.

<!-- SCREENSHOT PLACEHOLDER
ID: orientation-mental-model
Where it appears: What is Cellucid? → The Cellucid mental model
Capture:
  - Asset type: diagram (SVG preferred) OR a clean annotated screenshot
  - Content: dataset → embedding(s) → fields → views/snapshots → filters/highlights → analysis → exports/sessions
  - Goal: communicate scope + flow (not every UI detail)
Crop:
  - Include only the diagram (no browser chrome)
Redact:
  - Not applicable (use generic labels like `pbmc_demo`)
Annotations:
  - Keep it minimal; use a single directional flow with 6–8 nodes
Alt text:
  - Diagram of the Cellucid workflow from dataset to exports.
Caption:
  - The Cellucid mental model: dataset → embeddings → fields → views → highlights/filters → analysis → exports/sessions.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder diagram for the Cellucid mental model.
:width: 100%

The Cellucid mental model: dataset → embeddings → fields → views → highlights/filters → analysis → exports/sessions.
```

Key terms (you’ll see these throughout the docs):

- **Dataset**: the loaded data bundle (points + metadata + optional gene expression).
- **Embedding / dimension (1D/2D/3D)**: where points live in space (UMAP/tSNE/etc).
- **Field**: something you can color by (categorical/continuous metadata or gene expression).
- **View**: a panel showing the dataset with a specific dimension, camera, fields, filters, and highlights.
- **Live view vs snapshot view**: live is your “working” view; snapshots are “kept” panels for comparison.
- **Artifacts**: exports (figures) and saved sessions (`.cellucid-session`) for sharing/reproducibility.

---

## What Cellucid is not (important expectations)

Cellucid intentionally does **not** try to be everything:

- It is not a full analysis notebook replacement (it’s an interactive viewer with targeted analysis modules).
- It is not a data-cleaning tool; if your inputs are malformed (NaNs, wrong shapes, mismatched cell order), the correct fix is usually upstream in data preparation.
- It is not a genome browser; gene expression support is for coloring/analysis, not raw read inspection.
- It is not a secure multi-user database by default. For multi-user label workflows, use **Community Annotation** (GitHub-backed) rather than “sending a screenshot around”.

---

## Relationship to `cellucid-python` (how data gets into the browser)

Most non-demo workflows look like this:

1) **Prepare** a dataset in Python (usually an AnnData): export embeddings + fields + gene expression into a folder layout.
2) **Load** that folder (or a server that hosts it) in the Cellucid web app.
3) **Explore**: switch dimensions, color by fields, filter, highlight.
4) **Save / export**:
   - figures for papers,
   - a `.cellucid-session` bundle for reproducibility and sharing.

If you’re choosing a workflow, use the decision tree: `a_orientation/05_which_workflow_is_for_me_decision_tree`.

---

## State & persistence (what gets remembered)

Cellucid is a **stateful application**. Most UI actions modify state such as:

- camera position + navigation mode,
- active fields (what you are coloring by),
- filters/outlier thresholds,
- highlight groups/pages,
- view layout and per‑view contexts.

You have three common “persistence” layers:

1) **In‑memory state** (lost on refresh/crash)
2) **Saved session bundle** (`Save State` → `.cellucid-session`)  
   Intended to restore the app state later or share with collaborators.
3) **Small browser preferences** (e.g. theme/background) stored in `localStorage`  
   These are convenience settings; they are not a reproducibility mechanism.

If you care about reproducibility, prefer **session bundles + exported figures**, not “remembering what sliders you touched”.

---

## Troubleshooting (starter set)

If you’re blocked immediately:

- If you see “WebGL2 is required…” or the canvas is blank: go to `a_orientation/02_system_requirements`.
- If the UI loads but you can’t move the camera: go to `c_core_interactions/06_troubleshooting_core_interactions`.
- If you can move the camera but no data appears: go to `b_data_loading/index` (data loading section).

---

## Next steps

- New user: `a_orientation/03_quick_tour_60_seconds`
- Confused by terms: `a_orientation/04_ui_glossary_terminology`
- Ready to learn navigation + multiview: `c_core_interactions/index`
