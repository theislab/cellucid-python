# Public functions and classes (coverage map)

This page is a **single source of truth** for:
- what we consider the public `cellucid` Python API,
- where each item is documented,
- and what “done” means for documentation coverage.

Public API definition:
- `cellucid.__all__` is the canonical list of exported symbols.
- The CLI entry point `cellucid` is part of the public API surface.

---

## Public API map (what users can import)

| Category | Symbol | Kind | Primary docs page |
|---|---|---|---|
| Data preparation | {func}`~cellucid.prepare` | function | {doc}`api/prepare` |
| Data preparation | {func}`~cellucid.restore_original_order` | function | {doc}`api/prepare` |
| Jupyter | {func}`~cellucid.show_anndata` | function | {doc}`api/jupyter` |
| Jupyter | {class}`~cellucid.AnnDataViewer` | class | {doc}`api/viewers` |
| Server | {func}`~cellucid.serve_anndata` | function | {doc}`api/server` |
| Server | {class}`~cellucid.AnnDataServer` | class | {doc}`api/server` |
| Web cache | {func}`~cellucid.get_web_cache_dir` | function | {doc}`api/server` |
| Web cache | {func}`~cellucid.clear_web_cache` | function | {doc}`api/server` |
| Sessions | {class}`~cellucid.CellucidSessionBundle` | class | {doc}`api/sessions` |
| Sessions | {func}`~cellucid.apply_cellucid_session_to_anndata` | function | {doc}`api/sessions` |
| Subsampling | {func}`~cellucid.subsample_cells` | function | {doc}`api/subsample` |
| Vector fields | {func}`~cellucid.compute_transition_drift` | function | {doc}`api/vector_fields` |
| Vector fields | {func}`~cellucid.add_transition_drift_to_obsm` | function | {doc}`api/vector_fields` |
| Version | `cellucid.__version__` | str | {doc}`api/index` |

CLI surface:
- `cellucid` (entry point: `cellucid.cli:main`) → documented in {doc}`api/cli`

---

## Documentation coverage checklist (definition of “complete”)

For each **public function/class**, documentation is “complete” when it contains:

### 1) A clear “why/when” explanation (mixed audience)
- 1–2 sentence definition (what it does)
- 2–5 use cases (wet lab + computational)
- “When NOT to use it” (common footguns)

### 2) A runnable minimal example (fast path)
- copy/pasteable snippet
- “what success looks like” description

### 3) A detailed reference section (practical path)
- parameters + defaults + types
- return value / side effects
- what files are written/served (if applicable)
- environment/dependency requirements

### 4) Edge cases (mandatory)
- shape mismatches
- missing keys / missing metadata
- large dataset behavior
- NaN/Inf and constant-value behavior (where relevant)

### 5) Troubleshooting (mandatory)
At least 8–15 symptom→diagnosis→fix entries for major entry points (`prepare`, `show_anndata`, `serve_anndata`, sessions).

### 6) Cross-links
- Link to the relevant user-guide chapter(s) where the long-form workflows live.

---

## How to keep this page in sync (maintainers)

Whenever you change `cellucid.__all__` or add a new CLI command:
1. Add/update the row in the “Public API map” table above.
2. Ensure there is a dedicated doc page (or section) in `g_api_reference_coverage/api/`.
3. Add/expand troubleshooting entries for the new behavior.
