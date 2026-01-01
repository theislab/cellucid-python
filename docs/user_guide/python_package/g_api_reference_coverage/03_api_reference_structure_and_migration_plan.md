# API reference structure and migration plan

This page explains **where API reference content lives** and how we keep it organized.

Context:
- The historical API reference lived in `cellucid-python/docs/user_guide/python_package/api/`.
- The plan is to make `cellucid-python/docs/user_guide/python_package/g_api_reference_coverage/api/` the **canonical** home for API reference content.

Why:
- “API reference” is more than docstrings; it needs workflow context, edge cases, and troubleshooting.
- We want the API reference and the “coverage map” to live together so maintainers can keep it complete.

---

## Canonical layout (current)

```
cellucid-python/docs/user_guide/python_package/g_api_reference_coverage/
├── index.md
├── 01_public_functions_and_classes.md
├── 02_error_messages_and_exceptions_document_patterns.md
├── 03_api_reference_structure_and_migration_plan.md
└── api/
    ├── index.md
    ├── jupyter.md
    ├── server.md
    ├── export.md
    ├── viewers.md
    ├── adapters.md
    ├── sessions.md
    ├── vector_fields.md
    └── cli.md
```

---

## Migration strategy (how we avoid breaking links)

Recommended approach:
1. Treat `g_api_reference_coverage/api/` as canonical.
2. Keep the legacy `python_package/api/` pages as **thin wrappers** (link or include) for backward compatibility.
3. Update internal references and navigation to point to the canonical location.
4. Only delete legacy pages once external links have been audited/updated.

---

## Authoring rules (so pages stay consistent)

Each page in `g_api_reference_coverage/api/` should include:
- Fast path (copy/paste) for beginners
- Practical path (parameter decisions + workflow patterns)
- Deep path (mental model, environment constraints, internals if relevant)
- Edge cases (explicit list)
- Troubleshooting (symptom → diagnosis → fix)
- Screenshot placeholders where UI steps are described (optional but recommended)

See:
- `cellucid/markdown/DOCUMENTATION_MASTER_GUIDE.md`
- `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`

---

## When you add a new public API (maintainers)

Whenever you add/remove/rename a public symbol (or a CLI command):
1. Update `cellucid.__all__` (for Python API) and/or `pyproject.toml` scripts (for CLI).
2. Add/update a row in {doc}`01_public_functions_and_classes`.
3. Add/update the appropriate API page under `g_api_reference_coverage/api/`.
4. Add at least 3–5 troubleshooting entries for the new behavior.
