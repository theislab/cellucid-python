# API Reference Coverage (and the canonical API reference)

This section serves two audiences:

1) **End users** who want to use the Python package → start here:
- {doc}`api/index`

2) **Maintainers/contributors** who want to keep the documentation complete and consistent:
- {doc}`01_public_functions_and_classes`
- {doc}`02_error_messages_and_exceptions_document_patterns`
- {doc}`03_api_reference_structure_and_migration_plan`

The goal is that “API reference” does not just mean docstrings; it also means:
- copy/pasteable examples that actually match real workflows,
- edge cases that prevent silent misuse,
- and symptom→diagnosis→fix troubleshooting that works for both beginners and experts.

---

## What counts as “public API”?

For `cellucid-python`, the public surface area includes:
- everything exported from `cellucid` (see `cellucid.__all__`)
- the CLI entry point (`cellucid …`)

Internal modules (e.g. prefixed with `_`) may be referenced for explanation, but are not stability-guaranteed.

---

## Documentation conventions used in this section

- Each API page aims for *layered writing*:
  - **Fast path** (beginner-friendly)
  - **Practical path** (computational workflows)
  - **Deep path** (expert mental models + internals)
- Each API page includes:
  - edge cases
  - troubleshooting
  - (optional) screenshot placeholders with capture specs

See also:
- `cellucid/markdown/DOCUMENTATION_MASTER_GUIDE.md`
- `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`

```{toctree}
:maxdepth: 2
:hidden:

api/index
```

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
