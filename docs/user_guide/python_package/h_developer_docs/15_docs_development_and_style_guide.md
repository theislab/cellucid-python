# Docs development and style guide

This page documents how the Cellucid docs site is built and how to write documentation that works for mixed audiences (wet lab, computational, developer).

This repo uses:
- Sphinx
- MyST Markdown (`.md`)
- MyST-NB for notebooks (`.ipynb`)

---

## Where the docs live

The documentation site is in:
- `cellucid-python/docs/`

Key files/folders:

- `cellucid-python/docs/index.md`: docs home page
- `cellucid-python/docs/user_guide/`: main documentation “book”
  - `web_app/`: user guide for the web app (even though code is in `cellucid/`)
  - `python_package/`: user guide for the Python package + CLI

Developer docs you are reading now:
- `cellucid-python/docs/user_guide/python_package/h_developer_docs/`

---

## Build the docs locally (the one command)

From `cellucid-python/`:

```bash
make -C docs html
```

Open:
- `cellucid-python/docs/_build/html/index.html`

Clean rebuild:

```bash
make -C docs clean html
```

---

## Writing style: layered pages for mixed audiences

Cellucid docs intentionally use layered writing (see `cellucid/markdown/DOCUMENTATION_MASTER_GUIDE.md`):

1) **Fast path**: minimal steps, minimal theory, “what to click / what to run”.
2) **Practical path**: parameters, data requirements, performance.
3) **Deep path**: architecture, schemas, edge cases, debugging.

For developer docs:
- keep a short “why this exists” at the top,
- then provide a step-by-step workflow,
- then include an explicit “do not break these invariants” list,
- then a troubleshooting section.

---

## Page templates (recommended)

Most pages should contain:

- **Audience** (who should read what)
- **Prerequisites**
- **Fast path**
- **Deep details**
- **Edge cases**
- **Troubleshooting**
- **Related pages**

This consistency is what makes the docs usable for both beginners and experts.

---

## Screenshots and figures

Cellucid is a UI-heavy project. Screenshots are often the difference between:
- a wet-lab user succeeding vs giving up,
- and a maintainer understanding a bug report vs guessing.

### Where to store images

Recommended:
- `cellucid-python/docs/_static/screenshots/<topic>/<filename>.png`

### How to add a screenshot placeholder

Use the MyST `{figure}` directive + an HTML comment production spec.

Example:

````md
<!-- SCREENSHOT PLACEHOLDER
ID: unique-id
Capture:
  - UI location: ...
  - Action to reach state: ...
Alt text:
  - ...
Caption:
  - ...
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot.
:width: 100%

Caption text.
```
````

The placeholder file used throughout the docs is:
- `cellucid-python/docs/_static/screenshots/placeholder-screenshot.svg`

Detailed guidance for captions/alt text/redaction:
- `cellucid/markdown/DOCUMENTATION_SCREENSHOTS_AND_FIGURES_GUIDE.md`

---

## Notebooks in docs (MyST-NB)

Notebooks under `docs/` are rendered via MyST-NB.

Current config:
- notebook execution is disabled in the doc build (`nb_execution_mode = "off"`)

That means:
- notebooks are documentation artifacts, not CI-executed tests.

Writing expectations for notebooks (especially important for Cellucid):
- be verbose and explicit (copy/pasteable steps),
- explain edge cases and failure modes,
- include a big troubleshooting section per notebook section,
- avoid embedding private data or identifiers.

If you need something to be tested in CI, write a unit test under `tests/` instead.

---

## Linking and cross-references

Prefer `{doc}` links over raw URLs for internal pages:

```md
See: {doc}`../h_developer_docs/12_debugging_playbook`
```

This keeps links stable when the site is reorganized.

---

## Troubleshooting

### Symptom: “Sphinx build fails with MyST errors”

Common causes:
- unclosed code fences,
- incorrect directive indentation,
- broken `{doc}` targets.

Fix:
1) run `make -C docs clean html`
2) search the error message for the filename + line

### Symptom: “My notebook renders strangely”

Confirm:
- the notebook outputs are not enormous,
- the notebook does not rely on execution during build.

If you need “live execution”, that is a separate design choice and should be discussed before enabling it.
