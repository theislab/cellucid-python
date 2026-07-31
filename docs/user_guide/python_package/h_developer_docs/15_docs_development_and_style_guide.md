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

Cellucid docs use layered writing:

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

### How to add a verified screenshot

Use the MyST `{figure}` directive with a committed, visually reviewed image.

Example:

````md
```{figure} ../../../_static/screenshots/figure_export/style.png
:alt: Figure Export Style controls for background, typography, and text sizes.
:width: 224px

Style controls expose background, typography, and text sizing.
```
````

Before committing, reproduce the documented state, inspect the image at full
resolution, remove private identifiers, and confirm that the caption describes
only what is visibly present.

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

**Every reference to another page on this site is written as a `{doc}` role.**
Never name a page in bare backticks.

```md
See: {doc}`12_debugging_playbook`
See: {doc}`../../web_app/e_filtering/index`
```

This is not a style preference. A bare backtick is inline literal text, so
Sphinx does not resolve it, does not check that the target exists, and does not
turn it into a link. A `{doc}` role does all three, and the docs build runs with
`-W`, so a `{doc}` target that does not exist **fails the build**. Written in
backticks, the same wrong target renders as plain grey text and ships silently:
a reference to `k_sessions_sharing/index` survived in the highlighting section
long after the directory was renamed to `l_sessions_sharing/`, because nothing
was ever asked to resolve it.

Rules:

- **Targets are relative to the page doing the linking**, exactly as a file path
  would be: `06_troubleshooting_highlighting` for a sibling page,
  `../e_filtering/index` for another section, `../../python_package/index` for
  the other guide. Do not write the path from the `docs/` root or from a section
  root — Sphinx resolves it against the current document, so
  `` {doc}`a_orientation/04_ui_glossary_terminology` `` written *inside*
  `a_orientation/` points at a page that does not exist.
- **Use custom link text when the page title does not fit the sentence.** The
  default renders the target page's title, which is usually what you want:

  ```md
  Load it over HTTP with {doc}`Python server <04_server_tutorial>`.
  ```

- **Keep bare backticks for things that really are literals**: filenames and
  paths (`dataset_identity.json`), JSON keys and API parameters
  (`vector_fields`), CSS classes, code identifiers, and literal values
  (`local-user`). Converting one of those to `{doc}` is a regression that no
  test will catch — `vector_fields` happens to collide with a generated API page
  name, so it *would* resolve, and it would still be wrong.

Section `index.md` pages also list their children in a `{toctree}`; that is
navigation, and it does not remove the need for `{doc}` roles in prose.

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
