# Contributing to Cellucid (R package)

The `cellucid-r` package exports data to the same on-disk format used by the Cellucid viewer.

## Where the R contribution guide lives

`cellucid-r/CONTRIBUTING.md`, in the
[`cellucid-r`](https://github.com/theislab/cellucid-r) repository, is the
**authoritative** guide for working on the R package: the development
environment, `devtools::test()` / `devtools::check()`, how the hand-written
`man/` pages and `NAMESPACE` are maintained, the design constraints, how to
validate an export end-to-end, pull request rules, and troubleshooting. It
ships beside the code it describes and is the file GitHub shows a contributor
when they open an issue or a pull request, so it is kept self-sufficient.

**This page is not a copy of that file, and neither file is generated from the
other.** It owns the part of R contribution that only exists on this
documentation site:

| Question | Answer lives in |
|---|---|
| Which repository does my change belong in? | `cellucid-r/CONTRIBUTING.md` and {doc}`../contributing` |
| How do I set up R, run the tests, and open a pull request? | `cellucid-r/CONTRIBUTING.md` |
| How are `man/` and `NAMESPACE` maintained? | `cellucid-r/CONTRIBUTING.md` |
| How do I write or change an R page on this site? | this page, below |
| Which URL form and version marker do I use across repositories? | {doc}`web_app` |

The split is deliberate, not an oversight, and it is not a mirror that fell out
of date. The documentation build tree contains only `cellucid-python`
(`.readthedocs.yaml` clones no siblings and declares no submodules), so this
page **cannot** `{include}` a file from the `cellucid-r` repository the way
{doc}`python` includes `cellucid-python/CONTRIBUTING.md` from its own
repository. The choice was therefore between hand-copying content across a
repository boundary — which drifts, silently, because no build or test spans
both repositories — and giving each file one audience and cross-linking. Each
fact below has exactly one home.

One consequence is worth stating plainly, because a previous revision of this
page got it wrong in a way that would have broken a contributor's checkout: the
R package runs **no documentation generator**. `man/*.Rd` and `NAMESPACE` are
hand-written and authoritative, `R/` contains no roxygen comments, and
`NAMESPACE` carries a `useDynLib()` registration that no generator would emit.
Running `devtools::document()` there removes it. `cellucid-r/CONTRIBUTING.md`
is where that workflow is specified; do not infer it from anywhere else.

---

## I want to work on docs for the R package

The documentation site is built from `cellucid-python/docs/`, so an R
documentation change on this site is a pull request against `cellucid-python`,
not against `cellucid-r`. The R pages are {doc}`../user_guide/r_package/index`
(`cellucid-python/docs/user_guide/r_package/`).

The R package's own documentation — `man/`, `vignettes/`, and
`cellucid-r/README.md` — lives in the `cellucid-r` repository and is changed
there. A behaviour change usually needs both: the help page beside the code and
the page on this site.

When one of those pages points at another page on this site, write the
reference as a `{doc}` role, never in bare backticks — a backtick is inline
literal text that Sphinx will not resolve, so a wrong target renders as plain
text and survives the `-W` build. The rule, the relative-path form, and where
backticks are still correct are in
{doc}`../user_guide/python_package/h_developer_docs/15_docs_development_and_style_guide`.

The reverse also holds: `cellucid-r/CONTRIBUTING.md` is read on GitHub, where a
`{doc}` role is meaningless, so references from that file are written as
repository-relative paths with their `.md` extension. A path in that file and a
role on this page are not drift; they are the same reference rendered for two
readers.

Build the site before opening the pull request — warnings are errors:

```bash
sphinx-build -W --keep-going -b html docs docs/_build/html
```

---

## Conventions that span every repository

The canonical URL forms (`www` for links, the bare apex for frozen identifiers)
and the `CELLUCID_VERSION` marker convention apply to `cellucid-r` as well as
to this site. They are stated once, in
{doc}`web_app`, and enforced by
`cellucid-python/tests/test_canonical_url_contract.py`.
