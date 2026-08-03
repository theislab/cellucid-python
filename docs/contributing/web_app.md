# Contributing to Cellucid (Web App)

Cellucid is a web app (browser UI + state layer + WebGL renderer) hosted at `https://www.cellucid.com`.

## Where the web-app contribution guide lives

`cellucid/CONTRIBUTING.md`, in the [`cellucid`](https://github.com/theislab/cellucid)
repository, is the **authoritative** guide for working on the web app: which
repository a change belongs in, what to put in a bug report, how to run the app
locally, where the code lives, the coding and performance rules, how to run the
tests, and what a good pull request looks like. It ships beside the code it
describes and is the file GitHub shows a contributor when they open an issue or
a pull request, so it is kept self-sufficient.

**This page is not a copy of that file, and neither file is generated from the
other.** It owns the part of web-app contribution that only exists on this
documentation site:

| Question | Answer lives in |
|---|---|
| Which repository does my change belong in? | `cellucid/CONTRIBUTING.md` and {doc}`../contributing` |
| How do I report a UI or rendering bug? | `cellucid/CONTRIBUTING.md` |
| How do I run, test, and open a pull request against the web app? | `cellucid/CONTRIBUTING.md` |
| Where does the web app's code live, and what are the performance rules? | `cellucid/CONTRIBUTING.md` |
| How do I write or change a web-app page on this site? | this page, below |
| Which URL form and version marker do I use across repositories? | this page, below |

The split is deliberate, not an oversight, and it is not a mirror that fell out
of date. The documentation build tree contains only `cellucid-python`
(`.readthedocs.yaml` clones no siblings and declares no submodules), so this
page **cannot** `{include}` a file from the `cellucid` repository the way
{doc}`python` includes `cellucid-python/CONTRIBUTING.md` from its own
repository. The choice was therefore between hand-copying content across a
repository boundary — which drifts, silently, because no build or test spans
both repositories — and giving each file one audience and cross-linking. Each
fact below has exactly one home.

---

## I want to work on docs for the web app

The documentation site is built from `cellucid-python/docs/`, so a web-app
documentation change is a pull request against `cellucid-python`, not against
`cellucid`. The web-app pages are:

- end-user documentation: {doc}`/user_guide/web_app/index`
  (`cellucid-python/docs/user_guide/web_app/`)
- developer documentation: {doc}`/user_guide/web_app/p_developer_docs/index`
  (`cellucid-python/docs/user_guide/web_app/p_developer_docs/`)

If you change UI or behaviour in the `cellucid` repository, update the matching
page here in the same change, so the description users read on ReadTheDocs
matches the app they are running.

When one of those pages points at another page on this site, write the
reference as a `{doc}` role, never in bare backticks — a backtick is inline
literal text that Sphinx will not resolve, so a wrong target renders as plain
text and survives the `-W` build. The rule, the relative-path form, and where
backticks are still correct are in
{doc}`/user_guide/python_package/h_developer_docs/15_docs_development_and_style_guide`.

The reverse also holds: `cellucid/CONTRIBUTING.md` is read on GitHub, where a
`{doc}` role is meaningless, so references from that file are written as
repository-relative paths with their `.md` extension. A path on this page and a
role in that file are not drift; they are the same reference rendered for two
readers.

Build the site before opening the pull request — warnings are errors:

```bash
sphinx-build -W --keep-going -b html docs docs/_build/html
```

The web-app testing, CI, and release checklist is
{doc}`/user_guide/web_app/p_developer_docs/14_testing_ci_and_release_process`.

---

## Cross-repository conventions

These two conventions span every Cellucid repository. Both are enforced, so a
change that ignores them fails a check rather than drifting quietly.

### Canonical URL forms: `www` for links, bare apex for identifiers

`cellucid.com` is written two ways on purpose. They are not interchangeable and
neither is a typo for the other.

| Form | Use it for | Why |
|---|---|---|
| `https://www.cellucid.com` | Every **navigational** URL: documentation links, README links, `CITATION.cff` `url:`, package homepages, CORS origin the app is served from | This is the host `cellucid/CNAME` publishes and the host the app’s own `<link rel="canonical">` declares, so it is where a browser should land |
| `https://cellucid.com/…` | Every **frozen identifier**: JSON Schema `$id`, XML namespace, exported-figure provenance strings | Identifiers must never change once published. A schema `$id` is a name, not an address — it is compared as a string and is not required to resolve |

The identifiers currently frozen at the bare apex are:

- `https://cellucid.com/contracts/community-annotation/user-v1.schema.json`
- `https://cellucid.com/contracts/community-annotation/config-v1.schema.json`
- `https://cellucid.com/contracts/community-annotation/merges-v1.schema.json`
- `https://cellucid.com/ns#` (SVG export namespace)
- `https://cellucid.com` as the figure-export `generator` / `Website` provenance
  string embedded in every exported PNG and SVG

Rewriting any of those to the `www` form is a **breaking change**: it
invalidates already-published schema references and makes figures exported by
two builds disagree about their own origin. Adding a `www` prefix to them is a
bug, and so is using the bare apex for a link a reader is meant to click.

Separately, `https://cellucid.com` and `https://*.cellucid.com` appear in the
Python server’s CORS documentation as **origin patterns**. Those describe what
`CORSMixin._get_allowed_origin` matches (it accepts the apex and any
`.cellucid.com` subdomain over HTTPS), so they are neither links nor
identifiers and stay as written.

Enforcement: `cellucid-python/tests/test_canonical_url_contract.py` sweeps the
documentation site and the Python sources, requires every navigational URL to
use `www`, and positively pins each frozen identifier at the apex — so drift in
either direction fails.

### `CELLUCID_VERSION` markers and the cross-repository sweep

Every hardcoded release version that a version bump must touch carries a
`CELLUCID_VERSION` marker comment on its line, so `grep -rn CELLUCID_VERSION`
lists the full edit set. `DESCRIPTION` in `cellucid-r` is DCF and cannot carry a
comment, so it declares the convention with a
`Config/cellucid/version-marker: CELLUCID_VERSION` field instead.

Markers only help inside one checkout. To confirm that all six repositories
agree, run the sweep from a workspace that contains them side by side:

```bash
python cellucid-python/scripts/sweep_versions.py
```

It enumerates every version declaration in every repository, groups them by
value with file and line, and separates three things that must not be confused:
declarations that have to agree, deliberately independent versions (the npm
name-reservation placeholder at `0.0.1`), and artifact stamps recording which
exporter release produced a committed export. It exits non-zero only on real
drift.
