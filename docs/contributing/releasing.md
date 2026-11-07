# Releasing Cellucid

Every Cellucid package is released by one habit: **publish a GitHub Release
named `vX.Y.Z` in its repository, and that repository's workflow publishes
it.** The web app and the Python and R packages share one version, and the
npm package is released under it too.

| To release | Publish a GitHub Release `vX.Y.Z` in | The workflow | Publishes to |
|---|---|---|---|
| The Python package `cellucid` | `theislab/cellucid-python` | `pypi-publish.yml` ("Publish to PyPI") | PyPI, through trusted publishing |
| The R package `cellucid` | `theislab/cellucid-r` | `release.yaml` ("Publish to CRAN") | CRAN: attaches the checked tarball to the release and submits it; the maintainer confirms CRAN's email |
| The npm package `cellucid` (a name-reservation placeholder) | `theislab/cellucid` | `npm-publish.yml` ("Publish to npm") | npm, through trusted publishing, with provenance |

Each workflow runs on `release: published` and on nothing else. Pushing a tag
publishes nothing, none of them has a manual trigger, and none reads a stored
token. The web app itself is not published by a release: `www.cellucid.com`
serves `main` of `theislab/cellucid`.

## Before releasing

From a workspace that holds the repositories side by side:

1. **Versions agree.** `python cellucid-python/scripts/sweep_versions.py`
   ends with `VERDICT: OK`; every declaration, `cellucid/npm/package.json`
   included, names the version being released.
2. **CI is green on `main`** of each repository being released (`test.yml` and
   `docs-check.yml` in `cellucid-python`, `R-CMD-check.yaml` in `cellucid-r`,
   `test.yml` in `cellucid`), on the commit the release will tag.
3. **The Python release contract holds.** In `cellucid-python`,
   `python scripts/validate_release.py --tag vX.Y.Z` prints
   `Cellucid Python release contract is valid for vX.Y.Z.`
   ({doc}`/user_guide/python_package/h_developer_docs/14_release_process`).
4. **The R package checks clean.** `testthat::test_local()` and
   `R CMD check --as-cran` pass, and `cran-comments.md` is current
   ({doc}`/user_guide/r_package/h_developer_docs/03_release_process`).
5. **The changelogs are ready and undated.** The current
   `cellucid-python/CHANGELOG.md` heading is `## [X.Y.Z]` with no date, and
   the current `cellucid-r/NEWS.md` heading is `# cellucid X.Y.Z`.

## Order

1. **Python and R**, in either order: `theislab/cellucid-python`, then
   `theislab/cellucid-r` (or the reverse). They are peers; neither waits for
   the other.
2. **npm**, once both have published: `theislab/cellucid`. The placeholder
   points at the Python and R packages, so it goes last.

When a session format changes, deploy `cellucid-datasets` before the web app.

In each repository: **Releases** → **Draft a new release** → **Choose a tag**
→ type `vX.Y.Z` → **Create new tag**, target `main` → a title and the
version's changelog section as notes → leave **Set as a pre-release**
unticked → **Publish release**. Then watch the workflow under **Actions**.

Every workflow refuses, before it builds or publishes anything:

- a draft release;
- a release whose **pre-release** mark does not match the version (an
  `X.Y.Z` version is never published from a release marked pre-release; the
  Python and R packages declare only `X.Y.Z` versions);
- a tag other than `v` plus the package's version (`pyproject.toml`,
  `DESCRIPTION`, `npm/package.json`);
- a tagged commit that is not on `main`.

## Re-running a release

A run that failed for a passing reason (a runner, the network, CRAN's upload
form, a trusted publisher set up late) is retried from the run itself:
**Actions** → the run → **Re-run jobs** (or **Re-run failed jobs**). It runs
the same release and the same tag again.

Do not re-run a CRAN run whose submission step already succeeded: that submits
the package to CRAN a second time. Answer CRAN's confirmation email instead.

When the source has to change instead, delete the GitHub Release and its tag,
fix `main`, and publish the release again. Once a registry has accepted a
version (any PyPI file, the npm version), it can never be replaced: release
the corrected source under the next version.

## Dates

- `cellucid-python/CHANGELOG.md`: the current version's heading stays
  undated until its GitHub Release is published (`validate_release.py`
  refuses a date before publication). The version bump that opens the next
  version dates it, `## [X.Y.Z] - YYYY-MM-DD`, with the day the Python
  release was published:
  `gh release view vX.Y.Z --repo theislab/cellucid-python --json publishedAt`
  (UTC).
- `cellucid-r/NEWS.md`: headings carry no date; CRAN records the
  publication day.
- `CITATION.cff` files carry no `date-released:`.

## One-time setup

Done once per registry, by an owner, outside the repositories:

- **PyPI**: the trusted publisher of the `cellucid` project names owner
  `theislab`, repository `cellucid-python`, workflow `pypi-publish.yml` and
  environment `pypi`, and the repository has a protected `pypi` environment
  ({doc}`/user_guide/python_package/h_developer_docs/14_release_process`).
- **CRAN**: nothing beyond the maintainer in `cellucid-r/DESCRIPTION`, who
  receives and confirms CRAN's email for each submission; the workflow asks
  for `contents: write` itself to attach the tarball
  (`cellucid-r/publishing.md`).
- **npm**: on npmjs.com, signed in as an owner of the package,
  **https://www.npmjs.com/package/cellucid** → **Settings** →
  **Trusted Publisher** → **Select your publisher**: **GitHub Actions**, with
  these exact, case-sensitive fields:

  | Field | Value |
  |---|---|
  | Organization or user | `theislab` |
  | Repository | `cellucid` |
  | Workflow filename | `npm-publish.yml` |
  | Environment name | (empty) |

  npm requires a new trusted-publisher configuration to complete its first
  publish within 2 days, or it expires: set it up shortly before the first
  release that should reach npm. Afterwards, **Settings** → **Publishing
  access** → **Require two-factor authentication and disallow tokens**
  refuses any token-based publish (`cellucid/npm/publishing.md`).
