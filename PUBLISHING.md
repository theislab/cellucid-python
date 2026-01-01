# Publishing Cellucid (Beginner Guide)

This guide walks you through publishing the `cellucid` Python package to:
- PyPI (pip install)
- conda-forge (conda install)
- Bioconda (bioinformatics-focused conda channel)
- Read the Docs (documentation hosting)

## 0) What gets published where

- **PyPI**: the canonical Python release artifact (wheel + sdist). Most downstream ecosystems (including conda-forge) pull from PyPI source releases.
- **conda-forge**: a separate “feedstock” repo maintained by conda-forge; it builds conda packages from your PyPI sdist.
- **Bioconda**: a separate `bioconda-recipes` PR; Bioconda packages typically depend on conda-forge for most Python dependencies.
- **Read the Docs**: builds documentation from your repo (usually from a Git tag or branch).

## 0.1) Recommended release model (lowest friction)

This repo already includes GitHub Actions workflows that publish when you push a Git tag:
- **PyPI publish**: `.github/workflows/pypi-publish.yml` (runs on tags like `v0.1.0`)
- **RTD trigger**: `.github/workflows/readthedocs.yml` (runs on `main` and `v*` tags)

If you use this model, you typically do **not** run `twine upload` locally.

## 1) One-time setup

### 1.1 Accounts + access

- Create a **PyPI** account: https://pypi.org/account/register/
- Enable **2FA** on PyPI (recommended/commonly required).
- (Optional but recommended) Create a **TestPyPI** account: https://test.pypi.org/account/register/
- Ensure you have GitHub permissions to create releases/tags in the `cellucid-python` repo.
  - If you will use GitHub Actions for PyPI: you also need access to set repo secrets.

### 1.2 Local tools

Create a clean environment (recommended):

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -U build twine
```

### 1.3 GitHub repo secrets (for automated publishing)

If you want publishing to happen automatically on tag push:

1) Create a PyPI API token (PyPI → Account settings → API tokens).
2) In GitHub repo settings → **Secrets and variables** → **Actions**:
   - Add `PYPI_API_TOKEN`
   - Add `READTHEDOCS_TOKEN` (only if you want RTD triggered from GitHub Actions)

## 2) Release checklist (do this every time)

### 2.1 Decide the version

Update `cellucid-python/pyproject.toml`:
- `version = "..."` (PEP 440 format)

Examples:
- stable: `0.1.0`
- prerelease: `0.1.0a1`, `0.1.0b1`, `0.1.0rc1`

### 2.2 Update changelog

Update `cellucid-python/CHANGELOG.md` with:
- what changed
- any breaking changes
- migration notes (if needed)

### 2.3 Run tests locally

From `cellucid-python/`:

```bash
pytest
```

### 2.3.1 Commit and tag (if using GitHub Actions)

From the `cellucid-python/` repo:

```bash
git status
git add -A
git commit -m "Release v<VERSION>"
git tag v<VERSION>
git push origin main
git push origin v<VERSION>
```

### 2.4 Build the artifacts (wheel + sdist)

From `cellucid-python/`:

```bash
rm -rf dist/
python -m build
```

You should now have:
- `dist/cellucid-<version>-py3-none-any.whl`
- `dist/cellucid-<version>.tar.gz`

### 2.5 Sanity-check the built artifacts (recommended)

```bash
python -m twine check dist/*
```

## 3) Publish to PyPI (pip)

### 3.0 Publish via GitHub Actions (recommended)

If `PYPI_API_TOKEN` is configured (Section 1.3):
1) Push a tag like `v0.0.9` (Section 2.3.1). <!-- CELLUCID_VERSION -->
2) Watch the GitHub Actions run: **Actions → Publish to PyPI**.

That workflow builds and uploads both wheel and sdist.

### 3.1 Publish to TestPyPI first (recommended)

```bash
python -m twine upload --repository testpypi dist/*
```

Then test install in a fresh env:

```bash
python -m venv /tmp/cellucid-test
source /tmp/cellucid-test/bin/activate
python -m pip install -U pip
python -m pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ cellucid
python -c "import cellucid; print(cellucid.__version__)"
```

### 3.2 Publish to real PyPI

```bash
python -m twine upload dist/*
```

## 4) Publish to conda-forge (conda)

conda-forge publishing happens in a separate repo called a **feedstock**.

### 4.1 First-time: create a conda-forge feedstock

1) Make sure the new version is on PyPI (Section 3).
2) Open a PR to conda-forge’s staged-recipes:
   - https://github.com/conda-forge/staged-recipes

In that PR you add a `recipe/meta.yaml` for `cellucid`.

Common ways to generate a starting recipe:
- `grayskull pypi cellucid` (then edit)
- manually write `meta.yaml` (fine for simple Python packages)

Minimal recipe ingredients you’ll need:
- `package: name/version`
- `source: url` pointing to the PyPI sdist + its `sha256`
- `build: noarch: python`
- `requirements: host/run` matching `pyproject.toml` dependencies
- `test: imports` and/or `pytest` invocation
- `about: license`, `license_file`, `home`, `summary`

**Important (common conda-forge failure):** conda-forge build jobs typically run with **no network access**.

If your recipe uses plain `pip install .` under PEP 517 build isolation, pip may try to download build requirements from PyPI and fail.

In conda recipes, prefer:

```yaml
script: {{ PYTHON }} -m pip install . --no-deps --no-build-isolation -vv
```

and list build requirements in `requirements: host:`.

### 4.1.1 How to get the sdist URL + sha256

In the conda recipe you need the exact PyPI source tarball URL and its SHA-256 hash.

1) Download the sdist:

```bash
python -m pip download --no-binary :all: --no-deps cellucid==<VERSION>
```

2) Compute sha256:

```bash
shasum -a 256 cellucid-<VERSION>.tar.gz
```

3) Use that URL + sha256 in `meta.yaml`.

### 4.1.2 A minimal `meta.yaml` skeleton (for reference)

This is a *starting point* (you’ll likely tweak deps/tests):

```yaml
{% set name = "cellucid" %}
{% set version = "<VERSION>" %}

package:
  name: {{ name|lower }}
  version: {{ version }}

source:
  url: https://pypi.io/packages/source/{{ name[0] }}/{{ name }}/{{ name }}-{{ version }}.tar.gz
  sha256: <SHA256>

build:
  noarch: python
  script: {{ PYTHON }} -m pip install . --no-deps --no-build-isolation -vv

requirements:
  host:
    - python >=3.10
    - pip
    - setuptools
  run:
    - python >=3.10
    - numpy >=1.21
    - pandas >=1.4
    - scipy >=1.7
    - tqdm >=4.45
    - anndata >=0.8
    - ipython >=7.23
    - jupyter-server-proxy >=4.1

test:
  imports:
    - cellucid

about:
  home: https://github.com/theislab/cellucid-python
  license: BSD-3-Clause
  license_file: LICENSE
  summary: Interactive Single-Cell Data Visualization
```

### 4.1.3 Common conda-forge CI failure modes (and fixes)

- **PEP 517 build isolation tries to download from PyPI**: use `--no-build-isolation` and list build deps under `requirements: host:`.
- **Tests accidentally start the server / open a browser**: keep conda tests to `python -c "import cellucid"` or a small unit test set; don’t run `cellucid serve` as a recipe test.
- **Prerelease versions (`a1`, `b1`, `rc1`)**: if reviewers push back, publish a stable `0.x.y` first.

After the staged-recipes PR is merged:
- conda-forge creates `cellucid-feedstock`
- CI builds and uploads the package to conda-forge

### 4.2 Updating versions after the feedstock exists

After the feedstock exists, updates are usually automatic:
- conda-forge’s **autotick bot** opens PRs when it detects new PyPI releases.
- you (or maintainers) review/merge the bot PR.

## 5) Publish to Bioconda

Bioconda is a separate ecosystem:
- https://github.com/bioconda/bioconda-recipes

Typical workflow:
1) Ensure the version is on PyPI (Section 3).
2) Ensure runtime dependencies exist on conda-forge/bioconda.
3) Add or update a recipe in `bioconda-recipes` via PR.

For most Python tools, the pragmatic path is:
1) Get the package onto conda-forge first (Section 4).
2) Then make a Bioconda recipe that depends on the conda-forge package.

Bioconda has stricter policies and review expectations than many projects:
- correct license metadata
- correct dependency declarations
- tests that don’t require network access

Note: in many cases you should publish to **conda-forge first**, then Bioconda can depend on the conda-forge package.

## 6) Publish docs on Read the Docs (RTD)

If RTD is already working, you typically only need to:
- ensure the docs build for the new tag/version
- keep dependencies pinned in the RTD config

Common workflow:
1) Create a Git tag (e.g. `v0.1.0`) and push it.
2) In RTD project settings, enable building tags (if desired).
3) Confirm the build succeeds for the new tag.

If you use the GitHub Actions trigger (`.github/workflows/readthedocs.yml`), pushing a tag is enough.

## 7) Troubleshooting (high-signal)

### “Build fails with pyproject.toml validation errors”

- Fix invalid PEP 621 fields in `cellucid-python/pyproject.toml`.
- Re-run `python -m build`.

### “conda-forge/bioconda builds fail but pip works”

- conda builds from sdists in clean environments with stricter dependency resolution.
- ensure your dependencies are correctly declared and available in conda.
- ensure tests do not rely on network.

### “RTD build fails but local docs build works”

- RTD uses a clean environment; add missing doc deps to the RTD config / docs extras.
- pin incompatible Sphinx extensions.

## 8) Quick reference (minimal happy path)

From `cellucid-python/`:

```bash
pytest
python -m build
python -m twine check dist/*
python -m twine upload dist/*
```
