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

This section provides a complete, step-by-step guide for publishing `cellucid` to conda-forge.

### 4.0 Understanding conda-forge

**What is conda-forge?**
- conda-forge is a community-maintained collection of conda packages
- It's a separate ecosystem from PyPI - packages need their own "recipe"
- Once set up, updates happen mostly automatically

**How it works:**
1. You submit a **recipe** (`meta.yaml`) to conda-forge's `staged-recipes` repo
2. After review and merge, conda-forge creates a **feedstock** repo for your package
3. The feedstock automatically builds and publishes your package
4. Future updates are handled by a bot that detects new PyPI releases

**Prerequisites:**
- Your package must already be published on PyPI (Section 3)
- You need a GitHub account
- Basic familiarity with Git and pull requests

---

### 4.1 Step-by-step: First-time submission

#### Step 1: Ensure your package is on PyPI

Before starting, verify your package is available on PyPI:

```bash
pip install cellucid  # Should work
```

#### Step 2: Generate the recipe with grayskull

We've already generated a `meta.yaml` for you. It's located at:

```
scripts/publishing/meta.yaml
```

If you need to regenerate it (e.g., for a new version), run:

```bash
# Install grayskull if you haven't
pip install grayskull

# Generate recipe (run from any directory)
grayskull pypi cellucid --strict-conda-forge
```

This creates a `cellucid/meta.yaml` file with the correct structure.

#### Step 3: Fork the staged-recipes repository

1. Go to: https://github.com/conda-forge/staged-recipes
2. Click **Fork** (top-right corner)
3. This creates a copy at `https://github.com/YOUR-USERNAME/staged-recipes`

#### Step 4: Clone your fork locally

```bash
git clone https://github.com/YOUR-USERNAME/staged-recipes.git
cd staged-recipes
```

#### Step 5: Create a new branch

```bash
git checkout -b add-cellucid
```

#### Step 6: Add your recipe

Create the recipe directory and copy the meta.yaml:

```bash
# Create the recipe folder
mkdir -p recipes/cellucid

# Copy the meta.yaml (adjust the source path as needed)
cp /path/to/cellucid-python/scripts/publishing/meta.yaml recipes/cellucid/meta.yaml
```

Your directory structure should look like:

```
staged-recipes/
└── recipes/
    └── cellucid/
        └── meta.yaml
```

#### Step 7: Validate the recipe locally (optional but recommended)

```bash
# Install conda-build if you haven't
conda install conda-build

# Lint the recipe
conda-build --check -c conda-forge recipes/cellucid

# Optionally, try a full build (takes a while)
conda-build -c conda-forge recipes/cellucid
```

#### Step 8: Commit and push

```bash
git add recipes/cellucid/meta.yaml
git commit -m "Add recipe for cellucid"
git push origin add-cellucid
```

#### Step 9: Create the pull request

1. Go to your fork: `https://github.com/YOUR-USERNAME/staged-recipes`
2. You'll see a banner: "Compare & pull request" - click it
3. Fill in the PR template:
   - Title: `Add cellucid`
   - Description: Brief description of what cellucid does
4. Click **Create pull request**

#### Step 10: Wait for CI and review

- conda-forge CI will automatically build your package on Linux, macOS, and Windows
- A conda-forge maintainer will review your recipe
- Address any feedback by pushing more commits to your branch
- Once approved and merged, your feedstock is created!

---

### 4.2 The meta.yaml explained

Here's what each section of the `meta.yaml` does:

```yaml
{% set name = "cellucid" %}
{% set version = "0.0.9" %}
```
**Variables:** Jinja2 variables for reuse throughout the recipe.

```yaml
package:
  name: {{ name|lower }}
  version: {{ version }}
```
**Package info:** The conda package name and version.

```yaml
source:
  url: https://pypi.io/packages/source/c/cellucid/cellucid-0.0.9.tar.gz
  sha256: 486d0bdd015f8a77a92c76f8fd28a5fe105f1f542366c1f337d6c3a1053567fc
```
**Source:** Where to download the package. The SHA256 ensures integrity.

```yaml
build:
  entry_points:
    - cellucid = cellucid.cli:main
  noarch: python
  script: {{ PYTHON }} -m pip install . -vv --no-deps --no-build-isolation
  number: 0
```
**Build settings:**
- `entry_points`: CLI commands your package provides
- `noarch: python`: Package works on all platforms (pure Python)
- `script`: How to install (pip install with no network access flags)
- `number`: Build number (increment for same version rebuilds)

```yaml
requirements:
  host:
    - python >=3.10
    - setuptools >=68
    - pip
  run:
    - python >=3.10
    - numpy >=1.21
    # ... other dependencies
```
**Dependencies:**
- `host`: Needed at build time (includes build backend)
- `run`: Needed at runtime (your actual dependencies)

```yaml
test:
  imports:
    - cellucid
  commands:
    - pip check
    - cellucid --help
  requires:
    - pip
```
**Tests:** Verifies the package works after installation.

```yaml
about:
  home: https://github.com/theislab/cellucid-python
  summary: Interactive Single-Cell Data Visualization
  license: BSD-3-Clause
  license_file: LICENSE
```
**Metadata:** Package information displayed on conda-forge.

```yaml
extra:
  recipe-maintainers:
    - inecik
```
**Maintainers:** GitHub usernames of recipe maintainers (that's you!).

---

### 4.3 Updating the SHA256 for new versions

When you release a new version, you need to update the `version` and `sha256` in `meta.yaml`.

**Method 1: Use grayskull (easiest)**

```bash
grayskull pypi cellucid --strict-conda-forge
# Copy the sha256 from the generated meta.yaml
```

**Method 2: Calculate manually**

```bash
# Download the sdist from PyPI
pip download --no-binary :all: --no-deps cellucid==NEW_VERSION

# Calculate SHA256
shasum -a 256 cellucid-NEW_VERSION.tar.gz
# or on Linux:
sha256sum cellucid-NEW_VERSION.tar.gz
```

**Method 3: Get from PyPI website**

1. Go to: https://pypi.org/project/cellucid/#files
2. Click on the `.tar.gz` file
3. Click "view hashes" → copy the SHA256

---

### 4.4 After feedstock creation: Updating versions

Once your feedstock exists (e.g., `conda-forge/cellucid-feedstock`), updates are mostly automatic:

1. **Automatic (recommended):**
   - conda-forge's **regro-cf-autotick-bot** monitors PyPI
   - When you release a new version on PyPI, the bot opens a PR
   - Review the PR and merge it

2. **Manual updates:**
   - Fork the feedstock: `https://github.com/conda-forge/cellucid-feedstock`
   - Edit `recipe/meta.yaml` with new version and sha256
   - Open a PR to the feedstock

---

### 4.5 Common issues and solutions

| Issue | Solution |
|-------|----------|
| **Build fails with network error** | Ensure `--no-build-isolation` is in the script and all build deps are in `host:` |
| **"Package not found" for a dependency** | Check if the dependency exists on conda-forge. If not, it needs to be added first |
| **Tests fail with import error** | Check that all runtime dependencies are listed in `run:` |
| **CI timeout** | Simplify tests; avoid running the server or heavy computations |
| **Prerelease version rejected** | conda-forge prefers stable versions; publish `0.x.y` instead of `0.x.ya1` |
| **License not detected** | Ensure `LICENSE` file exists in your repo root |

---

### 4.6 Useful links

- **staged-recipes repo:** https://github.com/conda-forge/staged-recipes
- **conda-forge docs:** https://conda-forge.org/docs/maintainer/adding_pkgs.html
- **grayskull docs:** https://github.com/conda/grayskull
- **Recipe linter:** https://github.com/conda-forge/conda-smithy
- **Your feedstock (after creation):** https://github.com/conda-forge/cellucid-feedstock

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
