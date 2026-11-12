# Release process

This is the exact production release process for the `cellucid` Python package.
To release it, publish a GitHub Release named `vX.Y.Z` in
`theislab/cellucid-python`; the workflow `pypi-publish.yml` ("Publish to PyPI")
checks it and publishes it to PyPI. Pushing a tag alone publishes nothing.
{doc}`/contributing/releasing` covers all three Cellucid packages in order.

## Versioning policy

`cellucid` is currently in beta (`0.9.8`). <!-- CELLUCID_VERSION -->

- Package versions use three numeric components.
- Git tags use the exact package version with a `v` prefix, such as `v0.9.8`. <!-- CELLUCID_VERSION -->
- PyPI distributions, citation metadata, documentation, and the downstream
  recipe declare the same version.
- The normalized source distribution has one reproducible byte stream. The
  downstream recipe's `source.sha256` is recorded from the file PyPI
  published, when the recipe goes to conda-forge
  (`scripts/publishing/PUBLISHING.md`); no release step compares it.
- Wheel ZIP timestamps use the fixed ZIP-compatible epoch
  `SOURCE_DATE_EPOCH=315532800` (1980-01-01T00:00:00Z), so two builds of the
  same source and exact backend produce identical wheel bytes.
- PyPI files are immutable. A changed release candidate receives a new version.

## One-time trusted-publisher setup

The workflow uses PyPI trusted publishing through GitHub OIDC. It does not read
a long-lived PyPI token.

In the PyPI `cellucid` project settings, register this publisher:

| Setting | Exact value |
|---|---|
| Owner | `theislab` |
| Repository | `cellucid-python` |
| Workflow | `pypi-publish.yml` |
| Environment | `pypi` |

In the GitHub repository, create the `pypi` environment and restrict deployment
to `v*` tags and authorized maintainers. The workflow requests
`id-token: write` only in its final publication job.

```{warning}
Release preflight on 2026-07-26 queried the repository through the read-only
GitHub API and found **zero configured GitHub environments**. This does not make
the workflow mechanically unrunnable: GitHub documents that a workflow which
references a nonexistent environment creates it. That implicit environment
would not carry the reviewed release protections required by this project.

Publication is therefore blocked by release-security policy until an authorized
repository administrator explicitly creates and protects `pypi` *before* any
GitHub Release is published. An authorized PyPI project administrator must also
confirm or configure the exact trusted-publisher record above, including the
`pypi` environment claim. Do not publish a GitHub Release until both external
settings have been verified; repository source validation cannot prove either
setting. See GitHub's
[environment behavior and protection documentation](https://docs.github.com/en/actions/how-tos/deploy/configure-and-manage-deployments/manage-environments)
and PyPI's
[trusted-publisher claim model](https://docs.pypi.org/trusted-publishers/internals/).
```

## Prepare the release candidate

### 1. Align release metadata

Update all release surfaces:

- `pyproject.toml`
- `CHANGELOG.md`
- `CITATION.cff`
- `docs/conf.py`
- the installation, packaging, and release pages
- `scripts/publishing/meta.yaml`

**Dates.** The current version's heading stays undated (`## [X.Y.Z]`) until
its GitHub Release is published, and `CITATION.cff` never carries
`date-released:`: nobody knows the publication day before it happens, so
`validate_release.py` refuses both. The bump that opens the next version dates
the released heading with the day its GitHub Release was published, read with
`gh release view vX.Y.Z --repo theislab/cellucid-python --json publishedAt`
(UTC), as `## [X.Y.Z] - YYYY-MM-DD`.

Every declaration a version bump has to touch carries a `CELLUCID_VERSION`
marker comment (`grep -rn CELLUCID_VERSION`). From a workspace that holds the
Cellucid repositories side by side, check that they all agree:

```bash
python cellucid-python/scripts/sweep_versions.py
```

Run the repository validator after all release metadata except the source
digest is aligned:

```bash
python scripts/validate_release.py --tag v0.9.8  # CELLUCID_VERSION
```

The command must finish with:

```text
Cellucid Python release contract is valid for v0.9.8.
```

### 2. Run source and documentation gates

From the repository root, in an environment with `.[dev,docs]` installed (the
`dev` extra brings `mudata`, which the MuData tests need):

```bash
python -m ruff check src tests scripts/validate_release.py scripts/normalize_sdist.py
python -m mypy src/cellucid
python -m pytest
python -m sphinx -W --keep-going -b html docs docs/_build/html
```

### 3. Build, normalize, and inspect both distributions

Use a clean environment with Python 3.14:

```bash
python -m pip install --upgrade pip
python -m pip install "build==1.5.0" "twine==6.2.0"
SOURCE_DATE_EPOCH=315532800 python -m build
python scripts/normalize_sdist.py dist
python scripts/validate_release.py --tag v0.9.8 --sdist dist  # CELLUCID_VERSION
python -m twine check --strict dist/*
```

The validator requires one source distribution named for the version that
ships neither the downstream recipe nor the tests. They are held out by two
different mechanisms. `scripts/` is not a package and no rule adds
it, so the downstream recipe never enters the archive's file set at all;
`MANIFEST.in` is the one line `prune tests`, and that line is what removes the
repository tests, which setuptools would otherwise carry. `normalize_sdist.py`
canonicalizes member ordering, timestamps, owners, modes, and the gzip byte
stream; rebuilding and normalizing the same source with the pinned setuptools
83.0.0 backend must produce the same digest. Normalization rejects a source
distribution larger than one mebibyte, so the deterministic stored-block gzip
representation has an explicit network cost ceiling.
The fixed build epoch controls archive metadata only; it is not a release date
or scientific provenance timestamp.

Install the wheel and source distribution into separate clean environments.
For each artifact, verify:

```bash
python -m pip check
python -c "import cellucid; assert cellucid.__version__ == '0.9.8'"  # CELLUCID_VERSION
cellucid --version
cellucid --help
cellucid serve --help
```

### 4. Review the exact candidate

Before publication authorization:

- inspect the complete diff;
- confirm the working tree contains no generated build output;
- confirm CI is green on the release commit;
- confirm wheel and source-distribution metadata declare the same dependencies;
- confirm `validate_release.py --sdist` proves that the downstream recipe
  records the exact normalized candidate hash.

## Publish after explicit authorization

Commit and push the reviewed source first. Wait for branch CI and docs checks to
finish successfully on `main`. Then publish a GitHub Release named `v0.9.8` <!-- CELLUCID_VERSION -->
in `theislab/cellucid-python`:

1. **Releases** → **Draft a new release**.
2. **Choose a tag** → type the tag, `vX.Y.Z` (the version above) → **Create new
   tag** on publish, target `main`.
3. Title `cellucid X.Y.Z`; paste the version's `CHANGELOG.md` section as the
   notes.
4. Leave **Set as a pre-release** unticked (the package declares `X.Y.Z`
   versions only, so a pre-release is refused) and press **Publish release**.

The `Publish to PyPI` workflow then:

1. checks out complete Git history, validates `v0.9.8` against the published <!-- CELLUCID_VERSION -->
   release (`validate_release.py --release-event`: the release is not a draft,
   is not marked pre-release, and names this tag), and requires the tagged
   `GITHUB_SHA` to be an ancestor of the fetched `origin/main`;
2. runs Ruff, mypy, the complete pytest suite, and a strict documentation build;
3. builds one wheel and one source distribution;
4. normalizes the source distribution and requires its SHA-256 digest to equal
   the downstream recipe;
5. runs strict metadata checks;
6. installs both artifacts on Linux, macOS, and Windows with Python 3.14;
7. publishes those exact gated files through the `pypi` environment and OIDC.

The workflow runs on `release: published` only. There is no tag-push trigger,
no dispatch button, and no token-based publication path in the repository.

## Verify the published release

After the workflow succeeds, create a clean environment and install from PyPI:

```bash
python -m venv cellucid-release-check
source cellucid-release-check/bin/activate
python -m pip install --upgrade pip
python -m pip install "cellucid==0.9.8"  # CELLUCID_VERSION
python -m pip check
python -c "import cellucid; assert cellucid.__version__ == '0.9.8'"  # CELLUCID_VERSION
cellucid --version
```

On Windows PowerShell, activate with:

```powershell
.\cellucid-release-check\Scripts\Activate.ps1
```

Verify that PyPI lists exactly one wheel and one source distribution and that
their hashes match the workflow artifact. Then verify the tagged Read the Docs
build and the stable documentation link.

## Failure handling

- A tag/version mismatch stops before build or publication.
- A tag whose commit is not contained in `origin/main` stops before the source
  or artifact gates start.
- A source-distribution digest mismatch stops before artifact upload.
- A failed source, documentation, metadata, or installed-artifact gate prevents
  the publication job from starting.
- A gate that failed for a passing reason (a runner or network failure) is
  retried with **Re-run jobs** (or **Re-run failed jobs**) on that run in the
  Actions tab; it runs the same release again.
- If publication has not started and the source must change, delete the GitHub
  Release and its tag, correct the release commit on `main`, and publish the
  release again once every gate is green.
- If PyPI accepted any file, do not replace it. Diagnose the cause, increment
  the package version, rebuild, and run the complete process again.
