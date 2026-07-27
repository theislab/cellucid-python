# Release process

This is the exact production release process for the `cellucid` Python package.
Pushing a version tag starts PyPI publication, so the tag is created only after
the release candidate is reviewed and publication is explicitly authorized.

## Versioning policy

`cellucid` is currently in beta (`0.9.1`). <!-- CELLUCID_VERSION -->

- Package versions use three numeric components.
- Git tags use the exact package version with a `v` prefix, such as `v0.9.1`.
- PyPI distributions, citation metadata, documentation, and the downstream
  recipe declare the same version.
- The normalized source distribution has one reproducible byte stream, and its
  SHA-256 digest must equal `scripts/publishing/meta.yaml`.
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
to version tags and authorized maintainers. The workflow requests
`id-token: write` only in its final publication job.

```{warning}
Release preflight on 2026-07-26 queried the repository through the read-only
GitHub API and found **zero configured GitHub environments**. This does not make
the workflow mechanically unrunnable: GitHub documents that a workflow which
references a nonexistent environment creates it. That implicit environment
would not carry the reviewed release protections required by this project.

Publication is therefore blocked by release-security policy until an authorized
repository administrator explicitly creates and protects `pypi` *before* any
release tag is pushed. An authorized PyPI project administrator must also
confirm or configure the exact trusted-publisher record above, including the
`pypi` environment claim. Do not create a release tag until both external
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

Run the repository validator after all release metadata except the source
digest is aligned:

```bash
python scripts/validate_release.py --tag v0.9.1
```

The command must finish with:

```text
Cellucid Python release contract is valid for v0.9.1.
```

### 2. Run source and documentation gates

From the repository root:

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
python -c "import hashlib, pathlib; p=next(pathlib.Path('dist').glob('cellucid-*.tar.gz')); print(hashlib.sha256(p.read_bytes()).hexdigest())"
python scripts/validate_release.py --tag v0.9.1 --sdist dist
python -m twine check --strict dist/*
```

Record the printed SHA-256 digest as `source.sha256` in
`scripts/publishing/meta.yaml`, then rerun the validator with `--sdist`.
`MANIFEST.in` excludes the downstream recipe and repository tests from the
public source distribution, so recording the digest does not alter the
candidate archive. `normalize_sdist.py` canonicalizes member ordering,
timestamps, owners, modes, and the gzip byte stream; rebuilding and normalizing
the same source with the pinned setuptools 83.0.0 backend must produce the same
digest. Normalization rejects a source distribution larger than one mebibyte,
so the deterministic stored-block gzip representation has an explicit network
cost ceiling.
The fixed build epoch controls archive metadata only; it is not a release date
or scientific provenance timestamp.

Install the wheel and source distribution into separate clean environments.
For each artifact, verify:

```bash
python -m pip check
python -c "import cellucid; assert cellucid.__version__ == '0.9.1'"
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
finish successfully. Then create and push the exact annotated tag:

```bash
git tag -a v0.9.1 -m "Cellucid Python 0.9.1"
git push origin v0.9.1
```

The `Publish to PyPI` workflow:

1. checks out complete Git history, validates `v0.9.1`, and requires the tagged
   `GITHUB_SHA` to be an ancestor of the fetched `origin/main`;
2. runs Ruff, mypy, the complete pytest suite, and a strict documentation build;
3. builds one wheel and one source distribution;
4. normalizes the source distribution and requires its SHA-256 digest to equal
   the downstream recipe;
5. runs strict metadata checks;
6. installs both artifacts on Linux, macOS, and Windows with Python 3.14;
7. publishes those exact gated files through the `pypi` environment and OIDC.

There is no production dispatch button and no token-based publication path in
the repository.

## Verify the published release

After the workflow succeeds, create a clean environment and install from PyPI:

```bash
python -m venv cellucid-release-check
source cellucid-release-check/bin/activate
python -m pip install --upgrade pip
python -m pip install "cellucid==0.9.1"
python -m pip check
python -c "import cellucid; assert cellucid.__version__ == '0.9.1'"
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
- If publication has not started, correct the release commit and use the
  intended tag only after every gate is green.
- If PyPI accepted any file, do not replace it. Diagnose the cause, increment
  the package version, rebuild, and run the complete process again.
