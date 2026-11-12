# Publishing Cellucid Python

The canonical production procedure is
[`docs/user_guide/python_package/h_developer_docs/14_release_process.md`](../../docs/user_guide/python_package/h_developer_docs/14_release_process.md).

## PyPI

To release, publish a GitHub Release named `vX.Y.Z` in
`theislab/cellucid-python`; the workflow `.github/workflows/pypi-publish.yml`
("Publish to PyPI") publishes it to PyPI. It runs on `release: published`
only: pushing a tag publishes nothing, and a failed run is retried with
**Re-run jobs** on that run. The workflow requires:

- a published (not draft) release, not marked pre-release, naming the tag;
- exact tag/package/release-metadata agreement;
- proof that the exact tagged commit belongs to `origin/main`;
- Ruff, mypy, pytest, and strict Sphinx gates;
- deterministic source-distribution normalization and exact recipe-digest
  verification;
- strict wheel and source-distribution metadata checks;
- installation of both artifacts on Linux, macOS, and Windows;
- the `pypi` GitHub environment; and
- PyPI trusted-publisher identity through GitHub OIDC.

The repository has no tag-push trigger, dispatch event, token secret, or local
upload command. The order across the three Cellucid packages and the shared
checklist are in `docs/contributing/releasing.md`.

For release 0.9.8, validate the candidate with: <!-- CELLUCID_VERSION -->

```bash
python scripts/validate_release.py --tag v0.9.8  # CELLUCID_VERSION
python -m ruff check src tests scripts/validate_release.py scripts/normalize_sdist.py
python -m mypy src/cellucid
python -m pytest
python -m sphinx -W --keep-going -b html docs docs/_build/html
python -m pip install "build==1.5.0" "twine==6.2.0"
SOURCE_DATE_EPOCH=315532800 python -m build
python scripts/normalize_sdist.py dist
python scripts/validate_release.py --tag v0.9.8 --sdist dist  # CELLUCID_VERSION
python -m twine check --strict dist/*
```

## Downstream conda recipe

`scripts/publishing/meta.yaml` records the exact source release and runtime
requirements for a conda-forge recipe submission. Before opening a downstream
pull request:

1. verify the version and runtime requirements equal the Python package;
2. download the source distribution PyPI published for that version and
   record its SHA-256 digest as `source.sha256`;
3. run the recipe import, `pip check`, version, and CLI tests; and
4. submit the recipe through the conda-forge review process.

The digest is recorded from the published file, never checked before it: a
release run that required the candidate archive to match it failed whenever
the source had changed since it was last recorded (0.9.6's first PyPI run),
and nothing before the release could notice. The release run checks only the
recipe's identity and requirements.

The recipe is downstream metadata. It is never sent to PyPI and does not change
the Python wheel.

## Read the Docs

Branch and pull-request documentation is validated by `docs-check.yml`.
`readthedocs.yml` requests the hosted `latest` documentation build on every
push to `main`; a tag's own documentation version is built by Read the
Docs' project settings (its version automation), not by this workflow. A release is
verified only after the tagged documentation build and stable link both
succeed.
