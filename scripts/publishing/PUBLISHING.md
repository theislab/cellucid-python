# Publishing Cellucid Python

The canonical production procedure is
[`docs/user_guide/python_package/h_developer_docs/14_release_process.md`](../../docs/user_guide/python_package/h_developer_docs/14_release_process.md).

## PyPI

PyPI publication is performed only by
`.github/workflows/pypi-publish.yml` after an explicitly authorized version tag
is pushed. The workflow requires:

- exact tag/package/release-metadata agreement;
- proof that the exact tagged commit belongs to `origin/main`;
- Ruff, mypy, pytest, and strict Sphinx gates;
- deterministic source-distribution normalization and exact recipe-digest
  verification;
- strict wheel and source-distribution metadata checks;
- installation of both artifacts on Linux, macOS, and Windows;
- the `pypi` GitHub environment; and
- PyPI trusted-publisher identity through GitHub OIDC.

The repository has no production dispatch event, token secret, or local upload
command.

For release 0.9.1, validate the candidate with:

```bash
python scripts/validate_release.py --tag v0.9.1
python -m ruff check src tests scripts/validate_release.py scripts/normalize_sdist.py
python -m mypy src/cellucid
python -m pytest
python -m sphinx -W --keep-going -b html docs docs/_build/html
python -m pip install "build==1.5.0" "twine==6.2.0"
SOURCE_DATE_EPOCH=315532800 python -m build
python scripts/normalize_sdist.py dist
python scripts/validate_release.py --tag v0.9.1 --sdist dist
python -m twine check --strict dist/*
```

## Downstream conda recipe

`scripts/publishing/meta.yaml` records the exact source release and runtime
requirements for a conda-forge recipe submission. Before opening a downstream
pull request:

1. verify the version and runtime requirements equal the Python package;
2. require the normalized candidate digest to equal `source.sha256` before
   publication;
3. download the published source distribution;
4. require its SHA-256 digest to equal the already-gated `source.sha256`;
5. run the recipe import, `pip check`, version, and CLI tests; and
6. submit the recipe through the conda-forge review process.

The recipe is downstream metadata. It is never sent to PyPI and does not change
the Python wheel.

## Read the Docs

Branch and pull-request documentation is validated by `docs-check.yml`.
`readthedocs.yml` requests the hosted documentation build for the exact pushed
branch or tag. A release is verified only after the tagged documentation build
and stable link both succeed.
