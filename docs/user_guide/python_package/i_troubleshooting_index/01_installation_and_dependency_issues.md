# Installation and dependency issues

## `ModuleNotFoundError: No module named 'cellucid'`

Confirm that `pip` and Python refer to the same interpreter:

```bash
python -m pip show cellucid
python -c "import sys; print(sys.executable)"
python -c "import cellucid; print(cellucid.__version__)"
```

If `pip show` reports no installation, activate the intended environment and run
`python -m pip install cellucid`.

## Python version rejected

Cellucid supports Python 3.11 through 3.14. Check with `python --version`, then
create an environment in that range if necessary.

## CLI command not found

Run `python -m pip show cellucid`, then inspect the environment that owns
`python`. A CLI installed in another environment will not be on the active
environment's executable path.

For environment creation, editable installs, and the full verification sequence,
see {doc}`../a_landing_pages/02_installation`.
