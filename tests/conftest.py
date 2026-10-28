"""Configuration every test in this suite runs under."""

from __future__ import annotations

from collections.abc import Iterator

import pytest


@pytest.fixture(autouse=True, scope="session")
def _mudata_keeps_modality_columns_in_the_modality() -> Iterator[None]:
    """Hold every MuData the suite builds or reads to one behaviour.

    The ``mudata`` extra admits mudata 0.3 and 0.4, and pip installs 0.3 on
    Python 3.11 (mudata 0.4 needs Python 3.12). mudata 0.4 keeps each
    modality's obs and var columns in the modality. mudata 0.3, unless told
    otherwise, copies them into the global obs and var whenever a MuData is
    built or read, and warns ``FutureWarning`` each time that this is about to
    change. Setting ``pull_on_update=False``, the switch mudata 0.3 names in
    that warning, gives 0.3 the 0.4 behaviour, so a read-back is compared with
    one expectation and checked for warnings on every supported mudata.
    mudata 0.4 has that behaviour only, and no ``set_options``.
    """
    try:
        import mudata
    except ImportError:
        yield
        return
    set_options = getattr(mudata, "set_options", None)
    if set_options is None:
        yield
        return
    with set_options(pull_on_update=False):
        yield


@pytest.fixture(autouse=True)
def _a_scratch_folder_is_removed_with_nothing_mapped(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    """Fail any test in which a served store removes its scratch folder while
    a file in it is still memory-mapped.

    Windows removes no mapped file, so there the folder would stay behind;
    other systems remove it without complaint. Checked when the folder is
    removed (a store's close swallows its closers' errors, so it is recorded
    and asserted after the test), it fails on every system where Windows would.
    """
    import shutil
    from pathlib import Path

    from _mapped_files import mapped_files

    from cellucid._served_store import _SCRATCH_PREFIX

    remove = shutil.rmtree
    still_mapped: list[Path] = []

    def checked(path, *args, **kwargs):  # type: ignore[no-untyped-def]
        if Path(path).name.startswith(_SCRATCH_PREFIX):
            still_mapped.extend(mapped_files(Path(path)))
        return remove(path, *args, **kwargs)

    monkeypatch.setattr(shutil, "rmtree", checked)
    yield
    assert still_mapped == [], f"removed while still memory-mapped: {still_mapped}"
