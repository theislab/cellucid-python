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
