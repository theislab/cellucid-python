"""The gene sample, as the notebook side reads it.

A large dataset stored cell by cell opens every cell on the map in the
viewer and reads a gene for a sample of the open cells (the first of them in
their order). Two records cross into Python: a session's ``genes`` (how many
open cells, every cell's order when every cell was open, their digest) and
the ``color_by`` event's ``cells_read``. Both are checked exactly as the web
app checks them (``app/session/cell-subset-record.js``
``decodeGeneSampleRecord``, ``data/viewer-host-source.js`` ``notifyColorBy``);
R checks them the same way (``tests/testthat/test-gene-sample.R``).
"""

from __future__ import annotations

import pytest

from cellucid.jupyter._wire import _require_inbound_jupyter_event
from cellucid.session_bundle import _decode_genes, _validate_manifest

GENES = {"cells": 40, "order": {"method": "blocks", "seed": 0, "blockSize": 19}, "digest": "0123456789abcdef"}


def test_a_session_gene_sample_is_read_as_recorded() -> None:
    assert _decode_genes(GENES, 100) == {"cells": 40, "order": {"method": "blocks", "seed": 0, "blockSize": 19}, "digest": "0123456789abcdef"}
    assert _decode_genes({**GENES, "order": None}, 100)["order"] is None


@pytest.mark.parametrize(
    ("genes", "message"),
    [
        ({**GENES, "cells": 100}, "genes.cells must be a count of open cells below"),
        ({**GENES, "cells": 0}, "genes.cells must be a count of open cells below"),
        ({**GENES, "cells": True}, "genes.cells must be a count of open cells below"),
        ({**GENES, "extra": 1}, "genes must have exactly cells, order and digest"),
        ({**GENES, "order": {"method": "uniform", "seed": 0, "blockSize": 19}}, "genes.order must be"),
        ({**GENES, "order": {"method": "blocks", "seed": -1, "blockSize": 19}}, "genes.order.seed must be a uint32"),
        ({**GENES, "order": {"method": "blocks", "seed": 0, "blockSize": 0}}, "genes.order.blockSize must be a positive count"),
        ({**GENES, "digest": "ABC"}, "genes.digest must be 16 lowercase hex digits"),
    ],
)
def test_a_session_gene_sample_out_of_its_contract_is_refused(genes, message) -> None:
    with pytest.raises(ValueError, match=message):
        _decode_genes(genes, 100)


def test_a_manifest_names_the_gene_sample_beside_its_other_keys() -> None:
    manifest = {
        "createdAt": "2026-10-28T12:00:00.000Z",
        "datasetFingerprint": {"cellCount": 100},
        "source": None,
        "requires": {"obs": [], "var": []},
        "chunks": [],
        "genes": GENES,
        "unexpected": 1,
    }
    # The gene sample is a known key; any other is still refused.
    with pytest.raises(ValueError, match="unknown unexpected"):
        _validate_manifest(manifest)


def test_color_by_says_the_cells_a_gene_is_read_for() -> None:
    event = {"type": "color_by", "field": "ENSG00000168484", "source": "var", "symbol": "SFTPC", "cells_read": 40_000}
    assert _require_inbound_jupyter_event(event)[1]["cells_read"] == 40_000
    assert _require_inbound_jupyter_event({**event, "cells_read": None})[1]["cells_read"] is None
    with pytest.raises(ValueError, match="cells read only for a gene"):
        _require_inbound_jupyter_event({"type": "color_by", "field": "leiden", "source": "obs", "symbol": None, "cells_read": 3})
    with pytest.raises(ValueError, match="cells_read"):
        _require_inbound_jupyter_event({**event, "cells_read": -1})
    with pytest.raises(ValueError, match="missing cells_read"):
        _require_inbound_jupyter_event({key: value for key, value in event.items() if key != "cells_read"})
