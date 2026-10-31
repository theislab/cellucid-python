"""A served matrix's entries asked for while it is laid out gene-major.

The served store lays a cell-major expression matrix out gene-major in the
background (:class:`~cellucid.prepare_data._gene_major.GeneMajorLayout`, the
store's worker). A range asked for before the placing pass has reached it is
made from that pass's progress (:meth:`GeneMajorLayout.entries_now`): the
entries already placed are copied and only those of the rows not reached yet
are made. These tests stand the pass halfway and check that each range is
exactly what the pass places, and count the rows a request passed over and the
entries it copied -- never the time it took. cellucid-r does the same in C
(``tests/testthat/test-served-layout-progress.R``).
"""

from __future__ import annotations

import threading
from typing import Any

import numpy as np
import pytest
from scipy import sparse

from cellucid._source_arrays import sparse_source
from cellucid.prepare_data import _gene_major
from cellucid.prepare_data._gene_major import GeneMajorLayout


def _matrix(n_obs: int, n_vars: int, per_row: int, *, seed: int, canonical: bool = True) -> sparse.csr_matrix:
    """A cells x genes CSR matrix, ``per_row`` values a row; not canonical:
    each row's columns shuffled and one column stored twice."""
    rng = np.random.default_rng(seed)
    columns = np.stack([rng.choice(n_vars, per_row, replace=False) for _ in range(n_obs)])
    if canonical:
        columns.sort(axis=1)
    else:
        columns[:, -1] = columns[:, 0]
    values = rng.integers(1, 50, n_obs * per_row).astype(np.float32)
    indptr = np.arange(0, n_obs * per_row + 1, per_row, dtype=np.int64)
    matrix = sparse.csr_matrix((values, columns.reshape(-1).astype(np.int32), indptr), shape=(n_obs, n_vars))
    if not canonical:
        matrix.has_sorted_indices = False
        matrix.has_canonical_format = False
    return matrix


def _layout(matrix: sparse.csr_matrix, *, as_stored: bool, genes: np.ndarray | None = None, seed: int = 7) -> GeneMajorLayout:
    order = np.random.default_rng(seed).permutation(matrix.shape[0]).astype(np.int64)
    new_row = np.empty_like(order)
    new_row[order] = np.arange(order.size, dtype=np.int64)
    layout = GeneMajorLayout(sparse_source(matrix, None), order=order, new_row=new_row, label="X", genes=genes, as_stored=as_stored)
    layout.count()
    return layout


class _Halfway:
    """The layout's placing pass on a thread, standing before its block ``stand_at`` until let go."""

    def __init__(self, layout: GeneMajorLayout, stand_at: int, monkeypatch: pytest.MonkeyPatch) -> None:
        self.reached = threading.Event()
        self.go = threading.Event()
        blocks = {"placed": 0}
        place_block = GeneMajorLayout._place_block

        def standing(self_: GeneMajorLayout, *args: Any, **kwargs: Any) -> np.ndarray:
            if self_ is layout:
                if blocks["placed"] == stand_at:
                    self.reached.set()
                    assert self.go.wait(30)
                blocks["placed"] += 1
            return place_block(self_, *args, **kwargs)

        monkeypatch.setattr(GeneMajorLayout, "_place_block", standing)
        self.thread = threading.Thread(target=layout.place)
        self.thread.start()
        assert self.reached.wait(30)

    def finish(self) -> None:
        self.go.set()
        self.thread.join(30)


def _ranges(indptr: np.ndarray) -> list[tuple[int, int]]:
    nnz = int(indptr[-1])
    return [(0, nnz), (0, 1), (nnz - 1, nnz), (17, min(nnz, 950)), (int(indptr[3]), int(indptr[5])), (int(indptr[2]) + 3, int(indptr[6]) - 2)]


@pytest.mark.parametrize(
    ("as_stored", "canonical", "kept"),
    [
        (False, True, False),  # an object's snapshot, canonical: the scan
        (True, False, False),  # a file's matrix as stored, rows unsorted, a column twice: the scan
        (False, True, True),  # some genes kept: the scan through the column map
        (False, False, False),  # a snapshot the layout makes canonical: its rows not reached, line by line
    ],
)
def test_a_range_asked_halfway_is_made_from_the_progress_passing_over_only_the_rows_left(
    monkeypatch: pytest.MonkeyPatch, as_stored: bool, canonical: bool, kept: bool
) -> None:
    monkeypatch.setattr(_gene_major, "PASS_ENTRIES", 400)
    matrix = _matrix(1200, 31, 6, seed=1, canonical=canonical)
    genes = np.array([30, 2, 17, 3, 9, 11, 25], dtype=np.int64) if kept else None
    expected = _layout(matrix, as_stored=as_stored, genes=genes).place()
    layout = _layout(matrix, as_stored=as_stored, genes=genes)
    halfway = _Halfway(layout, stand_at=5, monkeypatch=monkeypatch)
    try:
        rows = layout._rows_placed
        assert 0 < rows < matrix.shape[0]
        passed = copied = 0
        for start, end in _ranges(expected.indptr):
            indices, data = layout.entries_now(start, end)
            np.testing.assert_array_equal(indices, expected.indices[start:end])
            np.testing.assert_array_equal(data, expected.data[start:end])
            passed += matrix.shape[0] - rows
            copied += int(np.count_nonzero(expected.indices[start:end] < rows))
            assert (layout.request_rows, layout.request_copied) == (passed, copied), (start, end)
    finally:
        halfway.finish()
    assert layout.result is not None
    np.testing.assert_array_equal(layout.result.indices, expected.indices)
    np.testing.assert_array_equal(layout.result.data, expected.data)
    # Laid out: a range is the pass's own output, nothing passed over.
    indices, _ = layout.entries_now(5, 800)
    np.testing.assert_array_equal(indices, expected.indices[5:800])
    assert (layout.request_rows, layout.request_copied) == (passed, copied)


def test_a_range_asked_before_the_pass_begins_passes_over_every_row() -> None:
    matrix = _matrix(500, 13, 4, seed=2)
    expected = _layout(matrix, as_stored=False).place()
    layout = _layout(matrix, as_stored=False)
    indices, data = layout.entries_now(10, 300)
    np.testing.assert_array_equal(indices, expected.indices[10:300])
    np.testing.assert_array_equal(data, expected.data[10:300])
    assert (layout.request_rows, layout.request_copied) == (500, 0)


def test_the_pass_is_waited_for_once_requests_have_spent_what_it_has_left(monkeypatch: pytest.MonkeyPatch) -> None:
    layout = _layout(_matrix(1000, 13, 4, seed=4), as_stored=False)
    clock = {"now": 100.0}
    monkeypatch.setattr(_gene_major.time, "monotonic", lambda: clock["now"])
    assert not layout.sooner()  # not begun
    # A quarter of the rows placed in 1 s: 3 s left.
    layout._placing_began, layout._rows_placed = 99.0, 250
    assert not layout.sooner()
    layout._request_last, layout._request_spent = 0.5, 2.9
    assert not layout.sooner()
    # Requests have spent as long as the pass has left: waited for.
    layout._request_spent = 3.0
    assert layout.sooner()
    # Or the last request alone took as long.
    layout._request_last, layout._request_spent = 3.2, 3.2
    assert layout.sooner()
    layout._request_last, layout._request_spent = 0.1, 0.1
    clock["now"] = 100.0 + 9.0  # 10 s for a quarter: 30 s left
    assert not layout.sooner()


def test_every_pass_steps_by_one_size_taken_from_the_machine_s_memory() -> None:
    from cellucid import _served_store
    from cellucid._scratch import machine_memory_bytes

    # One definition: the served store's passes step by the layout's.
    assert _served_store.PASS_ENTRIES == _gene_major.PASS_ENTRIES
    assert not hasattr(_served_store, "_PASS_ENTRIES")
    assert max(1, machine_memory_bytes() // 64 // _gene_major._BLOCK_BYTES_PER_VALUE) == _gene_major.PASS_ENTRIES
    assert _served_store.default_cache_bytes() == machine_memory_bytes() // 64


@pytest.mark.parametrize("as_stored", [False, True])
def test_the_layout_is_the_same_bytes_whatever_a_pass_steps_by(monkeypatch: pytest.MonkeyPatch, as_stored: bool) -> None:
    matrix = _matrix(900, 23, 5, seed=6, canonical=not as_stored)
    results = []
    for step in (37, 1000, _gene_major.PASS_ENTRIES):
        monkeypatch.setattr(_gene_major, "PASS_ENTRIES", step)
        results.append(_layout(matrix, as_stored=as_stored).place())
    for result in results[1:]:
        np.testing.assert_array_equal(result.indptr, results[0].indptr)
        np.testing.assert_array_equal(result.indices, results[0].indices)
        np.testing.assert_array_equal(result.data, results[0].data)
