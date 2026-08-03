"""A payload the server refuses has to say why, and the browser has to hear it.

The defect: ``AnnDataAdapter`` validates every gene column before sending a byte,
and ``_handle_var`` caught that specific ``ValueError`` alongside every other
exception and answered ``500 Internal server error`` with a twenty-one byte body.
Measured against a real 18,142,044-cell object, all five thousand genes answered
500 after about fifty seconds while continuous ``obs`` columns of identical size
answered 200 in three. The reason existed — it named the gene — and never left
the machine.

Four things are held here.

* The refusal counts. "One cell of eighteen million" and "every cell" were the
  same sentence before, and they want opposite responses.
* It is HTTP 422, not 500. The request was well formed and the server is
  working; the object it names cannot be drawn. A 500 sends the reader to the
  server logs to look for a crash that did not happen.
* The body is JSON the browser can present without parsing prose back apart.
* The same refusal covers continuous ``obs`` columns and the exporter, because it
  is the same defect wherever a continuous payload is published.
"""

from __future__ import annotations

import http.client
import json
from collections.abc import Iterator
from contextlib import contextmanager

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid.anndata_adapter import (
    GENE_EXPRESSION_CACHE_BYTES,
    AnnDataAdapter,
    LRUCache,
)
from cellucid.anndata_server import AnnDataServer
from cellucid.continuous_payload_diagnosis import (
    NON_FINITE_EXAMPLE_LIMIT,
    NonFinitePayloadError,
    diagnose_continuous_payload,
)


def _adata(
    *,
    gene_values: np.ndarray | None = None,
    obs_values: np.ndarray | None = None,
) -> ad.AnnData:
    cells = 4
    expression = (
        np.array([[1.0], [2.0], [3.0], [4.0]], dtype=np.float32)
        if gene_values is None
        else gene_values.reshape(cells, 1).astype(np.float32)
    )
    score = (
        np.array([0.25, 0.5, 0.75, 1.0], dtype=np.float32)
        if obs_values is None
        else obs_values.astype(np.float32)
    )
    adata = ad.AnnData(
        X=expression,
        obs=pd.DataFrame(
            {"score": score},
            index=pd.Index([f"cell-{i}" for i in range(cells)], dtype=object),
        ),
        var=pd.DataFrame(index=pd.Index(["Gene_A"], dtype=object)),
    )
    adata.obsm["X_umap_2d"] = np.arange(cells * 2, dtype=np.float32).reshape(
        cells, 2
    )
    return adata


def _adapter(adata: ad.AnnData) -> AnnDataAdapter:
    return AnnDataAdapter(
        adata,
        dataset_name="Diagnosis contract",
        dataset_id="diagnosis-contract",
    )


@contextmanager
def _running_server(adata: ad.AnnData) -> Iterator[tuple[str, int]]:
    server = AnnDataServer(
        adata,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Diagnosis contract",
        dataset_id="diagnosis-contract",
        serve_connectivity=False,
        serve_vector_fields=False,
    )
    server.start_background()
    assert server._server is not None
    host, port = server._server.server_address
    try:
        yield host, port
    finally:
        server.stop()


def _get(host: str, port: int, path: str) -> tuple[int, str, bytes]:
    connection = http.client.HTTPConnection(host, port, timeout=10)
    try:
        connection.request("GET", path)
        response = connection.getresponse()
        body = response.read()
        return response.status, response.getheader("Content-Type") or "", body
    finally:
        connection.close()


# ---------------------------------------------------------------------------
# The diagnosis itself
# ---------------------------------------------------------------------------


def test_a_publishable_column_costs_no_diagnosis() -> None:
    # None is what tells a caller "the values are fine" and lets every other
    # refusal — a wrong dtype, a wrong shape — keep its own message.
    assert (
        diagnose_continuous_payload(
            np.ones(8, dtype=np.float32), kind="gene", name="Gene_A"
        )
        is None
    )
    assert (
        diagnose_continuous_payload(
            np.array(["a", "b"]), kind="field", name="label"
        )
        is None
    )


def test_every_unpublishable_kind_is_counted_separately() -> None:
    values = np.ones(64, dtype=np.float64)
    values[1] = np.nan
    values[2] = np.inf
    values[3] = -np.inf
    values[4] = 1e300  # finite in float64, infinite in float32
    values[5] = 1e-320  # nonzero in float64, exactly zero in float32

    diagnosis = diagnose_continuous_payload(values, kind="gene", name="Gene_A")
    assert diagnosis is not None
    assert diagnosis.total == 64
    assert diagnosis.nan == 1
    assert diagnosis.positive_infinity == 1
    assert diagnosis.negative_infinity == 1
    assert diagnosis.float32_overflow == 1
    assert diagnosis.float32_underflow == 1
    assert diagnosis.offending == 5
    assert diagnosis.examples == (1, 2, 3, 4, 5)

    message = str(diagnosis)
    assert "of 64 cells" in message
    assert "1 NaN" in message
    assert "1 +Inf" in message
    assert "1 -Inf" in message
    assert "beyond the float32 range" in message
    assert "below the smallest float32" in message


def test_the_examples_are_bounded_and_the_count_is_not() -> None:
    values = np.full(10_000, np.nan)
    diagnosis = diagnose_continuous_payload(values, kind="gene", name="Gene_A")
    assert diagnosis is not None
    assert diagnosis.nan == 10_000
    assert len(diagnosis.examples) == NON_FINITE_EXAMPLE_LIMIT
    assert "10,000 NaN" in str(diagnosis)
    assert ", ..." in str(diagnosis)


def test_a_gene_and_a_field_are_given_different_remedies() -> None:
    values = np.array([1.0, np.inf])
    gene = str(diagnose_continuous_payload(values, kind="gene", name="Gene_A"))
    field = str(diagnose_continuous_payload(values, kind="field", name="score"))
    assert gene.startswith("Gene 'Gene_A' cannot be published")
    assert "adata.X" in gene and "adata.obs" not in gene
    assert field.startswith("Continuous obs field 'score' cannot be published")
    assert "adata.obs['score']" in field


def test_the_diagnosis_refuses_an_unknown_kind_or_empty_name() -> None:
    with pytest.raises(ValueError, match='"gene" or "field"'):
        diagnose_continuous_payload(np.ones(2), kind="other", name="x")
    with pytest.raises(ValueError, match="non-empty string"):
        diagnose_continuous_payload(np.ones(2), kind="gene", name="")


# ---------------------------------------------------------------------------
# The adapter raises it
# ---------------------------------------------------------------------------


def test_the_adapter_refuses_a_non_finite_gene_with_its_counts() -> None:
    adapter = _adapter(
        _adata(gene_values=np.array([1.0, np.nan, 3.0, np.inf]))
    )
    with pytest.raises(NonFinitePayloadError) as refusal:
        adapter.get_gene_expression("Gene_A")
    diagnosis = refusal.value.diagnosis
    assert diagnosis.kind == "gene"
    assert diagnosis.name == "Gene_A"
    assert diagnosis.nan == 1
    assert diagnosis.positive_infinity == 1
    assert diagnosis.examples == (1, 3)


def test_the_adapter_refuses_a_non_finite_obs_column_the_same_way() -> None:
    # Continuous obs columns are classified when the adapter is built, so this
    # one refuses to start the server rather than to answer one request. That is
    # the message with the least room to be vague.
    with pytest.raises(NonFinitePayloadError) as refusal:
        _adapter(_adata(obs_values=np.array([0.25, np.inf, 0.75, 1.0])))
    diagnosis = refusal.value.diagnosis
    assert diagnosis.kind == "field"
    assert diagnosis.name == "score"
    assert diagnosis.positive_infinity == 1
    assert diagnosis.examples == (1,)


def test_a_publishable_object_is_unaffected() -> None:
    adapter = _adapter(_adata())
    assert len(adapter.get_gene_expression("Gene_A")) == 4 * 4
    assert len(adapter.get_obs_continuous_values("score")) == 4 * 4


# ---------------------------------------------------------------------------
# The server answers 422 with a body the browser can present
# ---------------------------------------------------------------------------


def test_a_non_finite_gene_is_422_with_a_json_diagnosis() -> None:
    adata = _adata(gene_values=np.array([1.0, np.nan, 3.0, 4.0]))
    with _running_server(adata) as (host, port):
        status, content_type, body = _get(host, port, "/var/0.values.f32")

    assert status == 422, "a refused payload is not an internal server error"
    assert content_type == "application/json"
    payload = json.loads(body)
    assert payload["error"] == "non_finite_continuous_payload"
    assert payload["kind"] == "gene"
    assert payload["name"] == "Gene_A"
    assert payload["counts"]["nan"] == 1
    assert payload["counts"]["total"] == 4
    assert payload["counts"]["offending"] == 1
    assert payload["examples"] == [1]
    assert "Gene 'Gene_A' cannot be published" in payload["message"]
    # The message must carry the remedy, because the browser shows it verbatim.
    assert "np.nan_to_num" in payload["message"]


def test_a_non_finite_obs_column_refuses_to_start_the_server() -> None:
    # Continuous obs columns are classified while the adapter is built, so an
    # unusable one stops the server rather than one request. The message is the
    # whole diagnosis, because it is the only thing the operator will see.
    adata = _adata(obs_values=np.array([0.25, 0.5, -np.inf, 1.0]))
    with pytest.raises(NonFinitePayloadError) as refusal:
        with _running_server(adata):
            pass
    diagnosis = refusal.value.diagnosis
    assert diagnosis.kind == "field"
    assert diagnosis.name == "score"
    assert diagnosis.negative_infinity == 1
    assert diagnosis.examples == (2,)


def test_a_column_that_turns_unusable_while_serving_is_422() -> None:
    # `serve_anndata` holds the caller's object, and a notebook can assign to it
    # after the server is up. The per-request validation is the reason that is
    # answered rather than crashed on, and it answers with the same diagnosis.
    adata = _adata()
    with _running_server(adata) as (host, port):
        healthy_status, _, _ = _get(host, port, "/obs/0.values.f32")
        assert healthy_status == 200

        adata.obs["score"] = np.array([0.25, np.inf, 0.75, 1.0], dtype=np.float32)
        status, content_type, body = _get(host, port, "/obs/0.values.f32")

    assert status == 422
    assert content_type == "application/json"
    payload = json.loads(body)
    assert payload["kind"] == "field"
    assert payload["name"] == "score"
    assert payload["counts"]["positive_infinity"] == 1
    assert payload["examples"] == [1]


def test_a_publishable_object_still_serves_binary() -> None:
    with _running_server(_adata()) as (host, port):
        gene_status, gene_type, gene_body = _get(host, port, "/var/0.values.f32")
        obs_status, _, obs_body = _get(host, port, "/obs/0.values.f32")

    assert gene_status == 200
    assert gene_type == "application/octet-stream"
    assert np.frombuffer(gene_body, dtype="<f4").tolist() == [1.0, 2.0, 3.0, 4.0]
    assert obs_status == 200
    assert np.frombuffer(obs_body, dtype="<f4").tolist() == [0.25, 0.5, 0.75, 1.0]


def test_an_unknown_gene_is_still_a_404() -> None:
    # 422 means "this exists and cannot be drawn". It must not swallow "this
    # does not exist", which is a different thing for a reader to do about.
    with _running_server(_adata()) as (host, port):
        status, _, _ = _get(host, port, "/var/99.values.f32")
    assert status == 404


# ---------------------------------------------------------------------------
# The cache that holds served columns is bounded by bytes
# ---------------------------------------------------------------------------


def test_the_gene_cache_is_bounded_by_bytes_not_by_entries() -> None:
    # Counting entries is what made this dangerous: one hundred gene columns is
    # 29 MB on a 72,000-cell object and 7.3 GB on an 18,142,044-cell one, and the
    # second is exactly the size of object the cache exists for.
    assert GENE_EXPRESSION_CACHE_BYTES == 256 * 1024 * 1024

    cache = LRUCache(max_bytes=1000)
    column = np.zeros(100, dtype=np.float32)  # 400 bytes each
    cache.put("a", column)
    cache.put("b", column.copy())
    assert len(cache) == 2
    assert cache.resident_bytes == 800

    cache.put("c", column.copy())
    assert cache.resident_bytes <= 1000
    assert cache.get("a") is None, "the least recently used entry is evicted"
    assert cache.get("c") is not None

    # A single value larger than the whole budget is not cached, and is not an
    # error: a cache must not decide what the server may publish.
    cache.put("huge", np.zeros(1000, dtype=np.float32))
    assert cache.get("huge") is None
    assert cache.resident_bytes <= 1000

    cache.clear()
    assert len(cache) == 0
    assert cache.resident_bytes == 0

    # One 18.1M-cell gene column is 72.6 MB and still fits.
    assert 18_142_044 * 4 <= GENE_EXPRESSION_CACHE_BYTES


def test_a_served_gene_round_trips_through_the_cache() -> None:
    adapter = _adapter(_adata())
    first = adapter.get_gene_expression("Gene_A")
    second = adapter.get_gene_expression("Gene_A")
    assert first == second
    assert adapter._gene_expression_cache.resident_bytes == 4 * 4


def test_a_sparse_matrix_is_diagnosed_like_a_dense_one() -> None:
    dense = np.array([[1.0], [np.inf], [3.0], [4.0]], dtype=np.float32)
    adata = _adata()
    adata.X = sparse.csr_matrix(dense)
    adapter = _adapter(adata)
    with pytest.raises(NonFinitePayloadError) as refusal:
        adapter.get_gene_expression("Gene_A")
    assert refusal.value.diagnosis.positive_infinity == 1
    assert refusal.value.diagnosis.examples == (1,)
