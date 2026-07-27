from __future__ import annotations

import gzip
import hashlib
import inspect
import json
import os
import stat
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid._compression import deterministic_gzip_compress
from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.prepare_data import prepare

CANONICAL_GZIP_HEADER = bytes.fromhex("1f8b08000000000000ff")
GOLDEN_PAYLOAD_SHA256 = "24d39f0523c693080fc935a0de51d88def96ff9ab5492a07a629f35203a79c09"
FIXED_CREATED_AT = "2026-07-26T12:34:56Z"


def _assert_canonical_gzip(payload: bytes) -> None:
    assert payload[:10] == CANONICAL_GZIP_HEADER
    assert payload[3] & 0x08 == 0
    assert int.from_bytes(payload[4:8], byteorder="little") == 0
    assert payload[9] == 255
    assert gzip.decompress(payload)


def _prepared_fixture(output: Path, *, created_at: str | None = FIXED_CREATED_AT) -> None:
    embedding = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ],
        dtype=np.float32,
    )
    vectors = np.array(
        [
            [0.1, 0.2],
            [0.2, 0.1],
            [-0.1, 0.2],
            [0.2, -0.1],
        ],
        dtype=np.float32,
    )
    prepare(
        latent_space=embedding.copy(),
        obs=pd.DataFrame(
            {
                "score": np.array([0.25, 0.5, 0.75, 1.0], dtype=np.float32),
                "cluster": pd.Categorical(["a", "a", "b", "b"]),
            }
        ),
        var=pd.DataFrame(index=["GeneA", "GeneB"]),
        gene_expression=np.array(
            [
                [1.0, 2.0],
                [2.0, 3.0],
                [3.0, 4.0],
                [4.0, 5.0],
            ],
            dtype=np.float32,
        ),
        connectivities=sparse.csr_matrix(
            np.array(
                [
                    [0.0, 1.0, 0.0, 0.0],
                    [1.0, 0.0, 0.5, 0.0],
                    [0.0, 0.5, 0.0, 1.0],
                    [0.0, 0.0, 1.0, 0.0],
                ],
                dtype=np.float64,
            )
        ),
        X_umap_2d=embedding,
        vector_fields={"velocity_umap_2d": vectors},
        out_dir=output,
        dataset_name="Deterministic gzip",
        dataset_id="deterministic-gzip",
        created_at=created_at,
        obs_categorical_dtype="uint16",
        centroid_min_points=1,
        compression=6,
    )


def _snapshot(root: Path) -> dict[str, bytes]:
    return {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _direct_adapter() -> AnnDataAdapter:
    embedding = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ],
        dtype=np.float32,
    )
    adata = ad.AnnData(
        X=np.array(
            [
                [1.0, 2.0],
                [2.0, 3.0],
                [3.0, 4.0],
                [4.0, 5.0],
            ],
            dtype=np.float32,
        ),
        obs=pd.DataFrame(
            {
                "score": np.array([0.25, 0.5, 0.75, 1.0], dtype=np.float32),
                "cluster": pd.Categorical(["a", "a", "b", "b"]),
            },
            index=["cell-1", "cell-2", "cell-3", "cell-4"],
        ),
        var=pd.DataFrame(index=["GeneA", "GeneB"]),
    )
    adata.obsm["X_umap_2d"] = embedding
    adata.obsm["X_pca"] = embedding.copy()
    adata.obsm["velocity_umap_2d"] = np.array(
        [
            [0.1, 0.2],
            [0.2, 0.1],
            [-0.1, 0.2],
            [0.2, -0.1],
        ],
        dtype=np.float32,
    )
    adata.obsp["connectivities"] = sparse.csr_matrix(
        np.array(
            [
                [0.0, 1.0, 0.0, 0.0],
                [1.0, 0.0, 0.5, 0.0],
                [0.0, 0.5, 0.0, 1.0],
                [0.0, 0.0, 1.0, 0.0],
            ],
            dtype=np.float64,
        )
    )
    return AnnDataAdapter(
        adata,
        latent_key="X_pca",
        centroid_min_points=1,
        dataset_name="Deterministic direct AnnData",
        dataset_id="deterministic-direct-anndata",
    )


def test_canonical_gzip_bytes_are_time_and_platform_header_independent() -> None:
    first = deterministic_gzip_compress(b"cellucid-gzip-contract", compresslevel=6)
    second = deterministic_gzip_compress(b"cellucid-gzip-contract", compresslevel=6)

    assert first == second
    assert hashlib.sha256(first).hexdigest() == GOLDEN_PAYLOAD_SHA256
    _assert_canonical_gzip(first)
    assert gzip.decompress(first) == b"cellucid-gzip-contract"


def test_created_at_is_one_optional_keyword_only_public_parameter() -> None:
    parameter = inspect.signature(prepare).parameters["created_at"]

    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is None


def test_repeated_prepared_exports_have_identical_complete_trees(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    clock = {"value": 1_700_000_000.0}
    monkeypatch.setattr(gzip.time, "time", lambda: clock["value"])
    first = tmp_path / "first" / "export"
    second = tmp_path / "second" / "export"

    _prepared_fixture(first)
    clock["value"] = 1_800_000_000.0
    _prepared_fixture(second)

    first_tree = _snapshot(first)
    second_tree = _snapshot(second)
    assert first_tree == second_tree
    identity = json.loads(first_tree["dataset_identity.json"])
    assert identity["created_at"] == FIXED_CREATED_AT

    first_payloads = {
        relative_path: payload
        for relative_path, payload in first_tree.items()
        if relative_path.endswith(".gz")
    }
    assert len(first_payloads) >= 10
    for payload in first_payloads.values():
        _assert_canonical_gzip(payload)
    if os.name == "posix":
        for path in first.rglob("*.gz"):
            assert stat.S_IMODE(path.stat().st_mode) == 0o644


@pytest.mark.parametrize(
    ("created_at", "expected_error"),
    [
        (123, TypeError),
        ("", ValueError),
        ("2026-7-26T12:34:56Z", ValueError),
        ("2026-07-26 12:34:56Z", ValueError),
        ("2026-07-26T12:34:56+00:00", ValueError),
        ("2026-07-26T12:34:56.000Z", ValueError),
        ("2026-02-29T12:34:56Z", ValueError),
        ("2026-07-26T24:00:00Z", ValueError),
    ],
)
def test_invalid_created_at_rejects_before_filesystem_mutation(
    tmp_path: Path,
    created_at: object,
    expected_error: type[Exception],
) -> None:
    parent = tmp_path / "must-not-exist"
    target = parent / "export"

    with pytest.raises(expected_error):
        prepare(
            out_dir=target,
            dataset_name="Rejected timestamp",
            dataset_id="rejected-timestamp",
            created_at=created_at,  # type: ignore[arg-type]
            obs_categorical_dtype="uint16",
        )

    assert not parent.exists()


def test_every_direct_anndata_gzip_payload_is_byte_stable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    adapter = _direct_adapter()
    clock = {"value": 1_700_000_000.0}
    monkeypatch.setattr(gzip.time, "time", lambda: clock["value"])

    def compressed_payloads() -> list[bytes]:
        categorical, _categories, _missing = adapter.get_obs_categorical_codes(
            "cluster",
            compress=True,
        )
        sources, destinations, weights, _n_edges, _max_neighbors = adapter.get_connectivity_edges(
            compress=True
        )
        return [
            adapter.get_points_binary(2, compress=True),
            adapter.get_vector_field_binary(
                "velocity_umap",
                2,
                compress=True,
            ),
            adapter.get_obs_continuous_values("score", compress=True),
            categorical,
            adapter.get_obs_outlier_quantiles("cluster", compress=True),
            adapter.get_gene_expression("GeneA", compress=True),
            sources,
            destinations,
            weights,
        ]

    first = compressed_payloads()
    clock["value"] = 1_800_000_000.0
    second = compressed_payloads()

    assert first == second
    assert len(first) == 9
    for payload in first:
        _assert_canonical_gzip(payload)
