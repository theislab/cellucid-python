"""The diagnostic dataset list must accept exactly what the server serves.

``debug_connection()`` fetches ``/_cellucid/datasets`` and validates it with
:func:`cellucid.jupyter._wire._require_debug_dataset_list` before probing every
declared identity. The validator and the server hold two separate declarations
of the same entry shape, and nothing forced them to agree: a prepared export
that ships a published default session serves a five-key entry, the validator
accepted only three, and the report silently lost ``dataset_identity_probes``.

These tests tie the two declarations together in both directions -- what the
server emits is accepted, and the exactness that makes the validator worth
having is preserved for everything else.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

from cellucid import jupyter
from cellucid.jupyter._wire import (
    _PUBLISHED_STATE_FIELDS,
    _PUBLISHED_STATE_MANIFEST_NAME,
    _require_debug_dataset_list,
)
from cellucid.server._datasets import _list_exported_datasets
from cellucid.server._state import (
    _PUBLISHED_STATE_BUNDLE,
    _PUBLISHED_STATE_MANIFEST,
)

STATE_MANIFEST_BYTES = b'{ "states": ["default.cellucid-session"] }\n'
STATE_BUNDLE_BYTES = b"exact published state"


def _prepared_dataset(root: Path, *, dataset_id: str) -> Path:
    root.mkdir(parents=True)
    (root / "dataset_identity.json").write_text(
        json.dumps(
            {
                "version": 2,
                "id": dataset_id,
                "name": f"Prepared {dataset_id}",
                "description": "",
                "stats": {
                    "n_cells": 1,
                    "n_genes": 0,
                    "n_obs_fields": 0,
                    "n_categorical_fields": 0,
                    "n_continuous_fields": 0,
                    "has_connectivity": False,
                    "n_edges": None,
                },
                "embeddings": {
                    "available_dimensions": [2],
                    "default_dimension": 2,
                    "files": {"2d": "points_2d.bin"},
                },
                "obs_fields": [],
                "export_settings": {
                    "compression": None,
                    "var_quantization": None,
                    "obs_continuous_quantization": None,
                    "obs_categorical_dtype": "uint16",
                },
            }
        ),
        encoding="utf-8",
    )
    (root / "obs_manifest.json").write_text(
        json.dumps(
            {
                "_format": "compact_v1",
                "n_points": 1,
                "centroid_outlier_quantile": None,
                "latent_key": "latent_space",
                "compression": None,
                "_obsSchemas": {},
                "_continuousFields": [],
                "_categoricalFields": [],
            }
        ),
        encoding="utf-8",
    )
    (root / "points_2d.bin").write_bytes(b"\x00" * 8)
    return root


def _publish_state(dataset: Path) -> str:
    assert len(STATE_MANIFEST_BYTES) == 43
    (dataset / _PUBLISHED_STATE_MANIFEST).write_bytes(STATE_MANIFEST_BYTES)
    (dataset / _PUBLISHED_STATE_BUNDLE).write_bytes(STATE_BUNDLE_BYTES)
    return hashlib.sha256(STATE_BUNDLE_BYTES).hexdigest()


def _entry(**overrides: object) -> dict[str, object]:
    entry: dict[str, object] = {
        "id": "alpha",
        "path": "/alpha/",
        "name": "Alpha",
        "state_manifest": _PUBLISHED_STATE_MANIFEST_NAME,
        "state_sha256": "0" * 64,
    }
    entry.update(overrides)
    return entry


def test_wire_and_server_declare_the_same_published_state_pair() -> None:
    """One published-state shape, two modules, and now one check that ties them."""
    assert _PUBLISHED_STATE_MANIFEST_NAME == _PUBLISHED_STATE_MANIFEST
    assert sorted(_PUBLISHED_STATE_FIELDS) == ["state_manifest", "state_sha256"]


def test_diagnostics_accept_the_exact_catalog_the_prepared_server_serves(
    tmp_path: Path,
) -> None:
    """The served entry is the input; no hand-written stand-in can drift from it."""
    exports = tmp_path / "exports"
    exports.mkdir()
    _prepared_dataset(exports / "alpha", dataset_id="alpha")
    beta = _prepared_dataset(exports / "beta", dataset_id="beta")
    digest = _publish_state(beta)

    served = _list_exported_datasets(exports)
    assert served == [
        {"id": "alpha", "path": "/alpha/", "name": "Prepared alpha"},
        {
            "id": "beta",
            "path": "/beta/",
            "name": "Prepared beta",
            "state_manifest": _PUBLISHED_STATE_MANIFEST,
            "state_sha256": digest,
        },
    ]

    assert _require_debug_dataset_list({"datasets": served}) == served


def test_debug_connection_probes_identity_for_a_published_state_catalog(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The diagnostic's whole point is the identity probe; it must still run."""
    viewer = jupyter.BaseViewer(
        port=8765,
        client_server_url="http://127.0.0.1:8765",
        web_cache_dir=tmp_path / "cache",
    )
    responses: dict[str, object] = {
        "http://127.0.0.1:8765/_cellucid/health": {"status": "ok"},
        "http://127.0.0.1:8765/_cellucid/info": {"version": "0.9.1"},
        "http://127.0.0.1:8765/_cellucid/datasets": {
            "datasets": [
                {
                    "id": "beta",
                    "path": "/beta/",
                    "name": "Prepared beta",
                    "state_manifest": _PUBLISHED_STATE_MANIFEST_NAME,
                    "state_sha256": "b" * 64,
                }
            ]
        },
        "http://127.0.0.1:8765/beta/dataset_identity.json": {
            "id": "beta",
            "name": "Prepared beta",
        },
    }

    class Response:
        def __init__(self, payload: object) -> None:
            self._payload = json.dumps(payload).encode("utf-8")
            self.headers: dict[str, str] = {}

        def __enter__(self) -> Response:
            return self

        def __exit__(self, *_: object) -> None:
            return None

        def read(self) -> bytes:
            return self._payload

    def urlopen(request: Any, *, timeout: float) -> Response:
        url = request.full_url if hasattr(request, "full_url") else request
        if url not in responses:
            raise OSError(f"unmapped diagnostic URL: {url}")
        return Response(responses[url])

    monkeypatch.setattr("urllib.request.urlopen", urlopen)
    report = viewer.debug_connection()

    assert "server_datasets_error" not in report
    assert report["server_datasets"] == [
        {
            "id": "beta",
            "path": "/beta/",
            "name": "Prepared beta",
            "state_manifest": _PUBLISHED_STATE_MANIFEST_NAME,
            "state_sha256": "b" * 64,
        }
    ]
    assert report["dataset_identity_probes"]["beta"] == {
        "path": "/beta/",
        "url": "http://127.0.0.1:8765/beta/dataset_identity.json",
        "identity": {"id": "beta", "name": "Prepared beta"},
    }


def test_diagnostics_still_accept_a_catalog_without_a_published_state() -> None:
    payload = {"datasets": [{"id": "alpha", "path": "/alpha/", "name": "Alpha"}]}
    assert _require_debug_dataset_list(payload) == [
        {"id": "alpha", "path": "/alpha/", "name": "Alpha"}
    ]


@pytest.mark.parametrize(
    ("entry", "match"),
    [
        pytest.param(
            {"id": "alpha", "path": "/alpha/", "name": "Alpha", "state_manifest": "x.json"},
            "exactly id, name, and path",
            id="manifest-without-digest",
        ),
        pytest.param(
            {"id": "alpha", "path": "/alpha/", "name": "Alpha", "state_sha256": "0" * 64},
            "exactly id, name, and path",
            id="digest-without-manifest",
        ),
        pytest.param(
            _entry(unexpected="value"),
            "exactly id, name, and path",
            id="unknown-extra-key",
        ),
        pytest.param(
            _entry(state_manifest="states.json"),
            "state_manifest must be exactly",
            id="foreign-manifest-name",
        ),
        pytest.param(
            _entry(state_manifest=None),
            "state_manifest must be exactly",
            id="manifest-not-text",
        ),
        pytest.param(
            _entry(state_sha256="0" * 63),
            "state_sha256 must be exactly",
            id="digest-too-short",
        ),
        pytest.param(
            _entry(state_sha256="A" * 64),
            "state_sha256 must be exactly",
            id="digest-uppercase",
        ),
        pytest.param(
            _entry(state_sha256="../etc/passwd".ljust(64, "0")),
            "state_sha256 must be exactly",
            id="digest-not-hexadecimal",
        ),
        pytest.param(
            _entry(state_sha256=0),
            "state_sha256 must be exactly",
            id="digest-not-text",
        ),
    ],
)
def test_diagnostics_reject_every_malformed_published_state_declaration(
    entry: dict[str, object],
    match: str,
) -> None:
    with pytest.raises((TypeError, ValueError), match=match):
        _require_debug_dataset_list({"datasets": [entry]})
