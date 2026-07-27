from __future__ import annotations

import hashlib
import json
import urllib.error
import urllib.request
from pathlib import Path

import pytest

from cellucid.server import CellucidServer

STATE_MANIFEST_NAME = "state-snapshots.json"
STATE_BUNDLE_NAME = "default.cellucid-session"
STATE_MANIFEST_BYTES = b'{ "states": ["default.cellucid-session"] }\n'


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


def _publish_state(dataset: Path, bundle: bytes = b"exact published state") -> bytes:
    assert len(STATE_MANIFEST_BYTES) == 43
    (dataset / STATE_MANIFEST_NAME).write_bytes(STATE_MANIFEST_BYTES)
    (dataset / STATE_BUNDLE_NAME).write_bytes(bundle)
    return bundle


def _server(data_dir: Path) -> CellucidServer:
    return CellucidServer(
        data_dir,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
    )


def _request(
    url: str,
    *,
    method: str = "GET",
    headers: dict[str, str] | None = None,
) -> tuple[int, bytes, urllib.response.addinfourl]:
    request = urllib.request.Request(
        url,
        method=method,
        headers={} if headers is None else headers,
    )
    try:
        response = urllib.request.urlopen(request, timeout=5)
    except urllib.error.HTTPError as error:
        return error.code, error.read(), error
    with response:
        return response.status, response.read(), response


def _exact_size_json(value: object) -> bytes:
    compact = json.dumps(value, separators=(",", ":")).encode("utf-8")
    assert len(compact) <= 42
    return compact + (b" " * (42 - len(compact))) + b"\n"


def test_absent_published_state_preserves_catalog_and_artifact_inventory(
    tmp_path: Path,
) -> None:
    dataset = _prepared_dataset(tmp_path / "dataset", dataset_id="absent-state")
    (dataset / "unadvertised.cellucid-session").write_bytes(b"manual user state")
    server = _server(dataset)

    assert server._datasets == [
        {
            "id": "absent-state",
            "path": "/",
            "name": "Prepared absent-state",
        }
    ]
    assert STATE_MANIFEST_NAME not in server._artifact_inventory
    assert STATE_BUNDLE_NAME not in server._artifact_inventory
    assert "unadvertised.cellucid-session" not in server._artifact_inventory


@pytest.mark.parametrize("mode", ["direct", "multi"])
def test_exact_published_state_is_advertised_and_served_in_every_prepared_mode(
    tmp_path: Path,
    mode: str,
) -> None:
    if mode == "direct":
        data_dir = _prepared_dataset(tmp_path / "dataset", dataset_id="state-direct")
        dataset = data_dir
        request_prefix = ""
        public_path = "/"
    else:
        data_dir = tmp_path / "exports"
        data_dir.mkdir()
        dataset = _prepared_dataset(data_dir / "alpha", dataset_id="state-multi")
        request_prefix = "alpha/"
        public_path = "/alpha/"

    bundle = _publish_state(dataset)
    digest = hashlib.sha256(bundle).hexdigest()
    server = _server(data_dir)

    assert server._datasets == [
        {
            "id": f"state-{mode}",
            "path": public_path,
            "name": f"Prepared state-{mode}",
            "state_manifest": STATE_MANIFEST_NAME,
            "state_sha256": digest,
        }
    ]
    assert f"{request_prefix}{STATE_MANIFEST_NAME}" in server._artifact_inventory
    assert f"{request_prefix}{STATE_BUNDLE_NAME}" in server._artifact_inventory
    assert digest == digest.lower()
    assert len(digest) == 64

    server.start_background()
    try:
        status, catalog_body, _catalog_headers = _request(
            f"{server.url}/_cellucid/datasets"
        )
        assert status == 200
        assert json.loads(catalog_body) == {"datasets": server._datasets}

        status, manifest_body, manifest_headers = _request(
            f"{server.url}/{request_prefix}{STATE_MANIFEST_NAME}",
            headers={"Origin": "https://www.cellucid.com"},
        )
        assert status == 200
        assert manifest_body == STATE_MANIFEST_BYTES
        assert manifest_headers.headers["Content-Type"] == "application/json"
        assert (
            manifest_headers.headers["Access-Control-Allow-Origin"]
            == "https://www.cellucid.com"
        )

        status, bundle_body, bundle_headers = _request(
            f"{server.url}/{request_prefix}{STATE_BUNDLE_NAME}",
            headers={"Range": "bytes=2-7"},
        )
        assert status == 206
        assert bundle_body == bundle[2:8]
        assert bundle_headers.headers["Content-Type"] == "application/octet-stream"
        assert bundle_headers.headers["Content-Range"] == f"bytes 2-7/{len(bundle)}"
        assert bundle_headers.headers["Accept-Ranges"] == "bytes"

        status, head_body, head_headers = _request(
            f"{server.url}/{request_prefix}{STATE_BUNDLE_NAME}",
            method="HEAD",
        )
        assert status == 200
        assert head_body == b""
        assert int(head_headers.headers["Content-Length"]) == len(bundle)
    finally:
        server.stop()


def test_multi_dataset_catalog_augments_only_the_exact_opted_in_generation(
    tmp_path: Path,
) -> None:
    exports = tmp_path / "exports"
    exports.mkdir()
    published = _prepared_dataset(exports / "alpha", dataset_id="published")
    _prepared_dataset(exports / "beta", dataset_id="plain")
    bundle = _publish_state(published)

    server = _server(exports)

    assert server._datasets == [
        {
            "id": "published",
            "path": "/alpha/",
            "name": "Prepared published",
            "state_manifest": STATE_MANIFEST_NAME,
            "state_sha256": hashlib.sha256(bundle).hexdigest(),
        },
        {
            "id": "plain",
            "path": "/beta/",
            "name": "Prepared plain",
        },
    ]
    assert f"alpha/{STATE_MANIFEST_NAME}" in server._artifact_inventory
    assert f"alpha/{STATE_BUNDLE_NAME}" in server._artifact_inventory
    assert f"beta/{STATE_MANIFEST_NAME}" not in server._artifact_inventory
    assert f"beta/{STATE_BUNDLE_NAME}" not in server._artifact_inventory


@pytest.mark.parametrize(
    ("manifest_bytes", "message"),
    [
        (b'{"states":["default.cellucid-session"]}', "43 bytes"),
        (_exact_size_json([]), "JSON object"),
        (_exact_size_json({"states": [], "extra": 1}), "exact keys"),
        (_exact_size_json({"states": []}), "exactly"),
        (_exact_size_json({"states": ["x", "y"]}), "exactly"),
        (_exact_size_json({"states": ["other.cellucid-session"]}), "exactly"),
        (_exact_size_json({"states": [1]}), "exactly"),
        (b"{" + (b" " * 42), "UTF-8 JSON"),
        ((b"\xff" * 43), "UTF-8 JSON"),
    ],
)
def test_published_state_manifest_is_one_exact_closed_43_byte_object(
    tmp_path: Path,
    manifest_bytes: bytes,
    message: str,
) -> None:
    dataset = _prepared_dataset(tmp_path / "dataset", dataset_id="invalid-manifest")
    _publish_state(dataset)
    (dataset / STATE_MANIFEST_NAME).write_bytes(manifest_bytes)

    with pytest.raises((TypeError, ValueError), match=message):
        _server(dataset)


@pytest.mark.parametrize(
    ("present_name", "missing_name"),
    [
        (STATE_MANIFEST_NAME, STATE_BUNDLE_NAME),
        (STATE_BUNDLE_NAME, STATE_MANIFEST_NAME),
    ],
)
def test_published_state_rejects_a_partial_sidecar_pair(
    tmp_path: Path,
    present_name: str,
    missing_name: str,
) -> None:
    dataset = _prepared_dataset(tmp_path / "dataset", dataset_id="partial-state")
    payload = STATE_MANIFEST_BYTES if present_name == STATE_MANIFEST_NAME else b"bundle"
    (dataset / present_name).write_bytes(payload)

    with pytest.raises(ValueError, match=rf"{present_name}.*{missing_name}|exact pair"):
        _server(dataset)


@pytest.mark.parametrize(
    ("bundle", "message"),
    [
        (b"", "1 through 32768"),
        (b"x" * 32769, "1 through 32768"),
    ],
    ids=("empty", "one-byte-over-limit"),
)
def test_published_state_bundle_has_one_exact_bounded_transport(
    tmp_path: Path,
    bundle: bytes,
    message: str,
) -> None:
    dataset = _prepared_dataset(tmp_path / "dataset", dataset_id="bounded-state")
    _publish_state(dataset, bundle)

    with pytest.raises(ValueError, match=message):
        _server(dataset)


@pytest.mark.parametrize("size", [1, 32768])
def test_published_state_bundle_accepts_both_exact_size_boundaries(
    tmp_path: Path,
    size: int,
) -> None:
    dataset = _prepared_dataset(tmp_path / "dataset", dataset_id=f"boundary-{size}")
    bundle = _publish_state(dataset, b"x" * size)

    server = _server(dataset)

    assert server._datasets[0]["state_sha256"] == hashlib.sha256(bundle).hexdigest()
    assert server._artifact_inventory[STATE_BUNDLE_NAME].size == size


@pytest.mark.parametrize("sidecar_name", [STATE_MANIFEST_NAME, STATE_BUNDLE_NAME])
def test_published_state_sidecars_must_be_regular_non_symlink_in_root_files(
    tmp_path: Path,
    sidecar_name: str,
) -> None:
    dataset = _prepared_dataset(tmp_path / "dataset", dataset_id="symlink-state")
    _publish_state(dataset)
    sidecar = dataset / sidecar_name
    outside = tmp_path / f"outside-{sidecar_name}"
    outside.write_bytes(sidecar.read_bytes())
    sidecar.unlink()
    try:
        sidecar.symlink_to(outside)
    except OSError as error:
        pytest.skip(f"File symlinks are unavailable on this platform: {error}")

    with pytest.raises(ValueError, match="symbolic link"):
        _server(dataset)


def test_published_state_rejects_extra_session_files(tmp_path: Path) -> None:
    dataset = _prepared_dataset(tmp_path / "dataset", dataset_id="extra-state")
    _publish_state(dataset)
    (dataset / "other.cellucid-session").write_bytes(b"undeclared state")

    with pytest.raises(ValueError, match="extra.*session|exactly one"):
        _server(dataset)


def test_multi_dataset_mode_rejects_one_invalid_publication_atomically(
    tmp_path: Path,
) -> None:
    exports = tmp_path / "exports"
    exports.mkdir()
    first = _prepared_dataset(exports / "first", dataset_id="first")
    second = _prepared_dataset(exports / "second", dataset_id="second")
    _publish_state(first)
    (second / STATE_MANIFEST_NAME).write_bytes(STATE_MANIFEST_BYTES)

    with pytest.raises(ValueError, match="exact pair"):
        _server(exports)
