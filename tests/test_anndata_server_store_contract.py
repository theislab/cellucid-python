"""What an AnnData server publishes: the store's bytes, and how to read them.

The viewer reads the store itself (``cellucid/assets/js/data/served-anndata.js``
over ``formats/io/byte-source.js`` and ``formats/zarr/stores/http.js``), so the
server's whole data contract is:

- ``/_cellucid/datasets``: the dataset's id and name, its store's route, kind,
  layout and size, and the launch defaults;
- the store's bytes: one exact ``Range`` per request answered ``206`` with
  ``Content-Range``, ``HEAD`` with the length, an ``ETag``, no transfer
  encoding; a Zarr folder key by key, with 404 for an absent key and
  consolidated metadata built when the store has none;
- HTTP/1.1 keep-alive, answered by a thread per connection.

A file in priority order (``prepare()``'s layout) is served as it is, so the
files here are put in priority order first (:func:`_in_priority_order`); an
object, and a file in any other order, are served as a store in priority
order (``test_served_priority_store.py``).
"""

from __future__ import annotations

import http.client
import json
import os
import shutil
import threading
import zipfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid._anndata_files import (
    browser_undecodable,
    classify_anndata_path,
    read_anndata_contents,
)
from cellucid._file_versions import path_version
from cellucid.anndata_server import AnnDataServer
from cellucid.prepare_data._priority import ORIGINAL_POSITION_COLUMN, CellPriority, prefix_levels

FIXTURES = Path(__file__).parent / "fixtures" / "older_anndata"


def _adata() -> ad.AnnData:
    rng = np.random.default_rng(7)
    adata = ad.AnnData(
        X=sparse.random(12, 5, density=0.4, format="csr", dtype=np.float32, random_state=7),
        obs=pd.DataFrame(
            {
                "cell_type": pd.Categorical(["a", "b", "c"] * 4),
                "score": rng.random(12).astype(np.float32),
            },
            index=pd.Index([f"cell-{i}" for i in range(12)], dtype=object),
        ),
        var=pd.DataFrame(
            {"symbol": [f"S{i}" for i in range(5)]},
            index=pd.Index([f"g{i}" for i in range(5)], dtype=object),
        ),
    )
    adata.obsm["X_umap"] = rng.random((12, 2)).astype(np.float32)
    adata.obsm["X_pca"] = rng.random((12, 6)).astype(np.float32)
    adata.obsm["velocity_umap"] = rng.random((12, 2)).astype(np.float32)
    # Named for X_pca, but narrower than it: no field.
    adata.obsm["drift_pca"] = rng.random((12, 2)).astype(np.float32)
    adata.obsp["connectivities"] = sparse.csr_matrix(np.eye(12, dtype=np.float32))
    adata.layers["counts"] = adata.X.copy()
    adata.raw = adata.copy()
    return adata


def _in_priority_order(adata: ad.AnnData) -> ad.AnnData:
    """``adata`` with its cells in priority order and the records saying so (served as it is)."""
    priority = CellPriority.of(adata.n_obs, 0)
    assert priority is not None
    placed = adata[priority.order].copy()
    placed.obs[ORIGINAL_POSITION_COLUMN] = priority.original.astype(np.uint32)
    placed.uns["cellucid"] = {"priority": priority.record(prefix_levels(adata.n_obs))}
    return placed


@contextmanager
def _serving(data: object, **choices: object) -> Iterator[tuple[AnnDataServer, str, int]]:
    server = AnnDataServer(
        data,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Store contract",
        dataset_id="store-contract",
        **choices,
    )
    server.start_background()
    assert server._server is not None
    host, port = server._server.server_address[:2]
    try:
        yield server, host, port
    finally:
        server.stop()


def _request(
    host: str,
    port: int,
    method: str,
    path: str,
    headers: dict[str, str] | None = None,
) -> tuple[int, dict[str, str], bytes]:
    connection = http.client.HTTPConnection(host, port, timeout=10)
    try:
        connection.request(method, path, headers=headers or {})
        response = connection.getresponse()
        return response.status, dict(response.getheaders()), response.read()
    finally:
        connection.close()


def _json(host: str, port: int, path: str) -> dict:
    status, headers, body = _request(host, port, "GET", path)
    assert status == 200
    assert headers["Content-Type"] == "application/json"
    return json.loads(body)


# =============================================================================
# The datasets document
# =============================================================================


def test_the_catalog_names_the_dataset_its_store_and_nothing_private(tmp_path: Path) -> None:
    private = tmp_path / "private-name-of-a-local-file.h5ad"
    _adata().write_h5ad(private)
    with _serving(private) as (_server, host, port):
        catalog = _json(host, port, "/_cellucid/datasets")
        health = _json(host, port, "/_cellucid/health")
        info = _json(host, port, "/_cellucid/info")
    assert catalog == {
        "datasets": [
            {
                "id": "store-contract",
                "name": "Store contract",
                "anndata": {
                    # Not in priority order: served as a store in priority order.
                    "file": {
                        "url": "_cellucid/anndata/store-contract.zarr/",
                        "kind": "zarr",
                        "layout": "folder",
                    },
                    # Nothing was chosen at launch: the viewer takes the
                    # store's own uns['cellucid'] defaults, else its rules.
                    "defaults": {},
                },
            }
        ]
    }
    assert "private-name" not in json.dumps([catalog, health, info])
    assert str(tmp_path) not in json.dumps([catalog, health, info])
    assert set(health) == {"status", "type", "version", "format", "n_cells", "n_genes"}
    assert (health["type"], health["format"], health["n_cells"], health["n_genes"]) == (
        "anndata",
        "zarr",
        12,
        5,
    )
    assert set(info) == {"version", "type", "format", "host", "port", "n_cells", "n_genes", "runtime"}
    # The program serving the viewer, by a name the shared protocol declares.
    protocol = json.loads(
        (Path(__file__).parent / "fixtures" / "viewer_protocol.json").read_text(encoding="utf-8")
    )
    assert info["runtime"] == "python"
    assert info["runtime"] in protocol["runtimes"]


def test_every_launch_choice_is_checked_and_published_in_the_viewer_vocabulary() -> None:
    with _serving(
        _adata(),
        embedding_keys={"2d": ("X_pca", [3, 1]), "3d": "X_pca"},
        knn_key="connectivities",
        expression_source="layer:counts",
        gene_id_column="symbol",
        latent_key="X_pca",
        vector_field_default="velocity",
        obs_keys=["score", "cell_type"],
        normalize_embeddings=False,
        centroid_min_points=3,
    ) as (_server, host, port):
        defaults = _json(host, port, "/_cellucid/datasets")["datasets"][0]["anndata"]["defaults"]
    assert defaults == {
        "embedding2d": {"key": "X_pca", "columns": [3, 1]},
        "embedding3d": {"key": "X_pca", "columns": [0, 1, 2]},
        "knn": "connectivities",
        "expression": "layers/counts",
        "geneIdColumn": "symbol",
        "latentKey": "X_pca",
        "vectorFieldDefault": "velocity",
        "obsKeys": ["score", "cell_type"],
        "normalizeEmbeddings": False,
        "centroidMinPoints": 3,
    }


@pytest.mark.parametrize(
    ("expression", "published"),
    [("X", "X"), ("raw", "raw/X"), ("layer:counts", "layers/counts")],
)
def test_each_expression_source_has_one_viewer_spelling(expression: str, published: str) -> None:
    with _serving(_adata(), expression_source=expression) as (_server, host, port):
        defaults = _json(host, port, "/_cellucid/datasets")["datasets"][0]["anndata"]["defaults"]
    assert defaults["expression"] == published


@pytest.mark.parametrize(
    ("choices", "message"),
    [
        ({"knn_key": "distances"}, r"knn_key names obsp\['distances'\]"),
        ({"embedding_keys": {"2d": "X_tsne"}}, r"obsm\['X_tsne'\], which this store does not have"),
        ({"embedding_keys": {"3d": "X_umap"}}, "needs at least 3 columns"),
        ({"embedding_keys": {"2d": ("X_pca", [0, 0])}}, "different column"),
        ({"embedding_keys": {"2d": ("X_pca", [0, 9])}}, "does not exist"),
        ({"embedding_keys": {"4d": "X_pca"}}, "'1d', '2d', or '3d'"),
        ({"expression_source": "layer:log1p"}, "names a layer this store does not have"),
        ({"expression_source": "raw.X"}, "expression_source must be"),
        ({"gene_id_column": "gene_symbols"}, r"var\['gene_symbols'\]"),
        ({"latent_key": "X_scVI"}, r"latent_key names obsm\['X_scVI'\]"),
        ({"vector_field_default": "velocity_umap"}, "does not declare"),
        ({"vector_field_default": "drift"}, "does not declare"),
        ({"obs_keys": ["missing"]}, r"obs_keys names obs\['missing'\]"),
        ({"obs_keys": ["score", "score"]}, "more than once"),
        ({"obs_keys": []}, "at least one obs column"),
        ({"obs_keys": "score"}, "sequence of obs column names"),
        ({"centroid_min_points": 0}, "centroid_min_points"),
        ({"normalize_embeddings": 1}, "normalize_embeddings"),
    ],
)
def test_a_launch_choice_the_store_cannot_honour_is_refused_before_binding(
    choices: dict[str, object],
    message: str,
) -> None:
    with pytest.raises((ValueError, TypeError), match=message):
        AnnDataServer(
            _adata(),
            port=0,
            quiet=True,
            serve_web_ui=False,
            dataset_name="Refused",
            dataset_id="refused",
            **choices,
        )


def test_a_default_vector_field_must_be_drawn_on_an_embedding_shown_at_open() -> None:
    """The web app selects only a field bound to an embedding it shows when the
    store opens: a named slot, the default rule's, or a 1-3 column array that
    is not the latent space. A wide basis, or the latent space, is not shown."""
    adata = _adata()
    rng = np.random.default_rng(11)
    adata.obsm["drift_pca"] = rng.random((12, 6)).astype(np.float32)
    adata.obsm["X_tsne"] = rng.random((12, 2)).astype(np.float32)
    adata.obsm["drift_tsne"] = rng.random((12, 2)).astype(np.float32)

    def defaults(**choices: object) -> dict[str, object]:
        with _serving(adata, **choices) as (_server, host, port):
            catalog = _json(host, port, "/_cellucid/datasets")
        return catalog["datasets"][0]["anndata"]["defaults"]

    # drift is drawn on X_tsne, a narrow array shown beside the default UMAP.
    assert defaults(vector_field_default="drift")["vectorFieldDefault"] == "drift"
    # With X_tsne as the latent space, drift's only shown basis is gone.
    with pytest.raises(ValueError, match="the viewer does not show"):
        AnnDataServer(
            adata,
            port=0,
            quiet=True,
            serve_web_ui=False,
            dataset_name="Refused",
            dataset_id="refused",
            latent_key="X_tsne",
            vector_field_default="drift",
        )
    # Showing the wide X_pca draws drift there again.
    shown = defaults(
        latent_key="X_tsne",
        embedding_keys={"2d": ("X_pca", [0, 1])},
        vector_field_default="drift",
    )
    assert shown["vectorFieldDefault"] == "drift"


def test_a_gene_id_column_is_looked_up_on_the_axis_of_the_chosen_matrix() -> None:
    adata = _adata()
    adata.raw = ad.AnnData(
        X=adata.X.copy(),
        var=pd.DataFrame({"raw_symbol": list("abcde")}, index=adata.var_names),
    )
    with pytest.raises(ValueError, match=r"raw.var\['symbol'\]"):
        AnnDataServer(adata, port=0, quiet=True, serve_web_ui=False, dataset_name="x",
                      dataset_id="x", expression_source="raw", gene_id_column="symbol")
    with _serving(adata, expression_source="raw", gene_id_column="raw_symbol") as (_s, host, port):
        defaults = _json(host, port, "/_cellucid/datasets")["datasets"][0]["anndata"]["defaults"]
    assert defaults["geneIdColumn"] == "raw_symbol"


def _described(adata: ad.AnnData) -> ad.AnnData:
    adata.uns["cellucid"] = {
        "schema": 1,
        "id": "declared-id",
        "name": "Declared name",
        "defaults": {"knn": "connectivities"},
    }
    return adata


@pytest.mark.parametrize("where", ["file", "object"])
def test_a_store_that_names_itself_needs_no_identity_at_launch(tmp_path: Path, where: str) -> None:
    adata = _described(_adata())
    data: object = adata
    if where == "file":
        data = tmp_path / "described.h5ad"
        adata.write_h5ad(data)
    server = AnnDataServer(data, port=0, quiet=True, serve_web_ui=False)
    try:
        assert (server.dataset_id, server.dataset_name) == ("declared-id", "Declared name")
        # Served in priority order: a Zarr store named by the declared id.
        assert server.file_route == "_cellucid/anndata/declared-id.zarr/"
        # What the store declares is read by the viewer; launch options override it.
        assert server.defaults == {}
    finally:
        server.stop()
    named = AnnDataServer(data, port=0, quiet=True, serve_web_ui=False,
                          dataset_id="given-id", dataset_name="Given name", knn_key="connectivities")
    try:
        assert (named.dataset_id, named.dataset_name) == ("given-id", "Given name")
        assert named.defaults == {"knn": "connectivities"}
    finally:
        named.stop()


def test_a_store_without_identity_must_be_named(tmp_path: Path) -> None:
    source = tmp_path / "anonymous.h5ad"
    _adata().write_h5ad(source)
    with pytest.raises(ValueError, match=r"dataset_name and dataset_id are required.*uns\['cellucid'\]"):
        AnnDataServer(source, port=0, quiet=True, serve_web_ui=False)
    with pytest.raises(ValueError, match="dataset_id is required"):
        AnnDataServer(_adata(), port=0, quiet=True, serve_web_ui=False, dataset_name="Only a name")
    adata = _adata()
    adata.uns["cellucid"] = {"schema": 1, "id": "bad/id", "name": "Name"}
    with pytest.raises(ValueError, match=r"uns\['cellucid'\]\['id'\]"):
        AnnDataServer(adata, port=0, quiet=True, serve_web_ui=False)


# =============================================================================
# The bytes of a file
# =============================================================================


def test_a_file_is_read_with_exact_ranges(tmp_path: Path) -> None:
    source = tmp_path / "x.h5ad"
    _in_priority_order(_adata()).write_h5ad(source)
    data = source.read_bytes()
    size = len(data)
    route = "/_cellucid/anndata/store-contract.h5ad"
    with _serving(source) as (_server, host, port):
        status, headers, body = _request(host, port, "GET", route, {"Range": "bytes=8-23"})
        assert status == 206
        assert body == data[8:24]
        assert headers["Content-Range"] == f"bytes 8-23/{size}"
        assert headers["Content-Length"] == "16"
        assert headers["Accept-Ranges"] == "bytes"
        assert headers["Content-Type"] == "application/octet-stream"
        etag = headers["ETag"]
        assert etag.startswith('"') and etag.endswith('"')
        expose = {name.strip() for name in headers["Access-Control-Expose-Headers"].split(",")}
        assert {"Accept-Ranges", "Content-Length", "Content-Range", "ETag"} <= expose

        # The probe the viewer makes: HEAD for the size, then one byte.
        status, headers, body = _request(host, port, "HEAD", route)
        assert (status, body, headers["Content-Length"], headers["ETag"]) == (200, b"", str(size), etag)
        status, headers, body = _request(host, port, "GET", route, {"Range": "bytes=0-0"})
        assert (status, body, headers["Content-Range"]) == (206, data[:1], f"bytes 0-0/{size}")

        # An open end and a suffix; a range past the end is clipped.
        assert _request(host, port, "GET", route, {"Range": f"bytes={size - 4}-"})[2] == data[-4:]
        assert _request(host, port, "GET", route, {"Range": "bytes=-5"})[2] == data[-5:]
        assert _request(host, port, "GET", route, {"Range": f"bytes=0-{size + 99}"})[2] == data

        # No transfer encoding, even when asked: a range is of the stored bytes.
        status, headers, body = _request(host, port, "GET", route, {"Accept-Encoding": "gzip"})
        assert (status, body) == (200, data)
        assert "Content-Encoding" not in headers

        for refused in (f"bytes={size}-", "bytes=0-1,4-5", "bytes=5-2", "items=0-1", "bytes=-0"):
            status, headers, _ = _request(host, port, "GET", route, {"Range": refused})
            assert status == 416, refused
            assert headers["Content-Range"] == f"bytes */{size}"

        for absent in ("/_cellucid/anndata/x.h5ad", "/_cellucid/anndata/store-contract.h5ad/x",
                       "/store-contract.h5ad", "/_cellucid/anndata/"):
            assert _request(host, port, "GET", absent)[0] == 404, absent


def test_a_file_that_changes_on_disk_names_its_new_version(tmp_path: Path) -> None:
    """The dataset list states the validator the file's bytes are answered
    with -- the version the viewer files what it keeps of the file under
    (its byte cache, a kept gene index, a converted store) -- read from the
    file when the list is asked for: grown, rewritten in place with its size
    and modification time kept, or replaced by a copy moved over it, the file
    is a new version."""
    source = tmp_path / "x.h5ad"
    _in_priority_order(_adata()).write_h5ad(source)
    route = "/_cellucid/anndata/store-contract.h5ad"

    def listed() -> dict:
        return _json(host, port, "/_cellucid/datasets")["datasets"][0]["anndata"]["file"]

    with _serving(source) as (_server, host, port):
        first = listed()
        before = _request(host, port, "HEAD", route)[1]["ETag"]
        assert first["etag"] == before == path_version(source)[1]
        assert first["size"] == source.stat().st_size
        assert listed() == first
        with source.open("ab") as file:
            file.write(b"\0")
        after = _request(host, port, "HEAD", route)[1]
        assert after["ETag"] != before
        assert after["Content-Length"] == str(source.stat().st_size)
        assert listed() == {**first, "size": source.stat().st_size, "etag": after["ETag"]}

        # In place, the same size, the modification time put back.
        grown = listed()["etag"]
        status = source.stat()
        with source.open("r+b") as file:
            head = file.read(8)
            file.seek(0)
            file.write(head[::-1])
        os.utime(source, ns=(status.st_atime_ns, status.st_mtime_ns))
        assert source.stat().st_size == status.st_size and source.stat().st_mtime_ns == status.st_mtime_ns
        rewritten = listed()["etag"]
        assert rewritten != grown
        assert _request(host, port, "HEAD", route)[1]["ETag"] == rewritten

        # Replaced by a copy moved over it, its times copied with it.
        copy = tmp_path / "copy.h5ad"
        shutil.copy2(source, copy)
        os.replace(copy, source)
        assert listed()["etag"] not in (grown, rewritten)

        # Gone: not listed as if it were there.
        source.unlink()
        status, _, body = _request(host, port, "GET", "/_cellucid/datasets")
        assert (status, body) == (404, b"x.h5ad is no longer readable")


def test_the_listed_version_is_read_from_the_open_file_as_its_bytes_are(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The validator the dataset list states and the one the file's bytes are
    answered with are read the same way, from the file opened: on Windows a
    path's ``os.stat`` reports the creation time as ``st_ctime`` and an open
    file's ``os.fstat`` the change time (Python 3.12 on), so a tag read from
    the path would name another version than the bytes' tag, and the viewer
    would refuse the file. Here a path's status reports another status-change
    time, as Windows' does."""
    source = tmp_path / "x.h5ad"
    _in_priority_order(_adata()).write_h5ad(source)
    path_status = Path.stat

    def as_windows(self: Path, *args: object, **kwargs: object) -> os.stat_result:
        status = path_status(self, *args, **kwargs)  # type: ignore[arg-type]
        if self.name != source.name:
            return status
        fields = list(status)
        return os.stat_result(fields, {"st_ctime_ns": status.st_ctime_ns - 1, "st_mtime_ns": status.st_mtime_ns,
                                       "st_atime_ns": status.st_atime_ns})

    with _serving(source) as (_server, host, port):
        monkeypatch.setattr(Path, "stat", as_windows)
        assert source.stat().st_ctime_ns != os.stat(source).st_ctime_ns
        listed = _json(host, port, "/_cellucid/datasets")["datasets"][0]["anndata"]["file"]
        route = "/_cellucid/anndata/store-contract.h5ad"
        assert listed["etag"] == _request(host, port, "HEAD", route)[1]["ETag"]
        assert listed["etag"] == _request(host, port, "GET", route, {"Range": "bytes=0-7"})[1]["ETag"]
        assert listed["etag"] == path_version(source)[1]


def test_a_browser_reuses_the_bytes_it_kept_of_an_unchanged_file(tmp_path: Path) -> None:
    """Store bytes are cacheable by validation: a reload revalidates (304) and
    transfers nothing it has; a changed file is never mixed with kept bytes."""
    source = tmp_path / "x.h5ad"
    _in_priority_order(_adata()).write_h5ad(source)
    data = source.read_bytes()
    route = "/_cellucid/anndata/store-contract.h5ad"
    with _serving(source) as (_server, host, port):
        status, headers, _ = _request(host, port, "GET", route, {"Range": "bytes=0-15"})
        etag = headers["ETag"]
        assert headers["Cache-Control"] == "no-cache"
        # Revalidation of kept bytes: nothing is sent again.
        for method in ("GET", "HEAD"):
            status, headers, body = _request(
                host, port, method, route, {"Range": "bytes=0-15", "If-None-Match": etag}
            )
            assert (status, body, headers["ETag"]) == (304, b"", etag)
        assert _request(host, port, "GET", route, {"If-None-Match": '"other", ' + etag})[0] == 304
        # Completing kept bytes of the same version: the range.
        status, _, body = _request(host, port, "GET", route, {"Range": "bytes=16-31", "If-Range": etag})
        assert (status, body) == (206, data[16:32])
        with source.open("ab") as file:
            file.write(b"\0")
        changed = source.read_bytes()
        # The kept version is no longer current: the whole new file, never a range of it.
        status, headers, body = _request(host, port, "GET", route, {"Range": "bytes=16-31", "If-Range": etag})
        assert (status, body) == (200, changed)
        assert headers["ETag"] != etag
        assert _request(host, port, "GET", route, {"If-None-Match": etag})[0] == 200


def test_file_bytes_are_sent_by_the_kernel_without_delaying_small_answers() -> None:
    from cellucid.anndata_server import AnnDataRequestHandler

    assert AnnDataRequestHandler.disable_nagle_algorithm is True
    assert AnnDataRequestHandler.protocol_version == "HTTP/1.1"


def test_connections_are_kept_alive_and_answered_concurrently(tmp_path: Path) -> None:
    source = tmp_path / "x.h5ad"
    _in_priority_order(_adata()).write_h5ad(source)
    data = source.read_bytes()
    route = "/_cellucid/anndata/store-contract.h5ad"
    with _serving(source) as (_server, host, port):
        idle = http.client.HTTPConnection(host, port, timeout=10)
        idle.request("GET", route, headers={"Range": "bytes=0-3"})
        first = idle.getresponse()
        assert first.version == 11
        assert first.read() == data[:4]
        # The same connection serves the next request.
        idle.request("GET", route, headers={"Range": "bytes=4-7"})
        assert idle.getresponse().read() == data[4:8]
        # A kept-alive connection holds nothing up: others are answered at once.
        results: list[bytes] = []
        threads = [
            threading.Thread(
                target=lambda i=i: results.append(
                    _request(host, port, "GET", route, {"Range": f"bytes={i}-{i}"})[2]
                )
            )
            for i in range(8)
        ]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=10)
        idle.close()
    assert sorted(results) == sorted(data[i : i + 1] for i in range(8))


# =============================================================================
# Zarr stores
# =============================================================================


def _zarr_folder(tmp_path: Path, *, consolidated: bool) -> Path:
    folder = tmp_path / "local.zarr"
    _in_priority_order(_adata()).write_zarr(folder)
    if not consolidated:
        (folder / ".zmetadata").unlink()
    return folder


def test_a_zarr_folder_is_read_key_by_key_with_404_for_an_absent_key(tmp_path: Path) -> None:
    folder = _zarr_folder(tmp_path, consolidated=True)
    outside = tmp_path / "secret.txt"
    outside.write_text("not part of the store")
    (folder / "obs" / "escape").symlink_to(outside)
    route = "/_cellucid/anndata/store-contract.zarr"
    with _serving(folder) as (server, host, port):
        catalog = _json(host, port, "/_cellucid/datasets")["datasets"][0]["anndata"]["file"]
        assert catalog == {"url": "_cellucid/anndata/store-contract.zarr/", "kind": "zarr", "layout": "folder"}
        assert server._served.generated_metadata is None
        status, _, body = _request(host, port, "GET", f"{route}/.zmetadata")
        assert (status, body) == (200, (folder / ".zmetadata").read_bytes())
        chunk = next(path for path in (folder / "X" / "data").iterdir() if not path.name.startswith("."))
        key = chunk.relative_to(folder).as_posix()
        status, headers, body = _request(host, port, "GET", f"{route}/{key}", {"Range": "bytes=0-1"})
        assert (status, body, headers["Content-Range"]) == (206, chunk.read_bytes()[:2],
                                                            f"bytes 0-1/{chunk.stat().st_size}")
        for absent in ("X/data/99", "obs/escape", "../secret.txt", "obs/../../secret.txt", "", "obs",
                       "obs//cell_type"):
            assert _request(host, port, "GET", f"{route}/{absent}")[0] == 404, absent
        assert _request(host, port, "GET", route)[0] == 404


def test_a_zarr_v2_folder_without_consolidated_metadata_gets_it_built(tmp_path: Path) -> None:
    folder = _zarr_folder(tmp_path, consolidated=False)
    with _serving(folder) as (_server, host, port):
        status, headers, body = _request(host, port, "GET", "/_cellucid/anndata/store-contract.zarr/.zmetadata")
        assert status == 200
        assert headers["Content-Length"] == str(len(body))
        assert "ETag" in headers
    document = json.loads(body)
    assert document["zarr_consolidated_format"] == 1
    on_disk = {
        path.relative_to(folder).as_posix(): json.loads(path.read_text())
        for path in folder.rglob("*")
        if path.name in {".zgroup", ".zattrs", ".zarray"}
    }
    assert document["metadata"] == on_disk


def test_a_zarr_v3_folder_without_consolidated_metadata_gets_it_inline(tmp_path: Path) -> None:
    folder = tmp_path / "v3.zarr"
    with ad.settings.override(zarr_write_format=3):
        _in_priority_order(_adata()).write_zarr(folder)
    root = json.loads((folder / "zarr.json").read_text())
    root.pop("consolidated_metadata", None)
    (folder / "zarr.json").write_text(json.dumps(root))
    with _serving(folder) as (_server, host, port):
        served = json.loads(_request(host, port, "GET", "/_cellucid/anndata/store-contract.zarr/zarr.json")[2])
    nodes = served["consolidated_metadata"]["metadata"]
    assert served["consolidated_metadata"]["kind"] == "inline"
    assert {key: value for key, value in served.items() if key != "consolidated_metadata"} == root
    assert "X/data" in nodes and "obs" in nodes and "" not in nodes
    assert nodes["X/data"] == json.loads((folder / "X" / "data" / "zarr.json").read_text())


def test_a_zipped_zarr_store_is_one_file(tmp_path: Path) -> None:
    _zarr_folder(tmp_path, consolidated=True)
    archive = tmp_path / "local.zarr.zip"
    shutil.make_archive(str(archive.with_suffix("")), "zip", root_dir=tmp_path, base_dir="local.zarr")
    store = classify_anndata_path(archive)
    assert (store.kind, store.layout, store.zip_root) == ("zarr", "zip", "local.zarr")
    with _serving(archive) as (_server, host, port):
        file = _json(host, port, "/_cellucid/datasets")["datasets"][0]["anndata"]["file"]
        assert file == {"url": "_cellucid/anndata/store-contract.zarr.zip", "kind": "zarr",
                        "layout": "zip", "size": archive.stat().st_size, "etag": path_version(archive)[1]}
        status, _, body = _request(host, port, "GET", "/_cellucid/anndata/store-contract.zarr.zip",
                                   {"Range": "bytes=0-3"})
    assert (status, body) == (206, b"PK\x03\x04")
    assert read_anndata_contents(store).n_obs == 12


def test_a_zip_that_holds_no_zarr_store_is_refused(tmp_path: Path) -> None:
    archive = tmp_path / "notes.zip"
    with zipfile.ZipFile(archive, "w") as file:
        file.writestr("notes.json", "{}")
    with pytest.raises(ValueError, match="holds no Zarr store"):
        classify_anndata_path(archive)


# =============================================================================
# Objects, and stores the browser cannot decode
# =============================================================================


def test_an_hdf5_filter_the_browser_cannot_undo_is_decoded_on_the_kernel(tmp_path: Path) -> None:
    """A file in priority order with an array the browser cannot decode is
    served through the store in priority order, in the file's own order."""
    source = tmp_path / "scaleoffset.h5ad"
    placed = _in_priority_order(_adata())
    placed.write_h5ad(source)
    with h5py.File(source, "r+") as file:
        values = file["obsm/X_pca"][...]
        attrs = dict(file["obsm/X_pca"].attrs)
        del file["obsm/X_pca"]
        dataset = file["obsm"].create_dataset("X_pca", data=values, scaleoffset=6, chunks=(4, 6))
        dataset.attrs.update(attrs)
    with h5py.File(source, "r") as file:
        values = file["obsm/X_pca"][...]
    store = classify_anndata_path(source)
    assert browser_undecodable(store) == ["obsm/X_pca (scaleoffset)"]
    with _serving(source) as (server, host, port):
        assert server._served.priority is not None
        assert _json(host, port, "/_cellucid/datasets")["datasets"][0]["anndata"]["file"] == {
            "url": "_cellucid/anndata/store-contract.zarr/", "kind": "zarr", "layout": "folder",
        }
        metadata = json.loads(_request(host, port, "GET", "/_cellucid/anndata/store-contract.zarr/obsm/X_pca/.zarray")[2])
        assert (metadata["compressor"], metadata["filters"]) == (None, None)
        chunk = _request(host, port, "GET", "/_cellucid/anndata/store-contract.zarr/obsm/X_pca/0.0")[2]
        rows = metadata["chunks"][0]
        decoded = np.frombuffer(chunk, dtype=metadata["dtype"]).reshape(metadata["chunks"])[: placed.n_obs]
        np.testing.assert_array_equal(decoded, values[:rows])
        # The file's own order and records, as stored.
        positions = _request(host, port, "GET", f"/_cellucid/anndata/store-contract.zarr/obs/{ORIGINAL_POSITION_COLUMN}/0")[2]
        np.testing.assert_array_equal(
            np.frombuffer(positions, dtype=np.uint32)[: placed.n_obs], placed.obs[ORIGINAL_POSITION_COLUMN]
        )
    assert source.exists()


def test_zarr_codecs_the_browser_cannot_decode_are_named(tmp_path: Path) -> None:
    folder = _zarr_folder(tmp_path, consolidated=True)
    zarray = folder / "obsm" / "X_pca" / ".zarray"
    meta = json.loads(zarray.read_text())
    meta["compressor"] = {"id": "bz2", "level": 1}
    meta["filters"] = [{"id": "delta", "dtype": "<f4"}]
    zarray.write_text(json.dumps(meta))
    assert browser_undecodable(classify_anndata_path(folder)) == ["obsm/X_pca (bz2, delta)"]


def test_a_browser_decodable_store_needs_no_copy(tmp_path: Path) -> None:
    source = tmp_path / "plain.h5ad"
    adata = _adata()
    adata.write_h5ad(source, compression="gzip")
    assert browser_undecodable(classify_anndata_path(source)) == []
    lzf = tmp_path / "lzf.h5ad"
    adata.write_h5ad(lzf, compression="lzf")
    assert browser_undecodable(classify_anndata_path(lzf)) == []
    folder = _zarr_folder(tmp_path, consolidated=True)
    assert browser_undecodable(classify_anndata_path(folder)) == []


# =============================================================================
# What a store holds, from metadata alone
# =============================================================================


@pytest.mark.parametrize(
    "name",
    ["legacy-0.6.h5ad", "legacy-0.6-compound-obsm.h5ad", "legacy-0.7-csr.h5ad",
     "legacy-0.7-dense.h5ad", "legacy-h5sparse.h5ad", "legacy-0.7.zarr.zip"],
)
def test_older_layouts_are_described_like_their_current_twin(name: str) -> None:
    contents = read_anndata_contents(classify_anndata_path(FIXTURES / name))
    twin = read_anndata_contents(classify_anndata_path(FIXTURES / "legacy-twin-current.h5ad"))
    assert (contents.n_obs, contents.n_vars) == (twin.n_obs, twin.n_vars)
    assert contents.obs_columns == twin.obs_columns
    assert contents.obsm == twin.obsm
    assert contents.has_x
    assert contents.raw_var_columns == twin.raw_var_columns


def test_columns_seuratdisk_leaves_out_of_column_order_are_still_named() -> None:
    contents = read_anndata_contents(classify_anndata_path(FIXTURES / "legacy-seuratdisk.h5ad"))
    with h5py.File(FIXTURES / "legacy-seuratdisk.h5ad", "r") as file:
        listed = set(file["obs"].attrs["column-order"])
    assert listed < contents.obs_columns
    assert "batch" in contents.obs_columns


# =============================================================================
# Diagnostics read the same catalog
# =============================================================================


@pytest.mark.parametrize("layout", ["file", "folder"])
def test_debug_connection_accepts_the_catalog_and_probes_the_store(tmp_path: Path, layout: str) -> None:
    from cellucid.jupyter._viewer import _probe_served_store
    from cellucid.jupyter._wire import _require_debug_dataset_list

    data: object = _adata() if layout == "file" else _zarr_folder(tmp_path, consolidated=True)
    with _serving(data) as (server, host, port):
        catalog = _json(host, port, "/_cellucid/datasets")
        assert _require_debug_dataset_list(catalog) == catalog["datasets"]
        probe = _probe_served_store(f"http://{host}:{port}/{server.file_route}")
    assert "error" not in probe, probe
    assert probe["head"]["accept_ranges"] == "bytes"
    assert probe["range"]["status"] == 206
    assert probe["range"]["content_range"].startswith("bytes 0-0/")
