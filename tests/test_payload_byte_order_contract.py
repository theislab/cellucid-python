"""Every published payload is little-endian, whatever the host's byte order is.

This is the one format rule that no reader can check. The dtype names the
manifests publish -- ``float32``, ``uint16``, ``float64`` -- carry no
byte-order component, and the web app builds its typed arrays straight from the
bytes it receives, which is host order by definition in JavaScript. A payload
written in the wrong order therefore decodes into plausible coordinates and
plausible expression values that are simply not the ones the producer had, and
nothing anywhere in the chain can raise.

These tests have to be written against that fact. A test that writes an array
on a little-endian machine and asserts the bytes are little-endian proves
nothing: it passes just as happily when the writer emits native order. Each
test below therefore presents the writer with a **big-endian** array holding
the intended values -- the condition every array is in on a big-endian host --
and requires the published bytes to be unchanged.
"""

from __future__ import annotations

import ast
import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellucid import anndata_adapter as adapter_module
from cellucid._byte_order import (
    PAYLOAD_BYTE_ORDER,
    little_endian_payload,
    little_endian_payload_bytes,
    little_endian_payload_dtype,
    little_endian_payload_view,
)
from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.prepare_data import _generation as prepare_generation
from cellucid.prepare_data import prepare
from cellucid.prepare_data._binary_io import _write_binary

SOURCE_ROOT = Path(__file__).resolve().parents[1] / "src" / "cellucid"
BYTE_ORDER_MODULE = SOURCE_ROOT / "_byte_order.py"

# Every scalar type the export format publishes, with the dtype character the
# array-protocol string uses for it.
PUBLISHED_DTYPES = ["float32", "float64", "uint8", "uint16", "uint32"]

N_CELLS = 64
FIXED_CREATED_AT = "2026-01-01T00:00:00Z"


def _as_big_endian(array: np.ndarray) -> np.ndarray:
    """Return the same values in the byte order a big-endian host would hold.

    Single-byte types are returned unchanged, exactly as they would be on such
    a host: a ``uint8`` payload has no byte order to get wrong.
    """
    return array.astype(array.dtype.newbyteorder(">"))


def _reference_little_endian_bytes(values: list[float], dtype: str) -> bytes:
    """Build the expected bytes independently of the code under test."""
    return np.array(values, dtype=np.dtype(dtype).newbyteorder("<")).tobytes()


# ---------------------------------------------------------------------------
# The conversion itself
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("dtype", PUBLISHED_DTYPES)
def test_a_big_endian_array_is_converted_to_the_published_order(dtype: str) -> None:
    values = [1.0, 2.0, 3.0, 4.0]
    big_endian = np.array(values, dtype=np.dtype(dtype).newbyteorder(">"))

    published = little_endian_payload_bytes(big_endian)

    assert published == _reference_little_endian_bytes(values, dtype)
    # Positive behaviour, not merely a byte comparison: the reader's own
    # interpretation of those bytes recovers the values that went in.
    decoded = np.frombuffer(published, dtype=np.dtype(dtype).newbyteorder("<"))
    assert decoded.tolist() == values


@pytest.mark.parametrize("dtype", ["float32", "float64", "uint16", "uint32"])
def test_a_big_endian_payload_is_not_published_verbatim(dtype: str) -> None:
    """The fixture is only evidence if the two orders actually differ."""
    big_endian = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.dtype(dtype).newbyteorder(">"))

    assert little_endian_payload_bytes(big_endian) != big_endian.tobytes()


def test_a_single_byte_payload_has_no_byte_order_to_convert() -> None:
    codes = np.arange(8, dtype=np.uint8)

    assert little_endian_payload(codes) is codes
    assert little_endian_payload_bytes(codes) == codes.tobytes()


@pytest.mark.parametrize("dtype", PUBLISHED_DTYPES)
def test_the_published_storage_dtype_is_always_little_endian(dtype: str) -> None:
    storage = little_endian_payload_dtype(np.dtype(dtype).newbyteorder(">"))

    assert storage.str == np.dtype(dtype).newbyteorder("<").str
    assert storage.str[0] in {"<", "|"}


def test_an_already_published_payload_is_borrowed_rather_than_copied() -> None:
    """Conversion must cost nothing on the hosts anyone actually exports from."""
    contiguous = np.arange(12, dtype="<f4").reshape(3, 4)
    contiguous_view = little_endian_payload_view(contiguous)

    assert contiguous_view.obj is contiguous
    assert contiguous_view.tobytes() == contiguous.tobytes(order="C")


def test_a_fortran_ordered_payload_is_published_row_major() -> None:
    fortran_ordered = np.asfortranarray(np.arange(12, dtype="<f4").reshape(3, 4))
    fortran_view = little_endian_payload_view(fortran_ordered)

    assert fortran_view.obj is not fortran_ordered
    assert fortran_view.tobytes() == fortran_ordered.tobytes(order="C")


def test_a_big_endian_fortran_payload_is_published_row_major_and_little_endian() -> None:
    source = np.arange(12, dtype="<f4").reshape(3, 4)
    hostile = np.asfortranarray(_as_big_endian(source))

    assert little_endian_payload_view(hostile).tobytes() == source.tobytes(order="C")


# ---------------------------------------------------------------------------
# The export writer
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("compression", [None, 6])
@pytest.mark.parametrize("dtype", PUBLISHED_DTYPES)
def test_the_writer_publishes_a_big_endian_array_little_endian(
    tmp_path: Path,
    dtype: str,
    compression: int | None,
) -> None:
    values = [0.0, 1.0, 2.0, 3.0, 4.0]
    big_endian = np.array(values, dtype=np.dtype(dtype).newbyteorder(">"))

    written = _write_binary(tmp_path / "payload.bin", big_endian, compression)

    stored = written.read_bytes()
    published = gzip.decompress(stored) if compression is not None else stored
    assert published == _reference_little_endian_bytes(values, dtype)


def test_a_big_endian_host_publishes_a_byte_identical_export(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The whole export, produced twice, once as a big-endian host would.

    Every array reaching the writer is byte-swapped to big-endian while
    keeping its values, which is exactly the state every array is in on such a
    host. The two export directories must be identical file for file. With the
    writer emitting native order this fails on every multi-byte payload:
    points, quantized codes, outlier quantiles, gene values, and both
    connectivity endpoint arrays.
    """
    native_dir = tmp_path / "native"
    big_endian_dir = tmp_path / "big-endian"
    _prepare_full_export(native_dir)

    real_write_binary = prepare_generation._write_binary
    swapped_payloads: list[str] = []

    def write_from_a_big_endian_host(path, data, compression=None):
        array = np.asarray(data)
        swapped_payloads.append(Path(path).name)
        return real_write_binary(path, _as_big_endian(array), compression)

    monkeypatch.setattr(prepare_generation, "_write_binary", write_from_a_big_endian_host)
    _prepare_full_export(big_endian_dir)

    native_files = sorted(p.relative_to(native_dir) for p in native_dir.rglob("*") if p.is_file())
    big_endian_files = sorted(
        p.relative_to(big_endian_dir) for p in big_endian_dir.rglob("*") if p.is_file()
    )
    assert native_files == big_endian_files

    # The simulation must actually have reached every payload family, or the
    # comparison above proves nothing about the ones it missed.
    assert {name.split(".")[0] for name in swapped_payloads} >= {
        "points_1d",
        "points_2d",
        "points_3d",
        "0",
        "1",
    }
    assert "edges" in {name.split(".")[0] for name in swapped_payloads}

    differing = [
        name
        for name in native_files
        if (native_dir / name).read_bytes() != (big_endian_dir / name).read_bytes()
    ]
    assert differing == []


def _prepare_full_export(out_dir: Path) -> None:
    """One export exercising every binary payload the format defines."""
    rng = np.random.default_rng(20260801)
    embedding = rng.normal(size=(N_CELLS, 3))
    codes = np.arange(N_CELLS) % 4
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical.from_codes(
                codes,
                categories=[f"Group {index}" for index in range(4)],
            ),
            "score": rng.normal(size=N_CELLS),
        },
        index=[f"cell-{index}" for index in range(N_CELLS)],
    )
    connectivities = np.zeros((N_CELLS, N_CELLS), dtype=np.float64)
    for index in range(N_CELLS - 1):
        connectivities[index, index + 1] = 0.25 + index / 1000.0
        connectivities[index + 1, index] = 0.25 + index / 1000.0

    prepare(
        latent_space=embedding.copy(),
        obs=obs,
        X_umap_1d=embedding[:, :1],
        X_umap_2d=embedding[:, :2],
        X_umap_3d=embedding,
        var=pd.DataFrame(index=pd.Index(["gene_a", "gene_b", "gene_c"], dtype=str)),
        gene_expression=rng.random((N_CELLS, 3), dtype=np.float64),
        connectivities=connectivities,
        vector_fields={"velocity_umap_2d": rng.normal(size=(N_CELLS, 2))},
        out_dir=out_dir,
        dataset_id="byte-order-contract",
        dataset_name="Byte order contract",
        obs_categorical_dtype="uint16",
        obs_continuous_quantization=16,
        var_quantization=8,
        centroid_min_points=4,
        # Two exports are compared byte for byte, so they must not carry two
        # different timestamps. Left to default, dataset_identity.json differs
        # whenever the wall clock ticks a second between the two calls, which
        # failed about one run in six and said nothing about byte order.
        created_at=FIXED_CREATED_AT,
    )


# ---------------------------------------------------------------------------
# The direct AnnData server
# ---------------------------------------------------------------------------


def _served_adapter() -> AnnDataAdapter:
    rng = np.random.default_rng(4242)
    codes = np.arange(N_CELLS) % 4
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical.from_codes(
                codes,
                categories=[f"Group {index}" for index in range(4)],
            ),
            "score": rng.normal(size=N_CELLS).astype(np.float32),
        },
        index=[f"cell-{index}" for index in range(N_CELLS)],
    )
    adata = ad.AnnData(
        X=rng.random((N_CELLS, 2), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(["gene_a", "gene_b"], dtype=str)),
    )
    adata.obsm["X_umap_2d"] = rng.normal(size=(N_CELLS, 2))
    adata.obsm["velocity_umap_2d"] = rng.normal(size=(N_CELLS, 2))
    adata.obsm["X_latent"] = rng.normal(size=(N_CELLS, 2))
    connectivities = np.zeros((N_CELLS, N_CELLS), dtype=np.float64)
    for index in range(N_CELLS - 1):
        connectivities[index, index + 1] = 0.5
        connectivities[index + 1, index] = 0.5
    adata.obsp["connectivities"] = connectivities
    return AnnDataAdapter(
        adata,
        latent_key="X_latent",
        centroid_min_points=4,
        dataset_name="Byte order contract",
        dataset_id="byte-order-contract",
    )


def _served_payloads(adapter: AnnDataAdapter) -> dict[str, object]:
    return {
        "points": adapter.get_points_binary(2),
        "vectors": adapter.get_vector_field_binary("velocity_umap", 2),
        "continuous": adapter.get_obs_continuous_values("score"),
        "codes": adapter.get_obs_categorical_codes("cell_type"),
        "outliers": adapter.get_obs_outlier_quantiles("cell_type"),
        "expression": adapter.get_gene_expression("gene_a"),
        "connectivity": adapter.get_connectivity_edges(),
    }


def test_a_big_endian_host_serves_byte_identical_payloads(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The served payloads, produced twice, once as a big-endian host would.

    The count assertion is what makes this discriminating: it fails if any
    emission point stops routing through the byte-order rule, which is the
    only way a served payload could carry the host's order instead of the
    format's.
    """
    adapter = _served_adapter()
    try:
        native = _served_payloads(adapter)

        real_payload_bytes = adapter_module.little_endian_payload_bytes
        conversions: list[int] = []

        def serve_from_a_big_endian_host(data: np.ndarray) -> bytes:
            array = np.asarray(data)
            conversions.append(array.size)
            return real_payload_bytes(_as_big_endian(array))

        monkeypatch.setattr(
            adapter_module,
            "little_endian_payload_bytes",
            serve_from_a_big_endian_host,
        )

        big_endian_adapter = _served_adapter()
        try:
            big_endian = _served_payloads(big_endian_adapter)
        finally:
            big_endian_adapter.close()
    finally:
        adapter.close()

    assert big_endian == native
    # Seven getters, and the connectivity one emits three aligned arrays.
    assert len(conversions) == 9


# ---------------------------------------------------------------------------
# No second way to emit payload bytes
# ---------------------------------------------------------------------------


def _array_to_bytes_conversions(source: Path) -> list[str]:
    """Find every ``ndarray`` -> bytes conversion in one module."""
    tree = ast.parse(source.read_text(encoding="utf-8"), filename=str(source))
    found: list[str] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        if node.func.attr in {"tobytes", "tofile"}:
            found.append(f"{source.name}:{node.lineno}: .{node.func.attr}()")
    return found


def test_no_module_converts_an_array_to_bytes_outside_the_byte_order_rule() -> None:
    """A new payload path must not be able to reintroduce host byte order.

    The defect this guards against is not a wrong line, it is a *missing* one:
    ``tobytes()`` and ``tofile()`` publish native order and say nothing about
    it, so a future write path that calls either of them directly is wrong on
    a big-endian host and silent everywhere else.
    """
    offenders: list[str] = []
    for source in sorted(SOURCE_ROOT.rglob("*.py")):
        if source == BYTE_ORDER_MODULE:
            continue
        offenders.extend(_array_to_bytes_conversions(source))

    assert offenders == []


def test_the_published_byte_order_is_defined_in_exactly_one_module() -> None:
    """One definition of '<', so no call site can quietly pick a different one."""
    assert PAYLOAD_BYTE_ORDER == "<"

    sources_naming_the_primitive = sorted(
        source.relative_to(SOURCE_ROOT).as_posix()
        for source in SOURCE_ROOT.rglob("*.py")
        if "newbyteorder" in source.read_text(encoding="utf-8")
    )

    assert sources_naming_the_primitive == [BYTE_ORDER_MODULE.name]
