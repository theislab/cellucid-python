"""Exact contract for payload naming: a filename never carries dataset content.

Every scientific payload an export writes is named by its integer index or by a
fixed neutral name, and the manifest is the only place a field's identity lives.
That removes a whole class of defect at the source: a gene called ``HLA-DRB1/2``
or ``CON``, an obs column called ``% mito``, and two fields that differ only in
case all used to fail an export or collide at a path, because the identifier had
to survive being a filename on three operating systems.

What replaces the filename rule is an index rule, and it is stricter in the one
way that matters. Within an axis directory the indices are exactly ``0..N-1``,
each used once, and ``obs`` shares that one space across ``_continuousFields``
and ``_categoricalFields`` because both write into ``obs/``. Two fields holding
one index write into one file: the second silently overwrites the first, and the
viewer then draws one field's values under the other field's name. Nothing
downstream can detect that, so the writer asserts it about the manifest it has
just built, and re-derives every declared payload path to compare it against the
directory it actually wrote.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cellucid.prepare_data import _generation, prepare
from cellucid.prepare_data._binary_io import _write_binary
from cellucid.prepare_data._manifest import (
    _gene_names_from_compact_manifest,
    _identity_obs_fields_from_compact_manifest,
    _require_declared_payloads_on_disk,
    _require_dense_payload_indices,
)
from cellucid.server._artifacts import _declared_dataset_artifacts

# Identifiers that no filesystem could carry, and that no longer have to.
HOSTILE_OBS_KEYS = ["% mito", "cell type", "CON", "Field", "field"]
HOSTILE_GENE_NAMES = ["HLA-DRB1/2", "NUL.txt", "trailing.", "細胞", "Gene", "gene"]

INDEXED_PAYLOAD = re.compile(r"^(0|[1-9][0-9]*)\.(values|codes|outliers)\.(f32|u8|u16)(\.gz)?$")
INDEXED_VECTOR = re.compile(r"^(0|[1-9][0-9]*)_[123]d\.bin(\.gz)?$")


def _hostile_export(out_dir: Path, **overrides: object) -> None:
    n_cells = 6
    coordinates = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 0.25], [0.25, 0.75]],
        dtype=np.float32,
    )
    # Distinct within-category distances, so the generated outlier quantiles
    # vary and the quantized variant of this export is encodable.
    latent = np.array(
        [[0.0, 0.0], [3.0, 0.0], [0.5, 0.0], [7.0, 0.0], [1.5, 0.0], [9.0, 0.0]],
        dtype=np.float32,
    )
    connectivities = np.zeros((n_cells, n_cells), dtype=np.float64)
    for source, destination in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5)):
        connectivities[source, destination] = 0.5
        connectivities[destination, source] = 0.5

    obs = pd.DataFrame(
        {
            "% mito": np.linspace(0.0, 1.0, n_cells, dtype=np.float32),
            "cell type": pd.Categorical(["a", "a", "a", "b", "b", "b"]),
            "CON": pd.Categorical(["x", "y", "y", "y", "x", "x"]),
            "Field": np.linspace(1.0, 2.0, n_cells, dtype=np.float32),
            "field": np.linspace(3.0, 4.0, n_cells, dtype=np.float32),
        }
    )
    assert list(obs.columns) == HOSTILE_OBS_KEYS

    kwargs: dict[str, object] = {
        "latent_space": latent,
        "obs": obs,
        "var": pd.DataFrame(index=HOSTILE_GENE_NAMES),
        "gene_expression": np.arange(
            n_cells * len(HOSTILE_GENE_NAMES), dtype=np.float32
        ).reshape(n_cells, len(HOSTILE_GENE_NAMES)),
        "connectivities": connectivities,
        "X_umap_2d": coordinates,
        "vector_fields": {
            "Velocity_umap_2d": coordinates * 0.25,
            "velocity_umap_2d": coordinates * 0.5,
        },
        "vector_field_default": "Velocity_umap",
        "out_dir": out_dir,
        "dataset_id": "payload-index",
        "dataset_name": "Payload index",
        "created_at": "2026-01-01T00:00:00Z",
        "obs_categorical_dtype": "uint16",
        "centroid_min_points": 1,
    }
    kwargs.update(overrides)
    prepare(**kwargs)  # type: ignore[arg-type]


@pytest.fixture(scope="module")
def export(tmp_path_factory) -> Path:
    out_dir = tmp_path_factory.mktemp("payload_index") / "exports"
    _hostile_export(out_dir)
    return out_dir


def test_identifiers_that_can_never_be_filenames_export_without_complaint(
    export: Path,
) -> None:
    obs_manifest = json.loads((export / "obs_manifest.json").read_text(encoding="utf-8"))
    var_manifest = json.loads((export / "var_manifest.json").read_text(encoding="utf-8"))

    exported_obs_keys = [field[1] for field in obs_manifest["_continuousFields"]]
    exported_obs_keys += [field[1] for field in obs_manifest["_categoricalFields"]]
    assert sorted(exported_obs_keys) == sorted(HOSTILE_OBS_KEYS)
    assert [field[1] for field in var_manifest["fields"]] == HOSTILE_GENE_NAMES


def test_no_payload_filename_carries_any_identifier(export: Path) -> None:
    for directory, pattern in (("obs", INDEXED_PAYLOAD), ("var", INDEXED_PAYLOAD)):
        names = sorted(path.name for path in (export / directory).iterdir())
        assert names, directory
        for name in names:
            assert pattern.fullmatch(name), f"{directory}/{name}"

    for name in sorted(path.name for path in (export / "vectors").iterdir()):
        assert INDEXED_VECTOR.fullmatch(name), name

    assert sorted(path.name for path in (export / "connectivity").iterdir()) == [
        "edges.dst.bin",
        "edges.src.bin",
        "edges.weights.f64.bin",
    ]

    # Nothing anywhere in the tree spells any identifier, including the
    # sanitized forms an escaping scheme would have produced.
    every_relative_path = "\n".join(
        path.relative_to(export).as_posix() for path in export.rglob("*")
    )
    for identifier in (*HOSTILE_OBS_KEYS, *HOSTILE_GENE_NAMES, "Velocity_umap", "velocity_umap"):
        assert identifier not in every_relative_path
        assert identifier.replace(" ", "_").replace("/", "_") not in every_relative_path


def test_obs_indices_are_one_dense_space_shared_by_both_field_arrays(
    export: Path,
) -> None:
    manifest = json.loads((export / "obs_manifest.json").read_text(encoding="utf-8"))
    continuous = manifest["_continuousFields"]
    categorical = manifest["_categoricalFields"]

    # The space follows the obs column order, so it interleaves the two arrays.
    assert [field[:2] for field in continuous] == [[0, "% mito"], [3, "Field"], [4, "field"]]
    assert [field[:2] for field in categorical] == [[1, "cell type"], [2, "CON"]]

    indices = [field[0] for field in continuous] + [field[0] for field in categorical]
    assert sorted(indices) == list(range(len(indices)))


def test_var_and_vector_indices_are_dense_and_position_based(export: Path) -> None:
    var_manifest = json.loads((export / "var_manifest.json").read_text(encoding="utf-8"))
    assert [field[0] for field in var_manifest["fields"]] == list(
        range(len(HOSTILE_GENE_NAMES))
    )

    identity = json.loads((export / "dataset_identity.json").read_text(encoding="utf-8"))
    assert {
        field_id: field["files"]
        for field_id, field in identity["vector_fields"]["fields"].items()
    } == {
        "Velocity_umap": {"2d": "vectors/0_2d.bin"},
        "velocity_umap": {"2d": "vectors/1_2d.bin"},
    }


def test_manifest_path_patterns_substitute_the_index(export: Path) -> None:
    obs_manifest = json.loads((export / "obs_manifest.json").read_text(encoding="utf-8"))
    var_manifest = json.loads((export / "var_manifest.json").read_text(encoding="utf-8"))

    assert obs_manifest["_obsSchemas"]["continuous"]["pathPattern"] == "obs/{index}.values.f32"
    assert obs_manifest["_obsSchemas"]["categorical"] == {
        "codesPathPattern": "obs/{index}.codes.{ext}",
        "outlierPathPattern": "obs/{index}.outliers.f32",
        "outlierExt": "f32",
        "outlierDtype": "float32",
        "outlierQuantized": False,
    }
    assert var_manifest["_varSchema"]["pathPattern"] == "var/{index}.values.f32"


def test_every_declared_artifact_exists_and_every_file_is_declared(export: Path) -> None:
    """The prepared-server reader and the written tree agree exactly."""
    declared = _declared_dataset_artifacts(export, include_published_state=False)
    on_disk = {
        path.relative_to(export).as_posix() for path in export.rglob("*") if path.is_file()
    }

    assert declared == on_disk


def test_a_shared_payload_index_is_rejected_rather_than_silently_overwritten() -> None:
    """The failure this rule exists for: one file holding two fields' values."""
    with pytest.raises(RuntimeError, match=r"exactly 0\.\.1, each used once"):
        _require_dense_payload_indices([0, 0], axis="Observation")
    with pytest.raises(RuntimeError, match=r"exactly 0\.\.2, each used once"):
        _require_dense_payload_indices([0, 1, 3], axis="Gene")
    with pytest.raises(RuntimeError, match="must be a native integer"):
        _require_dense_payload_indices(["0"], axis="Gene")

    _require_dense_payload_indices([], axis="Vector field")
    _require_dense_payload_indices([2, 0, 1], axis="Observation")


def test_the_obs_manifest_validator_rejects_a_collided_index() -> None:
    manifest = {
        "_format": "compact_v1",
        "n_points": 3,
        "centroid_outlier_quantile": 0.95,
        "latent_key": "latent_space",
        "compression": None,
        "_obsSchemas": {},
        # One continuous and one categorical field claiming the same obs/ file.
        "_continuousFields": [[0, "score"]],
        "_categoricalFields": [[0, "group", ["a"], "uint8", 255, {}]],
    }

    with pytest.raises(RuntimeError, match="Observation payload indices"):
        _identity_obs_fields_from_compact_manifest(manifest)


def test_the_var_manifest_validator_rejects_a_collided_index() -> None:
    with pytest.raises(RuntimeError, match="Gene payload indices"):
        _gene_names_from_compact_manifest(
            {"fields": [[0, "TNMD"], [0, "CIITA"]]},
        )
    with pytest.raises(ValueError, match=r"\[index, name\]"):
        _gene_names_from_compact_manifest({"fields": [["TNMD"]]})


def test_a_manifest_that_does_not_describe_its_own_files_is_rejected(
    tmp_path: Path,
) -> None:
    directory = tmp_path / "var"
    directory.mkdir()
    (directory / "0.values.f32").write_bytes(b"")

    with pytest.raises(RuntimeError, match=r"Declared but absent: \['var/1.values.f32'\]"):
        _require_declared_payloads_on_disk(
            tmp_path,
            directory_name="var",
            declared={"var/0.values.f32", "var/1.values.f32"},
            axis="Gene",
        )

    with pytest.raises(RuntimeError, match=r"Written but undeclared: \['var/0.values.f32'\]"):
        _require_declared_payloads_on_disk(
            tmp_path,
            directory_name="var",
            declared=set(),
            axis="Gene",
        )

    # An axis directory declares files and nothing else, so a directory inside
    # one is refused by kind rather than compared against a declared name.
    (directory / "nested").mkdir()
    with pytest.raises(
        RuntimeError,
        match=r"Gene payload directory holds a non-file entry: .*/var/nested$",
    ):
        _require_declared_payloads_on_disk(
            tmp_path,
            directory_name="var",
            declared={"var/0.values.f32"},
            axis="Gene",
        )


def test_the_export_root_must_hold_exactly_what_the_export_declares(
    tmp_path: Path,
) -> None:
    """The root carries the point payloads, and they are declared by path there."""
    root = tmp_path / "export"
    (root / "obs").mkdir(parents=True)
    for name in ("dataset_identity.json", "obs_manifest.json", "points_2d.bin.gz"):
        (root / name).write_bytes(b"")

    def reconcile(directories: set[str] | None = None) -> None:
        _require_declared_payloads_on_disk(
            root,
            directory_name=None,
            declared={"dataset_identity.json", "obs_manifest.json", "points_2d.bin.gz"},
            axis="Export",
            declared_directories={"obs"} if directories is None else directories,
        )

    reconcile()

    # The mutation that reached disk once: the coordinates were written
    # uncompressed while dataset_identity.json declared the compressed name.
    (root / "points_2d.bin.gz").rename(root / "points_2d.bin")
    with pytest.raises(
        RuntimeError,
        match=(
            r"Declared but absent: \['points_2d\.bin\.gz'\]\. "
            r"Written but undeclared: \['points_2d\.bin'\]\."
        ),
    ):
        reconcile()
    (root / "points_2d.bin").rename(root / "points_2d.bin.gz")
    reconcile()

    # A payload directory the export never created is absent, not invisible.
    with pytest.raises(
        RuntimeError,
        match=r"Declared but absent: \['var'\]\. Written but undeclared: \[\]\.",
    ):
        reconcile(directories={"obs", "var"})

    # A directory the export does not declare is refused by kind, so it can
    # never be mistaken for a declared payload of the same name.
    (root / "vectors").mkdir()
    with pytest.raises(
        RuntimeError,
        match=r"Export payload directory holds a non-file entry: .*/vectors$",
    ):
        reconcile()
    (root / "vectors").rmdir()

    # A declared directory that reached disk as a regular file is reported on
    # both sides, because it is neither the directory declared nor a file
    # anything declared.
    (root / "var").write_bytes(b"")
    with pytest.raises(
        RuntimeError,
        match=r"Declared but absent: \['var'\]\. Written but undeclared: \['var'\]\.",
    ):
        reconcile(directories={"obs", "var"})
    (root / "var").unlink()

    (root / "scratch.tmp").write_bytes(b"")
    with pytest.raises(
        RuntimeError,
        match=r"Declared but absent: \[\]\. Written but undeclared: \['scratch\.tmp'\]\.",
    ):
        reconcile()


def _root_reconciliation_export(out_dir: Path) -> None:
    coordinates = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]], dtype=np.float32)
    prepare(
        latent_space=coordinates,
        obs=pd.DataFrame({"group": pd.Categorical(["a", "b", "a"])}),
        X_umap_2d=coordinates,
        out_dir=out_dir,
        dataset_id="root-reconciliation",
        dataset_name="Root reconciliation",
        created_at="2026-01-01T00:00:00Z",
        obs_categorical_dtype="uint16",
        centroid_min_points=1,
        compression=6,
    )


def test_a_generation_that_declares_the_points_it_wrote_publishes(tmp_path: Path) -> None:
    out_dir = tmp_path / "points"
    _root_reconciliation_export(out_dir)

    identity = json.loads((out_dir / "dataset_identity.json").read_text(encoding="utf-8"))
    assert identity["embeddings"]["files"] == {"2d": "points_2d.bin.gz"}
    assert (out_dir / "points_2d.bin.gz").is_file()


@pytest.mark.parametrize(
    ("mutation", "expected"),
    [
        (
            "uncompressed",
            r"Declared but absent: \['points_2d\.bin\.gz'\]\. "
            r"Written but undeclared: \['points_2d\.bin'\]\.",
        ),
        (
            "unwritten",
            r"Declared but absent: \['points_2d\.bin\.gz'\]\. Written but undeclared: \[\]\.",
        ),
        (
            "stray",
            r"Declared but absent: \[\]\. Written but undeclared: \['scratch\.tmp'\]\.",
        ),
    ],
)
def test_a_published_generation_cannot_disagree_with_its_own_export_root(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
    expected: str,
) -> None:
    """The failure this rule exists for: an export the browser cannot draw.

    The declared name and the written name of a point payload are two
    expressions of one compression setting, and until the root was reconciled
    a generation whose points disagreed published successfully and then failed
    at read time with no coordinates.
    """
    out_dir = tmp_path / "points"

    def patched_write_binary(
        path: Path,
        data: np.ndarray,
        compression: int | None = None,
    ) -> Path:
        if path.name.startswith("points_"):
            if mutation == "uncompressed":
                return _write_binary(path, data, None)
            if mutation == "unwritten":
                return Path(f"{path}.gz")
        written = _write_binary(path, data, compression)
        if mutation == "stray" and path.name.startswith("points_"):
            (path.parent / "scratch.tmp").write_bytes(b"")
        return written

    monkeypatch.setattr(_generation, "_write_binary", patched_write_binary)

    with pytest.raises(RuntimeError, match=expected):
        _root_reconciliation_export(out_dir)
    assert not out_dir.exists()


@pytest.mark.parametrize("compression", [None, 6])
def test_a_full_export_root_holds_exactly_the_artifacts_it_declares(
    tmp_path: Path,
    compression: int | None,
) -> None:
    """Genes, connectivity, vector fields, and two dimension levels, both codecs."""
    out_dir = tmp_path / f"full-{compression}"
    _hostile_export(
        out_dir,
        X_umap_1d=np.linspace(0.0, 1.0, 6, dtype=np.float32).reshape(6, 1),
        compression=compression,
    )

    identity = json.loads((out_dir / "dataset_identity.json").read_text(encoding="utf-8"))
    suffix = ".gz" if compression is not None else ""
    assert identity["embeddings"]["files"] == {
        "1d": f"points_1d.bin{suffix}",
        "2d": f"points_2d.bin{suffix}",
    }
    assert sorted(entry.name for entry in out_dir.iterdir()) == [
        "connectivity",
        "connectivity_manifest.json",
        "dataset_identity.json",
        "obs",
        "obs_manifest.json",
        f"points_1d.bin{suffix}",
        f"points_2d.bin{suffix}",
        "var",
        "var_manifest.json",
        "vectors",
    ]


def test_quantized_entries_keep_the_index_and_gain_their_bounds(tmp_path: Path) -> None:
    out_dir = tmp_path / "quantized"
    _hostile_export(
        out_dir,
        var_quantization=8,
        obs_continuous_quantization=8,
        compression=6,
    )

    var_manifest = json.loads((out_dir / "var_manifest.json").read_text(encoding="utf-8"))
    assert var_manifest["_varSchema"]["pathPattern"] == "var/{index}.values.u8.gz"
    for position, field in enumerate(var_manifest["fields"]):
        assert len(field) == 4
        assert field[0] == position
        assert field[1] == HOSTILE_GENE_NAMES[position]
        assert field[2] < field[3]

    obs_manifest = json.loads((out_dir / "obs_manifest.json").read_text(encoding="utf-8"))
    assert obs_manifest["_obsSchemas"]["continuous"]["pathPattern"] == "obs/{index}.values.u8.gz"
    assert (
        obs_manifest["_obsSchemas"]["categorical"]["outlierPathPattern"]
        == "obs/{index}.outliers.u8.gz"
    )
    for field in obs_manifest["_continuousFields"]:
        assert len(field) == 4
    for field in obs_manifest["_categoricalFields"]:
        assert len(field) == 8

    assert _declared_dataset_artifacts(out_dir, include_published_state=False) == {
        path.relative_to(out_dir).as_posix() for path in out_dir.rglob("*") if path.is_file()
    }
