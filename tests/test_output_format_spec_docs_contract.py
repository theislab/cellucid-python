"""Contract tests for the on-disk export format specification page.

The specification page is the written contract shared by the Python writer, the
R writer, the web app reader, and both dataset repositories. Nothing else in
this repository reads it, which is how it drifted away from what the exporter
actually writes. These tests pin the page to a real export produced by
:func:`cellucid.prepare` and to the exporter's own constants, so a format change
that is not documented fails here.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cellucid.connectivity_contract import connectivity_index_storage
from cellucid.prepare_data import prepare

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SPEC_PAGE = (
    REPOSITORY_ROOT
    / "docs"
    / "user_guide"
    / "python_package"
    / "c_data_preparation_api"
    / "09_output_format_specification_exports_directory.md"
)

INPUT_GENES = ["gene_a", "gene_b", "gene_c", "gene_d"]
EXPORTED_GENES = ["gene_c", "gene_a"]

OBS_SECTION = ("## Obs manifest:", "## Var manifest:")
VAR_SECTION = ("## Var manifest:", "## Connectivity manifest:")
CONNECTIVITY_SECTION = ("## Connectivity manifest:", "## Vector files:")
IDENTITY_SECTION = ("## Required file: `dataset_identity.json`", "## Points files:")
FAST_PATH_SECTION = ("### Minimal export folder (required)", "### Typical full export folder")


def _spec() -> str:
    return SPEC_PAGE.read_text(encoding="utf-8")


def _section(text: str, bounds: tuple[str, str]) -> str:
    start = text.index(bounds[0])
    return text[start : text.index(bounds[1], start)]


def _documented_top_level_keys(section: str) -> set[str]:
    """Return the keys of the ``Top-level keys`` bullet list of one section."""
    start = section.index("Top-level keys")
    end = section.index("\n#### ", start)
    return set(re.findall(r'^- `"([A-Za-z_]+)"', section[start:end], re.MULTILINE))


def _documented_json_object(section: str, marker: str) -> dict:
    for block in re.findall(r"```json\n(.*?)\n```", section, re.DOTALL):
        if marker in block:
            return json.loads(block)
    raise AssertionError(f"no JSON example containing {marker!r} in the section")


@pytest.fixture(scope="module")
def export(tmp_path_factory) -> Path:
    """One real export exercising every manifest the specification describes."""
    out_dir = tmp_path_factory.mktemp("spec_contract_export") / "exports"
    n_cells = 6
    coordinates = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 0.5], [0.25, 0.75]],
        dtype=np.float32,
    )
    connectivities = np.zeros((n_cells, n_cells), dtype=np.float64)
    for source, destination in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5)):
        connectivities[source, destination] = 0.5
        connectivities[destination, source] = 0.5

    prepare(
        latent_space=coordinates,
        obs=pd.DataFrame(
            {
                "pct_mito": np.linspace(0.0, 1.0, n_cells, dtype=np.float32),
                "leiden": pd.Categorical(["a", "b", "a", "b", "c", "c"]),
            }
        ),
        var=pd.DataFrame(index=INPUT_GENES),
        gene_expression=np.arange(n_cells * len(INPUT_GENES), dtype=np.float32).reshape(
            n_cells, len(INPUT_GENES)
        ),
        gene_identifiers=EXPORTED_GENES,
        connectivities=connectivities,
        X_umap_2d=coordinates,
        out_dir=out_dir,
        dataset_id="spec-contract",
        dataset_name="Spec contract",
        created_at="2026-01-01T00:00:00Z",
        obs_categorical_dtype="uint16",
    )
    return out_dir


def _manifest(export: Path, name: str) -> dict:
    return json.loads((export / name).read_text(encoding="utf-8"))


def test_obs_manifest_documents_every_emitted_top_level_key(export: Path) -> None:
    """CEL-AUDIT-0029: the page listed 6 of the 8 required obs manifest keys."""
    documented = _documented_top_level_keys(_section(_spec(), OBS_SECTION))
    emitted = set(_manifest(export, "obs_manifest.json"))

    assert documented == emitted
    assert {"centroid_outlier_quantile", "latent_key"} <= documented


def test_var_manifest_documents_every_emitted_top_level_key(export: Path) -> None:
    documented = _documented_top_level_keys(_section(_spec(), VAR_SECTION))

    assert documented == set(_manifest(export, "var_manifest.json"))


def test_connectivity_manifest_example_matches_the_emitted_key_set(export: Path) -> None:
    documented = _documented_json_object(
        _section(_spec(), CONNECTIVITY_SECTION),
        '"edge_pairs"',
    )

    assert set(documented) == set(_manifest(export, "connectivity_manifest.json"))


def test_identity_example_covers_every_emitted_key(export: Path) -> None:
    documented = _documented_json_object(
        _section(_spec(), IDENTITY_SECTION),
        '"cellucid_data_version"',
    )
    identity = _manifest(export, "dataset_identity.json")

    assert set(identity) <= set(documented)
    assert set(documented) - set(identity) <= {"source", "vector_fields"}
    assert set(documented["stats"]) == set(identity["stats"])
    assert set(documented["embeddings"]) == set(identity["embeddings"])
    assert set(documented["export_settings"]) == set(identity["export_settings"])


def test_identity_obs_fields_example_lists_continuous_before_categorical(
    export: Path,
) -> None:
    """The reader compares identity obs_fields to the obs manifest positionally."""
    documented = _documented_json_object(
        _section(_spec(), IDENTITY_SECTION),
        '"cellucid_data_version"',
    )
    documented_kinds = [field["kind"] for field in documented["obs_fields"]]
    emitted_kinds = [
        field["kind"] for field in _manifest(export, "dataset_identity.json")["obs_fields"]
    ]

    assert documented_kinds == sorted(documented_kinds, key="continuous".__ne__)
    assert emitted_kinds == sorted(emitted_kinds, key="continuous".__ne__)
    assert set(documented_kinds) == {"continuous", "category"}


def test_identity_n_genes_is_documented_as_the_exported_subset(export: Path) -> None:
    """CEL-AUDIT-0031: the page claimed n_genes was the input matrix gene count."""
    identity = _manifest(export, "dataset_identity.json")
    var_manifest = _manifest(export, "var_manifest.json")

    assert identity["stats"]["n_genes"] == len(EXPORTED_GENES)
    assert identity["stats"]["n_genes"] == len(var_manifest["fields"])
    assert identity["stats"]["n_genes"] < len(INPUT_GENES)

    spec = _spec()
    assert "reflects the input matrix gene count" not in spec
    assert "`n_genes` is the **exported** gene count" in spec


def test_minimal_export_folder_lists_every_always_required_entry() -> None:
    """CEL-AUDIT-0030: the page claimed identity plus one points file sufficed."""
    fast_path = _section(_spec(), FAST_PATH_SECTION)

    for required in ("dataset_identity.json", "obs_manifest.json", "obs/", "points_2d.bin"):
        assert required in fast_path
    assert "**not optional and never degrades gracefully**" in fast_path


def test_connectivity_index_width_documents_both_exporter_widths() -> None:
    """CEL-AUDIT-0032: the page documented uint16 only; shipped datasets use uint32."""
    _, small_dtype, small_bytes = connectivity_index_storage(65_536)
    _, large_dtype, large_bytes = connectivity_index_storage(65_537)
    assert (small_dtype, small_bytes) == ("uint16", 2)
    assert (large_dtype, large_bytes) == ("uint32", 4)

    section = _section(_spec(), CONNECTIVITY_SECTION)
    assert "| `1` … `65536` | `\"uint16\"` | `2` |" in section
    assert "| `65537` … `4294967296` | `\"uint32\"` | `4` |" in section
    assert "not fixed at `uint16`" in section


def test_connectivity_edge_rules_match_the_written_bytes(export: Path) -> None:
    manifest = _manifest(export, "connectivity_manifest.json")
    index_dtype = np.dtype(manifest["index_dtype"]).newbyteorder("<")
    sources = np.fromfile(export / manifest["sourcesPath"], dtype=index_dtype)
    destinations = np.fromfile(export / manifest["destinationsPath"], dtype=index_dtype)
    weights = np.fromfile(export / manifest["weightsPath"], dtype="<f8")

    assert sources.size == destinations.size == weights.size == manifest["n_edges"]
    assert np.all(sources < destinations)
    assert np.all(weights > 0)
    order = np.lexsort((destinations, sources))
    assert np.array_equal(order, np.arange(sources.size))
    assert len({(int(a), int(b)) for a, b in zip(sources, destinations, strict=True)}) == (
        sources.size
    )

    section = _section(_spec(), CONNECTIVITY_SECTION)
    assert "**`src[i] < dst[i]`**" in section
    assert "**finite and strictly positive**" in section
    assert "**strict lexicographic order** by `(src, dst)`" in section
    assert "weights are finite, non-negative" not in section


def test_exact_key_set_rule_is_stated_for_every_manifest() -> None:
    spec = _spec()

    assert "## Exact key sets (applies to every manifest)" in spec
    for manifest_name in (
        "`dataset_identity.json`",
        "`obs_manifest.json`",
        "`var_manifest.json`",
        "`connectivity_manifest.json`",
    ):
        assert manifest_name in _section(spec, ("## Exact key sets", "---"))
    assert "an **unknown extra** key is equally fatal" in spec
    assert "all eight are required, and no other key is allowed" in spec
    assert "all seven are required, and no other key is allowed" in spec
    assert "all twelve keys below and no others" in spec
