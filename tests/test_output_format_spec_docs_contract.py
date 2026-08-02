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

from cellucid._compression import deterministic_gzip_compress
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
VECTOR_SECTION = ("## Vector files:", "## Compression strategy")
DETERMINISM_SECTION = (
    "## Determinism and reproducibility",
    "## Multi-dataset export roots",
)
KEY_SET_SECTION = ("## Exact key sets (applies to every manifest)", "---")


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


@pytest.fixture(scope="module")
def vector_export(tmp_path_factory) -> Path:
    """One compressed export carrying vector fields and a name-only ``source``.

    The uncompressed export above cannot exercise the ``.gz`` producer paths,
    the gzip header the page pins, or the optional-``source`` rule, so those
    claims are checked against this second export instead.
    """
    out_dir = tmp_path_factory.mktemp("spec_contract_vectors") / "exports"
    n_cells = 6
    coordinates_2d = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 0.5], [0.25, 0.75]],
        dtype=np.float32,
    )
    coordinates_3d = np.hstack(
        [coordinates_2d, np.linspace(0.0, 1.0, n_cells, dtype=np.float32).reshape(-1, 1)]
    )

    prepare(
        latent_space=coordinates_2d,
        obs=pd.DataFrame(
            {
                "pct_mito": np.linspace(0.0, 1.0, n_cells, dtype=np.float32),
                "leiden": pd.Categorical(["a", "b", "a", "b", "c", "c"]),
            }
        ),
        X_umap_2d=coordinates_2d,
        X_umap_3d=coordinates_3d,
        vector_fields={
            "velocity_umap_2d": coordinates_2d * 0.1,
            "velocity_umap_3d": coordinates_3d * 0.1,
        },
        vector_field_default="velocity_umap",
        compression=6,
        out_dir=out_dir,
        dataset_id="spec-contract-vectors",
        dataset_name="Spec contract vectors",
        created_at="2026-01-01T00:00:00Z",
        obs_categorical_dtype="uint16",
        source_name="Spec contract source",
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


def test_identity_example_stats_counts_are_self_consistent(export: Path) -> None:
    """CEL-0012 sibling: the example declared 12 obs fields beside a 2-element array.

    The reader rejects both disagreements, so an example that shows them is an
    example nobody can load.
    """
    documented = _documented_json_object(
        _section(_spec(), IDENTITY_SECTION),
        '"cellucid_data_version"',
    )
    identity = _manifest(export, "dataset_identity.json")

    for payload in (documented, identity):
        stats = payload["stats"]
        assert stats["n_obs_fields"] == len(payload["obs_fields"])
        assert stats["n_obs_fields"] == (
            stats["n_categorical_fields"] + stats["n_continuous_fields"]
        )

    identity_section = _section(_spec(), IDENTITY_SECTION)
    assert "`stats.n_obs_fields` must equal `stats.n_categorical_fields +" in identity_section


def test_identity_source_documents_name_as_its_only_required_key(
    vector_export: Path,
) -> None:
    """CEL-0011: the page called every key of ``source`` required; two are optional."""
    identity = _manifest(vector_export, "dataset_identity.json")

    assert identity["source"] == {"name": "Spec contract source"}
    assert "filename" not in identity["source"]

    key_sets = _section(_spec(), KEY_SET_SECTION)
    assert "`name` is required while `url`" in key_sets
    assert "Only `dataset_identity.json` has optional keys at all" not in _spec()

    identity_section = _section(_spec(), IDENTITY_SECTION)
    assert "only `name` is required" in identity_section
    assert "`source.filename`" in identity_section
    assert "**not part of the export format**" in identity_section


def test_vector_field_label_and_basis_are_documented_as_derived(
    vector_export: Path,
) -> None:
    """CEL-0010: the example's ``"label": "Velocity (UMAP)"`` cannot be loaded."""
    identity = _manifest(vector_export, "dataset_identity.json")
    documented = _documented_json_object(
        _section(_spec(), IDENTITY_SECTION),
        '"cellucid_data_version"',
    )

    for payload in (identity, documented):
        fields = payload["vector_fields"]["fields"]
        compression = payload["export_settings"]["compression"]
        suffix = "" if compression is None else ".gz"
        for index, (field_id, field) in enumerate(fields.items()):
            assert field["label"] == field_id
            assert field["basis"] == "umap"
            assert field["default_dimension"] == max(field["available_dimensions"])
            assert field["files"] == {
                f"{dimension}d": f"vectors/{index}_{dimension}d.bin{suffix}"
                for dimension in field["available_dimensions"]
            }
        assert payload["vector_fields"]["default_field"] in fields

    for relative in identity["vector_fields"]["fields"]["velocity_umap"]["files"].values():
        assert (vector_export / relative).is_file()

    identity_section = _section(_spec(), IDENTITY_SECTION)
    documented_example = re.search(
        r"```json\n(.*?\"cellucid_data_version\".*?)\n```",
        identity_section,
        re.DOTALL,
    )
    assert documented_example is not None
    assert '"label": "Velocity (UMAP)"' not in documented_example.group(1)
    assert "**derived by the writer, not chosen by you**" in identity_section
    # The page must name the failure, not merely avoid demonstrating it.
    assert '`"label": "Velocity (UMAP)"` on a field called' in identity_section


def test_vector_files_section_documents_the_recomputed_producer_path() -> None:
    """The reader rebuilds every vector path; the page claimed it never has to."""
    section = _section(_spec(), VECTOR_SECTION)

    assert "so a reader never has to build one" not in _spec()
    assert "vectors/{index}_{dim}d.bin.gz   otherwise" in section
    assert "it rejects the dataset when the declared\nstring is anything else" in section
    assert "must be the **largest** entry of `available_dimensions`" in section


def test_categorical_code_width_caps_are_documented(tmp_path) -> None:
    """The reserved missing marker caps a field's category count; the page omitted it."""
    section = _section(_spec(), OBS_SECTION)
    assert "at most 255 categories" in section
    assert "at most 65,535 under `\"uint16\"`" in section

    n_cells = 256
    coordinates = np.stack(
        [
            np.linspace(0.0, 1.0, n_cells, dtype=np.float32),
            np.linspace(1.0, 0.0, n_cells, dtype=np.float32),
        ],
        axis=1,
    )
    with pytest.raises(ValueError, match="uint8 supports at most 255"):
        prepare(
            latent_space=coordinates,
            obs=pd.DataFrame({"crowded": [f"c{index}" for index in range(n_cells)]}),
            X_umap_2d=coordinates,
            centroid_min_points=1,
            out_dir=tmp_path / "exports",
            dataset_id="spec-contract-cap",
            dataset_name="Spec contract cap",
            created_at="2026-01-01T00:00:00Z",
            obs_categorical_dtype="uint8",
        )


def test_connectivity_working_set_ceiling_is_documented() -> None:
    """The reader refuses a conformant export above a 512 MiB edge working set."""
    section = _section(_spec(), CONNECTIVITY_SECTION)

    assert "536,870,912 bytes (512 MiB)" in section
    assert "| `\"uint32\"` | `44 × n_edges + 4 × n_cells` |" in section
    assert "| `\"uint16\"` | `48 × n_edges + 4 × n_cells` |" in section
    # The stated edge ceilings must follow from the stated budget and rates.
    limit = 512 * 1024 * 1024
    assert limit // 44 == 12_201_611
    assert (limit - 4 * 65_536) // 48 == 11_179_349
    assert "12.2 million edges" in section
    assert "11.1 million" in section
    # The largest shipped graph must still clear the ceiling the page quotes.
    assert "6,279,148 edges over 561,947 cells" in section
    assert limit > 44 * 6_279_148 + 4 * 561_947


def test_gzip_header_is_canonical_at_every_documented_level(
    vector_export: Path,
) -> None:
    """CEL-0012(a)/CEL-0009: the whole ten-byte header, not only ``OS``, is pinned.

    The page names every byte of the header and quotes the RFC 1952 2.3.1 ``XFL``
    rule the R writer was made to match, so each named byte is checked against a
    real export and against the encoder the exporter shares with the server.
    """
    header = (vector_export / "points_2d.bin.gz").read_bytes()[:10]

    assert header[:2] == b"\x1f\x8b"  # gzip magic
    assert header[2] == 0x08  # deflate
    assert header[3] == 0x00  # no optional header fields, FNAME included
    assert header[4:8] == b"\x00\x00\x00\x00"  # MTIME 0
    assert header[8] == 0x00  # XFL at the level 6 the fixture exports with
    assert header[9] == 0xFF  # OS: "unknown", not the build platform

    for level, extra_flags in ((1, 0x04), (6, 0x00), (9, 0x02)):
        member = deterministic_gzip_compress(b"spec contract", compresslevel=level)
        assert member[:8] == b"\x1f\x8b\x08\x00\x00\x00\x00\x00"
        assert member[8] == extra_flags
        assert member[9] == 0xFF

    section = _section(_spec(), DETERMINISM_SECTION)
    assert "ten-byte gzip member header is identical between the two writers, on" in section
    assert "every platform, unconditionally**" in section
    assert "`MTIME = 0`" in section
    assert "(`0x04`\nat level 1, `0x02` at level 9, `0x00` in between)" in section
    assert "`OS = 0xff`, the\n\"unknown\" code" in section


def test_cross_writer_differences_match_the_current_writers() -> None:
    """CEL-0012: the page promised byte-identical payloads including their ``.gz``.

    CEL-0009 then made the R header canonical, which retired the first of the six
    differences the retraction listed. This pins both halves: the five surviving
    differences with their numbering, and the absence of the obsolete ``OS``-byte
    claim, so neither the stale nor the narrowed version can drift back.
    """
    spec = _spec()
    section = _section(spec, DETERMINISM_SECTION)

    assert "They are also deterministic across the two writers." not in spec
    assert "Two encoder differences remain" not in spec
    assert (
        "**semantic equality of every export, plus byte\nidentity of the uncompressed "
        "binary payloads**"
    ) in section

    # The gzip header is no longer a cross-writer difference anywhere on the page.
    for retracted in (
        "Byte identity does **not** extend to the compressed `.gz` files",
        "gzip member header `OS` byte",
        "`gzfile()` leaves the header to zlib",
        "`0x03` on Unix, `0x0b` on Windows",
        "Difference 1 alone is enough",
    ):
        assert retracted not in spec

    # What replaced it: a guarantee narrowed to the deflate stream, with its
    # measurement, so the page cannot decay back into "the .gz files differ".
    assert "byte-identical across the writers whenever the\ntwo zlib builds agree" in section
    assert "verified for a 4 KiB float32 payload at all nine" in section
    assert "R linked against zlib 1.3.2 and" in section
    assert "Python against zlib 1.2.12" in section
    assert "4,495,576-byte payload compressed at level 6 to" in section
    assert "4,160,196 bytes under R and 4,163,855 bytes under Python" in section
    assert "example of build skew rather than a property of the packages" in section
    assert "chunked deflate equals one-shot deflate" in section
    assert "`deflateInit2` parameters differed" in section
    assert "compare the payloads either uncompressed or after decompression" in section

    # The table enumerates exactly the five surviving differences, in order.
    rows = re.findall(r"^\| (\d+) \| (.+?) \|", section, re.MULTILINE)
    assert rows == [
        ("1", "end of every JSON file"),
        ("2", "`dataset_identity.json` layout"),
        ("3", "separators in the three compact manifests"),
        ("4", "non-ASCII text in any JSON file"),
        ("5", "categorical centroid `position`"),
    ]
    for known_difference in (
        "one trailing `\\n` (`writeLines()`)",
        "`jsonlite::prettify(indent = 2)`",
        "`, ` and `: `",
        "escaped as `\\uXXXX`",
        "sum a category's coordinates in a different order",
    ):
        assert known_difference in section
    assert "Difference 5 is the only numeric one" in section


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
