from __future__ import annotations

import struct
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
WEB_GUIDE = REPOSITORY_ROOT / "docs" / "user_guide" / "web_app"
SCREENSHOTS = REPOSITORY_ROOT / "docs" / "_static" / "screenshots"


def _read(relative_path: str) -> str:
    return (WEB_GUIDE / relative_path).read_text(encoding="utf-8")


def _png_dimensions(path: Path) -> tuple[int, int]:
    header = path.read_bytes()[:24]
    assert header[:8] == b"\x89PNG\r\n\x1a\n"
    assert header[12:16] == b"IHDR"
    return struct.unpack(">II", header[16:24])


def _jpeg_dimensions(path: Path) -> tuple[int, int]:
    payload = path.read_bytes()
    assert payload[:2] == b"\xff\xd8"
    offset = 2
    start_of_frame_markers = {
        0xC0,
        0xC1,
        0xC2,
        0xC3,
        0xC5,
        0xC6,
        0xC7,
        0xC9,
        0xCA,
        0xCB,
        0xCD,
        0xCE,
        0xCF,
    }
    while offset < len(payload):
        assert payload[offset] == 0xFF
        while payload[offset] == 0xFF:
            offset += 1
        marker = payload[offset]
        offset += 1
        if marker in {0xD8, 0xD9}:
            continue
        segment_length = int.from_bytes(payload[offset : offset + 2], "big")
        assert segment_length >= 2
        if marker in start_of_frame_markers:
            height = int.from_bytes(payload[offset + 3 : offset + 5], "big")
            width = int.from_bytes(payload[offset + 5 : offset + 7], "big")
            return width, height
        offset += segment_length
    raise AssertionError("JPEG does not contain a start-of-frame marker")


def test_home_internal_button_uses_renderable_content() -> None:
    home = (REPOSITORY_ROOT / "docs" / "index.md").read_text(encoding="utf-8")
    tour_button = home.split(
        "```{button-ref} user_guide/web_app/a_orientation/03_quick_tour_60_seconds",
        1,
    )[1].split("```", 1)[0]

    assert "🚀 Take the 60-second tour" in tour_button
    assert "{octicon}" not in tour_button

    hero = home.split("```{figure} ", 1)[1].split("```", 1)[0]
    assert hero.startswith("_static/screenshots/web_app/suo-cell-type-close-up.jpg")
    assert ":width: 70%" in hero
    assert ":align: center" in hero
    assert "Suo atlas" in hero
    assert "Pancreas" not in hero
    assert _jpeg_dimensions(
        SCREENSHOTS / "web_app" / "suo-cell-type-close-up.jpg"
    ) == (1920, 1201)


def test_smoke_and_multiview_docs_publish_current_exact_limits() -> None:
    smoke = _read("c_core_interactions/03_render_modes_points_vs_volumetric_smoke.md")
    multiview = _read("c_core_interactions/04_view_layout_live_snapshots_small_multiples.md")

    for exact in (
        "`Grid density:`",
        "80 / `128³`",
        "`Ray quality:`",
        "75 / `High`",
        "`Render resolution:`",
        "15 / `0.51x`",
        "`Noise detail:`",
        "58 / `128³`",
        "approximately 300 ms debounce",
        "Volumetric smoke requires a single view. Clear snapshots first.",
    ):
        assert exact in smoke
    assert "currently visible" in smoke
    assert "A `.cellucid-session` saves `Render mode:`" in smoke

    for exact in (
        "at most **three** kept snapshots",
        "`Keep view`",
        "`Locked Cam`",
        "`Unlocked Cam`",
        "`Layout: Grid compare`",
        "`Layout: Edit selected view`",
        "`Orb`, `Pan`, or `Fly`",
        "starts multiview with linked cameras",
    ):
        assert exact in multiview
    assert "Highlight pages are a separate feature" in " ".join(multiview.split())
    assert "A matching-dataset `.cellucid-session` saves and restores" in multiview


def test_velocity_docs_use_real_pancreas_captures_and_current_controls() -> None:
    index = _read("i_vector_field_velocity/index.md")
    enable = _read("i_vector_field_velocity/02_enabling_overlay_and_selecting_field.md")
    core = _read("i_vector_field_velocity/03_core_parameters_document_exact_ui_labels.md")
    captures = _read("i_vector_field_velocity/08_screenshots.md")

    assert "real\nPancreas velocity field" in index
    assert "synthetic" not in (index + captures).lower()
    assert "off by default" in enable
    assert "does not silently choose the first field" in enable
    assert "Failed to load vector field." in enable
    assert "Loading or replacing a dataset always" in enable
    assert "For the same loaded dataset, a `.cellucid-session` captures" in enable

    for exact in (
        "1K–500K",
        "0.05×–5.00×",
        "0.1s–15.0s",
        "0.5–30",
        "0%–100%",
        "`Viridis`, `Plasma`, `Turbo`, `Cividis`, `Magma`, `Coolwarm`",
        "`Sync with LOD`",
    ):
        assert exact in core

    for filename in (
        "overlay-controls.png",
        "advanced-settings.png",
        "pancreas-velocity-3d.png",
    ):
        image = SCREENSHOTS / "vector_field_velocity" / filename
        assert _png_dimensions(image) == (1440, 1000)
        directive = f"../../../_static/screenshots/vector_field_velocity/{filename}"
        assert directive in captures
    assert captures.count(":width: 1440px") == 3
    assert "3,696 cells" in captures
    assert "2D **Planar**" in captures
    assert "3D `velocity_umap`" in captures


def test_every_reused_velocity_capture_keeps_real_pancreas_provenance() -> None:
    relative_paths = (
        "docs/user_guide/python_package/c_data_preparation_api/"
        "08_vector_fields_velocity_displacement.md",
        "docs/user_guide/python_package/d_viewing_apis/"
        "08_anndata_mode_show_anndata_and_serve_anndata.md",
        "docs/user_guide/python_package/f_notebooks_tutorials/"
        "33_vector_fields_and_velocity_overlay_end_to_end.md",
        "docs/user_guide/web_app/b_data_loading/10_standard_pancreas_dataset.md",
    )
    for relative_path in relative_paths:
        source = (REPOSITORY_ROOT / relative_path).read_text(encoding="utf-8")
        assert "synthetic 2D dataset" not in source
        assert "Build 2026-07-27.1" in source


def test_community_annotation_docs_distinguish_local_and_portable_sessions() -> None:
    index = _read("j_community_annotation/index.md")
    ui = _read("j_community_annotation/03_ui_reference.md")
    author = _read("j_community_annotation/02_author_guide.md")

    assert "Community Annotation** accordion starts collapsed" in index
    for exact in (
        "**Connect repo**",
        "**Continue with GitHub**",
        "**Add repo**",
        "**Reload**",
        "**Connect**",
        "**Pull latest**",
        "**Publish**",
    ):
        assert exact in index + ui
    assert "dataset id + `owner/repo` + branch + numeric GitHub user id" in index
    assert "not** part of an ordinary `.cellucid-session`" in index
    assert "intentionally excluded from" in ui
    assert "cellucid-consensus.json" in ui
    assert "consensus_<datasetId>.json" not in (index + ui + author)


def test_community_annotation_docs_publish_current_network_and_repository_contract() -> None:
    index = _read("j_community_annotation/index.md")
    annotator = _read("j_community_annotation/01_annotator_guide.md")
    author = _read("j_community_annotation/02_author_guide.md")
    ui = _read("j_community_annotation/03_ui_reference.md")
    combined = "\n".join((index, annotator, author, ui))
    normalized = " ".join(combined.split())
    author_normalized = " ".join(author.split())

    for exact in (
        "annotations/config.schema.json",
        "annotations/moderation/merges.schema.json",
        "rsync -a --exclude '.git/'",
        "fails the complete Pull without changing annotation state",
        "before selecting the raw-cache scope or downloading any user or moderation files",
        "exactly the same keys as `fieldsToAnnotate`",
        "There is no UI route override",
        "Branch protection alone does not select fork + Pull Request",
        "equal net vote and equal upvote count",
        "`<fieldKey>:<categoryLabel>` when the field key contains no colon",
        "`fk~${encodeURIComponent(fieldKey)}:<categoryLabel>`",
        "`%3A` in the prefixed component",
        "from urllib.parse import unquote",
        "`/auth/login`",
        "`/auth/callback`",
        "`/auth/user`",
        "`/auth/installations`",
        "`/auth/installation-repos`",
        "`/cap/lookup-cells`",
        "`/cap/search-datasets`",
        "`/api/repos/*`",
        "`ALLOWED_ORIGINS`, `GITHUB_CLIENT_ID`, and `GITHUB_CLIENT_SECRET`",
        '`{"status":"ok","service":"Cellucid GitHub Auth","contractVersion":1,',
        "1,000,000 decoded UTF-8 bytes",
        "10,000 active user files",
        "64,000,000 aggregate decoded UTF-8 bytes across those user files",
        "persistent raw-file cache requires IndexedDB",
        "Cellucid never switches to an in-memory cache",
        "Only the GitHub identity that created a merge record can edit or delete its note",
    ):
        assert exact in author_normalized

    for guide in (annotator, author, ui):
        guide_normalized = " ".join(guide.split())
        assert "configured Cellucid Worker" in guide_normalized
        assert "`/cap/lookup-cells`" in guide_normalized
        assert "GitHub bearer token and OAuth cookies are not sent to CAP" in guide_normalized
    assert "`/cap/search-datasets`" not in annotator
    assert "`/cap/search-datasets`" not in ui
    assert "`/cap/lookup-cells` or" not in annotator
    assert "`/cap/lookup-cells` or" not in ui
    assert "related datasets" not in annotator

    for guide in (index, annotator, author):
        guide_normalized = " ".join(guide.split())
        assert "fork + Pull Request only when the repository allows forking" in guide_normalized
        assert "does not switch to the other route" in guide_normalized
    assert "without either route, they cannot publish" in normalized

    ui_normalized = " ".join(ui.split())
    assert "**Name:**" in ui
    assert "direct browser request to `https://pub.orcid.org`" in ui_normalized
    assert "at least three Unicode code points" in ui_normalized
    assert "no surrounding whitespace and at most 240 Unicode code points" in ui_normalized
    assert "`/v3.0/{id}/person`" in ui_normalized
    assert "`/v3.0/expanded-search/`" in ui_normalized
    assert "at most eight results" in ui_normalized
    assert "60-second per-modal cache" in ui_normalized
    assert "omits credentials and referrer information" in ui_normalized
    assert "empty ORCID or an exact checksum-valid ORCID iD" in ui_normalized
    assert "name text is search input, not a valid saved ORCID iD" in ui_normalized
    assert "direct vote on the merged target wins" in ui_normalized
    assert "majority of that user's votes across the merged members" in ui_normalized

    assert 'encoded_field_key = field_key_part[3:]' in author
    assert '"%3a" in encoded_field_key.lower()' in author

    for stale in (
        "Cellucid will warn and ignore it",
        "If a field has no explicit entry in `annotatableSettings`, Cellucid uses defaults",
        "ties between top net-vote suggestions are always Disputed",
        "If the top net-vote score is tied across multiple suggestions",
        "`/api/*` proxy",
        "`GITHUB_APP_ID`",
        "`GITHUB_PRIVATE_KEY`",
        "some UIs may show delegated bundle votes",
        "Depending on the UI version",
        "If you do not have push access: Cellucid uses a fork + Pull Request flow",
        "CAP search sends queries to `https://celltype.info/graphql`",
        "same-origin Worker",
        "otherwise memory-only",
        "otherwise PR flow",
        "cached in-memory",
        "downloaded cache may already contain validated remote files",
        "Validated files are cached before the dataset-field check",
        "publish to a dedicated “staging” branch",
    ):
        assert stale not in normalized

    assert "raw-file cache stores file bodies in IndexedDB and its SHA index" in normalized
    assert "Raw pulled GitHub file bodies use IndexedDB" in normalized


def test_official_sample_state_docs_publish_integrity_verified_static_loading() -> None:
    index = _read("l_sessions_sharing/index.md")
    manual = _read("l_sessions_sharing/03_save_restore_ux.md")
    presets = _read("l_sessions_sharing/04_official_sample_states.md")
    manual_normalized = " ".join(manual.split())
    normalized = " ".join(presets.split())

    assert "{doc}`04_official_sample_states`" in index
    for dataset_id in ("suo", "garcia", "he", "pancreas"):
        assert f"exports/{dataset_id}/default.cellucid-session" in presets
        assert f"exports/{dataset_id}/state-snapshots.json" in presets

    for exact in (
        '{ "states": ["default.cellucid-session"] }',
        "`state_manifest` and `state_sha256` as a pair",
        "catalog's SHA-256",
        "Only after the scientific dataset has been published",
        "automatically applies the verified static starting view",
        "`cell_type` observation field",
        "3-D **Orbit** camera",
        "light theme and grid background",
        "no `cinematic/camera` chunk",
        "**Loop playback**, or **Autoplay**",
        "moving-camera playback remains off",
        "The official catalog uses this capability",
        "A custom `exportsBaseUrl` catalog may deliberately opt in",
        "same paired fields and exact sidecars",
        "only explicit catalog advertisement enables the state requests",
        "Jupyter, GitHub, local-user, and remote sources are not probed",
        "a GitHub repository (`github-repo`)",
        "a locally selected prepared export (`local-user`)",
        "a remote prepared-data server (`remote`)",
        "does not infer sidecar paths",
        "manual, identity-matched working-state file",
        "07_versioning_compatibility_and_dataset_identity",
        # A published state is regenerated whenever the reviewed view or the
        # session format changes, so a `state_sha256` copied into prose goes
        # stale silently. The page has to send readers to the catalog instead.
        "Read the current value out of the catalog",
        "all four generation directories",
    ):
        assert exact in normalized

    assert "applies the matching static state automatically" in index
    assert (
        "integrity-verified static starting state applies automatically"
        in manual_normalized
    )
    for stale in (
        "cellucid-datasets/sessions",
        "__close-up.cellucid-session",
        "never restored automatically",
        "Custom and non-advertised sources are not probed",
        # The four published presets were regenerated; this pinned the digest
        # they had before, and there is no fifth (kanemaru) generation.
        "5581569118eaa992",
        "all five generation directories",
    ):
        assert stale not in presets


def test_session_docs_publish_only_the_atomic_current_contract() -> None:
    relative_paths = (
        "l_sessions_sharing/01_session_mental_model.md",
        "l_sessions_sharing/02_what_gets_saved_and_restored.md",
        "l_sessions_sharing/03_save_restore_ux.md",
        "l_sessions_sharing/04_official_sample_states.md",
        "l_sessions_sharing/05_share_workflows_links_bundles_exports.md",
        "l_sessions_sharing/06_collaboration_best_practices.md",
        "l_sessions_sharing/07_versioning_compatibility_and_dataset_identity.md",
        "l_sessions_sharing/08_security_privacy_and_trust.md",
        "l_sessions_sharing/09_edge_cases.md",
        "l_sessions_sharing/10_troubleshooting_sessions.md",
        "l_sessions_sharing/12_reference.md",
        "p_developer_docs/10_sessions_persistence_and_serialization.md",
    )
    sources = {path: _read(path) for path in relative_paths}
    combined = "\n".join(sources.values())

    reference = sources["l_sessions_sharing/12_reference.md"]
    developer = sources["p_developer_docs/10_sessions_persistence_and_serialization.md"]
    manual = sources["l_sessions_sharing/03_save_restore_ux.md"]
    official = sources["l_sessions_sharing/04_official_sample_states.md"]
    reference_normalized = " ".join(reference.split())
    developer_normalized = " ".join(developer.split())

    for exact in (
        "`createdAt`, `datasetFingerprint`, and `chunks`",
        # The fingerprint carries five fields, not four: the viewer records a
        # `cellOrder` digest so a dataset republished re-sorted under an
        # unchanged identity cannot silently re-point saved selections.
        "`sourceType`, `datasetId`, `cellCount`, `varCount`, and `cellOrder`",
        "`cellOrder` has exactly `dimension` and `digest`",
        "16 lowercase hexadecimal characters",
        "All five fingerprint fields must match",
        "`analysis/cache-inventory`",
        "`cinematic/camera`",
        "one for every categorical user-defined field",
        "one membership payload for every advertised highlight group",
        "exact ordered match to `analysis/cache-inventory`",
        "**Session fully restored**",
        "whole restore",
        "rolls back",
    ):
        assert exact in reference_normalized

    for exact in (
        "one current `.cellucid-session` writer/reader contract",
        "accepts exactly one current schema",
        "complete document all-or-nothing",
        "preflights the complete",
        "before constructing the native decompressor",
        "Direct serializer APIs reject with the coded `AbortError`",
        "Do not introduce a reader for an older shape",
    ):
        assert exact in developer_normalized

    assert "strict manual loader rejects" in " ".join(manual.split())
    assert (
        "Downloading it and choosing **Load State** is unsupported"
        in " ".join(official.split())
    )

    # `07_...` is where every other session page sends the reader for "the exact
    # matching rules", so it has to state the whole five-field rule and give the
    # cell-order failures their own outcome rows.
    identity = " ".join(
        sources[
            "l_sessions_sharing/07_versioning_compatibility_and_dataset_identity.md"
        ].split()
    )
    troubleshooting = " ".join(
        sources["l_sessions_sharing/10_troubleshooting_sessions.md"].split()
    )
    for exact in (
        "The session stores five fields",
        "All five must equal the currently loaded dataset",
        "`cellOrder`: a record of *which ordering*",
        "exactly 1, 2, or 3",
        "A different embedding dimension is on screen",
        "stores its cells in a different order",
        "A session file saved before Cellucid recorded cell order",
        "accepts a fingerprint **with or without** `cellOrder`",
    ):
        assert exact in identity
    for exact in (
        "All five must match",
        "`cellOrder`, itself exactly a `dimension` (1, 2, or 3)",
    ):
        assert exact in troubleshooting

    for retired in (
        "Restoring only dataset-agnostic layout",
        "dataset-dependent chunks are skipped",
        "skips dataset-dependent chunks",
        "skipped dataset-dependent",
        "Session restored (eager stage complete)",
        "`dependsOn` (optional)",
        "skips unknown chunk contributors",
        "isolates chunk restore failures",
        "accept a partial restore",
        "may trigger a background lazy restore",
        # The fingerprint stopped being four fields when `cellOrder` was added;
        # these are the exact sentences that said otherwise.
        "All four fingerprint fields must match",
        "All four must equal",
        "The session stores four fields",
        "All four must match",
    ):
        assert retired not in combined
