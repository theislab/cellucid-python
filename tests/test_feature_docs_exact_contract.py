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


def test_home_page_leads_with_the_app_and_current_captures() -> None:
    home = (REPOSITORY_ROOT / "docs" / "index.md").read_text(encoding="utf-8")

    # `button-ref` drops a nested role, so the internal button is plain text.
    tour_button = home.split(
        "```{button-ref} user_guide/web_app/a_orientation/03_quick_tour_60_seconds",
        1,
    )[1].split("```", 1)[0]
    assert "Take the 60-second tour" in tour_button
    assert "{octicon}" not in tour_button

    # Open Cellucid is the first, filled action; the tour is the second.
    assert home.index("```{button-link} https://www.cellucid.com") < home.index(
        "```{button-ref} user_guide/web_app/a_orientation/03_quick_tour_60_seconds"
    )
    assert ":color: primary" in home.split("```{button-link}", 1)[1].split("```", 1)[0]

    # The hero is a light and a dark capture of the Suo sample, each at two
    # widths for `srcset`, shown by the theme the reader chose.
    hero = home.split('<figure class="cl-hero-shot">', 1)[1].split("</figure>", 1)[0]
    for theme in ("light", "dark"):
        image = hero.split(f'<img class="only-{theme}"', 1)[1].split(">", 1)[0]
        for width in (1280, 2400):
            source = f"_static/screenshots/home/hero-{theme}-{width}.webp"
            assert f"{source} {width}w" in image
            assert (SCREENSHOTS / "home" / f"hero-{theme}-{width}.webp").is_file()
    assert "Suo" in hero
    assert "561,947" in hero
    assert not (SCREENSHOTS / "web_app" / "suo-cell-type-close-up.jpg").exists()


def test_smoke_and_multiview_docs_publish_current_exact_limits() -> None:
    smoke = _read("c_core_interactions/03_render_modes_points_vs_volumetric_smoke.md")
    multiview = _read("c_core_interactions/04_view_layout_live_snapshots_small_multiples.md")

    for exact in (
        "`Grid density`",
        "80 / `128³`",
        "`Ray quality`",
        "75 / `High`",
        "`Render resolution`",
        "33 / `0.50×`",
        "`Noise detail`",
        "80 / `128³`",
        "65 / `1.00`",
        "40 / `2.0`",
        "15 / `0.3`",
        "50 / `1.00×`",
        "`Still`",
        "approximately 300 ms debounce",
        "Volumetric smoke requires a single view. Clear snapshots first.",
    ):
        assert exact in smoke
    assert "currently visible" in smoke
    assert "A `.cellucid-session` saves `Render mode`" in smoke

    for exact in (
        "any number of kept views",
        "GPU memory, not a count, is\nthe limit",
        "`Keep view`",
        "`View 2 shows`",
        "`Link cameras`",
        "`Layout: Grid compare`",
        "`Layout: Edit selected view`",
        "`Orbit`, `Pan` or `Fly`",
        "**Link cameras** is on by default",
        "Views that show different embeddings always keep their\nown cameras",
    ):
        assert exact in multiview
    assert "at most **three**" not in multiview
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
    assert "3D `velocity` payload" in captures
    assert "velocity_umap" not in captures


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
        # The capture's own record names the web commit it shows
        # (``captures.json``); a build stamp in prose only goes stale.
        assert "Pancreas sample" in source
        assert "Build 2026-" not in source


def test_community_annotation_docs_distinguish_local_and_portable_sessions() -> None:
    index = _read("j_community_annotation/index.md")
    ui = _read("j_community_annotation/03_ui_reference.md")
    author = _read("j_community_annotation/02_author_guide.md")

    assert "Community annotation** accordion starts collapsed" in index
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
        "`/api/repos/*`",
        "`ALLOWED_ORIGINS`, `GITHUB_CLIENT_ID`, and `GITHUB_CLIENT_SECRET`",
        '`{"status":"ok","service":"Cellucid GitHub Auth","contractVersion":1,',
        "1,000,000 decoded UTF-8 bytes",
        "A Pull has no limit on how many user files there are or on their total size",
        "streaming the cache one file at a time",
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
    assert "`/cap/search-datasets`" not in author
    assert "`/cap/lookup-cells` or" not in annotator
    assert "`/cap/lookup-cells` or" not in ui
    assert "related datasets" not in annotator

    for guide in (index, annotator, author):
        guide_normalized = " ".join(guide.split())
        assert "fork + Pull Request only when the repository allows forking" in guide_normalized
        assert "does not switch to the other route" in guide_normalized
    assert "without either route, they cannot publish" in normalized

    ui_normalized = " ".join(ui.split())
    assert "**Name**" in ui
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
    index_normalized = " ".join(index.split())
    normalized = " ".join(presets.split())

    assert "{doc}`04_official_sample_states`" in index
    for exact in (
        # The catalog names a dataset's default session by its `session`
        # entry, pinned by exact length and digest.
        '"session": {',
        '"path": "pancreas.cellucid-session"',
        "stops reading past the byte length the catalog declares",
        "verifies its SHA-256 against the catalog",
        "Only after the dataset has been published",
        "Cellucid applies the verified static starting view",
        "`cell_type` observation field",
        "3-D **Orbit** camera",
        "light theme and the **Studio (light)** background",
        "no `cinematic/camera` chunk",
        "**Loop playback**, or **Autoplay**",
        "moving-camera playback remains off",
        "All three are required",
        "Nothing else is probed",
        "a file or folder chosen with **Open…** or dropped from your computer",
        "a dataset served by `cellucid serve`, `serve_anndata()` or Jupyter",
        "a `session` key there is ignored",
        "does not infer session paths",
        "manual, identity-matched working-state file",
        "07_versioning_compatibility_and_dataset_identity",
        # A published session is regenerated whenever the reviewed view or the
        # session format changes, so a digest copied into prose goes stale
        # silently. The page has to send readers to the catalog instead.
        "Read the current values out of the catalog",
        "scripts/publish_session_preset.mjs",
    ):
        assert exact in normalized

    assert "the session's bytes are verified" in index_normalized
    assert "integrity-verified and applies automatically" in manual_normalized
    for stale in (
        "state-snapshots.json",
        "state_manifest",
        "state_sha256",
        "state_bytes",
        "exports-v3",
        "exportsBaseUrl",
        "prepared",
        "5581569118eaa992",
        "20769ea20075a212",
    ):
        assert stale not in presets


#: Names that only the prepared export format had (removed in 0.9.5): its
#: sample tree, its query key, its catalog fields, its file-picker button and
#: the screenshot host that served it.
PREPARED_FORMAT_NAMES = (
    "exports-v3",
    "exportsBaseUrl",
    "state-snapshots.json",
    "state_manifest",
    "state_sha256",
    "state_bytes",
    "user-data-browse-btn",
    "serve-exports",
)


def test_no_documentation_page_names_the_removed_prepared_format() -> None:
    """Every page, published or working material, names the AnnData catalog,
    ``?sampleCatalog=`` and the current controls. The screenshot checklist's
    index keeps the one table that maps the old names to the new ones."""
    docs = REPOSITORY_ROOT / "docs"
    mapping_table = WEB_GUIDE / "r_screenshot_checklist" / "index.md"
    failures = []
    for page in sorted(docs.rglob("*.md")):
        if "_build" in page.parts or page == mapping_table:
            continue
        text = page.read_text(encoding="utf-8")
        failures.extend(
            f"{page.relative_to(docs)}: {name}" for name in PREPARED_FORMAT_NAMES if name in text
        )
    assert failures == []


def test_session_docs_publish_only_the_atomic_current_contract() -> None:
    relative_paths = (
        "l_sessions_sharing/01_session_mental_model.md",
        "l_sessions_sharing/02_what_gets_saved_and_restored.md",
        "l_sessions_sharing/03_save_restore_ux.md",
        "l_sessions_sharing/04_official_sample_states.md",
        "l_sessions_sharing/05_share_workflows_links_sessions_data_files.md",
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
        # Session format 2: the version is part of the signature.
        "`CELLUCID_SESSION/2\\n`",
        "`createdAt`, `datasetFingerprint`, `source`, `requires`, and `chunks`",
        # The fingerprint carries six fields: the viewer records a `cellOrder`
        # digest so a dataset republished re-sorted under an unchanged identity
        # cannot silently re-point saved selections, the data mapping the
        # dataset was shown with, and the display name used in messages.
        "`datasetId`, `name`, `cellCount`, `varCount`, `cellOrder`, and `dataMapping`",
        "`cellOrder` has exactly `embedding`, `dimension`, and `digest`",
        "16 lowercase hexadecimal characters",
        "naming every missing key at once",
        "`ui/figure-export`",
        "`analysis/docked`",
        "`cinematic/camera`",
        "one for every categorical user-defined field",
        "one membership payload for every advertised highlight group",
        "The category this selection was saved relative to no longer holds the same cells.",
        "**Session fully restored**",
        "whole restore",
        "rolls back",
    ):
        assert exact in reference_normalized

    for exact in (
        "one current `.cellucid-session` writer/reader contract",
        "accepts exactly one session format",
        "complete document all-or-nothing",
        "preflights the complete",
        "before constructing the native decompressor",
        "Direct serializer APIs reject with the coded `AbortError`",
        "Do not introduce a reader for an older shape",
    ):
        assert exact in developer_normalized

    assert "strict manual loader rejects" in " ".join(manual.split())
    assert (
        "Downloading it and choosing **Load** is unsupported"
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
        "The session stores six fields",
        "The three scalars and `cellOrder` must equal the currently loaded dataset",
        # Every data-mapping group, not only the embeddings, is restored when
        # the open dataset can use the saved choice (`describeDataMappingMismatch`).
        "The data mapping need not equal what is on screen",
        "`cellOrder`: a record of *which ordering*",
        "`dataMapping`: which stored arrays the dataset was shown with",
        "exactly 1, 2, or 3",
        "Restoring is *verify, then apply*",
        "The session turns on other embeddings, or names other defaults, that this dataset can show",
        "This session shows obsm['X_tsne'] columns [0, 1], which this dataset does not have.",
        "The session names another graph, modality, gene names, expression matrix, latent space or vector bindings, ones this dataset can use",
        "stores its cells in a different order",
        "A file that does not start with `CELLUCID_SESSION/2`",
        "This file is not a Cellucid session file of format 2 (it does not start with CELLUCID_SESSION/2).",
        "The dataset lacks an observation field or gene the session uses",
        # The Python reader requires the same fingerprint and never verifies
        # the digest; the page must say both.
        "a file without `cellOrder.embedding` or `dataMapping` is refused",
        "Python validates `cellOrder`'s shape and never its digest",
    ):
        assert exact in identity
    for retired in (
        "accepts a fingerprint **with or without** `cellOrder`",
        "A different embedding dimension is on screen",
        "All six must equal the currently loaded dataset",
        "A session file saved before Cellucid recorded cell order",
        "Cellucid does not store a format-version field",
        "in a session format this release no longer reads",
        "The data mapping must too, with one exception: its embeddings",
    ):
        assert retired not in identity
    for exact in (
        "All but `name` must match",
        "`cellOrder`, itself exactly the `embedding` the digest was taken over",
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
        # Format 1 leftovers.
        "`analysis/cache-inventory`",
        "analysis/artifacts/bulk-gene",
        "five-chunk",
        "All five fingerprint fields must match",
    ):
        assert retired not in combined
