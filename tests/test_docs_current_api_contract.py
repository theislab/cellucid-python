from __future__ import annotations

import ast
import binascii
import inspect
import json
import re
import struct
import textwrap
from collections.abc import Iterator
from pathlib import Path

import h5py

import cellucid
from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.session_bundle import CellucidSessionBundle

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
README = REPOSITORY_ROOT / "README.md"
DOCS_ROOT = REPOSITORY_ROOT / "docs"
EXAMPLES_ROOT = REPOSITORY_ROOT / "examples"
TUTORIAL_H5AD = (
    DOCS_ROOT
    / "user_guide"
    / "python_package"
    / "f_notebooks_tutorials"
    / "datasets"
    / "endocrinogenesis_day15.5_raw.h5ad"
)
PYTHON_FENCE_PATTERN = re.compile(
    r"(?ms)^[ \t]*```(?P<info>[^\n]*)\n(?P<source>.*?)^[ \t]*```[ \t]*$"
)
PYTHON_FENCE_INFO = {
    "python",
    "py",
    "{code-block} python",
    "{code-block} py",
    "{code-cell} python",
    "{code-cell} py",
}
PUBLIC_CALL_NAMES = {name for name in cellucid.__all__ if callable(getattr(cellucid, name))} | {
    "from_file",
    "apply_to_anndata",
}
COMMENTED_CALL_START = re.compile(
    r"^\s*"
    r"(?:(?:await|return)\s+)?"
    r"(?:[A-Za-z_]\w*\s*=\s*)?"
    r"(?:[A-Za-z_]\w*\.)*"
    rf"(?:{'|'.join(sorted(map(re.escape, PUBLIC_CALL_NAMES), key=len, reverse=True))})"
    r"\s*\("
)
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def _validated_png_dimensions(path: Path) -> tuple[int, int]:
    payload = path.read_bytes()
    assert payload.startswith(PNG_SIGNATURE), path

    offset = len(PNG_SIGNATURE)
    dimensions: tuple[int, int] | None = None
    saw_image_data = False
    saw_end = False
    while offset < len(payload):
        assert offset + 12 <= len(payload), path
        chunk_length = struct.unpack(">I", payload[offset : offset + 4])[0]
        chunk_type = payload[offset + 4 : offset + 8]
        chunk_data_start = offset + 8
        chunk_data_end = chunk_data_start + chunk_length
        chunk_crc_end = chunk_data_end + 4
        assert chunk_crc_end <= len(payload), (path, chunk_type)

        chunk_data = payload[chunk_data_start:chunk_data_end]
        expected_crc = struct.unpack(">I", payload[chunk_data_end:chunk_crc_end])[0]
        actual_crc = binascii.crc32(chunk_type + chunk_data) & 0xFFFFFFFF
        assert actual_crc == expected_crc, (path, chunk_type)

        if chunk_type == b"IHDR":
            assert dimensions is None and offset == len(PNG_SIGNATURE), path
            assert chunk_length == 13, path
            width, height = struct.unpack(">II", chunk_data[:8])
            assert width > 0 and height > 0, path
            assert chunk_data[10:12] == b"\x00\x00", path
            assert chunk_data[12] in {0, 1}, path
            dimensions = (width, height)
        elif chunk_type == b"IDAT":
            assert dimensions is not None, path
            saw_image_data = True
        elif chunk_type == b"IEND":
            assert chunk_length == 0 and not saw_end, path
            saw_end = True
            offset = chunk_crc_end
            break

        offset = chunk_crc_end

    assert dimensions is not None and saw_image_data and saw_end, path
    assert offset == len(payload), path
    return dimensions


def _markdown_files() -> list[Path]:
    return [README, *sorted(DOCS_ROOT.rglob("*.md"))]


def _notebook_files() -> list[Path]:
    files = [
        *DOCS_ROOT.rglob("*.ipynb"),
        *EXAMPLES_ROOT.rglob("*.ipynb"),
    ]
    return sorted(path for path in files if "_build" not in path.parts)


def _markdown_python_sources() -> Iterator[tuple[str, str]]:
    for path in _markdown_files():
        text = path.read_text(encoding="utf-8")
        for match in PYTHON_FENCE_PATTERN.finditer(text):
            if match.group("info").strip().lower() not in PYTHON_FENCE_INFO:
                continue
            source = textwrap.dedent(match.group("source"))
            line = text[: match.start("source")].count("\n") + 1
            yield f"{path.relative_to(REPOSITORY_ROOT)}:{line}", source


def _notebook_python_sources() -> Iterator[tuple[str, str]]:
    for path in _notebook_files():
        notebook = json.loads(path.read_text(encoding="utf-8"))
        for cell_number, cell in enumerate(notebook.get("cells", []), start=1):
            if cell.get("cell_type") != "code":
                continue
            source = "".join(cell.get("source", []))
            # IPython line magics are valid notebook syntax but not Python AST.
            source = "\n".join(
                line for line in source.splitlines() if not line.lstrip().startswith("%")
            )
            if not source.strip():
                continue
            yield f"{path.relative_to(REPOSITORY_ROOT)}:cell {cell_number}", source


def _python_sources() -> Iterator[tuple[str, str]]:
    yield from _markdown_python_sources()
    yield from _notebook_python_sources()


def _commented_api_sources() -> Iterator[tuple[str, str]]:
    """Yield exact public-API calls intentionally shown as commented code."""
    for location, source in _python_sources():
        lines = source.splitlines()
        line_index = 0
        while line_index < len(lines):
            stripped = lines[line_index].lstrip()
            if not stripped.startswith("#"):
                line_index += 1
                continue

            first = stripped[1:].lstrip()
            if not COMMENTED_CALL_START.match(first):
                line_index += 1
                continue

            candidate_lines = [first]
            try:
                ast.parse(first)
            except SyntaxError:
                # A documented multiline call must begin with an open call.
                # This condition excludes prose such as
                # ``# prepare(...) writes an export``.
                if not first.rstrip().endswith("("):
                    line_index += 1
                    continue
                next_index = line_index + 1
                while next_index < len(lines):
                    continuation = lines[next_index].lstrip()
                    if not continuation.startswith("#"):
                        break
                    candidate_lines.append(continuation[1:].lstrip())
                    next_index += 1
                    try:
                        ast.parse("\n".join(candidate_lines))
                    except SyntaxError:
                        continue
                    break
                line_index = next_index
            else:
                line_index += 1

            yield (
                f"{location} (commented line {line_index + 1})",
                textwrap.dedent("\n".join(candidate_lines)),
            )


def _all_python_sources() -> Iterator[tuple[str, str]]:
    yield from _python_sources()
    yield from _commented_api_sources()


def _python_examples() -> Iterator[tuple[str, ast.Module]]:
    for location, source in _all_python_sources():
        try:
            tree = ast.parse(source)
        except SyntaxError as error:
            raise AssertionError(
                f"{location}: syntactically incomplete Python example: {error.msg}"
            ) from error
        yield location, tree


def _call_name(call: ast.Call) -> str | None:
    if isinstance(call.func, ast.Name):
        return call.func.id
    if isinstance(call.func, ast.Attribute):
        return call.func.attr
    return None


def test_python_examples_are_syntactically_complete() -> None:
    failures: list[str] = []
    for location, source in _all_python_sources():
        try:
            ast.parse(source)
        except SyntaxError as error:
            failures.append(
                f"{location}:{error.lineno}: {error.msg}\n"
                f"{error.text.rstrip() if error.text else source.rstrip()}"
            )

    assert not failures, "\n\n".join(failures)


def test_runnable_python_examples_match_current_public_signatures() -> None:
    public_objects = {
        name: getattr(cellucid, name)
        for name in cellucid.__all__
        if callable(getattr(cellucid, name))
    }
    public_objects["from_file"] = AnnDataAdapter.from_file
    public_objects["apply_to_anndata"] = CellucidSessionBundle.apply_to_anndata
    signatures = {
        name: inspect.signature(public_object) for name, public_object in public_objects.items()
    }
    failures: list[str] = []

    for location, tree in _python_examples():
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            name = _call_name(node)
            if name not in signatures:
                continue
            signature = signatures[name]
            keywords = {keyword.arg for keyword in node.keywords if keyword.arg is not None}
            has_expansion = any(keyword.arg is None for keyword in node.keywords)
            has_starred_args = any(isinstance(argument, ast.Starred) for argument in node.args)
            unknown = sorted(keywords - set(signature.parameters))
            if unknown:
                failures.append(f"{location}:{node.lineno}: {name} has removed keywords {unknown}")
            if has_expansion or has_starred_args:
                failures.append(
                    f"{location}:{node.lineno}: {name} uses dynamic argument expansion, "
                    "so its exact current call contract cannot be verified"
                )
                continue

            parameters = list(signature.parameters.values())
            if (
                isinstance(node.func, ast.Attribute)
                and parameters
                and parameters[0].name in {"self", "cls"}
            ):
                parameters = parameters[1:]
            positional_parameters = [
                parameter
                for parameter in parameters
                if parameter.kind
                in {
                    inspect.Parameter.POSITIONAL_ONLY,
                    inspect.Parameter.POSITIONAL_OR_KEYWORD,
                }
            ]
            has_varargs = any(
                parameter.kind is inspect.Parameter.VAR_POSITIONAL for parameter in parameters
            )
            if len(node.args) > len(positional_parameters) and not has_varargs:
                failures.append(
                    f"{location}:{node.lineno}: {name} has too many positional arguments"
                )
            consumed_positionals = {
                parameter.name for parameter in positional_parameters[: len(node.args)]
            }
            missing_positionals = sorted(
                parameter.name
                for parameter in positional_parameters
                if parameter.default is inspect.Parameter.empty
                and parameter.name not in consumed_positionals
                and parameter.name not in keywords
            )
            if missing_positionals:
                failures.append(
                    f"{location}:{node.lineno}: {name} misses required positional "
                    f"arguments {missing_positionals}"
                )

            required_keywords = {
                parameter.name
                for parameter in parameters
                if parameter.kind is inspect.Parameter.KEYWORD_ONLY
                and parameter.default is inspect.Parameter.empty
            }
            missing = sorted(required_keywords - keywords)
            if missing:
                failures.append(
                    f"{location}:{node.lineno}: {name} misses required keywords {missing}"
                )

    assert not failures, "\n".join(failures)


def test_cli_examples_use_only_the_exact_current_contract() -> None:
    failures: list[str] = []
    for path in [*_markdown_files(), *_notebook_files()]:
        text = path.read_text(encoding="utf-8")
        relative = path.relative_to(REPOSITORY_ROOT)
        if "--no-backed" in text or "--backed" in text:
            failures.append(f"{relative}: documents a removed backed-mode CLI option")

        logical_text = re.sub(r"\\\s*\n\s*", " ", text)
        for line_number, line in enumerate(logical_text.splitlines(), start=1):
            if "cellucid serve" not in line or not re.search(r"\.(?:h5ad|zarr)\b", line):
                continue
            if "--dataset-name" not in line or "--dataset-id" not in line:
                failures.append(
                    f"{relative}:{line_number}: direct AnnData CLI example lacks exact identity"
                )

    assert not failures, "\n".join(failures)


def test_docs_do_not_reintroduce_removed_or_invented_runtime_contracts() -> None:
    text_by_path = {
        path.relative_to(REPOSITORY_ROOT): path.read_text(encoding="utf-8")
        for path in _markdown_files()
    }
    all_docs_text = "\n".join(text_by_path.values())
    python_guide_text = "\n".join(
        text
        for path, text in text_by_path.items()
        if path == Path("README.md") or path.parts[:3] == ("docs", "user_guide", "python_package")
    )

    stale_claims = [
        "prefetching (one-time per web build)",
        "Viewer UI cache ready",
        "later entries overwrite earlier ones",
        "jupyter_server_proxy",
        "automatic proxy selection",
        "first available port",
        "ports 8765, 8766",
        "warn-and-skip or partial-apply mismatch policy",
        "mismatch policies",
        "not all implemented",
        "add a screenshot",
        "depends on UI version",
        "show_anndata(data, **kwargs)",
        "**adapter_kwargs",
        "Cellucid already catches exceptions inside hook callbacks",
        "Python should ignore unknown keys",
        "the server may still start",
        "the session was for a different dataset and was skipped",
        "apply_cellucid_session_to_anndata requires pandas",
        "For large `.h5ad`, prefer `.zarr`",
        "Clear UI cache and retry",
        "Retry `viewer.get_session_bundle()`",
        "skip unknown chunk contributors",
        "Pandas is a core `cellucid` dependency and is used to write the planned",
    ]
    for claim in stale_claims:
        assert claim not in python_guide_text

    retired_web_claims = [
        "bridge.html",
        "bridge.js",
        "fetchWithExportsBridge",
        "Browser **.zarr folder picker**",
        "Browser **.zarr** folder picker",
        "H5AD/Zarr/Prepared",
        "user selects folder/h5ad/zarr",
        "Load a `.zarr` Directory Directly",
        "Zarr directory selected, but it errors immediately",
        "Switch to Chrome/Edge",
        "Chrome is the most reliable for folder access",
        "`velocity_umap` if present, otherwise the first field",
        "synthetic rendering benchmark** (when available)",
    ]
    for claim in retired_web_claims:
        assert claim not in all_docs_text

    incomplete_identity_calls = [
        "show_anndata(adata)",
        'show_anndata("data.h5ad")',
        'show_anndata("data.zarr")',
        "show_anndata(<data.h5ad>)",
        "show_anndata(<data.zarr or AnnData>)",
        "serve_anndata(adata)",
        "serve_anndata(<data.h5ad>)",
        "serve_anndata(<data.zarr>)",
    ]
    for call in incomplete_identity_calls:
        assert call not in python_guide_text

    assert not re.search(
        r"serve_anndata\s*\([^)]*\bbacked\s*=\s*True",
        python_guide_text,
        re.DOTALL,
    )
    for nonexistent_command in ("cellucid prepare", "cellucid cache", "cellucid session"):
        assert nonexistent_command not in python_guide_text

    screenshot_placeholder = DOCS_ROOT / "_static" / "screenshots" / "placeholder-screenshot.svg"
    assert not screenshot_placeholder.exists()

    vector_page = text_by_path[
        Path(
            "docs/user_guide/python_package/d_viewing_apis/"
            "08_anndata_mode_show_anndata_and_serve_anndata.md"
        )
    ]
    assert "vector_field_default" in vector_page
    assert "rejects duplicate gene IDs" in vector_page
    cli_page = text_by_path[
        Path("docs/user_guide/python_package/d_viewing_apis/04_cli_cellucid_serve_quickstart.md")
    ]
    assert "--vector-field-default" in cli_page

    environment_page = text_by_path[
        Path(
            "docs/user_guide/python_package/h_developer_docs/"
            "06_configuration_env_vars_and_logging.md"
        )
    ]
    assert "no environment-variable configuration layer" in environment_page
    assert "`CELLUCID_CLIENT_SERVER_URL` are ignored" in environment_page

    installation_page = text_by_path[
        Path("docs/user_guide/python_package/a_landing_pages/02_installation.md")
    ]
    assert r".\.venv\Scripts\Activate.ps1" in installation_page
    assert r".\.venv\Scripts\activate.bat" in installation_page

    server_source = (REPOSITORY_ROOT / "src/cellucid/server.py").read_text(encoding="utf-8")
    anndata_server_source = (REPOSITORY_ROOT / "src/cellucid/anndata_server.py").read_text(
        encoding="utf-8"
    )
    for source in (server_source, anndata_server_source):
        assert "one-time per web build" not in source
        assert '"Viewer UI cache"' not in source
        assert '"Viewer UI generation"' in source
        assert '"establishing exact configured source"' in source


def test_documented_vector_declarations_require_dimension_suffixes() -> None:
    all_text = "\n".join(
        path.read_text(encoding="utf-8") for path in [*_markdown_files(), *_notebook_files()]
    )
    stale_contract_phrases = [
        "Implicit keys are also supported",
        "**Fallback (implicit key",
        "- Implicit: `<field>_<basis>`",
        "Style B: implicit dimension",
        "explicit_dim_suffix=True",
        "explicit_dim_suffix=False",
    ]
    for phrase in stale_contract_phrases:
        assert phrase not in all_text

    failures: list[str] = []
    for location, tree in _python_examples():
        for node in ast.walk(tree):
            dictionaries: list[ast.Dict] = []
            if isinstance(node, ast.Assign) and isinstance(node.value, ast.Dict):
                if any(
                    isinstance(target, ast.Name) and target.id == "vector_fields"
                    for target in node.targets
                ):
                    dictionaries.append(node.value)
            elif isinstance(node, ast.Call):
                for keyword in node.keywords:
                    if keyword.arg in {"obsm", "vector_fields"} and isinstance(
                        keyword.value, ast.Dict
                    ):
                        dictionaries.append(keyword.value)

            for dictionary in dictionaries:
                for key in dictionary.keys:
                    if (
                        isinstance(key, ast.Constant)
                        and isinstance(key.value, str)
                        and key.value != "X_umap"
                        and key.value.endswith("_umap")
                    ):
                        failures.append(
                            f"{location}:{key.lineno}: unsuffixed vector declaration {key.value!r}"
                        )

            if not isinstance(node, ast.Assign):
                continue
            for target in node.targets:
                if not isinstance(target, ast.Subscript):
                    continue
                if not (
                    isinstance(target.value, ast.Attribute)
                    and target.value.attr == "obsm"
                    and isinstance(target.slice, ast.Constant)
                    and isinstance(target.slice.value, str)
                ):
                    continue
                key = target.slice.value
                if key != "X_umap" and key.endswith("_umap"):
                    failures.append(
                        f"{location}:{target.lineno}: unsuffixed obsm vector declaration {key!r}"
                    )

    assert not failures, "\n".join(failures)


def test_documented_multi_field_vector_defaults_use_exact_runtime_field_ids() -> None:
    relative_paths = (
        "docs/user_guide/python_package/d_viewing_apis/"
        "08_anndata_mode_show_anndata_and_serve_anndata.md",
        "docs/user_guide/python_package/f_notebooks_tutorials/"
        "33_vector_fields_and_velocity_overlay_end_to_end.md",
        "docs/user_guide/web_app/b_data_loading/05_jupyter_tutorial.md",
    )
    page_text = {
        relative: (REPOSITORY_ROOT / relative).read_text(encoding="utf-8")
        for relative in relative_paths
    }
    combined = "\n".join(page_text.values())

    assert 'vector_field_default="velocity_umap"' in page_text[relative_paths[0]]
    assert 'vector_field_default="velocity_umap"' in page_text[relative_paths[1]]
    assert 'default_vector_field = "velocity_umap"' in page_text[relative_paths[1]]
    assert "`velocity_umap` or\n  `T_fwd_umap`" in page_text[relative_paths[1]]
    assert 'vector_field_default="velocity_umap"' in page_text[relative_paths[2]]
    assert 'vector_field_default="T_fwd_umap"' in page_text[relative_paths[2]]
    assert not re.search(
        r'\bvector_field_default\s*=\s*["\'](?:velocity|T_fwd)["\']',
        combined,
    )
    assert 'default_vector_field = "velocity"' not in combined


def test_documented_embedding_requirements_use_exact_current_keys() -> None:
    requirement_pages = [
        DOCS_ROOT
        / "user_guide"
        / "web_app"
        / "b_data_loading"
        / "03_browser_file_picker_tutorial.md",
        DOCS_ROOT
        / "user_guide"
        / "python_package"
        / "c_data_preparation_api"
        / "01_prepare_export_overview.md",
        DOCS_ROOT
        / "user_guide"
        / "python_package"
        / "c_data_preparation_api"
        / "02_input_requirements_global.md",
        DOCS_ROOT / "user_guide" / "python_package" / "c_data_preparation_api" / "index.md",
    ]
    exact_keys = ("X_umap_1d", "X_umap_2d", "X_umap_3d")

    failures: list[str] = []
    for path in requirement_pages:
        text = path.read_text(encoding="utf-8")
        relative = path.relative_to(REPOSITORY_ROOT)
        for key in exact_keys:
            if key not in text:
                failures.append(f"{relative}: omits current embedding key {key!r}")
        for unsupported_access in ("obsm['X_umap']", 'obsm["X_umap"]'):
            if unsupported_access in text:
                failures.append(
                    f"{relative}: presents unsupported embedding access {unsupported_access!r}"
                )

    overview = requirement_pages[1].read_text(encoding="utf-8")
    for required_minimum_artifact in (
        "points_1d.bin(.gz)",
        "obs_manifest.json",
        "dataset_identity.json",
    ):
        assert required_minimum_artifact in overview

    assert not failures, "\n".join(failures)


def test_every_documented_connectivity_payload_is_the_weighted_current_contract() -> None:
    connectivity_pages = {
        path.relative_to(REPOSITORY_ROOT): path.read_text(encoding="utf-8")
        for path in _markdown_files()
        if "edges.src.bin" in path.read_text(encoding="utf-8")
    }
    assert connectivity_pages
    failures: list[str] = []
    for path, text in connectivity_pages.items():
        if "edges.weights.f64.bin" not in text and "Float64 weight" not in text:
            failures.append(f"{path}: source indices are documented without weights")
        for stale_claim in (
            "unweighted undirected edge list",
            "writes two parallel arrays",
            "writes two binary arrays",
            "exactly symmetric, binary graph",
            "matrix must already be square, finite, binary",
            "Weighted graphs are rejected",
        ):
            if stale_claim in text:
                failures.append(f"{path}: contains stale claim {stale_claim!r}")
    assert not failures, "\n".join(failures)

    format_pages = "\n".join(connectivity_pages.values())
    for required_manifest_field in (
        "weightsPath",
        "weight_dtype",
        "weight_bytes",
    ):
        assert required_manifest_field in format_pages


def test_server_and_jupyter_docs_state_exact_current_security_and_loading() -> None:
    text_by_path = {
        path.relative_to(REPOSITORY_ROOT): path.read_text(encoding="utf-8")
        for path in _markdown_files()
    }
    all_text = "\n".join(text_by_path.values())

    stale_security_claims = [
        "There is no per-event cryptographic auth token",
        "does not currently enforce `viewerToken`",
        "not authenticated beyond `viewerId`",
        "viewerId-routed but not authenticated",
        "Authentication-ish",
    ]
    for claim in stale_security_claims:
        assert claim not in all_text

    stale_zarr_claims = [
        "lazy-ish (chunked access)",
        "inherently chunked/lazy",
        "scripted server startup with chunked access",
        "array access is chunked/lazy",
        "lazy loading from `.h5ad` backed mode or `.zarr`",
        "lazy loading** for backed `.h5ad` and `.zarr`",
    ]
    for claim in stale_zarr_claims:
        assert claim not in all_text

    prepared_server_page = text_by_path[
        Path(
            "docs/user_guide/python_package/d_viewing_apis/"
            "05_python_serve_and_serve_anndata_quickstart.md"
        )
    ]
    advanced_server_page = text_by_path[
        Path("docs/user_guide/python_package/d_viewing_apis/09_server_mode_advanced.md")
    ]
    for page in (prepared_server_page, advanced_server_page):
        assert "/?source=remote" in page
        assert "first paint" in page

    failures: list[str] = []
    bare_loopback_viewer = re.compile(
        r"http://127\.0\.0\.1:(?:<port>|\d+)/"
        r"(?=(?:\s|[`)\],.;:]|$))"
    )
    launch_text_by_path = dict(text_by_path)
    launch_text_by_path.update(
        {
            path.relative_to(REPOSITORY_ROOT): path.read_text(encoding="utf-8")
            for path in _notebook_files()
        }
    )
    for path, text in launch_text_by_path.items():
        for line_number, line in enumerate(text.splitlines(), start=1):
            if bare_loopback_viewer.search(line):
                failures.append(
                    f"{path}:{line_number}: prepared/direct server launch "
                    "must use its exact printed Viewer URL"
                )

    cli_source = (REPOSITORY_ROOT / "src" / "cellucid" / "cli.py").read_text(encoding="utf-8")
    assert not bare_loopback_viewer.search(cli_source)
    assert "open the exact Viewer URL printed by Cellucid" in cli_source

    session_route = re.compile(r"/_cellucid/session_bundle\?[^\s`]*")
    for path, text in text_by_path.items():
        for line_number, line in enumerate(text.splitlines(), start=1):
            for match in session_route.finditer(line):
                if "viewerToken=" not in match.group(0):
                    failures.append(f"{path}:{line_number}: session upload route omits viewerToken")
    assert not failures, "\n".join(failures)


def test_bundled_tutorial_h5ad_uses_the_exact_current_root_encoding() -> None:
    with h5py.File(TUTORIAL_H5AD, "r") as handle:
        assert handle.attrs["encoding-type"] == "anndata"
        assert handle.attrs["encoding-version"] == "0.1.0"
        assert "X_umap_2d" in handle["obsm"]
        assert "X_umap" not in handle["obsm"]


def test_docs_preserve_closed_jupyter_session_and_zarr_contracts() -> None:
    markdown_text = "\n".join(path.read_text(encoding="utf-8") for path in _markdown_files())
    pyproject_text = (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")

    assert "anndata>=0.12.19,<0.13" in pyproject_text
    assert "zarr>=3.1.4,<4" in pyproject_text
    assert "numcodecs>=0.16.3,<0.17" in pyproject_text
    assert "pip install zarr" not in markdown_text.lower()
    assert "frontend_console" not in markdown_text

    removed_custom_tutorial = (
        DOCS_ROOT
        / "user_guide"
        / "python_package"
        / "f_notebooks_tutorials"
        / "31_custom_message_schemas_and_frontend_extensions.md"
    )
    exact_schema_tutorial = removed_custom_tutorial.with_name(
        "31_exact_message_schemas_and_diagnostics.md"
    )
    assert not removed_custom_tutorial.exists()
    assert exact_schema_tutorial.is_file()

    event_page = (
        DOCS_ROOT
        / "user_guide"
        / "python_package"
        / "e_jupyter_hooks"
        / "07_frontend_to_python_events.md"
    ).read_text(encoding="utf-8")
    for event_type in (
        "selection",
        "hover",
        "click",
        "ready",
        "pong",
        "debug_snapshot",
        "session_bundle",
    ):
        assert f"### `{event_type}`" in event_page
    assert "Unknown event types" in event_page
    assert "rejected before any hook runs" in event_page

    session_reference = (
        DOCS_ROOT
        / "user_guide"
        / "python_package"
        / "g_api_reference_coverage"
        / "api"
        / "sessions.md"
    ).read_text(encoding="utf-8")
    assert (
        "An unknown chunk, unknown contributor, or mismatched contributor/chunk pair"
        in session_reference
    )
    assert "inplace=False` rejects a backed target" in session_reference


def test_embedding_docs_require_a_defined_normalization_range() -> None:
    page = (
        DOCS_ROOT
        / "user_guide"
        / "python_package"
        / "c_data_preparation_api"
        / "03_embeddings_and_coordinates.md"
    ).read_text(encoding="utf-8")
    assert "all points identical is rejected" in page
    assert "no finite normalization scale exists" in page
    assert "nonzero finite range on at least one axis" in page
    assert "collapsed point cloud" not in page
    assert "export is “valid”" not in page


def test_current_startup_screenshots_and_documented_defaults_are_exact() -> None:
    screenshot_directory = DOCS_ROOT / "_static" / "screenshots" / "web_app"
    for filename in ("welcome-startup.png", "startup-loaded-build.png"):
        assert _validated_png_dimensions(screenshot_directory / filename) == (
            1440,
            1000,
        )

    quick_tour = (
        DOCS_ROOT / "user_guide" / "web_app" / "a_orientation" / "03_quick_tour_60_seconds.md"
    ).read_text(encoding="utf-8")
    glossary = (
        DOCS_ROOT / "user_guide" / "web_app" / "a_orientation" / "04_ui_glossary_terminology.md"
    ).read_text(encoding="utf-8")
    camera_page = (
        DOCS_ROOT / "user_guide" / "web_app" / "c_core_interactions" / "07_screenshots.md"
    ).read_text(encoding="utf-8")
    annotation_index = (
        DOCS_ROOT / "user_guide" / "web_app" / "j_community_annotation" / "index.md"
    ).read_text(encoding="utf-8")
    combined = "\n".join((quick_tour, glossary))

    assert "../../../_static/screenshots/web_app/welcome-startup.png" in quick_tour
    assert "../../../_static/screenshots/web_app/startup-loaded-build.png" in combined
    assert "Every ordinary bundled or catalog startup opens the welcome overlay" in quick_tour
    assert "explicit Jupyter, remote-server, or GitHub-served startup skips" in quick_tour
    assert "Explicit Jupyter, remote-server, and GitHub-served startup" in glossary
    assert "Suo dataset remains the sample\ncatalog default" in quick_tour
    assert "1D or 2D starts in **Planar**" in quick_tour
    assert "3D starts in **Orbit**" in quick_tour
    assert "default light grid" in quick_tour
    assert "Build 2026-07-27.1" in combined
    assert "exact value in bug reports" in glossary
    assert "Camera Path** accordion starts collapsed" in camera_page
    assert "Community Annotation** accordion starts collapsed" in annotation_index
    assert "unless the user explicitly enables **Autoplay**" in camera_page


def test_standard_pancreas_sample_contract_is_documented_exactly() -> None:
    page = (
        DOCS_ROOT / "user_guide" / "web_app" / "b_data_loading" / "10_standard_pancreas_dataset.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(page.split())

    for exact_value in (
        "Pancreatic endocrinogenesis (scVelo)",
        "3,696 cells",
        "3,753 genes",
        "`clusters_coarse` (5 categories)",
        "`clusters` (8 categories)",
        "`cell_type` (8 categories)",
        "`S_score`",
        "`G2M_score`",
        "33,476 exact weighted connectivity edges",
        "`velocity_umap`",
        "f6cad69fd509be44c8453205309f4c3c3c37ba34",
        "27,998 ordered genes",
        "3,774-file, 2,477,205-byte",
        "268a36b0c62f1e872bc4ebc271283fe09b9036c2b8db2b7bbf2d7d737556c865",
        "GSE132188",
        "10.1242/dev.173849",
    ):
        assert exact_value in normalized

    assert "Suo remains the catalog default" in page
    assert "compact Session statistics round both counts to **4K**" in normalized
    assert "**Showing all 3,696 points**" in normalized
    assert "default 3D embedding with **Orbit**" in normalized
    assert "Selecting 1D or 2D" in page and "**Planar**" in page
    assert "No camera path is created or played automatically" in normalized
    assert (
        "four source-owned observation fields (`clusters_coarse`, `clusters`, "
        "`S_score`, and `G2M_score`)"
    ) in normalized
    assert "it is not a fifth source-owned field" in normalized
    assert "category order and per-cell codes exactly equal source-owned `clusters`" in normalized
    assert "neither is sliced, padded, or copied from the 2D coordinates" in normalized
    assert "analysis-only working copy" in normalized
    assert "sources/pancreas.json" in page


def test_data_loading_matrix_uses_the_exact_prepared_picker_label() -> None:
    page = (
        DOCS_ROOT / "user_guide" / "web_app" / "b_data_loading" / "index.md"
    ).read_text(encoding="utf-8")

    assert "| Browser Prepared picker |" in page
    assert "Browser folder picker" not in page


def test_remote_and_github_url_docs_require_exact_multi_dataset_selection() -> None:
    page = (
        DOCS_ROOT
        / "user_guide"
        / "web_app"
        / "p_developer_docs"
        / "04_configuration_env_vars_and_feature_flags.md"
    ).read_text(encoding="utf-8")

    assert "?remote=https://data.example.org&dataset=study-a" in page
    assert "?github=owner/repo/path&dataset=study-a" in page
    assert "manifest containing exactly one dataset selects that sole" in page
    assert "requires `dataset=<exact-dataset-id>`" in page
    assert "does not choose an arbitrary\n  catalog entry" in page
    assert "loads its first dataset" not in page


def test_prepared_server_docs_open_the_catalog_without_guessing_a_dataset() -> None:
    pages = [
        (
            DOCS_ROOT
            / "user_guide"
            / "python_package"
            / "d_viewing_apis"
            / "05_python_serve_and_serve_anndata_quickstart.md"
        ),
        (
            DOCS_ROOT
            / "user_guide"
            / "python_package"
            / "d_viewing_apis"
            / "09_server_mode_advanced.md"
        ),
        (
            DOCS_ROOT
            / "user_guide"
            / "python_package"
            / "g_api_reference_coverage"
            / "api"
            / "server.md"
        ),
    ]
    combined = "\n".join(page.read_text(encoding="utf-8") for page in pages)
    normalized = " ".join(combined.split())

    assert "`server.viewer_url` opens the served catalog" in normalized
    assert "catalog contains exactly one dataset" in normalized
    assert "multiple datasets" in normalized
    assert "exact dataset-id selection" in normalized
    assert "never chooses the first entry" in normalized
    assert "stray or partial subdirectory rejects the root" in normalized
    assert "loads its first dataset" not in combined


def test_debug_connection_docs_use_only_the_keyed_all_dataset_report() -> None:
    markdown = "\n".join(path.read_text(encoding="utf-8") for path in _markdown_files())
    assert markdown.count("dataset_identity_probes") >= 11
    assert "dataset_identity_url" not in markdown
    assert "dataset_identity_error" not in markdown
    assert "recent frontend console warnings" not in markdown
    assert "recent frontend console warnings/errors" not in markdown
    assert "every declared dataset identity under its exact id/path" in markdown


def test_every_documented_screenshot_uses_its_intrinsic_responsive_width() -> None:
    screenshot_root = (DOCS_ROOT / "_static" / "screenshots").resolve()
    referenced: set[Path] = set()
    declared_widths: dict[Path, set[str]] = {}
    figure_count = 0

    for page in sorted(DOCS_ROOT.rglob("*.md")):
        lines = page.read_text(encoding="utf-8").splitlines()
        for index, line in enumerate(lines):
            if not line.startswith("```{figure} "):
                continue
            target = line.removeprefix("```{figure} ")
            image = (page.parent / target).resolve()
            if image.suffix.lower() != ".png":
                continue
            assert image.is_relative_to(screenshot_root), (page, target)
            assert image.is_file(), (page, target)

            closing_index = lines.index("```", index + 1)
            options = lines[index + 1 : closing_index]
            alt_options = [
                option.removeprefix(":alt: ") for option in options if option.startswith(":alt: ")
            ]
            width_options = [
                option.removeprefix(":width: ")
                for option in options
                if option.startswith(":width: ")
            ]
            assert len(alt_options) == 1 and alt_options[0].strip()
            assert len(width_options) == 1

            intrinsic_width, _ = _validated_png_dimensions(image)
            expected_width = f"{intrinsic_width}px"
            assert width_options[0] == expected_width, (
                page.relative_to(REPOSITORY_ROOT),
                target,
                intrinsic_width,
                width_options[0],
            )

            referenced.add(image)
            declared_widths.setdefault(image, set()).add(width_options[0])
            figure_count += 1

    screenshots = {path.resolve() for path in screenshot_root.rglob("*.png")}
    assert len(screenshots) >= 38
    assert figure_count >= 97
    assert referenced == screenshots
    assert all(len(widths) == 1 for widths in declared_widths.values())

    stylesheet = (DOCS_ROOT / "_static" / "custom.css").read_text(encoding="utf-8")
    assert "figure img {" in stylesheet
    assert "height: auto;" in stylesheet
    assert "max-width: 100%;" in stylesheet
