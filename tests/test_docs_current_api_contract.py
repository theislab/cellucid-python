from __future__ import annotations

import ast
import inspect
import json
import re
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
        DOCS_ROOT
        / "user_guide"
        / "python_package"
        / "c_data_preparation_api"
        / "index.md",
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

    failures: list[str] = []
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

    assert "anndata>=0.11.4,<0.12" in pyproject_text
    assert "zarr>=2.18.3,<3" in pyproject_text
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
