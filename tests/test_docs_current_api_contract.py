from __future__ import annotations

import ast
import binascii
import hashlib
import inspect
import json
import re
import struct
import textwrap
from collections.abc import Iterator
from pathlib import Path

import h5py

import cellucid
from cellucid.data_mapping import EMBEDDING_KEY_PREFIX, parse_vector_field_key
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


def _without_doctest_output(source: str) -> str:
    """The code of a doctest-style example: its ``>>>``/``...`` lines, unprompted.

    An interactive transcript shows a call and what it printed; only the call
    is Python, and only the call is checked against the current API.
    """
    lines = source.splitlines()
    if not lines or not lines[0].startswith(">>>"):
        return source
    return "\n".join(
        line[4:] for line in lines if line.startswith((">>> ", "... ")) or line in (">>>", "...")
    )


def _markdown_python_sources() -> Iterator[tuple[str, str]]:
    for path in _markdown_files():
        text = path.read_text(encoding="utf-8")
        for match in PYTHON_FENCE_PATTERN.finditer(text):
            if match.group("info").strip().lower() not in PYTHON_FENCE_INFO:
                continue
            source = textwrap.dedent(match.group("source"))
            line = text[: match.start("source")].count("\n") + 1
            yield f"{path.relative_to(REPOSITORY_ROOT)}:{line}", _without_doctest_output(source)


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


_FENCE = re.compile(r"^(?P<ticks>`{3,})(?P<info>.*)$")

# A transcript that documents the *missing-identity failure* has to show the
# command without its identity flags -- that incomplete command is the whole
# subject of the example. The exemption is the failure the CLI itself prints for
# exactly that mistake, so it cannot be spelled over a command that is being
# recommended rather than diagnosed.
_MISSING_IDENTITY_REPORT = re.compile(
    r"^Error: (?:dataset_name and dataset_id are|dataset_(?:name|id) is) required\b"
)


def _fenced_block_of_each_line(text: str) -> list[int | None]:
    """Return, per line, the index of the fenced block it belongs to."""
    membership: list[int | None] = []
    open_ticks: str | None = None
    block_index = -1
    for line in text.splitlines():
        fence = _FENCE.match(line)
        if open_ticks is None:
            if fence is not None:
                open_ticks = fence["ticks"]
                block_index += 1
            membership.append(None)
            continue
        if fence is not None and len(fence["ticks"]) >= len(open_ticks) and not fence["info"].strip():
            open_ticks = None
            membership.append(None)
            continue
        membership.append(block_index)
    return membership


def test_cli_examples_use_only_the_exact_current_contract() -> None:
    failures: list[str] = []
    for path in [*_markdown_files(), *_notebook_files()]:
        text = path.read_text(encoding="utf-8")
        relative = path.relative_to(REPOSITORY_ROOT)
        if "--no-backed" in text or "--backed" in text:
            failures.append(f"{relative}: documents a removed backed-mode CLI option")

        logical_text = re.sub(r"\\\s*\n\s*", " ", text)
        lines = logical_text.splitlines()
        membership = _fenced_block_of_each_line(logical_text)
        diagnostic_blocks = {
            block
            for block, line in zip(membership, lines, strict=True)
            if block is not None and _MISSING_IDENTITY_REPORT.match(line.strip())
        }
        # A file written by ``cellucid.prepare()`` declares its own identity in
        # ``uns['cellucid']``, so a block that says it serves one may omit the
        # flags; it has to say so, in the block, where the reader copies from.
        prepared_blocks = {
            block
            for block, line in zip(membership, lines, strict=True)
            if block is not None and "cellucid.prepare()" in line
        }
        for line_number, line in enumerate(lines, start=1):
            # A command names its store after the subcommand; a line that only
            # mentions `cellucid serve` (a quoted message listing where a file
            # goes) is not a command.
            if "cellucid serve" not in line or not re.search(
                r"\.(?:h5ad|zarr)\b", line.split("cellucid serve", 1)[1]
            ):
                continue
            # The identity flags are optional for a store that declares
            # uns['cellucid']['id'] and ['name'], so prose may name the bare
            # command. A command in a code block is one a reader copies onto
            # their own file, which usually declares no identity, so it names
            # the dataset.
            if membership[line_number - 1] is None:
                continue
            if membership[line_number - 1] in diagnostic_blocks | prepared_blocks:
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
    assert "rejects duplicate gene names" not in vector_page
    assert "var_names_make_unique()" in vector_page
    assert "embedding_keys" in vector_page
    assert "knn_key" in vector_page
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
        Path("docs/user_guide/python_package/installation.md")
    ]
    assert r".\.venv\Scripts\Activate.ps1" in installation_page
    assert r".\.venv\Scripts\activate.bat" in installation_page

    anndata_server_source = (REPOSITORY_ROOT / "src/cellucid/anndata_server.py").read_text(
        encoding="utf-8"
    )
    assert "one-time per web build" not in anndata_server_source
    assert '"Viewer UI cache"' not in anndata_server_source
    assert '"Viewer UI generation"' in anndata_server_source
    assert '"establishing exact configured source"' in anndata_server_source


#: The prepared-export world Cellucid removed: its server, its viewer, its
#: on-disk format and its knobs. The Python guide, the README and the site home
#: describe a Cellucid-ready AnnData instead, so none of these may come back.
REMOVED_EXPORT_VOCABULARY = (
    "CellucidViewer",
    "CellucidServer",
    "generate_datasets_manifest",
    "var_quantization",
    "obs_continuous_quantization",
    "obs_categorical_dtype",
    "latent_space=",
    "out_dir=",
    "dataset_identity.json",
    "obs_manifest.json",
    "var_manifest.json",
    "connectivity_manifest.json",
    "edges.src.bin",
    "/?source=remote",
    "Embedding(",
    "export folder",
    "export directory",
    "exported directory",
    "prepared export",
)
#: A call of the removed ``show(...)``/``serve(...)``/``pack(...)`` helpers,
#: and not ``show_anndata(...)``/``serve_anndata(...)``.
REMOVED_EXPORT_CALL = re.compile(r"(?<![\w.])(?:cellucid\.)?(?:show|serve|pack)\(")


def test_python_guide_describes_no_prepared_export_contract() -> None:
    pages = [
        README,
        DOCS_ROOT / "index.md",
        *sorted((DOCS_ROOT / "user_guide" / "python_package").rglob("*.md")),
        *sorted((DOCS_ROOT / "user_guide" / "python_package").rglob("*.ipynb")),
    ]
    failures: list[str] = []
    for path in pages:
        if "_build" in path.parts:
            continue
        text = path.read_text(encoding="utf-8")
        relative = path.relative_to(REPOSITORY_ROOT)
        lowered = text.lower()
        for term in REMOVED_EXPORT_VOCABULARY:
            if term.lower() in lowered:
                failures.append(f"{relative}: describes removed {term!r}")
        for match in REMOVED_EXPORT_CALL.finditer(text):
            line = text[: match.start()].count("\n") + 1
            failures.append(f"{relative}:{line}: calls removed {match.group(0)!r}")
    assert not failures, "\n".join(failures)


def test_documented_vector_declarations_follow_the_vector_field_key_rule() -> None:
    """Every documented vector declaration binds by the vector-field key rule.

    A vector field is ``obsm['<field>_<basis>']`` over ``obsm['X_<basis>']``
    (:func:`cellucid.data_mapping.parse_vector_field_key`). A documented
    ``obsm`` dictionary must therefore name the embedding each of its vector
    keys is drawn over, and the removed per-dimension UMAP grammar is not
    described as the current rule.
    """
    all_text = "\n".join(
        path.read_text(encoding="utf-8") for path in [*_markdown_files(), *_notebook_files()]
    )
    stale_contract_phrases = [
        "Implicit keys are also supported",
        "**Fallback (implicit key",
        "Style B: implicit dimension",
        "explicit_dim_suffix=True",
        "explicit_dim_suffix=False",
        "<field>_umap_<d>d",
        "<field>_umap_<n>d",
        "<field>_umap_<1|2|3>d",
    ]
    for phrase in stale_contract_phrases:
        assert phrase not in all_text, phrase

    failures: list[str] = []
    for location, tree in _python_examples():
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            for keyword in node.keywords:
                # prepare() binds each vector field to the embedding it is
                # drawn over, so every value of its vector_fields= mapping is
                # a mapping of 'X_<basis>' embedding keys to arrays.
                if (
                    _call_name(node) == "prepare"
                    and keyword.arg == "vector_fields"
                    and isinstance(keyword.value, ast.Dict)
                ):
                    for value in keyword.value.values:
                        if not isinstance(value, ast.Dict):
                            failures.append(
                                f"{location}:{value.lineno}: prepare vector field is not "
                                "bound to embedding keys"
                            )
                            continue
                        for key in value.keys:
                            if (
                                isinstance(key, ast.Constant)
                                and isinstance(key.value, str)
                                and not key.value.startswith(EMBEDDING_KEY_PREFIX)
                            ):
                                failures.append(
                                    f"{location}:{key.lineno}: prepare binds a vector field "
                                    f"to {key.value!r}, which is not an 'X_<basis>' key"
                                )
                # An obsm dictionary that declares a vector-looking key next to
                # embeddings must name the embedding that key is drawn over.
                if keyword.arg != "obsm" or not isinstance(keyword.value, ast.Dict):
                    continue
                keys = [
                    key.value
                    for key in keyword.value.keys
                    if isinstance(key, ast.Constant) and isinstance(key.value, str)
                ]
                if not any(key.startswith(EMBEDDING_KEY_PREFIX) for key in keys):
                    continue
                for key in keys:
                    if key.startswith(("velocity", "T_fwd", "T_bwd")) and (
                        parse_vector_field_key(key, set(keys)) is None
                    ):
                        failures.append(
                            f"{location}:{keyword.value.lineno}: obsm vector {key!r} names no "
                            f"embedding among {sorted(keys)}"
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

    # A field is listed by the <field> part of '<field>_<basis>', so a launch
    # default names `velocity` or `T_fwd`, never a full obsm key.
    assert 'vector_field_default="velocity"' in page_text[relative_paths[0]]
    assert 'vector_field_default="velocity"' in page_text[relative_paths[1]]
    assert 'default_vector_field = "velocity"' in page_text[relative_paths[1]]
    assert 'vector_field_default="velocity"' in page_text[relative_paths[2]]
    assert 'vector_field_default="T_fwd"' in page_text[relative_paths[2]]
    assert not re.search(
        r'\bvector_field_default\s*=\s*["\'](?:velocity|T_fwd)_umap\w*["\']',
        combined,
    )
    assert 'default_vector_field = "velocity_umap"' not in combined


def test_documented_embedding_requirements_use_exact_current_keys() -> None:
    """Two paths name embeddings, and each page must state its own rule.

    ``prepare()`` takes every embedding as an array in ``embeddings={...}``,
    stored whole under its key, and names a dimension's default with
    ``embedding_defaults``. The readers -- the Python serve path and the
    browser's own file picker -- read ``obsm`` keys, and a plain ``X_umap`` at
    the dimension its column count states, so their pages must say so.
    """
    reader_pages = [
        DOCS_ROOT
        / "user_guide"
        / "web_app"
        / "b_data_loading"
        / "03_browser_file_picker_tutorial.md",
    ]
    preparation = DOCS_ROOT / "user_guide" / "python_package" / "c_data_preparation_api"
    preparation_pages = [
        preparation / "01_prepare_overview.md",
        preparation / "03_embeddings_and_coordinates.md",
        preparation / "index.md",
    ]

    failures: list[str] = []
    for path in preparation_pages:
        text = path.read_text(encoding="utf-8")
        relative = path.relative_to(REPOSITORY_ROOT)
        if "embeddings=" not in text and "embeddings={" not in text:
            failures.append(f"{relative}: omits the embeddings= argument")
        for removed in ("X_umap_1d=", "X_umap_2d=", "X_umap_3d=", "Embedding("):
            if removed in text:
                failures.append(f"{relative}: documents removed argument {removed!r}")

    for path in reader_pages:
        text = path.read_text(encoding="utf-8")
        relative = path.relative_to(REPOSITORY_ROOT)
        for key in ("X_umap_1d", "X_umap_2d", "X_umap_3d"):
            if key not in text:
                failures.append(f"{relative}: omits current embedding key {key!r}")
        if "X_umap" not in text or "column count" not in text:
            failures.append(
                f"{relative}: must state that a plain X_umap is read at the "
                "dimension its column count states"
            )

    overview = preparation_pages[0].read_text(encoding="utf-8")
    for container in ("`.h5ad`", "`.zarr`", "`.zarr.zip`", "uns['cellucid']"):
        assert container in overview
    embeddings_page = preparation_pages[1].read_text(encoding="utf-8")
    for exact in ("embedding_defaults", "latent_key", "normalize_embeddings", "every column"):
        assert exact in embeddings_page

    assert not failures, "\n".join(failures)


def test_connectivity_docs_state_the_stored_graph_contract() -> None:
    """``prepare(connectivities=...)`` stores one symmetric CSR float32 graph."""
    page = (
        DOCS_ROOT
        / "user_guide"
        / "python_package"
        / "c_data_preparation_api"
        / "07_connectivities_knn_graph.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(page.split())
    for exact in (
        "`obsp['connectivities']`, a CSR float32 matrix",
        "**exactly symmetric**",
        "zero on the diagonal",
        "Connectivity matrix must be exactly symmetric.",
    ):
        assert exact in normalized, exact


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

    python_server_page = text_by_path[
        Path(
            "docs/user_guide/python_package/d_viewing_apis/"
            "05_python_serve_anndata_quickstart.md"
        )
    ]
    advanced_server_page = text_by_path[
        Path("docs/user_guide/python_package/d_viewing_apis/09_server_mode_advanced.md")
    ]
    for page in (python_server_page, advanced_server_page):
        assert "anndata=true" in page
        assert "first paint" in page
        assert "source=remote" not in page

    failures: list[str] = []
    # ``?data=http://127.0.0.1:8000/`` names a static host for the web app to
    # read from, not a Cellucid server's launch, so it is not a Viewer URL.
    bare_loopback_viewer = re.compile(
        r"(?<!data=)http://127\.0\.0\.1:(?:<port>|\d+)/"
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
                    f"{path}:{line_number}: server launch "
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
        "command_error",
        "data_mapping",
    ):
        assert f"### `{event_type}`" in event_page
    assert "Unknown event types" in event_page
    assert "rejected before any hook runs" in event_page
    assert "one of nine closed event records" in event_page
    assert "@viewer.on_command_error" in event_page
    # A command that failed must not read as one that applied.
    assert "did not take effect in the viewer" in event_page

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


def test_embedding_docs_state_that_the_viewer_normalizes() -> None:
    """``prepare()`` stores coordinates as given; the browser normalizes them."""
    page = (
        DOCS_ROOT
        / "user_guide"
        / "python_package"
        / "c_data_preparation_api"
        / "03_embeddings_and_coordinates.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(page.split())
    assert "Nothing is normalized or cropped." in normalized
    assert "The viewer scales each embedding into its view" in normalized
    assert "normalizes each embedding on its own" not in normalized


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
    assert "default light studio floor" in quick_tour
    # The sidebar footer's build identity bumps on every web release — seven
    # times in the five days after ``2026-07-27.1``. Pinning a literal here made
    # the assertion enforce staleness: it kept passing while both pages, and the
    # screenshot they share, described a build nobody was running, and it turned
    # every web change into a documentation change. The durable claim is *where*
    # to read the value and that it is yours, not ours; this forbids a literal
    # coming back.
    assert not re.search(r"Build \d{4}-\d{2}-\d{2}\.\d+", combined), combined
    assert "changes with every web release" in " ".join(quick_tour.split())
    assert "read yours off your own footer" in " ".join(glossary.split())
    assert "exact value in bug reports" in glossary
    assert "Camera path** accordion starts collapsed" in camera_page
    assert "Community annotation** accordion starts collapsed" in annotation_index
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
        "one vector field, `velocity` (the `obsm` arrays `velocity_umap_1d`, "
        "`velocity_umap` and `velocity_umap_3d`",
        "| `X_umap_1d` | 1 | derived | yes (1D default) | `velocity_umap_1d` |",
        "| `X_umap` | 2 | source | yes (2D default) | `velocity_umap` |",
        "| `X_umap_3d` | 3 | derived | yes (3D default) | `velocity_umap_3d` |",
        "| `X_pca` | 50 | source | yes, components 1 and 2 as a 2D view, enabled by the starting view",
        "selector reports `velocity`",
        "f6cad69fd509be44c8453205309f4c3c3c37ba34",
        "27,998 ordered genes",
        "SHA-256 store digest covering every file of `anndata/pancreas.zarr`",
        "recorded in `sources/pancreas.json` rather than repeated here",
        "GSE132188",
        "10.1242/dev.173849",
    ):
        assert exact_value in normalized

    assert "Suo remains the catalog default" in page
    # The file stores the field under scVelo's names for it on each UMAP
    # (`<field>_<basis>`), which the viewer lists as the one field `velocity`.
    assert "`velocity_umap` field" not in page
    assert "velocity_umap_2d" not in page
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


def test_data_loading_matrix_uses_the_exact_local_picker_labels() -> None:
    page = (
        DOCS_ROOT / "user_guide" / "web_app" / "b_data_loading" / "index.md"
    ).read_text(encoding="utf-8")

    assert "| A Zarr store folder (`.zarr`), for example from `prepare()` | **Open… › Zarr folder…**, or drop the folder |" in page
    assert "| A portable `.zarr.zip` or `.zip` | **Open… › File…** |" in page
    assert "| A `.h5ad` of any size, and no Python | **Open… › File…** |" in page
    assert "Prepared" not in page
    assert "Browser folder picker" not in page
    for removed in ("**H5AD**", "**Zarr ZIP**", "**Files**", "Open from link", "**Connect**"):
        assert removed not in page, removed


def test_remote_and_github_url_docs_require_exact_multi_dataset_selection() -> None:
    page = (
        DOCS_ROOT
        / "user_guide"
        / "web_app"
        / "p_developer_docs"
        / "04_configuration_env_vars_and_feature_flags.md"
    ).read_text(encoding="utf-8")

    assert "?remote=https://data.example.org&dataset=study-a" in page
    assert "?data=gh:owner/repo/path&dataset=study-a" in page
    assert "dataset list containing exactly one dataset selects that sole" in page
    assert "requires `dataset=<exact-dataset-id>`" in page
    assert "does not choose an arbitrary\n  catalog entry" in page
    assert "loads its first dataset" not in page


def test_server_docs_open_the_one_served_dataset() -> None:
    pages = [
        (
            DOCS_ROOT
            / "user_guide"
            / "python_package"
            / "d_viewing_apis"
            / "05_python_serve_anndata_quickstart.md"
        ),
        (
            DOCS_ROOT
            / "user_guide"
            / "python_package"
            / "d_viewing_apis"
            / "09_server_mode_advanced.md"
        ),
    ]
    for page in pages:
        normalized = " ".join(page.read_text(encoding="utf-8").split())
        assert "`server.viewer_url`" in normalized, page
        assert "opens the one dataset the server lists at `/_cellucid/datasets`" in normalized, page
        assert "loads its first dataset" not in normalized, page


def test_debug_connection_docs_use_only_the_keyed_store_probe_report() -> None:
    markdown = "\n".join(path.read_text(encoding="utf-8") for path in _markdown_files())
    python_guide = "\n".join(
        path.read_text(encoding="utf-8")
        for path in sorted((DOCS_ROOT / "user_guide" / "python_package").rglob("*.md"))
    )
    assert python_guide.count("anndata_store_probes") >= 11
    assert "dataset_identity_probes" not in python_guide
    assert "dataset_identity_url" not in markdown
    assert "dataset_identity_error" not in markdown
    assert "recent frontend console warnings" not in markdown
    assert "recent frontend console warnings/errors" not in markdown


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


# Published images that ``docs/_tooling/screenshots/captures.json`` has no record
# for. They predate the capture tool — several were made by a script that was
# never committed — so ``node capture.mjs check`` cannot re-verify them, which is
# how an empty screenshot once shipped. Backfilling one means authoring a
# scenario and recapturing, not editing this number.
#
# **This number may only go down.** It is a ratchet: it stops the backlog growing
# while it is worked off, and reaching 0 is what closes the gap.
UNPROVENANCED_SCREENSHOT_BUDGET = 62

# Records whose ``sha256`` does not describe the file they name, so the published
# bytes have no provenance even though a record exists. Recapturing the scenario
# is the fix; correcting the digest in place would be inventing provenance for an
# image nobody can show came from that run.
DIVERGENT_SCREENSHOT_RECORDS: frozenset[str] = frozenset()


def test_screenshot_provenance_records_describe_the_published_images() -> None:
    """``captures.json`` must describe the bytes actually published, and the
    set of images it cannot account for must only ever shrink."""
    screenshot_root = (DOCS_ROOT / "_static" / "screenshots").resolve()
    manifest = json.loads(
        (DOCS_ROOT / "_tooling" / "screenshots" / "captures.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["version"] == 1

    recorded: set[Path] = set()
    divergent: set[str] = set()
    for capture_id, record in manifest["captures"].items():
        for field in ("topic", "name", "sha256", "web_repository_revision"):
            assert record.get(field), (capture_id, field)
        image = (screenshot_root / record["topic"] / f"{record['name']}.png").resolve()
        if record["topic"] == "_proof":
            # ``_proof`` scenarios demonstrate that the tool works and are never
            # emitted; an on-disk file for one would be an orphan PNG.
            assert not image.exists(), capture_id
            continue
        assert image.is_file(), (capture_id, image)
        recorded.add(image)
        digest = hashlib.sha256(image.read_bytes()).hexdigest()
        if digest != record["sha256"]:
            divergent.add(capture_id)

    assert divergent == set(DIVERGENT_SCREENSHOT_RECORDS), sorted(
        divergent.symmetric_difference(DIVERGENT_SCREENSHOT_RECORDS)
    )

    published = {path.resolve() for path in screenshot_root.rglob("*.png")}
    unprovenanced = published - recorded
    assert len(unprovenanced) <= UNPROVENANCED_SCREENSHOT_BUDGET, sorted(
        str(path.relative_to(screenshot_root)) for path in unprovenanced
    )


# Every page that enumerates the viewer's hook decorators, and the name of the
# heading or lead-in that makes the list read as complete. A hook added to
# ``AnnDataViewer`` has to reach all of them: the web-app Jupyter tutorial silently
# lost ``on_command_error`` because no check tied these lists to the class.
HOOK_ENUMERATION_PAGES = (
    "docs/user_guide/web_app/b_data_loading/05_jupyter_tutorial.md",
    "docs/user_guide/python_package/e_jupyter_hooks/"
    "01_overview_bidirectional_communication.md",
    "docs/user_guide/python_package/e_jupyter_hooks/16_reference.md",
    "docs/user_guide/python_package/b_concepts_mental_models/"
    "01_what_is_a_viewer_object.md",
    "docs/user_guide/python_package/h_developer_docs/"
    "10_jupyter_embedding_architecture.md",
)

# ``command_error`` is accepted by the wire protocol ahead of any web build that
# emits it, so every page naming the hook also has to say the event does not
# arrive yet. Without that, a reader registers a handler and reads its silence
# as a viewer that never refuses a command.
COMMAND_ERROR_AVAILABILITY = re.compile(
    r"(?:only the web build that emits|sent only by the web build that emits)"
    r" `command_error`",
    re.IGNORECASE,
)


def test_every_documented_hook_list_matches_the_viewer_decorator_surface() -> None:
    from cellucid.jupyter import AnnDataViewer

    hooks = {
        name
        for name, member in inspect.getmembers(AnnDataViewer, lambda value: isinstance(value, property))
        if name.startswith("on_")
    }
    assert hooks == {
        "on_selection",
        "on_hover",
        "on_click",
        "on_ready",
        "on_command_error",
        "on_data_mapping",
        "on_message",
    }, hooks

    failures: list[str] = []
    for relative in HOOK_ENUMERATION_PAGES:
        page = (REPOSITORY_ROOT / relative).read_text(encoding="utf-8")
        # The sentence is wrapped differently on every page, so match it against
        # the page with its line breaks collapsed rather than as written.
        unwrapped = " ".join(page.split())
        for hook in sorted(hooks):
            if hook not in page:
                failures.append(f"{relative}: hook list omits {hook}")
        if COMMAND_ERROR_AVAILABILITY.search(unwrapped) is None:
            failures.append(
                f"{relative}: names on_command_error without saying the event arrives "
                "only from the web build that emits command_error"
            )

    assert not failures, "\n".join(failures)


# Every page that enumerates the routes the Python servers answer, and so reads
# as complete. ``/_cellucid/protocol`` is how a web build discovers what this
# installation accepts, which makes an endpoint list that omits it actively
# misleading rather than merely short.
ENDPOINT_ENUMERATION_PAGES = (
    "docs/user_guide/python_package/h_developer_docs/"
    "09_server_mode_architecture_endpoints_and_security.md",
    "docs/user_guide/python_package/h_developer_docs/"
    "17_security_privacy_and_networking.md",
    "docs/user_guide/python_package/e_jupyter_hooks/16_reference.md",
    "docs/user_guide/python_package/e_jupyter_hooks/12_security_model.md",
    "docs/user_guide/python_package/d_viewing_apis/09_server_mode_advanced.md",
    "docs/user_guide/python_package/d_viewing_apis/"
    "04_cli_cellucid_serve_quickstart.md",
    "docs/user_guide/python_package/g_api_reference_coverage/api/server.md",
)

# The pages that print the route's body. A capability list written out in prose
# is a claim about what the server answers, so it is checked against the answer
# rather than proof-read.
PROTOCOL_BODY_PAGES = (
    "docs/user_guide/python_package/h_developer_docs/"
    "09_server_mode_architecture_endpoints_and_security.md",
    "docs/user_guide/python_package/h_developer_docs/"
    "11_hooks_events_protocol_and_schema.md",
)


def test_every_documented_endpoint_list_names_the_capability_route() -> None:
    failures = [
        f"{relative}: endpoint list omits /_cellucid/protocol"
        for relative in ENDPOINT_ENUMERATION_PAGES
        if "/_cellucid/protocol" not in (REPOSITORY_ROOT / relative).read_text(encoding="utf-8")
    ]

    assert not failures, "\n".join(failures)


def test_every_documented_capability_body_matches_the_served_document() -> None:
    from cellucid.jupyter._wire import _protocol_capability_document

    served = _protocol_capability_document()
    failures: list[str] = []
    for relative in PROTOCOL_BODY_PAGES:
        page = (REPOSITORY_ROOT / relative).read_text(encoding="utf-8")
        bodies = [
            match.group("source")
            for match in PYTHON_FENCE_PATTERN.finditer(page)
            if match.group("info").strip() == "json" and '"commands"' in match.group("source")
        ]
        if len(bodies) != 1:
            failures.append(f"{relative}: expected exactly one capability body, found {len(bodies)}")
            continue
        if json.loads(bodies[0]) != served:
            failures.append(f"{relative}: capability body does not match /_cellucid/protocol")

    assert not failures, "\n".join(failures)


PREPARATION_DOCS = DOCS_ROOT / "user_guide" / "python_package" / "c_data_preparation_api"
SPECIFICATION = PREPARATION_DOCS / "09_cellucid_ready_anndata_specification.md"


def _table_keys(section: str) -> list[str]:
    """The backticked first-column keys of the first Markdown table in ``section``."""
    keys: list[str] = []
    for line in section.splitlines():
        if not line.startswith("| `"):
            if keys:
                break
            continue
        keys.extend(re.findall(r"`([A-Za-z0-9]+)`", line.split("|")[1]))
    return keys


def test_prepare_reference_signature_is_the_runtime_signature() -> None:
    """The signature block in the API reference is ``inspect.signature(prepare)``."""
    page = (
        DOCS_ROOT / "user_guide" / "python_package" / "g_api_reference_coverage" / "api" / "prepare.md"
    ).read_text(encoding="utf-8")
    block = page.split("## Signature", 1)[1].split("```text", 1)[1].split("```", 1)[0]
    documented: list[tuple[str, str | None]] = []
    keyword_only = False
    for raw in block.splitlines():
        entry = raw.strip().rstrip(",")
        if not entry or entry.endswith("(") or entry.startswith(")"):
            continue
        if entry == "*":
            keyword_only = True
            continue
        name, _, default = entry.partition("=")
        documented.append((name, default or None))
        assert keyword_only or name == "out_path", name

    signature = inspect.signature(cellucid.prepare)
    runtime = [
        (name, None if parameter.default is inspect.Parameter.empty else repr(parameter.default))
        for name, parameter in signature.parameters.items()
    ]
    assert documented == runtime
    assert [
        name
        for name, parameter in signature.parameters.items()
        if parameter.kind is not inspect.Parameter.KEYWORD_ONLY
    ] == ["out_path"]


def test_specification_lists_exactly_the_uns_cellucid_keys_the_writer_and_viewer_use(
    tmp_path: Path,
) -> None:
    """The ``uns['cellucid']`` and ``defaults`` tables match what ``prepare()`` writes.

    The viewer's own key lists (``anndata-cellucid-metadata.js`` and
    ``requireAnnDataLaunchDefaults``) are restated here; a change to either is
    a change to the specification page and to this list together.
    """
    import anndata as ad
    import numpy as np

    text = SPECIFICATION.read_text(encoding="utf-8")
    top = _table_keys(text.split("## 2. `uns['cellucid']`", 1)[1])
    defaults = _table_keys(text.split("### `defaults`", 1)[1])
    assert top == ["schema", "id", "name", "description", "citation", "source", "defaults", "outlierQuantiles", "session"]
    assert defaults == [
        "embedding1d",
        "embedding2d",
        "embedding3d",
        "knn",
        "expression",
        "geneIdColumn",
        "latentKey",
        "vectorFieldDefault",
        "obsKeys",
        "normalizeEmbeddings",
        "centroidMinPoints",
    ]

    rng = np.random.default_rng(0)
    path = cellucid.prepare(
        tmp_path / "every-key.h5ad",
        format="h5ad",
        embeddings={"X_umap": rng.normal(size=(6, 2)), "X_umap_3d": rng.normal(size=(6, 3))},
        dataset_id="every-key",
        dataset_name="Every key",
        vector_fields={"velocity": {"X_umap": rng.normal(size=(6, 2))}},
        dataset_description="All of them.",
        citation_text="Doe 2026",
        citation_url="https://example.org",
        source_name="Synthetic",
        source_url="https://example.org/source",
        embedding_defaults={"1d": ("X_umap", [1]), "2d": "X_umap", "3d": "X_umap_3d"},
        latent_key="X_umap_3d",
        vector_field_default="velocity",
        normalize_embeddings=False,
        centroid_min_points=3,
    )
    written = ad.read_h5ad(path).uns["cellucid"]
    assert set(written) <= set(top) - {"session"}
    assert set(written["defaults"]) <= set(defaults)
    # Everything prepare() can write is written here, except the two defaults
    # it leaves to the viewer (knn, expression) and the four it needs obs/var for.
    assert set(written) == set(top) - {"session"}
    assert set(written["defaults"]) == set(defaults) - {
        "knn",
        "expression",
        "obsKeys",
        "geneIdColumn",
    }
