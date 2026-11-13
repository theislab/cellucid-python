from __future__ import annotations

import ast
import binascii
import hashlib
import inspect
import json
import os
import re
import shutil
import struct
import subprocess
import textwrap
from collections.abc import Iterator
from pathlib import Path

import h5py
import pytest
from _repository_source import CHILD_PROCESS_FENCE
from _web_root import node_executable

import cellucid
from cellucid.data_mapping import EMBEDDING_KEY_PREFIX, parse_vector_field_key

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
    from IPython.core.inputtransformer2 import TransformerManager

    # IPython syntax (``%pip``, also inside an ``if`` block) is valid notebook
    # code but not Python: IPython's own transformer turns it into the Python
    # it runs as. A notebook of another kernel (the R tutorials) holds no Python.
    to_python = TransformerManager()
    for path in _notebook_files():
        notebook = json.loads(path.read_text(encoding="utf-8"))
        language = notebook.get("metadata", {}).get("kernelspec", {}).get("language", "python")
        if language.lower() != "python":
            continue
        for cell_number, cell in enumerate(notebook.get("cells", []), start=1):
            if cell.get("cell_type") != "code":
                continue
            source = to_python.transform_cell("".join(cell.get("source", [])))
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
    signatures = {
        name: inspect.signature(public_object) for name, public_object in public_objects.items()
    }
    failures: list[str] = []

    for location, tree in _python_examples():
        for node in ast.walk(tree):
            # A name imported from the package is one the package exports.
            if isinstance(node, ast.ImportFrom) and node.module == "cellucid" and node.level == 0:
                failures.extend(
                    f"{location}:{node.lineno}: cellucid exports no {alias.name!r}"
                    for alias in node.names
                    if alias.name not in cellucid.__all__
                )
            if not isinstance(node, ast.Call):
                continue
            name = _call_name(node)
            # A call written cellucid.<name>(...) is a call of the package: the
            # name is one the package has, or the page documents a function
            # that does not exist.
            if (
                isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "cellucid"
                and name not in cellucid.__all__
            ):
                failures.append(f"{location}:{node.lineno}: cellucid.{name} is not a function of the package")
                continue
            # A method called on a viewer is the viewer's (viewer.pull(adata)
            # is AnnDataViewer.pull, not cellucid.pull(adata, source)).
            viewer_method = (
                isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "viewer"
                and callable(getattr(cellucid.AnnDataViewer, name, None))
            )
            if viewer_method:
                signature = inspect.signature(getattr(cellucid.AnnDataViewer, name))
            elif name in signatures:
                signature = signatures[name]
            else:
                continue
            keywords = {keyword.arg for keyword in node.keywords if keyword.arg is not None}
            has_expansion = any(keyword.arg is None for keyword in node.keywords)
            has_starred_args = any(isinstance(argument, ast.Starred) for argument in node.args)
            takes_any_keyword = any(
                parameter.kind is inspect.Parameter.VAR_KEYWORD for parameter in signature.parameters.values()
            )
            unknown = [] if takes_any_keyword else sorted(keywords - set(signature.parameters))
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


def test_documented_vector_declarations_follow_the_vector_field_key_rule() -> None:
    """Every documented vector declaration binds by the vector-field key rule.

    A vector field is ``obsm['<field>_<basis>']`` over ``obsm['X_<basis>']``
    (:func:`cellucid.data_mapping.parse_vector_field_key`). A documented
    ``obsm`` dictionary must therefore name the embedding each of its vector
    keys is drawn over, and a per-dimension UMAP grammar, which is not the
    rule, is not described as one.
    """
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


def test_server_and_jupyter_docs_state_exact_current_security_and_loading() -> None:
    text_by_path = {
        path.relative_to(REPOSITORY_ROOT): path.read_text(encoding="utf-8")
        for path in _markdown_files()
    }
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

# Published images that ``docs/_tooling/screenshots/captures.json`` has no record
# for, so ``node capture.mjs check`` cannot re-verify them (an empty screenshot
# would pass unseen). Every published image is a scenario capture: a new image
# gets a scenario and a record, never an allowance here.
UNPROVENANCED_SCREENSHOT_BUDGET = 0

# Records whose ``sha256`` does not describe the file they name, so the published
# bytes have no provenance even though a record exists. Recapturing the scenario
# is the fix; correcting the digest in place would be inventing provenance for an
# image nobody can show came from that run.
DIVERGENT_SCREENSHOT_RECORDS: frozenset[str] = frozenset()

# What a record holds, as ``capture.mjs`` writes it: the image (``sha256``, its
# sizes and format), how it was taken, the Cellucid version of the application
# it shows (``cellucid_version``) and the browser. ``quality`` is a WebP's.
SCREENSHOT_RECORD_FIELDS = frozenset({
    "id", "topic", "name", "describes", "dataset", "theme", "viewport", "device_scale_factor", "crop",
    "captured_css", "native_px", "emitted_px", "format", "cursor", "rings", "settled", "uses_timeout",
    "page_errors", "sha256", "cellucid_version", "browser", "captured_at",
})  # fmt: skip


def _version_tuple(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


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
        assert set(record) - {"quality"} == SCREENSHOT_RECORD_FIELDS, (
            capture_id,
            sorted(set(record) ^ SCREENSHOT_RECORD_FIELDS),
        )
        for field in ("topic", "name", "sha256", "cellucid_version"):
            assert record.get(field), (capture_id, field)
        # The release the image shows: a Cellucid version, and never one
        # newer than this repository's.
        assert re.fullmatch(r"\d+\.\d+\.\d+", record["cellucid_version"]), (capture_id, record["cellucid_version"])
        assert _version_tuple(record["cellucid_version"]) <= _version_tuple(cellucid.__version__), (
            capture_id,
            record["cellucid_version"],
        )
        extension = record.get("format", "png")
        assert extension in {"png", "webp"}, (capture_id, extension)
        image = (screenshot_root / record["topic"] / f"{record['name']}.{extension}").resolve()
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
    # WebP captures are new, so none of them is allowed to lack a record.
    webp = {path.resolve() for path in screenshot_root.rglob("*.webp")}
    assert webp <= recorded, sorted(
        str(path.relative_to(screenshot_root)) for path in webp - recorded
    )


def test_the_capture_tool_emits_only_from_a_web_tree_without_uncommitted_changes(tmp_path: Path) -> None:
    """``capture.mjs --emit`` writes the published images and their records,
    and the application it serves is the web repository's files on disk: it
    refuses, before anything starts, while that tree has uncommitted changes."""
    node = node_executable("the screenshot capture tool")
    git = shutil.which("git")
    if git is None:
        pytest.skip("git is not installed: the capture tool's check of the web tree cannot run.")
    tooling = DOCS_ROOT / "_tooling" / "screenshots"
    manifest_before = (tooling / "captures.json").read_bytes()
    scenario = next(iter(json.loads(manifest_before)["captures"]))

    web = tmp_path / "cellucid"
    web.mkdir()
    (web / "index.html").write_text('<meta name="cellucid-version" content="0.0.1" />\n', encoding="utf-8")
    # The first thing the tool asks the web tree once it goes on is its ports: this stand-in records that, and fails.
    started = tmp_path / "web-tree-asked"
    (web / "scripts").mkdir()
    (web / "scripts" / "browser-test-ports.mjs").write_text(
        "import { writeFileSync } from 'node:fs';\n"
        "export function resolveBrowserTestPorts() {\n"
        "  writeFileSync(process.env.CELLUCID_TEST_STARTED_MARKER, 'asked');\n"
        "  throw new Error('a stand-in holds no app');\n"
        "}\n",
        encoding="utf-8",
    )

    def run_git(*arguments: str) -> None:
        subprocess.run(
            [git, "-c", "user.name=test", "-c", "user.email=test@example.invalid", *arguments],
            cwd=web, check=True, capture_output=True,
        )

    run_git("init", "--quiet")
    run_git("add", "index.html", "scripts")
    run_git("commit", "--quiet", "-m", "a web tree")

    def emit() -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [node, str(tooling / "capture.mjs"), "shoot", scenario, "--emit"],
            env={**os.environ, "CELLUCID_WEB_ROOT": str(web), "CELLUCID_TEST_STARTED_MARKER": str(started)},
            capture_output=True, text=True, encoding="utf-8", check=False, timeout=CHILD_PROCESS_FENCE,
        )

    # A committed tree is not refused: the run goes on to the web tree, which this stand-in fails.
    committed = emit()
    assert committed.returncode != 0
    assert started.exists(), committed.stderr

    # A tree with uncommitted changes is refused before anything starts or is written, and the paths are named.
    started.unlink()
    (web / "index.html").write_text('<meta name="cellucid-version" content="0.0.2" />\n', encoding="utf-8")
    (web / "scratch.js").write_text("// not committed\n", encoding="utf-8")
    refused = emit()
    assert refused.returncode != 0
    assert not started.exists(), refused.stderr
    assert "index.html" in refused.stderr and "scratch.js" in refused.stderr
    assert (tooling / "captures.json").read_bytes() == manifest_before


def _webp_dimensions(path: Path) -> tuple[int, int]:
    payload = path.read_bytes()
    assert payload[:4] == b"RIFF" and payload[8:12] == b"WEBP", path
    chunk = payload[12:16]
    if chunk == b"VP8X":
        width = int.from_bytes(payload[24:27], "little") + 1
        height = int.from_bytes(payload[27:30], "little") + 1
        return width, height
    if chunk == b"VP8 ":
        assert payload[23:26] == b"\x9d\x01\x2a", path
        width, height = struct.unpack("<HH", payload[26:30])
        return width & 0x3FFF, height & 0x3FFF
    assert chunk == b"VP8L" and payload[20] == 0x2F, path
    bits = int.from_bytes(payload[21:25], "little")
    return (bits & 0x3FFF) + 1, ((bits >> 14) & 0x3FFF) + 1


HTML_IMAGE = re.compile(r"<img\b[^>]*>", re.DOTALL)
HTML_ATTRIBUTE = re.compile(r'([a-z]+)="([^"]*)"')


def test_every_page_placed_image_declares_its_intrinsic_size() -> None:
    """A page that places a WebP capture itself (the home page) states the
    file's own size, so the browser reserves its box before it loads, names
    only files that exist, and leaves no WebP capture unreferenced."""
    screenshot_root = (DOCS_ROOT / "_static" / "screenshots").resolve()
    referenced: set[Path] = set()
    for page in sorted(DOCS_ROOT.rglob("*.md")):
        if "_build" in page.parts:
            continue
        for tag in HTML_IMAGE.findall(page.read_text(encoding="utf-8")):
            attributes = dict(HTML_ATTRIBUTE.findall(tag))
            if "_static/screenshots/" not in attributes.get("src", ""):
                continue
            image = (page.parent / attributes["src"]).resolve()
            assert image.suffix == ".webp" and image.is_file(), (page, attributes["src"])
            assert attributes.get("alt", "").strip(), (page, attributes["src"])
            assert (int(attributes["width"]), int(attributes["height"])) == _webp_dimensions(
                image
            ), (page, attributes["src"])
            referenced.add(image)
            for candidate in attributes.get("srcset", "").split(","):
                if not candidate.strip():
                    continue
                source, descriptor = candidate.split()
                variant = (page.parent / source).resolve()
                assert variant.is_file(), (page, source)
                assert descriptor == f"{_webp_dimensions(variant)[0]}w", (page, source)
                referenced.add(variant)

    published = {path.resolve() for path in screenshot_root.rglob("*.webp")}
    assert referenced == published, sorted(
        str(path.relative_to(screenshot_root)) for path in published ^ referenced
    )


# Every page that enumerates the viewer's hook decorators, and the name of the
# heading or lead-in that makes the list read as complete. A hook added to
# ``AnnDataViewer`` has to reach all of them: with no check tying these lists to
# the class, a page silently lacks one.
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

# The web build reports a refused command as ``command_error``
# (``data/viewer-host-source.js`` ``_reportMessageFailure``) whenever the host's
# ``/_cellucid/protocol`` lists the event, and the Python and R hosts list it.
# A page that still hedges the hook as arriving "only from the web build that
# emits it" tells a reader to doubt a working hook.
COMMAND_ERROR_HEDGE = re.compile(
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
        "on_history_change",
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
        if COMMAND_ERROR_HEDGE.search(unwrapped) is not None:
            failures.append(
                f"{relative}: hedges on_command_error as arriving only from a web build "
                "that emits command_error; the current web build emits it"
            )

    assert not failures, "\n".join(failures)


# The pages that print the route's body. A capability list written out in prose
# is a claim about what the server answers, so it is checked against the answer
# rather than proof-read.
PROTOCOL_BODY_PAGES = (
    "docs/user_guide/python_package/h_developer_docs/"
    "09_server_mode_architecture_endpoints_and_security.md",
    "docs/user_guide/python_package/h_developer_docs/"
    "11_hooks_events_protocol_and_schema.md",
)


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
        keys.extend(re.findall(r"`([A-Za-z0-9_]+)`", line.split("|")[1]))
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


#: The records the packages keep in ``uns['cellucid']`` on the object a person
#: works with (an apply, a pull); ``prepare()`` never writes them.
PACKAGE_RECORD_KEYS = {"history_applied", "gene_sets", "pulled"}


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
    import pandas as pd

    text = SPECIFICATION.read_text(encoding="utf-8")
    top = _table_keys(text.split("## 2. `uns['cellucid']`", 1)[1])
    defaults = _table_keys(text.split("### `defaults`", 1)[1])
    assert top == [
        "schema", "id", "name", "description", "citation", "source", "defaults", "outlierQuantiles", "totalCounts",
        "priority", "history_applied", "gene_sets", "pulled",
    ]
    assert defaults == [
        "embedding1d",
        "embedding2d",
        "embedding3d",
        "knn",
        "modality",
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
        # Counts, so each cell's total is recorded (totalCounts).
        gene_expression=rng.poisson(2.0, size=(6, 3)).astype(np.float32),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
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
    assert set(written) <= set(top) - PACKAGE_RECORD_KEYS
    assert set(written["defaults"]) <= set(defaults)
    # Everything prepare() can write is written here, except the two defaults
    # it leaves to the viewer (knn, expression), the two it needs obs/var for,
    # and the modality, which only a MuData has.
    assert set(written) == set(top) - PACKAGE_RECORD_KEYS
    assert set(written["defaults"]) == set(defaults) - {
        "knn",
        "expression",
        "obsKeys",
        "geneIdColumn",
        "modality",
    }
