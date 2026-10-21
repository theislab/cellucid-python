"""Every size or count limit in the package is a real one, and says why.

A refusal constant -- a name with MAX, LIMIT, CAP, CAPACITY, BUDGET or CEILING
in it, or a byte size written as ``N * 1024 * 1024`` -- has to appear below
with its class and the reason it exists. A new one fails this test until
someone classifies it, and a limit made up for convenience has no class to
take: the classes are

- FORMAT: fixed by a file or wire format (a u32 length, a DNS label);
- PLAT: a real limit of a platform the package runs on or publishes to;
- HOST: a bound on untrusted input, tied to what that input may declare;
- PERF: a cache or scheduling size that decides what is kept, never what is
  refused;
- MSG: how much of something is quoted in a message.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCANNED = ("src/cellucid", "scripts")
NAME = re.compile(r"(^|_)(MAX|MAXIMUM|LIMIT|LIMITS|CAP|CAPACITY|BUDGET|CEILING)(_|$)")
DECLARATION = re.compile(r"^\s*(_?[A-Z][A-Z0-9_]*)\s*(?::[^=\n]+)?=(?!=)", re.M)
MIB_LITERAL = re.compile(r"^\s*(_?[A-Z][A-Z0-9_]*)\s*(?::[^=\n]+)?=\s*[^\n]*\b1024\s*\*\s*1024", re.M)
MIB_INLINE = re.compile(r"\b1024\s*\*\s*1024")
CLASSES = {"FORMAT", "PLAT", "HOST", "PERF", "MSG"}

ALLOWLIST: dict[str, tuple[str, str]] = {
    "scripts/normalize_sdist.py:MAX_NORMALIZED_SDIST_BYTES": (
        "PLAT",
        "PyPI's upload limit for one distribution file (100 MB).",
    ),
    "src/cellucid/_server_base.py:MAX_SELECTION_BODY_BYTES": (
        "FORMAT",
        "A selection is little-endian uint32 cell indices: at most 2**32 of 4 bytes.",
    ),
    "src/cellucid/_server_base.py:_MAXIMUM_HOST_LABEL_LENGTH": (
        "FORMAT",
        "A DNS label is at most 63 octets (RFC 1035).",
    ),
    "src/cellucid/_server_base.py:_MAXIMUM_HOST_NAME_LENGTH": (
        "FORMAT",
        "A DNS name is at most 253 characters in text form (RFC 1035).",
    ),
    "src/cellucid/_server_base.py:_WEB_BUILD_ID_MAX_LEN": (
        "HOST",
        "Bounds the build-id string read from a downloaded index.html "
        "(the web app writes a date and a short hash).",
    ),
    "src/cellucid/continuous_payload_diagnosis.py:NON_FINITE_EXAMPLE_LIMIT": (
        "MSG",
        "How many offending values a diagnosis quotes; every value is still checked.",
    ),
    "src/cellucid/jupyter/_base.py:_DEBUG_COMMAND_ERROR_LIMIT": (
        "MSG",
        "How many viewer command errors debug_connection() quotes.",
    ),
    "src/cellucid/jupyter/_wire.py:_COMMAND_ERROR_REASON_LIMIT": (
        "HOST",
        "Bounds a command_error reason string the viewer reports over the wire.",
    ),
    "src/cellucid/session_bundle.py:MAX_FRAME_BYTES": (
        "FORMAT",
        "A session bundle frame length is a little-endian u32.",
    ),
    "src/cellucid/session_bundle.py:MAX_SAFE_INTEGER": (
        "FORMAT",
        "The largest integer a JavaScript number holds exactly (2**53 - 1).",
    ),
    "src/cellucid/session_codecs.py:MAX_SAFE_INTEGER": (
        "FORMAT",
        "The largest integer a JavaScript number holds exactly (2**53 - 1).",
    ),
    "src/cellucid/session_codecs.py:MAX_UINT16": (
        "FORMAT",
        "The largest value of a uint16 code array.",
    ),
    "src/cellucid/session_codecs.py:MAX_UINT32": (
        "FORMAT",
        "The largest value of a uint32 code array.",
    ),
    "src/cellucid/_server_base.py:MAX_EVENT_BODY_BYTES": (
        "HOST",
        "One posted viewer event record; selections and sessions have their own routes.",
    ),
    "src/cellucid/server/_handler.py:_READ_CHUNK_BYTES": (
        "PERF",
        "The size of one read while a served artifact is streamed.",
    ),
    "src/cellucid/anndata_mapping.py:DEFAULT_MAPPING_CACHE_BYTES": (
        "PERF",
        "How many decoded mapping arrays stay cached; a larger one is read again, never refused.",
    ),
    "src/cellucid/anndata_adapter.py:GENE_EXPRESSION_CACHE_BYTES": (
        "PERF",
        "How many decoded gene columns stay cached; others are read again, never refused.",
    ),
}


def _scanned_limits() -> set[str]:
    found: set[str] = set()
    for directory in SCANNED:
        for path in sorted((ROOT / directory).rglob("*.py")):
            text = path.read_text(encoding="utf-8")
            relative = path.relative_to(ROOT).as_posix()
            for match in DECLARATION.finditer(text):
                if NAME.search(match.group(1).lstrip("_")):
                    found.add(f"{relative}:{match.group(1)}")
            for match in MIB_LITERAL.finditer(text):
                found.add(f"{relative}:{match.group(1)}")
            # A byte size written inline is a limit nobody named or classified.
            for number, line in enumerate(text.splitlines(), 1):
                if MIB_INLINE.search(line) and not MIB_LITERAL.match(line + "\n"):
                    found.add(f"{relative}:line {number}")
    return found


def test_every_size_or_count_limit_is_classified() -> None:
    unclassified = sorted(_scanned_limits() - ALLOWLIST.keys())
    assert unclassified == [], (
        "These limits have no class and reason in tests/test_size_limit_allowlist.py: "
        f"{unclassified}. A limit stays only if it is FORMAT, PLAT, HOST, PERF or MSG."
    )


def test_the_allowlist_names_only_limits_that_exist() -> None:
    stale = sorted(ALLOWLIST.keys() - _scanned_limits())
    assert stale == [], f"Remove these entries; the limits are gone: {stale}"


def test_every_entry_has_a_real_class_and_a_reason() -> None:
    for key, (kind, reason) in ALLOWLIST.items():
        assert kind in CLASSES, f"{key}: {kind!r} is not a limit class"
        assert reason.strip(), f"{key}: say why the limit exists"
