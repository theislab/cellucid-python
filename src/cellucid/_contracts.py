"""
Contracts every user-supplied identity, label, and output path must satisfy.

These rules decide what Cellucid accepts as a name it will store, draw, or
write to: a portable filename component, an identifier distinct within its
axis, text the viewer can render exactly as stored, a native scalar that was
never coerced, a canonical UTC timestamp, and a directory that is safe to
replace whole.

The module depends on nothing but the standard library, and that is
load-bearing. The export writer, the two HTTP servers, the notebook
integration, and the web-asset cache all enforce these same rules, but only the
export writer needs numpy, pandas, and scipy. Keeping the contracts here is
what lets the light modules import them without paying for the heavy ones.
"""

import os
import re
import unicodedata
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path

JsonScalar = str | bool | int | float
_PORTABLE_COMPONENT_PATTERN = re.compile(
    r"^[A-Za-z0-9][A-Za-z0-9._-]{0,179}$",
    flags=re.ASCII,
)
_CANONICAL_UTC_TIMESTAMP_PATTERN = re.compile(
    r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$",
    flags=re.ASCII,
)
# Text the viewer prints verbatim -- a category label, a dataset name, a
# description, a source line -- is rendered with the browser's default
# white-space handling, in the legend, the field selector, and every exported
# figure. Characters that carry no glyph therefore change the value without
# changing the picture: 'Liver ' and 'Liver' are two distinct categories that
# look identical on screen and in an exported SVG.
#
# All three classes below are rejected rather than repaired. Trimming would
# rewrite a scientific label the caller never asked to change, and would merge
# two distinct categories into one, silently moving cells between them.
_CONTROL_CHARACTER_PATTERN = re.compile(r"[\x00-\x1f\x7f-\x9f]")
# Zero-width characters with no meaning of their own: ZERO WIDTH SPACE, WORD
# JOINER, and the byte-order mark a spreadsheet leaves at the front of a
# UTF-8 CSV. U+200C and U+200D are deliberately absent: they join Indic,
# Persian, and emoji sequences, so banning them would reject real text.
_INVISIBLE_CHARACTER_PATTERN = re.compile("[\u200b\u2060\ufeff]")
# A surrogate code point is only meaningful as half of a pair. Python holds an
# unpaired one, and json.dumps escapes it happily, so it reaches the browser -
# where JSON.parse yields an unpaired code unit that draws as U+FFFD. The
# stored label and the drawn label differ, which is this rule's whole subject.
_SURROGATE_PATTERN = re.compile("[\ud800-\udfff]")
# Every Unicode whitespace character that is not already a control character.
_EDGE_WHITESPACE_CHARACTERS = (
    "\u0020\u00a0\u1680\u2000-\u200a\u2028\u2029\u202f\u205f\u3000"
)
_EDGE_WHITESPACE_PATTERN = re.compile(
    f"^[{_EDGE_WHITESPACE_CHARACTERS}]|[{_EDGE_WHITESPACE_CHARACTERS}]$"
)
_WHITESPACE_RUN_PATTERN = re.compile(f"[{_EDGE_WHITESPACE_CHARACTERS}]+")
_WINDOWS_RESERVED_COMPONENTS = {
    "CON",
    "PRN",
    "AUX",
    "NUL",
    *(f"COM{index}" for index in range(1, 10)),
    *(f"LPT{index}" for index in range(1, 10)),
}


def _require_portable_filename_component(
    name: object,
    *,
    label: str = "Field key",
) -> str:
    """Require one exact cross-platform ASCII filename component."""
    if not isinstance(name, str):
        raise TypeError(f"{label} must be a native string.")
    if not _PORTABLE_COMPONENT_PATTERN.fullmatch(name):
        raise ValueError(
            f"{label} {name!r} must be 1-180 ASCII bytes, start with an "
            "ASCII letter or digit, and otherwise contain only ASCII "
            "letters, digits, '.', '_', or '-'."
        )
    if name.endswith("."):
        raise ValueError(f"{label} {name!r} must not end with '.'.")
    windows_stem = name.split(".", 1)[0].upper()
    if windows_stem in _WINDOWS_RESERVED_COMPONENTS:
        raise ValueError(f"{label} {name!r} is a reserved Windows filename.")
    return name


def _require_unique_identifiers(
    values: Sequence[str],
    *,
    label: str,
) -> list[str]:
    """Require one distinct identifier per position.

    This is the identity rule, not the payload-path rule. It is what makes a
    lookup by identifier resolve exactly one row, so it holds over every
    identifier supplied, whether or not that row is among the exported ones.
    """
    seen: set[str] = set()
    unique: list[str] = []
    for value in values:
        if value in seen:
            raise ValueError(f"{label} key {value!r} is duplicated.")
        seen.add(value)
        unique.append(value)
    return unique


def _require_string_identifiers(
    values: Sequence[object],
    *,
    label: str,
) -> list[str]:
    identifiers: list[str] = []
    for index, value in enumerate(values):
        if not isinstance(value, str):
            raise TypeError(
                f"{label} identifier at position {index} must be a string, "
                f"got {type(value).__name__}."
            )
        if not value:
            raise ValueError(f"{label} identifier at position {index} must be non-empty.")
        identifiers.append(value)
    return identifiers


def _require_nonempty_string(value: object, *, label: str) -> str:
    """Require explicit text without deriving or coercing an identity."""
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a non-empty string.")
    if not value:
        raise ValueError(f"{label} must be a non-empty string.")
    return value


def _describe_character(character: str) -> str:
    """Name one character exactly, so an invisible defect can be read aloud."""
    codepoint = f"U+{ord(character):04X}"
    if 0xD800 <= ord(character) <= 0xDFFF:
        # Surrogates carry no Unicode name, and calling one a control
        # character would send the reader looking for the wrong thing.
        return f"{codepoint} (unpaired surrogate)"
    name = unicodedata.name(character, "")
    if name:
        return f"{codepoint} {name}"
    # The C0/C1 ranges carry no Unicode name at all, so say what they are.
    return f"{codepoint} (control character)"


def _display_text_defect(text: str) -> str | None:
    """Describe why one string cannot be shown verbatim, or None when it can.

    The rejected classes are the ones whose stored form and drawn form differ:
    a control character, an unpaired surrogate, a zero-width character, and
    whitespace at either edge. The first three occupy no glyph of their own or
    draw as a replacement character; the last is how 'Liver ' reaches a legend
    that reads exactly like the separate 'Liver' category beside it.

    The surrogate case is Python's alone. R strings cannot hold an unpaired
    surrogate, so its writer has nothing to check.
    """
    control = _CONTROL_CHARACTER_PATTERN.search(text)
    if control is not None:
        return f"contains {_describe_character(control.group())}"
    surrogate = _SURROGATE_PATTERN.search(text)
    if surrogate is not None:
        return f"contains {_describe_character(surrogate.group())}"
    invisible = _INVISIBLE_CHARACTER_PATTERN.search(text)
    if invisible is not None:
        return f"contains {_describe_character(invisible.group())}"
    edge = _EDGE_WHITESPACE_PATTERN.search(text)
    if edge is not None:
        position = "starts with" if edge.start() == 0 else "ends with"
        return f"{position} {_describe_character(edge.group())}"
    return None


def _require_display_text(value: object, *, label: str, allow_empty: bool) -> str:
    """Require one string the viewer can draw exactly as it is stored."""
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string.")
    if not value:
        if allow_empty:
            return value
        raise ValueError(f"{label} must be a non-empty string.")
    defect = _display_text_defect(value)
    if defect is not None:
        raise ValueError(
            f"{label} is displayed verbatim, so it must not carry characters that "
            f"have no glyph: {value!r} {defect}. Cellucid does not remove them for "
            "you, because that would change text you did not ask to change. Pass "
            "the exact text you want shown."
        )
    return value


def _require_field_identities(
    values: Sequence[object],
    *,
    label: str,
) -> list[str]:
    """Require one distinct, drawable identity per exported field on one axis.

    A payload filename is an integer index, so an identifier is never a path and
    carries no filename rule at all: no ASCII restriction, no case-collision
    rule, no reserved Windows name. What survives is what the identity is
    actually for. It names one field in the manifest, so it must be a non-empty
    string, distinct within its axis so a lookup resolves one field, and text
    the viewer can draw exactly as it is stored -- the same rule every category
    label obeys, because a gene name and a category label are shown in the same
    legend.
    """
    identifiers = _require_string_identifiers(values, label=label)
    for position, identifier in enumerate(identifiers):
        _require_display_text(
            identifier,
            label=f"{label} identifier at position {position}",
            allow_empty=False,
        )
    return _require_unique_identifiers(identifiers, label=label)


def _require_dataset_name(value: object, *, label: str = "dataset_name") -> str:
    """Require one exact unpadded, nonblank human-readable dataset name."""
    return _require_display_text(value, label=label, allow_empty=False)


def _require_dataset_id(value: object, *, label: str = "dataset_id") -> str:
    """Require the portable producer identity accepted by every data path."""
    return _require_portable_filename_component(value, label=label)


def _require_native_boolean(value: object, *, label: str) -> bool:
    """Require one native boolean without truth-value coercion."""
    if type(value) is not bool:
        raise TypeError(f"{label} must be exactly True or False.")
    return value


def _require_positive_native_integer(value: object, *, label: str) -> int:
    """Require one positive native integer without numeric coercion."""
    if type(value) is not int:
        raise TypeError(f"{label} must be a native integer.")
    if value < 1:
        raise ValueError(f"{label} must be a positive integer.")
    return value


def _require_optional_native_string(
    value: object,
    *,
    label: str,
    allow_empty: bool,
) -> str | None:
    """Validate optional identity text without omission or string coercion."""
    if value is None:
        return None
    if type(value) is not str:
        raise TypeError(f"{label} must be None or a native string.")
    if not allow_empty and not value:
        raise ValueError(f"{label} must be None or a non-empty string.")
    return _require_display_text(value, label=label, allow_empty=True)


def _resolve_created_at(value: object) -> str:
    """Return one exact UTC-seconds timestamp for dataset identity metadata."""
    if value is None:
        return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    if type(value) is not str:
        raise TypeError(
            "created_at must be None or a native string in 'YYYY-MM-DDTHH:MM:SSZ' UTC format."
        )
    if not _CANONICAL_UTC_TIMESTAMP_PATTERN.fullmatch(value):
        raise ValueError("created_at must use exact 'YYYY-MM-DDTHH:MM:SSZ' UTC format.")
    try:
        datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ")
    except ValueError as error:
        raise ValueError(
            "created_at must be a valid UTC calendar timestamp in 'YYYY-MM-DDTHH:MM:SSZ' format."
        ) from error
    return value


def _protected_export_directories() -> frozenset[Path]:
    """Return every directory that must never be claimed as an export target.

    None of these is a directory anyone sets aside for one dataset: each is
    either a filesystem landmark or the place a user keeps unrelated work.
    """
    protected = {Path.cwd().resolve()}
    try:
        home = Path.home().resolve()
    except RuntimeError:
        return frozenset(protected)
    protected.add(home)
    protected.add(home.parent)
    return frozenset(protected)


def _require_export_output_directory(value: object, *, label: str = "out_dir") -> Path:
    """Require one dedicated export directory that is safe to replace whole.

    Publishing a generation renames the entire existing target aside and then
    removes it, so every file the named directory holds is destroyed by a
    replacement. A directory nobody set aside for one dataset -- the
    filesystem root, the working directory, the home directory, or the
    directory holding every home -- is refused here, where the argument is
    first resolved and before any path is created, locked, written, or
    removed.

    A leading ``~`` and every symbolic link are resolved before the
    comparison, both on the whole path and on the parent alone, so neither an
    alias to a protected directory nor a symbolic link standing in for one is
    accepted as a target.
    """
    if not isinstance(value, str | os.PathLike):
        raise TypeError(f"{label} must be a native string or os.PathLike path.")
    expanded = Path(value).expanduser()
    resolved = expanded.resolve()
    protected = _protected_export_directories()
    if (
        expanded.name in {"", ".."}
        or resolved.parent == resolved
        or resolved in protected
        or expanded.parent.resolve() / expanded.name in protected
    ):
        raise ValueError(
            f"{label} must name a dedicated dataset output directory, not {resolved}. "
            "Publishing replaces the whole directory, so the filesystem root, the "
            "current working directory, the home directory, and the directory holding "
            "it are refused. Pass a child directory of your own, such as "
            "'./exports/my_dataset'."
        )
    return expanded
