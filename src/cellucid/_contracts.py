"""
Contracts every user-supplied identity and label must satisfy.

These rules decide what Cellucid accepts as a name it will store or draw: a
portable dataset id, text the viewer can render exactly as stored, a native
scalar that was never coerced, and the directories nothing may claim whole.

The module depends on nothing but the standard library, and that is
load-bearing. ``prepare()``, the HTTP server, the notebook integration and the
web-asset cache all enforce these same rules, but only ``prepare()`` needs
numpy, pandas and scipy. Keeping the contracts here is what lets the light
modules import them without paying for the heavy ones.
"""

import re
import unicodedata
from pathlib import Path

_PORTABLE_COMPONENT_PATTERN = re.compile(
    r"^[A-Za-z0-9][A-Za-z0-9._-]{0,179}$",
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


def _protected_directories() -> frozenset[Path]:
    """Return every directory that must never be claimed whole as a target.

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
