from __future__ import annotations

from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
WEB_GUIDE = REPOSITORY_ROOT / "docs" / "user_guide" / "web_app"


def _read(relative_path: str) -> str:
    return (WEB_GUIDE / relative_path).read_text(encoding="utf-8")


def test_compact_help_behavior_is_documented_at_each_user_entry_point() -> None:
    glossary = _read("a_orientation/04_ui_glossary_terminology.md")
    accessibility = _read(
        "o_accessibility_privacy_security/01_accessibility.md"
    )
    data_loading = _read("b_data_loading/index.md")

    assert "Compact help (`i`) buttons" in glossary
    assert "does not replace the control's visible action label" in glossary
    assert "{kbd}`Enter`" in glossary
    assert "{kbd}`Space`" in glossary
    assert "{kbd}`Escape`" in glossary

    for exact_semantic in (
        'aria-controls',
        'aria-haspopup="dialog"',
        "aria-expanded",
        'role="dialog"',
        "{kbd}`Shift+Tab`",
        "Only one help dialog can be open.",
    ):
        assert exact_semantic in accessibility

    assert "keeps every loading action and connection state visible" in data_loading
    assert "compact `i` buttons only open explanatory help" in data_loading
