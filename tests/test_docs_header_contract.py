"""Exact contract for the shared documentation header."""

from __future__ import annotations

import ast
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPOSITORY_ROOT / "docs"


def _literal_config(name: str) -> object:
    tree = ast.parse((DOCS_ROOT / "conf.py").read_text(encoding="utf-8"))
    for statement in tree.body:
        if not isinstance(statement, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == name for target in statement.targets):
            return ast.literal_eval(statement.value)
    raise AssertionError(f"{name} is not assigned in docs/conf.py")


def test_docs_start_light_and_only_expose_the_web_app_header_link() -> None:
    theme_options = _literal_config("html_theme_options")
    html_context = _literal_config("html_context")
    assert isinstance(theme_options, dict)
    assert isinstance(html_context, dict)

    assert html_context["default_mode"] == "light"
    assert theme_options["navbar_end"] == ["theme-switcher", "navbar-icon-links"]
    assert theme_options["secondary_sidebar_items"] == ["page-toc"]
    assert "use_edit_page_button" not in theme_options
    assert theme_options["icon_links"] == [
        {
            "name": "Cellucid App",
            "url": "https://cellucid.com",
            "icon": "fa-solid fa-globe",
        }
    ]

    assert "switcher" not in theme_options
    assert "show_version_warning_banner" not in theme_options
    assert set(html_context) == {"default_mode"}
    assert not (DOCS_ROOT / "_static" / "switcher.json").exists()

    stylesheet = (DOCS_ROOT / "_static" / "custom.css").read_text(encoding="utf-8")
    assert "version-switcher" not in stylesheet
