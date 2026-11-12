"""
Cellucid's own category colors: the color the viewer gives category ``i`` of
a field that brings no colors of its own (no Scanpy ``uns['<key>_colors']``,
no ``palette``), as ``#rrggbb``.

The web app's ``getCategoryColor`` (``cellucid/assets/js/data/palettes.js``):
its 72 designed colors (the category palette, then the color picker's), then
hues fanned out around the color wheel by the golden angle, each channel the
byte it is drawn in (``round(channel × 255)``, ``utils/hex-color.js``
``hexOfRgb01``). cellucid-r keeps the same rule (``R/category-colors.R``);
the shared fixture ``category_default_colors.json`` (written by the web app's
``scripts/generate-category-default-colors.mjs``) pins the first 1,024
colors byte for byte in all three repositories.

An apply writes them as a user field's colors when the session gives none
(:mod:`cellucid.anndata_session._colors`): what the viewer showed.
"""

from __future__ import annotations

import math

#: The designed colors, as drawn (``CATEGORY_PALETTE`` then ``COLOR_PICKER_PALETTE``).
_SEEDS = (
    "#ff9900", "#0077c8", "#ad50e0", "#ffe500", "#00b783", "#e61977", "#4cafee", "#d15f00",
    "#009e9e", "#7d57c2", "#fa8032", "#5c6cd1", "#008b6e", "#dd2e5b", "#74ade1", "#ecb100",
    "#9339a7", "#3497b7", "#cb4b2e", "#7080c0", "#00a891", "#9294e3", "#c29411", "#347789",
    "#1f77b4", "#4169e1", "#00bfff", "#87ceeb", "#0000cd", "#191970", "#2ca02c", "#008000",
    "#32cd32", "#90ee90", "#006400", "#2e8b57", "#d62728", "#ff0000", "#ff6347", "#ff7f50",
    "#dc143c", "#b22222", "#ff7f0e", "#ffa500", "#ffd700", "#ffff00", "#bdb76b", "#ffe4b5",
    "#9467bd", "#800080", "#8a2be2", "#ee82ee", "#e377c2", "#ff69b4", "#8c564b", "#d2691e",
    "#a0522d", "#8b4513", "#deb887", "#f5deb3", "#7f7f7f", "#a9a9a9", "#d3d3d3", "#696969",
    "#2f4f4f", "#000000", "#17becf", "#00ced1", "#40e0d0", "#00ffff", "#48d1cc", "#5f9ea0",
)

_GOLDEN_ANGLE_DEGREES = 137.50776405003785


def _clamp(value: float, low: float, high: float) -> float:
    return min(high, max(low, value))


def _byte(channel: float) -> int:
    # JavaScript's Math.round: halves round up.
    return math.floor(channel * 255 + 0.5)


def default_category_color(index: int) -> str:
    """The color Cellucid gives category ``index`` (0-based) of a field that
    brings none, as ``#rrggbb``."""
    if isinstance(index, bool) or not isinstance(index, int) or index < 0:
        raise ValueError("A category's place is a non-negative whole number.")
    if index < len(_SEEDS):
        return _SEEDS[index]
    offset = index - len(_SEEDS)
    hue = ((offset * _GOLDEN_ANGLE_DEGREES) % 360) / 360
    saturation = _clamp(0.62 + 0.1 * math.sin(offset * 0.8), 0.45, 0.82)
    lightness = _clamp(0.54 + 0.12 * math.cos(offset * 0.6), 0.38, 0.72)
    chroma = (1 - abs(2 * lightness - 1)) * saturation
    sector = (hue % 1) * 6
    second = chroma * (1 - abs((sector % 2) - 1))
    if sector < 1:
        rgb = (chroma, second, 0.0)
    elif sector < 2:
        rgb = (second, chroma, 0.0)
    elif sector < 3:
        rgb = (0.0, chroma, second)
    elif sector < 4:
        rgb = (0.0, second, chroma)
    elif sector < 5:
        rgb = (second, 0.0, chroma)
    else:
        rgb = (chroma, 0.0, second)
    light = lightness - chroma / 2
    return "#" + "".join(f"{_byte(channel + light):02x}" for channel in rgb)
