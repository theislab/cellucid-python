"""The noun a count takes: ``plural(1, "cell")`` is ``"cell"``, ``plural(2,
"cell")`` is ``"cells"``.

Every message that writes a count before a noun writes the noun with this,
never as ``f"{n} cells"``, which reads "1 cells" the day ``n`` is one
(``tests/test_count_noun_agreement.py``). The count is written by the caller,
so its own formatting (``{n:,}``) stays as it was.
"""

from __future__ import annotations


def plural(count: int, singular: str, plural_form: str | None = None) -> str:
    """Return ``singular`` for a count of one, else its plural.

    ``plural_form`` is for a plural that is not ``singular + "s"``
    (``plural(n, "entry", "entries")``).
    """
    if not singular:
        raise ValueError("A counted noun must be non-empty text.")
    if count == 1:
        return singular
    return plural_form if plural_form is not None else f"{singular}s"
