"""
The data mapping: which stored arrays Cellucid plots and reads.

A *data mapping* names, for one dataset, the embeddings it can be shown with
(any number, of any mix of 1D, 2D and 3D: the *enabled* representations), the
default representation of each dimension, the vector fields bound to those
embeddings, and -- for an AnnData source -- the neighbour graph, the gene-name
set, the expression matrix and the latent space. The same JSON shape is written into session fingerprints,
sent over the Jupyter wire and derived from URLs, so it is defined once, here,
and mirrored by the web app's ``assets/js/data/data-mapping.js``.

Everything in this module is a pure rule over metadata, with two exceptions
that read an AnnData object without reading its large arrays:
:func:`build_anndata_candidates` (shapes and dtypes only) and
:func:`resolve_gene_names` (one ``var`` column). The rules that must agree
byte-for-byte between Python, R and JavaScript -- making names unique, ranking
best-guess embeddings, validating embedding columns, canonical embedding ids,
and the "nothing to plot" message -- are pinned by the shared conformance
fixture ``tests/fixtures/data_mapping_conformance.json``.

This module imports only the standard library at import time, because the
session bundle reader depends on it and must stay light.
"""

from __future__ import annotations

import copy
import math
import numbers
import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

#: The three dimension slots, in the order every mapping and identity lists them.
DIMENSION_SLOTS: tuple[str, ...] = ("1d", "2d", "3d")

#: How a caller names the embedding of one slot when serving AnnData: an
#: ``obsm`` key (read at columns ``0 .. d-1``) or a ``(key, columns)`` pair
#: with 0-based columns.
EmbeddingKeySpec = str | tuple[str, Sequence[int]]

#: The two kinds of data mapping. Exports fix the KNN graph, the expression
#: matrix and the latent space at ``prepare()`` time, so those groups exist
#: only on the AnnData kind.
EXPORT_KIND = "export"
ANNDATA_KIND = "anndata"

_EXPORT_MAPPING_KEYS = frozenset({"kind", "embeddings", "geneNames", "vectorFields"})
_ANNDATA_MAPPING_KEYS = frozenset(
    {"kind", "embeddings", "knn", "geneNames", "expression", "latent", "vectorFields"}
)
_EMBEDDING_REF_KEYS = frozenset({"key", "columns"})
_EMBEDDINGS_GROUP_KEYS = frozenset({"enabled", "defaults"})
_EXPORT_VECTOR_BINDING_KEYS = frozenset({"field", "embedding"})
_ANNDATA_VECTOR_BINDING_KEYS = frozenset({"key", "embedding"})

#: ``sc.tl.umap()`` writes ``X_umap``; earlier Cellucid releases asked for one
#: key per dimension. Both are the default rule (design §1.4).
DIMENSIONAL_UMAP_KEYS: tuple[str, ...] = ("X_umap_1d", "X_umap_2d", "X_umap_3d")
UNSUFFIXED_UMAP_KEY = "X_umap"

#: The AnnData ``obsm`` grammar of a per-dimension UMAP vector field. The three
#: dimensional embedding keys match it too (``X_umap_2d`` reads as field
#: ``X_umap``), so test keys with :func:`is_vector_field_key`, which excludes them.
VECTOR_KEY_PATTERN = re.compile(r"^(?P<field>.+_umap)_(?P<dimension>[123])d$")


def is_vector_field_key(key: object) -> bool:
    """Whether an ``obsm`` key declares a UMAP vector field.

    ``X_umap_1d``, ``X_umap_2d`` and ``X_umap_3d`` are embeddings, not vectors,
    although they fit the vector grammar.
    """
    return (
        isinstance(key, str)
        and key not in DIMENSIONAL_UMAP_KEYS
        and VECTOR_KEY_PATTERN.fullmatch(key) is not None
    )

_GENE_NAMES_COLUMN_PREFIX = "column:"
_LAYER_PREFIX = "layer:"

#: Best-guess tiers, in the order a name is tested against them. A name is
#: lowercased, and the first tier whose substring it contains decides its rank;
#: a name that matches none ranks 7. Tier 8 (latent spaces) is tested last so
#: that, for example, ``X_umap_on_pca`` is still a UMAP. The name keeps its
#: leading ``x_``: every other needle matches the same names either way, and
#: stripping it would leave the ``x_emb`` needle unable to match ``X_emb``.
_BEST_GUESS_TIERS: tuple[tuple[int, tuple[str, ...]], ...] = (
    (0, ("umap",)),
    (1, ("tsne", "t-sne", "t_sne")),
    (2, ("draw_graph", "force")),
    (3, ("phate",)),
    (4, ("pacmap", "trimap")),
    (5, ("spatial",)),
    (6, ("diffmap",)),
    (8, ("pca", "scvi", "latent", "harmony", "x_emb")),
)
_BEST_GUESS_NO_MATCH_TIER = 7


class DataMappingError(ValueError):
    """A mapping, column list or name set that breaks one exact rule.

    ``code`` is the stable machine-readable reason the conformance fixture
    pins; the message is for people.
    """

    def __init__(self, code: str, message: str) -> None:
        super().__init__(message)
        self.code = code


# =============================================================================
# Making names unique (design §1.5)
# =============================================================================


@dataclass(frozen=True)
class UniqueNames:
    """The result of :func:`make_unique_names`.

    ``renamed`` lists ``(from, to)`` for every renamed position, in input
    order; ``filled_from_index`` counts the empty values replaced by the
    position's index name before uniqueness was applied. ``filled`` holds each
    position's name after that replacement and before uniqueness, and
    ``from_index`` says which positions were replaced, so a caller can report
    on a subset of the positions.
    """

    names: list[str]
    renamed: list[tuple[str, str]]
    filled_from_index: int
    filled: list[str]
    from_index: list[bool]


def _is_missing_name(value: object) -> bool:
    """Whether one name is empty: ``None``, ``""``, a NaN, or a pandas missing value."""
    if value is None:
        return True
    if isinstance(value, str):
        return value == ""
    if isinstance(value, float):
        return math.isnan(value)
    # pandas.NA and pandas.NaT, detected without importing pandas.
    return type(value).__name__ in {"NAType", "NaTType"}


def gene_name_text(value: object) -> str | None:
    """The text of one gene-name value from a ``var`` column, or None when missing.

    The one rule every path applies (the web app's ``geneNameText``): a string
    is itself; ``None``, NaN and pandas missing values are missing; a boolean is
    ``True``/``False``; an integer is its decimal digits; any other real number
    is ``str(float)``. Anything else cannot name a gene.
    """
    if _is_missing_name(value):
        return None
    if isinstance(value, str):
        return value
    if isinstance(value, bool) or type(value).__name__ == "bool_":
        return "True" if value else "False"
    if isinstance(value, numbers.Integral):
        return str(int(value))
    if isinstance(value, numbers.Real):
        number = float(value)
        return None if math.isnan(number) else str(number)
    raise DataMappingError(
        "name_not_string",
        f"A gene name must be text or a number, got {type(value).__name__}.",
    )


def column_gene_names(var: Any, column: str | None) -> UniqueNames:
    """The unique gene names of one ``var`` frame: its index, or one column.

    A column must be one the candidate catalog lists as usable (text,
    categories or integers); its values are turned into text with
    :func:`gene_name_text`, an empty one takes the row's index name (which
    must be text), and the list is made unique with
    :func:`make_unique_names`. Every path -- ``prepare()``, the served gene
    list, the mapping routes and the session reader -- names genes with this
    one function.
    """
    index = var.index.tolist()
    if column is None:
        return make_unique_names(index)
    if column not in var.columns:
        raise KeyError(f"var column {column!r} is absent. Available columns: {list(var.columns)}")
    verdict = next(item for item in _describe_var_columns(var) if item["key"] == str(column))
    if not verdict["usable"]:
        raise DataMappingError("gene_names_column", f"var column {column!r}: {verdict['reason']}.")
    return make_unique_names([gene_name_text(value) for value in var[column].tolist()], index=index)


def make_unique_names(
    names: Sequence[object],
    *,
    index: Sequence[object] | None = None,
    join: str = "-",
) -> UniqueNames:
    """Make one name list unique with anndata's ``make_index_unique`` rule.

    1. An empty name (``None``, ``""`` or NaN) is replaced by the index value
       at the same position first. ``index`` is then required, and the value
       it supplies must itself be a non-empty string.
    2. The list is walked in order. The first occurrence of a name keeps it.
    3. Each later occurrence of ``name`` becomes ``name{join}k``. ``k`` is a
       counter per original name that only counts up, and a candidate is
       skipped while it already exists in the filled list or has already been
       produced.

    A list that is already unique comes back unchanged, so a unique gene axis
    exports exactly the names it exported before this rule existed.
    """
    if not isinstance(join, str) or not join:
        raise TypeError("join must be a non-empty string.")
    if index is not None and len(index) != len(names):
        raise ValueError(
            f"index holds {len(index)} names but the name list holds {len(names)}."
        )

    filled: list[str] = []
    from_index: list[bool] = []
    filled_from_index = 0
    for position, value in enumerate(names):
        from_index.append(_is_missing_name(value))
        if from_index[-1]:
            if index is None:
                raise DataMappingError(
                    "missing_name",
                    f"Name at position {position} is empty and no index was supplied "
                    "to take its name from.",
                )
            index_name = index[position]
            if not isinstance(index_name, str) or not index_name:
                raise DataMappingError(
                    "missing_index_name",
                    f"Name at position {position} is empty, and its index name is not "
                    f"a non-empty string: {index_name!r}.",
                )
            filled.append(index_name)
            filled_from_index += 1
            continue
        if not isinstance(value, str):
            raise DataMappingError(
                "name_not_string",
                f"Name at position {position} must be a string, got {type(value).__name__}.",
            )
        filled.append(value)

    taken = set(filled)
    seen: set[str] = set()
    counters: dict[str, int] = {}
    result: list[str] = []
    renamed: list[tuple[str, str]] = []
    for value in filled:
        if value not in seen:
            seen.add(value)
            result.append(value)
            continue
        while True:
            counters[value] = counters.get(value, 0) + 1
            candidate = f"{value}{join}{counters[value]}"
            if candidate not in taken:
                break
        taken.add(candidate)
        result.append(candidate)
        renamed.append((value, candidate))
    return UniqueNames(
        names=result,
        renamed=renamed,
        filled_from_index=filled_from_index,
        filled=filled,
        from_index=from_index,
    )


def _plural(count: int, singular: str, plural: str) -> str:
    return singular if count == 1 else plural


def unique_names_note(
    *,
    filled_from_index: int,
    renamed: Sequence[Sequence[str]],
    source: str | None,
) -> str | None:
    """Return the one-line note that reports what making names unique changed.

    ``source`` is the ``var`` column the names came from, or None for the var
    index (which cannot take names from itself). The text is shared with the web
    app and pinned by the conformance fixture, so it quotes names with plain
    single quotes and never escapes them.
    """
    parts: list[str] = []
    if filled_from_index:
        if source is None:
            raise ValueError("The var index cannot fill names from itself.")
        parts.append(
            f"{filled_from_index:,} "
            + _plural(
                filled_from_index,
                f"gene has no '{source}' value and uses its var index name.",
                f"genes have no '{source}' value and use their var index name.",
            )
        )
    if renamed:
        first_from, first_to = renamed[0][0], renamed[0][1]
        parts.append(
            f"{len(renamed):,} "
            + _plural(
                len(renamed),
                "duplicate name was made unique",
                "duplicate names were made unique",
            )
            + f", e.g. '{first_from}' → '{first_to}'."
        )
    return " ".join(parts) if parts else None


# =============================================================================
# Embedding references, columns and canonical ids (design §1.1, §1.3)
# =============================================================================


def slot_dimension(slot: str) -> int:
    """Return the dimension a ``"1d"``/``"2d"``/``"3d"`` slot names."""
    if slot not in DIMENSION_SLOTS:
        raise ValueError(f"Embedding slot must be one of {list(DIMENSION_SLOTS)}, got {slot!r}.")
    return int(slot[0])


def dimension_slot(dimension: int) -> str:
    """Return the slot name of one dimension."""
    if type(dimension) is not int or dimension not in (1, 2, 3):
        raise ValueError(f"Dimension must be exactly 1, 2, or 3, got {dimension!r}.")
    return f"{dimension}d"


def default_columns(dimension: int) -> list[int]:
    """The columns a dimension reads when none are named: ``0 .. d-1``."""
    return list(range(dimension))


def validate_embedding_columns(
    columns: object,
    *,
    dimension: int,
    width: int,
) -> list[int]:
    """Validate the columns one embedding reads for one dimension.

    The checks run in a fixed order and the first failure is reported, so every
    implementation reports the same reason for the same input:

    ``too_narrow``
        the array has fewer than ``dimension`` columns;
    ``columns_not_list``
        ``columns`` is not a list;
    ``column_count``
        the list does not hold exactly ``dimension`` entries;
    ``column_not_integer``
        an entry is not an integer (booleans are not integers);
    ``column_out_of_range``
        an entry is negative or not below ``width``;
    ``column_duplicate``
        two entries name the same column.
    """
    if type(dimension) is not int or dimension not in (1, 2, 3):
        raise ValueError(f"dimension must be exactly 1, 2, or 3, got {dimension!r}.")
    if type(width) is not int or width < 0:
        raise ValueError(f"width must be a non-negative integer, got {width!r}.")
    if width < dimension:
        raise DataMappingError(
            "too_narrow",
            f"A {dimension}D embedding needs at least {dimension} columns; "
            f"this array has {width}.",
        )
    if not isinstance(columns, list | tuple):
        raise DataMappingError("columns_not_list", "Embedding columns must be a list.")
    if len(columns) != dimension:
        raise DataMappingError(
            "column_count",
            f"A {dimension}D embedding reads exactly {dimension} columns; got {len(columns)}.",
        )
    resolved: list[int] = []
    for column in columns:
        if type(column) is not int:
            raise DataMappingError(
                "column_not_integer",
                f"Embedding column {column!r} must be an integer.",
            )
        if not 0 <= column < width:
            raise DataMappingError(
                "column_out_of_range",
                f"Embedding column {column} does not exist; the array has columns "
                f"0 to {width - 1}.",
            )
        resolved.append(column)
    if len(set(resolved)) != len(resolved):
        raise DataMappingError("column_duplicate", "Each axis needs a different column.")
    return resolved


def canonical_embedding_id(key: str, columns: Sequence[int] | None, *, width: int | None) -> str:
    """Return the id one embedding reference is known by.

    An export item is known by its own id (``columns`` is None). An AnnData
    ``obsm`` array is known by its key when it is read whole -- its width equals
    the number of columns and they are ``0 .. d-1`` in order -- and otherwise by
    ``key[c0,c1]``. Because ``(key, columns)`` is unique, so is the id.
    """
    if not isinstance(key, str) or not key:
        raise ValueError("Embedding key must be a non-empty string.")
    if columns is None:
        return key
    resolved = list(columns)
    if width is not None and width == len(resolved) and resolved == default_columns(len(resolved)):
        return key
    return f"{key}[{','.join(str(column) for column in resolved)}]"


# =============================================================================
# Validation of the DataMapping JSON (design §1.1)
# =============================================================================


def _require_exact_keys(value: object, expected: frozenset[str], *, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise DataMappingError("not_object", f"{label} must be an object.")
    missing = sorted(expected - set(value))
    unknown = sorted(str(key) for key in set(value) - expected)
    if missing or unknown:
        details: list[str] = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if unknown:
            details.append("unknown " + ", ".join(unknown))
        raise DataMappingError("keys", f"{label} fields: " + "; ".join(details) + ".")
    return value


def _require_mapping_string(value: object, *, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise DataMappingError("not_string", f"{label} must be a non-empty string.")
    return value


def validate_embedding_ref(
    value: object,
    *,
    kind: str,
    dimension: int | None,
    label: str,
) -> dict[str, Any]:
    """Validate one ``{key, columns}`` reference and return a canonical copy.

    An export reference names an exported embedding id and carries
    ``columns: null``, because the columns were fixed at export. An AnnData
    reference names an ``obsm`` key and always lists its 0-based columns.
    """
    ref = _require_exact_keys(value, _EMBEDDING_REF_KEYS, label=label)
    key = _require_mapping_string(ref["key"], label=f"{label}.key")
    columns = ref["columns"]
    if kind == EXPORT_KIND:
        if columns is not None:
            raise DataMappingError(
                "columns_not_null",
                f"{label}.columns must be null: an exported embedding's columns were "
                "fixed when it was prepared.",
            )
        return {"key": key, "columns": None}
    if kind != ANNDATA_KIND:
        raise ValueError(f"Unknown data mapping kind {kind!r}.")
    if not isinstance(columns, list):
        raise DataMappingError("columns_not_list", f"{label}.columns must be a list.")
    if dimension is not None and len(columns) != dimension:
        raise DataMappingError(
            "column_count",
            f"{label}.columns must list exactly {dimension} columns.",
        )
    if not 1 <= len(columns) <= 3:
        raise DataMappingError("column_count", f"{label}.columns must list 1 to 3 columns.")
    for column in columns:
        if type(column) is not int or column < 0:
            raise DataMappingError(
                "column_not_integer",
                f"{label}.columns must hold non-negative integers.",
            )
    if len(set(columns)) != len(columns):
        raise DataMappingError("column_duplicate", f"{label}: each axis needs a different column.")
    return {"key": key, "columns": list(columns)}


def _validate_gene_names(value: object, *, label: str) -> str:
    if value == "index":
        return "index"
    if (
        isinstance(value, str)
        and value.startswith(_GENE_NAMES_COLUMN_PREFIX)
        and len(value) > len(_GENE_NAMES_COLUMN_PREFIX)
    ):
        return value
    raise DataMappingError(
        "gene_names",
        f"{label} must be \"index\" or \"column:<var column>\".",
    )


def _validate_expression(value: object, *, label: str) -> str | None:
    if value is None:
        return None
    if isinstance(value, str) and (
        value in ("X", "raw")
        or (value.startswith(_LAYER_PREFIX) and len(value) > len(_LAYER_PREFIX))
    ):
        return value
    raise DataMappingError(
        "expression",
        f"{label} must be \"X\", \"raw\", \"layer:<name>\" or null.",
    )


def _optional_mapping_string(value: object, *, label: str) -> str | None:
    if value is None:
        return None
    return _require_mapping_string(value, label=label)


def embedding_ref_sort_key(ref: Mapping[str, Any]) -> tuple[Any, ...]:
    """The canonical order of embedding references.

    By ``key`` in code-point order, then ``columns``: an export's ``null``
    first, then lists compared element by element.
    """
    columns = ref["columns"]
    return (ref["key"], 0 if columns is None else 1, tuple(columns or ()))


def _validate_embeddings_group(value: object, *, kind: str, label: str) -> dict[str, Any]:
    """Validate the ``embeddings`` group: ``{"enabled": [...], "defaults": {...}}``.

    ``enabled`` lists every representation the dataset can be shown with: at
    least one, in canonical order (:func:`embedding_ref_sort_key`), without
    repeats. ``defaults`` names at most one enabled representation per
    dimension slot, and at least one in all.
    """
    group = _require_exact_keys(value, _EMBEDDINGS_GROUP_KEYS, label=label)
    raw_enabled = group["enabled"]
    if not isinstance(raw_enabled, list):
        raise DataMappingError("not_list", f"{label}.enabled must be a list.")
    if not raw_enabled:
        raise DataMappingError("no_embedding", f"{label}.enabled must list at least one embedding.")
    enabled = [
        validate_embedding_ref(ref, kind=kind, dimension=None, label=f"{label}.enabled[{position}]")
        for position, ref in enumerate(raw_enabled)
    ]
    order = [embedding_ref_sort_key(ref) for ref in enabled]
    if order != sorted(order):
        raise DataMappingError(
            "embedding_order",
            f"{label}.enabled must be sorted by embedding key, then columns.",
        )
    if len(set(order)) != len(order):
        raise DataMappingError("embedding_duplicate", f"{label}.enabled repeats an embedding.")
    slots = _require_exact_keys(group["defaults"], frozenset(DIMENSION_SLOTS), label=f"{label}.defaults")
    defaults: dict[str, dict[str, Any] | None] = {}
    enabled_keys = set(order)
    for slot in DIMENSION_SLOTS:
        ref = slots[slot]
        if ref is None:
            defaults[slot] = None
            continue
        valid = validate_embedding_ref(
            ref,
            kind=kind,
            dimension=slot_dimension(slot),
            label=f"{label}.defaults.{slot}",
        )
        if embedding_ref_sort_key(valid) not in enabled_keys:
            raise DataMappingError(
                "default_not_enabled",
                f"{label}.defaults.{slot} must be one of the enabled embeddings.",
            )
        defaults[slot] = valid
    if all(ref is None for ref in defaults.values()):
        raise DataMappingError("no_embedding", f"{label}.defaults must name at least one default.")
    return {"enabled": enabled, "defaults": defaults}


def _binding_sort_key(binding: Mapping[str, Any], name_key: str) -> tuple[Any, ...]:
    embedding = binding["embedding"]
    columns = embedding["columns"]
    return (
        embedding["key"],
        0 if columns is None else 1,
        tuple(columns or ()),
        binding[name_key],
    )


def validate_data_mapping(value: object, *, label: str = "dataMapping") -> dict[str, Any]:
    """Validate one DataMapping and return a canonical deep copy.

    The mapping is a closed shape per ``kind``: an unknown or missing key is a
    mapping this version does not understand. At least one embedding must be
    enabled and at least one default named, since a mapping that plots nothing
    describes no session; every vector binding is bound to an enabled
    embedding.

    ``vectorFields`` must already be in canonical order -- by the bound
    embedding's ``key``, then its ``columns`` (an export's ``null`` first), then
    the field id or ``obsm`` key -- and must not repeat a binding.
    """
    if not isinstance(value, dict):
        raise DataMappingError("not_object", f"{label} must be an object.")
    kind = value.get("kind")
    if kind == EXPORT_KIND:
        expected = _EXPORT_MAPPING_KEYS
        binding_keys = _EXPORT_VECTOR_BINDING_KEYS
        name_key = "field"
    elif kind == ANNDATA_KIND:
        expected = _ANNDATA_MAPPING_KEYS
        binding_keys = _ANNDATA_VECTOR_BINDING_KEYS
        name_key = "key"
    else:
        raise DataMappingError("kind", f"{label}.kind must be \"export\" or \"anndata\".")
    mapping = _require_exact_keys(value, expected, label=label)

    embeddings = _validate_embeddings_group(mapping["embeddings"], kind=kind, label=f"{label}.embeddings")

    raw_bindings = mapping["vectorFields"]
    if not isinstance(raw_bindings, list):
        raise DataMappingError("not_list", f"{label}.vectorFields must be a list.")
    bindings: list[dict[str, Any]] = []
    for position, raw in enumerate(raw_bindings):
        binding_label = f"{label}.vectorFields[{position}]"
        binding = _require_exact_keys(raw, binding_keys, label=binding_label)
        bindings.append(
            {
                name_key: _require_mapping_string(
                    binding[name_key],
                    label=f"{binding_label}.{name_key}",
                ),
                "embedding": validate_embedding_ref(
                    binding["embedding"],
                    kind=kind,
                    dimension=None,
                    label=f"{binding_label}.embedding",
                ),
            }
        )
    sort_keys = [_binding_sort_key(binding, name_key) for binding in bindings]
    if sort_keys != sorted(sort_keys):
        raise DataMappingError(
            "vector_order",
            f"{label}.vectorFields must be sorted by embedding key, embedding columns, "
            f"then {name_key}.",
        )
    if len(set(sort_keys)) != len(sort_keys):
        raise DataMappingError("vector_duplicate", f"{label}.vectorFields repeats a binding.")
    enabled_keys = {embedding_ref_sort_key(ref) for ref in embeddings["enabled"]}
    for position, binding in enumerate(bindings):
        if embedding_ref_sort_key(binding["embedding"]) not in enabled_keys:
            raise DataMappingError(
                "vector_not_enabled",
                f"{label}.vectorFields[{position}] is bound to an embedding that is not enabled.",
            )

    canonical: dict[str, Any] = {"kind": kind, "embeddings": embeddings}
    if kind == ANNDATA_KIND:
        canonical["knn"] = _optional_mapping_string(mapping["knn"], label=f"{label}.knn")
    canonical["geneNames"] = _validate_gene_names(mapping["geneNames"], label=f"{label}.geneNames")
    if kind == ANNDATA_KIND:
        canonical["expression"] = _validate_expression(
            mapping["expression"],
            label=f"{label}.expression",
        )
        canonical["latent"] = _optional_mapping_string(
            mapping["latent"],
            label=f"{label}.latent",
        )
    canonical["vectorFields"] = bindings
    return copy.deepcopy(canonical)


def export_data_mapping(identity: Mapping[str, Any], *, gene_name_source: str | None) -> dict[str, Any]:
    """The data mapping a prepared export is shown with when it opens.

    Every declared embedding is enabled, in canonical order; the defaults are
    the identity's ``embeddings.defaults``; genes are named by the default
    gene-name set (``gene_name_source`` is its ``var`` column, None for the
    index); every declared vector-field file is bound to its embedding. The
    same rule holds for exports written by the Python and the R package, which
    declare embeddings identically. Mirrors the web app's
    ``resolveExportDataMapping``.
    """
    embeddings = identity["embeddings"]
    defaults = embeddings["defaults"]
    bindings = [
        {"field": field, "embedding": {"key": embedding_id, "columns": None}}
        for field, entry in identity.get("vector_fields", {}).get("fields", {}).items()
        for embedding_id in entry["files"]
    ]
    bindings.sort(key=lambda binding: _binding_sort_key(binding, "field"))
    return validate_data_mapping(
        {
            "kind": EXPORT_KIND,
            "embeddings": {
                "enabled": [
                    {"key": embedding_id, "columns": None}
                    for embedding_id in sorted(embeddings["items"])
                ],
                "defaults": {
                    slot: ({"key": defaults[slot], "columns": None} if slot in defaults else None)
                    for slot in DIMENSION_SLOTS
                },
            },
            "geneNames": "index" if gene_name_source is None else f"{_GENE_NAMES_COLUMN_PREFIX}{gene_name_source}",
            "vectorFields": bindings,
        },
        label="Export data mapping",
    )


def gene_names_column(gene_names: str) -> str | None:
    """Return the ``var`` column a ``geneNames`` value names, or None for the index."""
    if gene_names == "index":
        return None
    if gene_names.startswith(_GENE_NAMES_COLUMN_PREFIX):
        return gene_names[len(_GENE_NAMES_COLUMN_PREFIX) :]
    raise ValueError(f"Invalid geneNames value {gene_names!r}.")


def expression_layer(expression: str) -> str | None:
    """Return the layer an ``expression`` value names, or None when it is not a layer."""
    if expression.startswith(_LAYER_PREFIX):
        return expression[len(_LAYER_PREFIX) :]
    return None


# =============================================================================
# Default resolution and the deterministic best guess (design §1.4)
# =============================================================================


@dataclass(frozen=True)
class ObsmEntry:
    """What the embedding rules need to know about one ``obsm`` entry.

    ``width`` is the column count of a dense array (1 for a plain vector) and
    None for anything that is not a dense array. ``usable`` is the listing-time
    verdict of design §1.3 -- dense, numeric, one row per cell, at least one
    column -- and ``reason`` says why not.
    """

    key: str
    width: int | None
    usable: bool
    reason: str | None = None


def default_embedding_slots(entries: Sequence[ObsmEntry]) -> dict[str, dict[str, Any] | None]:
    """Apply today's default embedding rule and return the filled slots.

    If any of ``X_umap_1d``, ``X_umap_2d`` or ``X_umap_3d`` is present, those
    keys alone fill their slots, reading columns ``0 .. d-1``. Otherwise a bare
    ``X_umap`` whose width is 1, 2 or 3 fills the slot of that width. Slots
    without a default stay None. Presence is all this rule asks: a present key
    that later fails its load-time check is reported, not skipped, because the
    caller named that array on purpose.
    """
    by_key = {entry.key: entry for entry in entries}
    slots: dict[str, dict[str, Any] | None] = dict.fromkeys(DIMENSION_SLOTS)
    dimensional = [key for key in DIMENSIONAL_UMAP_KEYS if key in by_key]
    if dimensional:
        for key in dimensional:
            dimension = int(key[len("X_umap_")])
            slots[dimension_slot(dimension)] = {"key": key, "columns": default_columns(dimension)}
        return slots
    bare = by_key.get(UNSUFFIXED_UMAP_KEY)
    if bare is not None and bare.width in (1, 2, 3):
        dimension = bare.width
        slots[dimension_slot(dimension)] = {
            "key": UNSUFFIXED_UMAP_KEY,
            "columns": default_columns(dimension),
        }
    return slots


def best_guess_tier(key: str) -> int:
    """Rank one ``obsm`` key for the best guess; a lower tier is a better guess."""
    name = key.lower()
    for tier, needles in _BEST_GUESS_TIERS:
        if any(needle in name for needle in needles):
            return tier
    return _BEST_GUESS_NO_MATCH_TIER


def best_guess_order(keys: Iterable[str]) -> list[str]:
    """Order keys by (tier, key in code-point order)."""
    return sorted(keys, key=lambda key: (best_guess_tier(key), key))


@dataclass(frozen=True)
class EmbeddingGuess:
    """The outcome of :func:`best_guess_embedding_slots`."""

    slots: dict[str, dict[str, Any] | None]
    guessed: list[str]
    rejected: list[tuple[str, str]]


LoadCheck = Callable[[str, list[int]], str | None]


def best_guess_embedding_slots(
    entries: Sequence[ObsmEntry],
    *,
    latent_key: str | None,
    load_check: LoadCheck | None = None,
) -> EmbeddingGuess:
    """Fill embedding slots when the default rule filled none.

    Candidates are the usable entries, minus vector-field keys
    (``<field>_umap_<d>d``) and the latent key. For each of 3D, 2D and 1D the
    slot takes the best-ranked *narrow* candidate (width 1-3) whose width is
    exactly that dimension. When no narrow candidate fills any slot, only the
    2D slot is filled, from the best-ranked *wide* candidate (width > 3) read
    at columns ``[0, 1]``.

    ``load_check(key, columns)`` returns None when the array loads, or the
    reason it does not. A candidate that fails is recorded in ``rejected`` and
    the next one is tried, so the guess is deterministic and bounded by the
    candidate count.
    """
    candidates = [
        entry
        for entry in entries
        if entry.usable
        and entry.width is not None
        and entry.width >= 1
        and not is_vector_field_key(entry.key)
        and entry.key != latent_key
    ]
    width_by_key = {entry.key: entry.width for entry in candidates}
    rejected: list[tuple[str, str]] = []
    failed: set[str] = set()
    slots: dict[str, dict[str, Any] | None] = dict.fromkeys(DIMENSION_SLOTS)

    def try_candidate(key: str, columns: list[int]) -> bool:
        if key in failed:
            return False
        reason = load_check(key, columns) if load_check is not None else None
        if reason is None:
            return True
        failed.add(key)
        rejected.append((key, reason))
        return False

    for dimension in (3, 2, 1):
        for key in best_guess_order(k for k, w in width_by_key.items() if w == dimension):
            columns = default_columns(dimension)
            if try_candidate(key, columns):
                slots[dimension_slot(dimension)] = {"key": key, "columns": columns}
                break

    if all(ref is None for ref in slots.values()):
        for key in best_guess_order(k for k, w in width_by_key.items() if (w or 0) > 3):
            columns = [0, 1]
            if try_candidate(key, columns):
                slots["2d"] = {"key": key, "columns": columns}
                break

    guessed = [f"embeddings.{slot}" for slot in DIMENSION_SLOTS if slots[slot] is not None]
    return EmbeddingGuess(slots=slots, guessed=guessed, rejected=rejected)


def resolve_embedding_slots(
    entries: Sequence[ObsmEntry],
    *,
    latent_key: str | None,
    load_check: LoadCheck | None = None,
) -> EmbeddingGuess:
    """Resolve the initial embedding slots: the default rule, else the best guess.

    When the default rule fills any slot, its slots are the answer and nothing
    is guessed -- so an object carrying ``X_umap`` or ``X_umap_<d>d`` resolves
    exactly as it did before the best guess existed. Only when it fills none
    does :func:`best_guess_embedding_slots` run.
    """
    slots = default_embedding_slots(entries)
    if any(ref is not None for ref in slots.values()):
        return EmbeddingGuess(slots=slots, guessed=[], rejected=[])
    return best_guess_embedding_slots(entries, latent_key=latent_key, load_check=load_check)


def default_enabled_embeddings(
    entries: Sequence[ObsmEntry],
    slots: Mapping[str, Mapping[str, Any] | None],
    *,
    latent_key: str | None,
    rejected: Iterable[str] = (),
) -> list[dict[str, Any]]:
    """The representations an AnnData object is shown with when it opens.

    Every filled slot, plus every other usable ``obsm`` array that is 1 to 3
    columns wide, read at all of its columns -- minus vector-field keys
    (``<field>_umap_<d>d``), the latent key and the ``rejected`` keys (arrays
    the best guess found unloadable). A wider array is enabled only when it
    fills a slot; the user can enable other column selections. Enabling reads
    nothing: an array is read only when a view shows it. The result is in
    canonical order (:func:`embedding_ref_sort_key`).
    """
    skipped = set(rejected)
    refs: dict[tuple[Any, ...], dict[str, Any]] = {}
    for ref in slots.values():
        if ref is not None:
            copy_ref = {"key": ref["key"], "columns": list(ref["columns"])}
            refs[embedding_ref_sort_key(copy_ref)] = copy_ref
    for entry in entries:
        if (
            entry.usable
            and entry.width in (1, 2, 3)
            and not is_vector_field_key(entry.key)
            and entry.key != latent_key
            and entry.key not in skipped
        ):
            assert entry.width is not None
            ref = {"key": entry.key, "columns": default_columns(entry.width)}
            refs.setdefault(embedding_ref_sort_key(ref), ref)
    return [refs[order] for order in sorted(refs)]


def default_vector_bindings(
    keys: Iterable[str],
    slots: Mapping[str, Mapping[str, Any] | None],
) -> list[dict[str, Any]]:
    """Bind ``<field>_umap_<d>d`` keys by today's rule.

    A key binds to the slot-``d`` embedding only when that slot holds
    ``X_umap_<d>d`` or ``X_umap``: a UMAP vector drawn over another layout would
    point the wrong way. The result is in canonical binding order.
    """
    bindings: list[dict[str, Any]] = []
    for key in keys:
        if not is_vector_field_key(key):
            continue
        match = VECTOR_KEY_PATTERN.fullmatch(key)
        assert match is not None
        slot = dimension_slot(int(match.group("dimension")))
        ref = slots.get(slot)
        if ref is None or ref["key"] not in (f"X_umap_{slot}", UNSUFFIXED_UMAP_KEY):
            continue
        bindings.append({"key": key, "embedding": {"key": ref["key"], "columns": list(ref["columns"])}})
    bindings.sort(key=lambda binding: _binding_sort_key(binding, "key"))
    return bindings


def default_expression(
    *,
    has_x: bool,
    has_raw: bool,
    layers: Sequence[str],
) -> tuple[str | None, bool]:
    """Return the default expression source and whether it was guessed.

    ``X`` when it exists. Otherwise ``raw``; otherwise the only layer;
    otherwise the first layer whose name contains ``norm`` or ``log``, then the
    remaining layers in code-point order, with ``counts`` and ``raw_counts``
    last. None when there is nothing to read.
    """
    if has_x:
        return "X", False
    if has_raw:
        return "raw", True
    if not layers:
        return None, False
    if len(layers) == 1:
        return f"{_LAYER_PREFIX}{layers[0]}", True

    def rank(name: str) -> tuple[int, str]:
        lowered = name.lower()
        if "norm" in lowered or "log" in lowered:
            return (0, name)
        if name in ("counts", "raw_counts"):
            return (2, name)
        return (1, name)

    return f"{_LAYER_PREFIX}{min(layers, key=rank)}", True


# =============================================================================
# The "nothing to plot" message (design §7)
# =============================================================================


def describe_obsm_entry(
    *,
    key: str,
    storage: str,
    shape: Sequence[int] | None,
    defect: str | None,
) -> str:
    """Describe one ``obsm`` entry the way the "nothing to plot" message lists it."""
    if storage == "dataframe":
        text = "table"
    else:
        extent = " × ".join(f"{int(size):,}" for size in shape) if shape else "?"
        text = f"sparse {extent}" if storage in ("csr", "csc") else extent
    if defect:
        text = f"{text}, {defect}"
    return f"{key} ({text})"


def obsm_entry_defect(
    item: Mapping[str, Any],
    *,
    latent_key: str | None,
    rejected: Mapping[str, str],
) -> str | None:
    """Say why one described ``obsm`` entry was not plotted, for the message.

    A table or a sparse array already reads as what it is, so only a dense
    array carries a reason: its listing-time defect, the load-time defect that
    rejected it, or the role that kept it out of the candidates.
    """
    key = item["key"]
    if key in rejected:
        return rejected[key]
    if item["storage"] != "dense":
        return None
    if not item["usable"]:
        return str(item["reason"])
    if is_vector_field_key(key):
        return "a vector field"
    if key == latent_key:
        return "the latent space"
    return None


def nothing_plottable_message(*, n_obs: int, descriptions: Sequence[str]) -> str:
    """Return the message shown when no ``obsm`` array can be plotted."""
    holds = ", ".join(descriptions) if descriptions else "nothing"
    return (
        "Cellucid found nothing to plot. It needs an obsm array with one row per "
        f"cell ({n_obs:,}) and at least one numeric column. obsm holds: {holds}."
    )


# =============================================================================
# The AnnData candidate catalog (design §1.2) -- metadata only
# =============================================================================


def _storage_of(value: Any) -> str:
    """Classify one AnnData element without reading its values."""
    module = type(value).__module__
    name = type(value).__name__
    if name == "DataFrame" and module.startswith("pandas"):
        return "dataframe"
    fmt = getattr(value, "format", None)
    if module.startswith("scipy.sparse") or name.startswith(("CSR", "CSC", "csr_", "csc_")):
        if fmt == "csr":
            return "csr"
        if fmt == "csc":
            return "csc"
        return "other"
    if module.startswith("anndata") and fmt in ("csr", "csc"):
        # Backed sparse datasets (anndata's BaseCompressedSparseDataset).
        return str(fmt)
    if hasattr(value, "shape") and hasattr(value, "dtype"):
        return "dense"
    return "other"


def _dtype_name(value: Any) -> str:
    dtype = getattr(value, "dtype", None)
    if dtype is None:
        return "object"
    name = getattr(dtype, "name", str(dtype))
    if name == "category":
        return "categorical"
    if name in ("str", "string") or name.startswith("<U"):
        return "string"
    return str(name)


def _numeric_kind(value: Any) -> bool:
    dtype = getattr(value, "dtype", None)
    kind = getattr(dtype, "kind", None)
    return kind in ("i", "u", "f")


def _shape_of(value: Any) -> list[int] | None:
    shape = getattr(value, "shape", None)
    if not isinstance(shape, tuple):
        return None
    return [int(extent) for extent in shape]


def describe_obsm(obsm: Any, *, n_obs: int) -> list[dict[str, Any]]:
    """Describe every ``obsm`` entry with the listing-time verdict of §1.3."""
    described: list[dict[str, Any]] = []
    for key in obsm:
        value = obsm[key]
        storage = _storage_of(value)
        shape = _shape_of(value)
        reason: str | None = None
        if storage == "dataframe":
            reason = "obsm table, not an array"
        elif storage in ("csr", "csc"):
            reason = "sparse obsm arrays cannot be plotted"
        elif storage != "dense":
            reason = "not an array"
        elif shape is None or len(shape) not in (1, 2):
            reason = "not a one- or two-dimensional array"
        elif shape[0] != n_obs:
            reason = f"has {shape[0]:,} rows, the dataset has {n_obs:,} cells"
        elif not _numeric_kind(value):
            reason = "not a numeric array"
        elif len(shape) == 2 and shape[1] < 1:
            reason = "has no columns"
        described.append(
            {
                "key": str(key),
                "storage": storage,
                "shape": shape,
                "dtype": _dtype_name(value),
                "usable": reason is None,
                "reason": reason,
            }
        )
    return described


def obsm_entries(described: Sequence[Mapping[str, Any]]) -> list[ObsmEntry]:
    """Turn an ``obsm`` description into the entries the embedding rules read."""
    entries: list[ObsmEntry] = []
    for item in described:
        shape = item["shape"]
        width: int | None = None
        if item["storage"] == "dense" and shape is not None:
            width = 1 if len(shape) == 1 else shape[1] if len(shape) == 2 else None
        entries.append(
            ObsmEntry(key=item["key"], width=width, usable=item["usable"], reason=item["reason"])
        )
    return entries


def _describe_obsp(obsp: Any, *, n_obs: int) -> list[dict[str, Any]]:
    described: list[dict[str, Any]] = []
    for key in obsp:
        value = obsp[key]
        storage = _storage_of(value)
        shape = _shape_of(value)
        reason: str | None = None
        if storage not in ("dense", "csr", "csc"):
            reason = "not a matrix"
        elif shape is None or len(shape) != 2 or shape[0] != shape[1]:
            reason = "not square"
        elif shape[0] != n_obs:
            reason = f"has {shape[0]:,} rows, the dataset has {n_obs:,} cells"
        elif not _numeric_kind(value):
            reason = "not a numeric matrix"
        entry: dict[str, Any] = {
            "key": str(key),
            "storage": storage,
            "shape": shape,
            "dtype": _dtype_name(value),
            "nnz": int(value.nnz) if storage in ("csr", "csc") and hasattr(value, "nnz") else None,
            "usable": reason is None,
            "reason": reason,
        }
        described.append(entry)
    return described


def _describe_var_columns(var: Any) -> list[dict[str, Any]]:
    described: list[dict[str, Any]] = []
    for key in var.columns:
        column = var[key]
        dtype = _dtype_name(column)
        kind = getattr(getattr(column, "dtype", None), "kind", None)
        reason: str | None = None
        if kind == "f":
            reason = "floating-point values cannot name genes"
        elif kind == "b" or dtype == "bool":
            reason = "true/false values cannot name genes"
        described.append({"key": str(key), "dtype": dtype, "usable": reason is None, "reason": reason})
    return described


def _describe_matrix(value: Any) -> dict[str, Any]:
    return {"storage": _storage_of(value), "dtype": _dtype_name(value)}


def _describe_layers(adata: Any) -> list[dict[str, Any]]:
    """Every layer with its listing-time verdict, in the file's key order."""
    n_obs = int(adata.n_obs)
    n_vars = int(adata.n_vars)
    layers: list[dict[str, Any]] = []
    for name in adata.layers:
        value = adata.layers[name]
        shape = _shape_of(value)
        storage = _storage_of(value)
        reason: str | None = None
        if storage not in ("dense", "csr", "csc"):
            reason = "not a matrix"
        elif shape != [n_obs, n_vars]:
            reason = f"has shape {' x '.join(str(s) for s in shape or [])}, expected {n_obs} x {n_vars}"
        elif not _numeric_kind(value):
            reason = "not a numeric matrix"
        layers.append(
            {
                "key": str(name),
                "storage": storage,
                "dtype": _dtype_name(value),
                "usable": reason is None,
                "reason": reason,
            }
        )
    return layers


def default_expression_for(adata: Any) -> tuple[str | None, bool]:
    """The default expression source of one AnnData object, from metadata.

    :func:`default_expression` over what the object holds: ``X`` when present;
    otherwise ``raw``, the only usable layer, or the best-named one, guessed.
    """
    return default_expression(
        has_x=adata.X is not None,
        has_raw=getattr(adata, "raw", None) is not None,
        layers=[layer["key"] for layer in _describe_layers(adata) if layer["usable"]],
    )


def build_anndata_candidates(
    adata: Any,
    *,
    latent_key: str | None = None,
    knn_key: str | None = None,
) -> dict[str, Any]:
    """Build the candidate catalog of one AnnData object from metadata alone.

    Shapes, dtypes and storage kinds are read; no array values are. The
    ``defaults`` entry is the mapping §1.4 resolves from that metadata (the
    best guess, when it applies, uses listing-time usability only), and
    ``guessed`` names every group that was guessed rather than defaulted.
    """
    n_obs = int(adata.n_obs)
    n_vars = int(adata.n_vars)
    obsm = describe_obsm(adata.obsm, n_obs=n_obs)
    obsp = _describe_obsp(adata.obsp, n_obs=n_obs)

    layers = _describe_layers(adata)

    x_value = adata.X
    raw = getattr(adata, "raw", None)
    raw_entry: dict[str, Any] | None = None
    if raw is not None:
        raw_entry = {
            "n_vars": int(raw.n_vars),
            "X": _describe_matrix(raw.X),
            "var": {"columns": _describe_var_columns(raw.var)},
        }

    entries = obsm_entries(obsm)
    guess = resolve_embedding_slots(entries, latent_key=latent_key)
    slots = guess.slots
    guessed: list[str] = list(guess.guessed)

    expression, expression_guessed = default_expression_for(adata)
    if expression_guessed:
        guessed.append("expression")

    defaults: dict[str, Any] = {
        "kind": ANNDATA_KIND,
        "embeddings": {
            "enabled": default_enabled_embeddings(
                entries,
                slots,
                latent_key=latent_key,
                rejected=(key for key, _ in guess.rejected),
            ),
            "defaults": slots,
        },
        "knn": knn_key,
        "geneNames": "index",
        "expression": expression,
        "latent": latent_key,
        "vectorFields": default_vector_bindings(
            (entry.key for entry in entries if entry.usable), slots
        ),
    }
    return {
        "kind": ANNDATA_KIND,
        "n_obs": n_obs,
        "obsm": obsm,
        "obsp": obsp,
        "var": {"n_vars": n_vars, "columns": _describe_var_columns(adata.var)},
        "X": None if x_value is None else _describe_matrix(x_value),
        "layers": layers,
        "raw": raw_entry,
        "defaults": defaults,
        "guessed": guessed,
    }


# =============================================================================
# Resolving a mapping against an AnnData object (sessions)
# =============================================================================


def mapped_var_axis(adata: Any, mapping: Mapping[str, Any]) -> Any:
    """Return the ``var`` frame the mapping's gene axis lives on.

    ``raw`` has its own gene axis, so an AnnData mapping that reads expression
    from ``raw`` names its genes from ``raw.var``; everything else, and every
    export mapping, uses ``var``.
    """
    if mapping["kind"] == ANNDATA_KIND and mapping["expression"] == "raw":
        raw = getattr(adata, "raw", None)
        if raw is None:
            raise ValueError("The data mapping reads expression from raw, and adata.raw is absent.")
        return raw.var
    return adata.var


def mapped_var_count(adata: Any, mapping: Mapping[str, Any]) -> int:
    """Return how many genes the mapping exposes: the size of its gene axis.

    An AnnData mapping without an expression matrix exposes none.
    """
    if mapping["kind"] == ANNDATA_KIND:
        if mapping["expression"] is None:
            return 0
        if mapping["expression"] == "raw":
            raw = getattr(adata, "raw", None)
            if raw is None:
                raise ValueError(
                    "The data mapping reads expression from raw, and adata.raw is absent."
                )
            return int(raw.n_vars)
    return int(adata.n_vars)


def resolve_gene_names(adata: Any, mapping: Mapping[str, Any]) -> UniqueNames:
    """Resolve the gene names a mapping shows, with the shared make-unique rule.

    The names come from the mapping's gene axis: its index, or one of its
    columns, where an empty value takes the row's index name. Index names are
    stringified the way anndata stores them.
    """
    var = mapped_var_axis(adata, mapping)
    column = gene_names_column(mapping["geneNames"])
    if column is not None and column not in var.columns:
        raise KeyError(f"The data mapping names genes from var column {column!r}, which is absent.")
    return column_gene_names(var, column)


def mapped_expression_matrix(adata: Any, mapping: Mapping[str, Any]) -> Any:
    """Return the expression matrix a mapping reads, or raise when it has none."""
    if mapping["kind"] == EXPORT_KIND:
        source: str | None = "X"
    else:
        source = mapping["expression"]
    if source is None:
        raise ValueError("The data mapping reads no expression matrix.")
    if source == "X":
        matrix = getattr(adata, "X", None)
        if matrix is None:
            raise ValueError("The data mapping reads expression from X, and adata.X is absent.")
        return matrix
    if source == "raw":
        raw = getattr(adata, "raw", None)
        if raw is None or raw.X is None:
            raise ValueError("The data mapping reads expression from raw, and adata.raw is absent.")
        return raw.X
    layer = expression_layer(source)
    if layer is None or layer not in adata.layers:
        raise ValueError(f"The data mapping reads expression from {source!r}, which is absent.")
    return adata.layers[layer]


def require_mapping_matches_anndata(adata: Any, mapping: Mapping[str, Any]) -> None:
    """Require every key an AnnData mapping names to exist on ``adata``.

    Checks are metadata only: the keys exist and embedding columns fit the
    array width. An export mapping names exported embedding ids, which an
    AnnData object need not carry, so only its gene-name column is checked.
    """
    problems: list[str] = []
    gene_column = gene_names_column(mapping["geneNames"])
    if mapping["kind"] == ANNDATA_KIND:
        n_obs = int(adata.n_obs)
        described = {item["key"]: item for item in describe_obsm(adata.obsm, n_obs=n_obs)}
        refs: list[tuple[str, Mapping[str, Any]]] = [
            (f"embeddings.enabled[{position}]", ref)
            for position, ref in enumerate(mapping["embeddings"]["enabled"])
        ]
        refs.extend(
            (f"vectorFields embedding of {binding['key']!r}", binding["embedding"])
            for binding in mapping["vectorFields"]
        )
        for label, ref in refs:
            item = described.get(ref["key"])
            if item is None:
                problems.append(f"{label} uses obsm[{ref['key']!r}], which this dataset does not have")
                continue
            entry = obsm_entries([item])[0]
            if not entry.usable or entry.width is None:
                problems.append(f"{label} uses obsm[{ref['key']!r}], which {item['reason']}")
                continue
            try:
                validate_embedding_columns(
                    list(ref["columns"]),
                    dimension=len(ref["columns"]),
                    width=entry.width,
                )
            except DataMappingError as error:
                problems.append(f"{label} uses obsm[{ref['key']!r}]: {error}")
        for binding in mapping["vectorFields"]:
            if binding["key"] not in described:
                problems.append(
                    f"vector field obsm[{binding['key']!r}] is not in this dataset"
                )
        if mapping["knn"] is not None and mapping["knn"] not in adata.obsp:
            problems.append(f"knn uses obsp[{mapping['knn']!r}], which this dataset does not have")
        if mapping["latent"] is not None and mapping["latent"] not in adata.obsm:
            problems.append(
                f"latent uses obsm[{mapping['latent']!r}], which this dataset does not have"
            )
        expression = mapping["expression"]
        if expression == "X" and adata.X is None:
            problems.append("expression uses X, which this dataset does not have")
        elif expression == "raw" and getattr(adata, "raw", None) is None:
            problems.append("expression uses raw, which this dataset does not have")
        elif expression is not None:
            layer = expression_layer(expression)
            if layer is not None and layer not in adata.layers:
                problems.append(
                    f"expression uses layers[{layer!r}], which this dataset does not have"
                )
    if gene_column is not None:
        try:
            var = mapped_var_axis(adata, mapping)
        except ValueError:
            var = None
        if var is not None and gene_column not in var.columns:
            problems.append(
                f"geneNames uses var column {gene_column!r}, which this dataset does not have"
            )
    if problems:
        raise ValueError(
            "The session's data mapping does not match this AnnData: " + "; ".join(problems) + "."
        )
