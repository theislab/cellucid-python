"""
The data mapping: which stored arrays Cellucid plots and reads.

A *data mapping* names, for one dataset, the embeddings it can be shown with
(any number, of any mix of 1D, 2D and 3D: the *enabled* representations), the
default representation of each dimension, the vector fields bound to those
embeddings, the neighbour graph, the gene-name set, the expression matrix and
the latent space, each named by its key in one AnnData object. The same JSON
shape is written into session fingerprints and sent over the Jupyter wire, so
it is validated once, here, and mirrored by the web app's
``assets/js/data/data-mapping.js``.

This module holds the parts of those rules the Python package applies: the
DataMapping validation (sessions and the notebook wire), making names unique
and the note that reports it (``prepare()`` and sessions), the embedding-column
check and the default UMAP rule (the server and ``prepare()``), the
vector-field key rule (``prepare()``, the launch defaults and the transition
drift utilities), and resolving a session's mapping against an AnnData object.
Everything is a pure rule over metadata, except :func:`describe_obsm` (shapes
and dtypes only) and :func:`resolve_gene_names` (one ``var`` column). The rules
that must agree byte-for-byte between Python, R and JavaScript are pinned by
the shared conformance fixture ``tests/fixtures/data_mapping_conformance.json``.

This module imports only the standard library at import time, because the
session bundle reader depends on it and must stay light.
"""

from __future__ import annotations

import copy
import math
import numbers
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from ._plural import plural

#: The three dimension slots, in the order every mapping and identity lists them.
DIMENSION_SLOTS: tuple[str, ...] = ("1d", "2d", "3d")

#: How a caller names the embedding of one slot when serving AnnData: an
#: ``obsm`` key (read at columns ``0 .. d-1``) or a ``(key, columns)`` pair
#: with 0-based columns.
EmbeddingKeySpec = str | tuple[str, Sequence[int]]

#: The kind every data mapping declares: Cellucid reads AnnData, and a mapping
#: names ``obsm``, ``obsp``, layer and ``var`` keys of one AnnData object.
ANNDATA_KIND = "anndata"

_ANNDATA_MAPPING_KEYS = frozenset(
    {
        "kind",
        "embeddings",
        "knn",
        "geneNames",
        "expression",
        "log1p",
        "highlyVariable",
        "latent",
        "vectorFields",
    }
)
#: The one key a mapping holds only for a dataset of several modalities
#: (feature spaces over the same cells): which one gene names, expression, its
#: log1p and the highly variable genes are read from. A mapping of a dataset
#: with one feature space has no ``modality`` key at all.
MODALITY_KEY = "modality"
_LOG1P_ENTRY_KEYS = frozenset({"base"})
_EMBEDDING_REF_KEYS = frozenset({"key", "columns"})
_EMBEDDINGS_GROUP_KEYS = frozenset({"enabled", "defaults"})
_VECTOR_BINDING_KEYS = frozenset({"key", "embedding"})

#: ``X_umap_1d`` .. ``X_umap_3d`` hold one UMAP per dimension count; they are
#: the default rule first (design §1.4).
DIMENSIONAL_UMAP_KEYS: tuple[str, ...] = ("X_umap_1d", "X_umap_2d", "X_umap_3d")

#: The conventional names of one UMAP, in the order the default rule takes
#: them: Scanpy's (``sc.tl.umap()`` writes ``X_umap``, as does a Seurat
#: object's ``umap`` reduction), a SingleCellExperiment's (scater's
#: ``runUMAP()`` names its reducedDim ``UMAP``), then a Seurat object's
#: weighted-nearest-neighbour and Azimuth reference UMAPs (``wnn.umap``,
#: ``ref.umap``). An R object's reductions and reducedDims are
#: ``obsm['X_<name>']`` with the name's case kept.
UMAP_KEYS: tuple[str, ...] = ("X_umap", "X_UMAP", "X_wnn.umap", "X_ref.umap")

#: The prefix AnnData names a representation with: ``obsm['X_<basis>']`` is the
#: embedding (or latent space) of basis ``<basis>``, as ``sc.tl.umap()`` writes
#: ``X_umap``. A key with this prefix is never a vector field.
EMBEDDING_KEY_PREFIX = "X_"


def obsm_namespace(key: str) -> tuple[str, str]:
    """The namespace of one ``obsm`` key: ``(prefix, local)``.

    A MuData modality's arrays are named ``<modality>:<key>`` (the web app's
    MuData view, ``prepare()`` and ``cellucid_prepare()`` alike), each
    modality's names its own; a key without ``:`` is the dataset's own
    (prefix ``""``). The key rules (:func:`parse_vector_field_key`, the default
    and best-guess embeddings) apply within one namespace. The web app's
    ``obsmNamespace`` and cellucid-r apply the same rule.
    """
    colon = key.find(":")
    if 0 < colon < len(key) - 1:
        return key[: colon + 1], key[colon + 1 :]
    return "", key


@dataclass(frozen=True)
class VectorFieldKey:
    """One ``obsm`` key read as a vector field: ``<field>_<basis>`` on ``X_<basis>``.

    ``field`` keeps the key's namespace (``rna:velocity``); ``basis`` is the
    name's own, and the embedding is in the same namespace.
    """

    field: str
    basis: str

    @property
    def embedding_key(self) -> str:
        """The ``obsm`` key of the embedding the field is named for, in the field's namespace."""
        return obsm_namespace(self.field)[0] + EMBEDDING_KEY_PREFIX + self.basis


def parse_vector_field_key(key: object, obsm_keys: Collection[str]) -> VectorFieldKey | None:
    """Read one ``obsm`` key with the vector-field key rule, scVelo's convention.

    ``<field>_<basis>`` is the vector field ``<field>`` on the embedding
    ``X_<basis>`` (``velocity_umap`` on ``X_umap``, ``velocity_umap_3d`` on
    ``X_umap_3d``, ``T_fwd_branches`` on ``X_branches``), within one namespace
    (:func:`obsm_namespace`): a modality's ``rna:velocity_umap`` is the field
    ``rna:velocity`` on ``rna:X_umap``, never on the dataset's ``X_umap``.

    - A key whose own name (after its namespace) starts with ``X_`` is never a
      vector field: AnnData names representations that way.
    - Otherwise every split of the name at one ``_`` into a non-empty
      ``<field>`` and a non-empty ``<basis>`` is tried, and the split whose
      ``X_<basis>``, in the same namespace, is in ``obsm_keys`` wins; when
      several do, the longest basis wins. Matching is exact and
      case-sensitive.
    - A key with no such split is an ordinary ``obsm`` array (None).

    The rule reads names only; whether the field is drawn also depends on
    column counts (the web app's ``defaultVectorBindings``). The web app's
    ``parseVectorFieldKey`` and cellucid-r apply the same rule, pinned by the
    shared conformance fixture.
    """
    if not isinstance(key, str):
        return None
    prefix, local = obsm_namespace(key)
    if local.startswith(EMBEDDING_KEY_PREFIX):
        return None
    split = local.find("_")
    while split != -1:
        # The first qualifying split from the left has the longest basis.
        if 0 < split < len(local) - 1:
            basis = local[split + 1 :]
            if prefix + EMBEDDING_KEY_PREFIX + basis in obsm_keys:
                return VectorFieldKey(field=prefix + local[:split], basis=basis)
        split = local.find("_", split + 1)
    return None


def vector_field_keys(keys: Iterable[str]) -> dict[str, VectorFieldKey]:
    """The vector fields among one object's ``obsm`` keys, by key, in input order."""
    listed = list(keys)
    obsm_keys = frozenset(listed)
    found: dict[str, VectorFieldKey] = {}
    for key in listed:
        parsed = parse_vector_field_key(key, obsm_keys)
        if parsed is not None:
            found[key] = parsed
    return found

_GENE_NAMES_COLUMN_PREFIX = "column:"
_LAYER_PREFIX = "layer:"


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

    A column must hold text, categories or integers; its values are turned into text with
    :func:`gene_name_text`, an empty one takes the row's index name (which
    must be text), and the list is made unique with
    :func:`make_unique_names`. ``prepare()`` and the session reader name
    genes with this one function, the rule the web app's ``geneNameText``
    and cellucid-r apply too.
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

    A list that is already unique comes back unchanged.
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
# Embedding references and columns (design §1.1, §1.3)
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
            f"A {dimension}D embedding needs at least {dimension} {plural(dimension, 'column')}; "
            f"this array has {width}.",
        )
    if not isinstance(columns, list | tuple):
        raise DataMappingError("columns_not_list", "Embedding columns must be a list.")
    if len(columns) != dimension:
        raise DataMappingError(
            "column_count",
            f"A {dimension}D embedding reads exactly {dimension} {plural(dimension, 'column')}; got {len(columns)}.",
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
    dimension: int | None,
    label: str,
) -> dict[str, Any]:
    """Validate one ``{key, columns}`` reference and return a canonical copy.

    A reference names an ``obsm`` key and lists its 0-based columns.
    """
    ref = _require_exact_keys(value, _EMBEDDING_REF_KEYS, label=label)
    key = _require_mapping_string(ref["key"], label=f"{label}.key")
    columns = ref["columns"]
    if not isinstance(columns, list):
        raise DataMappingError("columns_not_list", f"{label}.columns must be a list.")
    if dimension is not None and len(columns) != dimension:
        raise DataMappingError(
            "column_count",
            f"{label}.columns must list exactly {dimension} {plural(dimension, 'column')}.",
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


def _is_expression_matrix(value: object) -> bool:
    return isinstance(value, str) and (
        value in ("X", "raw")
        or (value.startswith(_LAYER_PREFIX) and len(value) > len(_LAYER_PREFIX))
    )


def _validate_log1p(value: object, *, label: str) -> dict[str, dict[str, Any]]:
    """Validate the ``log1p`` group (the web app's ``validateLog1p``).

    For each expression matrix holding ``log1p`` values (``X``, ``raw``,
    ``layer:<name>``), its ``{"base": ...}``: None for the natural logarithm,
    else a positive number other than 1. A matrix not listed holds its values
    as stored. The result lists the matrices in code-point order.
    """
    if not isinstance(value, dict):
        raise DataMappingError("log1p", f"{label} must be an object of expression matrices.")
    result: dict[str, dict[str, Any]] = {}
    for matrix in sorted(value):
        if not _is_expression_matrix(matrix):
            raise DataMappingError(
                "log1p",
                f'{label} keys must be expression matrices: "X", "raw" or "layer:<name>".',
            )
        entry = _require_exact_keys(value[matrix], _LOG1P_ENTRY_KEYS, label=f"{label}.{matrix}")
        base = entry["base"]
        if base is not None and not (
            isinstance(base, int | float)
            and not isinstance(base, bool)
            and math.isfinite(base)
            and base > 0
            and base != 1
        ):
            raise DataMappingError(
                "log1p_base",
                f"{label}.{matrix}.base must be null (base e) or a positive number other than 1.",
            )
        result[matrix] = {"base": base}
    return result


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

    By ``key`` in code-point order, then ``columns`` compared element by
    element.
    """
    return (ref["key"], tuple(ref["columns"]))


def _validate_embeddings_group(value: object, *, label: str) -> dict[str, Any]:
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
        validate_embedding_ref(ref, dimension=None, label=f"{label}.enabled[{position}]")
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


def _binding_sort_key(binding: Mapping[str, Any]) -> tuple[Any, ...]:
    return (*embedding_ref_sort_key(binding["embedding"]), binding["key"])


def _validate_modality(value: object, *, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise DataMappingError(
            "modality",
            f"{label} must name one of the dataset's modalities; a dataset with one feature "
            "space has no modality key.",
        )
    return value


def validate_data_mapping(value: object, *, label: str = "dataMapping") -> dict[str, Any]:
    """Validate one DataMapping and return a canonical deep copy.

    The mapping is a closed shape: ``kind`` is ``"anndata"``, and an unknown
    or missing key is a mapping this version does not understand. At least one
    embedding must be enabled and at least one default named, since a mapping
    that plots nothing describes no session; every vector binding is bound to
    an enabled embedding.

    ``vectorFields`` must already be in canonical order -- by the bound
    embedding's ``key``, then its ``columns``, then the vector's ``obsm`` key --
    and must not repeat a binding.

    ``log1p`` names, per expression matrix holding ``log1p`` values, its base
    (None for e); ``highlyVariable`` the boolean ``var`` column of the mapped
    expression matrix's genes that marks the highly variable ones, or None.

    ``modality`` is present only for a dataset of several modalities (a
    non-empty string naming one); ``geneNames``, ``expression``, ``log1p`` and
    ``highlyVariable`` are then that modality's.

    The canonical key order is ``kind``, ``embeddings``, ``knn``, ``modality``
    (when present), ``geneNames``, ``expression``, ``log1p``,
    ``highlyVariable``, ``latent``, ``vectorFields``.
    """
    if not isinstance(value, dict):
        raise DataMappingError("not_object", f"{label} must be an object.")
    if value.get("kind") != ANNDATA_KIND:
        raise DataMappingError("kind", f"{label}.kind must be \"{ANNDATA_KIND}\".")
    has_modality = MODALITY_KEY in value
    mapping = _require_exact_keys(
        value,
        _ANNDATA_MAPPING_KEYS | {MODALITY_KEY} if has_modality else _ANNDATA_MAPPING_KEYS,
        label=label,
    )

    embeddings = _validate_embeddings_group(mapping["embeddings"], label=f"{label}.embeddings")

    raw_bindings = mapping["vectorFields"]
    if not isinstance(raw_bindings, list):
        raise DataMappingError("not_list", f"{label}.vectorFields must be a list.")
    bindings: list[dict[str, Any]] = []
    for position, raw in enumerate(raw_bindings):
        binding_label = f"{label}.vectorFields[{position}]"
        binding = _require_exact_keys(raw, _VECTOR_BINDING_KEYS, label=binding_label)
        bindings.append(
            {
                "key": _require_mapping_string(binding["key"], label=f"{binding_label}.key"),
                "embedding": validate_embedding_ref(
                    binding["embedding"],
                    dimension=None,
                    label=f"{binding_label}.embedding",
                ),
            }
        )
    sort_keys = [_binding_sort_key(binding) for binding in bindings]
    if sort_keys != sorted(sort_keys):
        raise DataMappingError(
            "vector_order",
            f"{label}.vectorFields must be sorted by embedding key, embedding columns, "
            "then key.",
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

    canonical: dict[str, Any] = {
        "kind": ANNDATA_KIND,
        "embeddings": embeddings,
        "knn": _optional_mapping_string(mapping["knn"], label=f"{label}.knn"),
        **(
            {"modality": _validate_modality(mapping[MODALITY_KEY], label=f"{label}.modality")}
            if has_modality
            else {}
        ),
        "geneNames": _validate_gene_names(mapping["geneNames"], label=f"{label}.geneNames"),
        "expression": _validate_expression(mapping["expression"], label=f"{label}.expression"),
        "log1p": _validate_log1p(mapping["log1p"], label=f"{label}.log1p"),
        "highlyVariable": _optional_mapping_string(
            mapping["highlyVariable"], label=f"{label}.highlyVariable"
        ),
        "latent": _optional_mapping_string(mapping["latent"], label=f"{label}.latent"),
        "vectorFields": bindings,
    }
    return copy.deepcopy(canonical)


#: The names a gene-expression modality goes by, lower-cased: MuData's ``rna``,
#: 10x's ``Gene Expression``, Seurat's ``RNA`` and ``SCT``.
_GENE_EXPRESSION_MODALITIES = ("rna", "gex", "gene expression", "gene_expression", "sct")


def default_modality(names: Sequence[str]) -> str | None:
    """The modality a dataset opens with, of its modalities in stored order.

    The first whose name is a gene-expression name (``rna``, ``gex``, ``gene
    expression``, ``gene_expression``, ``sct``, in any case), else the first;
    None when there are none. The web app's ``defaultModality``.
    """
    for name in names:
        if name.lower() in _GENE_EXPRESSION_MODALITIES:
            return name
    return names[0] if names else None


#: The neighbour graph AnnData's ``sc.pp.neighbors`` writes, the default first.
DEFAULT_GRAPH_KEY = "connectivities"

_SNN_GRAPH_SUFFIX = "_snn"


def default_graph(keys: Iterable[str], *, modality: str | None = None) -> str | None:
    """The neighbour graph a dataset opens with, of its ``obsp`` keys.

    Scanpy's ``connectivities``; else Seurat's shared-nearest-neighbour graph
    of the modality (assay) it opens with (``<modality>_snn``), else the first
    other shared-nearest-neighbour graph (``<assay>_snn``, or ``SNN``, in any
    case) in code-point order; else None. Seurat's other graph
    (``<assay>_nn``, directed) and Scanpy's ``distances`` are no undirected
    neighbour graph and are never the default. ``prepare()`` writes the graph
    this rule chooses as ``obsp['connectivities']``; the web app's
    ``defaultGraph`` and cellucid-r apply the same rule, pinned by the shared
    conformance fixture.
    """
    listed = list(keys)
    if DEFAULT_GRAPH_KEY in listed:
        return DEFAULT_GRAPH_KEY
    if modality is not None and f"{modality}{_SNN_GRAPH_SUFFIX}" in listed:
        return f"{modality}{_SNN_GRAPH_SUFFIX}"
    shared = sorted(
        key
        for key in listed
        if key.lower() == "snn" or (key.lower().endswith(_SNN_GRAPH_SUFFIX) and len(key) > len(_SNN_GRAPH_SUFFIX))
    )
    return shared[0] if shared else None


def require_mapping_modality(mapping: Mapping[str, Any], modalities: Sequence[str] | None) -> None:
    """Check a mapping's ``modality`` against the dataset's modalities.

    A dataset of several modalities needs one named, a dataset of one feature
    space (``modalities`` None) none, and the name must be one of them.
    Raises :class:`DataMappingError` with code ``modality_missing``,
    ``modality_unexpected`` or ``modality_unknown``.
    """
    named = MODALITY_KEY in mapping
    if modalities is None:
        if named:
            raise DataMappingError(
                "modality_unexpected",
                f"This mapping reads the modality {mapping[MODALITY_KEY]!r}, and this dataset "
                "has one feature space.",
            )
        return
    if not named:
        raise DataMappingError(
            "modality_missing",
            f"This dataset has several modalities ({', '.join(modalities)}); its mapping must "
            "name the one it reads.",
        )
    if mapping[MODALITY_KEY] not in modalities:
        raise DataMappingError(
            "modality_unknown",
            f"This mapping reads the modality {mapping[MODALITY_KEY]!r}, which this dataset "
            f"does not have ({', '.join(modalities)}).",
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
# The default embedding rule (design §1.4)
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


def default_embedding_slots(
    entries: Sequence[ObsmEntry], *, modality: str | None = None
) -> dict[str, dict[str, Any] | None]:
    """Apply the default embedding rule and return the filled slots.

    If any of ``X_umap_1d``, ``X_umap_2d`` or ``X_umap_3d`` is present, those
    keys alone fill their slots, reading columns ``0 .. d-1``. Otherwise the
    first of ``UMAP_KEYS`` (``X_umap``, ``X_UMAP``, ``X_wnn.umap``,
    ``X_ref.umap``) whose width is 1, 2 or 3 fills the slot of that width. A
    dataset of several modalities whose own names fill no slot is given the
    same names in the namespace of the modality it opens with (``modality``:
    ``rna:X_umap``, :func:`obsm_namespace`). Slots without a default stay None.
    Presence is all this rule asks: a present key that later fails its
    load-time check is reported, not skipped, because the caller named that
    array on purpose.
    """
    by_key = {entry.key: entry for entry in entries}
    for prefix in ("",) if modality is None else ("", f"{modality}:"):
        slots: dict[str, dict[str, Any] | None] = dict.fromkeys(DIMENSION_SLOTS)
        dimensional = [key for key in DIMENSIONAL_UMAP_KEYS if prefix + key in by_key]
        if dimensional:
            for key in dimensional:
                dimension = int(key[len("X_umap_")])
                slots[dimension_slot(dimension)] = {"key": prefix + key, "columns": default_columns(dimension)}
            return slots
        for key in UMAP_KEYS:
            single = by_key.get(prefix + key)
            if single is not None and single.width in (1, 2, 3):
                dimension = single.width
                slots[dimension_slot(dimension)] = {"key": prefix + key, "columns": default_columns(dimension)}
                return slots
    return dict.fromkeys(DIMENSION_SLOTS)


#: The web app's best-guess ranking of ``obsm`` keys (``bestGuessTier`` in
#: ``data/data-mapping.js``): the first needle a lower-cased key contains
#: decides its tier, a lower tier is a better guess, and a key matching none
#: ranks between ``diffmap`` and the latent spaces.
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


def best_guess_order(keys: Iterable[str], *, modality: str | None = None) -> list[str]:
    """Order ``obsm`` keys for the best guess: by tier, then namespace, then code point.

    A key ranks by its own name (``rna:X_umap`` is a UMAP); among keys of one
    tier the dataset's own come first, then the modality it opens with
    (``modality``), then the other modalities (:func:`obsm_namespace`).
    """

    def tier(key: str) -> int:
        name = obsm_namespace(key)[1].lower()
        for rank, needles in _BEST_GUESS_TIERS:
            if any(needle in name for needle in needles):
                return rank
        return _BEST_GUESS_NO_MATCH_TIER

    def namespace(key: str) -> int:
        prefix = obsm_namespace(key)[0]
        if prefix == "":
            return 0
        return 1 if modality is not None and prefix == f"{modality}:" else 2

    return sorted(keys, key=lambda key: (tier(key), namespace(key), key))


def enabled_embedding_keys(
    entries: Sequence[ObsmEntry],
    *,
    named_slots: Mapping[str, str],
    latent_key: str | None,
    modality: str | None = None,
) -> frozenset[str]:
    """The ``obsm`` keys the web app shows an AnnData store with when it opens.

    The slots are the ones a launch named (``named_slots``, slot to key) over
    the default rule (:func:`default_embedding_slots`); when a launch names
    none and the default rule fills none, the best guess fills them: for 3D,
    2D and 1D the first key of that exact width in :func:`best_guess_order`,
    else the first wider key in 2D. Every usable array 1 to 3 columns wide is
    shown beside them. Vector-field keys and the latent key are never guessed
    or shown beside. This is the web app's ``resolveEmbeddingSlots`` and
    ``defaultEnabledEmbeddings`` read from metadata alone: the browser also
    loads each guessed array and passes over one that does not load, which
    metadata cannot tell. ``modality``: the modality a dataset of several
    opens with (:func:`default_embedding_slots`, :func:`best_guess_order`).
    """
    vector_keys = vector_field_keys(entry.key for entry in entries)
    candidates = {
        entry.key: entry.width
        for entry in entries
        if entry.usable
        and entry.width is not None
        and entry.width >= 1
        and entry.key not in vector_keys
        and entry.key != latent_key
    }
    slots = default_embedding_slots(entries, modality=modality)
    if named_slots:
        for slot, key in named_slots.items():
            slots[slot] = {"key": key}
    elif all(slots[slot] is None for slot in DIMENSION_SLOTS):
        for dimension in (3, 2, 1):
            exact = best_guess_order(
                (key for key, width in candidates.items() if width == dimension), modality=modality
            )
            if exact:
                slots[dimension_slot(dimension)] = {"key": exact[0]}
        if all(slots[slot] is None for slot in DIMENSION_SLOTS):
            wide = best_guess_order((key for key, width in candidates.items() if width > 3), modality=modality)
            if wide:
                slots[dimension_slot(2)] = {"key": wide[0]}
    shown = {ref["key"] for ref in slots.values() if ref is not None}
    shown.update(key for key, width in candidates.items() if width <= 3)
    return frozenset(shown)


# =============================================================================
# Describing obsm and var -- metadata only
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
            reason = f"has {shape[0]:,} {plural(shape[0], 'row')}, the dataset has {n_obs:,} {plural(n_obs, 'cell')}"
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


# =============================================================================
# Resolving a mapping against an AnnData object (sessions)
# =============================================================================


def mapped_var_axis(adata: Any, mapping: Mapping[str, Any]) -> Any:
    """Return the ``var`` frame the mapping's gene axis lives on.

    ``raw`` has its own gene axis, so a mapping that reads expression from
    ``raw`` names its genes from ``raw.var``; everything else uses ``var``.
    """
    if mapping["expression"] == "raw":
        raw = getattr(adata, "raw", None)
        if raw is None:
            raise ValueError("The data mapping reads expression from raw, and adata.raw is absent.")
        return raw.var
    return adata.var


def mapped_var_count(adata: Any, mapping: Mapping[str, Any]) -> int:
    """Return how many genes the mapping exposes: the size of its gene axis.

    A mapping without an expression matrix exposes none.
    """
    if mapping["expression"] is None:
        return 0
    if mapping["expression"] == "raw":
        raw = getattr(adata, "raw", None)
        if raw is None:
            raise ValueError("The data mapping reads expression from raw, and adata.raw is absent.")
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
    source: str | None = mapping["expression"]
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
    array width. An AnnData object has one feature space, so its mapping
    names no modality.
    """
    require_mapping_modality(mapping, None)
    problems: list[str] = []
    gene_column = gene_names_column(mapping["geneNames"])
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
    for matrix in mapping["log1p"]:
        if matrix == "X" and adata.X is None:
            problems.append("log1p names X, which this dataset does not have")
        elif matrix == "raw" and getattr(adata, "raw", None) is None:
            problems.append("log1p names raw, which this dataset does not have")
        elif expression_layer(matrix) is not None and expression_layer(matrix) not in adata.layers:
            problems.append(
                f"log1p names layers[{expression_layer(matrix)!r}], which this dataset does not have"
            )
    highly_variable = mapping["highlyVariable"]
    if highly_variable is not None:
        try:
            var = mapped_var_axis(adata, mapping)
        except ValueError:
            var = None
        if var is not None and highly_variable not in var.columns:
            problems.append(
                f"highlyVariable uses var column {highly_variable!r}, which this dataset does not have"
            )
        elif var is not None and str(var[highly_variable].dtype) not in ("bool", "boolean"):
            problems.append(
                f"highlyVariable uses var column {highly_variable!r}, which is not a true/false column"
            )
    if problems:
        raise ValueError(
            "The session's data mapping does not match this dataset: " + "; ".join(problems) + "."
        )
