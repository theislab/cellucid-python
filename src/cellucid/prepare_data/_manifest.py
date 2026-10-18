"""
The compact_v2 on-disk manifest format.

One place names what the export format calls things -- the format version, the
two payload directories, the codes extension of each categorical dtype -- and
one place re-derives, from an emitted manifest alone, exactly which payload
files a reader will go looking for. The writer checks itself against that
derivation before it publishes, so a manifest can never describe a directory
that was not written. The same check reconciles the ``embeddings/`` directory,
whose declaration is ``dataset_identity.json``, and the export root itself.

Those checks all run on the payload in memory, which proves nothing about the
file a reader opens. ``json.dumps`` is not a total function into the JSON the
browser accepts -- a non-finite float is written as the bare token ``NaN`` or
``Infinity``, which ``JSON.parse`` refuses and ``json.loads`` does not, and a
dict key that is not a string is silently coerced into one -- so every manifest
is read back under a parser as strict as the reader's and required to parse to
exactly the payload that was validated.
"""

import json
import math
from collections.abc import Sequence
from pathlib import Path

from .._contracts import _require_unique_identifiers

DEFAULT_OBS_DIRNAME = "obs"
DEFAULT_VAR_DIRNAME = "var"

# The format version both compact manifests (obs and var) carry. compact_v2
# keys categorical centroids by embedding id and adds gene-name sets to var.
MANIFEST_FORMAT_VERSION = "compact_v2"


def _reject_non_standard_json_constant(name: str) -> object:
    """Refuse the three tokens Python writes and JSON does not define."""
    raise ValueError(f"{name} is not a JSON value")


def _json_kind(value: object) -> str:
    """Name what one payload node is in JSON.

    ``bool`` is checked before ``int`` because it is a subclass of it, and JSON
    tells ``true`` and ``1`` apart even though Python compares them equal.
    """
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "boolean"
    if isinstance(value, str):
        return "string"
    if isinstance(value, (int, float)):
        return "number"
    if isinstance(value, dict):
        return "object"
    if isinstance(value, (list, tuple)):
        return "array"
    return "unsupported"


def _first_json_difference(expected: object, actual: object, path: str = "$") -> str | None:
    """Return where a parsed manifest stops matching its payload, or None.

    Numbers compare by value rather than by Python type, because an exported
    count is one JSON number whether it was built as ``int`` or ``float``. A
    number that became a *string* is caught by the kind, which is the case that
    matters.
    """
    expected_kind = _json_kind(expected)
    actual_kind = _json_kind(actual)
    if expected_kind == "unsupported":
        return f"{path} is not a JSON value"
    if expected_kind != actual_kind:
        return f"{path} was written as {actual_kind} but validated as {expected_kind}"
    if expected_kind == "null":
        return None
    if expected_kind == "object":
        assert isinstance(expected, dict) and isinstance(actual, dict)
        expected_keys = list(expected)
        actual_keys = list(actual)
        if any(not isinstance(key, str) for key in expected_keys):
            return f"{path} has a key that is not a string"
        if expected_keys != actual_keys:
            return (
                f"{path} keys are {actual_keys} but were validated as {expected_keys}"
            )
        for key in expected_keys:
            difference = _first_json_difference(expected[key], actual[key], f"{path}.{key}")
            if difference is not None:
                return difference
        return None
    if expected_kind == "array":
        assert isinstance(expected, (list, tuple)) and isinstance(actual, (list, tuple))
        if len(expected) != len(actual):
            return (
                f"{path} holds {len(actual)} elements but was validated "
                f"with {len(expected)}"
            )
        for index, (expected_item, actual_item) in enumerate(zip(expected, actual, strict=True)):
            difference = _first_json_difference(expected_item, actual_item, f"{path}[{index}]")
            if difference is not None:
                return difference
        return None
    if expected_kind == "number":
        assert isinstance(expected, (int, float)) and isinstance(actual, (int, float))
        if not math.isfinite(expected) or expected != actual:
            return f"{path} is {actual!r} but was validated as {expected!r}"
        return None
    if expected != actual:
        return f"{path} is {actual!r} but was validated as {expected!r}"
    return None


def _require_manifest_reads_back(path: Path, payload: object) -> None:
    """Require the manifest just written to say what the payload said.

    Everything else in this module validates the payload; this is the only
    place that validates the file. It runs inside the export transaction, so a
    generation whose manifest lost a value on the way to disk is rejected here
    instead of publishing and failing in the browser.
    """
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError as error:
        raise RuntimeError(
            f"{path.name} does not read back as the manifest that was "
            f"validated: the file is not valid UTF-8."
        ) from error
    try:
        parsed = json.loads(text, parse_constant=_reject_non_standard_json_constant)
    except ValueError as error:
        raise RuntimeError(
            f"{path.name} does not read back as the manifest that was validated: {error}"
        ) from error
    difference = _first_json_difference(payload, parsed)
    if difference is not None:
        raise RuntimeError(
            f"{path.name} does not read back as the manifest that was "
            f"validated: {difference}."
        )


def _write_manifest(path: Path, payload: object, *, indent: int | None = None) -> None:
    """Write one manifest and prove the bytes carry what was validated."""
    path.write_text(json.dumps(payload, indent=indent), encoding="utf-8")
    _require_manifest_reads_back(path, payload)


def _require_dense_payload_indices(
    indices: Sequence[object],
    *,
    axis: str,
) -> None:
    """Require the payload indices of one axis to be exactly 0..N-1, each once.

    The index *is* the filename, so two fields holding one index write into one
    file: the second overwrites the first and the viewer then draws one field's
    values under the other field's name. Nothing downstream can detect that, so
    it is asserted here in the writer, against the manifest that was just built,
    rather than only in a test.
    """
    resolved: list[int] = []
    for position, index in enumerate(indices):
        if type(index) is not int:
            raise RuntimeError(
                f"{axis} payload index at position {position} must be a native "
                f"integer, got {index!r}."
            )
        resolved.append(index)
    if sorted(resolved) != list(range(len(resolved))):
        raise RuntimeError(
            f"{axis} payload indices must be exactly 0..{len(resolved) - 1}, each "
            f"used once; got {sorted(resolved)}."
        )


def _expand_payload_pattern(
    pattern: object,
    *,
    index: int,
    label: str,
    extension: str | None = None,
) -> str:
    """Substitute one payload index, and one codes extension, into a pattern."""
    if not isinstance(pattern, str) or not pattern:
        raise RuntimeError(f"{label} must be a non-empty path pattern.")
    expanded = pattern.replace("{index}", str(index))
    if extension is not None:
        expanded = expanded.replace("{ext}", extension)
    if "{" in expanded or "}" in expanded:
        raise RuntimeError(f"{label} retains an unsubstituted placeholder: {expanded!r}.")
    return expanded


def _require_declared_payloads_on_disk(
    out_dir: Path,
    *,
    directory_name: str | None,
    declared: set[str],
    axis: str,
    declared_directories: set[str] | None = None,
) -> None:
    """Require one export directory to hold exactly the payloads it declares.

    The manifest is the only index the viewer has: a payload it does not declare
    is invisible, and a payload it declares but that was never written fails the
    dataset at read time, in the browser, long after the export succeeded. Both
    are caught here by re-expanding the emitted path patterns and comparing them
    against the directory that was actually written.

    ``directory_name=None`` reconciles the export root itself rather than one
    axis directory. The root holds only the manifests and the axis
    directories, so ``declared_directories`` names the subdirectories that
    belong to the generation: a directory that is not one of them is refused by
    kind, which also keeps a declared payload from being satisfied by a
    directory of the same name.
    """
    directory = out_dir if directory_name is None else out_dir / directory_name
    prefix = "" if directory_name is None else f"{directory_name}/"
    expected_directories = set() if declared_directories is None else declared_directories
    on_disk: set[str] = set()
    on_disk_directories: set[str] = set()
    if directory.exists():
        for entry in sorted(directory.iterdir()):
            relative = f"{prefix}{entry.name}"
            if entry.is_dir():
                if relative not in expected_directories:
                    raise RuntimeError(f"{axis} payload directory holds a non-file entry: {entry}")
                on_disk_directories.add(relative)
                continue
            if not entry.is_file():
                raise RuntimeError(f"{axis} payload directory holds a non-file entry: {entry}")
            on_disk.add(relative)
    missing = sorted((declared - on_disk) | (expected_directories - on_disk_directories))
    undeclared = sorted(on_disk - declared)
    if missing or undeclared:
        raise RuntimeError(
            f"{axis} manifest does not describe the payloads that were written. "
            f"Declared but absent: {missing}. Written but undeclared: {undeclared}."
        )


def _identity_obs_fields_from_compact_manifest(
    manifest: dict,
) -> list[dict]:
    expected_keys = {
        "_format",
        "n_points",
        "centroid_outlier_quantile",
        "latent_key",
        "compression",
        "_obsSchemas",
        "_continuousFields",
        "_categoricalFields",
    }
    if not isinstance(manifest, dict) or set(manifest) != expected_keys:
        raise ValueError(
            f"obs manifest must contain exactly the current {MANIFEST_FORMAT_VERSION} fields."
        )
    if manifest["_format"] != MANIFEST_FORMAT_VERSION:
        raise ValueError(f"obs manifest must use the current {MANIFEST_FORMAT_VERSION} format.")
    # The schema map is an object in every generation, including one carrying
    # no observation field at all, where it is empty. The viewer refuses a
    # non-object here, so the kind is asserted rather than only the contents.
    if _json_kind(manifest["_obsSchemas"]) != "object":
        raise ValueError("obs manifest _obsSchemas must be a JSON object.")

    continuous_fields = manifest["_continuousFields"]
    categorical_fields = manifest["_categoricalFields"]
    if not isinstance(continuous_fields, list):
        raise ValueError("compact obs manifest _continuousFields must be a list.")
    if not isinstance(categorical_fields, list):
        raise ValueError("compact obs manifest _categoricalFields must be a list.")

    identity_fields: list[dict] = []
    manifest_keys: list[str] = []
    # Both arrays write into obs/, so their payload indices share one space.
    payload_indices: list[object] = []
    for field in continuous_fields:
        if (
            not isinstance(field, list)
            or len(field) not in (2, 4)
            or type(field[0]) is not int
            or not isinstance(field[1], str)
            or not field[1]
        ):
            raise ValueError(
                "compact continuous observation fields must be exact "
                "[index, key] or [index, key, minValue, maxValue] tuples."
            )
        payload_indices.append(field[0])
        key = field[1]
        manifest_keys.append(key)
        identity_fields.append({"key": key, "kind": "continuous"})

    for field in categorical_fields:
        if (
            not isinstance(field, list)
            or len(field) not in (6, 8)
            or type(field[0]) is not int
            or not isinstance(field[1], str)
            or not field[1]
            or not isinstance(field[2], list)
            or _json_kind(field[5]) != "object"
        ):
            raise ValueError(
                "compact categorical observation fields must be exact six- or "
                "eight-member tuples with a category array and a centroid object "
                "keyed by embedding id."
            )
        payload_indices.append(field[0])
        key = field[1]
        manifest_keys.append(key)
        identity_fields.append(
            {
                "key": key,
                "kind": "category",
                "n_categories": len(field[2]),
            }
        )

    _require_unique_identifiers(manifest_keys, label="Observation field")
    _require_dense_payload_indices(payload_indices, axis="Observation")
    return identity_fields


def _gene_names_from_compact_manifest(manifest: dict) -> list[str]:
    """Validate one emitted var manifest and return its exact default gene names."""
    expected_keys = {
        "_format",
        "n_points",
        "compression",
        "quantization",
        "_varSchema",
        "fields",
        "gene_name_sets",
    }
    if not isinstance(manifest, dict) or set(manifest) != expected_keys:
        raise ValueError(
            f"var manifest must contain exactly the current {MANIFEST_FORMAT_VERSION} fields."
        )
    if manifest["_format"] != MANIFEST_FORMAT_VERSION:
        raise ValueError(f"var manifest must use the current {MANIFEST_FORMAT_VERSION} format.")
    fields = manifest["fields"]
    if not isinstance(fields, list):
        raise ValueError("compact var manifest fields must be a list.")

    names: list[str] = []
    payload_indices: list[object] = []
    for field in fields:
        if (
            not isinstance(field, list)
            or len(field) not in (2, 4)
            or type(field[0]) is not int
            or not isinstance(field[1], str)
            or not field[1]
        ):
            raise ValueError(
                "compact var fields must be exact [index, name] or "
                "[index, name, minValue, maxValue] tuples."
            )
        payload_indices.append(field[0])
        names.append(field[1])

    _require_unique_identifiers(names, label="Gene")
    _require_dense_payload_indices(payload_indices, axis="Gene")
    _require_gene_name_sets(manifest["gene_name_sets"], default_names=names)
    return names


_GENE_NAME_SET_KEYS = ["source", "filled_from_index", "renamed", "names"]


def _require_gene_name_sets(sets: object, *, default_names: list[str]) -> None:
    """Validate the gene-name sets of one var manifest against its field names.

    Set 0 is the default: its names are the field names, so its ``names`` is
    null. Every later set lists one name per field, null where it equals the
    default name, and resolves to names that are distinct within the set. No
    two sets read the same source.
    """
    if not isinstance(sets, list) or not sets:
        raise ValueError("var manifest gene_name_sets must be a non-empty list.")
    sources: list[object] = []
    for position, name_set in enumerate(sets):
        label = f"var manifest gene_name_sets[{position}]"
        if not isinstance(name_set, dict) or list(name_set) != _GENE_NAME_SET_KEYS:
            raise ValueError(f"{label} must hold exactly the keys {_GENE_NAME_SET_KEYS}.")
        source = name_set["source"]
        if source is not None and (not isinstance(source, str) or not source):
            raise ValueError(f"{label}.source must be null or a non-empty column name.")
        if position > 0 and source is None:
            raise ValueError(f"{label}.source must name a var column.")
        sources.append(source)
        filled = name_set["filled_from_index"]
        if type(filled) is not int or not 0 <= filled <= len(default_names):
            raise ValueError(f"{label}.filled_from_index must count exported genes.")
        if source is None and filled != 0:
            raise ValueError(f"{label}: the var index cannot fill names from itself.")
        renamed = name_set["renamed"]
        if not isinstance(renamed, list) or any(
            not isinstance(pair, list)
            or len(pair) != 2
            or not all(isinstance(name, str) and name for name in pair)
            for pair in renamed
        ):
            raise ValueError(f"{label}.renamed must be a list of [from, to] name pairs.")
        names = name_set["names"]
        if position == 0:
            if names is not None:
                raise ValueError(f"{label}.names must be null: the field names are this set.")
            resolved = list(default_names)
        else:
            if not isinstance(names, list) or len(names) != len(default_names):
                raise ValueError(f"{label}.names must hold one entry per gene.")
            if any(name is not None and (not isinstance(name, str) or not name) for name in names):
                raise ValueError(f"{label}.names must hold non-empty strings or null.")
            if any(name == default for name, default in zip(names, default_names, strict=True)):
                raise ValueError(
                    f"{label}.names must be null wherever a name equals the default name."
                )
            resolved = [
                default if name is None else name
                for name, default in zip(names, default_names, strict=True)
            ]
        _require_unique_identifiers(resolved, label=f"{label} gene")
    if len(set(sources)) != len(sources):
        raise ValueError("var manifest gene_name_sets must read each source once.")


def _declared_embedding_payload_paths(identity: dict) -> set[str]:
    """Re-derive every embedding payload path ``dataset_identity.json`` declares.

    The declared ``file`` of the item at position ``index`` must be exactly
    ``embeddings/{index}_{dimension}d.bin``, plus ``.gz`` when the export is
    compressed, so the declaration is checked against that rule before it is
    compared with the directory.
    """
    embeddings = identity["embeddings"]
    compression = identity["export_settings"]["compression"]
    suffix = ".gz" if compression is not None else ""
    declared: set[str] = set()
    for index, (embedding_id, item) in enumerate(embeddings["items"].items()):
        expected = f"embeddings/{index}_{item['dimension']}d.bin{suffix}"
        if item["file"] != expected:
            raise RuntimeError(
                f"Embedding {embedding_id!r} declares file {item['file']!r}; "
                f"the format requires {expected!r}."
            )
        declared.add(expected)
    return declared


_CODES_EXTENSION_BY_DTYPE = {"uint8": "u8", "uint16": "u16"}


def _declared_obs_payload_paths(manifest: dict) -> set[str]:
    """Re-derive every obs payload path the emitted manifest points a reader at."""
    schemas = manifest["_obsSchemas"]
    continuous_schema = schemas.get("continuous") or {}
    categorical_schema = schemas.get("categorical") or {}
    declared: set[str] = set()
    for field in manifest["_continuousFields"]:
        declared.add(
            _expand_payload_pattern(
                continuous_schema.get("pathPattern"),
                index=field[0],
                label="obs continuous pathPattern",
            )
        )
    for field in manifest["_categoricalFields"]:
        extension = _CODES_EXTENSION_BY_DTYPE.get(field[3])
        if extension is None:
            raise RuntimeError(
                f"Categorical obs field {field[1]!r} declares an unknown codes dtype {field[3]!r}."
            )
        declared.add(
            _expand_payload_pattern(
                categorical_schema.get("codesPathPattern"),
                index=field[0],
                extension=extension,
                label="obs categorical codesPathPattern",
            )
        )
        declared.add(
            _expand_payload_pattern(
                categorical_schema.get("outlierPathPattern"),
                index=field[0],
                label="obs categorical outlierPathPattern",
            )
        )
    return declared


def _declared_var_payload_paths(manifest: dict) -> set[str]:
    """Re-derive every var payload path the emitted manifest points a reader at."""
    schema = manifest["_varSchema"]
    return {
        _expand_payload_pattern(
            schema.get("pathPattern"),
            index=field[0],
            label="var pathPattern",
        )
        for field in manifest["fields"]
    }
