"""
What a person made in the Cellucid viewer, brought into an AnnData in one call.

The viewer sends every result it holds (``viewer-protocol.json``, kind
``result``): highlight groups, user-defined categorical fields (labels,
annotations, clusters), continuous fields and scores, gene sets, the filter's
visible cells and the cells it opened. :func:`pull` writes them onto
``obs`` -- the global ``obs`` of a MuData -- and the gene sets into
``uns['cellucid']['gene_sets']``, from a live viewer
(:class:`~cellucid.AnnDataViewer`) or from what the web app exports with no
notebook at all: obs-column CSV files, GMT files and ``.cellucid-session``
files.

Column names, one scheme (documented in the Python API reference, "Round
trip"):

- a highlight group: ``<highlight_prefix><group id>``
  (``cellucid_highlight__highlight_3``), ``bool``;
- a categorical or continuous field, a score: ``<field_prefix><key>``;
- the visible cells (``visibility=True``): ``cellucid_visible``.

Every column a pull writes is recorded in ``uns['cellucid']['pulled']`` with
where it came from (``pulled_from``: ``"viewer"`` or the file's name) and the
result it holds (``kind``, ``id``). A later pull of the same result writes the
same column again, whatever suffix it took; a viewer pull removes a column
only when the viewer no longer holds its result (``prune``); a column another
source wrote, or Cellucid did not write, is a conflict (``on_conflict``). The
column rules -- dtype, missing values, a subsample's unopened cells as
``<NA>`` -- are :mod:`._columns`', shared with
:func:`~cellucid.apply_cellucid_session_to_anndata`.
"""

from __future__ import annotations

import copy
import os
import re
import warnings
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .._anndata_files import is_mudata
from .._plural import plural
from ._columns import categorical_column, continuous_column, membership_column, write_columns
from ._records import _ColumnPlan

#: The ``uns['cellucid']`` key recording what :func:`pull` wrote.
PULLED_KEY = "pulled"
#: The ``uns['cellucid']`` key holding the gene sets :func:`pull` brought.
GENE_SETS_KEY = "gene_sets"
#: The column of the filter's visible cells.
VISIBLE_COLUMN = "cellucid_visible"
#: The first column of an obs-column CSV: the cells' obs names.
OBS_NAME_COLUMN = "obs_name"
#: The first cell of an obs-column CSV's second row, which declares each
#: column's kind.
KIND_ROW = "kind"
#: What a column of an obs-column CSV holds (its kind row).
CSV_KINDS = ("highlight", "categorical", "continuous")
#: Where a pull read a result: the live viewer (a file's results carry the
#: file's name).
VIEWER_SOURCE = "viewer"

_CONFLICT_RULES = frozenset({"error", "replace", "suffix"})


@dataclass(frozen=True)
class ViewerResult:
    """One result the viewer holds (the protocol's ``result``).

    ``cells`` are original positions (rows of ``obs``), ascending, or None for
    every cell; ``values`` are a field's codes or numbers, one per cell named
    (or per cell of the dataset), or None.
    """

    kind: str
    id: str
    revision: int
    n_obs: int
    meta: dict[str, Any]
    cells: np.ndarray | None
    values: np.ndarray | None
    #: The ``obs`` column it is written to, when it names one itself (a column
    #: of an obs-column CSV); None follows the naming scheme.
    column: str | None = None
    #: Where it was read: ``"viewer"``, or the name of the file it came from.
    source: str = VIEWER_SOURCE
    #: A highlight's known cells, when the source says which (an obs-column
    #: CSV lists them); None: the cells the viewer opened.
    known: np.ndarray | None = None

    @classmethod
    def from_event(cls, payload: Mapping[str, Any]) -> ViewerResult:
        """The result of one decoded ``result`` event payload."""
        return cls(
            kind=str(payload["kind"]),
            id=str(payload["id"]),
            revision=int(payload["revision"]),
            n_obs=int(payload["nObs"]),
            meta=copy.deepcopy(dict(payload["meta"])),
            cells=payload["cells"],
            values=payload["values"],
        )


@dataclass(frozen=True)
class Selection:
    """The cells of the viewer's last selection: original positions, their obs
    names (when the viewer's dataset has them in this process), and the tool
    that made it."""

    cells: np.ndarray
    obs_names: pd.Index | None
    source: str

    def __len__(self) -> int:
        return int(self.cells.size)


@dataclass(frozen=True)
class PullSummary:
    """What :func:`pull` wrote.

    ``columns`` maps each ``obs`` column written to what it holds (``kind``,
    ``id``, ``pulled_from`` and the result's label, key or page);
    ``removed`` lists the columns of an earlier viewer pull whose result the
    viewer no longer holds; ``gene_sets`` maps each gene set brought to its
    genes, and ``removed_gene_sets`` names the gene sets of an earlier viewer
    pull the viewer no longer holds; ``skipped`` names results left out and
    why.
    """

    columns: dict[str, dict[str, Any]] = field(default_factory=dict)
    removed: list[str] = field(default_factory=list)
    gene_sets: dict[str, list[str]] = field(default_factory=dict)
    removed_gene_sets: list[str] = field(default_factory=list)
    skipped: dict[str, str] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# GMT
# ---------------------------------------------------------------------------


_GMT_BREAKS = re.compile(r"[\t\r\n]+")


def _clean_gmt_text(value: object) -> str:
    """Text a GMT cell can hold: each run of tabs and line breaks one space,
    trimmed (the web app's ``formatGmt`` and cellucid-r's rule)."""
    return _GMT_BREAKS.sub(" ", str(value)).strip()


def unique_gene_set_name(name: str, taken: set[str]) -> str:
    """``name``, or ``"name (2)"``, ``"name (3)"``... -- the first not in
    ``taken`` -- added to ``taken``. The one rule for a gene set name used
    twice, in the web app, cellucid-python and cellucid-r alike."""
    candidate = name
    number = 2
    while candidate in taken:
        candidate = f"{name} ({number})"
        number += 1
    taken.add(candidate)
    return candidate


def write_gmt(
    gene_sets: Mapping[str, Sequence[str]],
    path: str | os.PathLike[str],
    *,
    descriptions: Mapping[str, str] | None = None,
) -> Path:
    """Write gene sets as a GMT file (MSigDB): one set per line,
    ``name<TAB>description<TAB>gene...``.

    Parameters
    ----------
    gene_sets
        Each set's name and its genes, in order.
    path
        The file to write.
    descriptions
        A description per set name; empty when not given.

    Returns
    -------
    Path
        The file written. Each run of tabs and line breaks in a name,
        description or gene is written as one space; a name that then reads
        as an earlier one's is written ``"name (2)"`` (as :func:`read_gmt`
        names a repeated name); each gene once. The Cellucid web app reads it
        (Gene signature and Dot plot: Import), as do :func:`read_gmt` and
        cellucid-r's ``cellucid_read_gmt()``.
    """
    if not isinstance(gene_sets, Mapping):
        raise TypeError("gene_sets must map each set name to its genes")
    descriptions = {} if descriptions is None else dict(descriptions)
    lines: list[str] = []
    taken: set[str] = set()
    for name, genes in gene_sets.items():
        clean_name = _clean_gmt_text(name)
        if not clean_name:
            raise ValueError("A gene set needs a name")
        if isinstance(genes, str):
            raise TypeError(f"Gene set {name!r} must list its genes, not one string")
        clean_genes = dict.fromkeys(gene for gene in (_clean_gmt_text(gene) for gene in genes) if gene)
        lines.append("\t".join([
            unique_gene_set_name(clean_name, taken), _clean_gmt_text(descriptions.get(name, "")), *clean_genes,
        ]))
    target = Path(path).expanduser()
    target.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")
    return target


def read_gmt(path: str | os.PathLike[str]) -> dict[str, list[str]]:
    """Read a GMT file into ``{set name: [genes]}``, in file order.

    One rule with the web app and cellucid-r: cells are tab-separated and
    trimmed; blank lines are skipped; a line without a name is ``"Set <n>"``
    (its place among the sets); a name used again is ``"name (2)"``,
    ``"name (3)"``...; each gene is kept once.
    """
    sets: dict[str, list[str]] = {}
    taken: set[str] = set()
    for line in Path(path).expanduser().read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        cells = line.split("\t")
        name = unique_gene_set_name(cells[0].strip() or f"Set {len(sets) + 1}", taken)
        sets[name] = list(dict.fromkeys(gene.strip() for gene in cells[2:] if gene.strip()))
    return sets


# ---------------------------------------------------------------------------
# Results -> columns
# ---------------------------------------------------------------------------


def _uns_record(target: Any) -> dict[str, Any]:
    uns = target.uns
    cellucid = uns.get("cellucid")
    if cellucid is None:
        return {}
    if not isinstance(cellucid, Mapping):
        raise TypeError("uns['cellucid'] must be a dictionary if present")
    return dict(cellucid)


def _require_conflict_rule(value: object) -> str:
    if value not in _CONFLICT_RULES:
        raise ValueError(f"on_conflict must be one of {sorted(_CONFLICT_RULES)}, not {value!r}")
    return str(value)


#: Which result a column holds: where it was pulled from, its kind and id.
_ResultKey = tuple[str, str, str]


def _key_of(info: Mapping[str, Any]) -> _ResultKey | None:
    """The result a ``pulled`` record names, or None for a record that names
    none (it is then never replaced or removed by a pull)."""
    pulled_from, kind, result_id = info.get("pulled_from"), info.get("kind"), info.get("id")
    if not all(isinstance(value, str) for value in (pulled_from, kind, result_id)):
        return None
    return (str(pulled_from), str(kind), str(result_id))


@dataclass
class _Naming:
    """Column names, unique among the pull's own.

    A result an earlier pull wrote takes the column it wrote then (``ours``,
    by result) -- its own name, or that name with the suffix it took (unless
    ``on_conflict="replace"`` asks for the name itself) -- so a pull run
    again writes the same column. Otherwise its name is used when no
    column has it, or when a pull of the same source wrote it; any other
    existing column is a conflict, settled by ``on_conflict``.
    """

    on_conflict: str
    existing: set[str]
    ours: dict[_ResultKey, str]
    written_by: dict[str, str]
    used: set[str] = field(default_factory=set)
    conflicts: list[str] = field(default_factory=list)

    def _suffixed(self, wanted: str) -> str:
        number = 2
        while f"{wanted}_{number}" in self.existing or f"{wanted}_{number}" in self.used:
            number += 1
        return f"{wanted}_{number}"

    def name(self, wanted: str, key: _ResultKey) -> str:
        earlier = self.ours.get(key)
        if earlier is not None and earlier not in self.used and (
            earlier == wanted
            or (self.on_conflict != "replace" and re.fullmatch(rf"{re.escape(wanted)}_[0-9]+", earlier) is not None)
        ):
            name = earlier
        elif wanted in self.used:
            name = self._suffixed(wanted)
        elif wanted not in self.existing or self.written_by.get(wanted) == key[0] or self.on_conflict == "replace":
            name = wanted
        elif self.on_conflict == "suffix":
            name = self._suffixed(wanted)
        else:
            self.conflicts.append(wanted)
            name = wanted
        self.used.add(name)
        return name


def _open_cells(results: Iterable[ViewerResult]) -> np.ndarray | None:
    for result in results:
        if result.kind == "subset":
            return result.cells
    return None


@dataclass(frozen=True)
class _Candidate:
    """One column a pull would write: the result it holds, the name it wants."""

    key: _ResultKey
    wanted: str
    values: pd.Series
    info: dict[str, Any]


def _result_candidates(
    target: Any,
    results: Sequence[ViewerResult],
    *,
    highlights: bool,
    fields: bool,
    scores: bool,
    visibility: bool,
    include_pushed: bool,
    highlight_prefix: str,
    field_prefix: str,
) -> tuple[list[_Candidate], dict[str, str]]:
    """The columns of ``results`` on ``target`` before naming, and the results
    skipped (with why)."""
    n_obs = int(target.n_obs)
    index = target.obs_names
    for result in results:
        if result.n_obs != n_obs:
            raise ValueError(
                f"The viewer shows a dataset of {result.n_obs} {plural(result.n_obs, 'cell')} and this object has {n_obs}: "
                "pull into the object the viewer was opened on, with its cells in the same order."
            )
    open_cells = _open_cells(results)
    candidates: list[_Candidate] = []
    skipped: dict[str, str] = {}
    for result in sorted(results, key=lambda item: (item.kind, item.revision)):
        key = (result.source, result.kind, result.id)
        meta = result.meta
        provenance = {"kind": result.kind, "id": result.id, "pulled_from": result.source}
        if result.kind == "highlight":
            if not highlights:
                continue
            wanted = result.column or f"{highlight_prefix}{result.id}"
            values = membership_column(
                n_obs, result.cells, open_cells=open_cells if result.known is None else result.known, index=index
            )
            info = {**provenance, "label": meta["label"], "tool": meta["tool"], "page": meta["pageName"],
                    "enabled": meta["enabled"]}
        elif result.kind in {"categorical", "continuous"}:
            origin = meta["origin"]
            if origin == "notebook" and not include_pushed:
                skipped[f"{result.kind}:{result.id}"] = "pushed from a notebook (include_pushed=True brings it back)"
                continue
            if origin == "analysis" and not scores:
                continue
            if origin != "analysis" and not fields:
                continue
            wanted = result.column or f"{field_prefix}{meta['key']}"
            assert result.values is not None
            if result.kind == "categorical":
                values = categorical_column(
                    result.values,
                    meta["categories"],
                    cells=result.cells,
                    n_obs=n_obs,
                    missing_code=int(np.iinfo(result.values.dtype).max),
                    index=index,
                )
            else:
                values = continuous_column(result.values, cells=result.cells, n_obs=n_obs, index=index)
            info = {**provenance, "key": meta["key"], "origin": origin}
        elif result.kind == "visibility":
            if not visibility:
                continue
            wanted = VISIBLE_COLUMN
            values = membership_column(n_obs, result.cells, open_cells=open_cells, index=index)
            info = provenance
        else:
            continue
        candidates.append(_Candidate(key=key, wanted=wanted, values=values, info=info))
    return candidates, skipped


def _name_candidates(
    target: Any, candidates: Sequence[_Candidate], *, on_conflict: str
) -> tuple[list[_ColumnPlan], dict[str, Any]]:
    """The plans of ``candidates``, named (see :class:`_Naming`), and the
    ``pulled`` record they are named against."""
    record = _uns_record(target)
    previous = record.get(PULLED_KEY) or {}
    if not isinstance(previous, Mapping):
        raise TypeError("uns['cellucid']['pulled'] must be a dictionary")
    columns = set(map(str, target.obs.columns))
    ours: dict[_ResultKey, str] = {}
    written_by: dict[str, str] = {}
    for column, info in previous.items():
        key = _key_of(info) if isinstance(info, Mapping) else None
        if key is None or column not in columns:
            continue
        ours[key] = str(column)
        written_by[str(column)] = key[0]
    naming = _Naming(on_conflict=on_conflict, existing=columns, ours=ours, written_by=written_by)
    plans = [
        _ColumnPlan(name=naming.name(candidate.wanted, candidate.key), values=candidate.values, metadata=candidate.info)
        for candidate in candidates
    ]
    if naming.conflicts:
        raise ValueError(
            "These obs columns exist already and no pull from the same source wrote them: "
            + ", ".join(repr(name) for name in naming.conflicts)
            + ". Pass on_conflict='replace' to overwrite them, on_conflict='suffix' to write "
            "beside them, or another highlight_prefix/field_prefix."
        )
    return plans, dict(previous)


def _pruned(
    target: Any, previous: Mapping[str, Any], plans: Sequence[_ColumnPlan], held: set[_ResultKey]
) -> list[str]:
    """The columns of an earlier viewer pull to remove: those whose result the
    viewer no longer holds, and those of a result this pull wrote under
    another name. A column a file pull wrote is never removed."""
    written = {plan.name for plan in plans}
    written_keys = {_key_of(plan.metadata) for plan in plans}
    remove = []
    for column, info in previous.items():
        key = _key_of(info) if isinstance(info, Mapping) else None
        if key is None or key[0] != VIEWER_SOURCE or column in written or column not in target.obs.columns:
            continue
        if key not in held or key in written_keys:
            remove.append(str(column))
    return sorted(remove)


def _gene_sets_of(results: Sequence[ViewerResult]) -> list[tuple[_ResultKey, str, dict[str, Any]]]:
    """Each gene set: the result it is, a name unique among them
    (:func:`unique_gene_set_name`), and its entry: its genes by identifier
    (the var index: ``adata.var_names``), and, from a viewer whose data
    mapping names the genes by a var column, the name it shows for each
    (``symbols``)."""
    sets = []
    taken: set[str] = set()
    for result in results:
        if result.kind != "geneSet":
            continue
        meta = result.meta
        sets.append((
            (result.source, result.kind, result.id),
            unique_gene_set_name(meta["name"], taken),
            {"genes": list(meta["genes"]), "description": meta["description"], "origin": meta["origin"],
             **({"symbols": list(meta["symbols"])} if meta.get("symbols") is not None else {})},
        ))
    return sets


def _merge_gene_sets(
    existing: Mapping[str, Any],
    brought: Sequence[tuple[_ResultKey, str, dict[str, Any]]],
    *,
    held: set[_ResultKey] | None,
) -> tuple[dict[str, Any], list[str], dict[str, list[str]]]:
    """``uns['cellucid']['gene_sets']`` with ``brought``: a set brought again
    replaces its earlier entry; with ``held`` (a viewer pull that prunes), a
    set of an earlier viewer pull the viewer no longer holds goes. Names are
    unique (:func:`unique_gene_set_name`). Returns the sets, the names
    removed, and each set brought by its name."""
    brought_keys = {key for key, _, _ in brought}
    merged: dict[str, Any] = {}
    removed: list[str] = []
    for name, entry in existing.items():
        key = _key_of({**entry, "kind": "geneSet"}) if isinstance(entry, Mapping) else None
        if key is not None and key in brought_keys:
            continue
        if key is not None and held is not None and key[0] == VIEWER_SOURCE and key not in held:
            removed.append(str(name))
            continue
        merged[str(name)] = entry
    taken = set(merged)
    names: dict[str, list[str]] = {}
    for key, name, entry in brought:
        final = unique_gene_set_name(name, taken)
        merged[final] = {**entry, "pulled_from": key[0], "id": key[2]}
        names[final] = list(entry["genes"])
    return merged, sorted(removed), names


def _write_pull(
    target: Any,
    plans: Sequence[_ColumnPlan],
    previous: Mapping[str, Any],
    remove: Sequence[str],
    *,
    gene_sets: Sequence[tuple[_ResultKey, str, dict[str, Any]]] | None,
    held_gene_sets: set[_ResultKey] | None,
    skipped: dict[str, str],
) -> PullSummary:
    """Write the named columns and the gene sets (None: leave them) onto
    ``target`` at once, and record what was written."""
    write_columns(target, plans, remove=list(remove))
    record = _uns_record(target)
    pulled = {
        str(name): info for name, info in previous.items() if name not in remove and name in target.obs.columns
    }
    pulled.update({plan.name: copy.deepcopy(plan.metadata) for plan in plans})
    record[PULLED_KEY] = pulled
    removed_sets: list[str] = []
    brought: dict[str, list[str]] = {}
    if gene_sets is not None:
        record[GENE_SETS_KEY], removed_sets, brought = _merge_gene_sets(
            dict(record.get(GENE_SETS_KEY) or {}), gene_sets, held=held_gene_sets
        )
    target.uns["cellucid"] = record
    return PullSummary(
        columns={plan.name: copy.deepcopy(plan.metadata) for plan in plans},
        removed=list(remove),
        gene_sets=brought,
        removed_gene_sets=removed_sets,
        skipped=skipped,
    )


def _require_pull_options(target: Any, flags: Mapping[str, object], prefixes: Mapping[str, object], on_conflict: object) -> str:
    for label, flag in flags.items():
        if type(flag) is not bool:
            raise TypeError(f"{label} must be exactly True or False")
    for label, prefix in prefixes.items():
        if not isinstance(prefix, str):
            raise TypeError(f"{label} must be a string")
    rule = _require_conflict_rule(on_conflict)
    if not hasattr(target, "obs") or not hasattr(target, "uns") or not hasattr(target, "n_obs"):
        raise TypeError("pull writes into an AnnData or a MuData")
    if getattr(target, "isbacked", False) and not is_mudata(target):
        raise ValueError("A backed AnnData cannot take new obs columns; load it with adata.to_memory() first")
    return rule


def apply_results(
    target: Any,
    results: Sequence[ViewerResult],
    *,
    highlights: bool = True,
    fields: bool = True,
    scores: bool = True,
    visibility: bool = False,
    gene_sets: bool = True,
    include_pushed: bool = False,
    highlight_prefix: str = "cellucid_highlight__",
    field_prefix: str = "",
    on_conflict: str = "error",
    prune: bool = True,
) -> PullSummary:
    """Write a viewer's ``results`` onto ``target`` in place (see :func:`pull`)."""
    on_conflict = _require_pull_options(
        target,
        {"highlights": highlights, "fields": fields, "scores": scores, "visibility": visibility,
         "gene_sets": gene_sets, "include_pushed": include_pushed, "prune": prune},
        {"highlight_prefix": highlight_prefix, "field_prefix": field_prefix},
        on_conflict,
    )
    candidates, skipped = _result_candidates(
        target, results, highlights=highlights, fields=fields, scores=scores, visibility=visibility,
        include_pushed=include_pushed, highlight_prefix=highlight_prefix, field_prefix=field_prefix,
    )
    plans, previous = _name_candidates(target, candidates, on_conflict=on_conflict)
    held = {(result.source, result.kind, result.id) for result in results}
    remove = _pruned(target, previous, plans, held) if prune else []
    return _write_pull(
        target, plans, previous, remove,
        gene_sets=_gene_sets_of(results) if gene_sets else None,
        held_gene_sets=held if prune else None,
        skipped=skipped,
    )


# ---------------------------------------------------------------------------
# The web app's exports
# ---------------------------------------------------------------------------


def _csv_column_result(
    name: str,
    kind: str,
    column: pd.Series,
    positions: np.ndarray,
    n_obs: int,
    *,
    source: str,
    revision: int,
) -> ViewerResult:
    """One CSV column as the result its kind row declares: ``highlight``
    (``true``/``false``; empty is unknown), ``categorical`` (labels; empty is
    no category) or ``continuous`` (numbers; empty is no value). A cell the
    file does not list is unknown, or has no value."""
    text = column.astype(str).str.strip()
    present = (text != "").to_numpy()
    order = np.argsort(positions, kind="stable")
    if kind == "highlight":
        lowered = text.str.lower()
        unreadable = present & ~lowered.isin(["true", "false"]).to_numpy()
        if unreadable.any():
            raise ValueError(
                f"{source}: the highlight column {name!r} holds {text.iloc[int(np.flatnonzero(unreadable)[0])]!r}; "
                "a highlight column holds true, false or nothing"
            )
        members = positions[present & (lowered == "true").to_numpy()]
        return ViewerResult(
            kind="highlight", id=name, revision=revision, n_obs=n_obs,
            meta={"label": name, "tool": "annotation", "pageId": "csv", "pageName": source, "color": "#000000",
                  "enabled": True},
            cells=np.sort(members).astype(np.uint32), values=None, column=name, source=source,
            known=np.sort(positions[present]).astype(np.int64),
        )
    if kind == "continuous":
        numbers = pd.to_numeric(text.where(present, None), errors="coerce").to_numpy(dtype=np.float64)
        unreadable = present & np.isnan(numbers)
        if unreadable.any() or np.isinf(numbers).any():
            bad = int(np.flatnonzero(unreadable | np.isinf(numbers))[0])
            raise ValueError(
                f"{source}: the continuous column {name!r} holds {text.iloc[bad]!r}; "
                "a continuous column holds finite numbers or nothing"
            )
        return ViewerResult(
            kind="continuous", id=name, revision=revision, n_obs=n_obs, meta={"key": name, "origin": "user"},
            cells=positions[order].astype(np.uint32), values=numbers[order], column=name, source=source,
        )
    categories = list(dict.fromkeys(text[present].tolist()))
    lookup = {category: code for code, category in enumerate(categories)}
    dtype = _code_dtype(len(categories))
    missing = np.iinfo(dtype).max
    codes = np.array([lookup[value] if keep else missing for value, keep in zip(text, present, strict=True)],
                     dtype=dtype)
    return ViewerResult(
        kind="categorical", id=name, revision=revision, n_obs=n_obs,
        meta={"key": name, "categories": categories, "origin": "user"},
        cells=positions[order].astype(np.uint32), values=codes[order], column=name, source=source,
    )


def read_obs_csv(path: str | os.PathLike[str], obs_names: pd.Index) -> list[ViewerResult]:
    """Read an obs-column CSV the web app exports: a header row (``obs_name``,
    then each column's name), a kind row (``kind``, then each column's kind:
    ``highlight``, ``categorical`` or ``continuous``), then a row per cell.
    A cell the file does not list is unknown; a name this object does not
    have is refused."""
    path = Path(path).expanduser()
    with path.open(encoding="utf-8", newline="") as handle:
        frame = pd.read_csv(handle, dtype=str, keep_default_na=False, header=None)
    if frame.shape[0] < 2 or frame.shape[1] < 1 or frame.iloc[0, 0] != OBS_NAME_COLUMN or frame.iloc[1, 0] != KIND_ROW:
        raise ValueError(
            f"{path.name}: an obs-column CSV starts with a header row whose first cell is {OBS_NAME_COLUMN!r} "
            f"and a kind row whose first cell is {KIND_ROW!r} (then each column's kind: {', '.join(CSV_KINDS)})"
        )
    header = [str(value) for value in frame.iloc[0]]
    kinds = [str(value) for value in frame.iloc[1]]
    for name, kind in zip(header[1:], kinds[1:], strict=True):
        if kind not in CSV_KINDS:
            raise ValueError(f"{path.name}: the column {name!r} is of kind {kind!r}, not one of {', '.join(CSV_KINDS)}")
    if len(set(header)) != len(header):
        raise ValueError(f"{path.name}: a column is named twice")
    rows = frame.iloc[2:].reset_index(drop=True)
    names = rows[0]
    if names.duplicated().any():
        raise ValueError(f"{path.name}: the cell {names[names.duplicated()].iloc[0]!r} is listed twice")
    if not obs_names.is_unique:
        raise ValueError("This object's obs names are not unique; make them unique (obs_names_make_unique()) first")
    positions = obs_names.get_indexer(names)
    if (positions < 0).any():
        unknown = names[positions < 0]
        raise ValueError(
            f"{path.name}: {len(unknown)} of its cells are not in this object (first: {unknown.iloc[0]!r})"
        )
    return [
        _csv_column_result(
            header[number], kinds[number], rows[number], positions.astype(np.int64), len(obs_names),
            source=path.name, revision=number,
        )
        for number in range(1, len(header))
    ]


def results_of_files(paths: Sequence[Path], target: Any) -> tuple[list[ViewerResult], list[Path]]:
    """The results of the web app's exports among ``paths`` (CSV and GMT), and
    the session files among them, which are applied as sessions."""
    results: list[ViewerResult] = []
    sessions: list[Path] = []
    revision = 1
    for path in paths:
        suffix = path.name.lower()
        if suffix.endswith(".csv"):
            for result in read_obs_csv(path, target.obs_names):
                results.append(result)
        elif suffix.endswith(".gmt"):
            for name, genes in read_gmt(path).items():
                results.append(ViewerResult(
                    kind="geneSet", id=name, revision=revision, n_obs=int(target.n_obs),
                    meta={"name": name, "description": "", "genes": genes, "origin": "user"},
                    cells=np.empty(0, dtype=np.uint32), values=None, source=path.name,
                ))
                revision += 1
        elif suffix.endswith(".cellucid-session"):
            sessions.append(path)
        else:
            raise ValueError(
                f"{path.name}: pull reads the web app's exports -- .csv (obs columns), .gmt (gene sets) "
                "and .cellucid-session files"
            )
    return results, sessions


def export_paths(source: str | os.PathLike[str] | Sequence[str | os.PathLike[str]]) -> list[Path]:
    """The files of ``source``: one file, a folder's exports, or a list of them."""
    items = [source] if isinstance(source, str | os.PathLike) else list(source)
    paths: list[Path] = []
    for item in items:
        path = Path(item).expanduser()
        if path.is_dir():
            paths.extend(
                sorted(
                    child for child in path.iterdir()
                    if child.is_file() and child.name.lower().endswith((".csv", ".gmt", ".cellucid-session"))
                )
            )
        elif path.is_file():
            paths.append(path)
        else:
            raise FileNotFoundError(f"{path} does not exist")
    if not paths:
        raise ValueError(f"{source} holds no .csv, .gmt or .cellucid-session file to pull")
    return paths




# ---------------------------------------------------------------------------
# The one call
# ---------------------------------------------------------------------------


def pull(
    adata: Any,
    source: Any,
    *,
    expected_dataset_id: str | None = None,
    highlights: bool = True,
    fields: bool = True,
    scores: bool = True,
    visibility: bool = False,
    gene_sets: bool = True,
    include_pushed: bool = False,
    highlight_prefix: str = "cellucid_highlight__",
    field_prefix: str = "",
    on_conflict: str = "error",
    prune: bool = True,
) -> PullSummary:
    """Write what was made in Cellucid onto ``adata`` in one call, in place.

    Every highlight group (a ``bool`` column each), every user-defined
    categorical field -- labels, annotations, clusters -- (``category``),
    every continuous field and score (``float32``) go to ``adata.obs`` (a
    MuData's global ``obs``), and every gene set to
    ``adata.uns['cellucid']['gene_sets']``: its genes by identifier -- the
    var index, ``adata.var_names`` -- and, from a viewer whose data mapping
    names the genes by a var column (such as ``feature_name``), the names it
    shows for them (``symbols``).

    Parameters
    ----------
    adata : anndata.AnnData or mudata.MuData
        The object the viewer shows, with its cells in the same order. From
        a viewer, an object holding the viewer's cells in another order is
        refused (nothing is written), as is one read from a Cellucid-ready
        file in its stored order (``obs['cellucid_original_position']``;
        :func:`cellucid.restore_original_order` puts it back).
    source : AnnDataViewer, path, or list of paths
        A viewer (a notebook widget or a browser-tab viewer): what it holds,
        as it last sent it -- every result as it was made, and all of them
        whenever it opened its dataset, so a cell run after the work pulls it
        all, with no waiting -- or what the web app exports with no
        notebook: obs-column ``.csv`` files (a header row, ``obs_name``
        then each column's name; a kind row, ``kind`` then each column's
        kind, ``highlight``, ``categorical`` or ``continuous``; a row per
        cell; each column written as named), ``.gmt`` gene sets and
        ``.cellucid-session`` files -- one file, a list, or a folder of them.
    expected_dataset_id : str, optional
        For a session file: the dataset it must have been saved on. Defaults
        to ``adata.uns['cellucid']['id']``.
    highlights, fields, scores, visibility, gene_sets : bool
        Which results to write: highlight groups; user-defined fields; scores
        an analysis made (a gene signature's); the filter's visible cells
        (``cellucid_visible``, off by default); gene sets.
    include_pushed : bool
        Also write the fields this notebook pushed (``viewer.push``), which
        ``adata`` holds already.
    highlight_prefix, field_prefix : str
        Each highlight's column is ``<highlight_prefix><group id>``; each
        field's and score's ``<field_prefix><key>``.
    on_conflict : {"error", "replace", "suffix"}
        What a column name that ``adata.obs`` has already does, when no pull
        from the same source (the viewer, or the same file) wrote it: raise
        (nothing is written), overwrite it, or write ``<name>_2`` beside it.
        A result an earlier pull wrote is written to the same column again,
        suffix and all.
    prune : bool
        For a viewer: remove the columns, and the gene sets, an earlier
        viewer pull wrote whose result the viewer no longer holds (a
        highlight deleted since). A result left out of this pull
        (``highlights=False``, a pushed field) keeps its column, and files
        never remove anything.

    Returns
    -------
    PullSummary
        The columns written (``columns``), removed (``removed``), the gene
        sets brought (``gene_sets``) and removed (``removed_gene_sets``),
        and the results left out (``skipped``).

    Notes
    -----
    When the viewer opened a subsample, a cell it did not open is ``<NA>`` in
    a highlight column and missing in a field's: unknown, not "not
    selected"; a cell an obs-column CSV does not list is the same.
    ``adata.uns['cellucid']['pulled']`` records each column written, with
    where it came from (``pulled_from``) and the result it holds (``kind``,
    ``id``); each gene set in ``uns['cellucid']['gene_sets']`` the same.
    """
    from ..jupyter._viewer import AnnDataViewer

    if isinstance(source, AnnDataViewer):
        source._require_viewer_order(adata)
        return apply_results(
            adata,
            source.results,
            highlights=highlights,
            fields=fields,
            scores=scores,
            visibility=visibility,
            gene_sets=gene_sets,
            include_pushed=include_pushed,
            highlight_prefix=highlight_prefix,
            field_prefix=field_prefix,
            on_conflict=on_conflict,
            prune=prune,
        )
    from ._apply import _plan_session

    on_conflict = _require_pull_options(
        adata,
        {"highlights": highlights, "fields": fields, "scores": scores, "visibility": visibility,
         "gene_sets": gene_sets, "include_pushed": include_pushed, "prune": prune},
        {"highlight_prefix": highlight_prefix, "field_prefix": field_prefix},
        on_conflict,
    )
    paths = export_paths(source)
    results, sessions = results_of_files(paths, adata)
    candidates, skipped = _result_candidates(
        adata, results, highlights=highlights, fields=fields, scores=scores, visibility=visibility,
        include_pushed=include_pushed, highlight_prefix=highlight_prefix, field_prefix=field_prefix,
    )
    sets = _gene_sets_of(results)
    record = _uns_record(adata)
    for path in sessions:
        dataset_id = expected_dataset_id if expected_dataset_id is not None else record.get("id")
        if not isinstance(dataset_id, str) or not dataset_id:
            raise ValueError(
                f"{path.name}: pass expected_dataset_id= (the dataset id the session was saved on); "
                "this object names none in uns['cellucid']['id']"
            )
        # Planned against no existing column: its names are then settled with
        # every other column of this pull (on_conflict, reuse).
        session = _plan_session(
            path,
            adata,
            expected_dataset_id=dataset_id,
            add_highlights=highlights,
            highlights_prefix=highlight_prefix,
            add_user_defined_fields=fields,
            user_defined_prefix=field_prefix,
            include_deleted_user_defined_fields=False,
            existing_columns=set(),
        )
        for plan in session.plans:
            result_id = str(plan.metadata.get("group_id") or plan.metadata.get("field_id"))
            info = {**plan.metadata, "id": result_id, "pulled_from": path.name}
            candidates.append(_Candidate(key=(path.name, str(plan.metadata["kind"]), result_id), wanted=plan.name,
                                         values=plan.values, info=info))
        taken = {name for _, name, _ in sets}
        for entry in session.gene_sets:
            sets.append((
                (path.name, "geneSet", entry["name"]),
                unique_gene_set_name(entry["name"], taken),
                {"genes": list(entry["genes"]), "description": entry["description"], "origin": "notebook"},
            ))
    plans, previous = _name_candidates(adata, candidates, on_conflict=on_conflict)
    return _write_pull(
        adata, plans, previous, [], gene_sets=sets if gene_sets else None, held_gene_sets=None, skipped=skipped
    )


# ---------------------------------------------------------------------------
# The other way: a column and the dataset's obs names
# ---------------------------------------------------------------------------


def read_obs_names(data: Any) -> pd.Index | None:
    """The obs names of what a viewer shows, by original position (the order
    the viewer names cells in): an AnnData or MuData object's, or those an
    AnnData or MuData file holds (read alone, nothing else); None for a format
    only the browser reads (``.loom``, 10x, MatrixMarket, R objects).

    A file in priority order (``uns['cellucid']['priority']``, as
    ``prepare()`` writes it) stores its rows in another order: its names are
    put back in original order, as the viewer reads it -- when the record
    describes the file's cells (:func:`cellucid.prepare_data._priority.stored_order`,
    the rule the viewer reads it by); a record that does not is reported and
    the names are taken as stored. An object's names follow the same rule:
    one read from such a file and still in its stored order is named in its
    original order, as it is served (:mod:`cellucid._served_store`)."""
    if hasattr(data, "obs_names"):
        from ..prepare_data._priority import stored_order_of

        names = pd.Index(data.obs_names)
        positions = stored_order_of(data).positions if hasattr(data, "uns") and hasattr(data, "obs") else None
        if positions is None:
            return names
        named = np.empty(len(names), dtype=object)
        named[positions] = np.asarray(names, dtype=object)
        return pd.Index(named, dtype=object)
    from .._anndata_files import (
        ANNDATA_KINDS,
        MUDATA_KINDS,
        _open_root,
        _text,
        classify_anndata_path,
    )

    store = classify_anndata_path(data)
    if store.kind not in ANNDATA_KINDS | MUDATA_KINDS:
        return None
    from anndata.io import read_elem

    from ..prepare_data._priority import ORIGINAL_POSITION_COLUMN, stored_order

    with _open_root(store) as root:
        obs = root["obs"]
        names = np.asarray(read_elem(obs[_text(obs.attrs.get("_index", "_index"))]), dtype=object)
        cellucid_uns = root["uns"]["cellucid"] if "uns" in root and "cellucid" in root["uns"] else None
        record = read_elem(cellucid_uns["priority"]) if cellucid_uns is not None and "priority" in cellucid_uns else None
        column = read_elem(obs[ORIGINAL_POSITION_COLUMN]) if ORIGINAL_POSITION_COLUMN in obs else None
        order = stored_order(record, names.size, column)
        if order.notice is not None:
            warnings.warn(f"{data}: {order.notice}", UserWarning, stacklevel=2)
        if order.positions is not None:
            in_original_order = np.empty_like(names)
            in_original_order[order.positions] = names
            names = in_original_order
        return pd.Index(names, dtype=object)


def _code_dtype(count: int) -> np.dtype:
    """The narrowest code type whose largest value is free for "no category"."""
    for dtype in (np.uint8, np.uint16, np.uint32):
        if count < np.iinfo(dtype).max:
            return np.dtype(dtype)
    raise ValueError(f"{count} {plural(count, 'category', 'categories')} are more than a field holds")


def field_message(values: Any, key: str, *, n_obs: int, obs_names: pd.Index | None) -> dict[str, Any]:
    """The ``setField`` command of one column, aligned to the viewer's cells.

    A Series is aligned by its index: the viewer's obs names (every cell, or
    some, the rest without a value), or a ``RangeIndex`` of every cell in
    order. An array or list is one value per cell, in order. Categorical,
    string, object and bool values make a categorical (a bool's labels are
    ``"false"`` and ``"true"``, as cellucid-r pushes them); numbers a
    continuous field.
    """
    if isinstance(values, pd.Series):
        series = values
    elif isinstance(values, np.ndarray | list | tuple | pd.Index | pd.Categorical):
        series = pd.Series(values)
        if len(series) != n_obs:
            raise ValueError(
                f"The column has {len(series)} {plural(len(series), 'value')} and the viewer's dataset "
                f"{n_obs} {plural(n_obs, 'cell')}"
            )
    else:
        raise TypeError("push takes a pandas Series (e.g. adata.obs['leiden']), an array or a list")
    cells: np.ndarray | None = None
    index = series.index
    if obs_names is not None and index.equals(obs_names):
        pass
    elif obs_names is not None and not isinstance(index, pd.RangeIndex):
        if not obs_names.is_unique:
            raise ValueError("The viewer's obs names are not unique: push an array in obs order instead")
        positions = obs_names.get_indexer(index)
        if (positions < 0).any():
            raise ValueError(
                f"{int((positions < 0).sum())} of the column's index names are not cells of the viewer's dataset "
                f"(first: {index[positions < 0][0]!r})"
            )
        if pd.Index(positions).has_duplicates:
            raise ValueError("The column names a cell twice")
        order = np.argsort(positions, kind="stable")
        cells = positions[order].astype(np.int64)
        series = series.iloc[order]
        if cells.size == n_obs:
            cells = None
    elif len(series) != n_obs or (isinstance(index, pd.RangeIndex) and not index.equals(pd.RangeIndex(n_obs))):
        raise ValueError(
            f"The column has {len(series)} {plural(len(series), 'value')} for the viewer's {n_obs} "
            f"{plural(n_obs, 'cell')} and an index that names none "
            "of them; push a Series indexed by obs names, or one value per cell in order"
        )
    dtype = series.dtype
    if dtype.kind == "b" or isinstance(dtype, pd.BooleanDtype):
        # A boolean column's labels are "false" and "true", in every language.
        labels = series.map({False: "false", True: "true"})
        present = set(labels.dropna())
        series = labels.astype(pd.CategoricalDtype([label for label in ("false", "true") if label in present]))
        dtype = series.dtype
    if isinstance(dtype, pd.CategoricalDtype) or dtype.kind in "OUS" or pd.api.types.is_string_dtype(dtype):
        categorical = series if isinstance(dtype, pd.CategoricalDtype) else series.astype("category")
        categories = [str(category) for category in categorical.cat.categories]
        if len(set(categories)) != len(categories):
            raise ValueError("Two categories read as the same text; rename them before pushing")
        if not categories:
            raise ValueError("The column has no value to push")
        codes = categorical.cat.codes.to_numpy().astype(np.int64)
        code_dtype = _code_dtype(len(categories))
        out = np.where(codes < 0, np.iinfo(code_dtype).max, codes).astype(code_dtype)
        kind, field_categories, field_values = "categorical", categories, out
    elif dtype.kind in "iuf":
        numbers = series.to_numpy(dtype=np.float64, na_value=np.nan)
        if bool(np.isinf(numbers).any()):
            raise ValueError("The column holds an infinity; the viewer takes finite numbers (NaN for no value)")
        if bool((np.abs(numbers[np.isfinite(numbers)]) > np.finfo(np.float32).max).any()):
            raise ValueError(
                "The column holds a value beyond the float32 range the viewer takes (about 3.4e38); "
                "scale it first, or push its logarithm"
            )
        kind, field_categories, field_values = "continuous", None, numbers.astype(np.float32)
    else:
        raise TypeError(f"Values of dtype {dtype} cannot be pushed: push labels or numbers")
    return {
        "type": "setField",
        "requestId": "",
        "key": key,
        "kind": kind,
        "categories": field_categories,
        "nObs": n_obs,
        "cells": cells,
        "values": field_values,
    }
