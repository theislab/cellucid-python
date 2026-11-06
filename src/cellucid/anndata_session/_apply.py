"""
Applying one step of a session to an AnnData or a MuData: the plan, then the change.

:func:`apply_cellucid_session_to_anndata` (and
:meth:`cellucid.CellucidSession.apply`) reads a session file, takes the step
asked for (by tracker, branch and step), checks that the object is the
dataset the session was saved on, and plans every change the step's state
holds as rows (:class:`ApplyPlan`, the history model's ``ApplyPlan``): one
per renamed column, renamed or merged category, deleted field, user-defined
field, gene set, highlight group, and presentation domain the step changed.
A dry run returns the plan; :meth:`ApplyPlan.apply` (or ``dry_run=False``)
writes the rows to write, all at once, and records what it did in
``uns['cellucid']['history_applied']`` (:class:`ApplyReport`).

The rules are the same in cellucid-python, cellucid-r and the web app's
"Save AnnData with this state…" (``app/history/exports/anndata-copy.js``),
compared with one shared expected plan (``history_golden.expected.json``,
``apply.sequence``):

- **Rows**, in domain order. ``field-overlays``: renamed columns, then
  renamed categories, then deleted fields, then user-defined fields;
  ``gene-sets``: each set; ``highlights``: each group of every page; every
  other domain whose value at the step differs from the tracker's opening is
  one ``presentation`` row, not applied (filters, colours, views, the camera,
  analyses are how the data was shown, not the data).
- **Statuses**: ``write`` (the object lacks it), ``unchanged`` (the object
  holds exactly it already: applying again is a no-op, never a failure),
  ``replace`` (a column or gene set an earlier apply wrote, holding other
  values, with ``overwrite=True``), ``not-applied`` (with the reason: a
  column the person made is never replaced, a column or category the object
  lacks, a deletion, a presentation domain).
- **Provenance**: one record per apply that wrote something, appended to
  ``uns['cellucid']['history_applied']`` (JSON text: a list of the model's
  ``Provenance``): the session's SHA-256, the tracker, branch and step, when,
  the package and its version, the columns it owns (which a later apply
  replaces only with ``overwrite=True``), what it renamed, what it did not
  apply, and, with ``keep_session=True``, the session file itself (base64),
  so the history can be rebuilt from the object.
- **Gene sets** go to ``uns['cellucid']['gene_sets']`` by name (origin
  ``session``) and, one ``bool`` column each, to the ``var`` of the genes the
  data mapping reads: ``var`` when it reads ``X`` or a layer, the modality's
  ``var`` of a MuData (``mdata.mod[modality].var``; the global ``mdata.var``
  is not written). The columns are named over every set of the step, in its
  order: ``cellucid_gene_set__<name>``, the highlight columns' rule, with
  ``_2``, ``_3``... on a clash (``gene_set_columns``). Only a set written
  or replaced gets its column, ``True`` on each row whose gene (as the
  mapping names genes) the set lists, a gene the dataset lacks left out and
  warned: a new column goes after ``var``'s columns, an earlier apply's
  (``replace``) is replaced in place, and a column of that name beside a set
  no apply wrote is the dataset's own, kept and warned. A mapping that reads
  ``raw`` or no expression has no ``var`` to mark: the sets stay in ``uns``
  only, warned. A set the step no longer holds has no row: its record and
  its column stay. ``expected.apply.sequence[].var_columns`` pins the
  columns. The plan's rows and the provenance's ``columns`` (obs columns)
  are the same as without them.
- **The priority-order column** ``obs['cellucid_original_position']`` (each
  row's original position in a file stored in priority order) is never
  written, renamed or replaced: a plan that would, or that renames another
  column onto it, is refused before anything is written, dry runs included.
"""

from __future__ import annotations

import base64
import copy
import hashlib
import json
import warnings
from collections.abc import Mapping
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from .. import _history_model as model
from .. import _history_rules as _rules
from .._anndata_files import is_mudata
from .._plural import plural
from ..data_mapping import (
    MODALITY_KEY,
    mapped_var_count,
    require_mapping_matches_anndata,
    require_mapping_modality,
    resolve_gene_names,
)
from ..history import DOMAIN_NAMES, Step, Tracker
from ..session import CellucidSession, package_version
from ._cells import known_cells
from ._columns import write_columns
from ._fields import _plan_user_fields
from ._highlights import _plan_highlights, column_part, unique_column
from ._labels import _plan_labels
from ._mudata import MuDataSessionView
from ._primitives import _require_exact_keys, _require_list, _require_text
from ._records import _Candidate, _ColumnPlan
from ._schema import (
    CELLS_DOMAIN,
    FIELD_OVERLAYS_DOMAIN,
    GENE_SET_KEYS,
    GENE_SETS_DOMAIN,
    GENE_SETS_KEYS,
    HIGHLIGHTS_DOMAIN,
    OVERLAY_KEYS,
)

#: The ``uns['cellucid']`` key every apply appends its record to (JSON text).
HISTORY_APPLIED_KEY = "history_applied"
#: Where gene sets go (``uns['cellucid']['gene_sets']``), as a pull writes them.
GENE_SETS_KEY = "gene_sets"
#: The domains an apply writes; every other domain is presentation.
DATA_DOMAINS = (CELLS_DOMAIN, FIELD_OVERLAYS_DOMAIN, GENE_SETS_DOMAIN, HIGHLIGHTS_DOMAIN)
#: The attributes a plan and a report keep beside their fields.
_CONTEXT = "_context"
_ADATA = "_adata"
#: Why a presentation domain is not applied.
PRESENTATION_NOTE = "Presentation: not applied to the object."
#: The start of every gene set's ``var`` column.
GENE_SET_PREFIX = "cellucid_gene_set__"
#: Why the gene sets an apply writes have no ``var`` columns: the mapping reads ``raw.var`` or no genes.
GENE_SETS_UNS_ONLY_NOTE = (
    "The gene sets are kept in uns['cellucid']['gene_sets'] only: the data mapping reads its genes from another "
    "matrix than var, so they have no var columns."
)
#: Why an apply refuses a plan that would write, rename or replace the priority-order column.
ORIGINAL_POSITION_MESSAGE = (
    "obs['cellucid_original_position'] holds each row's original position in a file stored in priority order, "
    "so an apply never writes, renames or replaces it."
)


def owned_note(name: str, step: str) -> str:
    """Why a column or gene set an earlier apply wrote is not replaced without ``overwrite``."""
    return f'An earlier apply (step {step}) wrote "{name}" with other values; overwrite replaces it.'


def taken_set_note(name: str) -> str:
    return f'The object already has a gene set "{name}" that no session wrote; it is never replaced.'


def taken_column_note(name: str) -> str:
    return f'The object already has an obs column "{name}" that Cellucid did not write; it is never replaced.'


def own_var_column_note(name: str, column: str) -> str:
    """Why a gene set no apply wrote gets no ``var`` column: one of that name is the dataset's own."""
    return f'The gene set "{name}" has no var column: the var column "{column}" is the dataset\'s own.'


def missing_genes_note(name: str, missing: list[str]) -> str:
    """The genes a gene set names that the dataset lacks, left out of its ``var`` column (the first three)."""
    shown = ", ".join(missing[:3]) + (", …" if len(missing) > 3 else "")
    return (
        f'The gene set "{name}" names {len(missing)} {plural(len(missing), "gene")} the dataset does not have '
        f"({shown}), left out of its var column."
    )


def steps_only_message() -> str:
    return (
        "This session was shared as its steps only, without the values they leave behind, so it cannot be "
        "applied to an object; replay it in Cellucid on its dataset instead."
    )


def _require_session(session: Any) -> CellucidSession:
    if isinstance(session, CellucidSession):
        return session
    if isinstance(session, str | Path):
        return CellucidSession(session)
    raise TypeError("session must be a CellucidSession or the path of a .cellucid-session file")


# ---------------------------------------------------------------------------
# The dataset the session was saved on


def _check_fingerprint(fingerprint: dict[str, Any], adata: Any, *, open_cells: int) -> Any:
    """Check the object against the session's dataset before anything is planned.

    Returns what the session is checked and planned against: ``adata``
    itself, or for a MuData the cells and the modality its data mapping
    names, as the viewer showed them (:class:`MuDataSessionView`). The order
    is the viewer's: the dataset's id (when the object names one in
    ``uns['cellucid']['id']`` and the session names one), the cell count,
    then the data mapping -- every key it names must exist on the object --
    and the gene count of the axis that mapping reads.
    """
    record = getattr(adata, "uns", {}).get("cellucid") if isinstance(getattr(adata, "uns", None), Mapping) else None
    own_id = record.get("id") if isinstance(record, Mapping) else None
    dataset_id = fingerprint["datasetId"]
    if isinstance(own_id, str) and own_id and dataset_id is not None and own_id != dataset_id:
        raise ValueError(
            f"This session was saved on the dataset {dataset_id!r}, and this object is {own_id!r} "
            "(uns['cellucid']['id']). Apply it to the dataset it was saved on."
        )
    if fingerprint["cellCount"] != open_cells:
        raise ValueError(
            f"This session was saved on {fingerprint['cellCount']} {plural(fingerprint['cellCount'], 'cell')}, "
            f"and this object has {open_cells}. Apply it to the dataset it was saved on."
        )
    mapping = fingerprint["dataMapping"]
    shown = adata
    if is_mudata(adata):
        require_mapping_modality(mapping, [str(name) for name in adata.mod])
        shown = MuDataSessionView(adata, mapping[MODALITY_KEY])
        mapping = {key: value for key, value in mapping.items() if key != MODALITY_KEY}
    require_mapping_matches_anndata(shown, mapping)
    expected_var_count = mapped_var_count(shown, mapping)
    if fingerprint["varCount"] != expected_var_count:
        axis = (
            "the raw matrix's genes" if mapping["expression"] == "raw"
            else "genes on a mapping without expression" if mapping["expression"] is None
            else "genes"
        )
        raise ValueError(
            f"This session was saved on a dataset of {fingerprint['varCount']} {plural(fingerprint['varCount'], 'gene')}, "
            f"and this object has {expected_var_count} {axis}. Apply it to the dataset it was saved on."
        )
    return shown


def cell_order_digest(values: np.ndarray, positions: np.ndarray | None = None) -> str:
    """The cell-order digest of some cells' stored embedding values
    (``data/cell-order-digest.js``): two MurmurHash3 lanes per cell over its
    original position and its float32 values on up to three axes, summed,
    finalised with the cell count; 16 lower-case hex digits.

    ``values`` is one row per cell, 1 to 3 columns; ``positions`` each row's
    original position (None: row ``i`` is position ``i``).
    """
    array = np.ascontiguousarray(values, dtype=np.float32)
    if array.ndim == 1:
        array = array[:, None]
    cells, width = array.shape
    if not 1 <= width <= 3:
        raise ValueError("A cell-order digest is taken over 1, 2 or 3 columns.")
    bits = array.view(np.uint32).astype(np.uint64)
    nan = (bits & 0x7F800000) == 0x7F800000
    bits = np.where(nan & ((bits & 0x007FFFFF) != 0), np.uint64(0x7FC00000), bits)
    mask = np.uint64(0xFFFFFFFF)

    def mul(a: np.ndarray, b: int) -> np.ndarray:
        return (a * np.uint64(b)) & mask

    def rotl(a: np.ndarray, r: int) -> np.ndarray:
        return ((a << np.uint64(r)) | (a >> np.uint64(32 - r))) & mask

    def scramble(word: np.ndarray) -> np.ndarray:
        return mul(rotl(mul(word, 0xCC9E2D51), 15), 0x1B873593)

    def step(h: np.ndarray, k: np.ndarray) -> np.ndarray:
        return (mul(rotl(h ^ k, 13), 5) + np.uint64(0xE6546B64)) & mask

    def fmix(h: np.ndarray) -> np.ndarray:
        h = h ^ (h >> np.uint64(16))
        h = mul(h, 0x85EBCA6B)
        h = h ^ (h >> np.uint64(13))
        h = mul(h, 0xC2B2AE35)
        mixed: np.ndarray = h ^ (h >> np.uint64(16))
        return mixed

    zero = scramble(np.zeros(cells, dtype=np.uint64))
    words = [scramble(bits[:, axis]) if axis < width else zero for axis in range(3)]
    position = np.arange(cells, dtype=np.uint64) if positions is None else np.asarray(positions, dtype=np.uint64)
    scrambled = scramble(position & mask)
    lanes = []
    for seed in (0x9E3779B1, 0x85EBCA6B):
        h = np.full(cells, seed, dtype=np.uint64)
        h = step(h, scrambled)
        for word in words:
            h = step(h, word)
        lanes.append(int(fmix(h ^ np.uint64(16)).sum() & 0xFFFFFFFF))
    low, high = lanes

    def fmix_one(value: int) -> int:
        return int(fmix(np.array([value & 0xFFFFFFFF], dtype=np.uint64))[0])

    first = fmix_one(high ^ fmix_one(cells ^ 0x297A2D39))
    second = fmix_one(low ^ fmix_one(cells ^ 0x2C1B3C6D))
    return f"{first:08x}{second:08x}"


def _check_cell_order(fingerprint: dict[str, Any], shown: Any, kept: np.ndarray | None) -> None:
    """The cell-order digest of the embedding the session names, over the
    cells it was saved on, when the object has that embedding: a mismatch
    means the rows are other cells (or in another order), so every cell
    position the session holds would name the wrong cell."""
    order = fingerprint["cellOrder"]
    key = order["embedding"]["key"]
    columns = list(order["embedding"]["columns"])
    obsm = getattr(shown, "obsm", None)
    if obsm is None or key not in obsm:
        return
    values = np.asarray(obsm[key])
    if values.ndim != 2 or max(columns) >= values.shape[1]:
        return
    rows = values[:, columns] if kept is None else values[np.asarray(kept, dtype=np.int64)][:, columns]
    digest = cell_order_digest(rows, kept)
    if digest != order["digest"]:
        raise ValueError(
            f"This object's cells are not the cells the session was saved on, in the same order: the cell-order "
            f"digest of obsm[{key!r}] (columns {', '.join(map(str, columns))}) is {digest}, and the session's is "
            f"{order['digest']}. Put the cells in their original order, or apply the session to the object it was "
            "saved on."
        )


def required_fields(state: Mapping[str, Any]) -> dict[str, list[str]]:
    """The dataset fields a step's state uses (``collectRequiredFields``): every
    obs column and gene a view colours by or filters on, or a highlight group
    is bound to, by the key the dataset gives it; fields the session made
    itself are not listed. Sorted as JavaScript sorts text."""
    required: dict[str, set[str]] = {"obs": set(), "var": set()}

    def add(ref: Any) -> None:
        if isinstance(ref, Mapping) and isinstance(ref.get("key"), str) and ref.get("source") in required:
            required[ref["source"]].add(ref["key"])

    for domain in ("filters", "colors"):
        value = state.get(domain)
        for view in (value or {}).get("views", []) if isinstance(value, Mapping) else []:
            for entry in view or []:
                if isinstance(entry, Mapping):
                    add(entry.get("field"))
    active = state.get("active-field")
    for ref in (active or {}).get("views", []) if isinstance(active, Mapping) else []:
        add(ref)
    highlights = state.get("highlights")
    for page in (highlights or {}).get("pages", []) if isinstance(highlights, Mapping) else []:
        for group in page.get("groups", []):
            add(group.get("field"))

    def utf16(text: str) -> bytes:
        return text.encode("utf-16-be", "surrogatepass")

    return {source: sorted(keys, key=utf16) for source, keys in required.items()}


def _require_session_fields(shown: Any, requires: Mapping[str, list[str]], mapping: dict[str, Any], step: str) -> None:
    """Refuse an object that lacks a field the step uses, naming all of them."""
    obs_keys = set(shown.obs.columns)
    missing_obs = [key for key in requires["obs"] if key not in obs_keys]
    missing_var: list[str] = []
    if requires["var"]:
        gene_names = set(resolve_gene_names(shown, mapping).names)
        missing_var = [key for key in requires["var"] if key not in gene_names]
    if not missing_obs and not missing_var:
        return
    parts: list[str] = []
    if missing_obs:
        parts.append(f"the obs {plural(len(missing_obs), 'column')} " + ", ".join(repr(key) for key in missing_obs))
    if missing_var:
        parts.append(f"the {plural(len(missing_var), 'gene')} " + ", ".join(repr(key) for key in missing_var))
    raise ValueError(f"The state at step {step} uses " + " and ".join(parts) + ", which this dataset does not have.")


def _read_gene_sets(value: Any) -> list[dict[str, Any]]:
    """The gene sets of a ``gene-sets`` value."""
    record = _require_exact_keys(value, GENE_SETS_KEYS, label="gene-sets value")
    sets: list[dict[str, Any]] = []
    for index, raw in enumerate(_require_list(record["sets"], label="gene-sets value sets")):
        label = f"gene-sets value sets[{index}]"
        entry = _require_exact_keys(raw, GENE_SET_KEYS, label=label)
        name = _require_text(entry["name"], label=f"{label}.name")
        if not isinstance(entry["description"], str):
            raise ValueError(f"This session's {label}.description must be a text.")
        genes = [
            _require_text(gene, label=f"{label}.genes[{position}]")
            for position, gene in enumerate(_require_list(entry["genes"], label=f"{label}.genes"))
        ]
        sets.append({"name": name, "description": entry["description"], "genes": genes})
    return sets


# ---------------------------------------------------------------------------
# What earlier applies wrote


def history_applied(adata: Any) -> list[dict[str, Any]]:
    """Every apply recorded on an object (``uns['cellucid']['history_applied']``),
    oldest first, each the model's ``Provenance`` as JSON."""
    uns = getattr(adata, "uns", None)
    if not isinstance(uns, Mapping):
        raise TypeError("adata.uns must be a dictionary")
    record = uns.get("cellucid")
    if record is None:
        return []
    if not isinstance(record, Mapping):
        raise TypeError("adata.uns['cellucid'] must be a dictionary if present")
    text = record.get(HISTORY_APPLIED_KEY)
    if text is None:
        return []
    if not isinstance(text, str):
        raise ValueError("adata.uns['cellucid']['history_applied'] must be JSON text (a list of apply records).")
    try:
        entries = json.loads(text)
    except ValueError as error:
        raise ValueError("adata.uns['cellucid']['history_applied'] is not JSON.") from error
    if not isinstance(entries, list):
        raise ValueError("adata.uns['cellucid']['history_applied'] must be a list of apply records.")
    for index, entry in enumerate(entries):
        try:
            _rules.check_record("Provenance", entry, f"history_applied[{index}]")
        except _rules.HistoryModelError as error:
            raise ValueError(f"adata.uns['cellucid']['history_applied'] is not as an apply writes it: {error}.") from error
    return entries


def _same_column(existing: pd.Series, planned: pd.Series) -> bool:
    return bool(existing.reset_index(drop=True).equals(planned.reset_index(drop=True)))


# ---------------------------------------------------------------------------
# The plan


@dataclass(frozen=True)
class _VarColumn:
    """One gene set's ``bool`` column of the ``var`` its data mapping reads."""

    name: str
    values: np.ndarray


@dataclass(frozen=True)
class _Writes:
    """What applying a plan writes: obs columns (added, replaced, rebuilt with
    renamed categories, renamed), gene sets, and their ``var`` columns (on
    ``adata.var``, or on modality ``var_modality``'s of a MuData), with what
    was left out of those (``var_notes``, warned when written)."""

    columns: list[_ColumnPlan]
    replaced: list[_ColumnPlan]
    rename: dict[str, str]
    gene_sets: dict[str, dict[str, Any]]
    renamed: list[dict[str, Any]]
    var_columns: list[_VarColumn]
    var_modality: str | None
    var_notes: list[str]


@dataclass(frozen=True)
class ApplyPlan(model.ApplyPlan):
    """What applying one step to an object would change, written nowhere yet.

    Every field of the history model's ``ApplyPlan`` (``session_sha256``,
    ``tracker``, ``branch``, ``step``, ``requires``, ``rows``); ``table`` is
    the rows as a DataFrame; :meth:`apply` writes them into the object the
    plan was made for.
    """

    @property
    def table(self) -> pd.DataFrame:
        """One row per change (``kind``, ``domain``, ``name``, ``column``,
        ``from``, ``cells``, ``status``, ``note``), in domain order."""
        return pd.DataFrame(
            [[row.kind, row.domain, row.name, row.column, list(row.from_), row.cells, row.status, row.note] for row in self.rows],
            columns=["kind", "domain", "name", "column", "from", "cells", "status", "note"],
        )

    @property
    def writes_anything(self) -> bool:
        """Whether applying it would write anything (a row to write or replace)."""
        return any(row.status in ("write", "replace") for row in self.rows)

    def apply(self) -> ApplyReport:
        """Write the plan into the object it was made for, all at once, and
        record it in ``uns['cellucid']['history_applied']``.

        The plan is made again first: an object changed since (another apply,
        a column added) is refused rather than half-written.
        """
        context = cast(_PlanContext, getattr(self, _CONTEXT))
        again = _plan(context.session, context.adata, tracker=context.tracker, branch=context.branch, step=context.step,
                      overwrite=context.overwrite, keep_session=False)
        if again.rows != self.rows:
            raise ValueError("The object changed since this plan was made; make the plan again (dry_run=True) and apply that.")
        object.__setattr__(again, _CONTEXT, replace(cast(_PlanContext, getattr(again, _CONTEXT)), kept=context.kept))
        return _write(again)

    def __repr__(self) -> str:
        counts = {status: sum(1 for row in self.rows if row.status == status) for status in model.APPLY_ROW_STATUSES}
        summary = ", ".join(f"{count} {status}" for status, count in counts.items() if count)
        return f"ApplyPlan(step {self.step} of branch {self.branch} of tracker {self.tracker}: {summary or 'no change'})"


@dataclass(frozen=True)
class ApplyReport(model.ApplyReport):
    """What an apply did: its plan, whether it wrote anything (``applied``),
    and the record it wrote (``provenance``, None when nothing was written).
    ``adata`` is the object it was applied to."""

    @property
    def adata(self) -> Any:
        return getattr(self, _ADATA)


@dataclass(frozen=True)
class _PlanContext:
    session: CellucidSession
    adata: Any
    tracker: str | None
    branch: str | None
    step: str | None
    overwrite: bool
    # The session file's bytes, read once when the plan is made with
    # keep_session (so the file may go before the plan is applied).
    kept: bytes | None
    writes: _Writes


def _choose(session: CellucidSession, *, tracker: str | None, branch: str | None, step: str | None) -> tuple[Tracker, str, Step]:
    """The tracker, branch and step an apply takes: the step asked for (on the
    branch asked for, when both are), a branch's tip, or the tracker's head
    on its active branch."""
    found = session.tracker(tracker)
    if step is not None:
        chosen = found.step(step)
        if branch is None:
            return found, chosen.branch, chosen
        on = found.branch(branch)
        if chosen.id not in {item.id for item in on.path}:
            raise ValueError(
                f'Step {chosen.id} is not on the branch "{on.name}" (its path runs from the root to step {on.tip}); '
                f'it was made on "{found.branch(chosen.branch).name}".'
            )
        return found, on.id, chosen
    if branch is not None:
        on = found.branch(branch)
        return found, on.id, found.step(on.tip)
    return found, found.active_branch, found.step(found.head_id)


def _row(kind: str, domain: str, name: str, status: str, *, column: str | None = None, from_: list[str] | None = None,
         cells: int | None = None, note: str | None = None) -> dict[str, Any]:
    return {"kind": kind, "domain": domain, "name": name, "column": column, "from": list(from_ or []), "cells": cells,
            "status": status, "note": note}


def _session_bytes(session: CellucidSession) -> bytes:
    """The session file's bytes, which must still be the file that was read."""
    data = session.path.read_bytes()
    if hashlib.sha256(data).hexdigest() != session.sha256:
        raise ValueError(f"{session.path.name} changed since it was read; read it again and make the plan again.")
    return data


def _plan(
    session: CellucidSession,
    adata: Any,
    *,
    tracker: str | None,
    branch: str | None,
    step: str | None,
    overwrite: bool,
    keep_session: bool,
) -> ApplyPlan:
    """Check the object and plan every row, writing nothing."""
    if session.scope.content == "steps":
        raise ValueError(steps_only_message())
    chosen_tracker, branch_id, chosen = _choose(session, tracker=tracker, branch=branch, step=step)
    if not hasattr(adata, "obs") or not hasattr(adata, "var") or not hasattr(adata, "uns"):
        raise TypeError("adata must be an AnnData-compatible object")
    if type(getattr(adata, "isbacked", None)) is not bool:
        raise TypeError("adata.isbacked must be exactly True or False")

    # An AnnData read from a file in priority order holds its cells in the
    # stored order; a session names cells in their original order.
    from ..prepare_data._priority import stored_order_of

    stored = stored_order_of(adata)
    if stored.notice is not None:
        warnings.warn(stored.notice, UserWarning, stacklevel=4)
    if stored.positions is not None:
        raise ValueError(
            "adata holds its cells in a Cellucid-ready file's stored order (uns['cellucid']['priority']), "
            "and a session names cells in their original order. Put them back first: "
            "adata = cellucid.restore_original_order(adata)."
        )
    n_obs = int(adata.n_obs)
    saved_on = session.cells
    open_cells = n_obs
    kept = None
    if saved_on is not None:
        if n_obs != saved_on["spec"]["nObs"]:
            raise ValueError(
                f"This session was saved on {saved_on['spec']['k']} of {saved_on['spec']['nObs']} "
                f"{plural(saved_on['spec']['nObs'], 'cell')}, and this object has {n_obs}. Apply it to the dataset it was saved on."
            )
        kept = saved_on["kept"]
        open_cells = int(kept.size)
    fingerprint = session.dataset_fingerprint
    shown = _check_fingerprint(fingerprint, adata, open_cells=open_cells)
    _check_cell_order(fingerprint, shown, kept)
    mapping = {key: value for key, value in fingerprint["dataMapping"].items() if key != MODALITY_KEY}
    raw = chosen_tracker._raw_state(chosen.id)
    root = chosen_tracker._hashes(chosen_tracker.root_id)
    hashes = chosen_tracker._hashes(chosen.id)
    for domain in DATA_DOMAINS:
        if domain not in raw:
            raise ValueError(f"This session's state at step {chosen.id} holds no {domain} value.")
    plain = chosen.state()
    requires = required_fields(plain)
    _require_session_fields(shown, requires, mapping, chosen.id)

    applied = history_applied(adata)
    owners: dict[str, str] = {}
    for record in applied:
        for column in record["columns"]:
            owners[column] = record["step"]
    earlier_renames = [entry for record in applied for entry in record["renamed"]]

    known = known_cells(saved_on=saved_on, cells_value=raw[CELLS_DOMAIN], shown=shown, n_obs=n_obs, step=chosen.id)
    overlays = _require_exact_keys(raw[FIELD_OVERLAYS_DOMAIN], OVERLAY_KEYS, label="field-overlays value")
    modality_columns = set(map(str, shown.obs.columns)) - set(map(str, adata.obs.columns))
    labels = _plan_labels(overlays, adata.obs, modality_columns=modality_columns, earlier=earlier_renames)
    present_after = (set(map(str, adata.obs.columns)) - set(labels.rename)) | set(labels.rename.values())

    columns: list[_ColumnPlan] = []

    def settle(candidate: _Candidate) -> dict[str, Any]:
        """A column's row: written, already there, replaced, or why not."""
        name = candidate.name
        if candidate.values is None:
            return _row(candidate.kind, candidate.domain, name, "not-applied", cells=candidate.cells, note=candidate.note)
        if name not in present_after:
            columns.append(_ColumnPlan(name=name, values=candidate.values, metadata=candidate.metadata))
            return _row(candidate.kind, candidate.domain, name, "write", cells=candidate.cells)
        existing_name = next((old for old, new in labels.rename.items() if new == name), name)
        existing = adata.obs[existing_name]
        if name not in owners:
            return _row(candidate.kind, candidate.domain, name, "not-applied", cells=candidate.cells, note=taken_column_note(name))
        if _same_column(existing, candidate.values):
            return _row(candidate.kind, candidate.domain, name, "unchanged", cells=candidate.cells)
        if overwrite:
            columns.append(_ColumnPlan(name=name, values=candidate.values, metadata=candidate.metadata))
            return _row(candidate.kind, candidate.domain, name, "replace", cells=candidate.cells,
                        note=f"Written by an earlier apply (step {owners[name]}); replaced.")
        return _row(candidate.kind, candidate.domain, name, "not-applied", cells=candidate.cells, note=owned_note(name, owners[name]))

    field_rows = [settle(candidate) for candidate in _plan_user_fields(overlays["userFields"], shown, mapping=mapping, known=known)]
    step_sets = _read_gene_sets(raw[GENE_SETS_DOMAIN])
    gene_set_rows, gene_sets = _plan_gene_sets(step_sets, adata, session, overwrite=overwrite)
    var_columns, var_notes = _plan_gene_set_columns(step_sets, gene_set_rows, gene_sets, shown=shown, mapping=mapping)
    highlight_rows = [settle(candidate) for candidate in _plan_highlights(raw[HIGHLIGHTS_DOMAIN], n_obs=n_obs, index=adata.obs_names, known=known)]
    _refuse_original_position(columns, labels.replaced, labels.rename)
    presentation = {
        domain: _row("presentation", domain, DOMAIN_NAMES[domain], "not-applied", note=PRESENTATION_NOTE)
        for domain in model.DOMAINS
        if domain not in DATA_DOMAINS and hashes.get(domain) != root.get(domain)
    }
    by_domain: dict[str, list[dict[str, Any]]] = {
        FIELD_OVERLAYS_DOMAIN: [*labels.rows, *field_rows],
        GENE_SETS_DOMAIN: gene_set_rows,
        HIGHLIGHTS_DOMAIN: highlight_rows,
    }
    rows = [row for domain in model.DOMAINS for row in by_domain.get(domain, [presentation[domain]] if domain in presentation else [])]
    plan = ApplyPlan(
        session_sha256=session.sha256,
        tracker=chosen_tracker.id,
        branch=branch_id,
        step=chosen.id,
        requires=model.Requires(obs=tuple(requires["obs"]), var=tuple(requires["var"])),
        rows=tuple(model.build("ApplyPlanRow", row) for row in rows),
    )
    writes = _Writes(
        columns=columns, replaced=labels.replaced, rename=labels.rename, gene_sets=gene_sets, renamed=labels.renamed,
        var_columns=var_columns, var_modality=fingerprint["dataMapping"][MODALITY_KEY] if is_mudata(adata) else None,
        var_notes=var_notes,
    )
    object.__setattr__(plan, _CONTEXT, _PlanContext(
        session=session, adata=adata, tracker=tracker, branch=branch, step=step, overwrite=overwrite,
        kept=_session_bytes(session) if keep_session else None, writes=writes,
    ))
    return plan


def _plan_gene_sets(
    sets: list[dict[str, Any]], adata: Any, session: CellucidSession, *, overwrite: bool,
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    """A row per gene set of the step, and the entries to write."""
    record = adata.uns.get("cellucid")
    current = record.get(GENE_SETS_KEY) if isinstance(record, Mapping) else None
    if current is not None and not isinstance(current, Mapping):
        raise TypeError("adata.uns['cellucid']['gene_sets'] must be a dictionary if present")
    current = current or {}
    rows: list[dict[str, Any]] = []
    writes: dict[str, dict[str, Any]] = {}
    for entry in sets:
        name = entry["name"]
        written = {"genes": list(entry["genes"]), "description": entry["description"], "origin": "session",
                   "pulled_from": session.path.name, "id": name}
        held = current.get(name)
        if held is None:
            writes[name] = written
            rows.append(_row("gene-set", GENE_SETS_DOMAIN, name, "write"))
            continue
        same = isinstance(held, Mapping) and list(held.get("genes", [])) == written["genes"] and held.get("description") == written["description"]
        owned = isinstance(held, Mapping) and held.get("origin") == "session"
        if same:
            rows.append(_row("gene-set", GENE_SETS_DOMAIN, name, "unchanged"))
        elif not owned:
            rows.append(_row("gene-set", GENE_SETS_DOMAIN, name, "not-applied", note=taken_set_note(name)))
        elif overwrite:
            writes[name] = written
            rows.append(_row("gene-set", GENE_SETS_DOMAIN, name, "replace", note="Written by an earlier apply; replaced."))
        else:
            rows.append(_row("gene-set", GENE_SETS_DOMAIN, name, "not-applied",
                             note=f'An earlier apply wrote the gene set "{name}" with other genes; overwrite replaces it.'))
    return rows, writes


def gene_set_columns(names: list[str]) -> dict[str, str]:
    """Each gene set's ``var`` column, by name: ``cellucid_gene_set__<name>``
    (the highlight columns' rule, ``column_part``: each run
    of characters other than ASCII letters and digits one underscore, none at
    either end; ``gene_set`` when nothing is left), ``_2``, ``_3``... on a
    clash. Every set of the step is named, in its order, so a set's column is
    the same whichever of them an apply writes."""
    taken: set[str] = set()
    return {name: unique_column(f"{GENE_SET_PREFIX}{column_part(name, 'gene_set')}", taken) for name in names}


def _plan_gene_set_columns(
    sets: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    writes: dict[str, dict[str, Any]],
    *,
    shown: Any,
    mapping: Mapping[str, Any],
) -> tuple[list[_VarColumn], list[str]]:
    """The ``var`` column of each gene set an apply writes or replaces, and
    what is left out of them.

    The columns go to the ``var`` of the genes the data mapping reads (``X``
    or a layer: ``var``, a MuData modality's for a MuData); a mapping that
    reads ``raw`` or no expression has none, and the sets stay in ``uns``
    only. A column of a set written (no apply wrote it before) whose name
    ``var`` holds already is the dataset's own and is not written; a set
    replaced replaces its column. ``True`` on each row whose gene, as the
    mapping names genes, the set lists; a gene the dataset lacks is left out.
    """
    if not writes:
        return [], []
    if mapping["expression"] in (None, "raw"):
        return [], [GENE_SETS_UNS_ONLY_NOTE]
    names = resolve_gene_names(shown, mapping).names
    row_of = {gene: position for position, gene in enumerate(names)}
    held = set(shown.var.columns)
    status = {row["name"]: row["status"] for row in rows}
    column_of = gene_set_columns([entry["name"] for entry in sets])
    columns: list[_VarColumn] = []
    notes: list[str] = []
    for name, entry in writes.items():
        column = column_of[name]
        if column in held and status[name] != "replace":
            notes.append(own_var_column_note(name, column))
            continue
        values = np.zeros(len(names), dtype=bool)
        missing: list[str] = []
        for gene in entry["genes"]:
            position = row_of.get(gene)
            if position is None:
                missing.append(gene)
            else:
                values[position] = True
        if missing:
            notes.append(missing_genes_note(name, missing))
        columns.append(_VarColumn(name=column, values=values))
    return columns, notes


def _refuse_original_position(columns: list[_ColumnPlan], replaced: list[_ColumnPlan], rename: Mapping[str, str]) -> None:
    """Refuse a plan that would write, replace or rename
    ``obs['cellucid_original_position']`` (or rename a column onto it): the
    column a file stored in priority order keeps each row's original position
    in (:func:`~cellucid.restore_original_order` reads and drops it)."""
    from ..prepare_data._priority import ORIGINAL_POSITION_COLUMN

    touched = {plan.name for plan in (*columns, *replaced)} | set(rename) | set(rename.values())
    if ORIGINAL_POSITION_COLUMN in touched:
        raise ValueError(ORIGINAL_POSITION_MESSAGE)


def _iso_now() -> str:
    now = datetime.now(UTC)
    return now.strftime("%Y-%m-%dT%H:%M:%S.") + f"{now.microsecond // 1000:03d}Z"


def _write(plan: ApplyPlan) -> ApplyReport:
    """Write a plan made just now: every obs column at once, the gene sets'
    ``var`` columns, then ``uns``; then warn what the ``var`` columns left out."""
    context = cast(_PlanContext, getattr(plan, _CONTEXT))
    target = context.adata
    if not plan.writes_anything:
        report = ApplyReport(plan=plan, applied=False, provenance=None)
        object.__setattr__(report, _ADATA, target)
        return report
    writes = context.writes
    provenance = {
        "session_sha256": plan.session_sha256,
        "tracker": plan.tracker,
        "branch": plan.branch,
        "step": plan.step,
        "applied_at": _iso_now(),
        "cellucid_version": package_version(),
        "package": "python",
        "columns": [row.name for row in plan.rows if row.kind in ("user-field", "highlight") and row.status in ("write", "replace")],
        "renamed": copy.deepcopy(writes.renamed),
        "skipped": [{"kind": row.kind, "name": row.name, "note": row.note or ""} for row in plan.rows if row.status == "not-applied"],
        "session": None if context.kept is None else base64.b64encode(context.kept).decode("ascii"),
    }
    _rules.check_record("Provenance", provenance, "provenance")
    uns = target.uns
    existing = uns.get("cellucid")
    record = copy.deepcopy(dict(existing)) if isinstance(existing, Mapping) else {}
    if writes.gene_sets:
        sets = dict(record.get(GENE_SETS_KEY) or {})
        sets.update(writes.gene_sets)
        record[GENE_SETS_KEY] = sets
    entries = history_applied(target)
    record[HISTORY_APPLIED_KEY] = json.dumps([*entries, provenance], ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    # A new column goes after var's columns; a replaced one keeps its place.
    var_owner = target if writes.var_modality is None else target.mod[writes.var_modality]
    next_var = None
    if writes.var_columns:
        next_var = var_owner.var.copy(deep=False)
        for column in writes.var_columns:
            next_var[column.name] = column.values
    write_columns(target, writes.columns, replace=writes.replaced, rename=writes.rename)
    if next_var is not None:
        var_owner.var = next_var
    target.uns["cellucid"] = record
    for note in writes.var_notes:
        warnings.warn(note, UserWarning, stacklevel=3)
    report = ApplyReport(plan=plan, applied=True, provenance=model.build("Provenance", provenance))
    object.__setattr__(report, _ADATA, target)
    return report


def _require_options(**flags: object) -> None:
    for label, flag in flags.items():
        if type(flag) is not bool:
            raise TypeError(f"{label} must be exactly True or False")


def apply_cellucid_session_to_anndata(
    adata: Any,
    session: CellucidSession | str | Path,
    *,
    tracker: str | None = None,
    branch: str | None = None,
    step: str | None = None,
    dry_run: bool = True,
    overwrite: bool = False,
    keep_session: bool = False,
) -> ApplyPlan | ApplyReport:
    """What applying one step of a session to its AnnData would change, or the change made.

    A session file holds the history of every step a person, an agent or a
    notebook took (:class:`~cellucid.CellucidSession`); this applies the
    state one step left -- its labels, user-defined fields, gene sets and
    highlight groups -- to the object it was saved on. With ``dry_run=True``
    (the default) nothing is written: the returned :class:`ApplyPlan` lists
    every change, one row each, with what would happen to it; ``plan.apply()``
    then writes it. With ``dry_run=False`` the plan is applied at once and
    the :class:`ApplyReport` returned. Every check runs before anything is
    written, and everything is written at once, so a failed apply leaves the
    object as it was.

    What is checked (only what the chosen step needs):

    - the dataset: the cell count (of the subsample the session was saved
      on, when it was), the data mapping's keys and gene count, the dataset
      id when the object names one (``uns['cellucid']['id']``), and the cell
      order: when the object has the embedding the session names, the
      cell-order digest of its stored values over the session's cells must be
      the session's (the rows are those cells, in that order); a digest is a
      sum over the cells, so a mismatch says the cells differ, not which;
    - the fields the chosen step's state uses (the obs columns and genes its
      views colour by or filter on, and its highlight groups are bound to);
      the session's own ``requires`` is the shown state's.

    What a plan's rows are, in domain order: ``column-renamed`` (an obs
    column renamed in place), ``categories-renamed`` and
    ``categories-merged`` (a category renamed, or renamed into one the column
    has: the two become one), ``field-deleted`` (recorded, never applied:
    the column stays), ``user-field`` (a column per field the session made: a
    categorical, or ``float32`` values; a field deleted in the session is not
    applied), ``gene-set`` (``uns['cellucid']['gene_sets']``, and a ``bool``
    ``var`` column, below), ``highlight``
    (one ``bool`` column per group, ``cellucid_highlight__<page>__<group>``
    from the page's name and the group's label, ASCII letters and digits with
    every other run one underscore, ``_2``, ``_3``... when two read alike) and
    ``presentation`` (a domain the step changed that is how the data was shown
    -- filters, colours, views, the camera, analyses -- never applied).

    What each row's status is: ``write``; ``unchanged`` (the object holds it
    already: applying the same step again changes nothing and reports so);
    ``replace`` (a column or gene set an earlier apply wrote, holding other
    values, with ``overwrite=True``); ``not-applied`` with the reason in
    ``note`` -- a column or gene set the person made is never replaced.

    Each gene set written or replaced also gets a ``bool`` column in the
    ``var`` of the genes the data mapping reads (``adata.var``; for a MuData,
    the ``var`` of the modality the mapping names), named
    ``cellucid_gene_set__<name>`` by the highlight columns' rule over every
    set of the step (``_2``, ``_3``... on a clash), ``True`` on each row whose
    gene, as the mapping names genes, the set lists. A new column goes after
    ``var``'s columns; a replaced set's column is replaced in place. Warned
    (``UserWarning``) when written: genes the dataset lacks (left out of the
    column), a ``var`` column of that name beside a set no apply wrote (the
    dataset's own: kept, no column written), and a mapping that reads ``raw``
    or no expression (no ``var`` to mark: the sets stay in ``uns`` only).
    Nothing is removed: a gene set the step no longer holds keeps its record
    and its column.

    ``obs['cellucid_original_position']`` (each row's original position in a
    file stored in priority order) is never written, renamed or replaced: a
    plan that would is refused.

    A session saved on a subsample, or a step taken while the viewer had a
    subsample open, knows the state of those cells only: a cell it did not
    open is ``<NA>`` in a highlight column (the nullable ``boolean`` dtype)
    and missing in a field's -- unknown, not "not selected".

    Parameters
    ----------
    adata : anndata.AnnData or mudata.MuData
        The dataset the session was saved on, its cells in their original
        order (an AnnData read from a Cellucid-ready file in its stored
        priority order is refused: :func:`~cellucid.restore_original_order`
        puts it back). For a MuData the session is checked against what the
        viewer showed -- the cells with every modality placed among them,
        and the modality its data mapping names -- and labels and columns
        apply to the global ``mdata.obs``, gene sets' ``var`` columns to that
        modality's ``var``; ``mdata.uns`` records the apply. Written in place.
    session : CellucidSession, str, or Path
        The session, or the path of a ``.cellucid-session`` file (read and
        checked whole first). One shared as its steps only cannot be applied.
    tracker : str, optional
        A tracker by its name or its id. Default: the one the session shows.
    branch : str, optional
        A branch by its name or its id: its tip is applied, or ``step`` when
        it is on the branch's path.
    step : str, optional
        A step by its id. Default: the branch's tip, or with no branch the
        tracker's head (what the viewer showed).
    dry_run : bool
        Plan only (the default); False applies the plan at once.
    overwrite : bool
        Replace a column or gene set an earlier apply wrote with other values.
    keep_session : bool
        Keep the session file itself (base64) in the apply's record, so the
        history can be rebuilt from the object alone. A record holds no
        ``source`` (where the viewer opened the dataset, which may be a
        private link); the whole file kept this way does, unless it was
        shared without one.

    Returns
    -------
    ApplyPlan or ApplyReport
        The plan (``dry_run=True``), or the report of the apply: its plan,
        whether anything was written, and the record written to
        ``uns['cellucid']['history_applied']`` (JSON text; read it with
        :func:`cellucid.history_applied`).

    Raises
    ------
    TypeError
        An argument of the wrong type.
    ValueError
        A session file that is damaged, of another format, written by a newer
        Cellucid, or shared as its steps only; a tracker, branch or step it
        does not hold, or a step not on the branch named; an object that is
        not the dataset the session was saved on, or lacks a field the step
        uses; a plan that would write, rename or replace
        ``obs['cellucid_original_position']``; a value this version does not
        read.
    """
    _require_options(dry_run=dry_run, overwrite=overwrite, keep_session=keep_session)
    for label, value in (("tracker", tracker), ("branch", branch), ("step", step)):
        if value is not None and not isinstance(value, str):
            raise TypeError(f"{label} must be a string (a name or an id) or None")
    plan = _plan(_require_session(session), adata, tracker=tracker, branch=branch, step=step, overwrite=overwrite,
                 keep_session=keep_session)
    if dry_run:
        return plan
    return _write(plan)


def session_columns(session: Any, adata: Any) -> tuple[list[_Candidate], list[dict[str, Any]]]:
    """The highlight groups and user-defined fields a session shows (its
    active tracker's head) as columns of ``adata``, and its gene sets: what
    :func:`~cellucid.pull` writes from a session file, checked as an apply
    checks it."""
    found = _require_session(session)
    plan = _plan(found, adata, tracker=None, branch=None, step=None, overwrite=False, keep_session=False)
    chosen = found.tracker().step(plan.step)
    raw = found.tracker()._raw_state(chosen.id)
    saved_on = found.cells
    shown = adata if not is_mudata(adata) else MuDataSessionView(adata, found.dataset_fingerprint["dataMapping"][MODALITY_KEY])
    known = known_cells(saved_on=saved_on, cells_value=raw[CELLS_DOMAIN], shown=shown, n_obs=int(adata.n_obs), step=chosen.id)
    mapping = {key: value for key, value in found.dataset_fingerprint["dataMapping"].items() if key != MODALITY_KEY}
    overlays = _require_exact_keys(raw[FIELD_OVERLAYS_DOMAIN], OVERLAY_KEYS, label="field-overlays value")
    candidates = [
        *_plan_highlights(raw[HIGHLIGHTS_DOMAIN], n_obs=int(adata.n_obs), index=adata.obs_names, known=known),
        *_plan_user_fields(overlays["userFields"], shown, mapping=mapping, known=known),
    ]
    return [candidate for candidate in candidates if candidate.values is not None], _read_gene_sets(raw[GENE_SETS_DOMAIN])
