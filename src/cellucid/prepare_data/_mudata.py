"""
Writing a Cellucid-ready MuData: one standard MuData, each modality a Cellucid-ready AnnData.

A Cellucid-ready MuData (``.h5mu``, or a MuData Zarr store as a folder or a
``.zarr.zip``) is what ``mudata.write_h5mu`` / ``mudata.write_zarr`` would
write, laid out for the web app's range reads:

- the root is the MuData group (``encoding-type`` ``MuData``, ``axis`` 0:
  the modalities share cells): the cells' ``obs``, ``obsm`` (shared
  embeddings, vector fields and mudata's per-modality boolean masks),
  ``obsp`` (the shared neighbour graph), ``uns`` (``title`` and
  ``uns['cellucid']``), the global ``var`` (every modality's feature names),
  ``obsmap`` and ``varmap`` (where each modality's cells and features sit in
  the global ones, 1-based, 0 for absent), and ``mod`` (``mod-order``);
- each ``mod/<name>`` is a Cellucid-ready AnnData without ``uns['cellucid']``:
  gene-major (CSC) ``X``, its ``var`` and ``obs``, its own ``obsm`` and
  ``obsp['connectivities']``, every array chunked and compressed as
  :mod:`._store` lays out an AnnData.

The modalities are written one after the other, each converted only while it
is written, so memory holds one modality's converted matrix at a time beside
the MuData itself.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
from scipy import sparse

from .._console import console_print
from .._plural import plural
from ..data_mapping import default_graph
from ._inputs import (
    cell_matrix,
    require_distinct_element_keys,
    require_element_key,
    resolve_connectivities,
    resolve_expression,
    resolve_genes,
    resolve_obs,
    storable_columns,
    stored_expression,
)
from ._priority import ORIGINAL_POSITION_COLUMN, _require_positions, without_positions
from ._store import (
    WHOLE_CHUNK_ELEMENTS,
    Elements,
    _cell_chunks,
    _row_chunks,
    _vector_chunks,
    _Writer,
)


def _modality_error(name: str, error: Exception) -> Exception:
    """``error`` naming the modality, as a plain ``TypeError`` or ``ValueError``.

    Not ``type(error)(...)``: a subclass with its own constructor
    (``DataMappingError`` takes a code first) cannot be rebuilt from a message.
    """
    kind = TypeError if isinstance(error, TypeError) else ValueError
    return kind(f"Modality {name!r}: {error}")


@dataclass(frozen=True)
class Modality:
    """One modality: what is checked before writing, and how to convert its matrix while writing."""

    name: str
    obs_index: np.ndarray
    obs: dict[str, Any]
    var_index: np.ndarray
    var: dict[str, Any]
    obsm: dict[str, np.ndarray]
    obsp: dict[str, sparse.csr_matrix]
    obsmap: np.ndarray
    varmap: np.ndarray
    #: The modality's X, checked when called (:func:`._inputs.resolve_expression`)
    #: and laid out as it is written (:meth:`stored_matrix`), None without X.
    matrix: Callable[[], Any] | None
    #: Each row's original row in the modality, for a MuData in a verified
    #: stored order (its own ``obs['cellucid_original_position']``); None:
    #: each row its own.
    given: np.ndarray | None = None
    #: The modality row each written row holds, and the written row of each
    #: modality row (``_priority.CellPriority.modality``); None: as given.
    order: np.ndarray | None = None
    new_row: np.ndarray | None = None
    #: In priority order, the modality's rows of cells below each global
    #: prefix level (``_priority.CellPriority.modality``); None as given.
    prefix_levels: np.ndarray | None = None

    def stored_matrix(self) -> Any:
        """The modality's X as it is written: gene-major in its written row
        order with its prefix counts (sparse), or its rows in that order (dense)."""
        if self.matrix is None:
            return None
        try:
            return stored_expression(
                self.matrix(),
                genes=None,
                order=self.order,
                new_row=self.new_row,
                levels=self.prefix_levels,
                label=f"Modality {self.name!r} X",
            )
        except (TypeError, ValueError) as error:
            raise _modality_error(self.name, error) from error


def _unique_names(names: Any, *, what: str, fix: str) -> list[str]:
    listed = [str(name) for name in names]
    seen: set[str] = set()
    repeated: list[str] = []
    for name in listed:
        if name in seen and name not in repeated:
            repeated.append(name)
        seen.add(name)
    if repeated:
        shown = ", ".join(repr(name) for name in repeated[:5])
        more = f" and {len(repeated) - 5:,} more" if len(repeated) > 5 else ""
        raise ValueError(
            f"{what} repeat ({shown}{more}). MuData places each modality's cells among the "
            f"global ones by name, so a repeated name is ambiguous; make them unique: {fix}"
        )
    return listed


def _placement(values: Any, *, total: int, present: int, what: str, fix: str) -> np.ndarray:
    """``obsmap``/``varmap``: one uint32 per global entry, 1-based, 0 where absent."""
    placed = np.asarray(values)
    if placed.shape != (total,) or placed.dtype.kind not in "iu":
        raise ValueError(f"{what} must hold one integer per global entry ({total:,}); {fix}")
    if total and (placed.min() < 0 or placed.max() > present):
        raise ValueError(f"{what} points outside the modality's {present:,} {plural(present, 'entry', 'entries')}; {fix}")
    found = np.sort(placed[placed > 0])
    if not np.array_equal(found, np.arange(1, present + 1)):
        raise ValueError(f"{what} must place each of the modality's {present:,} {plural(present, 'entry', 'entries')} exactly once; {fix}")
    return placed.astype(np.uint32, copy=False)


def check_mudata(mdata: Any, selected: list[str] | None = None) -> list[str]:
    """Refuse a MuData Cellucid cannot show, saying what to do; return the modality names written.

    ``selected`` (``prepare(modalities=)``) names the modalities written, in
    order, and only those are checked; None is every one, in ``mdata.mod``
    order.
    """
    axis = int(getattr(mdata, "axis", 0))
    if axis != 0:
        raise ValueError(
            f"This MuData shares features, not cells (axis={axis}); Cellucid shows cells, so a "
            "Cellucid-ready MuData has modalities over shared cells (axis=0). Prepare each "
            "modality as an AnnData instead: prepare(..., data=mdata.mod['<name>'])."
        )
    names = [str(name) for name in mdata.mod] if selected is None else list(selected)
    if not names:
        raise ValueError("This MuData has no modalities (mdata.mod is empty); add one before prepare.")
    for name in names:
        require_element_key(name, "Modality")
    require_distinct_element_keys(names, "Modality")
    update = "mdata.obs_names_make_unique(); mdata.update()"
    _unique_names(mdata.obs_names, what="The MuData's cell names (mdata.obs_names)", fix=update)
    for name in names:
        modality = mdata.mod[name]
        if modality.n_obs == 0 or modality.n_vars == 0:
            empty = "cells" if modality.n_obs == 0 else "features"
            raise ValueError(
                f"Modality {name!r} has no {empty} ({modality.n_obs:,} x {modality.n_vars:,}). "
                f"Remove it before prepare: del mdata.mod[{name!r}]; mdata.update()."
            )
        _unique_names(
            modality.obs_names,
            what=f"The cell names of modality {name!r}",
            fix=f"mdata.mod[{name!r}].obs_names_make_unique(); mdata.update()",
        )
        if modality.X is None:
            raise ValueError(
                f"Modality {name!r} has no X. Give it the matrix the viewer reads, for example "
                f"mdata.mod[{name!r}].X = mdata.mod[{name!r}].layers['counts']."
            )
    return names


def modality_cell_arrays(modality: Any, name: str) -> dict[str, np.ndarray]:
    """A modality's ``obsm``: every numeric per-cell array, as float32, under its own key."""
    from ._data import is_cell_array

    arrays: dict[str, np.ndarray] = {}
    keys = [key for key, value in modality.obsm.items() if is_cell_array(value)]
    for key in keys:
        require_element_key(key, f"Modality {name!r} obsm key")
    require_distinct_element_keys(keys, f"Modality {name!r} obsm key")
    for key in keys:
        arrays[key] = cell_matrix(modality.obsm[key], f"Modality {name!r} obsm[{key!r}]")
    return arrays


def resolve_modality(mdata: Any, name: str, *, stored: bool = False) -> Modality:
    """Everything of one modality but its matrix, checked; the matrix is laid out when written.

    ``stored``: the MuData is in a verified stored order, so the modality's own
    ``obs['cellucid_original_position']`` is each row's original row in it
    (``Modality.given``), composed with the written order.
    """
    modality = mdata.mod[name]
    update = "rebuild them with mdata.update()"
    try:
        frame, note = storable_columns(modality.obs, f"mdata.mod[{name!r}].obs")
        if note is not None:
            console_print(note)
        obs_index, obs = resolve_obs(frame, modality.n_obs)
        given = None
        if stored and ORIGINAL_POSITION_COLUMN in obs:
            given = _require_positions(modality.obs[ORIGINAL_POSITION_COLUMN].to_numpy(), modality.n_obs)
        obs = without_positions(obs, f"mdata.mod[{name!r}].obs", verified=given is not None)
        frame, note = storable_columns(modality.var, f"mdata.mod[{name!r}].var")
        if note is not None:
            console_print(note)
        _, var_index, var = resolve_genes(modality.n_vars, frame, None, None)
        obsp = {}
        graph = default_graph(modality.obsp.keys(), modality=name)
        if graph is not None:
            obsp["connectivities"] = resolve_connectivities(modality.obsp[graph], modality.n_obs)
    except (TypeError, ValueError) as error:
        raise _modality_error(name, error) from error
    obsmap = _placement(
        mdata.obsmap[name],
        total=mdata.n_obs,
        present=modality.n_obs,
        what=f"mdata.obsmap[{name!r}]",
        fix=update,
    )
    varmap = _placement(
        mdata.varmap[name],
        total=mdata.n_vars,
        present=modality.n_vars,
        what=f"mdata.varmap[{name!r}]",
        fix=update,
    )

    def matrix() -> Any:
        try:
            return resolve_expression(modality.X, modality.n_obs)
        except (TypeError, ValueError) as error:
            raise _modality_error(name, error) from error

    return Modality(
        name=name,
        obs_index=obs_index,
        obs=obs,
        var_index=var_index,
        var=var,
        obsm=modality_cell_arrays(modality, name),
        obsp=obsp,
        obsmap=obsmap,
        varmap=varmap,
        matrix=matrix,
        given=given,
    )


@dataclass(frozen=True)
class MuDataLayout:
    """One Cellucid-ready MuData, checked; written by :func:`._store.write_store`."""

    obs_index: np.ndarray
    obs: dict[str, Any]
    obsm: dict[str, np.ndarray]
    #: mudata's per-modality boolean ``obsm`` masks, which cells each
    #: measured, as mudata keeps them (one column, ``(n_cells, 1)``).
    masks: dict[str, np.ndarray]
    obsp: dict[str, sparse.csr_matrix]
    uns: dict[str, Any]
    var_index: np.ndarray
    modalities: tuple[Modality, ...]

    def write(self, root: Any, kind: Literal["h5ad", "zarr"], zarr_format: int | None) -> None:
        writer = _Writer(root, kind, zarr_format)
        root.attrs.update({"encoding-type": "MuData", "encoding-version": "0.1.0", "axis": 0})
        writer.uns(root, "uns", self.uns)
        writer.dataframe(root, "obs", self.obs_index, self.obs, opening_categories=True, per_cell=True)
        writer.dataframe(root, "var", self.var_index, {})
        obsm = writer.dict_group(root, "obsm")
        for key, values in self.obsm.items():
            writer.array(obsm, key, values, _row_chunks(values.shape))
        for key, mask in self.masks.items():
            chunks = _row_chunks(mask.shape) if mask.ndim == 2 else _cell_chunks(len(mask))
            writer.array(obsm, key, mask, chunks)
        obsp = writer.dict_group(root, "obsp")
        for key, graph in self.obsp.items():
            writer.sparse_entries(writer.sparse_pointers(obsp, key, graph), graph, sliced=False)
        writer.dict_group(root, "varm")
        writer.dict_group(root, "varp")
        for group_key, placement in (("obsmap", "obsmap"), ("varmap", "varmap")):
            group = writer.dict_group(root, group_key)
            for modality in self.modalities:
                values = getattr(modality, placement)
                # A cell map has one entry per cell; a feature map one per feature.
                chunks = _cell_chunks(len(values)) if group_key == "obsmap" else _vector_chunks(len(values), WHOLE_CHUNK_ELEMENTS)
                writer.array(group, modality.name, values, chunks)
        mod = root.create_group("mod")
        mod.attrs["mod-order"] = [modality.name for modality in self.modalities]
        writer.flush()
        groups = {modality.name: mod.create_group(modality.name) for modality in self.modalities}
        for modality in self.modalities:
            x = modality.stored_matrix()
            _Writer(groups[modality.name], kind, zarr_format).write(
                Elements(
                    x=x,
                    obs_index=modality.obs_index,
                    obs=modality.obs,
                    var_index=modality.var_index,
                    var=modality.var,
                    obsm=modality.obsm,
                    obsp=modality.obsp,
                    uns={},
                )
            )
            del x


def global_var_index(mdata: Any) -> np.ndarray:
    return np.asarray([str(name) for name in mdata.var_names], dtype=object)


def placed_rows(values: np.ndarray, obsmap: np.ndarray) -> np.ndarray:
    """A modality's per-cell rows placed among the global cells, NaN where it measured none."""
    placed = np.full((len(obsmap), values.shape[1]), np.nan, dtype=values.dtype)
    present = obsmap > 0
    placed[present] = values[obsmap[present].astype(np.int64) - 1]
    return placed


def modality_obsm_widths(modalities: Mapping[str, Modality]) -> dict[str, int]:
    """The ``<modality>:<key>`` obsm arrays the viewer shows, by width."""
    return {
        f"{name}:{key}": int(values.shape[1])
        for name, modality in modalities.items()
        for key, values in modality.obsm.items()
    }
