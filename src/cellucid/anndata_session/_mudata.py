"""
A MuData as the viewer showed it, for applying a session to it.

The viewer shows a MuData's cells with the modalities placed among them
(``data/anndata/views/mudata.js``): the global ``obs`` plus every modality's
own ``obs`` columns as ``<modality>:<column>``, the global ``obsm`` plus every
modality's arrays as ``<modality>:<key>``, and the genes, matrices and
``var`` of the modality the session's data mapping names. Values of a cell a
modality did not measure are missing (NaN), never zero.

:class:`MuDataSessionView` is that picture with the attributes the session
checks and planners read from an AnnData, so a session is validated against
exactly what it was saved on. What the session adds -- highlight groups and
user-defined fields -- goes to the MuData's global ``obs``.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse


def _placement(mdata: Any, name: str) -> np.ndarray:
    """Each global cell's 0-based row in modality ``name``, -1 where it measured none."""
    return np.asarray(mdata.obsmap[name], dtype=np.int64) - 1


def _placed_column(values: pd.Series, rows: np.ndarray) -> Any:
    """One modality ``obs`` column placed among the global cells, missing where not measured."""
    present = rows >= 0
    if isinstance(values.dtype, pd.CategoricalDtype):
        codes = np.full(len(rows), -1, dtype=np.int64)
        codes[present] = values.cat.codes.to_numpy()[rows[present]]
        return pd.Categorical.from_codes(codes, dtype=values.dtype)
    array = values.to_numpy()
    if array.dtype.kind in "iubf":
        placed = np.full(len(rows), np.nan, dtype=np.float64 if array.dtype.kind != "f" else array.dtype)
        placed[present] = array[rows[present]]
        return placed
    out = np.full(len(rows), None, dtype=object)
    out[present] = array[rows[present]]
    return out


def _placed_rows(values: Any, rows: np.ndarray) -> np.ndarray:
    array = np.asarray(values.toarray() if sparse.issparse(values) else values, dtype=np.float64)
    shape = (len(rows),) if array.ndim == 1 else (len(rows), array.shape[1])
    placed = np.full(shape, np.nan, dtype=np.float64)
    present = rows >= 0
    placed[present] = array[rows[present]]
    return placed


class _PlacedRows:
    """A modality's per-cell array placed among the global cells, made only when read.

    A session is checked against the arrays' keys, widths and types alone; a
    placed float64 copy of every modality array was built up front and never
    read (two gigabytes for a few embeddings of two million cells).
    """

    def __init__(self, values: Any, rows: np.ndarray) -> None:
        self._values = values
        self._rows = rows
        self.shape = (len(rows),) if values.ndim == 1 else (len(rows), int(values.shape[1]))
        self.dtype = np.dtype(np.float64)

    def __array__(self, dtype: Any = None, copy: Any = None) -> np.ndarray:
        placed = _placed_rows(self._values, self._rows)
        return placed if dtype is None else placed.astype(dtype)


class _PlacedMatrix:
    """A modality's cells x genes matrix read with one row per global cell (NaN where not measured)."""

    def __init__(self, matrix: Any, rows: np.ndarray) -> None:
        self._matrix = matrix
        self._rows = rows
        self.shape = (len(rows), int(matrix.shape[1]))

    def __getitem__(self, key: Any) -> np.ndarray:
        if not (isinstance(key, tuple) and len(key) == 2 and key[0] == slice(None)):
            raise TypeError("A MuData modality's matrix is read one gene at a time: matrix[:, j].")
        column = self._matrix[:, key[1]]
        values = column.toarray().ravel() if sparse.issparse(column) else np.asarray(column).ravel()
        placed = np.full(len(self._rows), np.nan, dtype=np.float64)
        present = self._rows >= 0
        placed[present] = values[self._rows[present]]
        return placed


class _PlacedRaw:
    def __init__(self, raw: Any, rows: np.ndarray) -> None:
        self.var = raw.var
        self.n_vars = int(raw.n_vars)
        self.X = _PlacedMatrix(raw.X, rows)


class MuDataSessionView:
    """The MuData a session was saved on, as the viewer showed it with modality ``modality``."""

    def __init__(self, mdata: Any, modality: str) -> None:
        names = [str(name) for name in mdata.mod]
        if modality not in names:
            raise ValueError(
                f"The session reads the modality {modality!r}, which this MuData does not have. "
                f"Its modalities: {names}."
            )
        self.n_obs = int(mdata.n_obs)
        self.obs_names = mdata.obs_names
        self.uns = mdata.uns
        self.isbacked = bool(getattr(mdata, "isbacked", False))
        rows = {name: _placement(mdata, name) for name in names}

        obs = mdata.obs.copy()
        for name in names:
            frame = mdata.mod[name].obs
            for column in frame.columns:
                key = f"{name}:{column}"
                if key not in obs.columns:
                    obs[key] = pd.Series(_placed_column(frame[column], rows[name]), index=obs.index)
        self.obs = obs

        obsm: dict[str, Any] = dict(mdata.obsm.items())
        obsp: dict[str, Any] = dict(mdata.obsp.items())
        for name in names:
            modality_data = mdata.mod[name]
            for key, value in modality_data.obsm.items():
                if isinstance(value, pd.DataFrame):
                    continue
                if np.dtype(getattr(value, "dtype", object)).kind in "iuf":
                    obsm[f"{name}:{key}"] = _PlacedRows(value, rows[name])
            for key in modality_data.obsp:
                obsp[f"{name}:{key}"] = modality_data.obsp[key]
        self.obsm = obsm
        self.obsp = obsp

        chosen = mdata.mod[modality]
        own = rows[modality]
        self.var = chosen.var
        self.n_vars = int(chosen.n_vars)
        self.X = None if chosen.X is None else _PlacedMatrix(chosen.X, own)
        self.raw = None if chosen.raw is None else _PlacedRaw(chosen.raw, own)
        self.layers: Mapping[str, Any] = {
            key: _PlacedMatrix(value, own) for key, value in chosen.layers.items()
        }

