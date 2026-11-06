"""
Which cells a session's state at one step says something about.

Every value of a session names cells by their original position (their row
in the order the dataset gives them), whatever cells the viewer had open. A
viewer that had a subsample open knew the state of those cells only: a cell
it did not open is unknown -- ``<NA>`` in a highlight column, missing in a
field's -- never "not selected". Two records say which cells were open:

- the manifest's ``cells``: the subsample the session was saved on (the
  cells its fingerprint counts);
- the ``cells`` state domain at the step applied: ``open`` is null when every
  cell was open, else the subsample's specification (``data/subsample.js``),
  which chooses the same cells again (:func:`cellucid.subsample_cells`).

A step's state is applied on the cells open in both: a step taken before
the person opened a share, or after she opened every cell again, is applied
on the cells it was taken on, and never beyond the cells the session was
saved on.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .._plural import plural
from ..session import _validate_subsample_spec
from ._primitives import _require_exact_keys
from ._schema import CELLS_KEYS


def _spec_without_count(spec: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in spec.items() if key != "nObs"}


def _chosen_again(spec: dict[str, Any], shown: Any, n_obs: int, step: str) -> np.ndarray:
    """The cells a subsample specification chooses on ``shown``, as the web app chose them."""
    from ..subsample import subsample_cells

    method = spec["method"]
    try:
        if method == "uniform":
            return subsample_cells(n_obs, spec["k"], method="uniform", seed=spec["seed"])
        if method == "blocks":
            return subsample_cells(n_obs, spec["k"], method="blocks", seed=spec["seed"], block_size=spec["blockSize"])
        if method == "stratified":
            field = spec["field"]
            if field not in shown.obs.columns:
                raise KeyError(f"obs has no column {field!r}")
            return subsample_cells(
                n_obs, spec["k"], method="stratified", seed=spec["seed"], groupby=shown.obs[field],
                allocation=spec["allocation"], min_per_group=spec["minPerCategory"],
            )
        embedding = spec["embedding"]
        if embedding not in shown.obsm:
            raise KeyError(f"obsm has no key {embedding!r}")
        return subsample_cells(
            n_obs, spec["k"], method="spatial", seed=spec["seed"], basis=np.asarray(shown.obsm[embedding]),
            columns=spec["columns"], grid_size=spec["gridSize"],
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            f"The step {step} was taken on a {method} subsample of {spec['k']} of the {n_obs} "
            f"{plural(n_obs, 'cell')}, which cannot be chosen again on this object ({error}). "
            "Apply a step taken on other cells, or the object the session was saved on."
        ) from error


def known_cells(
    *,
    saved_on: dict[str, Any] | None,
    cells_value: Any,
    shown: Any,
    n_obs: int,
    step: str,
) -> np.ndarray | None:
    """The cells (ascending original positions) a step's state knows; None for every cell.

    ``saved_on`` is the session's decoded manifest ``cells`` (None: saved on
    every cell); ``cells_value`` the ``cells`` domain's value at the step;
    ``shown`` the object as the viewer showed it (an AnnData, or a MuData's
    :class:`~cellucid.anndata_session._mudata.MuDataSessionView`) with
    ``n_obs`` cells.
    """
    value = _require_exact_keys(cells_value, CELLS_KEYS, label=f"cells value at step {step}")
    saved = None if saved_on is None else np.asarray(saved_on["kept"], dtype=np.int64)
    if value["open"] is None:
        return saved
    open_spec = value["open"]
    if not isinstance(open_spec, dict) or "nObs" in open_spec:
        raise ValueError(f"This session's cells value at step {step} names its open cells in a form this version does not read.")
    spec = _validate_subsample_spec({**open_spec, "nObs": n_obs}, f"Session cells value at step {step}")
    if saved_on is not None and _spec_without_count(spec) == _spec_without_count(saved_on["spec"]):
        return saved
    chosen = np.asarray(_chosen_again(spec, shown, n_obs, step), dtype=np.int64)
    return chosen if saved is None else np.intersect1d(saved, chosen, assume_unique=True)
