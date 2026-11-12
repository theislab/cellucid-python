"""Integer columns that hold identifiers: whole numbers past 2^53.

No double holds an integer past 2^53 exactly, and such integers are
identifiers (a database id, a barcode hashed to an integer), never
measurements. The web app therefore reads an obs column of 64-bit integers
holding one as its exact decimal text, a text column (``values.js``
``readExactIntegers``; a categorical's integer labels likewise), and says
why. What this package hands the viewer or writes back from it reads such a
column the same way: the labels a session names are those texts.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

#: The largest integer a double holds exactly together with its neighbors (2^53 - 1):
#: JavaScript's Number.MAX_SAFE_INTEGER. The viewer reads numbers as doubles, so a
#: larger integer would reach it rounded; such a column is carried as text instead.
SAFE_INTEGER_MAX = 2**53 - 1

#: Why such a column is text, as the viewer says it (``UNSAFE_INTEGERS_AS_TEXT``).
UNSAFE_INTEGERS_AS_TEXT = (
    "Read as text: it holds whole numbers past 2^53, which are identifiers, not measurements, "
    "and no number holds them exactly."
)


def holds_unsafe_integers(values: Any) -> bool:
    """Whether an integer column (NumPy, or pandas' nullable ``Int64``/``UInt64``)
    holds a value past the safe integers (beyond +/-(2^53 - 1)); missing values
    are not data. Any other column holds none."""
    dtype = getattr(values, "dtype", None)
    if dtype is None:
        values = np.asarray(values)
        dtype = values.dtype
    numpy_dtype = getattr(dtype, "numpy_dtype", dtype)
    if not isinstance(numpy_dtype, np.dtype) or numpy_dtype.kind not in "iu":
        return False
    if isinstance(values, pd.Series):
        values = values.array
    if isinstance(values, pd.arrays.IntegerArray):
        present = values.to_numpy(dtype=numpy_dtype, na_value=0)[~np.asarray(values.isna())]
    else:
        present = np.asarray(values)
    if present.size == 0:
        return False
    if numpy_dtype.kind == "u":
        return bool(present.max() > SAFE_INTEGER_MAX)
    return bool(present.max() > SAFE_INTEGER_MAX or present.min() < -SAFE_INTEGER_MAX)


def as_exact_text(series: pd.Series) -> pd.Series:
    """An integer column as its exact decimal texts (None where missing), same index."""
    values = series.array if isinstance(series.array, pd.arrays.IntegerArray) else np.asarray(series)
    if isinstance(values, pd.arrays.IntegerArray):
        missing = np.asarray(values.isna())
        numbers = values.to_numpy(dtype=values.dtype.numpy_dtype, na_value=0)
    else:
        missing = np.zeros(len(values), dtype=bool)
        numbers = values
    text = np.array([None if gone else str(int(number)) for number, gone in zip(numbers, missing, strict=True)], dtype=object)
    return pd.Series(text, index=series.index, name=series.name)
