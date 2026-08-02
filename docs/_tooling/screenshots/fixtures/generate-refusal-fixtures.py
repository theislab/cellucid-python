#!/usr/bin/env python3
"""Regenerate the on-disk fixtures the data-loading refusal captures select.

Run from anywhere with the repository's own environment::

    python docs/_tooling/screenshots/fixtures/generate-refusal-fixtures.py

Why a committed file rather than one built at capture time. ``capture.mjs`` has
no dependencies of its own -- Playwright is resolved out of the web repository
at run time and nothing else is needed -- and making one scenario shell out to a
Python interpreter would end that. The generator is committed beside its output
instead, so the fixture is reproducible without giving the Node tool a Python
dependency it would carry for every other scenario too.

The other refusal fixture needs no generator and is built by the scenario that
uses it: a CSV renamed ``.h5ad`` is a string.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# A current-schema AnnData whose ``obsm`` holds exactly one key, and that key is
# not an embedding Cellucid accepts. This is the ordinary version of the
# mistake: PCA was computed and UMAP was not, so the file looks finished and is
# refused for the one thing it is missing. The refusal quotes the ``obsm`` keys
# that *are* present, so ``X_pca`` is what the figure's caption names.
N_CELLS = 120
N_GENES = 6
SEED = 0


def build() -> ad.AnnData:
    rng = np.random.default_rng(SEED)
    counts = rng.poisson(2.0, size=(N_CELLS, N_GENES)).astype(np.float32)
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(
                ["Alpha", "Beta", "Delta"] * (N_CELLS // 3)
            ),
            "total_counts": counts.sum(axis=1),
        },
        index=[f"cell-{index:04d}" for index in range(N_CELLS)],
    )
    var = pd.DataFrame(index=[f"gene-{index}" for index in range(N_GENES)])
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.obsm["X_pca"] = rng.normal(size=(N_CELLS, 2)).astype(np.float32)
    return adata


def main() -> None:
    target = Path(__file__).resolve().parent / "pca-only.h5ad"
    adata = build()
    assert list(adata.obsm) == ["X_pca"], list(adata.obsm)
    adata.write_h5ad(target, compression="gzip")
    print(f"wrote {target} ({target.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
