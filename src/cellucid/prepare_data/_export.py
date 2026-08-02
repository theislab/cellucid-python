"""
The public export entry point.

:func:`prepare` owns everything that happens around a generation rather than
inside it: it validates the arguments a caller supplies, refuses an output
directory that is not safe to replace whole, takes the writer lock, resolves
any interrupted previous run, and publishes -- or rolls back -- the generation
the writer builds.
"""

from collections.abc import Sequence
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
from scipy import sparse

from .._contracts import (
    _require_dataset_id,
    _require_dataset_name,
    _require_export_output_directory,
    _require_native_boolean,
    _require_optional_native_string,
    _require_positive_native_integer,
    _resolve_created_at,
)
from ._generation import _prepare_generation
from ._locking import _exclusive_export_generation
from ._transaction import (
    _begin_export_transaction,
    _publish_export_generation,
    _recover_export_transaction,
)


def prepare(
    latent_space: np.ndarray | sparse.spmatrix | None = None,
    obs: pd.DataFrame | None = None,
    var: pd.DataFrame | None = None,
    gene_expression: np.ndarray | sparse.spmatrix | None = None,
    var_gene_id_column: str | None = None,
    gene_identifiers: Sequence[str] | None = None,
    connectivities: np.ndarray | sparse.spmatrix | None = None,
    out_dir: Path | str = "./exports",
    obs_keys: Sequence[str] | None = None,
    centroid_outlier_quantile: float | None = 0.95,
    centroid_min_points: int = 10,
    force: bool = False,
    var_quantization: int | None = None,
    obs_continuous_quantization: int | None = None,
    compression: int | None = None,
    *,
    obs_categorical_dtype: Literal["uint8", "uint16"],
    dataset_name: str,
    dataset_id: str,
    created_at: str | None = None,
    dataset_description: str | None = None,
    source_name: str | None = None,
    source_url: str | None = None,
    source_citation: str | None = None,
    X_umap_1d: np.ndarray | None = None,
    X_umap_2d: np.ndarray | None = None,
    X_umap_3d: np.ndarray | None = None,
    vector_fields: dict[str, np.ndarray | sparse.spmatrix] | None = None,
    vector_field_default: str | None = None,
) -> None:
    """Build and atomically publish one complete canonical Cellucid export generation.

    ``out_dir`` defaults to ``./exports`` and is resolved against the current
    working directory at call time, so a later :func:`os.chdir` changes where a
    subsequent call writes. A leading ``~`` is expanded to the user's home
    directory.

    ``out_dir`` must name a directory dedicated to this one export, because
    publishing replaces the whole directory and ``force=True`` removes
    everything the previous one held. The filesystem root, the current working
    directory, the home directory, and the directory holding every home are
    refused before any path is created, locked, written, or removed, whether
    they are named directly, through ``~``, or through a symbolic link.

    ``created_at`` defaults to the current UTC time. Reproducible builders can
    pass an exact ``YYYY-MM-DDTHH:MM:SSZ`` UTC timestamp; it is validated and
    preserved byte-for-byte in ``dataset_identity.json``.

    A quantized continuous payload is published with a trailing
    ``minValue``/``maxValue``. A gene, continuous obs field, or generated
    outlier-quantile set that carries one value is published as a constant
    field: ``minValue == maxValue`` and every code ``0``, which the viewer
    decodes back to that exact value rather than to an approximation of it. A
    gene detected in no cell of a subset is the usual cause, and it is
    published like any other gene.

    A generated outlier-quantile set in which *every* quantile is missing is
    the one continuous payload with nothing to publish, because no category
    holds ``centroid_min_points`` cells. Those sets are all checked before any
    output path is created and a single failure names every one of them: lower
    ``centroid_min_points``, drop the field from ``obs_keys``, or pass
    ``obs_continuous_quantization=None`` to write float32 values and no
    manifest bounds.
    """
    force = _require_native_boolean(force, label="force")
    dataset_name = _require_dataset_name(
        dataset_name,
        label="dataset_name",
    )
    dataset_id = _require_dataset_id(
        dataset_id,
        label="dataset_id",
    )
    centroid_min_points = _require_positive_native_integer(
        centroid_min_points,
        label="centroid_min_points",
    )
    dataset_description = _require_optional_native_string(
        dataset_description,
        label="dataset_description",
        allow_empty=True,
    )
    source_name = _require_optional_native_string(
        source_name,
        label="source_name",
        allow_empty=False,
    )
    source_url = _require_optional_native_string(
        source_url,
        label="source_url",
        allow_empty=False,
    )
    source_citation = _require_optional_native_string(
        source_citation,
        label="source_citation",
        allow_empty=False,
    )
    created_at = _resolve_created_at(created_at)
    if source_name is None and (source_url is not None or source_citation is not None):
        raise ValueError(
            "source_name is required whenever source_url or source_citation is provided."
        )

    target_dir = _require_export_output_directory(out_dir, label="out_dir")
    target_dir.parent.mkdir(parents=True, exist_ok=True)
    with _exclusive_export_generation(target_dir):
        _recover_export_transaction(target_dir)
        if target_dir.is_symlink():
            raise ValueError(f"out_dir must not be a symbolic link: {target_dir}")
        if target_dir.exists():
            if not target_dir.is_dir():
                raise NotADirectoryError(f"Export path exists but is not a directory: {target_dir}")
            if not force and any(target_dir.iterdir()):
                raise FileExistsError(
                    f"Refusing to replace non-empty export directory {target_dir}. "
                    "Pass force=True to publish a complete replacement generation."
                )

        try:
            transaction_id, had_target, staging_dir, _backup_dir = _begin_export_transaction(
                target_dir
            )
            _prepare_generation(
                latent_space=latent_space,
                obs=obs,
                var=var,
                gene_expression=gene_expression,
                var_gene_id_column=var_gene_id_column,
                gene_identifiers=gene_identifiers,
                connectivities=connectivities,
                out_dir=staging_dir,
                obs_keys=obs_keys,
                centroid_outlier_quantile=centroid_outlier_quantile,
                centroid_min_points=centroid_min_points,
                var_quantization=var_quantization,
                obs_continuous_quantization=obs_continuous_quantization,
                obs_categorical_dtype=obs_categorical_dtype,
                compression=compression,
                _published_out_dir=target_dir,
                dataset_name=dataset_name,
                dataset_id=dataset_id,
                created_at=created_at,
                dataset_description=dataset_description,
                source_name=source_name,
                source_url=source_url,
                source_citation=source_citation,
                X_umap_1d=X_umap_1d,
                X_umap_2d=X_umap_2d,
                X_umap_3d=X_umap_3d,
                vector_fields=vector_fields,
                vector_field_default=vector_field_default,
            )
            _publish_export_generation(
                staging_dir,
                target_dir,
                transaction_id=transaction_id,
                had_target=had_target,
            )
        except BaseException:
            try:
                _recover_export_transaction(target_dir)
            except BaseException as cleanup_error:
                raise RuntimeError(
                    f"Failed to recover rejected export transaction for {target_dir}."
                ) from cleanup_error
            raise
