import json

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid.prepare_data import prepare

_DEFAULT_VAR = object()


def _prepare_gene_export(
    tmp_path,
    *,
    var=_DEFAULT_VAR,
    gene_identifiers=None,
    var_gene_id_column=None,
    force=False,
    sparse_expression=False,
):
    if var is _DEFAULT_VAR:
        var = pd.DataFrame(index=["gene_a", "gene_b", "gene_c"])

    gene_expression = np.array(
        [
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ],
        dtype=np.float32,
    )
    if sparse_expression:
        gene_expression = sparse.csr_matrix(gene_expression)

    prepare(
        latent_space=np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [0.0, 1.0],
            ],
            dtype=np.float32,
        ),
        obs=pd.DataFrame(index=["cell_a", "cell_b", "cell_c"]),
        var=var,
        gene_expression=gene_expression,
        gene_identifiers=gene_identifiers,
        var_gene_id_column=var_gene_id_column,
        X_umap_2d=np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [0.0, 1.0],
            ],
            dtype=np.float32,
        ),
        out_dir=tmp_path,
        dataset_id="gene-subset",
        dataset_name="Gene subset",
        obs_categorical_dtype="uint16",
        force=force,
    )


@pytest.mark.parametrize("sparse_expression", [False, True])
def test_identity_gene_count_matches_exported_var_subset(
    tmp_path,
    sparse_expression,
):
    _prepare_gene_export(
        tmp_path,
        gene_identifiers=["gene_c", "gene_a"],
        sparse_expression=sparse_expression,
    )

    var_manifest = json.loads((tmp_path / "var_manifest.json").read_text(encoding="utf-8"))
    identity = json.loads((tmp_path / "dataset_identity.json").read_text(encoding="utf-8"))

    assert var_manifest["fields"] == [["gene_c"], ["gene_a"]]
    assert identity["stats"]["n_genes"] == len(var_manifest["fields"]) == 2


@pytest.mark.parametrize(
    ("var", "var_gene_id_column", "error_label"),
    [
        (
            pd.DataFrame(index=["gene_a", 2, "gene_c"]),
            None,
            "var index",
        ),
        (
            pd.DataFrame(
                {"gene_id": ["gene_a", 2, "gene_c"]},
                index=["row_a", "row_b", "row_c"],
            ),
            "gene_id",
            "var column 'gene_id'",
        ),
    ],
)
def test_var_identifiers_must_be_strings(
    tmp_path,
    var,
    var_gene_id_column,
    error_label,
):
    with pytest.raises(
        TypeError,
        match=rf"{error_label} identifier at position 1 must be a string",
    ):
        _prepare_gene_export(
            tmp_path,
            var=var,
            var_gene_id_column=var_gene_id_column,
        )


def test_requested_gene_identifiers_must_be_strings(tmp_path):
    with pytest.raises(
        TypeError,
        match="Requested gene identifier at position 1 must be a string",
    ):
        _prepare_gene_export(
            tmp_path,
            gene_identifiers=["gene_a", 2],
        )


def test_gene_identifier_string_scalar_is_not_treated_as_a_sequence(tmp_path):
    with pytest.raises(
        TypeError,
        match="gene_identifiers must be a sequence of strings",
    ):
        _prepare_gene_export(
            tmp_path,
            gene_identifiers="gene_a",
        )


def test_duplicate_var_identifiers_are_rejected(tmp_path):
    with pytest.raises(
        ValueError,
        match="Gene key 'gene_a' is duplicated",
    ):
        _prepare_gene_export(
            tmp_path,
            var=pd.DataFrame(index=["gene_a", "gene_a", "gene_c"]),
        )


def test_duplicate_requested_gene_identifiers_are_rejected(tmp_path):
    with pytest.raises(
        ValueError,
        match="Requested gene key 'gene_a' is duplicated",
    ):
        _prepare_gene_export(
            tmp_path,
            gene_identifiers=["gene_a", "gene_a"],
        )


def test_missing_requested_gene_identifier_is_rejected(tmp_path):
    with pytest.raises(
        KeyError,
        match="gene_identifiers contain identifiers not present in var",
    ):
        _prepare_gene_export(
            tmp_path,
            gene_identifiers=["gene_a", "missing_gene"],
        )


def test_none_selects_var_index_and_literal_index_selects_the_column(tmp_path):
    var = pd.DataFrame(
        {"index": ["column_a", "column_b", "column_c"]},
        index=["row_a", "row_b", "row_c"],
    )

    index_export = tmp_path / "var-index"
    _prepare_gene_export(index_export, var=var, var_gene_id_column=None)
    index_manifest = json.loads(
        (index_export / "var_manifest.json").read_text(encoding="utf-8")
    )
    assert index_manifest["var_gene_id_column"] is None
    assert index_manifest["fields"] == [["row_a"], ["row_b"], ["row_c"]]

    column_export = tmp_path / "literal-index-column"
    _prepare_gene_export(column_export, var=var, var_gene_id_column="index")
    column_manifest = json.loads(
        (column_export / "var_manifest.json").read_text(encoding="utf-8")
    )
    assert column_manifest["var_gene_id_column"] == "index"
    assert column_manifest["fields"] == [["column_a"], ["column_b"], ["column_c"]]


def test_existing_var_manifest_must_match_requested_subset(tmp_path):
    _prepare_gene_export(
        tmp_path,
        gene_identifiers=["gene_a", "gene_b"],
    )

    with pytest.raises(FileExistsError, match="non-empty export directory"):
        _prepare_gene_export(
            tmp_path,
            gene_identifiers=["gene_b", "gene_c"],
        )
