"""Exact contract for the gene-name sets of ``var_manifest.json`` (compact_v2).

One expression payload per gene, any number of names for it. ``fields[i][1]``
is the default name of gene ``i`` -- the column ``var_gene_id_column`` names, or
``var.index`` -- and ``gene_name_sets`` describes where every name set came
from. Set 0 is the default set and carries ``names: null``; every column in
``var_gene_name_columns`` adds a set whose ``names`` lists one name per exported
gene, null wherever it equals the default name. In every set an empty value
takes that row's ``var.index`` name, and repeated names are made unique with
anndata's rule over the whole var, before any subset is taken. What was filled
and renamed is reported over the exported genes, in exported order.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cellucid.prepare_data import prepare
from cellucid.prepare_data._manifest import _gene_names_from_compact_manifest

N_CELLS = 3


def _export(
    out_dir: Path,
    var: pd.DataFrame,
    **overrides: object,
) -> dict:
    kwargs: dict[str, object] = {
        "latent_space": np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
        "obs": pd.DataFrame(index=["cell-0", "cell-1", "cell-2"]),
        "embeddings": {"X_umap": np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])},
        "var": var,
        "gene_expression": np.arange(N_CELLS * len(var), dtype=np.float32).reshape(
            N_CELLS, len(var)
        ),
        "out_dir": out_dir,
        "dataset_name": "Gene name sets",
        "dataset_id": "gene-name-sets",
        "created_at": "2026-01-01T00:00:00Z",
        "obs_categorical_dtype": "uint8",
        "centroid_min_points": 1,
    }
    kwargs.update(overrides)
    prepare(**kwargs)  # type: ignore[arg-type]
    return json.loads((out_dir / "var_manifest.json").read_text(encoding="utf-8"))


def _var() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_name": ["MATR3", "MATR3", None, "CD4", "CD8A"],
            "symbol": ["S1", "S1", "", "CD4", np.nan],
            "ensembl": ["ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5"],
        },
        index=["g1", "g2", "g3", "g4", "g5"],
    )


def test_the_var_manifest_has_the_exact_compact_v2_shape(tmp_path: Path) -> None:
    manifest = _export(tmp_path, _var())
    assert list(manifest) == [
        "_format",
        "n_points",
        "compression",
        "quantization",
        "_varSchema",
        "fields",
        "gene_name_sets",
    ]
    assert manifest["_format"] == "compact_v2"
    assert "var_gene_id_column" not in manifest
    assert manifest["gene_name_sets"] == [
        {"source": None, "filled_from_index": 0, "renamed": [], "names": None}
    ]
    assert list(manifest["gene_name_sets"][0]) == [
        "source",
        "filled_from_index",
        "renamed",
        "names",
    ]
    assert manifest["fields"] == [[0, "g1"], [1, "g2"], [2, "g3"], [3, "g4"], [4, "g5"]]


def test_alternative_sets_take_the_index_for_empty_values_and_are_made_unique(
    tmp_path: Path,
) -> None:
    manifest = _export(
        tmp_path,
        _var(),
        var_gene_id_column="gene_name",
        var_gene_name_columns=["symbol", "ensembl"],
    )

    # None takes the index name 'g3'; the second MATR3 is made unique.
    assert manifest["fields"] == [
        [0, "MATR3"],
        [1, "MATR3-1"],
        [2, "g3"],
        [3, "CD4"],
        [4, "CD8A"],
    ]
    assert manifest["gene_name_sets"] == [
        {
            "source": "gene_name",
            "filled_from_index": 1,
            "renamed": [["MATR3", "MATR3-1"]],
            "names": None,
        },
        {
            # "" and NaN both take the index; the second S1 is made unique.
            "source": "symbol",
            "filled_from_index": 2,
            "renamed": [["S1", "S1-1"]],
            # Null wherever the name equals the default name.
            "names": ["S1", "S1-1", None, None, "g5"],
        },
        {
            "source": "ensembl",
            "filled_from_index": 0,
            "renamed": [],
            "names": ["ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5"],
        },
    ]
    # Metadata only: one payload per gene, whatever the number of name sets.
    assert sorted(path.name for path in (tmp_path / "var").iterdir()) == [
        f"{index}.values.f32" for index in range(5)
    ]
    assert _gene_names_from_compact_manifest(manifest) == [field[1] for field in manifest["fields"]]


def test_uniqueness_spans_the_whole_var_and_reports_the_exported_genes_in_order(
    tmp_path: Path,
) -> None:
    var = pd.DataFrame(
        {"symbol": ["A", "B", "A", "", "A"]},
        index=["g1", "g2", "g3", "g4", "g5"],
    )
    manifest = _export(
        tmp_path,
        var,
        var_gene_name_columns=["symbol"],
        gene_identifiers=["g5", "g4", "g1"],
    )

    assert manifest["fields"] == [[0, "g5"], [1, "g4"], [2, "g1"]]
    # Over the whole var the symbols are A, B, A-1, g4, A-2: the fifth row is
    # 'A-2' even though the third row is not exported.
    assert manifest["gene_name_sets"][1] == {
        "source": "symbol",
        "filled_from_index": 1,
        "renamed": [["A", "A-2"]],
        "names": ["A-2", None, "A"],
    }


def test_a_default_column_that_repeats_a_name_is_addressed_by_its_unique_name(
    tmp_path: Path,
) -> None:
    manifest = _export(
        tmp_path,
        _var(),
        var_gene_id_column="gene_name",
        gene_identifiers=["MATR3-1", "g3"],
    )
    assert manifest["fields"] == [[0, "MATR3-1"], [1, "g3"]]
    assert manifest["gene_name_sets"] == [
        {
            "source": "gene_name",
            "filled_from_index": 1,
            "renamed": [["MATR3", "MATR3-1"]],
            "names": None,
        }
    ]
    values = np.fromfile(tmp_path / "var" / "0.values.f32", dtype="<f4")
    assert values.tolist() == [1.0, 6.0, 11.0]


def test_what_was_filled_and_renamed_is_printed(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _export(tmp_path, _var(), var_gene_id_column="gene_name", var_gene_name_columns=["symbol"])
    printed = capsys.readouterr().out
    assert (
        "Gene names from var column 'gene_name': 1 gene has no 'gene_name' value and uses "
        "its var index name. 1 duplicate name was made unique, e.g. 'MATR3' → 'MATR3-1'."
    ) in printed
    assert (
        "Gene names from var column 'symbol': 2 genes have no 'symbol' value and use their "
        "var index name. 1 duplicate name was made unique, e.g. 'S1' → 'S1-1'."
    ) in printed


@pytest.mark.parametrize(
    ("overrides", "error", "message"),
    [
        ({"var_gene_name_columns": ["missing"]}, KeyError, r"'missing' not found in var"),
        (
            {"var_gene_id_column": "gene_name", "var_gene_name_columns": ["gene_name"]},
            ValueError,
            r"lists 'gene_name', which var_gene_id_column already uses",
        ),
        (
            {"var_gene_name_columns": ["symbol", "symbol"]},
            ValueError,
            r"lists 'symbol' more than once",
        ),
        ({"var_gene_name_columns": "symbol"}, TypeError, r"not a string scalar"),
        ({"var_gene_name_columns": [7]}, TypeError, r"must be a native non-empty string"),
        ({"var_gene_name_columns": [""]}, ValueError, r"must be a native non-empty string"),
        (
            {"var_gene_name_columns": ["count"]},
            TypeError,
            r"var column 'count': floating-point values cannot name genes",
        ),
        ({"var_gene_name_columns": ["padded"]}, ValueError, r"displayed verbatim"),
        ({"var_gene_name_columns": [" symbol"]}, ValueError, r"displayed verbatim"),
        (
            {"var_gene_id_column": "count"},
            TypeError,
            r"var column 'count': floating-point values cannot name genes",
        ),
        ({"var_gene_id_column": "missing"}, KeyError, r"var_gene_id_column 'missing' not found"),
    ],
)
def test_invalid_gene_name_sets_are_rejected_before_output_mutation(
    tmp_path: Path,
    overrides: dict[str, object],
    error: type[Exception],
    message: str,
) -> None:
    var = _var()
    # Integers would name genes by their digits; a float column cannot.
    var["count"] = [1.5, 2.5, 3.5, 4.5, 5.5]
    var["padded"] = ["a", "b ", "c", "d", "e"]
    var[" symbol"] = ["a", "b", "c", "d", "e"]
    out_dir = tmp_path / "must-not-exist"
    with pytest.raises(error, match=message):
        _export(out_dir, var, **overrides)
    assert not out_dir.exists()


def test_an_undrawable_name_outside_the_export_is_not_checked(tmp_path: Path) -> None:
    var = _var()
    var["padded"] = ["a", "b", "c", "d", "e "]
    manifest = _export(
        tmp_path,
        var,
        var_gene_name_columns=["padded"],
        gene_identifiers=["g1", "g2"],
    )
    assert manifest["gene_name_sets"][1]["names"] == ["a", "b"]


def test_alternative_names_without_gene_expression_are_refused(tmp_path: Path) -> None:
    out_dir = tmp_path / "must-not-exist"
    with pytest.raises(ValueError, match=r"var_gene_name_columns.*no gene_expression"):
        prepare(
            latent_space=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
            obs=pd.DataFrame(index=["cell-0", "cell-1", "cell-2"]),
            embeddings={"X_umap": np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])},
            var=_var(),
            var_gene_name_columns=["symbol"],
            out_dir=out_dir,
            dataset_name="Gene name sets",
            dataset_id="gene-name-sets",
            obs_categorical_dtype="uint8",
        )
    assert not out_dir.exists()


def test_an_empty_index_name_cannot_stand_in_for_an_empty_value(tmp_path: Path) -> None:
    var = pd.DataFrame({"symbol": ["A", ""]}, index=["g1", ""])
    out_dir = tmp_path / "must-not-exist"
    with pytest.raises(
        ValueError,
        match=r"var column 'symbol': Name at position 1 is empty, and its index name is not a non-empty",
    ):
        _export(out_dir, var, var_gene_id_column="symbol")
    assert not out_dir.exists()


def _manifest(fields: list, sets: list) -> dict:
    return {
        "_format": "compact_v2",
        "n_points": 3,
        "compression": None,
        "quantization": None,
        "_varSchema": {},
        "fields": fields,
        "gene_name_sets": sets,
    }


def _set(source: object, names: object, *, filled: int = 0, renamed: list | None = None) -> dict:
    return {
        "source": source,
        "filled_from_index": filled,
        "renamed": [] if renamed is None else renamed,
        "names": names,
    }


@pytest.mark.parametrize(
    ("sets", "message"),
    [
        ([], r"gene_name_sets must be a non-empty list"),
        ([_set(None, ["a", "b"])], r"\[0\]\.names must be null"),
        ([_set(None, None, filled=1)], r"var index cannot fill names from itself"),
        ([_set(None, None), _set(None, ["x", "y"])], r"\[1\]\.source must name a var column"),
        ([_set(None, None), _set("c", ["x"])], r"one entry per gene"),
        ([_set(None, None), _set("c", ["A", "y"])], r"null wherever a name equals the default"),
        ([_set(None, None), _set("c", ["B", None])], r"gene key 'B' is duplicated"),
        ([_set(None, None), _set("c", ["x", ""])], r"non-empty strings or null"),
        ([_set("c", None), _set("c", ["x", "y"])], r"read each source once"),
        ([_set(None, None, renamed=[["A"]])], r"list of \[from, to\] name pairs"),
        ([_set(None, None, filled=-1)], r"filled_from_index must count exported genes"),
        ([{"source": None, "names": None}], r"must hold exactly the keys"),
    ],
)
def test_the_manifest_validator_rejects_inconsistent_name_sets(sets: list, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        _gene_names_from_compact_manifest(_manifest([[0, "A"], [1, "B"]], sets))


def test_the_manifest_validator_rejects_the_previous_format() -> None:
    manifest = _manifest([[0, "A"]], [_set(None, None)])
    manifest["_format"] = "compact_v1"
    with pytest.raises(ValueError, match=r"current compact_v2 format"):
        _gene_names_from_compact_manifest(manifest)
    manifest = _manifest([[0, "A"]], [_set(None, None)])
    manifest["var_gene_id_column"] = None
    with pytest.raises(ValueError, match=r"exactly the current compact_v2 fields"):
        _gene_names_from_compact_manifest(manifest)
