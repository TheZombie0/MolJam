"""Snapshot tests for score_database() top-level API."""

import pytest
from moljam import score_database, MoleculeDBScorer


def test_score_database_returns_3_tuple(tiny_db_path):
    result = score_database(
        tiny_db_path,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
        use_parallel=False,
    )
    assert isinstance(result, tuple)
    assert len(result) == 3


def test_score_database_types(tiny_db_path):
    final_score, report, scorer = score_database(
        tiny_db_path,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
        use_parallel=False,
    )
    assert isinstance(final_score, float)
    assert isinstance(report, str)
    assert isinstance(scorer, MoleculeDBScorer)


def test_final_score_in_range(tiny_db_path):
    final_score, _, _ = score_database(
        tiny_db_path,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
        use_parallel=False,
    )
    assert 0 <= final_score <= 100


def test_report_contains_header(tiny_db_path):
    _, report, _ = score_database(
        tiny_db_path,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
        use_parallel=False,
    )
    assert "# Molecular Database Quality Assessment Report" in report


def test_scores_dict_has_expected_keys(tiny_db_path):
    _, _, scorer = score_database(
        tiny_db_path,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
        use_parallel=False,
    )
    expected_top_keys = [
        "Structural Integrity",
        "Data Quality",
        "Chemical Space Coverage",
        "Data Distribution",
        "Total Score",
        "Normalized Score",
        "Final Adjusted Score",
    ]
    for key in expected_top_keys:
        assert key in scorer.scores, f"Missing key: {key}"


def test_score_database_with_experimental_info(tiny_db_path):
    final_score, report, scorer = score_database(
        tiny_db_path,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
        experimental_method_cols=["exp_method"],
        time_col="time_point",
        use_parallel=False,
        include_experimental_info=True,
    )
    assert isinstance(final_score, float)
    assert 0 <= final_score <= 100
    assert "Experimental Information Quality" in scorer.scores
