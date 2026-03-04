"""Tests for MCP-friendly API."""

import json

import pytest

from moljam import score_database, score_database_json
from moljam.scoring.models import DatabaseScoringReport


def test_score_database_json_normal(tiny_db_path):
    report = score_database_json(
        tiny_db_path,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
        use_parallel=False,
    )
    assert isinstance(report, DatabaseScoringReport)
    assert report.final_score is not None
    assert 0 <= report.final_score <= 100
    # Verify JSON serialization
    json_str = report.to_json()
    parsed = json.loads(json_str)
    assert parsed["final_score"] is not None
    assert "scores" in parsed
    assert "analysis_results" in parsed


def test_score_database_json_failure():
    report = score_database_json(
        "/nonexistent/path.csv",
        smiles_col="smiles",
        use_parallel=False,
    )
    assert isinstance(report, DatabaseScoringReport)
    assert report.final_score is None
    assert report.scores is None


def test_score_database_unchanged(tiny_db_path):
    """Regression: original score_database still returns 3-tuple."""
    result = score_database(
        tiny_db_path,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
        use_parallel=False,
    )
    assert isinstance(result, tuple)
    assert len(result) == 3
    final_score, report, scorer = result
    assert isinstance(final_score, float)
    assert isinstance(report, str)


def test_mcp_score_database(tiny_db_path):
    from moljam.mcp import mcp_score_database

    result = mcp_score_database(
        tiny_db_path,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
        use_parallel=False,
    )
    assert isinstance(result, dict)
    assert "final_score" in result
    # Verify it's JSON-serializable
    json.dumps(result, default=str)


def test_mcp_classify_columns(tiny_db_path):
    from moljam.mcp import mcp_classify_columns

    result = mcp_classify_columns(tiny_db_path, smiles_col="smiles")
    assert isinstance(result, dict)
    assert result["status"] == "success"
    classifications = result["classifications"]
    assert "useful" in classifications
    assert "excluded" in classifications
    assert "unknown" in classifications


def test_mcp_classify_columns_error():
    from moljam.mcp import mcp_classify_columns

    result = mcp_classify_columns("/nonexistent/path.csv")
    assert isinstance(result, dict)
    assert result["status"] == "error"
    assert "error" in result


def test_mcp_clean_database(tiny_db_path):
    from moljam.mcp import mcp_clean_database

    result = mcp_clean_database(
        tiny_db_path,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
        use_parallel=False,
    )
    assert isinstance(result, dict)
    assert result["status"] == "success"
    assert "cleaning_report" in result


def test_mcp_clean_database_error():
    from moljam.mcp import mcp_clean_database

    result = mcp_clean_database("/nonexistent/path.csv")
    assert isinstance(result, dict)
    assert result["status"] == "error"
