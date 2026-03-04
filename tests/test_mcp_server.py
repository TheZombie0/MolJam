"""Tests for MCP Server tools (direct function calls, no server startup)."""

import json

import pytest


def test_mcp_server_importable():
    """Module imports and FastMCP instance is accessible."""
    from moljam.mcp_server import mcp

    assert mcp.name == "MolJam"


# ── score_database ──────────────────────────────────────────────────────────


def test_score_database_success(tiny_db_path):
    from moljam.mcp_server import score_database

    result_str = score_database(
        db_path=tiny_db_path,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
    )
    result = json.loads(result_str)
    assert result["status"] == "success"
    assert 0 <= result["final_score"] <= 100
    assert "scores" in result
    assert "dataset_info" in result


def test_score_database_error():
    from moljam.mcp_server import score_database

    result_str = score_database(db_path="/nonexistent/path.csv")
    result = json.loads(result_str)
    assert result["status"] == "error"
    assert "error" in result


# ── classify_columns ────────────────────────────────────────────────────────


def test_classify_columns_success(tiny_db_path):
    from moljam.mcp_server import classify_columns

    result_str = classify_columns(db_path=tiny_db_path, smiles_col="smiles")
    result = json.loads(result_str)
    assert result["status"] == "success"
    assert "classifications" in result
    classifications = result["classifications"]
    assert "useful" in classifications
    assert "excluded" in classifications
    assert "unknown" in classifications


def test_classify_columns_error():
    from moljam.mcp_server import classify_columns

    result_str = classify_columns(db_path="/nonexistent/path.csv")
    result = json.loads(result_str)
    assert result["status"] == "error"
    assert "error" in result


# ── clean_database ──────────────────────────────────────────────────────────


def test_clean_database_success(tiny_db_path, tmp_path):
    from moljam.mcp_server import clean_database

    output = str(tmp_path / "cleaned.csv")
    result_str = clean_database(
        db_path=tiny_db_path,
        output_path=output,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
    )
    result = json.loads(result_str)
    assert result["status"] == "success"
    assert result["output_path"] == output
    assert "cleaning_report" in result
    assert "rows_before" in result
    assert "rows_after" in result
    assert result["rows_after"] <= result["rows_before"]

    # Verify output file was written
    import pandas as pd

    df = pd.read_csv(output)
    assert len(df) == result["rows_after"]


def test_clean_database_error():
    from moljam.mcp_server import clean_database

    result_str = clean_database(
        db_path="/nonexistent/path.csv",
        output_path="/tmp/wont_be_created.csv",
    )
    result = json.loads(result_str)
    assert result["status"] == "error"
    assert "error" in result


# ── no report_text in score_database ────────────────────────────────────────


def test_score_database_no_report_text(tiny_db_path):
    """score_database tool should NOT include verbose report_text."""
    from moljam.mcp_server import score_database

    result_str = score_database(
        db_path=tiny_db_path,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
    )
    result = json.loads(result_str)
    assert result["status"] == "success"
    assert "report_text" not in result
