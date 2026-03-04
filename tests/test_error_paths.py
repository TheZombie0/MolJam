"""Tests for error handling paths."""

import pytest
from moljam import score_database


def test_missing_smiles_column(tiny_db_path):
    result = score_database(
        tiny_db_path,
        smiles_col="nonexistent_column",
        use_parallel=False,
    )
    assert result == (None, None, None)


def test_missing_file_path():
    result = score_database(
        "/nonexistent/path/to/file.csv",
        smiles_col="smiles",
        use_parallel=False,
    )
    assert result == (None, None, None)
