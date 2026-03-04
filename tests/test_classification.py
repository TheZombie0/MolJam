"""Tests for HybridColumnClassifier."""

import pandas as pd
import pytest

from moljam.classification.hybrid import HybridColumnClassifier


@pytest.fixture
def classifier():
    return HybridColumnClassifier()


def test_classify_column_returns_dict(classifier):
    col_data = pd.Series([1, 0, 1, 0, 1, 0, 1, 0])
    result = classifier.classify_column("activity", col_data)
    assert isinstance(result, dict)


def test_classify_column_has_required_keys(classifier):
    col_data = pd.Series([1, 0, 1, 0, 1, 0, 1, 0])
    result = classifier.classify_column("activity", col_data)
    assert "category" in result
    assert "confidence" in result


def test_classify_column_category_values(classifier):
    col_data = pd.Series([1, 0, 1, 0, 1, 0, 1, 0])
    result = classifier.classify_column("activity", col_data)
    assert result["category"] in ("useful", "excluded", "unknown")


def test_classify_columns_returns_dict(classifier, tiny_db_df):
    result = classifier.classify_columns(tiny_db_df, smiles_col="smiles")
    assert isinstance(result, dict)
    assert "useful" in result
    assert "excluded" in result
    assert "unknown" in result
