"""Tests for typed data models."""

import json

import pytest

from moljam.scoring.models import (
    CategoryScore,
    DatabaseScoringReport,
    ScoringResult,
)


@pytest.fixture
def sample_scores_dict():
    return {
        "Structural Integrity": {
            "Valid SMILES": 9.5,
            "Representation Consistency": 10.0,
            "Stereochemistry Completeness": 8.0,
            "Total": 27.5,
            "Normalized Total": 22.92,
        },
        "Data Quality": {
            "Label Consistency": 7.0,
            "Data Consistency and Reliability": 8.5,
            "Total": 15.5,
            "Normalized Total": 19.38,
        },
        "Chemical Space Coverage": {
            "Chemical Diversity": 6.0,
            "Drug-likeness": 7.5,
            "Total": 13.5,
            "Normalized Total": 16.88,
        },
        "Data Distribution": {
            "Data Size": 4.0,
            "Data Balance and Distribution": 6.0,
            "Total": 10.0,
            "Normalized Total": 12.5,
        },
        "Total Score": 66.5,
        "Normalized Score": 71.67,
        "Final Adjusted Score": 70.17,
        "Metric Penalties": {"Data Quality::Label Consistency": 0.0},
    }


def test_from_scores_dict_roundtrip(sample_scores_dict):
    result = ScoringResult.from_scores_dict(sample_scores_dict)
    restored = result.to_dict()
    # Check key fields match
    assert restored["Total Score"] == sample_scores_dict["Total Score"]
    assert restored["Normalized Score"] == sample_scores_dict["Normalized Score"]
    assert restored["Final Adjusted Score"] == sample_scores_dict["Final Adjusted Score"]
    assert restored["Structural Integrity"]["Valid SMILES"] == 9.5
    assert restored["Structural Integrity"]["Total"] == 27.5


def test_to_json_valid(sample_scores_dict):
    result = ScoringResult.from_scores_dict(sample_scores_dict)
    json_str = result.to_json()
    parsed = json.loads(json_str)
    assert isinstance(parsed, dict)
    assert "Structural Integrity" in parsed


def test_with_experimental_info():
    scores = {
        "Structural Integrity": {"Total": 0, "Normalized Total": 0},
        "Data Quality": {"Total": 0, "Normalized Total": 0},
        "Chemical Space Coverage": {"Total": 0, "Normalized Total": 0},
        "Data Distribution": {"Total": 0, "Normalized Total": 0},
        "Experimental Information Quality": {
            "Time Label Availability": 5.0,
            "Useful Column Quality": 3.0,
            "Classification Confidence": 7.0,
            "Type Diversity": 6.67,
            "Total": 21.67,
            "Normalized Total": 10.83,
        },
        "Total Score": 21.67,
        "Normalized Score": 10.83,
        "Final Adjusted Score": 10.83,
    }
    result = ScoringResult.from_scores_dict(scores)
    assert result.experimental_info_quality is not None
    assert result.experimental_info_quality.metrics["Time Label Availability"] == 5.0
    d = result.to_dict()
    assert "Experimental Information Quality" in d


def test_without_experimental_info(sample_scores_dict):
    result = ScoringResult.from_scores_dict(sample_scores_dict)
    assert result.experimental_info_quality is None
    d = result.to_dict()
    assert "Experimental Information Quality" not in d


def test_database_scoring_report_to_json():
    sr = ScoringResult(total_score=50.0, normalized_score=60.0, final_adjusted_score=55.0)
    report = DatabaseScoringReport(
        final_score=55.0,
        report_text="# Test",
        scores=sr,
        analysis_results={"key": {"nested": "value"}},
        dataset_info={"rows": 100},
    )
    json_str = report.to_json()
    parsed = json.loads(json_str)
    assert parsed["final_score"] == 55.0
    assert parsed["scores"]["Total Score"] == 50.0


def test_scorer_to_result(scorer):
    scorer.run_all_checks()
    result = scorer.to_result()
    assert isinstance(result, ScoringResult)
    assert result.final_adjusted_score == scorer.scores["Final Adjusted Score"]
    json_str = result.to_json()
    parsed = json.loads(json_str)
    assert isinstance(parsed, dict)
