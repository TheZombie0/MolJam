"""Tests for individual MoleculeDBScorer check methods."""

import pandas as pd
import pytest


class TestValidateSmiles:
    def test_returns_float(self, scorer):
        score = scorer.validate_smiles()
        assert isinstance(score, float)

    def test_score_in_range(self, scorer):
        score = scorer.validate_smiles()
        assert 0 <= score <= 10

    def test_populates_scores_dict(self, scorer):
        scorer.validate_smiles()
        assert scorer.scores["Structural Integrity"]["Valid SMILES"] > 0

    def test_adds_to_completed_checks(self, scorer):
        scorer.validate_smiles()
        assert "validate_smiles" in scorer.completed_checks

    def test_sets_valid_df(self, scorer):
        scorer.validate_smiles()
        assert hasattr(scorer, "valid_df")
        assert not scorer.valid_df.empty

    def test_detects_invalid_smiles(self, scorer):
        scorer.validate_smiles()
        # tiny_db has 3 invalid SMILES
        assert len(scorer.invalid_indices) >= 1


class TestRepresentationConsistency:
    def test_returns_numeric(self, scorer):
        scorer.validate_smiles()
        score = scorer.check_representation_consistency()
        assert isinstance(score, (int, float))

    def test_score_in_range(self, scorer):
        scorer.validate_smiles()
        score = scorer.check_representation_consistency()
        assert 0 <= score <= 10

    def test_populates_scores_dict(self, scorer):
        scorer.validate_smiles()
        scorer.check_representation_consistency()
        assert scorer.scores["Structural Integrity"]["Representation Consistency"] >= 0

    def test_adds_to_completed_checks(self, scorer):
        scorer.validate_smiles()
        scorer.check_representation_consistency()
        assert "check_representation_consistency" in scorer.completed_checks


class TestStereochemistry:
    def test_returns_float(self, scorer):
        scorer.validate_smiles()
        score = scorer.check_stereochemistry()
        assert isinstance(score, float)

    def test_score_in_range(self, scorer):
        scorer.validate_smiles()
        score = scorer.check_stereochemistry()
        assert 0 <= score <= 10

    def test_adds_to_completed_checks(self, scorer):
        scorer.validate_smiles()
        scorer.check_stereochemistry()
        assert "check_stereochemistry" in scorer.completed_checks


class TestLabelConsistency:
    def test_returns_float(self, scorer):
        scorer.validate_smiles()
        score = scorer.check_label_consistency()
        assert isinstance(score, float)

    def test_score_in_range(self, scorer):
        scorer.validate_smiles()
        score = scorer.check_label_consistency()
        assert 0 <= score <= 10

    def test_adds_to_completed_checks(self, scorer):
        scorer.validate_smiles()
        scorer.check_label_consistency()
        assert "check_label_consistency" in scorer.completed_checks


class TestDataConsistency:
    def test_returns_float(self, scorer):
        scorer.validate_smiles()
        score = scorer.check_data_consistency_and_reliability()
        assert isinstance(score, float)

    def test_score_in_range(self, scorer):
        scorer.validate_smiles()
        score = scorer.check_data_consistency_and_reliability()
        assert 0 <= score <= 10

    def test_adds_to_completed_checks(self, scorer):
        scorer.validate_smiles()
        scorer.check_data_consistency_and_reliability()
        assert "check_data_consistency_and_reliability" in scorer.completed_checks


class TestChemicalDiversity:
    def test_returns_float(self, scorer):
        scorer.validate_smiles()
        score = scorer.analyze_chemical_diversity()
        assert isinstance(score, float)

    def test_score_in_range(self, scorer):
        scorer.validate_smiles()
        score = scorer.analyze_chemical_diversity()
        assert 0 <= score <= 10

    def test_adds_to_completed_checks(self, scorer):
        scorer.validate_smiles()
        scorer.analyze_chemical_diversity()
        assert "analyze_chemical_diversity" in scorer.completed_checks


class TestDruglikeness:
    def test_returns_float(self, scorer):
        scorer.validate_smiles()
        score = scorer.analyze_druglikeness()
        assert isinstance(score, float)

    def test_score_in_range(self, scorer):
        scorer.validate_smiles()
        score = scorer.analyze_druglikeness()
        assert 0 <= score <= 10

    def test_adds_to_completed_checks(self, scorer):
        scorer.validate_smiles()
        scorer.analyze_druglikeness()
        assert "analyze_druglikeness" in scorer.completed_checks


class TestDataSize:
    def test_returns_numeric(self, scorer):
        score = scorer.check_data_size()
        assert isinstance(score, (int, float))

    def test_score_in_range(self, scorer):
        score = scorer.check_data_size()
        assert 0 <= score <= 10

    def test_adds_to_completed_checks(self, scorer):
        scorer.check_data_size()
        assert "check_data_size" in scorer.completed_checks


class TestDataBalance:
    def test_returns_float(self, scorer):
        scorer.validate_smiles()
        score = scorer.analyze_data_balance_and_distribution()
        assert isinstance(score, float)

    def test_score_in_range(self, scorer):
        scorer.validate_smiles()
        score = scorer.analyze_data_balance_and_distribution()
        assert 0 <= score <= 10

    def test_adds_to_completed_checks(self, scorer):
        scorer.validate_smiles()
        scorer.analyze_data_balance_and_distribution()
        assert "analyze_data_balance_and_distribution" in scorer.completed_checks


class TestRunAllChecks:
    def test_returns_float(self, scorer):
        score = scorer.run_all_checks()
        assert isinstance(score, float)

    def test_score_in_range(self, scorer):
        score = scorer.run_all_checks()
        assert 0 <= score <= 100

    def test_all_checks_completed(self, scorer):
        scorer.run_all_checks()
        expected = {
            "validate_smiles",
            "check_representation_consistency",
            "check_stereochemistry",
            "check_label_consistency",
            "check_data_consistency_and_reliability",
            "analyze_chemical_diversity",
            "analyze_druglikeness",
            "check_data_size",
            "analyze_data_balance_and_distribution",
        }
        assert expected.issubset(scorer.completed_checks)


class TestRunAllChecksWithExperimental:
    def test_score_in_range(self, scorer_with_experimental):
        score = scorer_with_experimental.run_all_checks()
        assert 0 <= score <= 100

    def test_experimental_scores_populated(self, scorer_with_experimental):
        scorer_with_experimental.run_all_checks()
        exp_scores = scorer_with_experimental.scores["Experimental Information Quality"]
        assert "Time Label Availability" in exp_scores
        assert "Useful Column Quality" in exp_scores
        assert "Classification Confidence" in exp_scores
        assert "Type Diversity" in exp_scores


class TestCleanDatabase:
    def test_returns_tuple(self, scorer):
        scorer.run_all_checks()
        result = scorer.clean_database(verbose=False)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_returns_dataframe_and_dict(self, scorer):
        scorer.run_all_checks()
        cleaned_df, report = scorer.clean_database(verbose=False)
        assert isinstance(cleaned_df, pd.DataFrame)
        assert isinstance(report, dict)

    def test_cleaned_df_smaller_or_equal(self, scorer):
        scorer.run_all_checks()
        cleaned_df, _ = scorer.clean_database(verbose=False)
        assert len(cleaned_df) <= len(scorer.df)
