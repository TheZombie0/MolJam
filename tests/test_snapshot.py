"""Tests for ScorerSnapshot."""

import pytest

from moljam.scoring.snapshot import ScorerSnapshot


def test_from_scorer_creates_snapshot(scorer):
    scorer.run_all_checks()
    snap = scorer.to_snapshot()
    assert isinstance(snap, ScorerSnapshot)


def test_snapshot_has_all_viz_attributes(scorer):
    scorer.run_all_checks()
    snap = scorer.to_snapshot()
    assert snap.scores is scorer.scores
    assert snap.analysis_results is scorer.analysis_results
    assert snap.num_molecules == scorer.num_molecules
    assert snap.smiles_col == scorer.smiles_col
    assert snap.activity_cols == list(scorer.activity_cols)
    assert snap.class_cols == list(scorer.class_cols)
    assert snap.include_experimental_info == scorer.include_experimental_info


def test_snapshot_matches_scorer_validation_attrs(scorer):
    scorer.run_all_checks()
    snap = scorer.to_snapshot()
    assert snap.invalid_rate == scorer.invalid_rate
    assert snap.non_canonical_rate == scorer.non_canonical_rate
    assert snap.valid_rate == scorer.valid_rate
    assert snap.invalid_indices == list(scorer.invalid_indices)
    assert snap.non_canonical_indices == list(scorer.non_canonical_indices)
    assert len(snap.valid_df) == len(scorer.valid_df)


def test_snapshot_df_not_empty(scorer):
    scorer.run_all_checks()
    snap = scorer.to_snapshot()
    assert len(snap.df) > 0
    assert len(snap.valid_df) > 0
