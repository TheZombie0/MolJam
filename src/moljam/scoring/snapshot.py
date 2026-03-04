"""ScorerSnapshot: the contract between Scorer and Viz layer."""

from dataclasses import dataclass, field
from typing import Any, Dict, List

import pandas as pd


@dataclass
class ScorerSnapshot:
    """Viz layer reads scorer data exclusively through this contract."""

    scores: Dict[str, Any] = field(default_factory=dict)
    analysis_results: Dict[str, Any] = field(default_factory=dict)
    num_molecules: int = 0
    df: pd.DataFrame = field(default_factory=pd.DataFrame)
    valid_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    invalid_indices: List[int] = field(default_factory=list)
    non_canonical_indices: List[int] = field(default_factory=list)
    invalid_rate: float = 0.0
    non_canonical_rate: float = 0.0
    valid_rate: float = 0.0
    smiles_col: str = "smiles"
    activity_cols: List[str] = field(default_factory=list)
    class_cols: List[str] = field(default_factory=list)
    include_experimental_info: bool = False

    @classmethod
    def from_scorer(cls, scorer) -> "ScorerSnapshot":
        """Create a snapshot from an active scorer object."""
        return cls(
            scores=scorer.scores,
            analysis_results=scorer.analysis_results,
            num_molecules=scorer.num_molecules,
            df=scorer.df,
            valid_df=getattr(scorer, "valid_df", pd.DataFrame()),
            invalid_indices=list(getattr(scorer, "invalid_indices", [])),
            non_canonical_indices=list(getattr(scorer, "non_canonical_indices", [])),
            invalid_rate=getattr(scorer, "invalid_rate", 0.0),
            non_canonical_rate=getattr(scorer, "non_canonical_rate", 0.0),
            valid_rate=getattr(scorer, "valid_rate", 0.0),
            smiles_col=scorer.smiles_col,
            activity_cols=list(scorer.activity_cols),
            class_cols=list(scorer.class_cols),
            include_experimental_info=scorer.include_experimental_info,
        )
