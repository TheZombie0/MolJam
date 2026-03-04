"""ScoringContext: shared state container replacing implicit self attributes."""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set

import pandas as pd

from .config import DEFAULT_CONFIG, ScoringConfig


@dataclass
class ScoringContext:
    """Shared state container for scoring checks.

    Replaces the implicit ``self.*`` attribute soup in the mixin-based design.
    New checks can use this instead of inheriting from MoleculeDBScorer.
    """

    df: pd.DataFrame = field(default_factory=pd.DataFrame)
    config: ScoringConfig = field(default_factory=lambda: DEFAULT_CONFIG)
    smiles_col: str = "smiles"
    activity_cols: List[str] = field(default_factory=list)
    class_cols: List[str] = field(default_factory=list)
    experimental_method_cols: List[str] = field(default_factory=list)
    id_col: Optional[str] = None
    name_col: Optional[str] = None
    time_col: Optional[str] = None
    use_parallel: bool = True
    include_experimental_info: bool = False
    num_molecules: int = 0

    # Populated during scoring
    scores: Dict[str, Any] = field(default_factory=dict)
    analysis_results: Dict[str, Any] = field(default_factory=dict)
    completed_checks: Set[str] = field(default_factory=set)

    # Set by validate_smiles
    valid_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    valid_mols: list = field(default_factory=list)
    invalid_indices: list = field(default_factory=list)
    non_canonical_indices: list = field(default_factory=list)
    valid_rate: float = 0.0
    invalid_rate: float = 0.0
    non_canonical_rate: float = 0.0
