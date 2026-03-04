import warnings

import numpy as np
import pandas as pd
from rdkit.rdBase import DisableLog

from .checks import (
    AnnotationChecksMixin,
    ChemicalSpaceChecksMixin,
    CleaningMixin,
    DataQualityChecksMixin,
    DistributionChecksMixin,
    OrchestrationMixin,
    ReportMixin,
    ScoringHelpersMixin,
    StructuralChecksMixin,
)
from .config import DEFAULT_CONFIG, ScoringConfig
from .models import ScoringResult
from .snapshot import ScorerSnapshot


class MoleculeDBScorer(
    ScoringHelpersMixin,
    StructuralChecksMixin,
    DataQualityChecksMixin,
    AnnotationChecksMixin,
    ChemicalSpaceChecksMixin,
    DistributionChecksMixin,
    OrchestrationMixin,
    ReportMixin,
    CleaningMixin,
):
    def __init__(self, df, smiles_col='smiles', activity_cols=None, class_cols=None,
             experimental_method_cols=None, id_col=None, name_col=None, time_col=None,
             use_parallel=True, include_experimental_info=False, config=None):
        """
        Initialize the scoring system

        Parameters:
            df: pandas DataFrame - Molecule dataset
            smiles_col: str - SMILES column name
            activity_cols: list/str - Column name(s) for activity/affinity (continuous values)
            class_cols: list/str - Column name(s) for class labels (binary or multi-class)
            experimental_method_cols: list/str - Column name(s) for experimental methods
            id_col: str - Column name for molecule ID
            name_col: str - Column name for molecule name
            time_col: str - Column name for time label
            use_parallel: bool - Whether to use parallel processing (default: True)
            include_experimental_info: bool - Whether to include Experimental Information Quality (5th dimension) in scoring (default: False)
            config: ScoringConfig - Scoring configuration (default: DEFAULT_CONFIG)
        """
        # Suppress RDKit warnings
        warnings.filterwarnings('ignore', module='rdkit')
        DisableLog('rdApp.warning')

        self.config = config or DEFAULT_CONFIG
        self.df = df.copy()
        self.smiles_col = smiles_col
        self.use_parallel = use_parallel
        self.include_experimental_info = include_experimental_info

        # Convert single column names to lists if necessary
        self.activity_cols = [activity_cols] if isinstance(activity_cols, str) else activity_cols or []
        self.class_cols = [class_cols] if isinstance(class_cols, str) else class_cols or []
        self.experimental_method_cols = [experimental_method_cols] if isinstance(experimental_method_cols, str) else experimental_method_cols or []

        self.id_col = id_col
        self.name_col = name_col
        self.time_col = time_col

        # Verify if SMILES column exists
        if self.smiles_col not in self.df.columns:
            raise ValueError(f"SMILES column '{self.smiles_col}' does not exist in the dataset")

        # Dataset information
        self.num_molecules = len(df)
        print(f"Dataset information: {self.num_molecules} molecules")

        # Initialize score results dictionary with conditional structure
        self.scores = {
            "Structural Integrity": {
                "Valid SMILES": 0.0,
                "Representation Consistency": 0.0,
                "Stereochemistry Completeness": 0.0,
                "Total": 0.0,
                "Normalized Total": 0.0
            },
            "Data Quality": {  # Changed from "Data Consistency"
                "Label Consistency": 0.0,
                "Data Consistency and Reliability": 0.0,  # New combined metric
                "Total": 0.0,
                "Normalized Total": 0.0
            },
            "Chemical Space Coverage": {
                "Chemical Diversity": 0.0,
                "Drug-likeness": 0.0,
                "Total": 0.0,
                "Normalized Total": 0.0
            },
            "Data Distribution": {
                "Data Size": 0.0,
                "Data Balance and Distribution": 0.0,  # Combined metric
                "Total": 0.0,
                "Normalized Total": 0.0
            },
            "Total Score": 0.0,
            "Normalized Score": 0.0,
            "Final Adjusted Score": 0.0
        }

        # Only add Experimental Information Quality if include_experimental_info is True
        if self.include_experimental_info:
            self.scores["Experimental Information Quality"] = {
                "Time Label Availability": 0.0,
                "Useful Column Quality": 0.0,
                "Classification Confidence": 0.0,
                "Type Diversity": 0.0,
                "Total": 0.0,
                "Normalized Total": 0.0
            }

        # Initialize analysis results
        self.analysis_results = {}

        # Track completed checks
        self.completed_checks = set()

    def to_result(self) -> ScoringResult:
        """Return a typed, JSON-serializable snapshot of current scores."""
        return ScoringResult.from_scores_dict(self.scores)

    def to_snapshot(self) -> ScorerSnapshot:
        """Return a snapshot for the Viz layer."""
        return ScorerSnapshot.from_scorer(self)

