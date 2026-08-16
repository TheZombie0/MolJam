import time

import numpy as np
import pandas as pd

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
    RUNTIME_CATEGORY_METRICS = {
        "Structural Integrity": [
            "Valid SMILES",
            "Representation Consistency",
            "Stereochemistry Completeness",
        ],
        "Data Quality": [
            "Label Consistency",
            "Data Consistency and Reliability",
        ],
        "Experimental Information Quality": [
            "Time Label Availability",
            "Annotation Support Quality",
            "Type Diversity",
        ],
        "Chemical Space Coverage": [
            "Chemical Diversity",
            "Drug-likeness",
        ],
        "Data Distribution": [
            "Data Size",
            "Data Balance and Distribution",
        ],
    }

    def __init__(self, df, smiles_col='smiles', activity_cols=None, class_cols=None,
             experimental_method_cols=None, id_col=None, name_col=None, time_col=None,
             use_parallel=True, experimental_info=True,
             parent_form_backend='dimorphite_dl', parent_form_ph=7.4,
             chemaxon_executable='cxcalc', dimorphite_python=None,
             dimorphite_conda_env='dimorphite'):
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
            experimental_info: bool - Whether to include Experimental Information Quality in scoring (default: True)
            parent_form_backend: str - Backend for deriving parent form ('dimorphite_dl' or 'chemaxon')
            parent_form_ph: float - Reference pH for parent-form derivation
            chemaxon_executable: str - ChemAxon CLI executable for the optional backend
            dimorphite_python: str|None - External Python executable for the Dimorphite-DL environment
            dimorphite_conda_env: str|None - Conda environment name for the Dimorphite-DL backend
        """
        self.df = df.copy()
        self.smiles_col = smiles_col
        self.use_parallel = use_parallel
        self.experimental_info = experimental_info
        self.parent_form_backend = parent_form_backend
        self.parent_form_ph = float(parent_form_ph)
        self.chemaxon_executable = chemaxon_executable
        self.dimorphite_python = dimorphite_python
        self.dimorphite_conda_env = dimorphite_conda_env

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

        if self.experimental_info:
            self.scores["Experimental Information Quality"] = {
                "Time Label Availability": 0.0,
                "Annotation Support Quality": 0.0,
                "Type Diversity": 0.0,
                "Total": 0.0,
                "Normalized Total": 0.0
            }

        # Initialize analysis results
        self.analysis_results = {}

        # Track completed checks
        self.completed_checks = set()

        # Runtime capture state for run_all_checks().
        self.runtime_profile = self._build_empty_runtime_profile()
        self._runtime_capture_active = False
        self._runtime_capture_start = None
        self._runtime_category_allocations = {}

    def _active_runtime_categories(self):
        categories = {}
        for category, metrics in self.RUNTIME_CATEGORY_METRICS.items():
            if category == "Experimental Information Quality" and not self.experimental_info:
                continue
            categories[category] = list(metrics)
        return categories

    def _build_empty_runtime_profile(self):
        categories = self._active_runtime_categories()
        metric_order = [metric for metrics in categories.values() for metric in metrics]
        return {
            "scope": "run_all_checks",
            "category_order": list(categories.keys()),
            "metric_order": metric_order,
            "category_metrics": categories,
            "total_seconds": 0.0,
            "finalization_seconds": 0.0,
            "molecules_per_second": 0.0,
            "category_seconds": {category: 0.0 for category in categories},
            "category_percentages": {category: 0.0 for category in categories},
            "metric_seconds": {metric: 0.0 for metric in metric_order},
            "metric_percentages": {metric: 0.0 for metric in metric_order},
            "metric_records": [],
            "finalization_window": None,
        }

    def _reset_runtime_profile(self):
        self.runtime_profile = self._build_empty_runtime_profile()
        self._runtime_capture_active = True
        self._runtime_capture_start = time.perf_counter()
        self._runtime_category_allocations = {}
        return self._runtime_capture_start

    def _queue_runtime_category_allocations(self, category, metric_seconds):
        if not self._runtime_capture_active:
            return

        self._runtime_category_allocations[category] = [
            (metric, max(0.0, float(seconds)))
            for metric, seconds in metric_seconds
        ]

    def _relative_runtime_seconds(self, absolute_time):
        if self._runtime_capture_start is None:
            return 0.0
        return max(0.0, float(absolute_time) - float(self._runtime_capture_start))

    def _record_runtime_metric_window(self, category, metric, window_start, window_end):
        if not self._runtime_capture_active:
            return

        record = {
            "category": category,
            "metric": metric,
            "start_seconds": self._relative_runtime_seconds(window_start),
            "end_seconds": self._relative_runtime_seconds(window_end),
        }
        record["seconds"] = max(0.0, record["end_seconds"] - record["start_seconds"])
        self.runtime_profile["metric_records"].append(record)

    def _record_runtime_allocated_windows(self, category, window_start, window_end, allocations):
        if not allocations:
            return

        total_window = max(0.0, float(window_end) - float(window_start))
        allocation_sum = sum(seconds for _, seconds in allocations)
        if allocation_sum > 0:
            scaled = [seconds * total_window / allocation_sum for _, seconds in allocations]
        else:
            scaled = [total_window / len(allocations) for _ in allocations]

        cursor = float(window_start)
        for index, ((metric, _), seconds) in enumerate(zip(allocations, scaled)):
            segment_end = float(window_end) if index == len(allocations) - 1 else cursor + seconds
            self._record_runtime_metric_window(category, metric, cursor, segment_end)
            cursor = segment_end

    def _record_runtime_step(self, category, metric, window_start, window_end):
        if not self._runtime_capture_active:
            return

        allocations = self._runtime_category_allocations.pop(category, None)
        if allocations:
            self._record_runtime_allocated_windows(category, window_start, window_end, allocations)
            return

        if metric is not None:
            self._record_runtime_metric_window(category, metric, window_start, window_end)

    def _record_runtime_finalization_window(self, window_start, window_end):
        if not self._runtime_capture_active:
            return

        finalization_seconds = max(0.0, float(window_end) - float(window_start))
        self.runtime_profile["finalization_window"] = {
            "start_seconds": self._relative_runtime_seconds(window_start),
            "end_seconds": self._relative_runtime_seconds(window_end),
            "seconds": finalization_seconds,
        }
        self.runtime_profile["finalization_seconds"] = finalization_seconds

    def _finalize_runtime_profile(self, total_end_time):
        total_seconds = max(0.0, float(total_end_time) - float(self._runtime_capture_start))
        categories = self.runtime_profile["category_order"]
        metrics = self.runtime_profile["metric_order"]

        category_seconds = {category: 0.0 for category in categories}
        metric_seconds = {metric: 0.0 for metric in metrics}

        for record in self.runtime_profile["metric_records"]:
            category_seconds[record["category"]] += record["seconds"]
            metric_seconds[record["metric"]] += record["seconds"]

        category_total = sum(category_seconds.values())
        metric_total = sum(metric_seconds.values())

        if category_total > 0:
            category_percentages = {
                category: category_seconds[category] / category_total * 100.0
                for category in categories
            }
        else:
            category_percentages = {category: 0.0 for category in categories}

        if metric_total > 0:
            metric_percentages = {
                metric: metric_seconds[metric] / metric_total * 100.0
                for metric in metrics
            }
        else:
            metric_percentages = {metric: 0.0 for metric in metrics}

        self.runtime_profile["total_seconds"] = total_seconds
        self.runtime_profile["molecules_per_second"] = (
            self.num_molecules / total_seconds if total_seconds > 0 else 0.0
        )
        self.runtime_profile["category_seconds"] = category_seconds
        self.runtime_profile["category_percentages"] = category_percentages
        self.runtime_profile["metric_seconds"] = metric_seconds
        self.runtime_profile["metric_percentages"] = metric_percentages

        self.analysis_results["Runtime Profile"] = {
            "Timed scope": "run_all_checks",
            "Total runtime (seconds)": round(total_seconds, 4),
            "Finalization runtime (seconds)": round(self.runtime_profile["finalization_seconds"], 4),
            "Molecules per second": round(self.runtime_profile["molecules_per_second"], 4),
            "Category runtimes (seconds)": {
                category: round(seconds, 4)
                for category, seconds in category_seconds.items()
            },
            "Category percentages": {
                category: f"{percentage:.2f}%"
                for category, percentage in category_percentages.items()
            },
            "Metric runtimes (seconds)": {
                metric: round(seconds, 4)
                for metric, seconds in metric_seconds.items()
            },
            "Metric percentages": {
                metric: f"{percentage:.2f}%"
                for metric, percentage in metric_percentages.items()
            },
        }

        self._runtime_capture_active = False
        self._runtime_capture_start = None
        self._runtime_category_allocations = {}
