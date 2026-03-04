"""Typed data models for scoring results. JSON-serializable."""

import json
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class CategoryScore:
    """Score for a single scoring category (e.g., Structural Integrity)."""

    metrics: Dict[str, Optional[float]] = field(default_factory=dict)
    total: float = 0.0
    normalized_total: float = 0.0

    def to_dict(self) -> dict:
        return {
            **{k: v for k, v in self.metrics.items()},
            "Total": self.total,
            "Normalized Total": self.normalized_total,
        }

    @classmethod
    def from_scores_dict(cls, d: dict) -> "CategoryScore":
        metrics = {
            k: v
            for k, v in d.items()
            if k not in ("Total", "Normalized Total")
        }
        return cls(
            metrics=metrics,
            total=d.get("Total", 0.0),
            normalized_total=d.get("Normalized Total", 0.0),
        )


@dataclass
class ScoringResult:
    """Typed representation of all scoring results."""

    structural_integrity: CategoryScore = field(default_factory=CategoryScore)
    data_quality: CategoryScore = field(default_factory=CategoryScore)
    chemical_space_coverage: CategoryScore = field(default_factory=CategoryScore)
    data_distribution: CategoryScore = field(default_factory=CategoryScore)
    experimental_info_quality: Optional[CategoryScore] = None
    total_score: float = 0.0
    normalized_score: float = 0.0
    final_adjusted_score: float = 0.0
    metric_penalties: Dict[str, float] = field(default_factory=dict)

    @classmethod
    def from_scores_dict(cls, scores: dict) -> "ScoringResult":
        exp = None
        if "Experimental Information Quality" in scores:
            exp_data = scores["Experimental Information Quality"]
            if isinstance(exp_data, dict):
                exp = CategoryScore.from_scores_dict(exp_data)

        penalties = scores.get("Metric Penalties", {})
        if not isinstance(penalties, dict):
            penalties = {}

        return cls(
            structural_integrity=CategoryScore.from_scores_dict(
                scores.get("Structural Integrity", {})
            ),
            data_quality=CategoryScore.from_scores_dict(
                scores.get("Data Quality", {})
            ),
            chemical_space_coverage=CategoryScore.from_scores_dict(
                scores.get("Chemical Space Coverage", {})
            ),
            data_distribution=CategoryScore.from_scores_dict(
                scores.get("Data Distribution", {})
            ),
            experimental_info_quality=exp,
            total_score=scores.get("Total Score", 0.0),
            normalized_score=scores.get("Normalized Score", 0.0),
            final_adjusted_score=scores.get("Final Adjusted Score", 0.0),
            metric_penalties=penalties,
        )

    def to_dict(self) -> dict:
        d = {
            "Structural Integrity": self.structural_integrity.to_dict(),
            "Data Quality": self.data_quality.to_dict(),
            "Chemical Space Coverage": self.chemical_space_coverage.to_dict(),
            "Data Distribution": self.data_distribution.to_dict(),
            "Total Score": self.total_score,
            "Normalized Score": self.normalized_score,
            "Final Adjusted Score": self.final_adjusted_score,
            "Metric Penalties": self.metric_penalties,
        }
        if self.experimental_info_quality is not None:
            d["Experimental Information Quality"] = self.experimental_info_quality.to_dict()
        return d

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False, default=str, indent=2)


@dataclass
class DatabaseScoringReport:
    """MCP-friendly complete scoring report. Fully JSON-serializable."""

    final_score: Optional[float] = None
    report_text: Optional[str] = None
    scores: Optional[ScoringResult] = None
    analysis_results: Optional[Dict[str, Any]] = None
    dataset_info: Optional[Dict[str, Any]] = None

    def to_dict(self) -> dict:
        d: Dict[str, Any] = {
            "final_score": self.final_score,
            "report_text": self.report_text,
        }
        if self.scores is not None:
            d["scores"] = self.scores.to_dict()
        else:
            d["scores"] = None
        d["analysis_results"] = self._sanitize(self.analysis_results)
        d["dataset_info"] = self.dataset_info
        return d

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False, default=str, indent=2)

    @staticmethod
    def _sanitize(obj):
        """Make analysis_results JSON-safe by converting non-serializable types."""
        if obj is None:
            return None
        if isinstance(obj, dict):
            return {k: DatabaseScoringReport._sanitize(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [DatabaseScoringReport._sanitize(item) for item in obj]
        if isinstance(obj, bool):
            return obj
        if isinstance(obj, (int, float)):
            return obj
        if isinstance(obj, str):
            return obj
        if isinstance(obj, set):
            return sorted(list(obj)) if obj and isinstance(next(iter(obj)), (int, float, str)) else [DatabaseScoringReport._sanitize(item) for item in obj]
        # Handle numpy types that may not be Python subclasses
        try:
            import numpy as np
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.bool_):
                return bool(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
        except ImportError:
            pass
        return str(obj)
