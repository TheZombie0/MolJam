from .._common import *


class LowScorePenaltyMixin:
    def _get_penalty_metrics(self):
        metrics = {
            "Structural Integrity": [
                "Valid SMILES",
                "Representation Consistency",
                "Stereochemistry Completeness",
            ],
            "Data Quality": [
                "Label Consistency",
                "Data Consistency and Reliability",
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

        if getattr(self, "experimental_info", False):
            metrics["Experimental Information Quality"] = [
                "Time Label Availability",
                "Annotation Support Quality",
                "Type Diversity",
            ]

        return metrics

    def apply_low_score_penalty(self):
        """
        Apply penalty to the final score if any individual metrics have very low scores.
        This implements the "barrel effect" - the weakest component limits overall quality.
        Returns both total penalty and per-metric penalty breakdown.
        """
        individual_scores = []
        score_details = []

        for category, metric_names in self._get_penalty_metrics().items():
            category_scores = self.scores.get(category, {})
            for metric in metric_names:
                score = category_scores.get(metric)
                if score is None or not isinstance(score, (int, float, np.number)):
                    continue

                individual_scores.append(float(score))
                score_details.append({
                    'category': category,
                    'metric': metric,
                    'score': float(score)
                })

        if not individual_scores:
            return 0, {}

        very_low_threshold = 2.0
        low_threshold = 4.0
        medium_threshold = 6.0

        # Calculate per-metric penalty
        metric_penalties = {}
        for detail in score_details:
            score = detail['score']
            category = detail['category']
            metric = detail['metric']

            # Calculate individual metric penalty based on score range
            if score < very_low_threshold:
                metric_penalty = 1.5  # Very low score penalty
            elif score < low_threshold:
                metric_penalty = 0.5  # Low score penalty
            elif score < medium_threshold:
                metric_penalty = 0.1  # Medium score penalty
            else:
                metric_penalty = 0.0  # No penalty

            metric_penalties[f"{category}::{metric}"] = metric_penalty

        # Calculate counts for overall penalty adjustment
        very_low_count = sum(1 for s in individual_scores if s < very_low_threshold)
        low_count = sum(1 for s in individual_scores if very_low_threshold <= s < low_threshold)
        medium_count = sum(1 for s in individual_scores if low_threshold <= s < medium_threshold)

        # Base penalty from individual metrics
        base_penalty = sum(metric_penalties.values())

        # Additional penalty if many metrics are below threshold
        below_threshold_ratio = (very_low_count + low_count + medium_count) / len(individual_scores)
        additional_penalty = 0
        if below_threshold_ratio > 0.3:
            additional_penalty = (below_threshold_ratio - 0.3) * 1.2

        # Total penalty (capped at 30)
        total_penalty = min(base_penalty + additional_penalty, 30)

        # Distribute additional penalty proportionally to existing penalties
        if additional_penalty > 0 and base_penalty > 0:
            penalty_scale = (base_penalty + additional_penalty) / base_penalty
            # Scale all metric penalties proportionally
            for key in metric_penalties:
                metric_penalties[key] *= penalty_scale
            # Apply cap to total
            if sum(metric_penalties.values()) > 30:
                cap_scale = 30 / sum(metric_penalties.values())
                for key in metric_penalties:
                    metric_penalties[key] *= cap_scale

        sorted_scores = sorted(score_details, key=lambda x: x['score'])
        lowest_metrics = sorted_scores[:5]

        self.analysis_results['Low Score Penalty'] = {
            'Total metrics evaluated': len(individual_scores),
            'Very low scores (< 2.0)': very_low_count,
            'Low scores (2.0-4.0)': low_count,
            'Medium scores (4.0-6.0)': medium_count,
            'Below threshold ratio': f"{below_threshold_ratio:.2%}",
            'Total penalty applied': f"{total_penalty:.2f}",
            'Metric penalties': metric_penalties,  # New: per-metric breakdown
            'Lowest scoring metrics': [
                {
                    'category': m['category'],
                    'metric': m['metric'],
                    'score': f"{m['score']:.2f}/10"
                } for m in lowest_metrics
            ]
        }

        return total_penalty, metric_penalties
