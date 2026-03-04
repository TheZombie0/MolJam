from ..._logging import get_logger

logger = get_logger(__name__)


class LowScorePenaltyMixin:
    def apply_low_score_penalty(self):
        """
        Apply penalty to the final score if any individual metrics have very low scores.
        This implements the "barrel effect" - the weakest component limits overall quality.
        Returns both total penalty and per-metric penalty breakdown.
        """
        individual_scores = []
        score_details = []

        for category, metrics in self.scores.items():
            if category in ["Total Score", "Normalized Score", "Final Adjusted Score"]:
                continue

            if isinstance(metrics, dict):
                for metric, score in metrics.items():
                    if metric not in ["Total", "Normalized Total"] and score is not None:
                        individual_scores.append(score)
                        score_details.append({
                            'category': category,
                            'metric': metric,
                            'score': score
                        })

        if not individual_scores:
            return 0, {}

        pcfg = self.config.penalty
        very_low_threshold = pcfg.very_low_threshold
        low_threshold = pcfg.low_threshold
        medium_threshold = pcfg.medium_threshold

        # Calculate per-metric penalty
        metric_penalties = {}
        for detail in score_details:
            score = detail['score']
            category = detail['category']
            metric = detail['metric']

            # Calculate individual metric penalty based on score range
            if score < very_low_threshold:
                metric_penalty = pcfg.very_low_penalty
            elif score < low_threshold:
                metric_penalty = pcfg.low_penalty
            elif score < medium_threshold:
                metric_penalty = pcfg.medium_penalty
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
        if below_threshold_ratio > pcfg.ratio_cutoff:
            additional_penalty = (below_threshold_ratio - pcfg.ratio_cutoff) * pcfg.additional_multiplier

        # Total penalty (capped)
        total_penalty = min(base_penalty + additional_penalty, pcfg.total_cap)

        # Distribute additional penalty proportionally to existing penalties
        if additional_penalty > 0 and base_penalty > 0:
            penalty_scale = (base_penalty + additional_penalty) / base_penalty
            # Scale all metric penalties proportionally
            for key in metric_penalties:
                metric_penalties[key] *= penalty_scale
            # Apply cap to total
            if sum(metric_penalties.values()) > pcfg.total_cap:
                cap_scale = pcfg.total_cap / sum(metric_penalties.values())
                for key in metric_penalties:
                    metric_penalties[key] *= cap_scale

        sorted_scores = sorted(score_details, key=lambda x: x['score'])
        lowest_metrics = sorted_scores[:5]

        self.analysis_results['Low Score Penalty'] = {
            'Total metrics evaluated': len(individual_scores),
            f'Very low scores (< {very_low_threshold})': very_low_count,
            f'Low scores ({very_low_threshold}-{low_threshold})': low_count,
            f'Medium scores ({low_threshold}-{medium_threshold})': medium_count,
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

