from ..._logging import get_logger

logger = get_logger(__name__)


class TotalScoreCalculationMixin:
    def calculate_total_scores(self):
        """Calculate category scores and final total score with normalization"""
        self.run_missing_checks()

        # Structural Integrity (3 metrics now, each 10 points)
        structure_scores = [
            self.scores["Structural Integrity"]["Valid SMILES"],
            self.scores["Structural Integrity"]["Representation Consistency"],
            self.scores["Structural Integrity"]["Stereochemistry Completeness"]
        ]

        structure_scores = [score for score in structure_scores if score is not None]
        structure_total = sum(structure_scores)
        structure_max = len(structure_scores) * 10
        self.scores["Structural Integrity"]["Total"] = structure_total

        # Data Quality (2 metrics, each 10 points)
        quality_scores = [
            self.scores["Data Quality"]["Label Consistency"],
            self.scores["Data Quality"]["Data Consistency and Reliability"]
        ]

        quality_scores = [score for score in quality_scores if score is not None]
        quality_total = sum(quality_scores)
        quality_max = len(quality_scores) * 10
        self.scores["Data Quality"]["Total"] = quality_total

        # Experimental Information Quality - conditional based on include_experimental_info
        if self.include_experimental_info:
            # 4 metrics, each 10 points = 40 points total
            annotation_scores = [
                self.scores["Experimental Information Quality"]["Time Label Availability"],
                self.scores["Experimental Information Quality"]["Useful Column Quality"],
                self.scores["Experimental Information Quality"]["Classification Confidence"],
                self.scores["Experimental Information Quality"]["Type Diversity"]
            ]

            annotation_scores = [score for score in annotation_scores if score is not None]
            annotation_total = sum(annotation_scores)
            annotation_max = 40  # 4 metrics × 10 points each
            self.scores["Experimental Information Quality"]["Total"] = annotation_total
        else:
            # If not including annotation quality, don't access the key at all
            annotation_total = 0
            annotation_max = 0

        # Chemical Space Coverage (2 metrics, each 10 points)
        coverage_scores = [
            self.scores["Chemical Space Coverage"]["Chemical Diversity"],
            self.scores["Chemical Space Coverage"]["Drug-likeness"]
        ]

        coverage_scores = [score for score in coverage_scores if score is not None]
        coverage_total = sum(coverage_scores)
        coverage_max = len(coverage_scores) * 10
        self.scores["Chemical Space Coverage"]["Total"] = coverage_total

        # Data Distribution (2 metrics, each 10 points)
        distribution_scores = [
            self.scores["Data Distribution"]["Data Size"],
            self.scores["Data Distribution"]["Data Balance and Distribution"]
        ]

        distribution_scores = [score for score in distribution_scores if score is not None]
        distribution_total = sum(distribution_scores)
        distribution_max = len(distribution_scores) * 10
        self.scores["Data Distribution"]["Total"] = distribution_total

        # Calculate normalized scores for each category
        # Adjust normalization based on whether experimental info is included
        ncfg = self.config.normalization
        if self.include_experimental_info:
            pts = ncfg.points_per_module_5
            normalized_structure = (structure_total / structure_max * pts) if structure_max > 0 else 0
            normalized_quality = (quality_total / quality_max * pts) if quality_max > 0 else 0
            normalized_annotation = (annotation_total / annotation_max * pts) if annotation_max > 0 else 0
            normalized_coverage = (coverage_total / coverage_max * pts) if coverage_max > 0 else 0
            normalized_distribution = (distribution_total / distribution_max * pts) if distribution_max > 0 else 0
        else:
            pts = ncfg.points_per_module_4
            normalized_structure = (structure_total / structure_max * pts) if structure_max > 0 else 0
            normalized_quality = (quality_total / quality_max * pts) if quality_max > 0 else 0
            normalized_annotation = 0
            normalized_coverage = (coverage_total / coverage_max * pts) if coverage_max > 0 else 0
            normalized_distribution = (distribution_total / distribution_max * pts) if distribution_max > 0 else 0

        self.scores["Structural Integrity"]["Normalized Total"] = normalized_structure
        self.scores["Data Quality"]["Normalized Total"] = normalized_quality
        if self.include_experimental_info:
            self.scores["Experimental Information Quality"]["Normalized Total"] = normalized_annotation
        self.scores["Chemical Space Coverage"]["Normalized Total"] = normalized_coverage
        self.scores["Data Distribution"]["Normalized Total"] = normalized_distribution

        # Final total score (sum of all category totals) - raw score
        total_score = structure_total + quality_total + annotation_total + coverage_total + distribution_total
        total_max = structure_max + quality_max + annotation_max + coverage_max + distribution_max
        self.scores["Total Score"] = total_score

        # Calculate normalized score (0-100 scale)
        normalized_score = normalized_structure + normalized_quality + normalized_annotation + normalized_coverage + normalized_distribution
        self.scores["Normalized Score"] = normalized_score

        # Apply low score penalty
        penalty, metric_penalties = self.apply_low_score_penalty()
        final_adjusted_score = max(0, normalized_score - penalty)
        self.scores["Final Adjusted Score"] = final_adjusted_score
        self.scores["Metric Penalties"] = metric_penalties  # Store for visualization

        print("\n=============== Final Score ===============")
        print(f"Structural Integrity: {structure_total:.2f}/{structure_max} (Normalized: {normalized_structure:.2f}/{pts:.0f})")
        print(f"Data Quality: {quality_total:.2f}/{quality_max} (Normalized: {normalized_quality:.2f}/{pts:.0f})")
        if self.include_experimental_info:
            print(f"Experimental Information Quality: {annotation_total:.2f}/{annotation_max} (Normalized: {normalized_annotation:.2f}/{pts:.0f})")
        print(f"Chemical Space Coverage: {coverage_total:.2f}/{coverage_max} (Normalized: {normalized_coverage:.2f}/{pts:.0f})")
        print(f"Data Distribution: {distribution_total:.2f}/{distribution_max} (Normalized: {normalized_distribution:.2f}/{pts:.0f})")
        print(f"Total Score: {total_score:.2f}/{total_max}")
        print(f"Normalized Score: {normalized_score:.2f}/100")
        print(f"Low Score Penalty: -{penalty:.2f}")
        print(f"Final Adjusted Score: {final_adjusted_score:.2f}/100")

        return final_adjusted_score

