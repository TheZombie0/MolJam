from ..._logging import get_logger

logger = get_logger(__name__)


class TimeLabelAvailabilityChecksMixin:
    def check_time_label_availability(self):
        """Check if time label is available in the dataset"""
        # Only run if experimental info should be included
        if not self.include_experimental_info:
            print("Experimental information not included in scoring, skipping time label check")
            self.completed_checks.add('check_time_label_availability')
            return 0

        if self.time_col and self.time_col in self.df.columns:
            time_count = self.df[self.time_col].notna().sum()
            time_coverage = (time_count / self.num_molecules) * 100 if self.num_molecules > 0 else 0

            unique_times = self.df[self.time_col].dropna().unique()
            num_time_points = len(unique_times)

            missing_time_rate = 100 - time_coverage
            score = self.calculate_quality_score(missing_time_rate, max_score=10)

            self.analysis_results['Time Label Availability'] = {
                'Time label column': self.time_col,
                'Time coverage rate': f"{time_coverage:.2f}%",
                'Number of unique time points': num_time_points,
                'Example time points': list(unique_times[:10])
            }
        else:
            score = 0
            self.analysis_results['Time Label Availability'] = {
                'Time label column': self.time_col if self.time_col else 'None',
                'Time coverage rate': "0.00%",
                'Note': "No time label column provided or column does not exist"
            }

        self.scores["Experimental Information Quality"]["Time Label Availability"] = score
        print(f"Time label availability: score: {score:.2f}/10")

        self.completed_checks.add('check_time_label_availability')
        return score

