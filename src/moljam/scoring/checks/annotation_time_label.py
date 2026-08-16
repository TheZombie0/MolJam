from .._common import *


class TimeLabelAvailabilityChecksMixin:
    def check_time_label_availability(self):
        """Check if time label is available in the dataset"""
        if not self.experimental_info:
            print("Experimental information not included in scoring, skipping time label check")
            self.completed_checks.add('check_time_label_availability')
            return 0

        time_col = self.time_col if self.time_col in self.df.columns else None
        if time_col:
            time_count = self.df[time_col].notna().sum()
            time_coverage = time_count / self.num_molecules if self.num_molecules > 0 else 0

            unique_times = self.df[time_col].dropna().unique()
            num_time_points = len(unique_times)

            richness_component = min(1.0, np.log2(1 + num_time_points) / np.log2(1 + 4)) if num_time_points > 0 else 0.0
            quality_signal = 0.8 * time_coverage + 0.2 * richness_component
            error_rate = 100 * (1 - quality_signal)
            score = self.calculate_quality_score(error_rate, max_score=10)

            self.analysis_results['Time Label Availability'] = {
                'Time label column': time_col,
                'Time coverage rate': f"{time_coverage * 100:.2f}%",
                'Number of unique time points': num_time_points,
                'Example time points': list(unique_times[:10]),
                'Coverage component': round(time_coverage, 4),
                'Richness component': round(richness_component, 4),
                'Quality signal': round(quality_signal, 4),
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
