import time

from ..._logging import get_logger

logger = get_logger(__name__)


class OrchestrationRunnerMixin:
    def run_missing_checks(self):
        """Run any checks that haven't been completed yet"""
        # Structural Integrity checks
        if 'validate_smiles' not in self.completed_checks:
            self.validate_smiles()
        if 'check_representation_consistency' not in self.completed_checks:
            self.check_representation_consistency()
        if 'check_stereochemistry' not in self.completed_checks:
            self.check_stereochemistry()

        # Data Quality checks
        if 'check_label_consistency' not in self.completed_checks:
            self.check_label_consistency()
        if 'check_data_consistency_and_reliability' not in self.completed_checks:
            self.check_data_consistency_and_reliability()

        # Experimental Information Quality checks - conditional
        if self.include_experimental_info:
            if 'check_annotation_quality' not in self.completed_checks:
                self.check_annotation_quality()

        # Chemical Space Coverage checks
        if 'analyze_chemical_diversity' not in self.completed_checks:
            self.analyze_chemical_diversity()
        if 'analyze_druglikeness' not in self.completed_checks:
            self.analyze_druglikeness()

        # Data Distribution checks
        if 'check_data_size' not in self.completed_checks:
            self.check_data_size()
        if 'analyze_data_balance_and_distribution' not in self.completed_checks:
            self.analyze_data_balance_and_distribution()

    def run_all_checks(self):
        """Run all checks and calculate the total score"""
        print("\n=============== Starting Molecular Database Quality Assessment ===============")
        total_start_time = time.time()

        try:
            # Structural Integrity checks
            self.validate_smiles()
            self.check_representation_consistency()
            self.check_stereochemistry()

            # Data Quality checks
            self.check_label_consistency()
            self.check_data_consistency_and_reliability()

            # Experimental Information Quality checks - conditional
            if self.include_experimental_info:
                self.check_annotation_quality()

            # Chemical Space Coverage checks
            self.analyze_chemical_diversity()
            self.analyze_druglikeness()

            # Data Distribution checks
            self.check_data_size()
            self.analyze_data_balance_and_distribution()

            # Calculate total score
            final_score = self.calculate_total_scores()

            total_elapsed = time.time() - total_start_time
            print(f"\n=============== Assessment completed in {total_elapsed:.2f} seconds ===============")

            return final_score
        except Exception as e:
            print(f"Error during scoring process: {str(e)}")
            # Try to calculate scores from completed checks
            try:
                return self.calculate_total_scores()
            except Exception:
                print("Unable to calculate total score, returning 0")
                return 0

