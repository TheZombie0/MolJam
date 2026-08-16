from .._common import *


class OrchestrationRunnerMixin:
    _TIMED_METRIC_SEQUENCE = [
        ("Structural Integrity", "Valid SMILES", "validate_smiles"),
        ("Structural Integrity", "Representation Consistency", "check_representation_consistency"),
        ("Structural Integrity", "Stereochemistry Completeness", "check_stereochemistry"),
        ("Data Quality", "Label Consistency", "check_label_consistency"),
        ("Data Quality", "Data Consistency and Reliability", "check_data_consistency_and_reliability"),
        ("Experimental Information Quality", None, "check_annotation_quality"),
        ("Chemical Space Coverage", "Chemical Diversity", "analyze_chemical_diversity"),
        ("Chemical Space Coverage", "Drug-likeness", "analyze_druglikeness"),
        ("Data Distribution", "Data Size", "check_data_size"),
        ("Data Distribution", "Data Balance and Distribution", "analyze_data_balance_and_distribution"),
    ]

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
        if self.experimental_info:
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
        self._reset_runtime_profile()

        try:
            for category, metric, method_name in self._TIMED_METRIC_SEQUENCE:
                if category == "Experimental Information Quality" and not self.experimental_info:
                    continue

                step_start = time.perf_counter()
                getattr(self, method_name)()
                step_end = time.perf_counter()
                self._record_runtime_step(category, metric, step_start, step_end)

            finalization_start = time.perf_counter()
            final_score = self.calculate_total_scores()
            finalization_end = time.perf_counter()
            self._record_runtime_finalization_window(finalization_start, finalization_end)
            self._finalize_runtime_profile(finalization_end)

            total_elapsed = time.time() - total_start_time
            print(f"\n=============== Assessment completed in {total_elapsed:.2f} seconds ===============")

            return final_score
        except Exception as e:
            print(f"Error during scoring process: {str(e)}")
            if self._runtime_capture_active:
                self._finalize_runtime_profile(time.perf_counter())
            # Try to calculate scores from completed checks
            try:
                return self.calculate_total_scores()
            except Exception:
                print("Unable to calculate total score, returning 0")
                return 0
