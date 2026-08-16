from .._common import *


class DrugLikenessChecksMixin:
    def analyze_druglikeness(self):
        """Analyze drug-likeness using RDKit's QED"""
        if 'validate_smiles' not in self.completed_checks:
            self.validate_smiles()

        if len(self.valid_mols) == 0:
            print("No valid molecules, skipping drug-likeness analysis")
            score = 0
            self.scores["Chemical Space Coverage"]["Drug-likeness"] = score
            self.analysis_results['Drug-likeness'] = {
                'Average QED score': "N/A",
                'QED distribution': "Not calculated",
                'Note': "No valid molecules"
            }
            self.completed_checks.add('analyze_druglikeness')
            return score

        print("Starting drug-likeness analysis...")
        start_time = time.time()

        # Calculate QED in parallel or serial
        canonical_smiles_list = self.valid_df['canonical_smiles'].tolist()

        if self.use_parallel:
            n_workers = min(cpu_count(), 100)
            with Pool(n_workers) as pool:
                qed_results = pool.map(calculate_qed_value, canonical_smiles_list)
        else:
            # Serial processing
            qed_results = [calculate_qed_value(smiles) for smiles in canonical_smiles_list]

        # Filter out None values
        qed_values = [qed for qed in qed_results if qed is not None]

        if not qed_values:
            print("Failed to calculate QED for any molecules")
            score = 0
            self.scores["Chemical Space Coverage"]["Drug-likeness"] = score
            self.analysis_results['Drug-likeness'] = {
                'Average QED score': "N/A",
                'QED distribution': "Calculation failed",
                'Note': "Failed to calculate QED for any molecules"
            }
            self.completed_checks.add('analyze_druglikeness')
            return score

        # Calculate QED statistics
        avg_qed = np.mean(qed_values)
        median_qed = np.median(qed_values)
        min_qed = np.min(qed_values)
        max_qed = np.max(qed_values)

        # Calculate distribution bins
        bins = {
            'Very high (0.9-1.0)': sum(1 for qed in qed_values if 0.9 <= qed <= 1.0) / len(qed_values),
            'High (0.7-0.9)': sum(1 for qed in qed_values if 0.7 <= qed < 0.9) / len(qed_values),
            'Medium (0.5-0.7)': sum(1 for qed in qed_values if 0.5 <= qed < 0.7) / len(qed_values),
            'Low (0.3-0.5)': sum(1 for qed in qed_values if 0.3 <= qed < 0.5) / len(qed_values),
            'Very low (0.0-0.3)': sum(1 for qed in qed_values if 0.0 <= qed < 0.3) / len(qed_values)
        }

        # Calculate score
        low_druglikeness_rate = (1 - avg_qed) * 100
        score = self.calculate_quality_score(low_druglikeness_rate, max_score=10, threshold_low=20, threshold_high=50)

        self.analysis_results['Drug-likeness'] = {
            'Average QED score': f"{avg_qed:.4f}",
            'Median QED score': f"{median_qed:.4f}",
            'QED range': f"{min_qed:.4f} to {max_qed:.4f}",
            'QED distribution': {k: f"{v*100:.2f}%" for k, v in bins.items()}
        }

        self.scores["Chemical Space Coverage"]["Drug-likeness"] = score

        elapsed_time = time.time() - start_time
        print(f"Drug-likeness analysis completed in {elapsed_time:.2f} seconds")
        print(f"Drug-likeness analysis: Average QED = {avg_qed:.4f}, score: {score:.2f}/10")

        self.completed_checks.add('analyze_druglikeness')
        return score

