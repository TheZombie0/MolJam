from .._common import *


class ChemicalDiversityChecksMixin:
    CHEMICAL_DIVERSITY_SAMPLE_SIZE = 1000

    def analyze_chemical_diversity(self):
        """Analyze chemical diversity of the dataset"""
        if 'validate_smiles' not in self.completed_checks:
            self.validate_smiles()

        if len(self.valid_mols) == 0:
            print("No valid molecules, skipping chemical diversity analysis")
            score = 0
            self.scores["Chemical Space Coverage"]["Chemical Diversity"] = score
            self.analysis_results['Chemical Diversity'] = {
                'Average Tanimoto similarity': "N/A",
                'Similarity sample size': 0,
                'Scaffold count': 0,
                'Scaffold ratio': "0.0000",
                'Full scaffold ratio': "0.0000",
                'Scaffold score sample size': 0,
                'Scaffold frequency': {},
                'Top scaffolds': [],
                'Note': "No valid molecules",
            }
            self.completed_checks.add('analyze_chemical_diversity')
            return score

        print("Starting chemical diversity analysis...")
        start_time = time.time()
        canonical_smiles_list = self.valid_df['canonical_smiles'].tolist()

        diversity_score = 0
        avg_similarity = "Calculation failed"
        similarity_sample_size = 0
        sampled_indices = []

        try:
            from rdkit import DataStructs

            if self.use_parallel:
                n_workers = min(cpu_count(), 100)
                with Pool(n_workers) as pool:
                    fingerprint_results = pool.map(calculate_mol_fingerprint, canonical_smiles_list)
            else:
                fingerprint_results = [calculate_mol_fingerprint(smiles) for smiles in canonical_smiles_list]

            valid_fingerprint_entries = [
                (idx, fp) for idx, fp in enumerate(fingerprint_results) if fp is not None
            ]

            similarity_sample_size = min(self.CHEMICAL_DIVERSITY_SAMPLE_SIZE, len(valid_fingerprint_entries))
            if similarity_sample_size > 0:
                rng = np.random.default_rng(42)
                sampled_positions = rng.choice(
                    len(valid_fingerprint_entries), similarity_sample_size, replace=False
                )
                sampled_entries = [valid_fingerprint_entries[i] for i in sampled_positions]
                sampled_indices = [idx for idx, _ in sampled_entries]
                sample_fps = [fp for _, fp in sampled_entries]
            else:
                sample_fps = []

            if len(sample_fps) > 1:
                similarities = []
                for i in range(len(sample_fps)):
                    for j in range(i + 1, len(sample_fps)):
                        similarity = DataStructs.TanimotoSimilarity(sample_fps[i], sample_fps[j])
                        similarities.append(similarity)

                if similarities:
                    avg_similarity = np.mean(similarities)
                    similarity_error_rate = avg_similarity * 100
                    diversity_score = self.calculate_quality_score(
                        similarity_error_rate,
                        max_score=5,
                        threshold_low=30,
                        threshold_high=70,
                    )
                else:
                    diversity_score = 0
                    avg_similarity = "Could not calculate"
            else:
                diversity_score = 0
                avg_similarity = "Too few samples"
        except Exception as e:
            print(f"Warning: Failed to calculate molecular similarity: {str(e)}")
            diversity_score = 0
            avg_similarity = "Calculation failed"

        if self.use_parallel:
            with Pool(min(cpu_count(), 100)) as pool:
                scaffold_results = pool.map(calculate_mol_scaffold, canonical_smiles_list)
        else:
            scaffold_results = [calculate_mol_scaffold(smiles) for smiles in canonical_smiles_list]

        scaffold_frequency = {}
        for scaffold in scaffold_results:
            if scaffold is not None:
                scaffold_frequency[scaffold] = scaffold_frequency.get(scaffold, 0) + 1

        sorted_scaffolds = sorted(scaffold_frequency.items(), key=lambda x: x[1], reverse=True)
        top_scaffolds = []
        for scaffold_smiles, count in sorted_scaffolds[:20]:
            percentage = (count / len(self.valid_mols)) * 100
            top_scaffolds.append({
                'scaffold_smiles': scaffold_smiles,
                'count': count,
                'percentage': f"{percentage:.2f}%",
            })

        full_scaffolds = set(s for s in scaffold_results if s is not None)
        scaffold_count = len(full_scaffolds)
        full_scaffold_ratio = scaffold_count / len(self.valid_mols) if self.valid_mols else 0

        if not sampled_indices and scaffold_results:
            fallback_sample_size = min(self.CHEMICAL_DIVERSITY_SAMPLE_SIZE, len(scaffold_results))
            if fallback_sample_size > 0:
                rng = np.random.default_rng(42)
                sampled_indices = rng.choice(
                    len(scaffold_results), fallback_sample_size, replace=False
                ).tolist()

        sampled_scaffolds = {
            scaffold_results[idx] for idx in sampled_indices if scaffold_results[idx] is not None
        }
        sampled_scaffold_ratio = (
            len(sampled_scaffolds) / len(sampled_indices) if sampled_indices else 0
        )

        low_diversity_rate = (1 - sampled_scaffold_ratio) * 100
        scaffold_score = self.calculate_quality_score(
            low_diversity_rate,
            max_score=5,
            threshold_low=20,
            threshold_high=60,
        )

        if isinstance(diversity_score, (int, float)):
            total_score = max(0, min(10, diversity_score + scaffold_score))
        else:
            total_score = scaffold_score

        self.analysis_results['Chemical Diversity'] = {
            'Average Tanimoto similarity': f"{avg_similarity:.4f}" if isinstance(avg_similarity, float) else avg_similarity,
            'Similarity sample size': similarity_sample_size,
            'Scaffold count': scaffold_count,
            'Scaffold ratio': f"{sampled_scaffold_ratio:.4f}",
            'Full scaffold ratio': f"{full_scaffold_ratio:.4f}",
            'Scaffold score sample size': len(sampled_indices),
            'Scaffold frequency': scaffold_frequency,
            'Top scaffolds': top_scaffolds,
        }

        self.scores["Chemical Space Coverage"]["Chemical Diversity"] = total_score

        elapsed_time = time.time() - start_time
        print(f"Chemical diversity analysis completed in {elapsed_time:.2f} seconds")
        similarity_display = avg_similarity if isinstance(avg_similarity, str) else f"{avg_similarity:.4f}"
        print(
            "Chemical diversity: "
            f"Average similarity {similarity_display}, "
            f"sampled scaffold ratio {sampled_scaffold_ratio:.4f}, "
            f"full scaffold ratio {full_scaffold_ratio:.4f}, "
            f"score: {total_score:.2f}/10"
        )

        self.completed_checks.add('analyze_chemical_diversity')
        return total_score
