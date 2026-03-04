import time
from multiprocessing import Pool, cpu_count

import numpy as np

from ..chem import calculate_mol_fingerprint, calculate_mol_scaffold
from ..._logging import get_logger

logger = get_logger(__name__)


class ChemicalDiversityChecksMixin:
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
                'Scaffold count': 0,
                'Scaffold ratio': "0.0000",
                'Scaffold frequency': {},  # 新增
                'Top scaffolds': [],  # 新增
                'Note': "No valid molecules"
            }
            self.completed_checks.add('analyze_chemical_diversity')
            return score
    
        print("Starting chemical diversity analysis...")
        start_time = time.time()
    
        try:
            from rdkit import DataStructs

            # Get canonical SMILES for fingerprint calculation
            canonical_smiles_list = self.valid_df['canonical_smiles'].tolist()

            # Calculate fingerprints in parallel or serial
            if self.use_parallel:
                n_workers = min(cpu_count(), 100)
                with Pool(n_workers) as pool:
                    fps = pool.map(calculate_mol_fingerprint, canonical_smiles_list)
            else:
                # Serial processing
                fps = [calculate_mol_fingerprint(smiles) for smiles in canonical_smiles_list]
    
            # Filter out None values
            fps = [fp for fp in fps if fp is not None]
    
            # Randomly sample molecules for similarity calculation
            sample_size = min(5000, len(fps))
            if sample_size > 1:
                sample_indices = np.random.choice(len(fps), sample_size, replace=False)
                sample_fps = [fps[i] for i in sample_indices]
    
                # Calculate similarities
                similarities = []
                for i in range(sample_size):
                    for j in range(i+1, sample_size):
                        similarity = DataStructs.TanimotoSimilarity(sample_fps[i], sample_fps[j])
                        similarities.append(similarity)
    
                if similarities:
                    avg_similarity = np.mean(similarities)
                    similarity_error_rate = avg_similarity * 100
                    diversity_score = self.calculate_quality_score(similarity_error_rate, max_score=5, threshold_low=30,    threshold_high=70)
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
    
        # Calculate Bemis-Murcko scaffold diversity in parallel or serial
        if self.use_parallel:
            with Pool(min(cpu_count(), 100)) as pool:
                scaffold_results = pool.map(calculate_mol_scaffold, canonical_smiles_list)
        else:
            # Serial processing
            scaffold_results = [calculate_mol_scaffold(smiles) for smiles in canonical_smiles_list]
    
        # 新增：统计骨架频率
        scaffold_frequency = {}
        for scaffold in scaffold_results:
            if scaffold is not None:
                if scaffold not in scaffold_frequency:
                    scaffold_frequency[scaffold] = 0
                scaffold_frequency[scaffold] += 1
        
        # 排序并获取top scaffolds
        sorted_scaffolds = sorted(scaffold_frequency.items(), key=lambda x: x[1], reverse=True)
        top_scaffolds = []
        for scaffold_smiles, count in sorted_scaffolds[:20]:  # 保存前20个
            percentage = (count / len(self.valid_mols)) * 100
            top_scaffolds.append({
                'scaffold_smiles': scaffold_smiles,
                'count': count,
                'percentage': f"{percentage:.2f}%"
            })
        
        # Create set for unique scaffolds
        scaffolds = set(s for s in scaffold_results if s is not None)
        scaffold_count = len(scaffolds)
        scaffold_ratio = scaffold_count / len(self.valid_mols) if self.valid_mols else 0
    
        # Scaffold ratio score
        low_diversity_rate = (1 - scaffold_ratio) * 100
        scaffold_score = self.calculate_quality_score(low_diversity_rate, max_score=5, threshold_low=20, threshold_high=60)
    
        # Total score
        if isinstance(diversity_score, (int, float)):
            total_score = max(0, min(10, diversity_score + scaffold_score))
        else:
            total_score = scaffold_score
    
        self.analysis_results['Chemical Diversity'] = {
            'Average Tanimoto similarity': f"{avg_similarity:.4f}" if isinstance(avg_similarity, float) else avg_similarity,
            'Scaffold count': scaffold_count,
            'Scaffold ratio': f"{scaffold_ratio:.4f}",
            'Scaffold frequency': scaffold_frequency,  # 新增：完整的骨架频率字典
            'Top scaffolds': top_scaffolds  # 新增：前20个高频骨架
        }
    
        self.scores["Chemical Space Coverage"]["Chemical Diversity"] = total_score
    
        elapsed_time = time.time() - start_time
        print(f"Chemical diversity analysis completed in {elapsed_time:.2f} seconds")
        print(f"Chemical diversity: Average similarity {avg_similarity if isinstance(avg_similarity, str) else f'{avg_similarity:.4f}   '}, scaffold ratio {scaffold_ratio:.4f}, score: {total_score:.2f}/10")
    
        self.completed_checks.add('analyze_chemical_diversity')
        return total_score

