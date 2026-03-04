import time
from multiprocessing import Pool, cpu_count

from ..chem import check_mol_chirality
from ..._logging import get_logger

logger = get_logger(__name__)


class StereochemistryChecksMixin:
    def check_stereochemistry(self):
        """Check completeness of stereochemistry"""
        if 'validate_smiles' not in self.completed_checks:
            self.validate_smiles()

        if len(self.valid_mols) == 0:
            print("No valid molecules, skipping stereochemistry check")
            score = 0
            self.scores["Structural Integrity"]["Stereochemistry Completeness"] = score
            self.analysis_results['Stereochemistry Completeness'] = {
                'Total chiral centers': 0,
                'Undefined chiral centers': 0,
                'Molecules with chiral centers': 0,
                'Molecules with undefined chirality': 0,
                'Undefined chirality ratio': "0.00%",
                'Molecules with undefined chirality ratio': "0.00%",
                'Undefined chirality count distribution': {},  # 新增
                'Example molecules with undefined chirality': [],
                'Note': "No valid molecules"
            }
            self.completed_checks.add('check_stereochemistry')
            return score

        print("Starting stereochemistry check...")
        start_time = time.time()

        # Prepare data for parallel processing
        mol_data = [(idx, self.valid_df.iloc[idx]['canonical_smiles'])
                    for idx in range(len(self.valid_mols))]

        # Process molecules in parallel or serial
        if self.use_parallel:
            n_workers = min(cpu_count(), 100)
            with Pool(n_workers) as pool:
                results = pool.map(check_mol_chirality, mol_data)
        else:
            # Serial processing
            results = [check_mol_chirality(data) for data in mol_data]

        # Aggregate results
        total_valid_mols = len(self.valid_mols)
        total_chiral_centers = 0
        undefined_chiral_centers = 0
        mols_with_chiral_centers = 0
        mols_with_undefined_chiral = 0
        undefined_chiral_examples = []

        # 新增：统计不同未定义手性数量的分子数
        undefined_count_distribution = {
            '1 undefined': 0,
            '2 undefined': 0,
            '3+ undefined': 0
        }

        for idx, num_centers, num_undefined, has_undefined in results:
            total_chiral_centers += num_centers
            undefined_chiral_centers += num_undefined

            if num_centers > 0:
                mols_with_chiral_centers += 1

            if has_undefined:
                mols_with_undefined_chiral += 1

                # 统计未定义手性数量分布
                if num_undefined == 1:
                    undefined_count_distribution['1 undefined'] += 1
                elif num_undefined == 2:
                    undefined_count_distribution['2 undefined'] += 1
                elif num_undefined >= 3:
                    undefined_count_distribution['3+ undefined'] += 1

                if len(undefined_chiral_examples) < 10:
                    undefined_chiral_examples.append({
                        'smiles': self.valid_df.iloc[idx][self.smiles_col],
                        'canonical_smiles': self.valid_df.iloc[idx]['canonical_smiles'],
                        'undefined_centers': num_undefined,
                        'total_centers': num_centers,
                        'original_index': self.valid_df.iloc[idx]['original_index'] + 1  # 1-based
                    })

        if mols_with_chiral_centers > 0:
            undefined_chiral_rate = (undefined_chiral_centers / total_chiral_centers) * 100
            undefined_mol_rate = (mols_with_undefined_chiral / total_valid_mols) * 100
        else:
            undefined_chiral_rate = 0
            undefined_mol_rate = 0

        score = self.calculate_quality_score(undefined_mol_rate, max_score=10)

        self.analysis_results['Stereochemistry Completeness'] = {
            'Total chiral centers': total_chiral_centers,
            'Undefined chiral centers': undefined_chiral_centers,
            'Molecules with chiral centers': mols_with_chiral_centers,
            'Molecules with undefined chirality': mols_with_undefined_chiral,
            'Undefined chirality ratio': f"{undefined_chiral_rate:.2f}%",
            'Molecules with undefined chirality ratio': f"{undefined_mol_rate:.2f}%",
            'Undefined chirality count distribution': undefined_count_distribution,  # 新增
            'Example molecules with undefined chirality': undefined_chiral_examples
        }

        self.scores["Structural Integrity"]["Stereochemistry Completeness"] = score

        elapsed_time = time.time() - start_time
        print(f"Stereochemistry check completed in {elapsed_time:.2f} seconds")
        print(f"Stereochemistry completeness: {undefined_mol_rate:.2f}% molecules with undefined chirality, score: {score:.2f}/10")

        self.completed_checks.add('check_stereochemistry')
        return score

