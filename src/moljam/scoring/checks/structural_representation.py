import time
from collections import defaultdict
from multiprocessing import Pool, cpu_count

from ..chem import process_mol_for_consistency
from ..._logging import get_logger

logger = get_logger(__name__)


class RepresentationConsistencyMixin:
    def check_representation_consistency(self):
        """
        Check for the same molecule existing in different acid/base/salt forms
        """
        if 'validate_smiles' not in self.completed_checks:
            self.validate_smiles()

        if len(self.valid_mols) == 0:
            print("No valid molecules, skipping representation consistency check")
            score = 0
            self.scores["Structural Integrity"]["Representation Consistency"] = score
            self.analysis_results['Representation Consistency'] = {
                'Molecules with multiple forms': 0,
                'Total redundant molecules': 0,
                'Redundancy rate': "0.00%",
                'Example groups': [],
                'Note': "No valid molecules"
            }
            self.completed_checks.add('check_representation_consistency')
            return score

        print("Starting representation consistency check...")
        start_time = time.time()

        # Prepare data for parallel processing - 使用原始SMILES而不是canonical_smiles
        mol_data = [(idx, self.valid_df.iloc[idx][self.smiles_col])
                    for idx in range(len(self.valid_df))]

        # Process molecules in parallel or serial
        if self.use_parallel:
            n_workers = min(cpu_count(), 100)
            with Pool(n_workers) as pool:
                results = pool.map(process_mol_for_consistency, mol_data)
        else:
            # Serial processing
            results = [process_mol_for_consistency(data) for data in mol_data]

        # Group molecules by their neutral form
        neutral_groups = defaultdict(list)
        for idx, original_smiles, neutral_smiles, formal_charge in results:
            neutral_groups[neutral_smiles].append({
                'index': idx,
                'original_smiles': original_smiles,
                'formal_charge': formal_charge
            })

        # Find groups with multiple representations
        redundant_groups = []
        total_redundant_molecules = 0

        for neutral_smiles, molecules in neutral_groups.items():
            if len(molecules) <= 1:
                continue
            
            original_smiles_set = set(mol_info['original_smiles'] for mol_info in molecules)
            charge_states = set(mol_info['formal_charge'] for mol_info in molecules)

            if len(original_smiles_set) > 1 or len(charge_states) > 1:
                redundant_groups.append({
                    'neutral_form': neutral_smiles,
                    'molecule_count': len(molecules),
                    'different_smiles': list(original_smiles_set),
                    'charge_states': list(charge_states),
                    'molecules': molecules
                })
                total_redundant_molecules += len(molecules) - 1

        # Use score_count_based_issues instead of score_low_count_issues
        redundant_count = len(redundant_groups)
        score = self.score_count_based_issues(
            redundant_count,
            len(self.valid_mols),
            max_score=10,
            severity='medium'  # Medium severity for representation consistency
        )

        # Calculate redundancy rate
        redundancy_rate = (total_redundant_molecules / len(self.valid_mols)) * 100 if self.valid_mols else 0

        # Prepare example groups
        example_groups = []
        for group in redundant_groups[:5]:
            example_groups.append({
                'neutral_form': group['neutral_form'],
                'variants': group['different_smiles'],
                'count': group['molecule_count']
            })

        self.analysis_results['Representation Consistency'] = {
            'Molecules with multiple forms': len(redundant_groups),
            'Total redundant molecules': total_redundant_molecules,
            'Redundancy rate': f"{redundancy_rate:.2f}%",
            'Example groups': example_groups
        }

        self.scores["Structural Integrity"]["Representation Consistency"] = score

        elapsed_time = time.time() - start_time
        print(f"Representation consistency check completed in {elapsed_time:.2f} seconds")
        print(f"Representation consistency: {len(redundant_groups)} molecules with multiple forms")
        print(f"Redundancy rate: {redundancy_rate:.2f}%, score: {score:.2f}/10")

        self.completed_checks.add('check_representation_consistency')
        return score

