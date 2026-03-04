import time
from multiprocessing import Pool, cpu_count

import pandas as pd
from rdkit import Chem

from ..chem import process_single_smiles
from ..._logging import get_logger

logger = get_logger(__name__)


class SmilesValidationMixin:
    def validate_smiles(self):
        """Validate SMILES validity and record canonical level (without scoring)"""
        print("Starting SMILES validation...")
        start_time = time.time()

        # Prepare data for parallel processing
        smiles_data = [(idx, smiles) for idx, smiles in enumerate(self.df[self.smiles_col])]

        if self.use_parallel:
            # Use multiprocessing to validate SMILES in parallel
            n_workers = min(cpu_count(), 100)
            with Pool(n_workers) as pool:
                results = pool.map(process_single_smiles, smiles_data)
        else:
            # Serial processing
            results = [process_single_smiles(data) for data in smiles_data]

        # Process results
        valid_mols = []
        invalid_indices = []
        non_canonical_indices = []
        mol_dict = {}  # Store mol objects by index
        canonical_smiles_list = []

        for result in results:
            idx, mol, canonical_smiles, is_invalid, is_non_canonical = result
            if is_invalid:
                invalid_indices.append(idx)
            else:
                valid_mols.append(mol)
                mol_dict[idx] = mol
                canonical_smiles_list.append(canonical_smiles)
                if is_non_canonical:
                    non_canonical_indices.append(idx)

        self.valid_mols = valid_mols
        self.invalid_indices = invalid_indices
        self.non_canonical_indices = non_canonical_indices

        # Calculate various rates
        self.valid_rate = (len(valid_mols) / self.num_molecules) * 100 if self.num_molecules > 0 else 0
        self.invalid_rate = (len(invalid_indices) / self.num_molecules) * 100 if self.num_molecules > 0 else 0
        self.non_canonical_rate = (len(non_canonical_indices) / self.num_molecules) * 100 if self.num_molecules > 0 else 0

        # Create dataframe for valid molecules with original indices
        valid_indices = [i for i in range(self.num_molecules) if i not in invalid_indices]
        self.valid_df = self.df.iloc[valid_indices].copy() if valid_indices else pd.DataFrame()

        if not self.valid_df.empty:
            # Recreate mols in correct order
            ordered_mols = []
            ordered_canonical_smiles = []
            for idx in valid_indices:
                ordered_mols.append(mol_dict[idx])
                ordered_canonical_smiles.append(Chem.MolToSmiles(mol_dict[idx]))

            self.valid_df['rdkit_mol'] = ordered_mols
            self.valid_df['canonical_smiles'] = ordered_canonical_smiles
            self.valid_df['original_index'] = valid_indices  # Store original indices
            self.valid_mols = ordered_mols  # Update to maintain order

        # Use score_count_based_issues instead of score_low_count_issues
        invalid_count = len(invalid_indices)
        score = self.score_count_based_issues(
            invalid_count, 
            self.num_molecules, 
            max_score=10,
            severity='high'  # High severity for invalid SMILES
        )

        self.analysis_results['Valid SMILES'] = {
            'Valid molecules': len(valid_mols),
            'Invalid molecules': len(invalid_indices),
            'Non-standardized molecules': len(non_canonical_indices),
            'Valid rate': f"{self.valid_rate:.2f}%",
            'Invalid rate': f"{self.invalid_rate:.2f}%",
            'Non-standardized rate': f"{self.non_canonical_rate:.2f}%",
            'Validity score': f"{score:.2f}/10"
        }

        self.scores["Structural Integrity"]["Valid SMILES"] = score

        elapsed_time = time.time() - start_time
        print(f"SMILES validation completed in {elapsed_time:.2f} seconds")
        print(f"SMILES validity: {self.valid_rate:.2f}% valid, {self.non_canonical_rate:.2f}% non-standardized")
        print(f"Invalid count: {invalid_count}, Score: {score:.2f}/10")

        self.completed_checks.add('validate_smiles')
        return score

