from .._common import *


class SmilesValidationMixin:
    def validate_smiles(self):
        """Validate SMILES and build the canonical/standardized/parent pipeline."""
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
        pipeline_outputs = {}

        for result in results:
            (
                idx,
                mol,
                canonical_smiles,
                is_invalid,
                is_non_canonical,
                standardized_smiles,
                observed_parent_smiles,
                parent_smiles,
                standardization_comment,
                parent_comment,
                removed_salts,
                removed_solvents,
                duplicate_parent_fragments,
                parent_fallback,
            ) = result
            if is_invalid:
                invalid_indices.append(idx)
            else:
                valid_mols.append(mol)
                mol_dict[idx] = mol
                pipeline_outputs[idx] = {
                    'canonical_smiles': canonical_smiles,
                    'standardized_smiles': standardized_smiles,
                    'observed_parent_smiles': observed_parent_smiles,
                    'parent_smiles': parent_smiles,
                    'standardization_comment': standardization_comment,
                    'parent_comment': parent_comment,
                    'removed_salts': removed_salts,
                    'removed_solvents': removed_solvents,
                    'duplicate_parent_fragments': duplicate_parent_fragments,
                    'parent_fallback': parent_fallback,
                }
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
            ordered_standardized_smiles = []
            ordered_observed_parent_smiles = []
            ordered_parent_smiles = []
            ordered_standardization_comments = []
            ordered_parent_comments = []
            ordered_removed_salts = []
            ordered_removed_solvents = []
            ordered_duplicate_parent_fragments = []
            ordered_parent_fallback = []
            for idx in valid_indices:
                ordered_mols.append(mol_dict[idx])
                ordered_canonical_smiles.append(pipeline_outputs[idx]['canonical_smiles'])
                ordered_standardized_smiles.append(pipeline_outputs[idx]['standardized_smiles'])
                ordered_observed_parent_smiles.append(pipeline_outputs[idx]['observed_parent_smiles'])
                ordered_parent_smiles.append(pipeline_outputs[idx]['parent_smiles'])
                ordered_standardization_comments.append(pipeline_outputs[idx]['standardization_comment'])
                ordered_parent_comments.append(pipeline_outputs[idx]['parent_comment'])
                ordered_removed_salts.append(pipeline_outputs[idx]['removed_salts'])
                ordered_removed_solvents.append(pipeline_outputs[idx]['removed_solvents'])
                ordered_duplicate_parent_fragments.append(pipeline_outputs[idx]['duplicate_parent_fragments'])
                ordered_parent_fallback.append(pipeline_outputs[idx]['parent_fallback'])

            self.valid_df['rdkit_mol'] = ordered_mols
            self.valid_df['canonical_smiles'] = ordered_canonical_smiles
            self.valid_df['standardized_smiles'] = ordered_standardized_smiles
            self.valid_df['observed_parent_smiles'] = ordered_observed_parent_smiles
            self.valid_df['parent_smiles'] = ordered_parent_smiles
            self.valid_df['standardization_comment'] = ordered_standardization_comments
            self.valid_df['parent_comment'] = ordered_parent_comments
            self.valid_df['removed_salts'] = ordered_removed_salts
            self.valid_df['removed_solvents'] = ordered_removed_solvents
            self.valid_df['duplicate_parent_fragments'] = ordered_duplicate_parent_fragments
            self.valid_df['parent_fallback'] = ordered_parent_fallback
            self.valid_df['original_index'] = valid_indices  # Store original indices
            self.valid_mols = ordered_mols  # Update to maintain order

            standardized_changed = int(
                (self.valid_df['standardized_smiles'] != self.valid_df['canonical_smiles']).sum()
            )
            parent_changed = int(
                (self.valid_df['parent_smiles'] != self.valid_df['standardized_smiles']).sum()
            )
            parent_fallback_count = int(self.valid_df['parent_fallback'].sum())
        else:
            standardized_changed = 0
            parent_changed = 0
            parent_fallback_count = 0

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
        self.analysis_results['Processing Pipeline'] = {
            'Standardized molecules changed vs canonical': standardized_changed,
            'Parent molecules changed vs standardized': parent_changed,
            'Parent fallback count': parent_fallback_count,
            'Default final result': 'parent_smiles',
        }

        self.scores["Structural Integrity"]["Valid SMILES"] = score

        elapsed_time = time.time() - start_time
        print(f"SMILES validation completed in {elapsed_time:.2f} seconds")
        print(f"SMILES validity: {self.valid_rate:.2f}% valid, {self.non_canonical_rate:.2f}% non-standardized")
        print(f"Invalid count: {invalid_count}, Score: {score:.2f}/10")

        self.completed_checks.add('validate_smiles')
        return score
