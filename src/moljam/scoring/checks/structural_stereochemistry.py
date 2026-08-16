from .._common import *


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
                'Total stereogenic double bonds': 0,
                'Undefined stereogenic double bonds': 0,
                'Molecules with stereogenic double bonds': 0,
                'Molecules with undefined double bond stereochemistry': 0,
                'Undefined chirality ratio': "0.00%",
                'Molecules with undefined chirality ratio': "0.00%",
                'Undefined double bond ratio': "0.00%",
                'Molecules with undefined double bond ratio': "0.00%",
                'Undefined chirality count distribution': {},
                'Undefined chirality exact count distribution': {},
                'Undefined double bond count distribution': {},
                'Undefined double bond exact count distribution': {},
                'Example molecules with undefined chirality': [],
                'Example molecules with undefined double bond stereochemistry': [],
                'Note': "No valid molecules"
            }
            self.completed_checks.add('check_stereochemistry')
            return score

        print("Starting stereochemistry check...")
        start_time = time.time()

        # Prepare data for parallel processing
        mol_data = [
            (
                idx,
                self.valid_df.iloc[idx]['canonical_smiles'],
                self.valid_df.iloc[idx].get('observed_parent_smiles', self.valid_df.iloc[idx]['canonical_smiles']),
            )
            for idx in range(len(self.valid_mols))
        ]

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
        total_stereogenic_double_bonds = 0
        undefined_double_bonds = 0
        mols_with_double_bonds = 0
        mols_with_undefined_double_bonds = 0
        undefined_double_bond_examples = []

        undefined_chirality_count_distribution = {
            '1 undefined': 0,
            '2 undefined': 0,
            '3+ undefined': 0
        }
        undefined_chirality_exact_count_distribution = {}
        undefined_double_bond_count_distribution = {
            '1 undefined': 0,
            '2 undefined': 0,
            '3+ undefined': 0
        }
        undefined_double_bond_exact_count_distribution = {}

        for result in results:
            idx = result['idx']
            detail_num_centers = result['detail_num_chiral_centers']
            detail_num_undefined = result['detail_undefined_chiral_centers']
            detail_chirality_basis_smiles = result['detail_chirality_basis_smiles']
            detail_chirality_basis_source = result['detail_chirality_basis_source']
            detail_num_double_bonds = result['detail_num_stereogenic_double_bonds']
            detail_num_undefined_double_bonds = result['detail_undefined_double_bonds']
            detail_double_bond_basis_smiles = result['detail_double_bond_basis_smiles']
            detail_double_bond_basis_source = result['detail_double_bond_basis_source']
            record_has_undefined = result['record_has_undefined_chiral']
            record_has_undefined_double_bond = result['record_has_undefined_double_bond']

            total_chiral_centers += detail_num_centers
            undefined_chiral_centers += detail_num_undefined
            total_stereogenic_double_bonds += detail_num_double_bonds
            undefined_double_bonds += detail_num_undefined_double_bonds

            if detail_num_centers > 0:
                mols_with_chiral_centers += 1
            if detail_num_double_bonds > 0:
                mols_with_double_bonds += 1

            if record_has_undefined:
                mols_with_undefined_chiral += 1

                undefined_chirality_exact_count_distribution[detail_num_undefined] = (
                    undefined_chirality_exact_count_distribution.get(detail_num_undefined, 0) + 1
                )
                if detail_num_undefined == 1:
                    undefined_chirality_count_distribution['1 undefined'] += 1
                elif detail_num_undefined == 2:
                    undefined_chirality_count_distribution['2 undefined'] += 1
                elif detail_num_undefined >= 3:
                    undefined_chirality_count_distribution['3+ undefined'] += 1

                if len(undefined_chiral_examples) < 10:
                    undefined_chiral_examples.append({
                        'smiles': self.valid_df.iloc[idx][self.smiles_col],
                        'canonical_smiles': self.valid_df.iloc[idx]['canonical_smiles'],
                        'analysis_smiles': detail_chirality_basis_smiles,
                        'analysis_source': detail_chirality_basis_source,
                        'undefined_centers': detail_num_undefined,
                        'total_centers': detail_num_centers,
                        'original_index': self.valid_df.iloc[idx]['original_index'] + 1  # 1-based
                    })

            if record_has_undefined_double_bond:
                mols_with_undefined_double_bonds += 1
                undefined_double_bond_exact_count_distribution[detail_num_undefined_double_bonds] = (
                    undefined_double_bond_exact_count_distribution.get(detail_num_undefined_double_bonds, 0) + 1
                )
                if detail_num_undefined_double_bonds == 1:
                    undefined_double_bond_count_distribution['1 undefined'] += 1
                elif detail_num_undefined_double_bonds == 2:
                    undefined_double_bond_count_distribution['2 undefined'] += 1
                elif detail_num_undefined_double_bonds >= 3:
                    undefined_double_bond_count_distribution['3+ undefined'] += 1
                if len(undefined_double_bond_examples) < 10:
                    undefined_double_bond_examples.append({
                        'smiles': self.valid_df.iloc[idx][self.smiles_col],
                        'canonical_smiles': self.valid_df.iloc[idx]['canonical_smiles'],
                        'analysis_smiles': detail_double_bond_basis_smiles,
                        'analysis_source': detail_double_bond_basis_source,
                        'undefined_double_bonds': detail_num_undefined_double_bonds,
                        'total_stereogenic_double_bonds': detail_num_double_bonds,
                        'original_index': self.valid_df.iloc[idx]['original_index'] + 1,  # 1-based
                    })

        if mols_with_chiral_centers > 0:
            undefined_chiral_rate = (undefined_chiral_centers / total_chiral_centers) * 100
            undefined_mol_rate = (mols_with_undefined_chiral / total_valid_mols) * 100
        else:
            undefined_chiral_rate = 0
            undefined_mol_rate = 0

        if total_stereogenic_double_bonds > 0:
            undefined_double_bond_rate = (undefined_double_bonds / total_stereogenic_double_bonds) * 100
            undefined_double_bond_mol_rate = (mols_with_undefined_double_bonds / total_valid_mols) * 100
        else:
            undefined_double_bond_rate = 0
            undefined_double_bond_mol_rate = 0

        score = self.calculate_quality_score(undefined_mol_rate, max_score=10)

        self.analysis_results['Stereochemistry Completeness'] = {
            'Total chiral centers': total_chiral_centers,
            'Undefined chiral centers': undefined_chiral_centers,
            'Molecules with chiral centers': mols_with_chiral_centers,
            'Molecules with undefined chirality': mols_with_undefined_chiral,
            'Total stereogenic double bonds': total_stereogenic_double_bonds,
            'Undefined stereogenic double bonds': undefined_double_bonds,
            'Molecules with stereogenic double bonds': mols_with_double_bonds,
            'Molecules with undefined double bond stereochemistry': mols_with_undefined_double_bonds,
            'Undefined chirality ratio': f"{undefined_chiral_rate:.2f}%",
            'Molecules with undefined chirality ratio': f"{undefined_mol_rate:.2f}%",
            'Undefined double bond ratio': f"{undefined_double_bond_rate:.2f}%",
            'Molecules with undefined double bond ratio': f"{undefined_double_bond_mol_rate:.2f}%",
            'Undefined chirality count distribution': undefined_chirality_count_distribution,
            'Undefined chirality exact count distribution': dict(
                sorted(undefined_chirality_exact_count_distribution.items())
            ),
            'Undefined double bond count distribution': undefined_double_bond_count_distribution,
            'Undefined double bond exact count distribution': dict(
                sorted(undefined_double_bond_exact_count_distribution.items())
            ),
            'Example molecules with undefined chirality': undefined_chiral_examples,
            'Example molecules with undefined double bond stereochemistry': undefined_double_bond_examples,
        }

        self.scores["Structural Integrity"]["Stereochemistry Completeness"] = score

        elapsed_time = time.time() - start_time
        print(f"Stereochemistry check completed in {elapsed_time:.2f} seconds")
        print(f"Stereochemistry completeness: {undefined_mol_rate:.2f}% molecules with undefined chirality, score: {score:.2f}/10")

        self.completed_checks.add('check_stereochemistry')
        return score
