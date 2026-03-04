import pandas as pd
from rdkit import Chem

from ..._logging import get_logger

logger = get_logger(__name__)


class CleaningCoreMixin:
    def clean_database(self, remove_invalid_smiles=True, 
                   remove_undefined_stereochemistry=True,
                   remove_conflicting_labels=True,
                   remove_consistent_duplicates=True,
                   verbose=True):
        """
        Clean the molecular database based on quality issues identified during scoring.
        
        Parameters:
            remove_invalid_smiles: bool - Remove molecules with invalid SMILES
            remove_undefined_stereochemistry: bool - Remove molecules with undefined chirality
            remove_conflicting_labels: bool - Remove molecules with conflicting binary labels
            remove_consistent_duplicates: bool - Remove structural duplicates with consistent activity data
            verbose: bool - Print detailed cleaning information
            
        Returns:
            cleaned_df: pandas DataFrame - Cleaned dataset
            cleaning_report: dict - Details about what was removed
        """
        # Ensure all necessary checks have been run
        if 'validate_smiles' not in self.completed_checks:
            self.validate_smiles()
        
        # Start with original dataframe
        cleaned_df = self.df.copy()
        cleaning_report = {}
        indices_to_remove = set()
        
        # 1. Remove invalid SMILES
        if remove_invalid_smiles and hasattr(self, 'invalid_indices'):
            invalid_count = len(self.invalid_indices)
            indices_to_remove.update(self.invalid_indices)
            cleaning_report['Invalid SMILES'] = {
                'removed_count': invalid_count,
                'removed_indices': list(self.invalid_indices[:10]) if len(self.invalid_indices) > 10 else list(self.invalid_indices),
                'total_invalid': invalid_count
            }
            if verbose:
                print(f"Removing {invalid_count} molecules with invalid SMILES")
        
        # 2. Remove molecules with undefined stereochemistry
        if remove_undefined_stereochemistry:
            if 'check_stereochemistry' not in self.completed_checks:
                self.check_stereochemistry()
            
            undefined_chirality_indices = []
            if hasattr(self, 'valid_df') and not self.valid_df.empty:
                # Re-check chirality for valid molecules
                for idx, row in self.valid_df.iterrows():
                    try:
                        mol = Chem.MolFromSmiles(row['canonical_smiles'])
                        if mol:
                            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                            if any(flag == '?' for _, flag in chiral_centers):
                                # Get original index from valid_df
                                original_idx = row['original_index']
                                undefined_chirality_indices.append(original_idx)
                    except Exception:
                        pass
                    
            if undefined_chirality_indices:
                indices_to_remove.update(undefined_chirality_indices)
                cleaning_report['Undefined Stereochemistry'] = {
                    'removed_count': len(undefined_chirality_indices),
                    'removed_indices': undefined_chirality_indices[:10] if len(undefined_chirality_indices) > 10 else undefined_chirality_indices
                }
                if verbose:
                    print(f"Removing {len(undefined_chirality_indices)} molecules with undefined stereochemistry")
        
        # 3. Remove molecules with conflicting labels AND duplicate molecules with identical labels
        if remove_conflicting_labels and self.class_cols:
            if 'check_label_consistency' not in self.completed_checks:
                self.check_label_consistency()

            conflicting_label_indices = []
            duplicate_label_indices = []

            if hasattr(self, 'valid_df') and not self.valid_df.empty:
                grouped = self.valid_df.groupby('canonical_smiles')

                for smiles, group in grouped:
                    if len(group) > 1:
                        # Collect all label values for each molecule in group
                        label_tuples = []
                        for idx, row in group.iterrows():
                            mol_labels = []
                            for class_col in self.class_cols:
                                if class_col in group.columns:
                                    value = row[class_col]
                                    mol_labels.append(value if pd.notna(value) else None)
                            label_tuples.append(tuple(mol_labels))

                        # Check if all labels are identical
                        unique_label_tuples = set(label_tuples)

                        if len(unique_label_tuples) > 1:
                            # Has conflicting labels - remove all
                            original_indices = group['original_index'].tolist()
                            conflicting_label_indices.extend(original_indices)
                        else:
                            # All labels are identical - keep only first
                            original_indices = group['original_index'].tolist()
                            duplicate_label_indices.extend(original_indices[1:])

            # Add conflicting label indices to removal set
            if conflicting_label_indices:
                indices_to_remove.update(conflicting_label_indices)
                if verbose:
                    print(f"Removing {len(conflicting_label_indices)} molecules with conflicting labels")

            # Add duplicate label indices to removal set
            if duplicate_label_indices:
                indices_to_remove.update(duplicate_label_indices)
                if verbose:
                    print(f"Removing {len(duplicate_label_indices)} duplicate molecules with identical labels")

            # Combined report
            if conflicting_label_indices or duplicate_label_indices:
                cleaning_report['Label Handling'] = {
                    'conflicting_labels_removed': len(conflicting_label_indices),
                    'duplicate_labels_removed': len(duplicate_label_indices),
                    'total_removed': len(conflicting_label_indices) + len(duplicate_label_indices),
                    'conflicting_indices': conflicting_label_indices[:5] if len(conflicting_label_indices) > 5 else conflicting_label_indices,
                    'duplicate_indices': duplicate_label_indices[:5] if len(duplicate_label_indices) > 5 else duplicate_label_indices
                }
        
        # 4. Remove structural duplicates with identical activity data
        if remove_consistent_duplicates and self.activity_cols:
            if 'check_data_consistency_and_reliability' not in self.completed_checks:
                self.check_data_consistency_and_reliability()
            
            duplicate_indices_to_remove = []
            if hasattr(self, 'valid_df') and not self.valid_df.empty:
                grouped = self.valid_df.groupby('canonical_smiles')
                
                for smiles, group in grouped:
                    if len(group) > 1:
                        # Check if activity data is exactly identical
                        activity_data_identical = True
                        
                        # Get activity data for all molecules in this group
                        activity_values_list = []
                        for idx, row in group.iterrows():
                            mol_activity_values = []
                            for activity_col in self.activity_cols:
                                if activity_col in group.columns:
                                    value = row[activity_col]
                                    # Use None for NaN values to distinguish from 0
                                    mol_activity_values.append(value if pd.notna(value) else None)
                            activity_values_list.append(tuple(mol_activity_values))
                        
                        # Find molecules with exactly identical activity values
                        seen_values = {}
                        for i, values in enumerate(activity_values_list):
                            original_idx = group.iloc[i]['original_index']
                            
                            if values in seen_values:
                                # This is a duplicate with identical values
                                duplicate_indices_to_remove.append(original_idx)
                            else:
                                # First occurrence of this value combination
                                seen_values[values] = original_idx
            
            if duplicate_indices_to_remove:
                indices_to_remove.update(duplicate_indices_to_remove)
                cleaning_report['Identical Duplicates'] = {
                    'removed_count': len(duplicate_indices_to_remove),
                    'removed_indices': duplicate_indices_to_remove[:10] if len(duplicate_indices_to_remove) > 10 else duplicate_indices_to_remove,
                    'note': 'Removed exact duplicates with identical activity values'
                }
                if verbose:
                    print(f"Removing {len(duplicate_indices_to_remove)} duplicate molecules with identical activity data")
        
        # Apply all removals
        if indices_to_remove:
            indices_to_keep = [i for i in range(len(cleaned_df)) if i not in indices_to_remove]
            cleaned_df = cleaned_df.iloc[indices_to_keep].reset_index(drop=True)
            
            cleaning_report['Summary'] = {
                'original_size': len(self.df),
                'removed_total': len(indices_to_remove),
                'final_size': len(cleaned_df),
                'retention_rate': f"{(len(cleaned_df) / len(self.df) * 100):.2f}%"
            }
            
            if verbose:
                print(f"\n=== Cleaning Summary ===")
                print(f"Original dataset: {len(self.df)} molecules")
                print(f"Total removed: {len(indices_to_remove)} molecules")
                print(f"Cleaned dataset: {len(cleaned_df)} molecules")
                print(f"Retention rate: {(len(cleaned_df) / len(self.df) * 100):.2f}%")
        else:
            cleaning_report['Summary'] = {
                'original_size': len(self.df),
                'removed_total': 0,
                'final_size': len(cleaned_df),
                'retention_rate': '100.00%',
                'note': 'No molecules removed'
            }
            if verbose:
                print("No molecules removed during cleaning")
        
        return cleaned_df, cleaning_report

