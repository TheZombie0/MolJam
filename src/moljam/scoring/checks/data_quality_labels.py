from .._common import *


class LabelConsistencyChecksMixin:
    def check_label_consistency(self):
        """Check label consistency, especially contradicting labels"""
        if not self.class_cols:
            print("No classification label columns provided, skipping label consistency check")
            self.scores["Data Quality"]["Label Consistency"] = None
            self.analysis_results['Label Consistency'] = {
                'Note': "No classification label columns provided"
            }
            self.completed_checks.add('check_label_consistency')
            return None

        if 'validate_smiles' not in self.completed_checks:
            self.validate_smiles()

        if len(self.valid_mols) == 0 or self.valid_df.empty:
            print("No valid molecules, skipping label consistency check")
            score = 0
            self.scores["Data Quality"]["Label Consistency"] = score
            self.analysis_results['Label Consistency'] = {
                'Molecules with contradictory labels': 0,
                'Contradictory label rate': "0.00%",
                'Note': "No valid molecules"
            }
            self.completed_checks.add('check_label_consistency')
            return score

        print("Starting label consistency check...")
        start_time = time.time()

        grouped = self.valid_df.groupby('canonical_smiles')
        column_results = {}
        contradictory_examples = []  # 新增：存储冲突示例

        for class_col in self.class_cols:
            if class_col not in self.valid_df.columns:
                column_results[class_col] = {
                    'Molecules with contradictory labels': 0,
                    'Contradictory label rate': "0.00%",
                    'Note': f"Column {class_col} does not exist in the valid dataset"
                }
                continue

            contradictory_labels = 0
            col_contradictory_examples = []  # 每列的冲突示例

            for smiles, group in grouped:
                if len(group) > 1:
                    valid_labels = [l for l in group[class_col].values if not pd.isnull(l)]
                    unique_labels = np.unique(valid_labels)
                    if len(unique_labels) > 1:
                        contradictory_labels += 1
                        # 记录冲突示例（限制数量）
                        # if len(col_contradictory_examples) < 10:
                        original_indices = group['original_index'].tolist()
                        original_smiles_list = group[self.smiles_col].tolist()
                        label_values = group[class_col].tolist()
                        col_contradictory_examples.append({
                            'canonical_smiles': smiles,
                            'original_smiles_with_indices': [
                                {
                                    'smiles': orig_smiles, 
                                    'index': idx + 1,  # 1-based indexing
                                    'label': label
                                }
                                for orig_smiles, idx, label in zip(original_smiles_list, original_indices, label_values)
                            ],
                            'conflicting_labels': list(unique_labels),
                            'count': len(group)
                        })

            if len(grouped) > 0:
                contradiction_rate = (contradictory_labels / len(grouped)) * 100
            else:
                contradiction_rate = 0

            column_results[class_col] = {
                'Molecules with contradictory labels': contradictory_labels,
                'Contradictory label rate': f"{contradiction_rate:.2f}%",
                'Contradictory examples': col_contradictory_examples  # 新增
            }

        # Calculate overall score
        if self.class_cols and len(grouped) > 0:
            unique_contradictory = set()
            overall_contradictory_examples = []  # 总体冲突示例

            for class_col in self.class_cols:
                if class_col not in self.valid_df.columns:
                    continue

                for smiles, group in grouped:
                    if len(group) > 1:
                        valid_labels = [l for l in group[class_col].values if not pd.isnull(l)]
                        unique_labels = np.unique(valid_labels)
                        if len(unique_labels) > 1:
                            if smiles not in unique_contradictory:
                                unique_contradictory.add(smiles)
                                # 记录总体冲突示例
                                if len(overall_contradictory_examples) < 20:
                                    original_indices = group['original_index'].tolist()
                                    original_smiles_list = group[self.smiles_col].tolist()

                                    # 收集所有class_cols的标签
                                    all_labels = {}
                                    for col in self.class_cols:
                                        if col in group.columns:
                                            all_labels[col] = group[col].tolist()

                                    overall_contradictory_examples.append({
                                        'canonical_smiles': smiles,
                                        'original_smiles_with_indices': [
                                            {'smiles': orig_smiles, 'index': idx + 1}
                                            for orig_smiles, idx in zip(original_smiles_list, original_indices)
                                        ],
                                        'all_labels': all_labels,
                                        'count': len(group)
                                    })

            overall_contradiction_rate = (len(unique_contradictory) / len(grouped)) * 100
        else:
            unique_contradictory = set()
            overall_contradiction_rate = 0
            overall_contradictory_examples = []

        # Use custom scoring function for contradictions
        contradictory_count = len(unique_contradictory)
        score = self.score_count_based_issues(
            contradictory_count,
            len(grouped),
            max_score=10,
            severity='high'  # High severity for label contradictions
        )

        self.analysis_results['Label Consistency'] = {
            'Total molecules with contradictory labels': contradictory_count,
            'Overall contradictory label rate': f"{overall_contradiction_rate:.2f}%",
            'Example indexing note': (
                "Indices in examples are 1-based row numbers in the original input table "
                "(header row excluded for CSV files)"
            ),
            'Example SMILES note': (
                "original_smiles_with_indices preserves the raw input SMILES before "
                "canonicalization, so entries in one contradictory group can have different "
                "SMILES strings while sharing the same canonical_smiles"
            ),
            'Details by column': column_results,
            'Overall contradictory examples': overall_contradictory_examples,  # 新增
            'Score': f"{score:.2f}/10"
        }

        self.scores["Data Quality"]["Label Consistency"] = score

        elapsed_time = time.time() - start_time
        print(f"Label consistency check completed in {elapsed_time:.2f} seconds")
        print(f"Label consistency: {contradictory_count} molecules with contradictory labels, score: {score:.2f}/10")

        self.completed_checks.add('check_label_consistency')
        return score
