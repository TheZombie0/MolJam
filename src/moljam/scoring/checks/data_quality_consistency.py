import time

import numpy as np

from ..._logging import get_logger

logger = get_logger(__name__)


class DataConsistencyReliabilityMixin:
    def check_data_consistency_and_reliability(self):
        """
        Combined check for structural duplication, activity data consistency, 
        data variability, and outlier detection
        """
        if 'validate_smiles' not in self.completed_checks:
            self.validate_smiles()

        print("Starting data consistency and reliability check...")
        start_time = time.time()

        # Initialize sub-scores
        duplication_score = 10
        activity_consistency_score = 10
        variability_score = 10
        outlier_score = 10
        
        sub_results = {}

        # 1. Structural Duplication Check
        if len(self.valid_mols) > 0 and not self.valid_df.empty:
            grouped = self.valid_df.groupby('canonical_smiles')
            duplicate_count = 0
            duplicate_examples = []
            
            for canonical_smiles, group in grouped:
                if len(group) > 1:
                    duplicate_count += len(group) - 1
                    if len(duplicate_examples) < 20:
                        original_indices = group['original_index'].tolist()
                        original_smiles_list = group[self.smiles_col].tolist()
                        duplicate_examples.append({
                            'canonical_smiles': canonical_smiles,
                            'original_smiles_with_indices': [
                                {'smiles': smiles, 'index': idx + 1}  # +1 for 1-based indexing
                                for smiles, idx in zip(original_smiles_list, original_indices)
                            ],
                            'count': len(group)
                        })
            
            duplicate_rate = (duplicate_count / len(self.valid_df)) * 100 if len(self.valid_df) > 0 else 0
            
            # Use custom scoring function for duplicates
            duplication_score = self.score_count_based_issues(
                duplicate_count,
                len(self.valid_df),
                max_score=10,
                severity='medium'
            )
            
            sub_results['Structural Duplication'] = {
                'Duplicate molecules': duplicate_count,
                'Duplication rate': f"{duplicate_rate:.2f}%",
                'Duplicate SMILES examples': duplicate_examples,
                'Score': f"{duplication_score:.2f}/10"
            }
        else:
            sub_results['Structural Duplication'] = {
                'Note': "No valid molecules for duplication check"
            }

        # 2. Activity Data Consistency Check
        if self.activity_cols and len(self.valid_mols) > 0 and not self.valid_df.empty:
            grouped = self.valid_df.groupby('canonical_smiles')
            inconsistent_data_count = 0
            activity_inconsistencies = {}

            for col in self.activity_cols:
                if col in self.valid_df.columns:
                    inconsistent_for_column = 0
                    for smiles, group in grouped:
                        if len(group) > 1:
                            values = group[col].dropna().values
                            if len(values) > 1:
                                try:
                                    mean = np.mean(values)
                                    std = np.std(values)
                                    cv = std / abs(mean) if mean != 0 else 0
                                    if cv > 0.1:
                                        inconsistent_for_column += 1
                                        inconsistent_data_count += 1
                                except Exception:
                                    pass
                    activity_inconsistencies[col] = inconsistent_for_column

            activity_consistency_score = max(0, 10 - inconsistent_data_count)
            sub_results['Activity Data Consistency'] = {
                'Molecules with inconsistent activity data': inconsistent_data_count,
                'Activity data inconsistencies by column': activity_inconsistencies
            }
        else:
            sub_results['Activity Data Consistency'] = {
                'Note': "No activity data columns or valid molecules"
            }

        # 3. Data Variability Check (keep original implementation)
        if self.activity_cols and len(self.valid_mols) > 0 and not self.valid_df.empty:
            column_results = {}
            all_cv_values = []
            
            for activity_col in self.activity_cols:
                if activity_col in self.valid_df.columns:
                    grouped = self.valid_df.groupby('canonical_smiles')
                    cv_values = []
                    
                    for smiles, group in grouped:
                        if len(group) > 1:
                            try:
                                values = group[activity_col].dropna().values
                                if len(values) > 1:
                                    mean = np.mean(values)
                                    std = np.std(values)
                                    cv = std / abs(mean) if mean != 0 else 0
                                    cv_values.append(cv)
                                    all_cv_values.append(cv)
                            except Exception:
                                pass
                    
                    if cv_values:
                        avg_cv = np.mean(cv_values)
                        max_cv = np.max(cv_values)
                        column_results[activity_col] = {
                            'Replicate measurements': len(cv_values),
                            'Average coefficient of variation': f"{avg_cv:.4f}",
                            'Maximum coefficient of variation': f"{max_cv:.4f}"
                        }
            
            if all_cv_values:
                overall_avg_cv = np.mean(all_cv_values)
                cv_rate = overall_avg_cv * 100
                variability_score = self.calculate_quality_score(cv_rate, max_score=10, threshold_low=10, threshold_high=30)
                sub_results['Data Variability'] = {
                    'Total replicate measurements': len(all_cv_values),
                    'Overall average coefficient of variation': f"{overall_avg_cv:.4f}",
                    'Details by column': column_results
                }
            else:
                sub_results['Data Variability'] = {
                    'Note': "No replicate measurements found"
                }

        # 4. Outlier Detection (keep original implementation)
        if self.activity_cols:
            column_results = {}
            total_outliers = 0
            total_valid_data = 0
            
            for activity_col in self.activity_cols:
                if activity_col in self.df.columns:
                    activity_data = self.df[activity_col].dropna()
                    
                    if len(activity_data) > 0:
                        try:
                            q1 = activity_data.quantile(0.25)
                            q3 = activity_data.quantile(0.75)
                            iqr = q3 - q1
                            
                            lower_bound = q1 - 1.5 * iqr
                            upper_bound = q3 + 1.5 * iqr
                            
                            outliers = activity_data[(activity_data < lower_bound) | (activity_data > upper_bound)]
                            outlier_count = len(outliers)
                            outlier_rate = (outlier_count / len(activity_data)) * 100
                            
                            column_results[activity_col] = {
                                'Outlier count': outlier_count,
                                'Outlier ratio': f"{outlier_rate:.2f}%",
                                'Interquartile range': f"{q1:.4f} - {q3:.4f}",
                                'Outlier range': f"< {lower_bound:.4f} or > {upper_bound:.4f}"
                            }
                            
                            total_outliers += outlier_count
                            total_valid_data += len(activity_data)
                        except Exception:
                            pass
            
            if total_valid_data > 0:
                overall_outlier_rate = (total_outliers / total_valid_data) * 100
                outlier_score = self.calculate_quality_score(overall_outlier_rate, max_score=10, threshold_low=5, threshold_high=15)
                sub_results['Outlier Detection'] = {
                    'Total outlier count': total_outliers,
                    'Total data points': total_valid_data,
                    'Overall outlier ratio': f"{overall_outlier_rate:.2f}%",
                    'Details by column': column_results
                }

        # Calculate combined score (average of all sub-scores that were calculated)
        scores_to_average = []
        if 'Structural Duplication' in sub_results and 'Note' not in sub_results['Structural Duplication']:
            scores_to_average.append(duplication_score)
        if 'Activity Data Consistency' in sub_results and 'Note' not in sub_results['Activity Data Consistency']:
            scores_to_average.append(activity_consistency_score)
        if 'Data Variability' in sub_results and 'Note' not in sub_results['Data Variability']:
            scores_to_average.append(variability_score)
        if 'Outlier Detection' in sub_results:
            scores_to_average.append(outlier_score)
        
        if scores_to_average:
            total_score = np.mean(scores_to_average)
        else:
            total_score = 0

        self.analysis_results['Data Consistency and Reliability'] = sub_results
        self.scores["Data Quality"]["Data Consistency and Reliability"] = total_score
        
        elapsed_time = time.time() - start_time
        print(f"Data consistency and reliability check completed in {elapsed_time:.2f} seconds")
        print(f"Data consistency and reliability: combined score: {total_score:.2f}/10")
        
        self.completed_checks.add('check_data_consistency_and_reliability')
        return total_score

