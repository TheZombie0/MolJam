from .._common import *


class DataBalanceDistributionMixin:
    LOG_SCALE_ACTIVITY_MARKERS = (
        "pchembl",
        "pic50",
        "pec50",
        "pki",
        "pkd",
    )

    def _is_log_scale_activity_column(self, column_name):
        column_lower = str(column_name).lower()
        return any(marker in column_lower for marker in self.LOG_SCALE_ACTIVITY_MARKERS)

    def _compute_range_adequacy(self, column_name, activity_data, std):
        q05 = float(activity_data.quantile(0.05))
        q95 = float(activity_data.quantile(0.95))
        effective_span = max(0.0, q95 - q05)

        if effective_span <= 0:
            return 0.0, "No effective range", q05, q95, effective_span

        if self._is_log_scale_activity_column(column_name):
            range_adequacy = min(1.0, effective_span / 3.0)
            return range_adequacy, "log_scale_value_span_q95_minus_q05", q05, q95, effective_span

        if q05 > 0 and q95 > 0:
            log_span = np.log10(q95) - np.log10(q05)
            range_adequacy = min(1.0, log_span / 3.0)
            return range_adequacy, "positive_value_log10_span_q95_over_q05", q05, q95, log_span

        if std > 0:
            range_adequacy = min(1.0, effective_span / (4 * std))
        else:
            range_adequacy = 0.0
        return range_adequacy, "robust_linear_span_over_4std", q05, q95, effective_span

    def analyze_data_balance_and_distribution(self):
        """Combined analysis of continuous and categorical data distribution"""
        sub_results = {}
        sub_scores = []
        
        print("Starting data balance and distribution analysis...")
        start_time = time.time()
        
        # 1. Continuous Data Distribution Analysis
        if self.activity_cols:
            continuous_results = {}
            continuous_scores = []
            
            for activity_col in self.activity_cols:
                if activity_col not in self.df.columns:
                    continuous_results[activity_col] = {
                        'Note': f"Column {activity_col} not found in dataset"
                    }
                    continue
                    
                activity_data = self.df[activity_col].dropna()
                
                if len(activity_data) <= 1:
                    continuous_results[activity_col] = {
                        'Note': f"Insufficient data in column {activity_col} for distribution analysis"
                    }
                    continue
                    
                try:
                    mean = np.mean(activity_data)
                    median = np.median(activity_data)
                    std = np.std(activity_data)
                    min_val = np.min(activity_data)
                    max_val = np.max(activity_data)
                    
                    skewness = ((mean - median) / std) if std > 0 else 0
                    
                    hist, bin_edges = np.histogram(activity_data, bins=10)
                    bin_occupancy = sum(1 for count in hist if count > 0) / len(hist)

                    range_adequacy, range_method, q05, q95, effective_span = self._compute_range_adequacy(
                        activity_col,
                        activity_data,
                        std,
                    )
                    
                    skewness_error_rate = abs(skewness) * 50
                    skewness_score = self.calculate_quality_score(skewness_error_rate, max_score=3.33, threshold_low=25, threshold_high=75)
                    
                    occupancy_error_rate = (1 - bin_occupancy) * 100
                    occupancy_score = self.calculate_quality_score(occupancy_error_rate, max_score=3.33, threshold_low=20, threshold_high=60)
                    
                    range_error_rate = (1 - range_adequacy) * 100
                    range_score = self.calculate_quality_score(range_error_rate, max_score=3.34, threshold_low=20, threshold_high=60)
                    
                    column_score = skewness_score + occupancy_score + range_score
                    continuous_scores.append(column_score)
                    
                    continuous_results[activity_col] = {
                        'Mean': f"{mean:.4f}",
                        'Median': f"{median:.4f}",
                        'Standard deviation': f"{std:.4f}",
                        'Range': f"{min_val:.4f} to {max_val:.4f}",
                        'Robust range (q05-q95)': f"{q05:.4f} to {q95:.4f}",
                        'Skewness': f"{skewness:.4f}",
                        'Bin occupancy': f"{bin_occupancy*100:.2f}%",
                        'Range adequacy': f"{range_adequacy:.4f}",
                        'Range method': range_method,
                        'Effective range span': f"{effective_span:.4f}",
                        'Score': f"{column_score:.2f}/10"
                    }
                    
                except Exception as e:
                    print(f"Warning: Failed to analyze distribution for {activity_col}: {str(e)}")
                    continuous_results[activity_col] = {
                        'Note': f"Distribution analysis failed: {str(e)}"
                    }
            
            if continuous_scores:
                continuous_total_score = np.mean(continuous_scores)
                sub_scores.append(continuous_total_score)
            
            sub_results['Continuous Data Distribution'] = {
                'Columns analyzed': len(continuous_scores),
                'Average score': f"{continuous_total_score:.2f}/10" if continuous_scores else "N/A",
                'Details by column': continuous_results
            }
        
        # 2. Categorical Data Distribution Analysis
        if self.class_cols:
            categorical_results = {}
            categorical_scores = []
            
            for class_col in self.class_cols:
                if class_col not in self.df.columns:
                    categorical_results[class_col] = {
                        'Note': f"Column {class_col} not found in dataset"
                    }
                    continue
                    
                class_data = self.df[class_col].dropna()
                
                if len(class_data) == 0:
                    categorical_results[class_col] = {
                        'Note': f"No valid data in column {class_col}"
                    }
                    continue
                    
                try:
                    class_counts = class_data.value_counts().to_dict()
                    num_classes = len(class_counts)
                    total = sum(class_counts.values())
                    
                    if num_classes == 2:
                        percentages = [count / total * 100 for count in class_counts.values()]
                        min_percentage = min(percentages)
                        max_percentage = max(percentages)
                        
                        if 30 <= min_percentage <= 50:
                            balance_score = 10
                        elif 20 <= min_percentage < 30:
                            penalty_rate = (30 - min_percentage) / 10 * 30
                            balance_score = self.calculate_quality_score(penalty_rate, max_score=10, threshold_low=10, threshold_high=40)
                        else:
                            if min_percentage < 20:
                                imbalance_rate = (20 - min_percentage) / 20 * 100
                            else:
                                imbalance_rate = 100
                            balance_score = self.calculate_quality_score(imbalance_rate, max_score=10, threshold_low=20, threshold_high=60)
                        
                        categorical_results[class_col] = {
                            'Number of classes': num_classes,
                            'Class distribution': {str(k): f"{v/total*100:.2f}%" for k, v in class_counts.items()},
                            'Binary classification': True,
                            'Distribution range': f"{min_percentage:.2f}% - {max_percentage:.2f}%",
                            'Balance score': f"{balance_score:.2f}/10"
                        }
                    else:
                        ideal_fraction = 1 / num_classes if num_classes > 0 else 0
                        deviations = []
                        for count in class_counts.values():
                            actual_fraction = count / total
                            deviation = abs(actual_fraction - ideal_fraction)
                            deviations.append(deviation)
                        
                        avg_deviation = np.mean(deviations) if deviations else 0
                        max_possible_deviation = 1 - ideal_fraction if ideal_fraction < 1 else 0
                        if max_possible_deviation > 0:
                            deviation_error_rate = (avg_deviation / max_possible_deviation) * 100
                            balance_score = self.calculate_quality_score(deviation_error_rate, max_score=10)
                        else:
                            balance_score = 10
                        
                        categorical_results[class_col] = {
                            'Number of classes': num_classes,
                            'Class distribution': {str(k): f"{v/total*100:.2f}%" for k, v in class_counts.items()},
                            'Average deviation from ideal': f"{avg_deviation:.4f}",
                            'Balance score': f"{balance_score:.2f}/10"
                        }
                    
                    categorical_scores.append(balance_score)
                    
                except Exception as e:
                    print(f"Warning: Failed to analyze distribution for {class_col}: {str(e)}")
                    categorical_results[class_col] = {
                        'Note': f"Distribution analysis failed: {str(e)}"
                    }
            
            if categorical_scores:
                categorical_total_score = np.mean(categorical_scores)
                sub_scores.append(categorical_total_score)
            
            sub_results['Categorical Data Distribution'] = {
                'Columns analyzed': len(categorical_scores),
                'Average score': f"{categorical_total_score:.2f}/10" if categorical_scores else "N/A",
                'Details by column': categorical_results
            }
        
        # Calculate combined score
        if sub_scores:
            total_score = np.mean(sub_scores)
        else:
            total_score = 0
        
        self.analysis_results['Data Balance and Distribution'] = sub_results
        self.scores["Data Distribution"]["Data Balance and Distribution"] = total_score
        
        elapsed_time = time.time() - start_time
        print(f"Data balance and distribution analysis completed in {elapsed_time:.2f} seconds")
        print(f"Data balance and distribution: combined score: {total_score:.2f}/10")
        
        self.completed_checks.add('analyze_data_balance_and_distribution')
        return total_score
