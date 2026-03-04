import pandas as pd

from ..._logging import get_logger

logger = get_logger(__name__)


class ExperimentalMethodChecksMixin:
    def check_experimental_methods(self):
        """Check completeness and consistency of experimental methods"""
        # Only run if experimental info should be included
        if not self.include_experimental_info:
            print("Experimental information not included in scoring, skipping experimental method check")
            self.completed_checks.add('check_experimental_methods')
            return 0

        if not self.experimental_method_cols:
            print("No experimental method columns provided, skipping experimental method check")
            score = 0
            self.scores["Experimental Information Quality"]["Experimental Method Completeness"] = score
            self.analysis_results['Experimental Method Completeness'] = {
                'Method coverage rate': "0.00%",
                'Number of different methods': 0,
                'Method list': [],
                'Note': "No experimental method columns provided"
            }
            self.completed_checks.add('check_experimental_methods')
            return score
            
        method_results = {}
        total_coverage = 0
        total_methods = 0
        
        for col in self.experimental_method_cols:
            if col not in self.df.columns:
                method_results[col] = {
                    'Method coverage rate': "0.00%",
                    'Number of different methods': 0,
                    'Method list': [],
                    'Note': f"Column {col} does not exist in the dataset"
                }
                continue
                
            method_count = self.df[col].notna().sum()
            method_coverage = (method_count / self.num_molecules) * 100 if self.num_molecules > 0 else 0
            
            unique_methods = self.df[col].unique()
            num_methods = sum(1 for m in unique_methods if pd.notna(m))
            
            method_results[col] = {
                'Method coverage rate': f"{method_coverage:.2f}%",
                'Number of different methods': num_methods,
                'Method list': [m for m in unique_methods if pd.notna(m)]
            }
            
            total_coverage += method_coverage
            total_methods += num_methods
        
        avg_coverage = total_coverage / len(self.experimental_method_cols) if self.experimental_method_cols else 0
        
        method_consistency_score = 5
        consistency_by_method = {}
        
        if self.activity_cols and total_methods > 1:
            for activity_col in self.activity_cols:
                if activity_col in self.df.columns:
                    for method_col in self.experimental_method_cols:
                        if method_col in self.df.columns:
                            try:
                                method_stats = self.df.groupby(method_col)[activity_col].agg(['mean', 'std']).dropna()
                                
                                if len(method_stats) > 1:
                                    overall_mean = method_stats['mean'].mean()
                                    method_variation = method_stats['mean'].std() / abs(overall_mean) if overall_mean != 0 else 0
                                    
                                    consistency_by_method[f"{activity_col} by {method_col}"] = {
                                        'Between-method variation coefficient': f"{method_variation:.4f}"
                                    }
                                    
                                    method_variation_rate = method_variation * 100
                                    method_consistency_score = min(method_consistency_score, 
                                                          self.calculate_quality_score(method_variation_rate, max_score=5))
                            except Exception:
                                pass
        
        missing_coverage_rate = 100 - avg_coverage
        coverage_score = self.calculate_quality_score(missing_coverage_rate, max_score=5)
        total_score = coverage_score + method_consistency_score
        
        self.analysis_results['Experimental Method Completeness'] = {
            'Average method coverage rate': f"{avg_coverage:.2f}%",
            'Total number of different methods': total_methods,
            'Between-method consistency': method_consistency_score,
            'Details by column': method_results,
            'Method consistency details': consistency_by_method
        }
        
        self.scores["Experimental Information Quality"]["Experimental Method Completeness"] = total_score
        print(f"Experimental method analysis: {avg_coverage:.2f}% coverage, {total_methods} methods, score: {total_score:.2f}/10")
        
        self.completed_checks.add('check_experimental_methods')
        return total_score

