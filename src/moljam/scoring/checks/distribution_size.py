import numpy as np

from ..._logging import get_logger

logger = get_logger(__name__)


class DataSizeChecksMixin:
    def check_data_size(self):
        """Check if the data size is adequate for machine learning using log scale"""
        # Use log10 to evaluate data size
        if self.num_molecules > 0:
            log_size = np.log10(self.num_molecules)
            
            # Define thresholds on log scale
            # log10(100) = 2, log10(1000) = 3, log10(10000) = 4, log10(100000) = 5
            if log_size >= 5:  # 100,000+
                score = 10
                size_category = 'Excellent'
            elif log_size >= 4:  # 10,000-99,999
                # Linear interpolation between 8 and 10
                score = 8 + (log_size - 4) * 2
                size_category = 'Very Good'
            elif log_size >= 3.5:  # ~3,162-9,999
                # Linear interpolation between 6 and 8
                score = 6 + (log_size - 3.5) * 4
                size_category = 'Good'
            elif log_size >= 3:  # 1,000-3,161
                # Linear interpolation between 4 and 6
                score = 4 + (log_size - 3) * 4
                size_category = 'Adequate'
            elif log_size >= 2.5:  # ~316-999
                # Linear interpolation between 2 and 4
                score = 2 + (log_size - 2.5) * 4
                size_category = 'Minimal'
            elif log_size >= 2:  # 100-315
                # Linear interpolation between 0 and 2
                score = (log_size - 2) * 4
                size_category = 'Poor'
            else:  # <100
                score = 0
                size_category = 'Insufficient'
        else:
            log_size = 0
            score = 0
            size_category = 'No data'
        
        self.analysis_results['Data Size'] = {
            'Total molecules': self.num_molecules,
            'Log10(size)': f"{log_size:.2f}",
            'Size category': size_category,
            'Score calculation': 'Based on log10 scale',
            'Score breakdown': {
                'Excellent (10 points)': '>= 100,000 molecules (log10 >= 5)',
                'Very Good (8-10 points)': '10,000-99,999 molecules (log10 = 4-5)',
                'Good (6-8 points)': '~3,162-9,999 molecules (log10 = 3.5-4)',
                'Adequate (4-6 points)': '1,000-3,161 molecules (log10 = 3-3.5)',
                'Minimal (2-4 points)': '~316-999 molecules (log10 = 2.5-3)',
                'Poor (0-2 points)': '100-315 molecules (log10 = 2-2.5)',
                'Insufficient (0 points)': '< 100 molecules (log10 < 2)'
            }
        }
        
        self.scores["Data Distribution"]["Data Size"] = score
        print(f"Data size: {self.num_molecules} molecules (log10 = {log_size:.2f}, {size_category}), score: {score:.2f}/10")
        
        self.completed_checks.add('check_data_size')
        return score

