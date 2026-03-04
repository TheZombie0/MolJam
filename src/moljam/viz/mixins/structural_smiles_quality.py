from .._common import *
from .. import plotting


class SmilesQualityPlotMixin:
    def plot_invalid_smiles_comparison(self, mode='both', figsize=(14, 8), save_path=None):
        """
        Plot invalid SMILES count and ratio comparison (updated both mode)
        """
        if not self.scoring_results:
            print("No scoring results available.")
            return
            
        db_names = []
        invalid_counts = []
        invalid_ratios = []
        
        for db_name, results in self.scoring_results.items():
            snap = results['snapshot']
            
            if hasattr(snap, 'invalid_indices'):
                invalid_count = len(snap.invalid_indices)
                invalid_ratio = snap.invalid_rate
            else:
                invalid_count = 0
                invalid_ratio = 0
                
            db_names.append(db_name)
            invalid_counts.append(invalid_count)
            invalid_ratios.append(invalid_ratio)
        
        if mode == 'both':
            # Create two subplots side by side
            fig, (ax1, ax2) = plotting.plt.subplots(1, 2, figsize=figsize)
            
            # Left plot: Count
            bars1 = ax1.bar(db_names, invalid_counts,
                           color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax1.set_xlabel('Database', fontsize=21)
            ax1.set_ylabel('Invalid SMILES Count', fontsize=21)
            ax1.set_title('Invalid SMILES Count', fontsize=22)
            ax1.tick_params(axis='x', rotation=45, labelsize=18)
            plotting.plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
            
            # Add value labels
            for bar, count in zip(bars1, invalid_counts):
                if int(count)>0:
                    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(invalid_counts)*0.01,
                            str(count), ha='center', va='bottom', fontsize=9)
            
            # Right plot: Ratio
            bars2 = ax2.bar(db_names, invalid_ratios,
                           color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax2.set_xlabel('Database', fontsize=21)
            ax2.set_ylabel('Invalid SMILES Ratio (%)', fontsize=21)
            ax2.set_title('Invalid SMILES Ratio', fontsize=22)
            plotting.plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
            ax2.tick_params(axis='x', rotation=45, labelsize=18)
            
            # Add value labels
            for bar, ratio in zip(bars2, invalid_ratios):
                if float(ratio)>0:
                    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(invalid_ratios)*0.01,
                            f'{ratio:.1f}%', ha='center', va='bottom', fontsize=9)
            
            plotting.plt.suptitle('Invalid SMILES Comparison Across Databases', fontsize=23)
            
        elif mode == 'count':
            fig, ax = plotting.plt.subplots(figsize=figsize)
            bars = ax.bar(db_names, invalid_counts,
                          color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax.set_xlabel('Database', fontsize=21)
            ax.set_ylabel('Invalid SMILES Count', fontsize=21)
            ax.set_title('Invalid SMILES Count Comparison', fontsize=22)
            
            for bar, count in zip(bars, invalid_counts):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(invalid_counts)*0.01,
                       str(count), ha='center', va='bottom', fontsize=9)
            
            plotting.plt.xticks(rotation=45, ha='right')
            ax.tick_params(labelsize=18)
            
        else:  # mode == 'ratio'
            fig, ax = plotting.plt.subplots(figsize=figsize)
            bars = ax.bar(db_names, invalid_ratios,
                          color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax.set_xlabel('Database', fontsize=21)
            ax.set_ylabel('Invalid SMILES Ratio (%)', fontsize=21)
            ax.set_title('Invalid SMILES Ratio Comparison', fontsize=22)
            
            for bar, ratio in zip(bars, invalid_ratios):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(invalid_ratios)*0.01,
                       f'{ratio:.1f}%', ha='center', va='bottom', fontsize=9)
            
            plotting.plt.xticks(rotation=45, ha='right')
            ax.tick_params(labelsize=18)
        
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.comparison_dir, f'invalid_smiles_{mode}.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

    def plot_non_standardized_smiles_comparison(self, mode='both', figsize=(14, 8), save_path=None):
        """
        Plot non-standardized SMILES count and ratio comparison (updated both mode)
        """
        if not self.scoring_results:
            print("No scoring results available.")
            return
            
        db_names = []
        non_canonical_counts = []
        non_canonical_ratios = []
        
        for db_name, results in self.scoring_results.items():
            snap = results['snapshot']
            
            if hasattr(snap, 'non_canonical_indices'):
                non_canonical_count = len(snap.non_canonical_indices)
                non_canonical_ratio = snap.non_canonical_rate
            else:
                non_canonical_count = 0
                non_canonical_ratio = 0
                
            db_names.append(db_name)
            non_canonical_counts.append(non_canonical_count)
            non_canonical_ratios.append(non_canonical_ratio)
        
        if mode == 'both':
            # Create two subplots side by side
            fig, (ax1, ax2) = plotting.plt.subplots(1, 2, figsize=figsize)
            
            # Left plot: Count
            bars1 = ax1.bar(db_names, non_canonical_counts,
                           color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax1.set_xlabel('Database', fontsize=21)
            ax1.set_ylabel('Non-standardized SMILES Count', fontsize=21)
            ax1.set_title('Non-standardized SMILES Count', fontsize=22)
            plotting.plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
            ax1.tick_params(labelsize=18)
            
            # Add value labels
            for bar, count in zip(bars1, non_canonical_counts):
                ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(non_canonical_counts)*0.01,
                        str(count), ha='center', va='bottom', fontsize=9)
            
            # Right plot: Ratio
            bars2 = ax2.bar(db_names, non_canonical_ratios,
                           color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax2.set_xlabel('Database', fontsize=21)
            ax2.set_ylabel('Non-standardized SMILES Ratio (%)', fontsize=21)
            ax2.set_title('Non-standardized SMILES Ratio', fontsize=22)
            plotting.plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
            ax2.tick_params(labelsize=18)
            
            # Add value labels
            for bar, ratio in zip(bars2, non_canonical_ratios):
                ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(non_canonical_ratios)*0.01,
                        f'{ratio:.1f}%', ha='center', va='bottom', fontsize=9)
            
            plotting.plt.suptitle('Non-standardized SMILES Comparison Across Databases', fontsize=23)
            
        elif mode == 'count':
            fig, ax = plotting.plt.subplots(figsize=figsize)
            bars = ax.bar(db_names, non_canonical_counts,
                          color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax.set_xlabel('Database', fontsize=21)
            ax.set_ylabel('Non-standardized SMILES Count', fontsize=21)
            ax.set_title('Non-standardized SMILES Count Comparison', fontsize=22)
            
            for bar, count in zip(bars, non_canonical_counts):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(non_canonical_counts)*0.01,
                       str(count), ha='center', va='bottom', fontsize=9)
            
            plotting.plt.xticks(rotation=45, ha='right')
            ax.tick_params(labelsize=18)
            
        else:  # mode == 'ratio'
            fig, ax = plotting.plt.subplots(figsize=figsize)
            bars = ax.bar(db_names, non_canonical_ratios,
                          color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax.set_xlabel('Database', fontsize=21)
            ax.set_ylabel('Non-standardized SMILES Ratio (%)', fontsize=21)
            ax.set_title('Non-standardized SMILES Ratio Comparison', fontsize=22)
            
            for bar, ratio in zip(bars, non_canonical_ratios):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(non_canonical_ratios)*0.01,
                       f'{ratio:.1f}%', ha='center', va='bottom', fontsize=9)
            
            plotting.plt.xticks(rotation=45, ha='right')
            ax.tick_params(labelsize=18)
        
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.comparison_dir, f'non_standardized_smiles_{mode}.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

