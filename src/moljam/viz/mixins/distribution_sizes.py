from .._common import *
from .. import plotting


class DatabaseSizePlotMixin:
    def plot_database_size_comparison(self, figsize=(14, 8), save_path=None):
        """Plot comparison of database sizes"""
        if not self.scoring_results:
            print("No scoring results available.")
            return
        
        db_names = list(self.scoring_results.keys())
        sizes = []
        
        for db_name, results in self.scoring_results.items():
            scorer = results['scorer']
            sizes.append(scorer.num_molecules)
        
        fig, ax = plotting.plt.subplots(figsize=figsize)
        
        bars = ax.bar(db_names, sizes,
                     color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
        
        ax.set_xlabel('Database', fontsize=21)
        ax.set_ylabel('Number of Molecules', fontsize=21)
        ax.set_title('Database Size Comparison', fontsize=23)
        plotting.plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        ax.tick_params(labelsize=18)
        
        # Add value labels
        for bar, size in zip(bars, sizes):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(sizes)*0.01,
                   f'{size:,}', ha='center', va='bottom', fontsize=9)
        
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.distribution_dir, 'database_size_comparison.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()
