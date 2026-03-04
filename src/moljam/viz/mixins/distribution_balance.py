from .._common import *
from .. import plotting


class DistributionBalancePlotMixin:
    def plot_data_imbalance_lorenz(self, figsize=(16, 16), save_path=None):
        """Plot Lorenz curve for data imbalance"""
        if not self.scoring_results:
            print("No scoring results available.")
            return
        
        fig, ax = plotting.plt.subplots(figsize=figsize)
        
        for db_idx, (db_name, results) in enumerate(self.scoring_results.items()):
            snap = results['snapshot']
            class_cols = results['class_cols']

            if class_cols and class_cols[0] in snap.df.columns:
                # Get class distribution
                class_counts = snap.df[class_cols[0]].value_counts().sort_values()
                
                if len(class_counts) > 1:
                    # Calculate cumulative proportions
                    cumsum = class_counts.cumsum() / class_counts.sum()
                    x = np.arange(len(cumsum)) / (len(cumsum) - 1)
                    
                    # Plot Lorenz curve
                    ax.plot(x, cumsum.values,
                           color=self.colors[db_idx % len(self.colors)],
                           linewidth=2, label=db_name, alpha=0.8)
        
        # Add line of equality
        # ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Perfect Balance')
        
        ax.set_xlabel('Cumulative Proportion of Classes', fontsize=21)
        ax.set_ylabel('Cumulative Proportion of Samples', fontsize=21)
        ax.set_title('Data Imbalance Lorenz Curves', fontsize=23)
        # ax.legend(loc='best', fontsize=16)
        ax.tick_params(labelsize=18)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.distribution_dir, 'data_imbalance_lorenz.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

    def plot_class_entropy_comparison(self, figsize=(12, 6), save_path=None):
        """Plot entropy of class distributions"""
        if not self.scoring_results:
            print("No scoring results available.")
            return
        
        db_names = []
        entropies = []
        normalized_entropies = []
        
        for db_name, results in self.scoring_results.items():
            snap = results['snapshot']
            class_cols = results['class_cols']

            if class_cols and class_cols[0] in snap.df.columns:
                # Get class distribution
                class_counts = snap.df[class_cols[0]].value_counts()
                
                if len(class_counts) > 1:
                    # Calculate entropy
                    proportions = class_counts / class_counts.sum()
                    entropy = -np.sum(proportions * np.log2(proportions + 1e-10))
                    max_entropy = np.log2(len(class_counts))
                    normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0
                    
                    db_names.append(db_name)
                    entropies.append(entropy)
                    normalized_entropies.append(normalized_entropy)
        
        if not db_names:
            print("No class distribution data available.")
            return
        
        fig, (ax1, ax2) = plotting.plt.subplots(1, 2, figsize=figsize)
        
        # Raw entropy
        bars1 = ax1.bar(db_names, entropies,
                       color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
        ax1.set_xlabel('Database', fontsize=21)
        ax1.set_ylabel('Entropy (bits)', fontsize=21)
        ax1.set_title('Class Distribution Entropy', fontsize=23)
        plotting.plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
        ax1.tick_params(labelsize=18)
        
        for bar, ent in zip(bars1, entropies):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(entropies)*0.01,
                    f'{ent:.2f}', ha='center', va='bottom', fontsize=9)
        
        # Normalized entropy
        bars2 = ax2.bar(db_names, normalized_entropies,
                       color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
        ax2.set_xlabel('Database', fontsize=21)
        ax2.set_ylabel('Normalized Entropy', fontsize=21)
        ax2.set_title('Normalized Class Distribution Entropy', fontsize=22)
        plotting.plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
        ax2.set_ylim(0, 1)
        ax2.tick_params(labelsize=18)
        
        for bar, ent in zip(bars2, normalized_entropies):
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f'{ent:.3f}', ha='center', va='bottom', fontsize=9)
        
        # Add reference line for perfect balance
        ax2.axhline(y=1.0, color='r', linestyle='--', alpha=0.5, label='Perfect Balance')
        ax2.legend()
        
        plotting.plt.suptitle('Class Distribution Entropy Analysis', fontsize=23)
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.distribution_dir, 'class_entropy_comparison.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

