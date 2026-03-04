from .._common import *
from .. import plotting


class QualityPlotMixin:
    def plot_data_quality_metrics(self, figsize=(20, 10), save_path=None):
        """
        Plot a comparison of Label Consistency (unique contradictory counts) and 
        Structural Duplication (duplicate counts) for different databases
        """
        if not self.scoring_results:
            print("No scoring results available. Please add databases first.")
            return
            
        # Collect data for each database
        db_names = []
        contradictory_counts = []
        duplicate_counts = []
        
        for db_name, results in self.scoring_results.items():
            snap = results['snapshot']

            # Get label consistency data
            label_consistency_data = snap.analysis_results.get('Label Consistency', {})
            contradictory_count = label_consistency_data.get('Total molecules with contradictory labels', 0)
            
            # Get structural duplication data - in Data Consistency and Reliability
            data_consistency = snap.analysis_results.get('Data Consistency and Reliability', {})
            structural_duplication = data_consistency.get('Structural Duplication', {})
            duplicate_count = structural_duplication.get('Duplicate molecules', 0)
            
            db_names.append(db_name)
            contradictory_counts.append(contradictory_count)
            duplicate_counts.append(duplicate_count)
        
        # Create figure
        fig, (ax1, ax2) = plotting.plt.subplots(1, 2, figsize=figsize)
        
        # Plot Label Consistency (Contradictory Labels)
        ax1.bar(db_names, contradictory_counts, color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
        ax1.set_title('Label Consistency Issues', fontsize=22)
        ax1.set_xlabel('Database', fontsize=21)
        ax1.set_ylabel('Number of Molecules with Contradictory Labels', fontsize=21)
        ax1.tick_params(axis='x', rotation=45, labelsize=18)
        
        # Add value labels
        for i, v in enumerate(contradictory_counts):
            ax1.text(i, v + 0.1, str(v), ha='center', va='bottom')
        
        # Plot Structural Duplication
        ax2.bar(db_names, duplicate_counts, color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
        ax2.set_title('Structural Duplication Issues', fontsize=21)
        ax2.set_xlabel('Database', fontsize=18)
        ax2.set_ylabel('Number of Duplicate Molecules', fontsize=18)
        ax2.tick_params(axis='x', rotation=45, labelsize=18)
        
        # Add value labels
        for i, v in enumerate(duplicate_counts):
            ax2.text(i, v + 0.1, str(v), ha='center', va='bottom')
        
        # Overall title
        fig.suptitle('Data Quality Metrics Comparison', fontsize=23)
        
        plotting.plt.tight_layout()
        
        if save_path:
            plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
    def plot_quality_heatmap(self, figsize=(14, 10), save_path=None):
        """Plot heatmap of all quality metrics"""
        if not self.scoring_results:
            print("No scoring results available.")
            return

        # Determine which categories are actually present across all databases
        available_categories = {}
        for category, metrics in self.categories.items():
            # Check if this category exists in at least one database
            category_exists = False
            for db_name, results in self.scoring_results.items():
                snap = results['snapshot']
                if category in snap.scores:
                    category_exists = True
                    break

            if category_exists:
                available_categories[category] = metrics

        # Collect all metrics from available categories only
        db_names = list(self.scoring_results.keys())
        all_metrics = []
        metric_names = []

        for category, metrics in available_categories.items():
            for metric in metrics:
                metric_names.append(f"{category[:3]}_{metric[:15]}")
                all_metrics.append(metric)

        # Create matrix - use NaN for None values instead of 0
        score_matrix = np.full((len(db_names), len(all_metrics)), np.nan)

        for i, db_name in enumerate(db_names):
            snap = self.scoring_results[db_name]['snapshot']
            j = 0
            for category, metrics in available_categories.items():
                for metric in metrics:
                    score = snap.scores[category].get(metric, 0)
                    if score is None:
                        # Keep as np.nan for label conflicts
                        score_matrix[i, j] = np.nan
                    else:
                        score_matrix[i, j] = score
                    j += 1

        # Plot heatmap with masked array to handle NaN values
        fig, ax = plotting.plt.subplots(figsize=figsize)

        # Create a masked array to properly handle NaN values
        masked_matrix = np.ma.masked_invalid(score_matrix)

        # Use RdYlGn colormap but set NaN color to light gray
        cmap = plotting.plt.cm.RdYlGn.copy()
        cmap.set_bad(color='lightgray')

        im = ax.imshow(masked_matrix, cmap=cmap, aspect='auto', vmin=0, vmax=10)

        # Set ticks
        ax.set_xticks(np.arange(len(metric_names)))
        ax.set_yticks(np.arange(len(db_names)))
        ax.set_xticklabels(metric_names, rotation=45, ha='right', fontsize=18)
        ax.set_yticklabels(db_names, fontsize=18)

        # Add colorbar
        cbar = plotting.plt.colorbar(im, ax=ax)
        cbar.set_label('Score (0-10)', rotation=270, labelpad=20, fontsize=18)

        # Add text annotations - display "NA" for NaN values
        for i in range(len(db_names)):
            for j in range(len(all_metrics)):
                if np.isnan(score_matrix[i, j]):
                    # Display NA for label conflicts
                    text = ax.text(j, i, 'NA',
                                  ha="center", va="center", color="black", fontsize=14, fontweight='bold')
                else:
                    text = ax.text(j, i, f'{score_matrix[i, j]:.1f}',
                                  ha="center", va="center", color="black", fontsize=14)

        ax.set_title('Quality Metrics Heatmap', fontsize=23)
        plotting.plt.tight_layout()

        if save_path is None:
            save_path = os.path.join(self.comprehensive_dir, 'quality_metrics_heatmap.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()
    def plot_quality_waterfall(self, figsize=(14, 8), save_path=None):
        """Plot waterfall chart of quality scores"""
        if not self.scoring_results:
            print("No scoring results available.")
            return
        
        # Collect final scores
        db_scores = []
        for db_name, results in self.scoring_results.items():
            snap = results['snapshot']
            final_score = snap.scores.get('Final Adjusted Score', 0)
            db_scores.append((db_name, final_score))
        
        # Sort by score
        db_scores.sort(key=lambda x: x[1], reverse=True)
        
        fig, ax = plotting.plt.subplots(figsize=figsize)
        
        # Create waterfall
        x = np.arange(len(db_scores))
        names, scores = zip(*db_scores)
        
        # Color based on score
        colors = []
        for score in scores:
            if score >= 80:
                colors.append('green')
            elif score >= 60:
                colors.append('yellow')
            elif score >= 40:
                colors.append('orange')
            else:
                colors.append('red')
        
        bars = ax.bar(x, scores, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
        
        # Add value labels
        for bar, score in zip(bars, scores):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                   f'{score:.1f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        # Add reference lines
        ax.axhline(y=80, color='green', linestyle='--', alpha=0.5, label='Excellent (≥80)')
        ax.axhline(y=60, color='yellow', linestyle='--', alpha=0.5, label='Good (≥60)')
        ax.axhline(y=40, color='orange', linestyle='--', alpha=0.5, label='Fair (≥40)')
        
        ax.set_xlabel('Database', fontsize=12)
        ax.set_ylabel('Final Adjusted Score', fontsize=12)
        ax.set_title('Database Quality Scores Waterfall', fontsize=16, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=45, ha='right')
        ax.set_ylim(0, 100)
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3, axis='y')
        
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.comprehensive_dir, 'quality_waterfall.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()
    def plot_before_after_cleaning(self, figsize=(14, 8), save_path=None):
        """Plot comparison before and after data cleaning"""
        # This would require running the cleaning function
        # Placeholder for now
        print("Before/after cleaning comparison requires cleaned datasets.")
        pass

