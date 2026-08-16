from .._common import *
from .. import plotting
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch


QUALITY_HEATMAP_CATEGORY_BAR_COLOR_SCHEMES = {
    "soft": {
        "Structural Integrity": "#adcef0",
        "Data Quality": "#afe2c3",
        "Experimental Information Quality": "#d6bae8",
        "Chemical Space Coverage": "#fae293",
        "Data Distribution": "#fdc9a1",
    },
    "deep": {
        "Structural Integrity": "#41719c",
        "Data Quality": "#548235",
        "Experimental Information Quality": "#7030a0",
        "Chemical Space Coverage": "#bf9000",
        "Data Distribution": "#c55a11",
    },
}

QUALITY_HEATMAP_TOP_AXIS_LABELS = {
    "Structural Integrity": "Structural\nIntegrity",
    "Data Quality": "Data\nQuality",
    "Experimental Information Quality": "Experimental\nInformation\nQuality",
    "Chemical Space Coverage": "Chemical Space\nCoverage",
    "Data Distribution": "Data\nDistribution",
}

QUALITY_HEATMAP_BAR_REFERENCE_DB_COUNT = 11
QUALITY_HEATMAP_BAR_REFERENCE_HEIGHT = 0.22
QUALITY_HEATMAP_BAR_REFERENCE_GAP = 0.03
QUALITY_HEATMAP_BAR_REFERENCE_ROUNDING = 0.10


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
            scorer = results['scorer']
            
            # Get label consistency data
            label_consistency_data = scorer.analysis_results.get('Label Consistency', {})
            contradictory_count = label_consistency_data.get('Total molecules with contradictory labels', 0)
            
            # Get structural duplication data - in Data Consistency and Reliability
            data_consistency = scorer.analysis_results.get('Data Consistency and Reliability', {})
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
            ax1.text(i, v + 0.1, str(v), ha='center', va='bottom', fontsize=9)
        
        # Plot Structural Duplication
        ax2.bar(db_names, duplicate_counts, color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
        ax2.set_title('Structural Duplication Issues', fontsize=22)
        ax2.set_xlabel('Database', fontsize=21)
        ax2.set_ylabel('Number of Duplicate Molecules', fontsize=21)
        ax2.tick_params(axis='x', rotation=45, labelsize=18)
        
        # Add value labels
        for i, v in enumerate(duplicate_counts):
            ax2.text(i, v + 0.1, str(v), ha='center', va='bottom', fontsize=9)
        
        # Overall title
        fig.suptitle('Data Quality Metrics Comparison', fontsize=23)
        
        plotting.plt.tight_layout()
        
        if save_path:
            plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
    def plot_quality_heatmap(
        self,
        figsize=(14, 10),
        save_path=None,
        category_bar_color_scheme="soft",
        top_label_size=20,
        title_pad=34,
    ):
        """Plot heatmap of all quality metrics"""
        if not self.scoring_results:
            print("No scoring results available.")
            return

        if category_bar_color_scheme not in QUALITY_HEATMAP_CATEGORY_BAR_COLOR_SCHEMES:
            valid = ", ".join(sorted(QUALITY_HEATMAP_CATEGORY_BAR_COLOR_SCHEMES))
            raise ValueError(
                f"Unknown category_bar_color_scheme '{category_bar_color_scheme}'. "
                f"Available: {valid}"
            )

        # Determine which categories are actually present across all databases
        available_categories = {}
        for category, metrics in self.categories.items():
            # Check if this category exists in at least one database
            category_exists = False
            for db_name, results in self.scoring_results.items():
                scorer = results['scorer']
                if category in scorer.scores:
                    category_exists = True
                    break

            if category_exists:
                available_categories[category] = metrics

        # Collect all metrics from available categories only
        db_names = list(self.scoring_results.keys())
        all_metrics = []
        metric_labels = []
        category_spans = []
        start = 0

        for category, metrics in available_categories.items():
            labels = self.categories_labels.get(category, metrics)
            for metric in metrics:
                all_metrics.append(metric)
            metric_labels.extend(labels[:len(metrics)])
            end = start + len(metrics)
            category_spans.append((category, start, end))
            start = end

        # Create matrix - use NaN for None values instead of 0
        score_matrix = np.full((len(db_names), len(all_metrics)), np.nan)

        for i, db_name in enumerate(db_names):
            scorer = self.scoring_results[db_name]['scorer']
            j = 0
            for category, metrics in available_categories.items():
                category_scores = scorer.scores.get(category, {})
                for metric in metrics:
                    score = category_scores.get(metric)
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

        # Use red-white-green colormap but set NaN color to light gray
        cmap = LinearSegmentedColormap.from_list(
            "RdWhGnThreeColor",
            ["#b2182b", "#ffffff", "#1a9850"],
        )
        cmap.set_bad(color='lightgray')

        im = ax.imshow(masked_matrix, cmap=cmap, aspect='auto', vmin=0, vmax=10)

        # Set ticks
        ax.set_xticks(np.arange(len(metric_labels)))
        ax.set_yticks(np.arange(len(db_names)))
        ax.set_xticklabels(
            metric_labels,
            rotation=50,
            ha='right',
            rotation_mode='anchor',
            fontsize=18,
            linespacing=1.0,
        )
        ax.set_yticklabels(db_names, fontsize=18)
        ax.tick_params(axis='x', length=7, width=1, direction='out', pad=2)
        ax.tick_params(axis='y', length=7, width=1, direction='out')

        for tick_label in ax.get_xticklabels():
            tick_label.set_color('black')
            tick_label.set_fontweight('normal')

        top_ax = ax.secondary_xaxis('top')
        top_ax.set_xticks([(cat_start + cat_end - 1) / 2 for _, cat_start, cat_end in category_spans])
        top_ax.set_xticklabels(
            [QUALITY_HEATMAP_TOP_AXIS_LABELS[category] for category, _, _ in category_spans],
            fontsize=top_label_size,
            fontweight='normal',
        )
        top_ax.tick_params(length=0, pad=10)
        for spine in top_ax.spines.values():
            spine.set_visible(False)

        for top_label in top_ax.get_xticklabels():
            top_label.set_color('black')
            top_label.set_fontweight('normal')

        bar_scale = max(len(db_names), 1) / QUALITY_HEATMAP_BAR_REFERENCE_DB_COUNT
        bar_height = QUALITY_HEATMAP_BAR_REFERENCE_HEIGHT * bar_scale
        bar_gap = QUALITY_HEATMAP_BAR_REFERENCE_GAP * bar_scale
        bar_rounding = QUALITY_HEATMAP_BAR_REFERENCE_ROUNDING * bar_scale
        bar_y = -0.5 - bar_gap - bar_height

        ax.spines['top'].set_visible(True)
        ax.spines['top'].set_linewidth(ax.spines['left'].get_linewidth())
        category_bar_colors = QUALITY_HEATMAP_CATEGORY_BAR_COLOR_SCHEMES[category_bar_color_scheme]
        for category, cat_start, cat_end in category_spans:
            ax.add_patch(
                FancyBboxPatch(
                    (cat_start - 0.46, bar_y),
                    (cat_end - cat_start) - 0.08,
                    bar_height,
                    boxstyle=f"round,pad=0,rounding_size={bar_rounding:.6f}",
                    facecolor=category_bar_colors[category],
                    edgecolor='none',
                    linewidth=0,
                    clip_on=False,
                )
            )

        # Add colorbar
        cbar = plotting.plt.colorbar(im, ax=ax, pad=0.02)
        cbar.set_label('Score (0-10)', rotation=270, labelpad=20, fontsize=20)
        cbar.ax.tick_params(labelsize=18)

        # Add text annotations - display "NA" for NaN values
        for i in range(len(db_names)):
            for j in range(len(all_metrics)):
                if np.isnan(score_matrix[i, j]):
                    # Display NA for label conflicts
                    text = ax.text(j, i, 'NA',
                                  ha="center", va="center", color="black", fontsize=14, fontweight='normal')
                else:
                    text = ax.text(j, i, f'{score_matrix[i, j]:.1f}',
                                  ha="center", va="center", color="black", fontsize=14)

        ax.set_title('Quality Metrics Heatmap', fontsize=26, pad=title_pad)
        plotting.plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.965))

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
            scorer = results['scorer']
            final_score = scorer.scores.get('Final Adjusted Score', 0)
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
                   f'{score:.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
        
        # Add reference lines
        ax.axhline(y=80, color='green', linestyle='--', alpha=0.5, label='Excellent (≥80)')
        ax.axhline(y=60, color='yellow', linestyle='--', alpha=0.5, label='Good (≥60)')
        ax.axhline(y=40, color='orange', linestyle='--', alpha=0.5, label='Fair (≥40)')
        
        ax.set_xlabel('Database', fontsize=21)
        ax.set_ylabel('Final Adjusted Score', fontsize=21)
        ax.set_title('Database Quality Scores Waterfall', fontsize=22)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=45, ha='right')
        ax.set_ylim(0, 100)
        ax.tick_params(labelsize=18)
        ax.legend(loc='best', fontsize=16)
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
