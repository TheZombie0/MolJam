from .._common import *
from .. import plotting


class AnnotationQualityPlotMixin:
    def plot_annotation_quality_breakdown(self, db_name, figsize=(16, 8), save_path=None):
        """
        Plot Experimental Information Quality breakdown for a specific database.
        Shows: 3 sub-metrics scores and column classification distribution.
        """
        if db_name not in self.scoring_results:
            print(f"Database '{db_name}' not found in scoring results")
            return

        scorer = self.scoring_results[db_name]['scorer']

        if 'Experimental Information Quality' not in scorer.scores:
            print(f"Experimental Information Quality not available for '{db_name}' (experimental_info=False)")
            return

        fig, axes = plotting.plt.subplots(1, 2, figsize=figsize)

        # Left plot: Sub-metrics bar chart
        ax1 = axes[0]
        metrics = ['Time Label\nAvailability', 'Annotation Support\nQuality', 'Type\nDiversity']
        scores = [
            scorer.scores['Experimental Information Quality'].get('Time Label Availability', 0),
            scorer.scores['Experimental Information Quality'].get('Annotation Support Quality', 0),
            scorer.scores['Experimental Information Quality'].get('Type Diversity', 0)
        ]

        colors = ['#45B7D1', '#6BC5D8', '#8CD3E5']
        bars = ax1.barh(metrics, scores, color=colors, edgecolor='white', height=0.6)

        # Add score labels
        for bar, score in zip(bars, scores):
            ax1.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
                     f'{score:.2f}/10', va='center', fontsize=11, fontweight='bold')

        ax1.set_xlim(0, 12)
        ax1.set_xlabel('Score', fontsize=12)
        ax1.set_title(f'Experimental Information Quality Breakdown\n(Total: {sum(scores):.2f}/30)', fontsize=14, fontweight='bold')
        ax1.axvline(x=10, color='gray', linestyle='--', alpha=0.5, label='Max Score')

        # Right plot: Column classification pie chart
        ax2 = axes[1]
        if 'Experimental Information Quality' in scorer.analysis_results:
            analysis = scorer.analysis_results['Experimental Information Quality']
            accepted = analysis.get('Accepted columns', analysis.get('Support columns', analysis.get('Useful columns', 0)))
            excluded = analysis.get('Excluded columns', 0)
            derived = analysis.get('Derived/Predicted columns', 0)
            unknown = analysis.get('Unknown columns', 0)

            pie_data = [
                ('Accepted', accepted, '#2ecc71'),
                ('Excluded', excluded, '#e74c3c'),
                ('Derived/Predicted', derived, '#f39c12'),
                ('Unknown', unknown, '#95a5a6'),
            ]
            pie_data = [item for item in pie_data if item[1] > 0]

            if pie_data:
                sizes = [item[1] for item in pie_data]
                labels = [f'{item[0]}\n({item[1]})' for item in pie_data]
                pie_colors = [item[2] for item in pie_data]
                explode = [0.05 if item[0] == 'Accepted' else 0 for item in pie_data]

                ax2.pie(
                    sizes,
                    labels=labels,
                    colors=pie_colors,
                    explode=explode,
                    autopct='%1.1f%%',
                    startangle=90,
                    textprops={'fontsize': 11},
                )
                ax2.set_title('Column Classification Distribution', fontsize=14, fontweight='bold')

                # Add types found info
                types_found = analysis.get('Types found', [])
                if types_found:
                    types_text = f"Types found: {', '.join(types_found)}"
                    ax2.text(0, -1.3, types_text, ha='center', fontsize=10,
                             style='italic', transform=ax2.transAxes)
            else:
                ax2.text(0.5, 0.5, 'No column data available', ha='center', va='center',
                         fontsize=14, transform=ax2.transAxes)
                ax2.set_title('Column Classification Distribution', fontsize=14, fontweight='bold')
        else:
            ax2.text(0.5, 0.5, 'Analysis results not available', ha='center', va='center',
                     fontsize=14, transform=ax2.transAxes)

        plotting.plt.suptitle(f'{db_name} - Experimental Information Quality Analysis', fontsize=16, fontweight='bold', y=1.02)
        plotting.plt.tight_layout()

        if save_path:
            plotting.plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"Saved annotation quality breakdown to {save_path}")
        else:
            # Save to database-specific directory
            db_dir = self.scoring_results[db_name]['output_dir']
            save_path = os.path.join(db_dir, f'{db_name}_annotation_quality_breakdown.png')
            plotting.plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"Saved annotation quality breakdown to {save_path}")

        plotting.plt.close()
