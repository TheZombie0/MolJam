from .._common import *
from .. import plotting


class UsefulColumnsPlotMixin:
    def plot_useful_columns_analysis(self, db_name, figsize=(18, 10), save_path=None):
        """
        Detailed analysis of annotation support columns for a specific database.
        Shows: column types, coverage rates, and confidence scores.
        """
        if db_name not in self.scoring_results:
            print(f"Database '{db_name}' not found in scoring results")
            return

        scorer = self.scoring_results[db_name]['scorer']

        if 'Experimental Information Quality' not in scorer.analysis_results:
            print(f"Experimental Information Quality analysis not available for '{db_name}'")
            return

        analysis = scorer.analysis_results['Experimental Information Quality']
        useful_details = analysis.get(
            'Accepted column details',
            analysis.get('Support column details', analysis.get('Useful column details', []))
        )

        if not useful_details:
            print(f"No useful column details available for '{db_name}'")
            return

        fig, axes = plotting.plt.subplots(2, 2, figsize=figsize)

        # 1. Top-left: Type distribution bar chart
        ax1 = axes[0, 0]
        role_counts = analysis.get('Role counts', {})
        type_names = ['time', 'activity', 'label', 'experimental_context']
        type_values = [role_counts.get(name, 0) for name in type_names]
        type_colors = ['#9b59b6', '#3498db', '#e74c3c', '#f39c12']

        bars = ax1.bar(type_names, type_values, color=type_colors, edgecolor='white')
        for bar, val in zip(bars, type_values):
            if val > 0:
                ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                         str(val), ha='center', va='bottom', fontweight='bold')

        ax1.set_ylabel('Number of Columns', fontsize=11)
        ax1.set_title('Accepted Annotation Columns by Type', fontsize=13, fontweight='bold')
        ax1.set_ylim(0, max(type_values) * 1.2 if max(type_values) > 0 else 1)

        # 2. Top-right: Coverage and confidence scatter
        ax2 = axes[0, 1]
        col_names = [col['name'][:15] + '...' if len(col['name']) > 15 else col['name']
                     for col in useful_details]
        confidences = [col.get('confidence', 0) for col in useful_details]

        # Calculate coverage for each column
        df = scorer.df
        coverages = []
        for col in useful_details:
            col_name = col['name']
            if col_name in df.columns:
                coverage = df[col_name].notna().sum() / len(df)
                coverages.append(coverage)
            else:
                coverages.append(0)

        scatter = ax2.scatter(coverages, confidences, c=range(len(col_names)),
                              cmap='viridis', s=100, alpha=0.7, edgecolors='white')
        ax2.set_xlabel('Data Coverage', fontsize=11)
        ax2.set_ylabel('Column Confidence', fontsize=11)
        ax2.set_title('Coverage vs Confidence', fontsize=13, fontweight='bold')
        ax2.set_xlim(-0.05, 1.05)
        ax2.set_ylim(-0.05, 1.05)
        ax2.axhline(y=0.8, color='green', linestyle='--', alpha=0.5, label='High Confidence')
        ax2.axvline(x=0.9, color='blue', linestyle='--', alpha=0.5, label='High Coverage')
        ax2.legend(loc='lower left', fontsize=9)

        # 3. Bottom-left: Column coverage bar chart
        ax3 = axes[1, 0]
        if len(col_names) > 10:
            col_names_display = col_names[:10]
            coverages_display = coverages[:10]
        else:
            col_names_display = col_names
            coverages_display = coverages

        colors = ['#2ecc71' if c >= 0.9 else '#f39c12' if c >= 0.7 else '#e74c3c'
                  for c in coverages_display]
        bars = ax3.barh(col_names_display, coverages_display, color=colors, edgecolor='white')
        ax3.set_xlabel('Coverage Rate', fontsize=11)
        ax3.set_title('Accepted Annotation Column Coverage (Top 10)', fontsize=13, fontweight='bold')
        ax3.set_xlim(0, 1.1)
        ax3.axvline(x=0.9, color='green', linestyle='--', alpha=0.5)

        # 4. Bottom-right: Summary statistics
        ax4 = axes[1, 1]
        ax4.axis('off')

        summary_text = f"""
        Summary Statistics
        ══════════════════════════════════

        Total Columns Analyzed: {analysis.get('Total columns analyzed', 'N/A')}

        Accepted Columns: {analysis.get('Accepted columns', analysis.get('Support columns', analysis.get('Useful columns', 'N/A')))}
        Excluded Columns: {analysis.get('Excluded columns', 'N/A')}
        Derived Columns: {analysis.get('Derived/Predicted columns', 'N/A')}
        Unknown Columns: {analysis.get('Unknown columns', 'N/A')}

        Average Coverage: {analysis.get('Average accepted coverage', 'N/A')}
        Average Confidence: {analysis.get('Average accepted confidence', 'N/A')}
        Agreement Rate: {analysis.get('Agreement rate', 'N/A')}
        Conflict Rate: {analysis.get('Conflict rate', 'N/A')}

        Types Found: {', '.join(analysis.get('Types found', [])) or 'None'}
        """

        ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, fontsize=12,
                 verticalalignment='top', fontfamily='monospace',
                 bbox=dict(boxstyle='round', facecolor='#f8f9fa', edgecolor='#dee2e6'))

        plotting.plt.suptitle(f'{db_name} - Annotation Support Analysis', fontsize=16, fontweight='bold', y=1.02)
        plotting.plt.tight_layout()

        if save_path:
            plotting.plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"Saved useful columns analysis to {save_path}")
        else:
            db_dir = self.scoring_results[db_name]['output_dir']
            save_path = os.path.join(db_dir, f'{db_name}_useful_columns_analysis.png')
            plotting.plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"Saved useful columns analysis to {save_path}")

        plotting.plt.close()
