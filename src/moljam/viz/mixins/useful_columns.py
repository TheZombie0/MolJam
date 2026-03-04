from .._common import *
from .. import plotting


class UsefulColumnsPlotMixin:
    def plot_useful_columns_analysis(self, db_name, figsize=(18, 10), save_path=None):
        """
        Detailed analysis of useful columns for a specific database.
        Shows: column types, coverage rates, and confidence scores.
        """
        if db_name not in self.scoring_results:
            print(f"Database '{db_name}' not found in scoring results")
            return

        snap = self.scoring_results[db_name]['snapshot']

        if 'Experimental Information Quality' not in snap.analysis_results:
            print(f"Experimental Information Quality analysis not available for '{db_name}'")
            return

        analysis = snap.analysis_results['Experimental Information Quality']
        useful_details = analysis.get('Useful column details', [])

        if not useful_details:
            print(f"No useful column details available for '{db_name}'")
            return

        fig, axes = plotting.plt.subplots(2, 2, figsize=figsize)

        # 1. Top-left: Type distribution bar chart
        ax1 = axes[0, 0]
        types_found = analysis.get('Types found', [])
        type_counts = {'activity': 0, 'label': 0, 'experimental': 0}

        type_keywords = {
            'activity': ['Activity', 'Measurement', 'Predicted', 'Calculated'],
            'label': ['label', 'categories', 'Classification'],
            'experimental': ['Experimental', 'condition']
        }

        for col in useful_details:
            reason = col.get('reason', '')
            for type_name, keywords in type_keywords.items():
                if any(kw.lower() in reason.lower() for kw in keywords):
                    type_counts[type_name] += 1
                    break

        type_names = list(type_counts.keys())
        type_values = list(type_counts.values())
        type_colors = ['#3498db', '#e74c3c', '#f39c12']

        bars = ax1.bar(type_names, type_values, color=type_colors, edgecolor='white')
        for bar, val in zip(bars, type_values):
            if val > 0:
                ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                         str(val), ha='center', va='bottom', fontweight='bold')

        ax1.set_ylabel('Number of Columns', fontsize=11)
        ax1.set_title('Useful Columns by Type', fontsize=13, fontweight='bold')
        ax1.set_ylim(0, max(type_values) * 1.2 if max(type_values) > 0 else 1)

        # 2. Top-right: Coverage and confidence scatter
        ax2 = axes[0, 1]
        col_names = [col['name'][:15] + '...' if len(col['name']) > 15 else col['name']
                     for col in useful_details]
        confidences = [col.get('confidence', 0) for col in useful_details]

        # Calculate coverage for each column
        df = snap.df
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
        ax2.set_ylabel('Classification Confidence', fontsize=11)
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
        ax3.set_title('Useful Columns Coverage (Top 10)', fontsize=13, fontweight='bold')
        ax3.set_xlim(0, 1.1)
        ax3.axvline(x=0.9, color='green', linestyle='--', alpha=0.5)

        # 4. Bottom-right: Summary statistics
        ax4 = axes[1, 1]
        ax4.axis('off')

        summary_text = f"""
        Summary Statistics
        ══════════════════════════════════

        Total Columns Analyzed: {analysis.get('Total columns analyzed', 'N/A')}

        Useful Columns: {analysis.get('Useful columns', 'N/A')}
        Excluded Columns: {analysis.get('Excluded columns', 'N/A')}
        Unknown Columns: {analysis.get('Unknown columns', 'N/A')}

        Average Coverage: {analysis.get('Average coverage of useful columns', 'N/A')}
        Useful Column Ratio: {analysis.get('Useful column ratio', 'N/A')}
        Average Confidence: {analysis.get('Average classification confidence', 'N/A')}

        Types Found: {', '.join(analysis.get('Types found', [])) or 'None'}
        """

        ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, fontsize=12,
                 verticalalignment='top', fontfamily='monospace',
                 bbox=dict(boxstyle='round', facecolor='#f8f9fa', edgecolor='#dee2e6'))

        plotting.plt.suptitle(f'{db_name} - Useful Columns Analysis', fontsize=16, fontweight='bold', y=1.02)
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

