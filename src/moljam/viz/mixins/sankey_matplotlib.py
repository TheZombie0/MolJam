from .._common import *
from .. import plotting


class SankeyMatplotlibFallbackMixin:
    def plot_sankey_for_database_matplotlib_fallback(self, db_name, figsize=(16, 12), save_path=None):
        """
        Plot Sankey diagram for a single database showing score and penalty flows

        Args:
            db_name: Database name
            figsize: Figure size
            save_path: Path to save the figure
        """
        if db_name not in self.scoring_results:
            print(f"Database {db_name} not found in scoring results.")
            return

        scorer = self.scoring_results[db_name]['scorer']

        # Get metric penalties
        metric_penalties = scorer.scores.get('Metric Penalties', {})

        # Prepare data structures
        flows_data = []

        # Define category colors (base colors for each category)
        category_colors = {
            'Structural Integrity': '#3498db',  # Blue
            'Data Quality': '#2ecc71',  # Green
            'Experimental Information Quality': '#f39c12',  # Orange
            'Chemical Space Coverage': '#9b59b6',  # Purple
            'Data Distribution': '#e74c3c'  # Red
        }

        # Collect all metrics data
        for category, metrics in self.categories.items():
            if category not in scorer.scores:
                continue

            base_color = category_colors.get(category, '#95a5a6')

            for idx, metric in enumerate(metrics):
                score = scorer.scores[category].get(metric, 0)
                if score is None:
                    score = 0

                # Calculate deduction (max 10 - actual score)
                deduction = 10 - score

                # Get penalty for this metric
                metric_key = f"{category}::{metric}"
                penalty = metric_penalties.get(metric_key, 0)

                # Calculate alpha based on metric index in category
                alpha = 1.0 - (idx * 0.25)  # Decrease alpha for each metric
                alpha = max(0.3, alpha)  # Minimum alpha

                flows_data.append({
                    'category': category,
                    'metric': metric,
                    'score': score,
                    'deduction': deduction,
                    'penalty': penalty,
                    'color': base_color,
                    'alpha': alpha
                })

        # Create figure
        fig, ax = plotting.plt.subplots(figsize=figsize)

        # Use plotly-style Sankey since matplotlib's Sankey is limited
        # We'll create a custom visualization using patches

        # Calculate positions
        left_x = 0.1
        right_x = 0.9

        # Group by category for positioning
        category_groups = {}
        for flow in flows_data:
            cat = flow['category']
            if cat not in category_groups:
                category_groups[cat] = []
            category_groups[cat].append(flow)

        # Calculate vertical positions
        total_metrics = len(flows_data)
        y_positions = {}
        current_y = 0.95

        for cat, flows in category_groups.items():
            for flow in flows:
                y_positions[f"{cat}::{flow['metric']}"] = current_y
                current_y -= (0.9 / total_metrics)

        # Right side nodes
        score_y = 0.7
        penalty_y = 0.3

        # Draw left side nodes (metrics)
        for flow in flows_data:
            metric_key = f"{flow['category']}::{flow['metric']}"
            y = y_positions[metric_key]

            # Draw metric box
            from matplotlib.patches import FancyBboxPatch
            from matplotlib import colors as mcolors

            color_rgba = mcolors.to_rgba(flow['color'], flow['alpha'])

            box = FancyBboxPatch((left_x - 0.05, y - 0.01), 0.05, 0.02,
                                boxstyle="round,pad=0.002",
                                facecolor=color_rgba,
                                edgecolor='black', linewidth=0.5)
            ax.add_patch(box)

            # Add metric label
            metric_label = flow['metric'][:25] + '...' if len(flow['metric']) > 25 else flow['metric']
            ax.text(left_x - 0.06, y, metric_label,
                   ha='center', va='center', fontsize=16)

        # Draw right side nodes
        # Score node
        score_box = FancyBboxPatch((right_x, score_y - 0.05), 0.05, 0.1,
                                  boxstyle="round,pad=0.005",
                                  facecolor='#27ae60', edgecolor='black', linewidth=1)
        ax.add_patch(score_box)
        ax.text(right_x + 0.06, score_y, 'Actual\nScores',
               ha='left', va='center', fontsize=16, fontweight='bold')

        # Penalty node
        penalty_box = FancyBboxPatch((right_x, penalty_y - 0.05), 0.05, 0.1,
                                    boxstyle="round,pad=0.005",
                                    facecolor='#c0392b', edgecolor='black', linewidth=1)
        ax.add_patch(penalty_box)
        ax.text(right_x + 0.06, penalty_y, 'Penalty\nDeductions',
               ha='left', va='center', fontsize=16, fontweight='bold')

        # Draw flows
        from matplotlib.patches import FancyArrowPatch
        from matplotlib.patches import PathPatch
        from matplotlib.path import Path

        for flow in flows_data:
            metric_key = f"{flow['category']}::{flow['metric']}"
            y = y_positions[metric_key]

            color_rgba = mcolors.to_rgba(flow['color'], flow['alpha'] * 0.3)

            # Flow to score (actual score value)
            if flow['score'] > 0:
                # Bezier curve to score node
                verts = [
                    (left_x, y),  # Start at metric
                    (left_x + 0.2, y),  # Control point 1
                    (right_x - 0.2, score_y),  # Control point 2
                    (right_x, score_y)  # End at score node
                ]
                codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]

                path = Path(verts, codes)
                # Width based on score value
                linewidth = max(0.5, flow['score'] * 0.8)
                patch = PathPatch(path, facecolor='none', edgecolor=color_rgba,
                                linewidth=linewidth, alpha=0.6)
                ax.add_patch(patch)

            # Flow to penalty (if penalty > 0)
            if flow['penalty'] > 0:
                verts = [
                    (left_x, y),
                    (left_x + 0.2, y),
                    (right_x - 0.2, penalty_y),
                    (right_x, penalty_y)
                ]
                codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]

                path = Path(verts, codes)
                linewidth = max(0.5, flow['penalty'] * 2)
                penalty_color = mcolors.to_rgba('#e74c3c', 0.4)
                patch = PathPatch(path, facecolor='none', edgecolor=penalty_color,
                                linewidth=linewidth, alpha=0.7)
                ax.add_patch(patch)

        # Add category legend
        legend_elements = []
        for cat, color in category_colors.items():
            legend_elements.append(plotting.mpatches.Patch(facecolor=color, edgecolor='black',
                                                 label=cat, alpha=0.7))

        ax.legend(handles=legend_elements, loc='upper center',
                 bbox_to_anchor=(0.5, -0.05), ncol=3, frameon=True,
                 fancybox=True, fontsize=9)

        # Set axis properties
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')

        # Add title
        total_score = scorer.scores.get('Normalized Score', 0)
        final_score = scorer.scores.get('Final Adjusted Score', 0)
        penalty_total = scorer.analysis_results.get('Low Score Penalty', {}).get('Total penalty applied', '0')

        ax.text(0.5, 0.98, f'Score Flow Diagram - {db_name}',
               ha='right', va='top', fontsize=16, fontweight='bold')
        ax.text(0.5, 0.94, f'Normalized Score: {total_score:.2f}/100  |  Penalty: -{penalty_total}  |  Final: {final_score:.2f}/100',
               ha='right', va='top', fontsize=11)

        plotting.plt.tight_layout()

        if save_path is None:
            save_path = os.path.join(self.scoring_results[db_name]['output_dir'],
                                   f'sankey_score_flow.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

        # print(f"Sankey diagram saved to: {save_path}")

