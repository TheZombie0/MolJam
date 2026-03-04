from .._common import *
from .. import plotting


class RadarPlotMixin:
    def plot_radar(self, figsize=(11, 11), save_path=None):
        """Plot a radar chart showing normalized scores for each category"""
        if not self.scoring_results:
            print("No scoring results available. Please add databases first.")
            return

        fig, ax = plotting.plt.subplots(figsize=figsize, subplot_kw=dict(projection='polar'))

        # Determine which categories are actually present across all databases
        available_categories = []
        for category in self.categories_labels.keys():
            # Check if this category exists in at least one database
            category_exists = False
            for db_name, results in self.scoring_results.items():
                snap = results['snapshot']
                if category in snap.scores:
                    category_exists = True
                    break

            if category_exists:
                available_categories.append(category)

        # Categories for radar plot (4 or 5 dimensions depending on available categories)
        categories = available_categories
        n_categories = len(categories)

        # Calculate angles for each category
        angles = [n / n_categories * 2 * pi for n in range(n_categories)]
        angles += angles[:1]  # Complete the circle

        # Plot data for each database
        for db_idx, (db_name, results) in enumerate(self.scoring_results.items()):
            snap = results['snapshot']

            # Calculate normalized scores for each category (0-100 scale)
            category_scores = []
            for category in categories:
                # Get total score for the category
                if category in snap.scores:
                    total = snap.scores[category].get("Normalized Total", 0)
                    category_scores.append(total)
                else:
                    category_scores.append(0)

            # Complete the circle
            category_scores += category_scores[:1]

            # Plot
            color = self.colors[db_idx % len(self.colors)]
            ax.plot(angles, category_scores, 'o-', linewidth=2,
                   label=f'{db_name} ({snap.scores["Normalized Score"]:.1f}/100)',
                   color=color)
            ax.fill(angles, category_scores, alpha=0.15, color=color)

        # Customize the radar chart
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories, size=12)

        # Adjust label positions to avoid overlap
        for label, angle in zip(ax.get_xticklabels(), angles[:-1]):
            if angle >= 0 and angle <  pi/2:
                label.set_horizontalalignment('left')
            elif  angle >= pi/2 and angle < 3*pi/2:
                label.set_horizontalalignment('right')
            else:
                label.set_horizontalalignment('left')

        # Set y-axis (radial) properties dynamically based on number of categories
        # If 4 categories: max 25 points per module
        # If 5 categories: max 20 points per module
        if n_categories == 4:
            max_val = 25
            ax.set_ylim(0, max_val)
            ax.set_yticks([5, 10, 15, 20, 25])
            ax.set_yticklabels(['5', '10', '15', '20', '25'], size=10)
        else:  # 5 categories or default
            max_val = 20
            ax.set_ylim(0, max_val)
            ax.set_yticks([4, 8, 12, 16, 20])
            ax.set_yticklabels(['4', '8', '12', '16', '20'], size=10)

        # Add grid
        ax.grid(True, linestyle='--', alpha=0.7)

        # Title and legend
        plotting.plt.title('Molecular Database Quality Assessment\nRadar Chart',
                 size=16, pad=30)
        # Move legend outside the plot area to avoid overlapping
        plotting.plt.legend(loc='center left', bbox_to_anchor=(1.30, 1.17),
                  frameon=True, fancybox=True)

        plotting.plt.tight_layout()

        if save_path:
            plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        #

