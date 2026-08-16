from .._common import *
from .. import plotting


class RadarPlotMixin:
    radar_category_labels = {
        "Structural Integrity": "Structural Integrity",
        "Data Quality": "Label Consistency",
        "Experimental Information Quality": "Experimental Information\nQuality",
        "Chemical Space Coverage": "Chemical Space Coverage",
        "Data Distribution": "Size and Distribution",
    }

    def _resolve_radar_database(self, db_selector):
        if isinstance(db_selector, str):
            if db_selector not in self.scoring_results:
                raise KeyError(f"Database '{db_selector}' not found in scoring results.")
            return db_selector, self.scoring_results[db_selector]

        if isinstance(db_selector, int):
            db_items = list(self.scoring_results.items())
            if db_selector < 1 or db_selector > len(db_items):
                available_names = [name for name, _ in db_items]
                raise IndexError(
                    f"Database index {db_selector} is out of range. Expected 1 to {len(db_items)}. "
                    f"Currently available databases: {available_names}. "
                    "Note: self.scoring_results is keyed by database name. If you added multiple "
                    "datasets with the same name, later entries overwrite earlier ones. Use unique "
                    "names such as 'BACE_raw'/'BACE_cleaned' and 'HIV_raw'/'HIV_cleaned'."
                )
            return db_items[db_selector - 1]

        raise TypeError("Database selector must be a database name or a 1-based integer index.")

    def _get_radar_categories(self, db_entries):
        available_categories = []
        for category in self.categories_labels.keys():
            if any(category in results['scorer'].scores for _, results in db_entries):
                available_categories.append(category)
        return available_categories

    def _get_radar_scores(self, scorer, categories):
        return [scorer.scores.get(category, {}).get("Normalized Total", 0) for category in categories]

    def _get_radar_category_labels(self, categories):
        return [self.radar_category_labels.get(category, category) for category in categories]

    def _get_radar_entry_colors(self, db_names):
        db_order = list(self.scoring_results.keys())
        return [
            self.colors[db_order.index(db_name) % len(self.colors)]
            for db_name in db_names
        ]

    def _normalize_radar_comparison_selectors(self, db_selectors):
        if len(db_selectors) == 1 and isinstance(db_selectors[0], (list, tuple)):
            db_selectors = tuple(db_selectors[0])

        if len(db_selectors) < 2:
            raise ValueError("Please select at least two databases for radar comparison.")

        resolved_entries = [self._resolve_radar_database(selector) for selector in db_selectors]
        resolved_names = [db_name for db_name, _ in resolved_entries]

        if len(set(resolved_names)) != len(resolved_names):
            raise ValueError("Please select different databases for radar comparison.")

        return resolved_entries

    def _style_radar_axes(self, ax, angles, categories, title, legend_bbox):
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(self._get_radar_category_labels(categories), size=16)

        for label, angle in zip(ax.get_xticklabels(), angles[:-1]):
            if angle >= 0 and angle < pi / 2:
                label.set_horizontalalignment('left')
            elif angle >= pi / 2 and angle < 3 * pi / 2:
                label.set_horizontalalignment('right')
            else:
                label.set_horizontalalignment('left')

        if len(categories) == 4:
            ax.set_ylim(0, 25)
            ax.set_yticks([5, 10, 15, 20, 25])
            ax.set_yticklabels(['5', '10', '15', '20', '25'], size=14)
        else:
            ax.set_ylim(0, 20)
            ax.set_yticks([4, 8, 12, 16, 20])
            ax.set_yticklabels(['4', '8', '12', '16', '20'], size=14)

        ax.grid(True, linestyle='--', alpha=0.7)
        plotting.plt.title(title, size=20, pad=30)
        plotting.plt.legend(
            loc='best',
            frameon=True,
            fancybox=True,
            fontsize=16,
        )

    def _plot_radar_entries(
        self,
        db_entries,
        figsize=(11, 11),
        save_path=None,
        include_average=False,
        colors=None,
        fill_alpha=0.15,
        line_alpha=1.0,
        title='Molecular Database Quality Assessment\nRadar Chart',
        legend_bbox=(1.30, 1.17),
    ):
        if not db_entries:
            print("No scoring results available. Please add databases first.")
            return

        categories = self._get_radar_categories(db_entries)
        if not categories:
            print("No radar categories available for the selected databases.")
            return

        fig, ax = plotting.plt.subplots(figsize=figsize, subplot_kw=dict(projection='polar'))

        n_categories = len(categories)
        angles = [n / n_categories * 2 * pi for n in range(n_categories)]
        angles += angles[:1]

        if colors is None:
            colors = [self.colors[idx % len(self.colors)] for idx in range(len(db_entries))]

        category_score_rows = []
        for db_idx, (db_name, results) in enumerate(db_entries):
            scorer = results['scorer']
            category_scores = self._get_radar_scores(scorer, categories)
            category_score_rows.append(category_scores.copy())

            category_scores += category_scores[:1]
            color = colors[db_idx % len(colors)]
            total_score = scorer.scores.get(
                "Final Adjusted Score",
                scorer.scores.get("Normalized Score", 0),
            )

            ax.plot(
                angles,
                category_scores,
                'o-',
                linewidth=2,
                label=f'{db_name} ({total_score:.1f}/100)',
                color=color,
                alpha=line_alpha,
            )
            ax.fill(angles, category_scores, color=color, alpha=fill_alpha)

        if include_average and category_score_rows:
            avg_scores = np.mean(category_score_rows, axis=0).tolist()
            avg_scores += avg_scores[:1]
            ax.plot(
                angles,
                avg_scores,
                linestyle='--',
                linewidth=2.4,
                color='#444444',
                label='Average',
            )

        self._style_radar_axes(ax, angles, categories, title, legend_bbox)
        plotting.plt.tight_layout()

        if save_path:
            plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')

        return fig, ax

    def plot_radar(self, figsize=(11, 11), save_path=None):
        """Plot a radar chart showing normalized scores for each category with an average line."""
        self._plot_radar_entries(
            list(self.scoring_results.items()),
            figsize=figsize,
            save_path=save_path,
            include_average=True,
            colors=self._get_radar_entry_colors(list(self.scoring_results.keys())),
        )

    def plot_radar_no_avg(self, figsize=(11, 11), save_path=None):
        """Plot a radar chart showing normalized scores for each category without the average line."""
        self._plot_radar_entries(
            list(self.scoring_results.items()),
            figsize=figsize,
            save_path=save_path,
            include_average=False,
            colors=self._get_radar_entry_colors(list(self.scoring_results.keys())),
        )

    def plot_radar_comparison(self, *db_selectors, figsize=(11, 11), save_path=None):
        """Plot a radar comparison for two or more databases.

        Selectors can be database names or 1-based integer indices.
        Supports both plot_radar_comparison(1, 2, 3) and
        plot_radar_comparison([1, 2, 3]).
        """
        if not self.scoring_results:
            print("No scoring results available. Please add databases first.")
            return

        comparison_entries = self._normalize_radar_comparison_selectors(db_selectors)
        comparison_names = [db_name for db_name, _ in comparison_entries]
        comparison_colors = self._get_radar_entry_colors(comparison_names)
        title = 'Molecular Database Quality Assessment'
        if len(comparison_names) == 2:
            title = f'{title}\nRadar Comparison: {comparison_names[0]} vs {comparison_names[1]}'
        self._plot_radar_entries(
            comparison_entries,
            figsize=figsize,
            save_path=save_path,
            include_average=False,
            colors=comparison_colors,
            title=title,
            legend_bbox=(1.24, 1.12),
        )
