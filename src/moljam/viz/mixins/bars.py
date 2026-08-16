from .._common import *
from .. import plotting


class BarsPlotMixin:
    def plot_bars_single_canvas(self, show_values=True, figsize=(20, 8), save_path=None):
        """Plot all scores in a single canvas with grouped bars"""
        if not self.scoring_results:
            print("No scoring results available. Please add databases first.")
            return

        fig, ax = plotting.plt.subplots(figsize=figsize)

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

        # Prepare data
        db_names = list(self.scoring_results.keys())
        n_dbs = len(db_names)

        # Fixed parameters
        bar_width = 0.8 / max(n_dbs, 1)  # Ensure consistent bar width
        group_width = bar_width * n_dbs  # Width of each metric group
        spacing_between_metrics = 0.2
        spacing_between_categories = 0.8

        # Calculate positions
        all_metrics_labels = []
        all_metrics_raw = []
        category_positions = []
        metric_centers = []

        current_pos = 0

        for cat_idx, (category, metrics) in enumerate(available_categories.items()):
            if cat_idx > 0:
                current_pos += spacing_between_categories
            
            category_start = current_pos
            
            for i, metric in enumerate(metrics):
                all_metrics_labels.append(self.categories_labels[category][i])
                all_metrics_raw.append(metric)
                
                # Calculate center position for this metric group
                metric_center = current_pos + group_width / 2
                metric_centers.append(metric_center)
                
                current_pos += group_width + spacing_between_metrics
            
            current_pos -= spacing_between_metrics  # Remove last spacing
            category_positions.append((category_start, current_pos))
            
        # Plot bars for each database
        for db_idx, (db_name, results) in enumerate(self.scoring_results.items()):
            scorer = results['scorer']
            metric_scores = []
            metric_positions = []
            
            # Collect scores for all metrics
            for metric in all_metrics_raw:
                for category, metrics in available_categories.items():
                    if metric in metrics:
                        score = scorer.scores.get(category, {}).get(metric, 0)
                        # Handle None values by treating them as 0
                        if score is None:
                            score = 0
                        metric_scores.append(score)
                        break
            
            # Calculate bar positions
            for center in metric_centers:
                pos = center - group_width/2 + db_idx * bar_width
                metric_positions.append(pos)
            
            # Plot bars for this database
            bars = ax.bar(metric_positions, metric_scores, bar_width, 
                          label=db_name, color=self.colors[db_idx % len(self.colors)])
            
            # Add value labels if requested
            if show_values:
                for bar in bars:
                    height = bar.get_height()
                    if height > 0:  # Only show label if there's a value
                        ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                               f'{height:.1f}', ha='center', va='bottom', fontsize=8)
        
        # Set x-axis
        ax.set_xticks(metric_centers)
        ax.set_xticklabels(all_metrics_labels, rotation=45, ha='right')
        
        # Add category labels and separators
        for i, (category, (start, end)) in enumerate(zip(available_categories.keys(), category_positions)):
            # Add vertical lines between categories
            if i > 0:
                ax.axvline(x=start - spacing_between_categories/2, color='gray', 
                          linestyle='--', alpha=0.5)
            
            # Add category label
            center = (start + end) / 2
            ax.text(center, -0.25, category, ha='right', va='top', 
                   fontsize=14, transform=ax.get_xaxis_transform())
        
        # Customize plot
        ax.set_ylabel('Score', fontsize=21)
        ax.set_title('Molecular Database Quality Assessment Scores', fontsize=23, pad=20)
        ax.tick_params(labelsize=14)
        ax.set_ylim(0, 11)
        ax.legend(loc='upper left',bbox_to_anchor=(1.05, 0.8), frameon=True, fancybox=True, fontsize=14)
        
        # Remove grid
        ax.grid(False)
        
        # Adjust x-axis limits to remove empty space
        if metric_centers:
            ax.set_xlim(-0.5, max(metric_centers) + group_width/2 + 0.5)
        
        # Adjust layout
        plotting.plt.tight_layout()
        plotting.plt.subplots_adjust(bottom=0.15)
        
        if save_path:
            plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        #
    def plot_bars_multi_canvas(self, show_values=True, figsize=(20, 12), save_path=None):
        """Plot scores in multiple canvases (one for each category)"""
        if not self.scoring_results:
            print("No scoring results available. Please add databases first.")
            return

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

        n_categories = len(available_categories)

        # Create dynamic subplot grid based on number of categories
        if n_categories == 4:
            # 2x2 grid for 4 categories
            fig, axes = plotting.plt.subplots(2, 2, figsize=figsize)
            axes = axes.flatten()
        else:  # 5 categories or default
            # 2x3 grid for 5 categories
            fig, axes = plotting.plt.subplots(2, 3, figsize=figsize)
            axes = axes.flatten()

        db_names = list(self.scoring_results.keys())
        n_dbs = len(db_names)
        bar_width = 0.6 / max(n_dbs, 1)  # Fixed bar width for all subplots

        for cat_idx, (category, metrics) in enumerate(available_categories.items()):
            ax = axes[cat_idx]

            # Prepare data for this category
            metric_positions = np.arange(len(metrics))

            for db_idx, (db_name, results) in enumerate(self.scoring_results.items()):
                scorer = results['scorer']
                scores = []
                for metric in metrics:
                    score = scorer.scores.get(category, {}).get(metric, 0)
                    # Handle None values by treating them as 0
                    if score is None:
                        score = 0
                    scores.append(score)
                metrics_label = [self.categories_labels[category][i] for i in range(len(metrics))]

                # Calculate positions for bars
                positions = metric_positions + db_idx * bar_width - (n_dbs - 1) * bar_width / 2

                bars = ax.bar(positions, scores, bar_width,
                             label=db_name, color=self.colors[db_idx % len(self.colors)])

                # Add value labels if requested
                if show_values:
                    for bar in bars:
                        height = bar.get_height()
                        if height > 0:
                            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                                   f'{height:.1f}', ha='center', va='bottom', fontsize=9)

            # Customize subplot
            ax.set_title(category, fontsize=22, pad=10)
            ax.set_xticks(metric_positions)
            ax.set_xticklabels(metrics_label, rotation=45, ha='right')
            ax.tick_params(labelsize=14)
            ax.set_ylabel('Score', fontsize=21)
            ax.set_ylim(0, 11)
            ax.grid(False)

        # Handle legend placement based on number of categories
        if n_categories == 4:
            # For 4 categories, use the entire figure for legend
            # Create legend handles
            handles = []
            for db_idx, db_name in enumerate(db_names):
                handles.append(plotting.mpatches.Patch(color=self.colors[db_idx % len(self.colors)], label=db_name))

            # Place legend below the subplots
            fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.05),
                      ncol=min(len(db_names), 4), frameon=True, fancybox=True, fontsize=14)
        else:  # 5 categories or more
            # Use the last subplot for legend as before
            ax_legend = axes[n_categories]
            ax_legend.axis('off')  # Turn off the axis

            # Create legend handles
            handles = []
            for db_idx, db_name in enumerate(db_names):
                handles.append(plotting.mpatches.Patch(color=self.colors[db_idx % len(self.colors)], label=db_name))

            # Place legend in the center of the last subplot
            ax_legend.legend(handles=handles, loc='center', frameon=True, fancybox=True, fontsize=14)

        # Overall title
        fig.suptitle('Molecular Database Quality Assessment Scores by Category',
                    fontsize=23)

        plotting.plt.tight_layout()

        if save_path:
            plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        #
