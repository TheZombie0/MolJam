from .._common import *
from .. import plotting


class ActivityDistributionPlotMixin:
    def plot_activity_boxplot_comparison(self, figsize=(14, 8), save_path=None):
        """Plot boxplot comparison of activity values"""
        if not self.scoring_results:
            print("No scoring results available.")
            return
        
        # Collect activity data
        activity_data = []
        db_names = []
        
        for db_name, results in self.scoring_results.items():
            scorer = results['scorer']
            activity_cols = results['activity_cols']
            
            if activity_cols and activity_cols[0] in scorer.df.columns:
                data = scorer.df[activity_cols[0]].dropna().values
                if len(data) > 0:
                    activity_data.append(data)
                    db_names.append(db_name)
        
        if not activity_data:
            print("No activity data available.")
            return
        
        fig, ax = plotting.plt.subplots(figsize=figsize)
        
        # Create boxplot
        bp = ax.boxplot(activity_data, patch_artist=True, labels=db_names, showfliers=False)
        
        # Color boxes
        for i, box in enumerate(bp['boxes']):
            box.set_facecolor(self.colors[i % len(self.colors)])
            box.set_alpha(0.7)
        
        ax.set_xlabel('Database', fontsize=21)
        ax.set_ylabel('Activity Value', fontsize=21)
        ax.set_title('Activity Value Distribution Comparison', fontsize=23)
        plotting.plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        ax.tick_params(labelsize=18)
        ax.grid(True, alpha=0.3)
        
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.distribution_dir, 'activity_boxplot_comparison.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

    def plot_activity_raincloud(self, activity_col=None, figsize=(14, 8), save_path=None, show_points=False):
        """Plot raincloud plot of activity values."""
        if not self.scoring_results:
            print("No scoring results available.")
            return
        
        # Prepare data for raincloud plot
        all_data = []
        db_names = []
        
        for db_name, results in self.scoring_results.items():
            scorer = results['scorer']
            activity_cols = results['activity_cols']
            
            # Use specified activity column or first activity column
            col = activity_col
            if col is None and activity_cols:
                col = activity_cols[0]
            
            if col and col in scorer.df.columns:
                data = scorer.df[col].dropna().values
                if len(data) > 0:
                    # Add to combined data
                    for val in data:
                        all_data.append({'Database': db_name, 'Activity': val})
                    db_names.append(db_name)
        
        if not all_data:
            print("No activity data available.")
            return
        
        # Convert to DataFrame
        df = pd.DataFrame(all_data)
        
        # Create figure
        fig, ax = plotting.plt.subplots(figsize=figsize)
        
        # Use seaborn for violin plot (density)
        plotting.sns.violinplot(x='Database', y='Activity', data=df,
                      palette=[self.colors[i % len(self.colors)] for i in range(len(db_names))],
                      inner=None, cut=0, ax=ax, alpha=0.5)
        
        # Add boxplot on top
        plotting.sns.boxplot(x='Database', y='Activity', data=df,
                    width=0.2, ax=ax, color='white',
                    boxprops=dict(alpha=0.7), showfliers=False)
        
        if show_points:
            # Add stripplot for individual points (subsample if too many)
            if len(df) > 5000:
                df_sample = df.sample(5000, random_state=42)
            else:
                df_sample = df

            plotting.sns.stripplot(x='Database', y='Activity', data=df_sample,
                          color='black', size=2, alpha=0.3, ax=ax, jitter=True)
        
        ax.set_xlabel('Database', fontsize=21)
        ax.set_ylabel('Activity Value', fontsize=21)
        ax.set_title('Activity Value Raincloud Plot', fontsize=23)
        plotting.plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        ax.tick_params(labelsize=18)
        ax.grid(True, alpha=0.3)
        
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.distribution_dir, 'activity_raincloud.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()
