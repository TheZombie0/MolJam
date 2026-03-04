from .._common import *
from .. import plotting


class ClassDistributionPlotMixin:
    def plot_class_distribution(self, db_name, class_col=None, figsize=(10, 8), save_path=None):
        """Plot class distribution for a specific database"""
        if db_name not in self.scoring_results:
            print(f"Database {db_name} not found in scoring results.")
            return
        
        snap = self.scoring_results[db_name]['snapshot']
        class_cols = self.scoring_results[db_name]['class_cols']

        # Use specified class column or first class column
        col = class_col
        if col is None and class_cols:
            col = class_cols[0]

        if not col or col not in snap.df.columns:
            print(f"Class column '{col}' not found in {db_name}.")
            return
        
        # Get class distribution
        class_counts = snap.df[col].value_counts()
        
        fig, ax = plotting.plt.subplots(figsize=figsize)
        
        # Create pie chart for binary/multi-class distribution
        if len(class_counts) <= 10:  # Reasonable number of classes for pie chart
            colors = [self.colors[i % len(self.colors)] for i in range(len(class_counts))]
            wedges, texts, autotexts = ax.pie(class_counts.values, labels=class_counts.index,
                                             autopct='%1.1f%%', startangle=90,
                                             colors=colors)
            
            # Improve text formatting
            for text in texts:
                text.set_fontsize(18)
            for autotext in autotexts:
                autotext.set_fontsize(18)
                autotext.set_color('white')
                
            ax.set_title(f'Class Distribution - {db_name}\n({col})', fontsize=23)
            
        else:  # Too many classes, use bar chart
            # Limit to top 20 classes
            top_classes = class_counts.head(20)
            
            bars = ax.bar(range(len(top_classes)), top_classes.values,
                         color=self.colors[0], alpha=0.7)
            
            ax.set_xlabel('Class', fontsize=21)
            ax.set_ylabel('Count', fontsize=21)
            ax.set_title(f'Top 20 Class Distribution - {db_name}\n({col})', fontsize=23)
            ax.set_xticks(range(len(top_classes)))
            ax.set_xticklabels([str(c)[:10] for c in top_classes.index], rotation=45, ha='right')
            ax.tick_params(labelsize=18)
        
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.scoring_results[db_name]['output_dir'],
                                   f'class_distribution_{col}.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

