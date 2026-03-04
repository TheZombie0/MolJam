from .._common import *
from .. import plotting


class ActivityConsistencyPlotMixin:
    def plot_activity_kde_distribution(self, figsize=(14, 8), save_path=None):
        """Plot KDE distribution of activity values for all databases"""
        if not self.scoring_results:
            print("No scoring results available.")
            return
            
        fig, ax = plotting.plt.subplots(figsize=figsize)
        
        for db_idx, (db_name, results) in enumerate(self.scoring_results.items()):
            snap = results['snapshot']
            activity_cols = results['activity_cols']

            if activity_cols and activity_cols[0] in snap.df.columns:
                activity_data = snap.df[activity_cols[0]].dropna()
                
                if len(activity_data) > 1:
                    # Calculate KDE
                    from scipy.stats import gaussian_kde
                    kde = gaussian_kde(activity_data)
                    
                    # Create range for plotting
                    x_range = np.linspace(activity_data.min(), activity_data.max(), 200)
                    density = kde(x_range)
                    
                    # Plot
                    ax.plot(x_range, density, 
                           color=self.colors[db_idx % len(self.colors)],
                           linewidth=2, label=f'{db_name} (n={len(activity_data)})',
                           alpha=0.8)
                    
                    # Fill area under curve
                    ax.fill_between(x_range, density, alpha=0.2,
                                   color=self.colors[db_idx % len(self.colors)])
        
        ax.set_xlabel('Activity Value', fontsize=21)
        ax.set_ylabel('Density', fontsize=21)
        ax.set_title('Activity Value Distribution (KDE)', fontsize=23)
        ax.legend(loc='best', fontsize=16)
        ax.grid(True, alpha=0.3)
        
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.quality_dir, 'activity_kde_distribution.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

    def plot_duplicate_activity_consistency(self, figsize=(14, 8), save_path=None):
        """Plot consistency of activity values for duplicate molecules"""
        if not self.scoring_results:
            print("No scoring results available.")
            return
            
        fig, (ax1, ax2) = plotting.plt.subplots(1, 2, figsize=figsize)
        
        db_names = []
        avg_cvs = []
        max_cvs = []
        
        for db_name, results in self.scoring_results.items():
            snap = results['snapshot']
            activity_cols = results['activity_cols']

            if activity_cols and hasattr(snap, 'valid_df') and not snap.valid_df.empty:
                grouped = snap.valid_df.groupby('canonical_smiles')
                cv_values = []
                
                for smiles, group in grouped:
                    if len(group) > 1 and activity_cols[0] in group.columns:
                        values = group[activity_cols[0]].dropna().values
                        if len(values) > 1:
                            mean = np.mean(values)
                            std = np.std(values)
                            cv = std / abs(mean) if mean != 0 else 0
                            cv_values.append(cv)
                
                if cv_values:
                    db_names.append(db_name)
                    avg_cvs.append(np.mean(cv_values))
                    max_cvs.append(np.max(cv_values))
        
        # Plot average CV
        bars1 = ax1.bar(db_names, avg_cvs,
                       color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
        ax1.set_xlabel('Database', fontsize=21)
        ax1.set_ylabel('Average Coefficient of Variation', fontsize=21)
        ax1.set_title('Average Activity Consistency for Duplicates', fontsize=22)
        plotting.plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
        ax1.tick_params(labelsize=18)
        
        for bar, cv in zip(bars1, avg_cvs):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(avg_cvs)*0.01,
                    f'{cv:.3f}', ha='center', va='bottom', fontsize=9)
        
        # Plot max CV
        bars2 = ax2.bar(db_names, max_cvs,
                       color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
        ax2.set_xlabel('Database', fontsize=21)
        ax2.set_ylabel('Maximum Coefficient of Variation', fontsize=21)
        ax2.set_title('Maximum Activity Inconsistency for Duplicates', fontsize=22)
        plotting.plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
        ax2.tick_params(labelsize=18)
        
        for bar, cv in zip(bars2, max_cvs):
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(max_cvs)*0.01,
                    f'{cv:.3f}', ha='center', va='bottom', fontsize=9)
        
        plotting.plt.suptitle('Activity Data Consistency for Duplicate Molecules', fontsize=23)
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.quality_dir, 'duplicate_activity_consistency.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

    def plot_label_conflict_heatmap(self, figsize=(12, 10), save_path=None):
        """Plot heatmap of label conflicts for databases with multiple class columns"""
        if not self.scoring_results:
            print("No scoring results available.")
            return
            
        # Find databases with multiple class columns
        multi_class_dbs = {}
        for db_name, results in self.scoring_results.items():
            class_cols = results['class_cols']
            if class_cols and len(class_cols) >= 1:
                multi_class_dbs[db_name] = class_cols[:10]  # Limit to 10 columns for visibility
        
        if not multi_class_dbs:
            print("No databases with multiple class columns found.")
            return
        
        for db_name, class_cols in multi_class_dbs.items():
            snap = self.scoring_results[db_name]['snapshot']

            # Create conflict matrix
            n_cols = len(class_cols)
            conflict_matrix = np.zeros((n_cols, n_cols))

            if hasattr(snap, 'valid_df') and not snap.valid_df.empty:
                grouped = snap.valid_df.groupby('canonical_smiles')
                
                for smiles, group in grouped:
                    if len(group) > 1:
                        for i, col1 in enumerate(class_cols):
                            for j, col2 in enumerate(class_cols):
                                if i != j and col1 in group.columns and col2 in group.columns:
                                    vals1 = group[col1].dropna().unique()
                                    vals2 = group[col2].dropna().unique()
                                    if len(vals1) > 1 or len(vals2) > 1:
                                        conflict_matrix[i, j] += 1
            
            # Plot heatmap
            fig, ax = plotting.plt.subplots(figsize=figsize)
            
            im = ax.imshow(conflict_matrix, cmap='YlOrRd', aspect='auto')
            
            # Set ticks
            ax.set_xticks(np.arange(n_cols))
            ax.set_yticks(np.arange(n_cols))
            ax.set_xticklabels([col[:20] for col in class_cols], rotation=45, ha='right', fontsize=18)
            ax.set_yticklabels([col[:20] for col in class_cols], fontsize=18)
            
            # Add colorbar
            cbar = plotting.plt.colorbar(im, ax=ax)
            cbar.set_label('Number of Conflicting Molecules', rotation=270, labelpad=20)
            
            # Add text annotations
            for i in range(n_cols):
                for j in range(n_cols):
                    text = ax.text(j, i, int(conflict_matrix[i, j]),
                                  ha="center", va="center", color="black", fontsize=14)
            
            ax.set_title(f'Label Conflict Heatmap - {db_name}', fontsize=23)
            plotting.plt.tight_layout()
            
            save_path_db = save_path or os.path.join(self.quality_dir, f'label_conflict_heatmap_{db_name}.png')
            plotting.plt.savefig(save_path_db, dpi=500, bbox_inches='tight')
            plotting.plt.close()

