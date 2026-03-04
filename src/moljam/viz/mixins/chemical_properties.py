from .._common import *
from .. import plotting
from ..chem import calculate_properties_parallel, calculate_qed_parallel


class MolecularPropertiesPlotMixin:
    def _compute_molecular_properties(self, db_name):
        """Compute and cache molecular properties for a database"""
        if db_name in self.property_cache:
            return self.property_cache[db_name]

        snap = self.scoring_results[db_name]['snapshot']

        if not hasattr(snap, 'valid_df') or snap.valid_df.empty:
            return None

        # Use parallel processing if enabled
        canonical_smiles = snap.valid_df['canonical_smiles'].tolist()

        if self.use_parallel:
            n_workers = min(cpu_count(), 100)
            with Pool(n_workers) as pool:
                properties_list = pool.map(calculate_properties_parallel, canonical_smiles)
        else:
            # Serial processing
            properties_list = [calculate_properties_parallel(smiles) for smiles in canonical_smiles]

        # Filter out None values and create DataFrame
        valid_properties = [p for p in properties_list if p is not None]
        if valid_properties:
            properties_df = pd.DataFrame(valid_properties)
            self.property_cache[db_name] = properties_df
            return properties_df

        return None

    def plot_molecular_properties_distribution(self, figsize=(20, 13), save_path=None):
        """Plot distribution of molecular properties"""
        if not self.scoring_results:
            print("No scoring results available.")
            return
        
        properties = ['MW', 'LogP', 'HBA', 'HBD', 'TPSA', 'RotBonds']
        property_names = ['Molecular Weight', 'LogP', 'H-Bond Acceptors', 
                         'H-Bond Donors', 'TPSA', 'Rotatable Bonds']
        
        fig, axes = plotting.plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()
        
        for prop_idx, (prop, prop_name) in enumerate(zip(properties, property_names)):
            ax = axes[prop_idx]
            
            for db_idx, db_name in enumerate(self.scoring_results.keys()):
                properties_df = self._compute_molecular_properties(db_name)
                
                if properties_df is not None and prop in properties_df.columns:
                    data = properties_df[prop].dropna()
                    
                    if len(data) > 0:
                        # Create violin plot
                        parts = ax.violinplot([data], positions=[db_idx], widths=0.7,
                                             showmeans=True, showmedians=True)
                        
                        # Color the violin
                        for pc in parts['bodies']:
                            pc.set_facecolor(self.colors[db_idx % len(self.colors)])
                            pc.set_alpha(0.6)
            
            # ax.set_xlabel('Database', fontsize=10)
            ax.set_ylabel(prop_name, fontsize=21)
            ax.set_title(f'{prop_name} Distribution', fontsize=22)
            ax.set_xticks(range(len(self.scoring_results)))
            ax.set_xticklabels(list(self.scoring_results.keys()), rotation=45, ha='right')
            ax.tick_params(labelsize=18)
        
        plotting.plt.suptitle('Molecular Properties Distribution', fontsize=23)
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.coverage_dir, 'molecular_properties_distribution.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

    def plot_lipinski_violations(self, figsize=(16, 8), save_path=None):
        """Plot Lipinski rule violations"""
        if not self.scoring_results:
            print("No scoring results available.")
            return
        
        fig, (ax1, ax2) = plotting.plt.subplots(1, 2, figsize=figsize)
        
        db_names = []
        violation_counts = {0: [], 1: [], 2: [], 3: [], 4: []}
        
        for db_name in self.scoring_results.keys():
            properties_df = self._compute_molecular_properties(db_name)
            
            if properties_df is not None and 'Lipinski' in properties_df.columns:
                db_names.append(db_name)
                
                # Count violations (4 - number of rules passed)
                violations = 4 - properties_df['Lipinski']
                
                for i in range(5):
                    count = (violations == i).sum()
                    violation_counts[i].append(count)
        
        if not db_names:
            print("No Lipinski data available.")
            return
        
        # Stacked bar chart for counts
        x = np.arange(len(db_names))
        width = 0.6
        
        bottoms = np.zeros(len(db_names))
        colors_violations = ['green', 'yellow', 'orange', 'red', 'darkred']
        
        for i in range(5):
            bars = ax1.bar(x, violation_counts[i], width, bottom=bottoms,
                          color=colors_violations[i], label=f'{i} violations',
                          alpha=0.8)
            bottoms += violation_counts[i]
        
        ax1.set_xlabel('Database', fontsize=21)
        ax1.set_ylabel('Number of Molecules', fontsize=21)
        ax1.set_title('Lipinski Rule Violations (Count)', fontsize=22)
        ax1.set_xticks(x)
        ax1.set_xticklabels(db_names, rotation=45, ha='right')
        ax1.tick_params(labelsize=18)
        ax1.legend(fontsize=16)
        
        # Percentage plot
        totals = bottoms
        percentages = {i: [] for i in range(5)}
        
        for j, total in enumerate(totals):
            if total > 0:
                for i in range(5):
                    percentages[i].append(violation_counts[i][j] / total * 100)
            else:
                for i in range(5):
                    percentages[i].append(0)
        
        bottoms_pct = np.zeros(len(db_names))
        
        for i in range(5):
            bars = ax2.bar(x, percentages[i], width, bottom=bottoms_pct,
                          color=colors_violations[i], label=f'{i} violations',
                          alpha=0.8)
            bottoms_pct += percentages[i]
        
        ax2.set_xlabel('Database', fontsize=21)
        ax2.set_ylabel('Percentage (%)', fontsize=21)
        ax2.set_title('Lipinski Rule Violations (Percentage)', fontsize=22)
        ax2.set_xticks(x)
        ax2.set_xticklabels(db_names, rotation=45, ha='right')
        ax2.tick_params(labelsize=18)
        ax2.legend(fontsize=16)
        
        plotting.plt.suptitle('Lipinski Rule of Five Analysis', fontsize=23)
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.coverage_dir, 'lipinski_violations.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

    def plot_qed_distribution(self, figsize=(12, 6), save_path=None):
        """
        绘制多个数据库的QED分数分布图（只显示拟合曲线）
        使用并行处理加速QED计算
        """
        if not self.scoring_results:
            print("No scoring results available.")
            return

        fig, ax = plotting.plt.subplots(figsize=figsize)

        for db_idx, (db_name, results) in enumerate(self.scoring_results.items()):
            snap = results['snapshot']

            # 获取QED数据
            if hasattr(snap, 'valid_df') and not snap.valid_df.empty:
                print(f"  Calculating QED for {db_name} ({len(snap.valid_df)} molecules)...")

                # 使用并行处理计算QED值
                canonical_smiles = snap.valid_df['canonical_smiles'].tolist()

                if self.use_parallel:
                    n_workers = min(cpu_count(), 100)
                    with Pool(n_workers) as pool:
                        qed_results = pool.map(calculate_qed_parallel, canonical_smiles)
                else:
                    # Serial processing
                    qed_results = [calculate_qed_parallel(smiles) for smiles in canonical_smiles]

                # 过滤掉None值
                qed_values = [qed for qed in qed_results if qed is not None]

                if qed_values:
                    print(f"    Successfully calculated QED for {len(qed_values)} molecules")

                    # 使用核密度估计拟合分布
                    from scipy.stats import gaussian_kde

                    kde = gaussian_kde(qed_values)
                    x_range = np.linspace(0, 1, 200)
                    density = kde(x_range)

                    # 绘制分布曲线
                    ax.plot(x_range, density,
                           color=self.colors[db_idx % len(self.colors)],
                           linewidth=2, label=f'{db_name} (n={len(qed_values)})')

        ax.set_xlabel('QED Score', fontsize=12)
        ax.set_ylabel('Density', fontsize=12)
        ax.set_title('QED Score Distribution Comparison', fontsize=16, fontweight='bold')
        ax.set_xlim(0, 1)
        ax.legend()
        ax.grid(True, alpha=0.3)

        plotting.plt.tight_layout()

        if save_path is None:
            save_path = os.path.join(self.comparison_dir, 'qed_distribution.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

