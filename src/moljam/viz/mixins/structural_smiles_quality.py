from .._common import *
from .. import plotting
from ..structural_breakdowns import (
    INVALID_SMILES_NORMALIZE_FILENAMES,
    build_invalid_smiles_plot_table,
    build_invalid_smiles_source_level_summaries,
    empty_invalid_smiles_details_frame,
    inspect_smiles_rows,
    inspections_to_frame,
    normalize_smiles_rows,
    render_invalid_smiles_stacked_horizontal_figure,
)


class SmilesQualityPlotMixin:
    def _build_invalid_smiles_error_breakdown_tables(self):
        source_payloads = {}
        detail_frames = []

        for db_name, results in self.scoring_results.items():
            scorer = results["scorer"]
            if "validate_smiles" not in scorer.completed_checks:
                scorer.validate_smiles()

            filtered_rows, _, _ = normalize_smiles_rows(scorer.df, scorer.smiles_col)
            invalid_indices = set(getattr(scorer, "invalid_indices", []))
            invalid_rows = filtered_rows.loc[filtered_rows["row_index"].isin(invalid_indices), ["row_index", "smiles"]].copy()

            if invalid_rows.empty:
                details_df = empty_invalid_smiles_details_frame()
            else:
                details_df = inspections_to_frame(inspect_smiles_rows(invalid_rows))
                details_df["source"] = db_name

            source_payloads[db_name] = {
                "total_count": int(len(filtered_rows)),
                "details_df": details_df,
            }
            detail_frames.append(details_df)

        source_level_df, per_source_category_df = build_invalid_smiles_source_level_summaries(
            source_payloads
        )
        combined_details_df = (
            pd.concat(detail_frames, ignore_index=True)
            if detail_frames and any(not frame.empty for frame in detail_frames)
            else empty_invalid_smiles_details_frame()
        )
        return combined_details_df, source_level_df, per_source_category_df

    def plot_invalid_smiles_error_categories(
        self,
        normalize_by="invalid",
        save_path=None,
        save_svg_path=None,
        title="Invalid SMILES error frequencies",
        panel_label=None,
    ):
        """Plot RDKit-derived invalid SMILES error categories across databases."""
        if not self.scoring_results:
            print("No scoring results available.")
            return

        _, source_level_df, per_source_category_df = self._build_invalid_smiles_error_breakdown_tables()
        plot_df = build_invalid_smiles_plot_table(
            source_level_df=source_level_df,
            per_source_category_df=per_source_category_df,
            normalize_by=normalize_by,
        )

        if save_path is None:
            save_path = os.path.join(
                self.comparison_dir,
                f"invalid_smiles_error_frequencies_{INVALID_SMILES_NORMALIZE_FILENAMES[normalize_by]}.png",
            )

        save_dir = os.path.dirname(save_path)
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
        if save_svg_path:
            save_svg_dir = os.path.dirname(save_svg_path)
            if save_svg_dir:
                os.makedirs(save_svg_dir, exist_ok=True)

        render_invalid_smiles_stacked_horizontal_figure(
            figure_png_path=save_path,
            figure_svg_path=save_svg_path,
            plot_df=plot_df,
            normalize_by=normalize_by,
            title=title,
            panel_label=panel_label,
        )
        print(f"Saved invalid SMILES error categories plot to: {save_path}")

    def plot_invalid_smiles_comparison(self, mode='both', figsize=(14, 8), save_path=None):
        """
        Plot invalid SMILES count and ratio comparison (updated both mode)
        """
        if not self.scoring_results:
            print("No scoring results available.")
            return
            
        db_names = []
        invalid_counts = []
        invalid_ratios = []
        
        for db_name, results in self.scoring_results.items():
            scorer = results['scorer']
            
            if hasattr(scorer, 'invalid_indices'):
                invalid_count = len(scorer.invalid_indices)
                invalid_ratio = scorer.invalid_rate
            else:
                invalid_count = 0
                invalid_ratio = 0
                
            db_names.append(db_name)
            invalid_counts.append(invalid_count)
            invalid_ratios.append(invalid_ratio)
        
        if mode == 'both':
            # Create two subplots side by side
            fig, (ax1, ax2) = plotting.plt.subplots(1, 2, figsize=figsize)
            
            # Left plot: Count
            bars1 = ax1.bar(db_names, invalid_counts,
                           color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax1.set_xlabel('Database', fontsize=21)
            ax1.set_ylabel('Invalid SMILES Count', fontsize=21)
            ax1.set_title('Invalid SMILES Count', fontsize=22)
            ax1.tick_params(axis='x', rotation=45, labelsize=18)
            plotting.plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
            
            # Add value labels
            for bar, count in zip(bars1, invalid_counts):
                if int(count)>0:
                    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(invalid_counts)*0.01,
                            str(count), ha='center', va='bottom', fontsize=9)
            
            # Right plot: Ratio
            bars2 = ax2.bar(db_names, invalid_ratios,
                           color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax2.set_xlabel('Database', fontsize=21)
            ax2.set_ylabel('Invalid SMILES Ratio (%)', fontsize=21)
            ax2.set_title('Invalid SMILES Ratio', fontsize=22)
            plotting.plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
            ax2.tick_params(axis='x', rotation=45, labelsize=18)
            
            # Add value labels
            for bar, ratio in zip(bars2, invalid_ratios):
                if float(ratio)>0:
                    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(invalid_ratios)*0.01,
                            f'{ratio:.1f}%', ha='center', va='bottom', fontsize=9)
            
            plotting.plt.suptitle('Invalid SMILES Comparison Across Databases', fontsize=23)
            
        elif mode == 'count':
            fig, ax = plotting.plt.subplots(figsize=figsize)
            bars = ax.bar(db_names, invalid_counts,
                          color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax.set_xlabel('Database', fontsize=21)
            ax.set_ylabel('Invalid SMILES Count', fontsize=21)
            ax.set_title('Invalid SMILES Count Comparison', fontsize=22)
            
            for bar, count in zip(bars, invalid_counts):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(invalid_counts)*0.01,
                       str(count), ha='center', va='bottom', fontsize=9)
            
            plotting.plt.xticks(rotation=45, ha='right')
            ax.tick_params(labelsize=18)
            
        else:  # mode == 'ratio'
            fig, ax = plotting.plt.subplots(figsize=figsize)
            bars = ax.bar(db_names, invalid_ratios,
                          color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax.set_xlabel('Database', fontsize=21)
            ax.set_ylabel('Invalid SMILES Ratio (%)', fontsize=21)
            ax.set_title('Invalid SMILES Ratio Comparison', fontsize=22)
            
            for bar, ratio in zip(bars, invalid_ratios):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(invalid_ratios)*0.01,
                       f'{ratio:.1f}%', ha='center', va='bottom', fontsize=9)
            
            plotting.plt.xticks(rotation=45, ha='right')
            ax.tick_params(labelsize=18)
        
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.comparison_dir, f'invalid_smiles_{mode}.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

    def plot_non_standardized_smiles_comparison(self, mode='both', figsize=(14, 8), save_path=None):
        """
        Plot non-standardized SMILES count and ratio comparison (updated both mode)
        """
        if not self.scoring_results:
            print("No scoring results available.")
            return
            
        db_names = []
        non_canonical_counts = []
        non_canonical_ratios = []
        
        for db_name, results in self.scoring_results.items():
            scorer = results['scorer']
            
            if hasattr(scorer, 'non_canonical_indices'):
                non_canonical_count = len(scorer.non_canonical_indices)
                non_canonical_ratio = scorer.non_canonical_rate
            else:
                non_canonical_count = 0
                non_canonical_ratio = 0
                
            db_names.append(db_name)
            non_canonical_counts.append(non_canonical_count)
            non_canonical_ratios.append(non_canonical_ratio)
        
        if mode == 'both':
            # Create two subplots side by side
            fig, (ax1, ax2) = plotting.plt.subplots(1, 2, figsize=figsize)
            
            # Left plot: Count
            bars1 = ax1.bar(db_names, non_canonical_counts,
                           color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax1.set_xlabel('Database', fontsize=21)
            ax1.set_ylabel('Non-standardized SMILES Count', fontsize=21)
            ax1.set_title('Non-standardized SMILES Count', fontsize=22)
            plotting.plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
            ax1.tick_params(labelsize=18)
            
            # Add value labels
            for bar, count in zip(bars1, non_canonical_counts):
                ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(non_canonical_counts)*0.01,
                        str(count), ha='center', va='bottom', fontsize=9)
            
            # Right plot: Ratio
            bars2 = ax2.bar(db_names, non_canonical_ratios,
                           color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax2.set_xlabel('Database', fontsize=21)
            ax2.set_ylabel('Non-standardized SMILES Ratio (%)', fontsize=21)
            ax2.set_title('Non-standardized SMILES Ratio', fontsize=22)
            plotting.plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
            ax2.tick_params(labelsize=18)
            
            # Add value labels
            for bar, ratio in zip(bars2, non_canonical_ratios):
                ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(non_canonical_ratios)*0.01,
                        f'{ratio:.1f}%', ha='center', va='bottom', fontsize=9)
            
            plotting.plt.suptitle('Non-standardized SMILES Comparison Across Databases', fontsize=23)
            
        elif mode == 'count':
            fig, ax = plotting.plt.subplots(figsize=figsize)
            bars = ax.bar(db_names, non_canonical_counts,
                          color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax.set_xlabel('Database', fontsize=21)
            ax.set_ylabel('Non-standardized SMILES Count', fontsize=21)
            ax.set_title('Non-standardized SMILES Count Comparison', fontsize=22)
            
            for bar, count in zip(bars, non_canonical_counts):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(non_canonical_counts)*0.01,
                       str(count), ha='center', va='bottom', fontsize=9)
            
            plotting.plt.xticks(rotation=45, ha='right')
            ax.tick_params(labelsize=18)
            
        else:  # mode == 'ratio'
            fig, ax = plotting.plt.subplots(figsize=figsize)
            bars = ax.bar(db_names, non_canonical_ratios,
                          color=[self.colors[i % len(self.colors)] for i in range(len(db_names))])
            
            ax.set_xlabel('Database', fontsize=21)
            ax.set_ylabel('Non-standardized SMILES Ratio (%)', fontsize=21)
            ax.set_title('Non-standardized SMILES Ratio Comparison', fontsize=22)
            
            for bar, ratio in zip(bars, non_canonical_ratios):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(non_canonical_ratios)*0.01,
                       f'{ratio:.1f}%', ha='center', va='bottom', fontsize=9)
            
            plotting.plt.xticks(rotation=45, ha='right')
            ax.tick_params(labelsize=18)
        
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.comparison_dir, f'non_standardized_smiles_{mode}.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()
