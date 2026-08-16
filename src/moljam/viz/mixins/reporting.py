from .._common import *
from .. import plotting
import tempfile


class ReportingMixin:
    def _append_saved_plot_to_pdf(self, pdf, title, plot_func, save_path):
        plotting.plt.close('all')
        plot_func(save_path=save_path)

        if not os.path.exists(save_path):
            print(f"    Skipped {title}: plot was not generated")
            return

        image = plotting.plt.imread(save_path)
        fig, ax = plotting.plt.subplots(figsize=(11, 8.5))
        ax.imshow(image)
        ax.set_title(title, fontsize=18, pad=12)
        ax.axis('off')
        pdf.savefig(fig, bbox_inches='tight')
        plotting.plt.close(fig)

    def generate_comprehensive_report(self, output_pdf='comprehensive_report.pdf'):
        """Generate a comprehensive PDF report with all visualizations"""
        pdf_path = os.path.join(self.comprehensive_dir, output_pdf)
        
        with plotting.PdfPages(pdf_path) as pdf, tempfile.TemporaryDirectory() as temp_dir:
            # Title page
            fig = plotting.plt.figure(figsize=(11, 8.5))
            fig.text(0.5, 0.6, 'Molecular Database Quality Assessment', 
                    ha='right', fontsize=24, fontweight='bold')
            fig.text(0.5, 0.5, 'Comprehensive Report', 
                    ha='right', fontsize=18)
            fig.text(0.5, 0.4, f'Generated: {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M")}', 
                    ha='right', fontsize=12)
            fig.text(0.5, 0.3, f'Number of Databases: {len(self.scoring_results)}', 
                    ha='right', fontsize=12)
            pdf.savefig(fig, bbox_inches='tight')
            plotting.plt.close()
            
            # Generate all plots and save to PDF
            print("Generating comprehensive report...")
            
            # List of all plot functions to call
            plot_functions = [
                ('Invalid SMILES Comparison', lambda save_path: self.plot_invalid_smiles_comparison('both', save_path=save_path)),
                ('Invalid SMILES Error Categories', lambda save_path: self.plot_invalid_smiles_error_categories('invalid', save_path=save_path)),
                ('Non-standardized SMILES Comparison', lambda save_path: self.plot_non_standardized_smiles_comparison('both', save_path=save_path)),
                ('Representation Consistency Categories', lambda save_path: self.plot_representation_consistency_categories('issue-only', save_path=save_path)),
                ('Undefined Chirality Count Distribution', lambda save_path: self.plot_undefined_chirality_comparison('count', use_log=False, save_path=save_path)),
                ('Undefined Chirality Count Distribution (Log Scale)', lambda save_path: self.plot_undefined_chirality_comparison('count', use_log=True, save_path=save_path)),
                ('Undefined Chirality Ratio Distribution', lambda save_path: self.plot_undefined_chirality_comparison('ratio', use_log=False, save_path=save_path)),
                ('Undefined Double-Bond E/Z Count Distribution', lambda save_path: self.plot_undefined_double_bond_comparison('count', use_log=False, save_path=save_path)),
                ('Undefined Double-Bond E/Z Count Distribution (Log Scale)', lambda save_path: self.plot_undefined_double_bond_comparison('count', use_log=True, save_path=save_path)),
                ('Undefined Double-Bond E/Z Ratio Distribution', lambda save_path: self.plot_undefined_double_bond_comparison('ratio', use_log=False, save_path=save_path)),
                ('Activity KDE Distribution', self.plot_activity_kde_distribution),
                ('Duplicate Activity Consistency', self.plot_duplicate_activity_consistency),
                ('Molecular Properties Distribution', self.plot_molecular_properties_distribution),
                ('Lipinski Violations', self.plot_lipinski_violations),
                ('Data Imbalance Lorenz Curves', self.plot_data_imbalance_lorenz),
                ('Class Entropy Comparison', self.plot_class_entropy_comparison),
                ('Activity Boxplot Comparison', self.plot_activity_boxplot_comparison),
                ('Quality Metrics Heatmap', self.plot_quality_heatmap),
                ('Quality Scores Waterfall', self.plot_quality_waterfall),
                ('Runtime Total Line', self.plot_total_runtime_line),
                ('Runtime Category Percentages', self.plot_category_runtime_percentage_bars),
                ('Runtime Category Composition (100% Stacked)', self.plot_category_runtime_stacked_percentage_bars),
                ('Runtime Metric Heatmap', self.plot_metric_runtime_heatmap),
            ]
            
            for idx, (title, plot_func) in enumerate(plot_functions):
                try:
                    print(f"  Adding {title}...")
                    save_path = os.path.join(temp_dir, f"report_plot_{idx:02d}.png")
                    self._append_saved_plot_to_pdf(pdf, title, plot_func, save_path)
                except Exception as e:
                    print(f"    Failed to generate {title}: {str(e)}")
            
            print(f"Comprehensive report saved to: {pdf_path}")
        
        return pdf_path
    
    # Keep all original methods from the original script...
    # (Include all the methods from the original script that aren't being modified)
    def generate_all_plots(self, show_values=True, save_prefix=None, plot_chirality_for=None,
                       top_scaffolds_n=10, activity_col=None, class_col=None,
                       tsne_samples=5000, include_tsne=True, include_sankey=False):
        """
        Generate all types of plots

        Args:
            show_values: Whether to show values on bar charts
            save_prefix: Prefix for saved files
            plot_chirality_for: List of databases to plot chiral molecules for
            top_scaffolds_n: Number of top scaffolds to show per database
            activity_col: Activity column to plot
            class_col: Class column to plot
            tsne_samples: Number of samples per database for t-SNE visualization.
                         Defaults to 5000. If None, use all molecules.
            include_tsne: Whether to include t-SNE plots (computationally expensive)
            include_sankey: Whether to include Sankey plots. Defaults to False.
        """

        print("\n" + "="*60)
        print("GENERATING COMPREHENSIVE VISUALIZATION REPORT")
        print("="*60)

        # Original comparison plots
        print("\n[1/7] Generating basic comparison plots...")

        # Single canvas bar chart
        save_path = os.path.join(self.comparison_dir, f'{save_prefix}_single_canvas.png') if save_prefix else None
        self.plot_bars_single_canvas(show_values=show_values, save_path=save_path)

        # Multi canvas bar chart
        save_path = os.path.join(self.comparison_dir, f'{save_prefix}_multi_canvas.png') if save_prefix else None
        self.plot_bars_multi_canvas(show_values=show_values, save_path=save_path)

        # Radar chart
        save_path = os.path.join(self.comparison_dir, f'{save_prefix}_radar.png') if save_prefix else None
        self.plot_radar(save_path=save_path)

        # Data quality metrics
        save_path = os.path.join(self.comparison_dir, f'{save_prefix}_data_quality_metrics.png') if save_prefix else None
        self.plot_data_quality_metrics(save_path=save_path)

        # Data quality comparison plots
        print("\n[2/7] Generating data quality comparison plots...")

        # Invalid SMILES comparison
        print("  - Invalid SMILES comparison...")
        self.plot_invalid_smiles_comparison(mode='both')
        self.plot_invalid_smiles_comparison(mode='count')
        self.plot_invalid_smiles_comparison(mode='ratio')
        self.plot_invalid_smiles_error_categories(normalize_by='invalid')
        self.plot_invalid_smiles_error_categories(normalize_by='all')

        # Representation consistency categories
        print("  - Representation consistency categories...")
        self.plot_representation_consistency_categories(normalize_by='issue-only')
        self.plot_representation_consistency_categories(normalize_by='all')

        # Non-standardized SMILES comparison
        print("  - Non-standardized SMILES comparison...")
        self.plot_non_standardized_smiles_comparison(mode='both')
        self.plot_non_standardized_smiles_comparison(mode='count')
        self.plot_non_standardized_smiles_comparison(mode='ratio')

        # Undefined stereochemistry comparison
        print("  - Undefined chirality comparison...")
        self.plot_undefined_chirality_comparison(mode='count', use_log=False)
        self.plot_undefined_chirality_comparison(mode='count', use_log=True)
        self.plot_undefined_chirality_comparison(mode='ratio', use_log=False)
        print("  - Undefined double-bond E/Z comparison...")
        self.plot_undefined_double_bond_comparison(mode='count', use_log=False)
        self.plot_undefined_double_bond_comparison(mode='count', use_log=True)
        self.plot_undefined_double_bond_comparison(mode='ratio', use_log=False)

        # Activity and label consistency
        print("  - Activity consistency analysis...")
        self.plot_duplicate_activity_consistency()
        self.plot_label_conflict_heatmap()

        # Chemical property plots
        print("\n[3/7] Generating chemical property analysis plots...")

        # QED distribution
        print("  - QED distribution...")
        self.plot_qed_distribution()

        # Molecular properties distribution
        print("  - Molecular properties distribution...")
        self.plot_molecular_properties_distribution()

        # Lipinski violations
        print("  - Lipinski rule analysis...")
        self.plot_lipinski_violations()

        # Chemical space visualization (optional due to computational cost)
        if include_tsne:
            print("  - Combined chemical space t-SNE for all databases...")
            self.plot_all_databases_combined_tsne(n_samples_per_db=tsne_samples)

            print("  - Chemical space t-SNE for each database (this may take a while)...")
            self.plot_chemical_space_tsne(n_samples=tsne_samples)

        # Distribution analysis plots
        print("\n[4/7] Generating distribution analysis plots...")

        # Database size comparison
        print("  - Database size comparison...")
        self.plot_database_size_comparison()

        # Activity distribution plots
        print("  - Activity distribution analysis...")
        self.plot_activity_kde_distribution()
        self.plot_activity_boxplot_comparison()
        self.plot_activity_raincloud(activity_col=activity_col)

        # Class distribution analysis
        print("  - Class distribution analysis...")
        self.plot_data_imbalance_lorenz()
        self.plot_class_entropy_comparison()

        # Comprehensive quality assessment
        print("\n[5/7] Generating comprehensive quality assessment plots...")

        # Quality heatmap
        print("  - Quality metrics heatmap...")
        self.plot_quality_heatmap()

        # Quality waterfall
        print("  - Quality scores waterfall...")
        self.plot_quality_waterfall()

        print("\n[6/7] Generating runtime analysis plots...")
        self.generate_runtime_analysis_plots()

        # Individual database plots
        print("\n[7/7] Generating individual database plots...")

        for db_idx, db_name in enumerate(self.scoring_results.keys()):
            print(f"\n  [{db_idx+1}/{len(self.scoring_results)}] Processing {db_name}...")

            # Top scaffolds
            print(f"    - Top {top_scaffolds_n} scaffolds...")
            self.plot_top_scaffolds(db_name, top_n=top_scaffolds_n)

            # Class distribution (if available)
            if self.scoring_results[db_name]['class_cols']:
                print(f"    - Class distribution...")
                self.plot_class_distribution(db_name, class_col=class_col)

            # Undefined chirality molecules (if specified)
            if plot_chirality_for and db_name in plot_chirality_for:
                print(f"    - Undefined chirality molecules...")
                save_path = os.path.join(self.scoring_results[db_name]['output_dir'],
                                       'undefined_chirality_molecules.png')
                if hasattr(self, 'plot_undefined_chirality_molecules'):
                    self.plot_undefined_chirality_molecules(db_name, save_path=save_path)

            print(f"    - Undefined stereochemistry distributions...")
            self.plot_undefined_chirality_distribution(
                db_name,
                save_path=os.path.join(
                    self.scoring_results[db_name]['output_dir'],
                    'undefined_chirality_distribution.png',
                ),
            )
            self.plot_undefined_double_bond_distribution(
                db_name,
                save_path=os.path.join(
                    self.scoring_results[db_name]['output_dir'],
                    'undefined_double_bond_distribution.png',
                ),
            )

            print(f"    - Representation consistency molecules...")
            self.plot_representation_consistency_molecules(db_name)

            # Structural duplication molecules
            print(f"    - Structural duplication molecules...")
            self.plot_structural_duplication_molecules(
                db_name,
                save_path=os.path.join(self.scoring_results[db_name]['output_dir'], 'structural_duplication_molecules.png')
            )

            # Structural duplication molecules - individual images
            print(f"    - Structural duplication molecules (individual)...")
            self.plot_structural_duplication_molecules_individual(db_name)

            # Contradictory label molecules
            print(f"    - Contradictory label molecules...")
            self.plot_contradictory_label_molecules(
                db_name,
                save_path=os.path.join(self.scoring_results[db_name]['output_dir'], 'contradictory_label_molecules.png')
            )

            # Contradictory label molecules - individual images
            print(f"    - Contradictory label molecules (individual)...")
            self.plot_contradictory_label_molecules_individual(db_name)

            if include_sankey:
                self.plot_sankey_for_database(db_name)
            self.plot_undefined_chirality_molecules(db_name, save_path=os.path.join(self.scoring_results[db_name]['output_dir'],
                                       'undefined_chirality_molecules.png'))

        print("\n" + "="*60)
        print("ALL PLOTS GENERATED SUCCESSFULLY!")
        print("="*60)

        # Print summary of output locations
        print("\nOutput locations:")
        print(f"  - Comparison plots: {self.comparison_dir}")
        print(f"  - Quality metrics: {self.quality_dir}")
        print(f"  - Chemical coverage: {self.coverage_dir}")
        print(f"  - Data distribution: {self.distribution_dir}")
        print(f"  - Runtime metrics: {self.runtime_dir}")
        print(f"  - Comprehensive analysis: {self.comprehensive_dir}")

        for db_name in self.scoring_results.keys():
            print(f"  - {db_name} specific plots: {self.scoring_results[db_name]['output_dir']}")
