from .._common import *
from .. import plotting


class ReportingMixin:
    def generate_comprehensive_report(self, output_pdf='comprehensive_report.pdf'):
        """Generate a comprehensive PDF report with all visualizations"""
        pdf_path = os.path.join(self.comprehensive_dir, output_pdf)
        
        with plotting.PdfPages(pdf_path) as pdf:
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
                ('Invalid SMILES Comparison', lambda: self.plot_invalid_smiles_comparison('both')),
                ('Non-standardized SMILES Comparison', lambda: self.plot_non_standardized_smiles_comparison('both')),
                ('Undefined Chirality (Normal)', lambda: self.plot_undefined_chirality_comparison('both', use_log=False)),
                ('Undefined Chirality (Log Scale)', lambda: self.plot_undefined_chirality_comparison('both', use_log=True)),
                ('Activity KDE Distribution', self.plot_activity_kde_distribution),
                ('Duplicate Activity Consistency', self.plot_duplicate_activity_consistency),
                ('Molecular Properties Distribution', self.plot_molecular_properties_distribution),
                ('Lipinski Violations', self.plot_lipinski_violations),
                ('Data Imbalance Lorenz Curves', self.plot_data_imbalance_lorenz),
                ('Class Entropy Comparison', self.plot_class_entropy_comparison),
                ('Activity Boxplot Comparison', self.plot_activity_boxplot_comparison),
                ('Quality Metrics Heatmap', self.plot_quality_heatmap),
                ('Quality Scores Waterfall', self.plot_quality_waterfall),
            ]
            
            for title, plot_func in plot_functions:
                try:
                    print(f"  Adding {title}...")
                    plot_func()
                    # Get the latest figure and add to PDF
                    fig = plotting.plt.gcf()
                    if fig.get_axes():
                        pdf.savefig(fig, bbox_inches='tight')
                    plotting.plt.close()
                except Exception as e:
                    print(f"    Failed to generate {title}: {str(e)}")
            
            print(f"Comprehensive report saved to: {pdf_path}")
        
        return pdf_path
    
    # Keep all original methods from the original script...
    # (Include all the methods from the original script that aren't being modified)
    def generate_all_plots(self, show_values=True, save_prefix=None, plot_chirality_for=None,
                       top_scaffolds_n=10, activity_col=None, class_col=None,
                       tsne_samples=None, include_tsne=True):
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
                         If None, use all molecules (default: None)
            include_tsne: Whether to include t-SNE plots (computationally expensive)
        """

        print("\n" + "="*60)
        print("GENERATING COMPREHENSIVE VISUALIZATION REPORT")
        print("="*60)

        # Original comparison plots
        print("\n[1/6] Generating basic comparison plots...")

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
        print("\n[2/6] Generating data quality comparison plots...")

        # Invalid SMILES comparison
        print("  - Invalid SMILES comparison...")
        self.plot_invalid_smiles_comparison(mode='both')
        self.plot_invalid_smiles_comparison(mode='count')
        self.plot_invalid_smiles_comparison(mode='ratio')

        # Non-standardized SMILES comparison
        print("  - Non-standardized SMILES comparison...")
        self.plot_non_standardized_smiles_comparison(mode='both')
        self.plot_non_standardized_smiles_comparison(mode='count')
        self.plot_non_standardized_smiles_comparison(mode='ratio')

        # Undefined chirality comparison
        print("  - Undefined chirality comparison...")
        self.plot_undefined_chirality_comparison(mode='both', use_log=False)
        self.plot_undefined_chirality_comparison(mode='both', use_log=True)
        self.plot_undefined_chirality_comparison(mode='both', use_log=True, add_inset=True)

        # Activity and label consistency
        print("  - Activity consistency analysis...")
        self.plot_duplicate_activity_consistency()
        self.plot_label_conflict_heatmap()

        # Chemical property plots
        print("\n[3/6] Generating chemical property analysis plots...")

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
        print("\n[4/6] Generating distribution analysis plots...")

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
        print("\n[5/6] Generating comprehensive quality assessment plots...")

        # Quality heatmap
        print("  - Quality metrics heatmap...")
        self.plot_quality_heatmap()

        # Quality waterfall
        print("  - Quality scores waterfall...")
        self.plot_quality_waterfall()

        # Individual database plots
        print("\n[6/6] Generating individual database plots...")

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
        print(f"  - Comprehensive analysis: {self.comprehensive_dir}")

        for db_name in self.scoring_results.keys():
            print(f"  - {db_name} specific plots: {self.scoring_results[db_name]['output_dir']}")

