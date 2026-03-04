import warnings

from rdkit.rdBase import DisableLog

from ...scoring.api import score_database
from .._common import *
from .. import plotting


class VisualizerCoreMixin:
    def __init__(self, output_dir='plots_bbbp_test', use_parallel=True):
        """Initialize the visualizer with default settings

        Parameters:
            output_dir: Output directory for plots
            use_parallel: Whether to use parallel processing (default: True)
        """
        # Suppress RDKit warnings
        warnings.filterwarnings('ignore', module='rdkit')
        DisableLog('rdApp.warning')
        DisableLog('rdApp.error')

        plotting.ensure_plotting_imports()

        # Define color palette for different databases
        self.colors = plotting.sns.color_palette("husl", n_colors=11)
        self.use_parallel = use_parallel

        # Create output directory structure
        self.output_dir = output_dir
        self.comparison_dir = os.path.join(output_dir, 'comparison')
        self.quality_dir = os.path.join(output_dir, 'quality_metrics')
        self.coverage_dir = os.path.join(output_dir, 'chemical_coverage')
        self.distribution_dir = os.path.join(output_dir, 'data_distribution')
        self.comprehensive_dir = os.path.join(output_dir, 'comprehensive')
        
        # Define the 5 main categories and their sub-metrics
        self.categories = {
            "Structural Integrity": [
                "Valid SMILES",
                "Representation Consistency",
                "Stereochemistry Completeness"
            ],
            "Data Quality": [
                "Label Consistency",
                "Data Consistency and Reliability"
            ],
            "Experimental Information Quality": [
                "Time Label Availability",
                "Useful Column Quality",
                "Classification Confidence",
                "Type Diversity"
            ],
            "Chemical Space Coverage": [
                "Chemical Diversity",
                "Drug-likeness"
            ],
            "Data Distribution": [
                "Data Size",
                "Data Balance and Distribution"
            ]
        }

        self.categories_labels = {
            "Structural Integrity": [
                "Valid\nSMILES",
                "Representation\nConsistency",
                "Stereochemistry\nCompleteness"
            ],
            "Data Quality": [
                "Label\nConsistency",
                "Data Consistency\nand Reliability"
            ],
            "Experimental Information Quality": [
                "Time Label\nAvailability",
                "Useful Column\nQuality",
                "Classification\nConfidence",
                "Type\nDiversity"
            ],
            "Chemical Space Coverage": [
                "Chemical\nDiversity",
                "Drug-likeness"
            ],
            "Data Distribution": [
                "Data\nSize",
                "Data Balance and\nDistribution"
            ]
        }
        # Store scoring results
        self.scoring_results = {}

        # Cache for computed properties
        self.property_cache = {}

        # Store axis ranges for consistency between t-SNE plots
        self.tsne_axis_ranges = {'xlim': None, 'ylim': None}
    def _ensure_base_dirs(self):
        for dir_path in [
            self.comparison_dir,
            self.quality_dir,
            self.coverage_dir,
            self.distribution_dir,
            self.comprehensive_dir,
        ]:
            os.makedirs(dir_path, exist_ok=True)
    def add_database(self, name, file_path, smiles_col='smiles', activity_cols=None,
                     class_cols=None, experimental_method_cols=None, id_col=None, name_col=None, time_col=None, include_experimental_info=False):
        """Add a database and calculate its scores"""
        print(f"\nScoring database: {name}")
        self._ensure_base_dirs()
        score, report, scorer = score_database(
            file_path,
            smiles_col=smiles_col,
            activity_cols=activity_cols,
            class_cols=class_cols,
            experimental_method_cols=experimental_method_cols,
            id_col=id_col,
            name_col=name_col,
            time_col=time_col,
            include_experimental_info=include_experimental_info
        )
        
        if scorer is not None:
            # Create separate directory for each database
            db_dir = os.path.join(self.output_dir, name)
            os.makedirs(db_dir, exist_ok=True)
            
            self.scoring_results[name] = {
                'scorer': scorer,
                'snapshot': scorer.to_snapshot(),
                'total_score': score,
                'report': report,
                'output_dir': db_dir,
                'activity_cols': activity_cols,
                'class_cols': class_cols
            }
            print(f"Successfully scored {name}")
        else:
            print(f"Failed to score {name}")

