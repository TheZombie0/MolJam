from ...scoring.api import score_database
from .._common import *
from .. import plotting


class VisualizerCoreMixin:
    def __init__(
        self,
        output_dir='plots_bbbp_test',
        use_parallel=True,
        color_palette='husl',
        color_count=11,
    ):
        """Initialize the visualizer with default settings

        Parameters:
            output_dir: Output directory for plots
            use_parallel: Whether to use parallel processing (default: True)
            color_palette: Palette configuration passed to seaborn.color_palette
            color_count: Number of colors to generate for database plots
        """
        plotting.ensure_plotting_imports()

        if color_count <= 0:
            raise ValueError("color_count must be a positive integer")

        # Define color palette for different databases
        self.color_palette = color_palette
        self.color_count = color_count
        self.colors = plotting.sns.color_palette(self.color_palette, n_colors=self.color_count)
        self.use_parallel = use_parallel

        # Create output directory structure
        self.output_dir = output_dir
        self.comparison_dir = os.path.join(output_dir, 'comparison')
        self.quality_dir = os.path.join(output_dir, 'quality_metrics')
        self.coverage_dir = os.path.join(output_dir, 'chemical_coverage')
        self.distribution_dir = os.path.join(output_dir, 'data_distribution')
        self.runtime_dir = os.path.join(output_dir, 'runtime_metrics')
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
                "Annotation Support Quality",
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
                "Annotation Support\nQuality",
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
            self.runtime_dir,
            self.comprehensive_dir,
        ]:
            os.makedirs(dir_path, exist_ok=True)
    def add_database(self, name, file_path, smiles_col='smiles', activity_cols=None,
                     class_cols=None, experimental_method_cols=None, id_col=None, name_col=None,
                     time_col=None, experimental_info=True, parent_form_backend='dimorphite_dl',
                     parent_form_ph=7.4, chemaxon_executable='cxcalc',
                     dimorphite_python=None, dimorphite_conda_env='dimorphite'):
        """Add a database and calculate its scores"""
        print(f"\nScoring database: {name}")
        self._ensure_base_dirs()

        if name in self.scoring_results:
            raise ValueError(
                f"Duplicate database name '{name}'. Database names are used as unique keys in "
                "self.scoring_results, so adding this name again would overwrite the previous result. "
                f"Existing databases: {list(self.scoring_results.keys())}. Please use distinct names "
                f"such as '{name}_raw' and '{name}_cleaned'."
            )

        score, report, scorer = score_database(
            file_path,
            smiles_col=smiles_col,
            activity_cols=activity_cols,
            class_cols=class_cols,
            experimental_method_cols=experimental_method_cols,
            id_col=id_col,
            name_col=name_col,
            time_col=time_col,
            use_parallel=self.use_parallel,
            experimental_info=experimental_info,
            parent_form_backend=parent_form_backend,
            parent_form_ph=parent_form_ph,
            chemaxon_executable=chemaxon_executable,
            dimorphite_python=dimorphite_python,
            dimorphite_conda_env=dimorphite_conda_env,
        )
        
        if scorer is not None:
            # Create separate directory for each database
            db_dir = os.path.join(self.output_dir, name)
            os.makedirs(db_dir, exist_ok=True)
            
            self.scoring_results[name] = {
                'scorer': scorer,
                'total_score': score,
                'report': report,
                'output_dir': db_dir,
                'activity_cols': activity_cols,
                'class_cols': class_cols
            }
            print(f"Successfully scored {name}")
        else:
            print(f"Failed to score {name}")
