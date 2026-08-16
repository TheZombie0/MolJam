# MolJam

MolJam is a Python toolkit for molecular database quality scoring, cleaning, and visualization. It is designed for CSV-style molecular datasets that contain SMILES strings and optional activity, label, experimental-method, identifier, name, and time columns.

The project focuses on one practical question: before a molecular dataset is used for modeling, benchmarking, or downstream chemical analysis, how reliable is the database itself?

## Main Features

- Molecular database quality scoring on a normalized 0-100 scale.
- Five assessment dimensions: structural integrity, data quality, experimental information quality, chemical-space coverage, and data distribution.
- Detailed Markdown report generation for each scored database.
- Dataset cleaning based on detected quality issues.
- Multi-database visualization for quality scores, structural issues, molecular properties, chemical space, label distribution, activity distribution, and runtime behavior.
- Column classification utilities that combine keyword and statistical signals to identify useful annotation columns.

## Project Structure

```text
MolJam_final/
├── README.md
└── src/
    └── moljam/
        ├── classification/   # Column-role classifiers
        ├── scoring/          # Scoring, chemistry processing, cleaning, reports
        └── viz/              # Plotting and visualization tools
```

Important files:

- `README.md`: project overview and usage guide.
- `src/moljam/scoring/scorer.py`: core `MoleculeDBScorer` class.
- `src/moljam/scoring/api.py`: convenience `score_database(...)` function.
- `src/moljam/viz/visualizer.py`: `MoleculeDBVisualizer` class for multi-database plots.

## Installation

MolJam requires Python 3.8 or newer.

RDKit is required for scoring and chemistry-related visualization. The recommended installation method is Conda:

```bash
conda create -n moljam python=3.10
conda activate moljam
conda install -c conda-forge rdkit
```

Install MolJam in editable mode from the project directory:

```bash
cd MolJam_final
pip install -e .
```

Install optional dependencies when needed:

```bash
pip install -e ".[scoring]"
pip install -e ".[viz]"
pip install -e ".[ml]"
```

For full local development, a common setup is:

```bash
pip install -e ".[scoring,viz,ml]"
```

Notes:

- `numpy` and `pandas` are base dependencies.
- `scipy` is used by scoring extras.
- `matplotlib`, `seaborn`, `pyecharts`, and `snapshot-selenium` are used by visualization extras.
- `scikit-learn` is used by machine-learning and dimensionality-reduction utilities.
- RDKit is intentionally installed outside `pip` in many environments because Conda packages are usually more reliable.

## Input Data

The primary input is a CSV file or a pandas `DataFrame`.

Required column:

- `smiles`: SMILES string column by default. You can use another name by passing `smiles_col`.

Optional columns:

- `activity_cols`: continuous activity or affinity values, such as `IC50`, `pIC50`, `Ki`, `logS`, or `activity`.
- `class_cols`: binary or multi-class labels, such as `active`, `inactive`, `toxicity`, or `label`.
- `experimental_method_cols`: assay or experimental method descriptors.
- `id_col`: compound identifier.
- `name_col`: compound name.
- `time_col`: time, year, release date, or version column.

Example minimal table:

```text
smiles,compound_id,pIC50,active,assay_type,year
CCO,CMPD_001,5.2,0,binding,2021
CCN,CMPD_002,7.1,1,functional,2022
```

## Scoring System

When `experimental_info=True`, MolJam computes 12 metrics across 5 categories:

| Category | Metrics |
| --- | --- |
| Structural Integrity | Valid SMILES; Representation Consistency; Stereochemistry Completeness |
| Data Quality | Label Consistency; Data Consistency and Reliability |
| Experimental Information Quality | Time Label Availability; Annotation Support Quality; Type Diversity |
| Chemical Space Coverage | Chemical Diversity; Drug-likeness |
| Data Distribution | Data Size; Data Balance and Distribution |

Each metric is scored on a 0-10 scale. Category totals are normalized so the final score is reported on a 0-100 scale.

When `experimental_info=True`, each of the five categories contributes up to 20 points. When `experimental_info=False`, the experimental-information category is skipped and the remaining four categories each contribute up to 25 points.

MolJam also applies a low-score penalty to produce the final adjusted score:

```text
Final Adjusted Score = max(0, Normalized Score - Low Score Penalty)
```

For more detail about the metric pre-calculations, see `metric_formulas_zh.md`.

## Quick Start: Score a CSV File

```python
from pathlib import Path

from moljam import score_database

score, report, scorer = score_database(
    "example_database.csv",
    smiles_col="smiles",
    activity_cols=["pIC50"],
    class_cols=["active"],
    experimental_method_cols=["assay_type"],
    id_col="compound_id",
    name_col=None,
    time_col="year",
    experimental_info=True,
)

print(f"Final adjusted score: {score:.2f}/100")
Path("quality_report.md").write_text(report, encoding="utf-8")
```

`score_database(...)` returns:

- `score`: final adjusted score on a 0-100 scale.
- `report`: Markdown-formatted detailed report.
- `scorer`: the `MoleculeDBScorer` object containing scores, intermediate analysis results, runtime profile, and cleaning methods.

## Direct Use With pandas

```python
import pandas as pd

from moljam import MoleculeDBScorer

df = pd.read_csv("example_database.csv")

scorer = MoleculeDBScorer(
    df,
    smiles_col="smiles",
    activity_cols=["pIC50"],
    class_cols=["active"],
    experimental_method_cols=["assay_type"],
    id_col="compound_id",
    time_col="year",
    use_parallel=True,
)

final_score = scorer.run_all_checks()
detailed_report = scorer.get_detailed_report()

print(final_score)
print(scorer.scores)
print(scorer.analysis_results.keys())
```

## Cleaning a Database

After scoring, the same scorer can remove records with common quality problems:

```python
cleaned_df, cleaning_report = scorer.clean_database(
    remove_invalid_smiles=True,
    remove_undefined_stereochemistry=True,
    remove_conflicting_labels=True,
    remove_consistent_duplicates=True,
)

cleaned_df.to_csv("example_database_cleaned.csv", index=False)
print(cleaning_report)
```

You can also clean and save in one step:

```python
cleaned_df, cleaning_report = scorer.save_cleaned_database(
    "example_database_cleaned.csv"
)
```

The cleaning report records which issue types were removed and gives representative row indices.

## Visualization

Use `MoleculeDBVisualizer` when comparing one or more databases.

```python
from moljam.viz import MoleculeDBVisualizer

viz = MoleculeDBVisualizer(output_dir="plots", use_parallel=True)

viz.add_database(
    "raw",
    "example_database.csv",
    smiles_col="smiles",
    activity_cols=["pIC50"],
    class_cols=["active"],
    experimental_method_cols=["assay_type"],
    id_col="compound_id",
    time_col="year",
)

viz.add_database(
    "cleaned",
    "example_database_cleaned.csv",
    smiles_col="smiles",
    activity_cols=["pIC50"],
    class_cols=["active"],
    experimental_method_cols=["assay_type"],
    id_col="compound_id",
    time_col="year",
)

viz.plot_quality_heatmap(save_path="plots/quality_heatmap.png")
viz.plot_radar_comparison("raw", "cleaned", save_path="plots/radar_comparison.png")
viz.generate_all_plots(save_prefix="example", include_tsne=False)
```

Common visualization outputs include:

- score bar charts and radar charts;
- quality heatmaps and waterfall plots;
- invalid SMILES and non-standardized SMILES comparisons;
- representation-consistency issue summaries;
- undefined chirality and undefined double-bond E/Z distributions;
- activity and class distribution plots;
- molecular property, QED, Lipinski, scaffold, and t-SNE plots;
- runtime line, bar, heatmap, Pareto, treemap, and Gantt-style plots.