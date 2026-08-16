import pandas as pd

from .scorer import MoleculeDBScorer


def score_database(
    db_path,
    smiles_col="smiles",
    activity_cols=None,
    class_cols=None,
    experimental_method_cols=None,
    id_col=None,
    name_col=None,
    time_col=None,
    use_parallel=True,
    experimental_info=True,
    parent_form_backend="dimorphite_dl",
    parent_form_ph=7.4,
    chemaxon_executable="cxcalc",
    dimorphite_python=None,
    dimorphite_conda_env="dimorphite",
):
    """
    Score a given molecular database.

    Returns:
        final_adjusted_score, detailed_report, scorer_object
    """
    try:
        df = pd.read_csv(db_path)
        print(f"Successfully loaded database: {db_path}")
        print(f"Dataset size: {df.shape[0]} rows, {df.shape[1]} columns")
        print(f"Column names: {df.columns.tolist()}")

        # Verify SMILES column exists
        if smiles_col not in df.columns:
            print(f"Error: SMILES column '{smiles_col}' does not exist in the dataset")
            return None, None, None

        # Verify specified columns exist and warn if not
        if activity_cols:
            activity_cols = [activity_cols] if isinstance(activity_cols, str) else activity_cols
            missing_cols = [col for col in activity_cols if col not in df.columns]
            if missing_cols:
                print(f"Warning: Activity column(s) {missing_cols} not found in the dataset")
                # Remove missing columns
                activity_cols = [col for col in activity_cols if col in df.columns]

        if class_cols:
            class_cols = [class_cols] if isinstance(class_cols, str) else class_cols
            missing_cols = [col for col in class_cols if col not in df.columns]
            if missing_cols:
                print(f"Warning: Class column(s) {missing_cols} not found in the dataset")
                # Remove missing columns
                class_cols = [col for col in class_cols if col in df.columns]

        if experimental_method_cols:
            experimental_method_cols = (
                [experimental_method_cols]
                if isinstance(experimental_method_cols, str)
                else experimental_method_cols
            )
            missing_cols = [col for col in experimental_method_cols if col not in df.columns]
            if missing_cols:
                print(
                    f"Warning: Experimental method column(s) {missing_cols} not found in the dataset"
                )
                # Remove missing columns
                experimental_method_cols = [
                    col for col in experimental_method_cols if col in df.columns
                ]

        if id_col is not None and id_col not in df.columns:
            print(f"Warning: ID column '{id_col}' not found in the dataset")
            id_col = None

        if name_col is not None and name_col not in df.columns:
            print(f"Warning: Name column '{name_col}' not found in the dataset")
            name_col = None

        if time_col is not None and time_col not in df.columns:
            print(f"Warning: Time column '{time_col}' not found in the dataset")
            time_col = None

        scorer = MoleculeDBScorer(
            df,
            smiles_col=smiles_col,
            activity_cols=activity_cols,
            class_cols=class_cols,
            experimental_method_cols=experimental_method_cols,
            id_col=id_col,
            name_col=name_col,
            time_col=time_col,
            use_parallel=use_parallel,
            experimental_info=experimental_info,
            parent_form_backend=parent_form_backend,
            parent_form_ph=parent_form_ph,
            chemaxon_executable=chemaxon_executable,
            dimorphite_python=dimorphite_python,
            dimorphite_conda_env=dimorphite_conda_env,
        )

        final_adjusted_score = scorer.run_all_checks()
        report = scorer.get_detailed_report()
        print(report)
        return final_adjusted_score, report, scorer

    except Exception as e:
        print(f"Scoring failed: {str(e)}")
        import traceback

        traceback.print_exc()
        return None, None, None
