"""MCP tool entry points. Each function returns a JSON-serializable dict."""

import traceback
from typing import List, Optional

from .scoring.api import score_database, score_database_json


def mcp_score_database(
    db_path: str,
    smiles_col: str = "smiles",
    activity_cols: Optional[List[str]] = None,
    class_cols: Optional[List[str]] = None,
    experimental_method_cols: Optional[List[str]] = None,
    id_col: Optional[str] = None,
    name_col: Optional[str] = None,
    time_col: Optional[str] = None,
    use_parallel: bool = True,
    include_experimental_info: bool = False,
) -> dict:
    """Score a molecular database. Returns JSON-serializable dict."""
    report = score_database_json(
        db_path,
        smiles_col=smiles_col,
        activity_cols=activity_cols,
        class_cols=class_cols,
        experimental_method_cols=experimental_method_cols,
        id_col=id_col,
        name_col=name_col,
        time_col=time_col,
        use_parallel=use_parallel,
        include_experimental_info=include_experimental_info,
    )
    return report.to_dict()


def mcp_classify_columns(
    db_path: str,
    smiles_col: str = "smiles",
) -> dict:
    """Classify CSV columns as useful/excluded/unknown. Returns JSON-serializable dict."""
    try:
        import pandas as pd

        from .classification.hybrid import HybridColumnClassifier

        df = pd.read_csv(db_path)
        classifier = HybridColumnClassifier()
        raw = classifier.classify_columns(df, smiles_col=smiles_col)
        # Convert tuple entries to dicts
        result = {}
        for category, entries in raw.items():
            result[category] = [
                {
                    "column": entry[0],
                    "reason": entry[1],
                    "confidence": entry[2],
                }
                for entry in entries
            ]
        return {"status": "success", "classifications": result}
    except Exception as e:
        return {"status": "error", "error": str(e), "traceback": traceback.format_exc()}


def mcp_clean_database(
    db_path: str,
    output_path: Optional[str] = None,
    smiles_col: str = "smiles",
    activity_cols: Optional[List[str]] = None,
    class_cols: Optional[List[str]] = None,
    experimental_method_cols: Optional[List[str]] = None,
    id_col: Optional[str] = None,
    name_col: Optional[str] = None,
    time_col: Optional[str] = None,
    use_parallel: bool = True,
    include_experimental_info: bool = False,
    remove_invalid_smiles: bool = True,
    remove_undefined_stereochemistry: bool = True,
    remove_conflicting_labels: bool = True,
    remove_consistent_duplicates: bool = True,
) -> dict:
    """Clean a molecular database. Returns JSON-serializable cleaning report."""
    try:
        final_score, report_text, scorer = score_database(
            db_path,
            smiles_col=smiles_col,
            activity_cols=activity_cols,
            class_cols=class_cols,
            experimental_method_cols=experimental_method_cols,
            id_col=id_col,
            name_col=name_col,
            time_col=time_col,
            use_parallel=use_parallel,
            include_experimental_info=include_experimental_info,
        )
        if scorer is None:
            return {"status": "error", "error": "Scoring failed. Cannot clean database."}

        cleaned_df, cleaning_report = scorer.clean_database(
            remove_invalid_smiles=remove_invalid_smiles,
            remove_undefined_stereochemistry=remove_undefined_stereochemistry,
            remove_conflicting_labels=remove_conflicting_labels,
            remove_consistent_duplicates=remove_consistent_duplicates,
            verbose=False,
        )

        if output_path:
            cleaned_df.to_csv(output_path, index=False)

        return {
            "status": "success",
            "cleaning_report": cleaning_report,
            "output_path": output_path,
            "original_score": final_score,
        }
    except Exception as e:
        return {"status": "error", "error": str(e), "traceback": traceback.format_exc()}
