"""MCP Server for MolJam. Exposes scoring, classification, and cleaning tools via FastMCP."""

import json
from typing import Optional

from mcp.server.fastmcp import FastMCP

mcp = FastMCP(
    "MolJam",
    instructions="MolJam 是分子数据库质量评分工具包。可以评分、分类列、清洗分子数据集。",
)


@mcp.tool()
def score_database(
    db_path: str,
    smiles_col: str = "smiles",
    activity_cols: list[str] | None = None,
    class_cols: list[str] | None = None,
    experimental_method_cols: list[str] | None = None,
    id_col: str | None = None,
    name_col: str | None = None,
    time_col: str | None = None,
    include_experimental_info: bool = False,
) -> str:
    """Score a molecular database for quality across 4-5 dimensions.

    Evaluates structural integrity, data quality, chemical space coverage,
    data distribution, and optionally experimental information quality.
    Returns a JSON object with the final score (0-100), per-category scores,
    dataset info, and a concise analysis summary.

    Args:
        db_path: Path to a CSV file containing the molecular database.
        smiles_col: Name of the column containing SMILES strings.
        activity_cols: Column names for activity/potency values.
        class_cols: Column names for classification labels.
        experimental_method_cols: Column names for experimental method info.
        id_col: Column name for molecule identifiers.
        name_col: Column name for molecule names.
        time_col: Column name for time information.
        include_experimental_info: Whether to score experimental info quality.
    """
    try:
        from .scoring.api import score_database_json

        report = score_database_json(
            db_path,
            smiles_col=smiles_col,
            activity_cols=activity_cols,
            class_cols=class_cols,
            experimental_method_cols=experimental_method_cols,
            id_col=id_col,
            name_col=name_col,
            time_col=time_col,
            use_parallel=True,
            include_experimental_info=include_experimental_info,
        )

        if report.final_score is None:
            return json.dumps(
                {"status": "error", "error": "Scoring failed. Check that db_path exists and smiles_col is correct."},
                ensure_ascii=False,
            )

        # Build concise result (no report_text — too verbose for agents)
        result = {
            "status": "success",
            "final_score": report.final_score,
            "scores": report.scores.to_dict() if report.scores else None,
            "dataset_info": report.dataset_info,
        }

        # Add concise analysis summary instead of full analysis_results
        if report.analysis_results:
            summary = {}
            for key, value in report.analysis_results.items():
                if isinstance(value, dict) and len(str(value)) > 500:
                    # Trim large nested dicts to key counts
                    summary[key] = {k: v for k, v in list(value.items())[:10]}
                    if len(value) > 10:
                        summary[key]["_truncated"] = f"{len(value)} items total"
                else:
                    summary[key] = value
            result["analysis_summary"] = summary

        return json.dumps(result, ensure_ascii=False, default=str, indent=2)

    except Exception as e:
        return json.dumps({"status": "error", "error": str(e)}, ensure_ascii=False)


@mcp.tool()
def classify_columns(
    db_path: str,
    smiles_col: str = "smiles",
) -> str:
    """Classify columns in a molecular CSV as useful, excluded, or unknown.

    Analyzes each column in the CSV to determine if it contains useful data
    (activity values, class labels, identifiers, etc.), should be excluded,
    or has unknown purpose. Helps users understand their dataset structure.

    Args:
        db_path: Path to a CSV file containing the molecular database.
        smiles_col: Name of the column containing SMILES strings.
    """
    try:
        import pandas as pd

        from .classification.hybrid import HybridColumnClassifier

        df = pd.read_csv(db_path)
        classifier = HybridColumnClassifier()
        raw = classifier.classify_columns(df, smiles_col=smiles_col)

        classifications = {}
        for category, entries in raw.items():
            classifications[category] = [
                {
                    "column": entry[0],
                    "reason": entry[1],
                    "confidence": entry[2],
                }
                for entry in entries
            ]

        return json.dumps(
            {"status": "success", "classifications": classifications},
            ensure_ascii=False,
            default=str,
            indent=2,
        )

    except Exception as e:
        return json.dumps({"status": "error", "error": str(e)}, ensure_ascii=False)


@mcp.tool()
def clean_database(
    db_path: str,
    output_path: str,
    smiles_col: str = "smiles",
    activity_cols: list[str] | None = None,
    class_cols: list[str] | None = None,
    remove_invalid_smiles: bool = True,
    remove_undefined_stereochemistry: bool = True,
    remove_conflicting_labels: bool = True,
    remove_consistent_duplicates: bool = True,
) -> str:
    """Clean a molecular database by removing quality issues.

    Runs scoring first to identify issues, then removes rows with problems
    such as invalid SMILES, undefined stereochemistry, conflicting labels,
    and consistent duplicates. Saves the cleaned dataset to output_path.

    Args:
        db_path: Path to the input CSV file.
        output_path: Path where the cleaned CSV will be saved (required).
        smiles_col: Name of the column containing SMILES strings.
        activity_cols: Column names for activity/potency values.
        class_cols: Column names for classification labels.
        remove_invalid_smiles: Remove rows with invalid SMILES.
        remove_undefined_stereochemistry: Remove rows with undefined stereo.
        remove_conflicting_labels: Remove rows with conflicting activity labels.
        remove_consistent_duplicates: Remove exact duplicate molecules.
    """
    try:
        from .scoring.api import score_database as _score_database

        final_score, _report_text, scorer = _score_database(
            db_path,
            smiles_col=smiles_col,
            activity_cols=activity_cols,
            class_cols=class_cols,
            use_parallel=True,
        )

        if scorer is None:
            return json.dumps(
                {"status": "error", "error": "Scoring failed. Cannot clean database."},
                ensure_ascii=False,
            )

        cleaned_df, cleaning_report = scorer.clean_database(
            remove_invalid_smiles=remove_invalid_smiles,
            remove_undefined_stereochemistry=remove_undefined_stereochemistry,
            remove_conflicting_labels=remove_conflicting_labels,
            remove_consistent_duplicates=remove_consistent_duplicates,
            verbose=False,
        )

        cleaned_df.to_csv(output_path, index=False)

        return json.dumps(
            {
                "status": "success",
                "cleaning_report": cleaning_report,
                "output_path": output_path,
                "original_score": final_score,
                "rows_before": len(scorer.df),
                "rows_after": len(cleaned_df),
            },
            ensure_ascii=False,
            default=str,
            indent=2,
        )

    except Exception as e:
        return json.dumps({"status": "error", "error": str(e)}, ensure_ascii=False)


def main():
    """Entry point for MCP server (stdio transport)."""
    mcp.run()


if __name__ == "__main__":
    main()
