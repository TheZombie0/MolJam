from .._common import *
from ..structural_breakdowns import (
    build_representation_category_summary,
    empty_representation_details_frame,
    summarize_representation_groups,
)


class QualityTableMixin:
    QUALITY_ISSUE_STAT_COLUMNS = [
        "dataset",
        "number",
        "invalid_smiles_count",
        "invalid_smiles_rate",
        "undefined_chirality_count",
        "undefined_chirality_rate",
        "undefined_double_bond_count",
        "undefined_double_bond_rate",
        "structural_duplication_count",
        "structural_duplication_rate",
        "label_inconsistency_count",
        "label_inconsistency_rate",
    ]

    QUALITY_ISSUE_DISPLAY_COLUMNS = [
        "Dataset",
        "Number",
        "Invalid SMILES",
        "Chiral Centers",
        "Stereogenic Double Bond",
        "Structural Duplication",
        "Label Inconsistency",
    ]

    REPRESENTATION_ISSUE_STAT_COLUMNS = [
        "dataset",
        "representation_issue_count",
        "representation_issue_rate",
        "salt_count",
        "salt_rate",
        "acid_adduct_count",
        "acid_adduct_rate",
        "solvent_stripping_count",
        "solvent_stripping_rate",
        "protonated_count",
        "protonated_rate",
        "deprotonated_count",
        "deprotonated_rate",
        "duplicate_component_count",
        "duplicate_component_rate",
        "other_non_parent_form_count",
        "other_non_parent_form_rate",
    ]

    REPRESENTATION_ISSUE_DISPLAY_COLUMNS = [
        "Dataset",
        "Representation Issues",
        "Salt",
        "Acid Adduct",
        "Solvent Stripping",
        "Protonated",
        "Deprotonated",
        "Duplicate Component",
        "Other Non-parent Form",
    ]

    REPRESENTATION_NON_PARENT_CATEGORY_COLUMNS = [
        ("salt", "Salt", "salt_count", "salt_rate"),
        ("acid adduct", "Acid Adduct", "acid_adduct_count", "acid_adduct_rate"),
        (
            "solvent stripping",
            "Solvent Stripping",
            "solvent_stripping_count",
            "solvent_stripping_rate",
        ),
        ("protonated", "Protonated", "protonated_count", "protonated_rate"),
        ("deprotonated", "Deprotonated", "deprotonated_count", "deprotonated_rate"),
        (
            "duplicate-component",
            "Duplicate Component",
            "duplicate_component_count",
            "duplicate_component_rate",
        ),
        (
            "other non-parent form",
            "Other Non-parent Form",
            "other_non_parent_form_count",
            "other_non_parent_form_rate",
        ),
    ]

    @staticmethod
    def _coerce_int(value):
        if value is None or pd.isna(value):
            return None
        return int(value)

    @staticmethod
    def _parse_percent(value):
        if value is None or pd.isna(value):
            return None
        if isinstance(value, str):
            normalized = value.strip()
            if not normalized:
                return None
            if normalized.endswith("%"):
                normalized = normalized[:-1]
            return float(normalized)
        return float(value)

    @classmethod
    def _extract_issue_metric(cls, section, count_key, rate_key):
        if not isinstance(section, dict):
            return None, None
        count = cls._coerce_int(section.get(count_key))
        rate = cls._parse_percent(section.get(rate_key))
        if count is None and rate is None:
            return None, None
        return count, rate

    @staticmethod
    def _format_count_with_rate(count, rate):
        if count is None or rate is None:
            return "NA"
        return f"{int(count)} ({float(rate):.2f}%)"

    def _representation_summary_tables(self):
        source_summaries = []
        detail_rows = []

        for db_name, results in self.scoring_results.items():
            scorer = results["scorer"]
            completed_checks = getattr(scorer, "completed_checks", set())
            if (
                "check_representation_consistency" not in completed_checks
                and hasattr(scorer, "check_representation_consistency")
            ):
                scorer.check_representation_consistency()

            analysis_results = getattr(scorer, "analysis_results", {})
            valid_smiles = analysis_results.get("Valid SMILES", {})

            input_count = self._coerce_int(getattr(scorer, "num_molecules", None))
            if input_count is None and hasattr(scorer, "df"):
                input_count = len(scorer.df)
            if input_count is None:
                input_count = self._coerce_int(valid_smiles.get("Valid molecules"))
            if input_count is None:
                input_count = 0

            valid_count = None
            if hasattr(scorer, "valid_mols"):
                valid_count = len(scorer.valid_mols)
            if valid_count is None:
                valid_count = self._coerce_int(valid_smiles.get("Valid molecules"))
            if valid_count is None and hasattr(scorer, "valid_df"):
                valid_count = len(scorer.valid_df)
            if valid_count is None:
                valid_count = 0

            invalid_count = self._coerce_int(valid_smiles.get("Invalid molecules"))
            if invalid_count is None and hasattr(scorer, "invalid_indices"):
                invalid_count = len(scorer.invalid_indices)
            if invalid_count is None:
                invalid_count = max(0, input_count - valid_count)

            source_summary, source_detail_rows = summarize_representation_groups(
                db_name,
                input_count=input_count,
                valid_count=valid_count,
                invalid_count=invalid_count,
                groups=list(getattr(scorer, "representation_consistency_groups", [])),
            )
            source_summaries.append(source_summary)
            detail_rows.extend(source_detail_rows)

        source_summary_df = pd.DataFrame(source_summaries)
        details_df = (
            pd.DataFrame(detail_rows)
            if detail_rows
            else empty_representation_details_frame()
        )
        category_summary_df = build_representation_category_summary(details_df, source_summary_df)
        return source_summary_df, category_summary_df

    def build_quality_issue_summary_statistics(self):
        """Build raw statistics for the manuscript-quality issue summary table."""
        if not self.scoring_results:
            print("No scoring results available. Please add databases first.")
            return pd.DataFrame(columns=self.QUALITY_ISSUE_STAT_COLUMNS)

        rows = []
        for db_name, results in self.scoring_results.items():
            scorer = results["scorer"]
            analysis_results = getattr(scorer, "analysis_results", {})

            data_size = analysis_results.get("Data Size", {})
            valid_smiles = analysis_results.get("Valid SMILES", {})
            stereochemistry = analysis_results.get("Stereochemistry Completeness", {})
            data_consistency = analysis_results.get("Data Consistency and Reliability", {})
            structural_duplication = data_consistency.get("Structural Duplication", {})
            label_consistency = analysis_results.get("Label Consistency", {})

            number = self._coerce_int(data_size.get("Total molecules"))
            if number is None:
                number = self._coerce_int(getattr(scorer, "num_molecules", None))
            if number is None and hasattr(scorer, "df"):
                number = len(scorer.df)

            invalid_smiles_count, invalid_smiles_rate = self._extract_issue_metric(
                valid_smiles,
                "Invalid molecules",
                "Invalid rate",
            )
            undefined_chirality_count, undefined_chirality_rate = self._extract_issue_metric(
                stereochemistry,
                "Molecules with undefined chirality",
                "Molecules with undefined chirality ratio",
            )
            undefined_double_bond_count, undefined_double_bond_rate = self._extract_issue_metric(
                stereochemistry,
                "Molecules with undefined double bond stereochemistry",
                "Molecules with undefined double bond ratio",
            )
            structural_duplication_count, structural_duplication_rate = self._extract_issue_metric(
                structural_duplication,
                "Duplicate molecules",
                "Duplication rate",
            )
            label_inconsistency_count, label_inconsistency_rate = self._extract_issue_metric(
                label_consistency,
                "Total molecules with contradictory labels",
                "Overall contradictory label rate",
            )

            rows.append(
                {
                    "dataset": db_name,
                    "number": number,
                    "invalid_smiles_count": invalid_smiles_count,
                    "invalid_smiles_rate": invalid_smiles_rate,
                    "undefined_chirality_count": undefined_chirality_count,
                    "undefined_chirality_rate": undefined_chirality_rate,
                    "undefined_double_bond_count": undefined_double_bond_count,
                    "undefined_double_bond_rate": undefined_double_bond_rate,
                    "structural_duplication_count": structural_duplication_count,
                    "structural_duplication_rate": structural_duplication_rate,
                    "label_inconsistency_count": label_inconsistency_count,
                    "label_inconsistency_rate": label_inconsistency_rate,
                }
            )

        return pd.DataFrame(rows, columns=self.QUALITY_ISSUE_STAT_COLUMNS, dtype=object)

    def build_quality_issue_summary_table(self):
        """Build the formatted database-quality issue summary table."""
        stats_df = self.build_quality_issue_summary_statistics()
        if stats_df.empty:
            return pd.DataFrame(columns=self.QUALITY_ISSUE_DISPLAY_COLUMNS)

        display_rows = []
        for row in stats_df.itertuples(index=False):
            display_rows.append(
                {
                    "Dataset": row.dataset,
                    "Number": "NA" if row.number is None else str(int(row.number)),
                    "Invalid SMILES": self._format_count_with_rate(
                        row.invalid_smiles_count,
                        row.invalid_smiles_rate,
                    ),
                    "Chiral Centers": self._format_count_with_rate(
                        row.undefined_chirality_count,
                        row.undefined_chirality_rate,
                    ),
                    "Stereogenic Double Bond": self._format_count_with_rate(
                        row.undefined_double_bond_count,
                        row.undefined_double_bond_rate,
                    ),
                    "Structural Duplication": self._format_count_with_rate(
                        row.structural_duplication_count,
                        row.structural_duplication_rate,
                    ),
                    "Label Inconsistency": self._format_count_with_rate(
                        row.label_inconsistency_count,
                        row.label_inconsistency_rate,
                    ),
                }
            )

        return pd.DataFrame(display_rows, columns=self.QUALITY_ISSUE_DISPLAY_COLUMNS)

    def build_representation_issue_summary_statistics(self):
        """Build raw statistics for the representation-issue summary table."""
        if not self.scoring_results:
            print("No scoring results available. Please add databases first.")
            return pd.DataFrame(columns=self.REPRESENTATION_ISSUE_STAT_COLUMNS)

        source_summary_df, category_summary_df = self._representation_summary_tables()
        rows = []

        for source_row in source_summary_df.itertuples(index=False):
            source_categories = category_summary_df.loc[
                category_summary_df["source"] == source_row.source
            ].set_index("category")

            row = {
                "dataset": source_row.source,
                "representation_issue_count": int(source_row.issue_molecule_count),
                "representation_issue_rate": float(source_row.issue_fraction_valid) * 100.0,
            }

            for category, _label, count_key, rate_key in self.REPRESENTATION_NON_PARENT_CATEGORY_COLUMNS:
                if category in source_categories.index:
                    row[count_key] = int(source_categories.loc[category, "count"])
                    row[rate_key] = float(source_categories.loc[category, "fraction_of_issue_only"]) * 100.0
                else:
                    row[count_key] = 0
                    row[rate_key] = 0.0

            rows.append(row)

        return pd.DataFrame(rows, columns=self.REPRESENTATION_ISSUE_STAT_COLUMNS, dtype=object)

    def build_representation_issue_summary_table(self):
        """Build the formatted representation-issue summary table."""
        stats_df = self.build_representation_issue_summary_statistics()
        if stats_df.empty:
            return pd.DataFrame(columns=self.REPRESENTATION_ISSUE_DISPLAY_COLUMNS)

        display_rows = []
        for row in stats_df.itertuples(index=False):
            display_row = {
                "Dataset": row.dataset,
                "Representation Issues": self._format_count_with_rate(
                    row.representation_issue_count,
                    row.representation_issue_rate,
                ),
            }

            for _category, label, count_key, rate_key in self.REPRESENTATION_NON_PARENT_CATEGORY_COLUMNS:
                display_row[label] = self._format_count_with_rate(
                    getattr(row, count_key),
                    getattr(row, rate_key),
                )

            display_rows.append(display_row)

        return pd.DataFrame(display_rows, columns=self.REPRESENTATION_ISSUE_DISPLAY_COLUMNS)

    def export_quality_issue_summary_table(self, save_path=None, caption=None, label=None):
        """Export the formatted quality issue summary table as a LaTeX booktabs table."""
        self._ensure_base_dirs()

        if save_path is None:
            save_path = os.path.join(
                self.comparison_dir,
                "quality_issue_summary_table.tex",
            )

        output_dir = os.path.dirname(save_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        table_df = self.build_quality_issue_summary_table()
        latex = table_df.to_latex(
            index=False,
            escape=True,
            column_format="l" + "c" * (len(table_df.columns) - 1),
            caption=caption,
            label=label,
        )

        with open(save_path, "w", encoding="utf-8") as handle:
            handle.write(latex)

        print(f"Saved quality issue summary table to: {save_path}")
        return save_path

    def export_representation_issue_summary_table(self, save_path=None, caption=None, label=None):
        """Export the formatted representation-issue summary table as a LaTeX booktabs table."""
        self._ensure_base_dirs()

        if save_path is None:
            save_path = os.path.join(
                self.comparison_dir,
                "representation_issue_summary_table.tex",
            )

        output_dir = os.path.dirname(save_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        table_df = self.build_representation_issue_summary_table()
        latex = table_df.to_latex(
            index=False,
            escape=True,
            column_format="l" + "c" * (len(table_df.columns) - 1),
            caption=caption,
            label=label,
        )

        with open(save_path, "w", encoding="utf-8") as handle:
            handle.write(latex)

        print(f"Saved representation issue summary table to: {save_path}")
        return save_path
