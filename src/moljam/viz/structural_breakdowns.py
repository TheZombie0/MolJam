from __future__ import annotations

import io
import logging
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase

from . import plotting


DEFAULT_SMILES_COLUMNS = (
    "smiles",
    "SMILES",
    "canonical_smiles",
    "original_smiles",
    "parent_smiles",
)


INVALID_SMILES_SOURCE_NOTE = (
    "Category mapping based on the six RDKit-derived error groups cited in "
    "Skinnider 2024 (Nat. Mach. Intell.) and implemented in the "
    "UnCorrupt SMILES repository (analysis/figure.py::cat_errors)."
)

INVALID_SMILES_CATEGORY_ORDER = [
    "aromaticity error",
    "unclosed ring",
    "parentheses error",
    "valence error",
    "syntax error",
    "bond exists",
]

INVALID_SMILES_CATEGORY_LABELS = {
    "aromaticity error": "Aromaticity\nError",
    "unclosed ring": "Unclosed\nRing",
    "parentheses error": "Parentheses\nError",
    "valence error": "Valence\nError",
    "syntax error": "Syntax\nError",
    "bond exists": "Bond Already\nExists",
}

INVALID_SMILES_CATEGORY_COLORS = {
    "aromaticity error": "#7B66FF",
    "unclosed ring": "#FF6B6B",
    "parentheses error": "#FFB84D",
    "valence error": "#2E86AB",
    "syntax error": "#6C757D",
    "bond exists": "#2F9E44",
}

INVALID_SMILES_NORMALIZE_LABELS = {
    "invalid": "Composition of invalid strings (%)",
    "all": "Proportion of all strings (%)",
}

INVALID_SMILES_NORMALIZE_FILENAMES = {
    "invalid": "by_invalid",
    "all": "by_all",
}

REPRESENTATION_CATEGORY_ORDER = [
    "parent form",
    "salt",
    "acid adduct",
    "solvent stripping",
    "protonated",
    "deprotonated",
    "duplicate-component",
    "other non-parent form",
]

REPRESENTATION_CATEGORY_LABELS = {
    "parent form": "Parent Form",
    "salt": "Salt",
    "acid adduct": "Acid Adduct",
    "solvent stripping": "Solvent Stripping",
    "protonated": "Protonated",
    "deprotonated": "Deprotonated",
    "duplicate-component": "Duplicate Component",
    "other non-parent form": "Other Non-parent Form",
}

REPRESENTATION_CATEGORY_COLORS = {
    "parent form": "#9AA5B1",
    "salt": "#FF6B6B",
    "acid adduct": "#C77DFF",
    "solvent stripping": "#4D96FF",
    "protonated": "#2BB673",
    "deprotonated": "#00B8D9",
    "duplicate-component": "#FFB84D",
    "other non-parent form": "#8D99AE",
}

REPRESENTATION_NORMALIZE_LABELS = {
    "issue-only": "Composition of issue-group molecules (%)",
    "all": "Proportion of valid molecules (%)",
}

REPRESENTATION_NORMALIZE_FILENAMES = {
    "issue-only": "by_issue_only",
    "all": "by_all",
}

TEXT_DARK = "#1F1F1F"
TEXT_MUTED = "#5C534B"
BAR_EDGE = "#FFFFFF"


@dataclass
class SmilesInspection:
    row_index: int
    smiles: str
    is_valid: bool
    category: str | None
    rdkit_log: str


class RDKitLogCapture:
    """Route RDKit C++ logs into a temporary in-memory Python logger."""

    def __enter__(self):
        self.logger = logging.getLogger("rdkit")
        self.previous_handlers = list(self.logger.handlers)
        self.previous_level = self.logger.level
        self.previous_propagate = self.logger.propagate

        self.stream = io.StringIO()
        self.handler = logging.StreamHandler(self.stream)
        self.handler.setLevel(logging.ERROR)
        self.logger.handlers = [self.handler]
        self.logger.setLevel(logging.ERROR)
        self.logger.propagate = False
        rdBase.LogToPythonLogger()
        return self

    def inspect(self, smiles: str):
        self.stream.seek(0)
        self.stream.truncate(0)
        mol = Chem.MolFromSmiles(smiles)
        self.handler.flush()
        return mol, self.stream.getvalue()

    def __exit__(self, exc_type, exc, tb):
        self.logger.handlers = self.previous_handlers
        self.logger.setLevel(self.previous_level)
        self.logger.propagate = self.previous_propagate
        rdBase.LogToCppStreams()
        return False


def detect_smiles_column(df: pd.DataFrame, requested: str | None) -> str:
    if requested:
        if requested not in df.columns:
            raise ValueError(f"SMILES column '{requested}' not found. Available columns: {list(df.columns)}")
        return requested

    for candidate in DEFAULT_SMILES_COLUMNS:
        if candidate in df.columns:
            return candidate

    if len(df.columns) == 1:
        return df.columns[0]

    raise ValueError("Could not infer the SMILES column. Pass --smiles-column explicitly.")


def normalize_smiles_rows(
    raw_df: pd.DataFrame,
    smiles_column: str | None = None,
) -> tuple[pd.DataFrame, int, str]:
    selected_column = detect_smiles_column(raw_df, smiles_column)
    valid_string_mask = raw_df[selected_column].apply(
        lambda value: isinstance(value, str) and bool(value.strip())
    )
    excluded_count = int((~valid_string_mask).sum())
    filtered_df = raw_df.loc[valid_string_mask, [selected_column]].copy()
    filtered_df.rename(columns={selected_column: "smiles"}, inplace=True)
    filtered_df.insert(0, "row_index", filtered_df.index.astype(int))
    filtered_df["smiles"] = filtered_df["smiles"].str.strip()
    return filtered_df.reset_index(drop=True), excluded_count, selected_column


def load_smiles_table(input_path: Path, smiles_column: str | None) -> tuple[pd.DataFrame, int, str]:
    suffix = input_path.suffix.lower()
    if suffix in {".txt", ".smi"}:
        smiles = [line.strip() for line in input_path.read_text().splitlines() if line.strip()]
        raw_df = pd.DataFrame({"smiles": smiles})
        return normalize_smiles_rows(raw_df, "smiles")

    separator = "\t" if suffix in {".tsv", ".tab"} else ","
    raw_df = pd.read_csv(input_path, sep=separator)
    return normalize_smiles_rows(raw_df, smiles_column)


def normalize_rdkit_log(log_text: str) -> str:
    cleaned_lines = []
    for line in log_text.splitlines():
        line = re.sub(r"^\[\d{2}:\d{2}:\d{2}\]\s*", "", line).strip()
        if line:
            cleaned_lines.append(line)
    return " | ".join(cleaned_lines)


def classify_rdkit_error(log_text: str) -> str:
    normalized = normalize_rdkit_log(log_text).lower()

    if not normalized:
        return "syntax error"
    if "kekulize" in normalized or "marked aromatic" in normalized:
        return "aromaticity error"
    if "unclosed ring" in normalized:
        return "unclosed ring"
    if "parentheses" in normalized:
        return "parentheses error"
    if "valence" in normalized:
        return "valence error"
    if (
        "ring closure" in normalized
        or "duplicates bond" in normalized
        or "duplicated ring closure" in normalized
        or "bond already exists" in normalized
        or "bond exists" in normalized
    ):
        return "bond exists"
    if "syntax error" in normalized or "failed parsing" in normalized or "parse error" in normalized:
        return "syntax error"
    return "syntax error"


def inspect_smiles_rows(rows: pd.DataFrame) -> list[SmilesInspection]:
    inspections: list[SmilesInspection] = []
    with RDKitLogCapture() as capture:
        for row in rows.itertuples(index=False):
            smiles = row.smiles if isinstance(row.smiles, str) else ""
            smiles = smiles.strip()
            if not smiles:
                inspections.append(
                    SmilesInspection(
                        row_index=int(row.row_index),
                        smiles="" if row.smiles is None else str(row.smiles),
                        is_valid=False,
                        category="syntax error",
                        rdkit_log="",
                    )
                )
                continue

            mol, log_text = capture.inspect(smiles)
            category = None if mol is not None else classify_rdkit_error(log_text)
            inspections.append(
                SmilesInspection(
                    row_index=int(row.row_index),
                    smiles=smiles,
                    is_valid=mol is not None,
                    category=category,
                    rdkit_log=normalize_rdkit_log(log_text),
                )
            )
    return inspections


def inspections_to_frame(inspections: Iterable[SmilesInspection]) -> pd.DataFrame:
    return pd.DataFrame(asdict(inspection) for inspection in inspections)


def empty_invalid_smiles_details_frame() -> pd.DataFrame:
    return pd.DataFrame(
        columns=["source", "row_index", "smiles", "is_valid", "category", "rdkit_log"]
    )


def summarize_invalid_smiles_results(
    details_df: pd.DataFrame,
    *,
    total_count: int | None = None,
) -> tuple[pd.DataFrame, dict[str, float]]:
    total_count = int(total_count if total_count is not None else len(details_df))
    if "is_valid" in details_df.columns:
        invalid_df = details_df.loc[~details_df["is_valid"]].copy()
    else:
        invalid_df = details_df.copy()

    invalid_count = int(len(invalid_df))
    total_count = max(total_count, invalid_count)
    valid_count = total_count - invalid_count

    summary_rows = []
    for category in INVALID_SMILES_CATEGORY_ORDER:
        category_count = int((invalid_df["category"] == category).sum()) if not invalid_df.empty else 0
        summary_rows.append(
            {
                "category": category,
                "display_label": INVALID_SMILES_CATEGORY_LABELS[category].replace("\n", " "),
                "count": category_count,
                "fraction_of_invalid": (category_count / invalid_count) if invalid_count else 0.0,
                "fraction_of_all": (category_count / total_count) if total_count else 0.0,
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    metadata = {
        "total_count": total_count,
        "valid_count": valid_count,
        "invalid_count": invalid_count,
        "valid_fraction": (valid_count / total_count) if total_count else 0.0,
        "invalid_fraction": (invalid_count / total_count) if total_count else 0.0,
    }
    return summary_df, metadata


def build_invalid_smiles_source_level_summaries(
    source_payloads: Mapping[str, Mapping[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    source_rows = []
    category_rows = []

    for source_name, payload in source_payloads.items():
        total_count = int(payload.get("total_count", 0))
        details_df = payload.get("details_df")
        if not isinstance(details_df, pd.DataFrame):
            details_df = empty_invalid_smiles_details_frame()

        summary_df, metadata = summarize_invalid_smiles_results(details_df, total_count=total_count)
        source_rows.append(
            {
                "source": source_name,
                "total_count": metadata["total_count"],
                "valid_count": metadata["valid_count"],
                "invalid_count": metadata["invalid_count"],
                "valid_fraction": metadata["valid_fraction"],
                "invalid_fraction": metadata["invalid_fraction"],
            }
        )

        source_summary = summary_df.copy()
        source_summary["source"] = source_name
        source_summary["total_count"] = metadata["total_count"]
        source_summary["invalid_count"] = metadata["invalid_count"]
        category_rows.append(source_summary)

    source_level_df = pd.DataFrame(source_rows)
    per_source_category_df = (
        pd.concat(category_rows, ignore_index=True) if category_rows else pd.DataFrame()
    )
    return source_level_df, per_source_category_df


def build_invalid_smiles_plot_table(
    source_level_df: pd.DataFrame,
    per_source_category_df: pd.DataFrame,
    normalize_by: str,
) -> pd.DataFrame:
    fraction_column = f"fraction_of_{normalize_by}"
    rows = []

    for source_row in source_level_df.itertuples(index=False):
        category_subset = per_source_category_df.loc[
            per_source_category_df["source"] == source_row.source
        ].set_index("category")

        row = {
            "source": source_row.source,
            "total_count": int(source_row.total_count),
            "invalid_count": int(source_row.invalid_count),
            "invalid_fraction": float(source_row.invalid_fraction),
        }

        stack_total = 0.0
        for category in INVALID_SMILES_CATEGORY_ORDER:
            value = float(category_subset.loc[category, fraction_column]) if category in category_subset.index else 0.0
            row[category] = value
            stack_total += value

        row["stack_total"] = stack_total
        rows.append(row)

    plot_df = pd.DataFrame(rows)
    if not plot_df.empty:
        plot_df["source"] = pd.Categorical(
            plot_df["source"],
            categories=source_level_df["source"].tolist(),
            ordered=True,
        )
        plot_df = plot_df.sort_values("source").reset_index(drop=True)
        plot_df["source"] = plot_df["source"].astype(str)
    return plot_df


def render_invalid_smiles_stacked_horizontal_figure(
    figure_png_path: Path | str,
    plot_df: pd.DataFrame,
    normalize_by: str,
    title: str,
    panel_label: str | None = None,
    figure_svg_path: Path | str | None = None,
):
    plotting.ensure_plotting_imports()
    from matplotlib import transforms

    n_sources = len(plot_df)
    fig_height = max(3.2, 1.4 + 0.75 * max(n_sources, 1))
    fig, ax = plotting.plt.subplots(figsize=(8.8, fig_height))

    y_positions = np.arange(n_sources)
    left = np.zeros(n_sources, dtype=float)

    for category in INVALID_SMILES_CATEGORY_ORDER:
        widths = plot_df[category].to_numpy(dtype=float) * 100.0 if n_sources else np.array([])
        ax.barh(
            y_positions,
            widths,
            left=left,
            height=0.72,
            color=INVALID_SMILES_CATEGORY_COLORS[category],
            edgecolor=BAR_EDGE,
            linewidth=1.0,
            label=INVALID_SMILES_CATEGORY_LABELS[category].replace("\n", " "),
        )
        left += widths

    ax.set_xlim(0, 100)
    ax.set_xlabel(INVALID_SMILES_NORMALIZE_LABELS[normalize_by], fontsize=11)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(plot_df["source"].tolist(), fontsize=10)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle="--", alpha=0.22, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.tick_params(axis="x", labelsize=10, width=1.0, length=4)
    ax.tick_params(axis="y", width=0, length=0)

    handles, labels = ax.get_legend_handles_labels()
    legend = fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.90),
        ncol=3,
        frameon=False,
        fontsize=9,
        handlelength=1.3,
        columnspacing=1.2,
    )
    for text in legend.get_texts():
        text.set_color(TEXT_DARK)

    if panel_label:
        fig.text(
            0.02,
            0.89,
            panel_label,
            fontsize=16,
            fontweight="bold",
            va="top",
            ha="left",
        )
    fig.suptitle(title, fontsize=12.5, y=0.98, x=0.5, ha="center")

    annotation_transform = transforms.blended_transform_factory(ax.transAxes, ax.transData)
    for index, row in plot_df.iterrows():
        ax.text(
            1.01,
            y_positions[index],
            f"{row['invalid_count']}/{row['total_count']} invalid ({row['invalid_fraction'] * 100:.1f}%)",
            transform=annotation_transform,
            fontsize=8.8,
            color=TEXT_MUTED,
            va="center",
            ha="left",
            clip_on=False,
        )

    fig.subplots_adjust(left=0.12, right=0.80, bottom=0.18, top=0.72)
    fig.savefig(figure_png_path, dpi=260, bbox_inches="tight")
    if figure_svg_path is not None:
        fig.savefig(figure_svg_path, bbox_inches="tight")
    plotting.plt.close(fig)


def choose_primary_representation_category(representation_tags: list[str]) -> str:
    for tag in representation_tags:
        if tag in REPRESENTATION_CATEGORY_ORDER:
            return tag
    return "other non-parent form"


def summarize_representation_groups(
    source_name: str,
    *,
    input_count: int,
    valid_count: int,
    invalid_count: int,
    groups: list[dict],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    issue_group_count = len(groups)
    issue_molecule_count = int(sum(group["molecule_count"] for group in groups))

    source_summary = {
        "source": source_name,
        "input_count": int(input_count),
        "valid_count": int(valid_count),
        "invalid_count": int(invalid_count),
        "issue_group_count": int(issue_group_count),
        "issue_molecule_count": int(issue_molecule_count),
        "issue_fraction_valid": (issue_molecule_count / valid_count) if valid_count else 0.0,
    }

    detail_rows = []
    for group_index, group in enumerate(groups, start=1):
        for variant_index, variant in enumerate(group["variants"], start=1):
            tags = list(variant.get("representation_tags") or [])
            primary_category = choose_primary_representation_category(tags)
            detail_rows.append(
                {
                    "source": source_name,
                    "group_index": group_index,
                    "variant_index": variant_index,
                    "group_parent_smiles": group["group_parent_smiles"],
                    "parent_form": group["parent_form"],
                    "group_tags": " | ".join(group.get("group_tags") or []),
                    "canonical_smiles": variant["canonical_smiles"],
                    "observed_parent_smiles": variant["observed_parent_smiles"],
                    "representation_tags": " | ".join(tags),
                    "primary_category": primary_category,
                    "count": int(variant["count"]),
                    "distinct_representations": int(group["distinct_representations"]),
                    "group_molecule_count": int(group["molecule_count"]),
                    "parent_form_backend_used": group["parent_form_backend_used"],
                    "parent_form_comment": group["parent_form_comment"],
                    "removed_salts_display": " | ".join(variant.get("removed_salts_display") or []),
                    "removed_solvents_display": " | ".join(variant.get("removed_solvents_display") or []),
                    "duplicate_parent_fragments_display": " | ".join(
                        variant.get("duplicate_parent_fragments_display") or []
                    ),
                    "multi_tag_variant": len(tags) > 1,
                }
            )

    return source_summary, detail_rows


def empty_representation_details_frame() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "source",
            "group_index",
            "variant_index",
            "group_parent_smiles",
            "parent_form",
            "group_tags",
            "canonical_smiles",
            "observed_parent_smiles",
            "representation_tags",
            "primary_category",
            "count",
            "distinct_representations",
            "group_molecule_count",
            "parent_form_backend_used",
            "parent_form_comment",
            "removed_salts_display",
            "removed_solvents_display",
            "duplicate_parent_fragments_display",
            "multi_tag_variant",
        ]
    )


def build_representation_category_summary(
    details_df: pd.DataFrame,
    source_summary_df: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for source_row in source_summary_df.itertuples(index=False):
        source_name = source_row.source
        source_details = details_df.loc[details_df["source"] == source_name]
        category_counts = (
            source_details.groupby("primary_category", sort=False)["count"].sum().to_dict()
            if not source_details.empty
            else {}
        )

        for category in REPRESENTATION_CATEGORY_ORDER:
            count = int(category_counts.get(category, 0))
            rows.append(
                {
                    "source": source_name,
                    "category": category,
                    "count": count,
                    "fraction_of_issue_only": (
                        count / source_row.issue_molecule_count if source_row.issue_molecule_count else 0.0
                    ),
                    "fraction_of_all": (count / source_row.valid_count) if source_row.valid_count else 0.0,
                }
            )

    category_df = pd.DataFrame(rows)
    if not category_df.empty:
        category_df["source"] = pd.Categorical(
            category_df["source"],
            categories=source_summary_df["source"].tolist(),
            ordered=True,
        )
        category_df["category"] = pd.Categorical(
            category_df["category"],
            categories=REPRESENTATION_CATEGORY_ORDER,
            ordered=True,
        )
        category_df = category_df.sort_values(["source", "category"]).reset_index(drop=True)
        category_df["source"] = category_df["source"].astype(str)
        category_df["category"] = category_df["category"].astype(str)
    return category_df


def build_representation_plot_table(
    source_summary_df: pd.DataFrame,
    category_summary_df: pd.DataFrame,
    normalize_by: str,
) -> pd.DataFrame:
    fraction_column = "fraction_of_issue_only" if normalize_by == "issue-only" else "fraction_of_all"
    rows = []

    for source_row in source_summary_df.itertuples(index=False):
        source_categories = category_summary_df.loc[
            category_summary_df["source"] == source_row.source
        ].set_index("category")
        row = {
            "source": source_row.source,
            "input_count": int(source_row.input_count),
            "valid_count": int(source_row.valid_count),
            "invalid_count": int(source_row.invalid_count),
            "issue_group_count": int(source_row.issue_group_count),
            "issue_molecule_count": int(source_row.issue_molecule_count),
            "issue_fraction_valid": float(source_row.issue_fraction_valid),
        }

        stack_total = 0.0
        for category in REPRESENTATION_CATEGORY_ORDER:
            value = float(source_categories.loc[category, fraction_column]) if category in source_categories.index else 0.0
            row[category] = value
            stack_total += value
        row["stack_total"] = stack_total
        rows.append(row)

    plot_df = pd.DataFrame(rows)
    if not plot_df.empty:
        plot_df["source"] = pd.Categorical(
            plot_df["source"],
            categories=source_summary_df["source"].tolist(),
            ordered=True,
        )
        plot_df = plot_df.sort_values("source").reset_index(drop=True)
        plot_df["source"] = plot_df["source"].astype(str)
    return plot_df


def render_representation_stacked_horizontal_figure(
    figure_png_path: Path | str,
    plot_df: pd.DataFrame,
    normalize_by: str,
    title: str,
    panel_label: str | None = None,
    figure_svg_path: Path | str | None = None,
):
    plotting.ensure_plotting_imports()
    from matplotlib import transforms

    n_sources = len(plot_df)
    fig_height = max(3.2, 1.4 + 0.75 * max(n_sources, 1))
    fig, ax = plotting.plt.subplots(figsize=(9.8, fig_height))

    y_positions = list(range(n_sources))
    left = [0.0 for _ in range(n_sources)]

    for category in REPRESENTATION_CATEGORY_ORDER:
        widths = plot_df[category].astype(float).mul(100.0).tolist() if n_sources else []
        ax.barh(
            y_positions,
            widths,
            left=left,
            height=0.72,
            color=REPRESENTATION_CATEGORY_COLORS[category],
            edgecolor=BAR_EDGE,
            linewidth=1.0,
            label=REPRESENTATION_CATEGORY_LABELS[category],
        )
        left = [base + width for base, width in zip(left, widths)]

    ax.set_xlim(0, 100)
    ax.set_xlabel(REPRESENTATION_NORMALIZE_LABELS[normalize_by], fontsize=11)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(plot_df["source"].tolist(), fontsize=10)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle="--", alpha=0.22, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.tick_params(axis="x", labelsize=10, width=1.0, length=4)
    ax.tick_params(axis="y", width=0, length=0)

    handles, labels = ax.get_legend_handles_labels()
    legend = fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.90),
        ncol=4,
        frameon=False,
        fontsize=9,
        handlelength=1.3,
        columnspacing=1.2,
    )
    for text in legend.get_texts():
        text.set_color(TEXT_DARK)

    if panel_label:
        fig.text(
            0.02,
            0.89,
            panel_label,
            fontsize=16,
            fontweight="bold",
            va="top",
            ha="left",
        )
    fig.suptitle(title, fontsize=12.5, y=0.98, x=0.5, ha="center")

    annotation_transform = transforms.blended_transform_factory(ax.transAxes, ax.transData)
    for index, row in plot_df.iterrows():
        ax.text(
            1.01,
            y_positions[index],
            (
                f"{row['issue_group_count']} groups; "
                f"{row['issue_molecule_count']}/{row['valid_count']} valid molecules "
                f"({row['issue_fraction_valid'] * 100:.1f}%)"
            ),
            transform=annotation_transform,
            fontsize=8.8,
            color=TEXT_MUTED,
            va="center",
            ha="left",
            clip_on=False,
        )

    fig.subplots_adjust(left=0.16, right=0.80, bottom=0.18, top=0.72)
    fig.savefig(figure_png_path, dpi=260, bbox_inches="tight")
    if figure_svg_path is not None:
        fig.savefig(figure_svg_path, bbox_inches="tight")
    plotting.plt.close(fig)
