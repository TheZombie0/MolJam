import io
import os

from .._common import *
from .. import plotting
from ..structural_breakdowns import (
    REPRESENTATION_NORMALIZE_FILENAMES,
    build_representation_category_summary,
    build_representation_plot_table,
    empty_representation_details_frame,
    render_representation_stacked_horizontal_figure,
    summarize_representation_groups,
)


class RepresentationConsistencyPlotMixin:
    _CATEGORY_DIRS = {
        'salt': 'salt',
        'acid adduct': 'acid_adduct',
        'solvent stripping': 'solvent_stripping',
        'protonated': 'protonated',
        'deprotonated': 'deprotonated',
        'duplicate-component': 'duplicate_component',
        'other non-parent form': 'other_non_parent_form',
    }

    @staticmethod
    def _draw_representation_image(smiles, width=600, height=420):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        from rdkit.Chem import AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
        from PIL import Image

        AllChem.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        options = drawer.drawOptions()
        options.baseFontSize = 0.9
        drawer.SetLineWidth(3)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return Image.open(io.BytesIO(drawer.GetDrawingText()))

    @staticmethod
    def _variant_rows(variant):
        example_rows = []
        for example in variant.get('examples', []):
            row_number = example.get('original_index')
            if row_number is not None:
                example_rows.append(f"#{row_number}")
        return ", ".join(example_rows) if example_rows else "N/A"

    @classmethod
    def _variant_summary_lines(cls, variant):
        lines = [
            f"Tags: {' + '.join(variant['representation_tags'])}",
            f"Count: {variant['count']}",
            f"Rows: {cls._variant_rows(variant)}",
        ]
        if variant.get('removed_salts_display'):
            lines.append("Removed salts: " + ", ".join(variant['removed_salts_display']))
        if variant.get('removed_solvents_display'):
            lines.append("Removed solvents: " + ", ".join(variant['removed_solvents_display']))
        if variant.get('duplicate_parent_fragments_display'):
            lines.append("Duplicate components: " + ", ".join(variant['duplicate_parent_fragments_display']))
        return "\n".join(lines)

    @staticmethod
    def _parent_form_summary_lines(group):
        lines = [
            f"Backend: {group['parent_form_backend_used']}",
            f"pH: {group['parent_form_ph']:.1f}",
        ]
        candidates = group.get('parent_form_candidates') or []
        if len(candidates) > 1:
            lines.append(f"Candidates: {len(candidates)}")
        return "\n".join(lines)

    def _build_representation_consistency_category_tables(self):
        source_summaries = []
        detail_rows = []

        for db_name, results in self.scoring_results.items():
            scorer = results["scorer"]
            if "check_representation_consistency" not in scorer.completed_checks:
                scorer.check_representation_consistency()

            source_summary, source_detail_rows = summarize_representation_groups(
                db_name,
                input_count=len(getattr(scorer, "df", [])),
                valid_count=len(getattr(scorer, "valid_mols", [])),
                invalid_count=len(getattr(scorer, "invalid_indices", [])),
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
        return details_df, source_summary_df, category_summary_df

    def plot_representation_consistency_categories(
        self,
        normalize_by="issue-only",
        save_path=None,
        save_svg_path=None,
        title="Representation Consistency categories",
        panel_label=None,
    ):
        """Plot representation-consistency categories across databases."""
        if not self.scoring_results:
            print("No scoring results available.")
            return

        _, source_summary_df, category_summary_df = self._build_representation_consistency_category_tables()
        plot_df = build_representation_plot_table(
            source_summary_df=source_summary_df,
            category_summary_df=category_summary_df,
            normalize_by=normalize_by,
        )

        if save_path is None:
            save_path = os.path.join(
                self.comparison_dir,
                f"representation_consistency_categories_{REPRESENTATION_NORMALIZE_FILENAMES[normalize_by]}.png",
            )

        save_dir = os.path.dirname(save_path)
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
        if save_svg_path:
            save_svg_dir = os.path.dirname(save_svg_path)
            if save_svg_dir:
                os.makedirs(save_svg_dir, exist_ok=True)

        render_representation_stacked_horizontal_figure(
            figure_png_path=save_path,
            figure_svg_path=save_svg_path,
            plot_df=plot_df,
            normalize_by=normalize_by,
            title=title,
            panel_label=panel_label,
        )
        print(f"Saved representation consistency categories plot to: {save_path}")

    def plot_representation_consistency_molecules(self, db_name, output_dir=None):
        """Save one image per representation-consistency group into categorized folders."""
        if db_name not in self.scoring_results:
            print(f"Database {db_name} not found in scoring results.")
            return

        plotting.ensure_plotting_imports()
        scorer = self.scoring_results[db_name]['scorer']

        if 'check_representation_consistency' not in scorer.completed_checks:
            scorer.check_representation_consistency()

        groups = getattr(scorer, 'representation_consistency_groups', [])
        if not groups:
            print(f"No representation consistency issues found in {db_name}")
            return

        if output_dir is None:
            output_dir = os.path.join(
                self.scoring_results[db_name]['output_dir'],
                'representation_consistency',
            )
        os.makedirs(output_dir, exist_ok=True)

        for group_idx, group in enumerate(groups, start=1):
            n_variants = len(group['variants'])
            n_columns = n_variants + 1
            fig = plotting.plt.figure(figsize=(max(10, n_columns * 3.3), 4.8))
            gs = fig.add_gridspec(
                1,
                n_columns,
                left=0.06,
                right=0.94,
                top=0.70,
                bottom=0.18,
                wspace=0.28,
            )

            parent_ax = fig.add_subplot(gs[0, 0])
            parent_image = self._draw_representation_image(group['parent_form'])
            if parent_image is None:
                parent_ax.text(
                    0.5,
                    0.5,
                    'Failed to parse parent form',
                    ha='center',
                    va='center',
                    transform=parent_ax.transAxes,
                )
                parent_ax.axis('off')
            else:
                parent_ax.imshow(parent_image)
                parent_ax.axis('off')

            parent_ax.set_title(
                "Parent form",
                fontsize=10,
                pad=8,
                bbox=dict(boxstyle="round", facecolor="#d9ecff", alpha=0.95),
            )
            parent_ax.text(
                0.5,
                -0.10,
                self._parent_form_summary_lines(group),
                ha='center',
                va='top',
                transform=parent_ax.transAxes,
                fontsize=8,
            )

            for col_idx, variant in enumerate(group['variants'], start=1):
                ax = fig.add_subplot(gs[0, col_idx])
                image = self._draw_representation_image(variant['canonical_smiles'])
                if image is None:
                    ax.text(
                        0.5,
                        0.5,
                        'Failed to parse molecule',
                        ha='center',
                        va='center',
                        transform=ax.transAxes,
                    )
                    ax.axis('off')
                else:
                    ax.imshow(image)
                    ax.axis('off')

                ax.set_title(
                    " + ".join(variant['representation_tags']),
                    fontsize=10,
                    pad=8,
                    bbox=dict(boxstyle="round", facecolor="#f2f2f2", alpha=0.9),
                )
                ax.text(
                    0.5,
                    -0.10,
                    self._variant_summary_lines(variant),
                    ha='center',
                    va='top',
                    transform=ax.transAxes,
                    fontsize=8,
                )

            fig.text(
                0.5,
                0.95,
                f"Database: {db_name}",
                ha='center',
                va='top',
                fontsize=12,
                bbox=dict(boxstyle="round", facecolor="#d9ecff", alpha=0.9),
            )
            fig.text(
                0.5,
                0.90,
                f"Parent form: {group['parent_form'][:120]}{'...' if len(group['parent_form']) > 120 else ''}",
                ha='center',
                va='top',
                fontsize=10,
            )
            fig.text(
                0.5,
                0.86,
                (
                    f"Distinct representations: {group['distinct_representations']} | "
                    f"Rows in group: {group['molecule_count']} | "
                    f"Categories: {', '.join(group['group_tags'])} | "
                    f"Grouping parent: {group['group_parent_smiles'][:80]}{'...' if len(group['group_parent_smiles']) > 80 else ''}"
                ),
                ha='center',
                va='top',
                fontsize=9,
            )

            image_buffer = io.BytesIO()
            plotting.plt.savefig(image_buffer, dpi=300, bbox_inches='tight')
            plotting.plt.close(fig)
            image_bytes = image_buffer.getvalue()
            image_buffer.close()

            file_name = f"group_{group_idx:03d}.png"
            category_tags = [
                tag for tag in group['group_tags']
                if tag in self._CATEGORY_DIRS and tag != 'parent form'
            ]
            for tag in category_tags:
                category_dir = os.path.join(output_dir, self._CATEGORY_DIRS[tag])
                os.makedirs(category_dir, exist_ok=True)
                save_path = os.path.join(category_dir, file_name)
                with open(save_path, 'wb') as handle:
                    handle.write(image_bytes)

        print(f"Saved {len(groups)} representation consistency group images to: {output_dir}")
