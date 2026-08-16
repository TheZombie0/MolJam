import math

from .._common import *
from .. import plotting


class ChiralityPlotMixin:
    _STEREOCHEMISTRY_COMPARISON_FIGSIZE = (14, 8)
    _STEREOCHEMISTRY_AXIS_LABEL_FONTSIZE = 21
    _STEREOCHEMISTRY_SUBTITLE_FONTSIZE = 22
    _STEREOCHEMISTRY_SUPTITLE_FONTSIZE = 23
    _STEREOCHEMISTRY_TICK_FONTSIZE = 18
    _STEREOCHEMISTRY_ANNOTATION_FONTSIZE = 9
    _STEREOCHEMISTRY_RATIO_ANNOTATION_THRESHOLD = 0.05

    _STEREOCHEMISTRY_PLOT_CONFIGS = {
        'chirality': {
            'title': 'Undefined Chirality',
            'db_title': 'Undefined Chirality Distribution in {db_name}',
            'comparison_title': 'Undefined Chirality Distribution Across Databases',
            'exact_distribution_key': 'Undefined chirality exact count distribution',
            'legacy_distribution_key': 'Undefined chirality count distribution',
            'total_key': 'Molecules with undefined chirality',
            'filename_stem': 'undefined_chirality_distribution',
            'x_label': 'Undefined chiral centers per molecule',
            'count_ylabel': 'Number of molecules',
            'ratio_ylabel': 'Molecules / database (%)',
            'empty_message': 'No molecules with undefined chirality found in {db_name}',
            'empty_comparison_message': 'No molecules with undefined chirality were found across the selected databases.',
        },
        'double_bond': {
            'title': 'Undefined Double-Bond E/Z Stereochemistry',
            'db_title': 'Undefined Double-Bond E/Z Distribution in {db_name}',
            'comparison_title': 'Undefined Double-Bond E/Z Distribution Across Databases',
            'exact_distribution_key': 'Undefined double bond exact count distribution',
            'legacy_distribution_key': 'Undefined double bond count distribution',
            'total_key': 'Molecules with undefined double bond stereochemistry',
            'filename_stem': 'undefined_double_bond_distribution',
            'x_label': 'Undefined stereogenic double bonds per molecule',
            'count_ylabel': 'Number of molecules',
            'ratio_ylabel': 'Molecules / database (%)',
            'empty_message': 'No molecules with undefined double-bond E/Z stereochemistry found in {db_name}',
            'empty_comparison_message': 'No molecules with undefined double-bond E/Z stereochemistry were found across the selected databases.',
        },
    }

    def plot_undefined_chirality_molecules(self, db_name, figsize=(16, 16), save_path=None):
        """Plot molecules with undefined chirality centers."""
        plotting.ensure_plotting_imports()

        if db_name not in self.scoring_results:
            print(f"Database {db_name} not found in scoring results.")
            return

        scorer = self.scoring_results[db_name]['scorer']

        if 'Stereochemistry Completeness' not in scorer.analysis_results:
            print(f"No stereochemistry analysis found for {db_name}")
            return

        examples = scorer.analysis_results['Stereochemistry Completeness'].get(
            'Example molecules with undefined chirality',
            [],
        )

        if not examples:
            print(f"No molecules with undefined chirality found in {db_name}")
            return

        n_examples = len(examples)
        n_cols = min(3, n_examples)
        n_rows = (n_examples + n_cols - 1) // n_cols

        from matplotlib.lines import Line2D

        fig = plotting.plt.figure(figsize=figsize)
        gs = plotting.GridSpec(n_rows, n_cols, figure=fig, hspace=0.4, wspace=0.3)

        for idx, mol_info in enumerate(examples):
            row = idx // n_cols
            col = idx % n_cols
            ax = fig.add_subplot(gs[row, col])

            smiles = mol_info.get('analysis_smiles') or mol_info['canonical_smiles']
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                ax.text(0.5, 0.5, 'Failed to parse molecule', ha='center', va='center', transform=ax.transAxes)
                ax.axis('off')
                continue

            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            undefined_atoms = [atom_idx for atom_idx, flag in chiral_centers if flag == '?']

            from rdkit.Chem import AllChem
            from rdkit.Chem.Draw import rdMolDraw2D
            from PIL import Image
            import io

            AllChem.Compute2DCoords(mol)

            drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
            options = drawer.drawOptions()
            options.baseFontSize = 3
            drawer.SetLineWidth(4)

            highlight_atom_colors = {atom: (1.0, 0.0, 0.0) for atom in undefined_atoms}
            highlight_atom_radii = {atom: 0.8 for atom in undefined_atoms}

            drawer.DrawMolecule(
                mol,
                highlightAtoms=undefined_atoms,
                highlightAtomColors=highlight_atom_colors,
                highlightAtomRadii=highlight_atom_radii,
            )
            drawer.FinishDrawing()

            img = Image.open(io.BytesIO(drawer.GetDrawingText()))
            ax.imshow(img)
            ax.axis('off')

            title = f"Molecule {idx + 1}\n"
            title += f"Original: {mol_info['smiles'][:30]}{'...' if len(mol_info['smiles']) > 30 else ''}\n"
            title += (
                f"Canonical: {mol_info['canonical_smiles'][:30]}"
                f"{'...' if len(mol_info['canonical_smiles']) > 30 else ''}\n"
            )
            if mol_info.get('analysis_smiles') and mol_info['analysis_smiles'] != mol_info['canonical_smiles']:
                title += (
                    f"Count basis: {mol_info['analysis_smiles'][:30]}"
                    f"{'...' if len(mol_info['analysis_smiles']) > 30 else ''}\n"
                )
            title += f"Undefined centers: {mol_info['undefined_centers']}/{mol_info['total_centers']}"
            ax.set_title(title, fontsize=12, pad=10)

        fig.suptitle(f'Molecules with Undefined Chirality in {db_name}', fontsize=23)

        legend_elements = [
            Line2D(
                [0],
                [0],
                marker='o',
                color='w',
                markerfacecolor='r',
                markersize=10,
                label='Undefined chiral center',
            )
        ]
        fig.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, 0.2), ncol=1, fontsize=16)

        plotting.plt.tight_layout()

        if save_path:
            plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
            plotting.plt.close(fig)

    def _get_stereochemistry_plot_config(self, feature):
        if feature not in self._STEREOCHEMISTRY_PLOT_CONFIGS:
            raise ValueError(f"Unsupported stereochemistry feature: {feature}")
        return self._STEREOCHEMISTRY_PLOT_CONFIGS[feature]

    def _coerce_exact_distribution(self, exact_distribution, legacy_distribution=None):
        source_distribution = exact_distribution or legacy_distribution or {}
        normalized = {}

        for raw_key, raw_value in source_distribution.items():
            if raw_value in (None, ''):
                continue

            if isinstance(raw_key, (int, np.integer)):
                count = int(raw_key)
            else:
                digits = ''.join(ch for ch in str(raw_key) if ch.isdigit())
                if not digits:
                    continue
                count = int(digits)

            if count <= 0:
                continue

            normalized[count] = normalized.get(count, 0) + int(raw_value)

        return dict(sorted(normalized.items()))

    def _extract_stereochemistry_distribution(self, db_name, feature):
        if db_name not in self.scoring_results:
            raise KeyError(f"Database {db_name} not found in scoring results.")

        config = self._get_stereochemistry_plot_config(feature)
        scorer = self.scoring_results[db_name]['scorer']
        stereo_data = scorer.analysis_results.get('Stereochemistry Completeness', {})

        distribution = self._coerce_exact_distribution(
            stereo_data.get(config['exact_distribution_key'], {}),
            stereo_data.get(config['legacy_distribution_key'], {}),
        )
        total_undefined = int(stereo_data.get(config['total_key'], sum(distribution.values())))
        total_molecules = int(getattr(scorer, 'num_molecules', len(getattr(scorer, 'df', []))))
        ratio = (total_undefined / total_molecules * 100.0) if total_molecules > 0 else 0.0

        return {
            'db_name': db_name,
            'distribution': distribution,
            'total_undefined': total_undefined,
            'total_molecules': total_molecules,
            'ratio': ratio,
        }

    def _build_stereochemistry_distribution_frames(self, feature):
        payloads = [
            self._extract_stereochemistry_distribution(db_name, feature)
            for db_name in self.scoring_results
        ]

        all_counts = sorted(
            {
                count
                for payload in payloads
                for count in payload['distribution']
            }
        )

        count_df = pd.DataFrame(
            0,
            index=[payload['db_name'] for payload in payloads],
            columns=all_counts,
            dtype=float,
        )

        for payload in payloads:
            for count, value in payload['distribution'].items():
                count_df.loc[payload['db_name'], count] = int(value)

        total_molecules = pd.Series(
            {payload['db_name']: payload['total_molecules'] for payload in payloads},
            dtype=float,
        )
        total_undefined = pd.Series(
            {payload['db_name']: payload['total_undefined'] for payload in payloads},
            dtype=float,
        )
        ratio_df = count_df.div(total_molecules.replace(0, np.nan), axis=0).fillna(0.0) * 100.0

        return payloads, count_df, ratio_df, total_molecules, total_undefined

    def _build_stereochemistry_comparison_tick_labels(self, columns, max_labels=20):
        if len(columns) == 0:
            return []

        label_step = max(1, math.ceil(len(columns) / max_labels))
        return [
            str(column) if idx % label_step == 0 or idx == len(columns) - 1 else ''
            for idx, column in enumerate(columns)
        ]

    def _build_stereochemistry_annotations(self, data_frame, value_kind):
        annotations = data_frame.copy().astype(object)

        for row_idx in range(data_frame.shape[0]):
            for col_idx in range(data_frame.shape[1]):
                value = float(data_frame.iat[row_idx, col_idx])
                if value_kind == 'count':
                    annotations.iat[row_idx, col_idx] = '' if value == 0.0 else str(int(value))
                else:
                    annotations.iat[row_idx, col_idx] = (
                        ''
                        if value < self._STEREOCHEMISTRY_RATIO_ANNOTATION_THRESHOLD
                        else f"{value:.2f}%"
                    )

        return annotations

    def _resolve_stereochemistry_comparison_save_path(self, config, mode, use_log=False, save_path=None):
        if save_path is not None:
            return save_path

        suffix = f"_{mode}"
        if mode == 'count' and use_log:
            suffix += "_log"
        return os.path.join(self.comparison_dir, f"{config['filename_stem']}{suffix}.png")

    def _expand_stereochemistry_comparison_save_paths(self, config, use_log=False, save_path=None):
        if save_path is None:
            return {
                'count': self._resolve_stereochemistry_comparison_save_path(
                    config,
                    mode='count',
                    use_log=use_log,
                    save_path=None,
                ),
                'ratio': self._resolve_stereochemistry_comparison_save_path(
                    config,
                    mode='ratio',
                    use_log=False,
                    save_path=None,
                ),
            }

        save_root, save_ext = os.path.splitext(save_path)
        save_ext = save_ext or '.png'
        return {
            'count': f"{save_root}_count{'_log' if use_log else ''}{save_ext}",
            'ratio': f"{save_root}_ratio{save_ext}",
        }

    def _render_stereochemistry_comparison_heatmap(
        self,
        data_frame,
        annotations,
        config,
        subtitle,
        cmap,
        cbar_label,
        save_path,
        figsize=None,
        use_log=False,
    ):
        figure_size = figsize or self._STEREOCHEMISTRY_COMPARISON_FIGSIZE
        plot_frame = np.log10(data_frame + 1.0) if use_log else data_frame

        fig, ax = plotting.plt.subplots(figsize=figure_size)
        heatmap = plotting.sns.heatmap(
            plot_frame,
            annot=annotations,
            fmt='',
            cmap=cmap,
            linewidths=0.4,
            linecolor='white',
            cbar_kws={'label': cbar_label},
            annot_kws={'fontsize': self._STEREOCHEMISTRY_ANNOTATION_FONTSIZE},
            ax=ax,
        )

        colorbar = heatmap.collections[0].colorbar
        colorbar.ax.tick_params(labelsize=self._STEREOCHEMISTRY_TICK_FONTSIZE)
        colorbar.set_label(cbar_label, fontsize=self._STEREOCHEMISTRY_AXIS_LABEL_FONTSIZE)

        fig.suptitle(config['comparison_title'], fontsize=self._STEREOCHEMISTRY_SUPTITLE_FONTSIZE, y=0.98)
        ax.set_title(subtitle, fontsize=self._STEREOCHEMISTRY_SUBTITLE_FONTSIZE, pad=10)
        ax.set_xlabel(config['x_label'], fontsize=self._STEREOCHEMISTRY_AXIS_LABEL_FONTSIZE)
        ax.set_ylabel('Database', fontsize=self._STEREOCHEMISTRY_AXIS_LABEL_FONTSIZE)

        ax.set_yticklabels([str(label) for label in data_frame.index], rotation=0)
        ax.tick_params(axis='y', labelsize=self._STEREOCHEMISTRY_TICK_FONTSIZE)
        ax.tick_params(axis='x', labelsize=self._STEREOCHEMISTRY_TICK_FONTSIZE, rotation=45)
        ax.set_xticklabels(
            self._build_stereochemistry_comparison_tick_labels(data_frame.columns),
            rotation=45,
            ha='right',
        )

        fig.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close(fig)

    def _build_distribution_palette(self, n_colors):
        if n_colors <= 0:
            return []
        return plotting.sns.color_palette('crest', n_colors=max(n_colors, 3))[-n_colors:]

    def _annotate_bar_values(self, ax, bars, labels):
        max_height = max((bar.get_height() for bar in bars), default=0.0)
        offset = max(max_height * 0.02, 0.05)

        for bar, label in zip(bars, labels):
            height = bar.get_height()
            if height <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                height + offset,
                label,
                ha='center',
                va='bottom',
                fontsize=9,
            )

    def _render_empty_distribution_figure(self, title, message, save_path, figsize=None):
        figure_size = figsize or (8, 4.5)
        fig, ax = plotting.plt.subplots(figsize=figure_size)
        ax.axis('off')
        ax.text(0.5, 0.5, message, ha='center', va='center', fontsize=12, transform=ax.transAxes)
        fig.suptitle(title, fontsize=16)
        fig.tight_layout()
        fig.savefig(save_path, dpi=350, bbox_inches='tight')
        plotting.plt.close(fig)

    def _plot_stereochemistry_distribution(self, db_name, feature, figsize=None, save_path=None):
        plotting.ensure_plotting_imports()

        try:
            payload = self._extract_stereochemistry_distribution(db_name, feature)
        except KeyError:
            print(f"Database {db_name} not found in scoring results.")
            return

        config = self._get_stereochemistry_plot_config(feature)
        if save_path is None:
            output_dir = self.scoring_results[db_name].get('output_dir', self.comparison_dir)
            save_path = os.path.join(output_dir, f"{config['filename_stem']}.png")

        if not payload['distribution']:
            self._render_empty_distribution_figure(
                config['db_title'].format(db_name=db_name),
                config['empty_message'].format(db_name=db_name),
                save_path,
                figsize=figsize,
            )
            return

        counts = list(payload['distribution'].keys())
        molecule_counts = np.array([payload['distribution'][count] for count in counts], dtype=float)
        molecule_ratios = (
            molecule_counts / payload['total_molecules'] * 100.0
            if payload['total_molecules'] > 0
            else np.zeros_like(molecule_counts)
        )

        if figsize is None:
            figsize = (max(7.5, 4.2 + 0.8 * len(counts)), 8.2)

        fig, (ax_count, ax_ratio) = plotting.plt.subplots(2, 1, figsize=figsize, sharex=True)
        x = np.arange(len(counts))
        palette = self._build_distribution_palette(len(counts))

        count_bars = ax_count.bar(x, molecule_counts, color=palette, edgecolor='white', linewidth=1.1)
        ratio_bars = ax_ratio.bar(x, molecule_ratios, color=palette, edgecolor='white', linewidth=1.1)

        self._annotate_bar_values(ax_count, count_bars, [str(int(value)) for value in molecule_counts])
        self._annotate_bar_values(ax_ratio, ratio_bars, [f"{value:.2f}%" for value in molecule_ratios])

        for ax, ylabel in ((ax_count, config['count_ylabel']), (ax_ratio, config['ratio_ylabel'])):
            ax.set_xticks(x)
            ax.set_xticklabels(counts)
            ax.set_xlabel(config['x_label'], fontsize=11)
            ax.set_ylabel(ylabel, fontsize=11)
            ax.grid(axis='y', linestyle='--', alpha=0.25)
            ax.set_axisbelow(True)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        ax_count.set_title('Count', fontsize=13)
        ax_ratio.set_title('Percentage of database', fontsize=13)

        if np.any(molecule_counts > 0):
            ax_count.set_ylim(0, molecule_counts.max() * 1.18 + 0.5)
        if np.any(molecule_ratios > 0):
            ax_ratio.set_ylim(0, molecule_ratios.max() * 1.18 + 0.15)

        fig.suptitle(config['db_title'].format(db_name=db_name), fontsize=16)
        fig.text(
            0.5,
            0.015,
            f"{payload['total_undefined']} molecules affected "
            f"({payload['ratio']:.2f}% of {payload['total_molecules']} total molecules)",
            ha='center',
            va='bottom',
            fontsize=10,
        )
        fig.tight_layout(rect=[0, 0.04, 1, 0.95], h_pad=1.5)
        fig.savefig(save_path, dpi=350, bbox_inches='tight')
        plotting.plt.close(fig)

    def _plot_stereochemistry_comparison(self, feature, mode='both', use_log=False, add_inset=False,
                                         figsize=None, save_path=None):
        plotting.ensure_plotting_imports()

        if not self.scoring_results:
            print("No scoring results available.")
            return

        if mode not in {'both', 'count', 'ratio'}:
            raise ValueError("mode must be one of: 'both', 'count', 'ratio'")

        config = self._get_stereochemistry_plot_config(feature)
        _, count_df, ratio_df, _, _ = self._build_stereochemistry_distribution_frames(feature)

        if mode == 'both':
            save_paths = self._expand_stereochemistry_comparison_save_paths(
                config,
                use_log=use_log,
                save_path=save_path,
            )
            self._plot_stereochemistry_comparison(
                feature,
                mode='count',
                use_log=use_log,
                add_inset=add_inset,
                figsize=figsize,
                save_path=save_paths['count'],
            )
            self._plot_stereochemistry_comparison(
                feature,
                mode='ratio',
                use_log=False,
                add_inset=add_inset,
                figsize=figsize,
                save_path=save_paths['ratio'],
            )
            return

        save_path = self._resolve_stereochemistry_comparison_save_path(
            config,
            mode=mode,
            use_log=use_log,
            save_path=save_path,
        )

        if count_df.empty or count_df.shape[1] == 0 or float(count_df.to_numpy().sum()) == 0.0:
            self._render_empty_distribution_figure(
                config['comparison_title'],
                config['empty_comparison_message'],
                save_path,
                figsize=figsize or self._STEREOCHEMISTRY_COMPARISON_FIGSIZE,
            )
            return

        if mode == 'count':
            self._render_stereochemistry_comparison_heatmap(
                data_frame=count_df,
                annotations=self._build_stereochemistry_annotations(count_df, value_kind='count'),
                config=config,
                subtitle='Count (log color scale)' if use_log else 'Count',
                cmap='YlOrRd',
                cbar_label='log10(molecule count + 1)' if use_log else 'Number of molecules',
                save_path=save_path,
                figsize=figsize,
                use_log=use_log,
            )
            return

        self._render_stereochemistry_comparison_heatmap(
            data_frame=ratio_df,
            annotations=self._build_stereochemistry_annotations(ratio_df, value_kind='ratio'),
            config=config,
            subtitle='Percentage of database',
            cmap='YlGnBu',
            cbar_label='Molecules / database (%)',
            save_path=save_path,
            figsize=figsize,
            use_log=False,
        )

    def plot_undefined_chirality_distribution(self, db_name, figsize=None, save_path=None):
        """Plot the exact undefined-chirality distribution for one database."""
        self._plot_stereochemistry_distribution(db_name, 'chirality', figsize=figsize, save_path=save_path)

    def plot_undefined_double_bond_distribution(self, db_name, figsize=None, save_path=None):
        """Plot the exact undefined-double-bond distribution for one database."""
        self._plot_stereochemistry_distribution(db_name, 'double_bond', figsize=figsize, save_path=save_path)

    def plot_undefined_chirality_comparison(self, mode='both', use_log=False, add_inset=False,
                                            figsize=None, save_path=None):
        """
        Plot the exact undefined-chirality distribution across databases.

        Args:
            mode: 'both', 'count', or 'ratio'. 'both' saves separate count and ratio figures.
            use_log: Apply a log-color transform to the count heatmap.
            add_inset: Kept for backward compatibility; ignored for matrix plots.
        """
        self._plot_stereochemistry_comparison(
            'chirality',
            mode=mode,
            use_log=use_log,
            add_inset=add_inset,
            figsize=figsize,
            save_path=save_path,
        )

    def plot_undefined_double_bond_comparison(self, mode='both', use_log=False, add_inset=False,
                                              figsize=None, save_path=None):
        """
        Plot the exact undefined-double-bond distribution across databases.

        Args:
            mode: 'both', 'count', or 'ratio'. 'both' saves separate count and ratio figures.
            use_log: Apply a log-color transform to the count heatmap.
            add_inset: Kept for backward compatibility; ignored for matrix plots.
        """
        self._plot_stereochemistry_comparison(
            'double_bond',
            mode=mode,
            use_log=use_log,
            add_inset=add_inset,
            figsize=figsize,
            save_path=save_path,
        )
