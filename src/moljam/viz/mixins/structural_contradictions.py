from .._common import *
from .. import plotting


class ContradictoryLabelPlotMixin:
    def plot_contradictory_label_molecules(self, db_name, save_path=None):
        """Plot molecules with contradictory labels - each label column gets a separate figure with multiple rows"""
        if db_name not in self.scoring_results:
            print(f"Database {db_name} not found in scoring results.")
            return

        scorer = self.scoring_results[db_name]['scorer']

        # Get contradictory examples from analysis results
        if 'Label Consistency' not in scorer.analysis_results:
            print(f"No label consistency analysis found for {db_name}")
            return

        label_consistency = scorer.analysis_results['Label Consistency']
        details_by_column = label_consistency.get('Details by column', {})

        if not details_by_column:
            print(f"No label consistency details found in {db_name}")
            return

        # Process each label column separately (creates one figure per column)
        for column_name, column_data in details_by_column.items():
            contradictory_examples = column_data.get('Contradictory examples', [])

            if not contradictory_examples:
                print(f"No contradictory examples found for column {column_name} in {db_name}")
                continue

            # Fixed parameters for layout
            row_width = 4   # Fixed total width per row in inches
            row_height = 3  # Fixed total height per row in inches
            max_rows_per_figure = 20
            n_molecules_total = len(contradictory_examples)
            column_safe = column_name.replace('/', '_').replace(' ', '_')

            for page_idx, start_idx in enumerate(range(0, n_molecules_total, max_rows_per_figure), start=1):
                page_examples = contradictory_examples[start_idx:start_idx + max_rows_per_figure]
                n_rows = len(page_examples)
                total_height = n_rows * row_height

                fig = plotting.plt.figure(figsize=(row_width, total_height))
                current_y_position = 0
                row_specs = []

                for local_idx, example in enumerate(page_examples):
                    original_smiles_with_indices = example['original_smiles_with_indices']
                    n_instances_in_group = len(original_smiles_with_indices)

                    row_specs.append({
                        'global_idx': start_idx + local_idx,
                        'n_instances': n_instances_in_group,
                        'instances': original_smiles_with_indices,
                        'canonical_smiles': example['canonical_smiles'],
                        'conflicting_labels': example['conflicting_labels'],
                        'count': example['count'],
                        'y_start': current_y_position / total_height,
                        'y_height': row_height / total_height
                    })

                    current_y_position += row_height

                for row_info in row_specs:
                    n_instances_in_group = row_info['n_instances']
                    y_start = row_info['y_start']
                    y_height = row_info['y_height']

                    gs_row = fig.add_gridspec(1, n_instances_in_group,
                                             top=y_start + y_height,
                                             bottom=y_start)

                    for col_idx, instance in enumerate(row_info['instances']):
                        ax = fig.add_subplot(gs_row[0, col_idx])

                        smiles = instance['smiles']
                        index = instance['index']
                        label_value = instance['label']
                        mol = Chem.MolFromSmiles(smiles)

                        if mol is None:
                            ax.text(0.5, 0.5, 'Failed to parse molecule',
                                   ha='center', va='center', transform=ax.transAxes, fontsize=8)
                            ax.axis('off')
                        else:
                            from rdkit.Chem import AllChem
                            from rdkit.Chem.Draw import rdMolDraw2D
                            from PIL import Image
                            import io

                            AllChem.Compute2DCoords(mol)
                            drawer = rdMolDraw2D.MolDraw2DCairo(400, 400)
                            options = drawer.drawOptions()
                            options.baseFontSize = 2
                            drawer.SetLineWidth(4)
                            drawer.DrawMolecule(mol)
                            drawer.FinishDrawing()

                            img_data = drawer.GetDrawingText()
                            img = Image.open(io.BytesIO(img_data))
                            ax.imshow(img)
                            ax.axis('off')

                    row_title_y = y_start + y_height
                    display_mol_idx = n_molecules_total - row_info['global_idx']
                    row_title_text = f"Molecule {display_mol_idx} (Count: {row_info['count']})\nCanonical: {row_info['canonical_smiles'][:40]}{'...' if len(row_info['canonical_smiles']) > 40 else ''}"

                    molecule_labels = []
                    for mol_idx_inner, instance in enumerate(row_info['instances']):
                        smiles = instance['smiles']
                        index = instance['index']
                        label_value = instance['label']
                        molecule_labels.append(f"Molecule{mol_idx_inner+1}: {smiles}   Label: {label_value}   #{index}")
                    merged_labels = "\n".join(molecule_labels)

                    fig.text(0.5, row_title_y, row_title_text,
                            ha='center', va='top', fontsize=8,
                            bbox=dict(boxstyle="round", facecolor="lightcoral", alpha=0.7))

                    labels_y = row_title_y - 0.011
                    fig.text(0.5, labels_y, merged_labels,
                            ha='center', va='top', fontsize=8)

                title = f'Contradictory Labels in {db_name} - Column: {column_name}'
                if n_molecules_total > max_rows_per_figure:
                    title += f' (Part {page_idx})'
                fig.suptitle(title, fontsize=13, y=1.01)

                if save_path:
                    base_path = save_path.replace('.png', '')
                    column_save_path = f"{base_path}_{column_safe}"
                else:
                    column_save_path = os.path.join(
                        self.scoring_results[db_name]['output_dir'],
                        f'contradictory_labels_{column_safe}'
                    )

                if n_molecules_total > max_rows_per_figure:
                    column_save_path = f"{column_save_path}_part{page_idx}"

                column_save_path = f"{column_save_path}.png"
                plotting.plt.savefig(column_save_path, dpi=500, bbox_inches='tight')
                print(f"Saved contradictory labels plot for column '{column_name}' to: {column_save_path}")
                plotting.plt.close()

        #

    def plot_contradictory_label_molecules_individual(self, db_name, output_dir=None):
        """Plot each contradictory label molecule as a separate image file"""
        if db_name not in self.scoring_results:
            print(f"Database {db_name} not found in scoring results.")
            return

        scorer = self.scoring_results[db_name]['scorer']

        # Get contradictory examples from analysis results
        if 'Label Consistency' not in scorer.analysis_results:
            print(f"No label consistency analysis found for {db_name}")
            return

        label_consistency = scorer.analysis_results['Label Consistency']
        details_by_column = label_consistency.get('Details by column', {})

        if not details_by_column:
            print(f"No label consistency details found in {db_name}")
            return

        print(f"Generating individual contradictory label images for {db_name}...")

        # Process each label column separately
        total_images = 0
        for column_name, column_data in details_by_column.items():
            contradictory_examples = column_data.get('Contradictory examples', [])

            if not contradictory_examples:
                continue

            # Create output directory for this column
            column_safe = column_name.replace('/', '_').replace(' ', '_')
            if output_dir is None:
                col_output_dir = os.path.join(
                    self.scoring_results[db_name]['output_dir'],
                    'contradictory_labels_individual',
                    column_safe
                )
            else:
                col_output_dir = os.path.join(output_dir, column_safe)

            os.makedirs(col_output_dir, exist_ok=True)

            # Process each molecule in this column
            for mol_idx, example in enumerate(contradictory_examples):
                canonical_smiles = example['canonical_smiles']
                original_smiles_with_indices = example['original_smiles_with_indices']
                conflicting_labels = example['conflicting_labels']
                count = example['count']
                n_instances_in_group = len(original_smiles_with_indices)

                # Create figure for this molecule
                fig_width = n_instances_in_group * 2  # 2 inches per instance
                fig_height = 3  # Fixed height
                fig = plotting.plt.figure(figsize=(fig_width, fig_height))

                # Create plotting.GridSpec for this row
                gs = fig.add_gridspec(1, n_instances_in_group, left=0.05, right=0.95, top=0.75, bottom=0.05)

                # Process each instance
                for col_idx, instance in enumerate(original_smiles_with_indices):
                    ax = fig.add_subplot(gs[0, col_idx])

                    smiles = instance['smiles']
                    index = instance['index']
                    label_value = instance['label']
                    mol = Chem.MolFromSmiles(smiles)

                    if mol is None:
                        ax.text(0.5, 0.5, 'Failed to parse molecule',
                               ha='center', va='center', transform=ax.transAxes, fontsize=8)
                        ax.axis('off')
                    else:
                        # Generate 2D coordinates
                        from rdkit.Chem import AllChem
                        AllChem.Compute2DCoords(mol)

                        # Draw molecule
                        from rdkit.Chem.Draw import rdMolDraw2D
                        from PIL import Image
                        import io

                        drawer = rdMolDraw2D.MolDraw2DCairo(400, 400)
                        options = drawer.drawOptions()
                        options.baseFontSize = 2
                        drawer.SetLineWidth(4)

                        drawer.DrawMolecule(mol)
                        drawer.FinishDrawing()

                        img_data = drawer.GetDrawingText()
                        img = Image.open(io.BytesIO(img_data))

                        ax.imshow(img)
                        ax.axis('off')

                # Add title and labels at the top
                title_text = f"Canonical: {canonical_smiles[:60]}{'...' if len(canonical_smiles) > 60 else ''}"
                fig.text(0.5, 0.92, title_text, ha='center', va='top', fontsize=9,
                        bbox=dict(boxstyle="round", facecolor="lightcoral", alpha=0.7))

                # Add molecule labels
                molecule_labels = []
                for inst_idx, instance in enumerate(original_smiles_with_indices):
                    smiles = instance['smiles']
                    index = instance['index']
                    label_value = instance['label']
                    molecule_labels.append(f"Mol{inst_idx+1}: {smiles[:25]}{'...' if len(smiles) > 25 else ''}  Label: {label_value}  #{index}")
                merged_labels = "\n".join(molecule_labels)

                fig.text(0.5, 0.85, merged_labels, ha='center', va='top', fontsize=7)

                # Save figure
                save_path = os.path.join(col_output_dir, f'molecule_{mol_idx + 1}.png')
                plotting.plt.savefig(save_path, dpi=300, bbox_inches='tight')
                plotting.plt.close()
                total_images += 1

            print(f"  Saved {len(contradictory_examples)} images for column '{column_name}' to: {col_output_dir}")

        print(f"  Total: {total_images} contradictory label images generated")
