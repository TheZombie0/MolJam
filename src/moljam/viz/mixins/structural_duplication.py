from .._common import *
from .. import plotting


class StructuralDuplicationPlotMixin:
    def plot_structural_duplication_molecules(self, db_name, save_path=None):
        """Plot molecules with structural duplication - each row shows one duplicate group"""
        if db_name not in self.scoring_results:
            print(f"Database {db_name} not found in scoring results.")
            return

        scorer = self.scoring_results[db_name]['scorer']

        # Get structural duplication examples from analysis results
        if 'Data Consistency and Reliability' not in scorer.analysis_results:
            print(f"No data consistency analysis found for {db_name}")
            return

        structural_duplication = scorer.analysis_results['Data Consistency and Reliability'].get('Structural Duplication', {})
        duplicate_examples = structural_duplication.get('Duplicate SMILES examples', [])

        if not duplicate_examples:
            print(f"No structural duplication found in {db_name}")
            return

        # Fixed parameters for layout
        row_width = 4   # Fixed total width per row in inches
        row_height = 3  # Fixed total height per row in inches

        # Calculate total height based on number of groups
        n_groups = len(duplicate_examples)
        total_height = n_groups * row_height

        # Create single figure with all groups as rows
        fig = plotting.plt.figure(figsize=(row_width, total_height))

        # Create a list to track row positions and dimensions
        current_y_position = 0
        row_specs = []

        for group_idx, example in enumerate(duplicate_examples):
            canonical_smiles = example['canonical_smiles']
            original_smiles_with_indices = example['original_smiles_with_indices']
            count = example['count']

            # Get number of molecules in this group
            n_mols_in_group = len(original_smiles_with_indices)

            # Store this row's specification
            row_specs.append({
                'group_idx': group_idx,
                'n_mols': n_mols_in_group,
                'examples': original_smiles_with_indices,
                'canonical_smiles': canonical_smiles,
                'count': count,
                'y_start': current_y_position / total_height,
                'y_height': row_height / total_height
            })

            current_y_position += row_height

        # Now render each row with its own plotting.GridSpec
        for row_info in row_specs:
            group_idx = row_info['group_idx']
            n_mols_in_group = row_info['n_mols']
            y_start = row_info['y_start']
            y_height = row_info['y_height']

            # Create plotting.GridSpec for this row with proper spacing
            # top: leave space for row title (0.25 of row height)
            # bottom: leave space for subplot title (0.25 of row height)
            gs_row = fig.add_gridspec(1, n_mols_in_group,
                                     top=y_start + y_height,
                                     bottom=y_start)

            # Process each duplicate in this group (each becomes a column in this row)
            for col_idx, dup_info in enumerate(row_info['examples']):
                ax = fig.add_subplot(gs_row[0, col_idx])

                # Use the original SMILES for drawing
                smiles = dup_info['smiles']
                index = dup_info['index']  # Already 1-based from optimized_scorer.py
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

                    # Create drawer
                    drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
                    options = drawer.drawOptions()
                    options.baseFontSize = 6
                    drawer.SetLineWidth(4)

                    # Draw molecule
                    drawer.DrawMolecule(mol)
                    drawer.FinishDrawing()

                    # Convert to image
                    img_data = drawer.GetDrawingText()
                    img = Image.open(io.BytesIO(img_data))

                    # Display in subplot
                    ax.imshow(img)
                    ax.axis('off')

                # Add text information above the image (with reduced fontsize and padding)
                # title = f"Original: {smiles[:35]}{'...' if len(smiles) > 35 else ''}\n#{index}"
                # ax.set_title(title, fontsize=8) #, pad=2)

            # Add row title for this group (positioned above the plotting.GridSpec area)
            row_title_y = y_start + y_height
            display_group_idx = n_groups - group_idx  # Display in reverse order (1 to n_groups)
            row_title_text = f"Group {display_group_idx} (Count: {row_info['count']})\nCanonical: {row_info['canonical_smiles'][:40]}{'...' if len(row_info['canonical_smiles']) > 40 else ''}"

            # Create merged labels for all molecules in this row
            molecule_labels = []
            for mol_idx, dup_info in enumerate(row_info['examples']):
                smiles = dup_info['smiles']
                index = dup_info['index']
                molecule_labels.append(f"Molecule{mol_idx+1}:{smiles}   #{index}")
            merged_labels = "\n".join(molecule_labels)

            # Add main title
            fig.text(0.5, row_title_y, row_title_text,
                    ha='center', va='top', fontsize=8,
                    bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.7))

            # Add merged molecule labels below the main title
            labels_y = row_title_y - 0.01  # Position below the main title
            fig.text(0.5, labels_y, merged_labels,
                    ha='center', va='top', fontsize=8)
                    # bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

        # Add main title with proper spacing
        fig.suptitle(f'Structural Duplication in {db_name}', fontsize=16, y=1.01)

        # Use subplots_adjust to prevent overlap
        # plotting.plt.subplots_adjust(top=0.98, bottom=0.02, left=0.05, right=0.95, hspace=0.0)

        if save_path is None:
            save_path = os.path.join(self.scoring_results[db_name]['output_dir'],
                                   f'structural_duplication_molecules.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

    def plot_structural_duplication_molecules_individual(self, db_name, output_dir=None):
        """Plot each structural duplication group as a separate image file"""
        if db_name not in self.scoring_results:
            print(f"Database {db_name} not found in scoring results.")
            return

        scorer = self.scoring_results[db_name]['scorer']

        # Get structural duplication examples from analysis results
        if 'Data Consistency and Reliability' not in scorer.analysis_results:
            print(f"No data consistency analysis found for {db_name}")
            return

        structural_duplication = scorer.analysis_results['Data Consistency and Reliability'].get('Structural Duplication', {})
        duplicate_examples = structural_duplication.get('Duplicate SMILES examples', [])

        if not duplicate_examples:
            print(f"No structural duplication found in {db_name}")
            return

        # Create output directory
        if output_dir is None:
            output_dir = os.path.join(self.scoring_results[db_name]['output_dir'], 'structural_duplication_individual')
        os.makedirs(output_dir, exist_ok=True)

        print(f"Generating individual structural duplication images for {db_name}...")

        # Process each group
        for group_idx, example in enumerate(duplicate_examples):
            canonical_smiles = example['canonical_smiles']
            original_smiles_with_indices = example['original_smiles_with_indices']
            count = example['count']
            n_mols_in_group = len(original_smiles_with_indices)

            # Create figure for this group
            fig_width = n_mols_in_group * 2  # 2 inches per molecule
            fig_height = 3  # Fixed height
            fig = plotting.plt.figure(figsize=(fig_width, fig_height))

            # Create plotting.GridSpec for this row
            gs = fig.add_gridspec(1, n_mols_in_group, left=0.05, right=0.95, top=0.75, bottom=0.05)

            # Process each duplicate in this group
            for col_idx, dup_info in enumerate(original_smiles_with_indices):
                ax = fig.add_subplot(gs[0, col_idx])

                smiles = dup_info['smiles']
                index = dup_info['index']
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

                    drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
                    options = drawer.drawOptions()
                    options.baseFontSize = 6
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
                    bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.7))

            # Add molecule labels
            molecule_labels = []
            for mol_idx, dup_info in enumerate(original_smiles_with_indices):
                smiles = dup_info['smiles']
                index = dup_info['index']
                molecule_labels.append(f"Mol{mol_idx+1}: {smiles[:30]}{'...' if len(smiles) > 30 else ''}  #{index}")
            merged_labels = "\n".join(molecule_labels)

            fig.text(0.5, 0.85, merged_labels, ha='center', va='top', fontsize=7)

            # Save figure
            save_path = os.path.join(output_dir, f'group_{group_idx + 1}.png')
            plotting.plt.savefig(save_path, dpi=300, bbox_inches='tight')
            plotting.plt.close()

        print(f"  Saved {len(duplicate_examples)} structural duplication images to: {output_dir}")

