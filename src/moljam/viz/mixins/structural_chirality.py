from .._common import *
from .. import plotting


class ChiralityPlotMixin:
    def plot_undefined_chirality_molecules(self, db_name, figsize=(16, 16), save_path=None):
        """Plot molecules with undefined chirality centers"""
        if db_name not in self.scoring_results:
            print(f"Database {db_name} not found in scoring results.")
            return
            
        snap = self.scoring_results[db_name]['snapshot']

        # Get undefined chirality examples from analysis results
        if 'Stereochemistry Completeness' not in snap.analysis_results:
            print(f"No stereochemistry analysis found for {db_name}")
            return
            
        examples = snap.analysis_results['Stereochemistry Completeness'].get('Example molecules with undefined chirality', [])
        
        if not examples:
            print(f"No molecules with undefined chirality found in {db_name}")
            return
            
        # Create figure with subplots
        n_examples = len(examples)
        n_cols = min(3, n_examples)
        n_rows = (n_examples + n_cols - 1) // n_cols
        
        fig = plotting.plt.figure(figsize=figsize)
        gs = plotting.GridSpec(n_rows, n_cols, figure=fig, hspace=0.4, wspace=0.3)
        
        for idx, mol_info in enumerate(examples):
            row = idx // n_cols
            col = idx % n_cols
            ax = fig.add_subplot(gs[row, col])
            
            # Parse SMILES and create molecule
            smiles = mol_info['canonical_smiles']
            mol = Chem.MolFromSmiles(smiles)
            
            if mol is None:
                ax.text(0.5, 0.5, 'Failed to parse molecule', 
                       ha='right', va='center', transform=ax.transAxes)
                ax.axis('off')
                continue
            
            # Find undefined chiral centers
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            undefined_atoms = [atom_idx for atom_idx, flag in chiral_centers if flag == '?']
            
            # Generate 2D coordinates
            from rdkit.Chem import AllChem
            AllChem.Compute2DCoords(mol)
            
            # Draw molecule with highlighted undefined chiral centers
            from rdkit.Chem.Draw import rdMolDraw2D
            from PIL import Image
            import io
            
            # Create drawer
            drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
            options = drawer.drawOptions()
            options.baseFontSize = 3
            # drawer.SetDrawOptions = options
            drawer.SetLineWidth(4)
            
            # Set highlighting
            highlight_atoms = undefined_atoms
            highlight_atom_colors = {atom: (1.0, 0.0, 0.0) for atom in highlight_atoms}  # Red
            highlight_atom_radii = {atom: 0.8 for atom in highlight_atoms}
            
            # Draw options
            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms,
                               highlightAtomColors=highlight_atom_colors,
                               highlightAtomRadii=highlight_atom_radii)
            drawer.FinishDrawing()
            
            # Convert to image
            img_data = drawer.GetDrawingText()
            img = Image.open(io.BytesIO(img_data))
            
            # Display in subplot
            ax.imshow(img)
            ax.axis('off')
            
            # Add text information
            title = f"Molecule {idx + 1}\n"
            title += f"Original: {mol_info['smiles'][:30]}{'...' if len(mol_info['smiles']) > 30 else ''}\n"
            title += f"Canonical: {mol_info['canonical_smiles'][:30]}{'...' if len(mol_info['canonical_smiles']) > 30 else ''}\n"
            title += f"Undefined centers: {mol_info['undefined_centers']}/{mol_info['total_centers']}"
            
            ax.set_title(title, fontsize=12, pad=10)
        
        # Main title
        fig.suptitle(f'Molecules with Undefined Chirality in {db_name}', 
                    fontsize=23)
        
        # Add legend
        from matplotlib.lines import Line2D
        legend_elements = [Line2D([0], [0], marker='o', color='w', 
                                 markerfacecolor='r', markersize=10, 
                                 label='Undefined chiral center')]
        fig.legend(handles=legend_elements, loc='lower center', 
                  bbox_to_anchor=(0.5, 0.2), ncol=1, fontsize=16)
        
        plotting.plt.tight_layout()
        
        if save_path:
            plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        #

    def plot_undefined_chirality_comparison(self, mode='both', use_log=False, add_inset=False, 
                                           figsize=(14, 8), save_path=None):
        """
        Plot undefined chirality comparison with options for log scale and inset
        
        Args:
            mode: 'both', 'count', or 'ratio'
            use_log: Use log scale for count plot
            add_inset: Add inset plot for small values
        """
        if not self.scoring_results:
            print("No scoring results available.")
            return
            
        db_names = []
        undefined_data = []
        
        for db_name, results in self.scoring_results.items():
            snap = results['snapshot']

            stereo_data = snap.analysis_results.get('Stereochemistry Completeness', {})
            
            undefined_dist = stereo_data.get('Undefined chirality count distribution', {
                '1 undefined': 0,
                '2 undefined': 0,
                '3+ undefined': 0
            })
            
            total_undefined = stereo_data.get('Molecules with undefined chirality', 0)
            total_molecules = snap.num_molecules
            undefined_ratio = (total_undefined / total_molecules * 100) if total_molecules > 0 else 0
            
            db_names.append(db_name)
            undefined_data.append({
                '1 undefined': undefined_dist.get('1 undefined', 0),
                '2 undefined': undefined_dist.get('2 undefined', 0),
                '3+ undefined': undefined_dist.get('3+ undefined', 0),
                'total': total_undefined,
                'ratio': undefined_ratio
            })
        
        x = np.arange(len(db_names))
        
        if mode == 'both':
            fig, (ax1, ax2) = plotting.plt.subplots(1, 2, figsize=figsize)
            
            # Left plot: Count (stacked bar)
            width = 0.6
            bottoms = np.zeros(len(db_names))

            patterns = ['xx', '..', '//']
            labels = ['3+ undefined', '2 undefined', '1 undefined']

            for i, label in enumerate(labels):
                values = [d[label] for d in undefined_data]
                bars = ax1.bar(x, values, width, bottom=bottoms,
                              color=[self.colors[j % len(self.colors)] for j in range(len(db_names))],
                              alpha=0.8, label=label)

                for bar in bars:
                    bar.set_hatch(patterns[i])

                # Remove individual value labels (we'll add total at the top instead)

                bottoms += values

            # Add total count and percentage labels at the top of each stacked bar
            for i, db_name in enumerate(db_names):
                total_count = undefined_data[i]['total']
                total_ratio = undefined_data[i]['ratio']
                if total_count > 0:
                    ax1.text(i, bottoms[i] + max(bottoms) * 0.02,
                            f'{total_count}',
                            ha='center', va='bottom', fontsize=10)
            
            ax1.set_xlabel('Database', fontsize=21)
            ax1.set_ylabel('Number of Molecules' + (' (log scale)' if use_log else ''), fontsize=21)
            ax1.set_title('Undefined Chirality Count', fontsize=22)
            ax1.set_xticks(x)
            ax1.set_xticklabels(db_names, rotation=45, ha='right')
            ax1.tick_params(labelsize=18)
            ax1.legend()
            
            if use_log:
                ax1.set_yscale('log')
                ax1.set_ylim(bottom=0.1)
            
            # Add inset if requested
            if add_inset and not use_log:
                # Find databases with small values
                total_counts = [d['total'] for d in undefined_data]
                threshold = np.percentile(total_counts, 50)
                
                # Create inset
                axins = plotting.inset_axes(ax1, width="40%", height="40%", loc='upper right')
                
                # Plot only small values in inset
                small_indices = [i for i, count in enumerate(total_counts) if count <= threshold]
                if small_indices:
                    x_small = np.arange(len(small_indices))
                    bottoms_small = np.zeros(len(small_indices))
                    
                    for i, label in enumerate(labels):
                        values_small = [undefined_data[idx][label] for idx in small_indices]
                        bars_small = axins.bar(x_small, values_small, 0.6, bottom=bottoms_small,
                                              color=[self.colors[idx % len(self.colors)] for idx in small_indices],
                                              alpha=0.8)
                        for bar in bars_small:
                            bar.set_hatch(patterns[i])
                        bottoms_small += values_small
                    
                    axins.set_xticks(x_small)
                    axins.set_xticklabels([db_names[idx] for idx in small_indices], 
                                          rotation=45, ha='right', fontsize=8)
                    axins.set_title('Detail View', fontsize=10)
                    axins.tick_params(labelsize=8)
            
            # Right plot: Ratio (stacked bar showing percentage relative to total database molecules)
            width = 0.6
            bottoms_ratio = np.zeros(len(db_names))

            # Calculate percentage of each type relative to total database molecules
            for i, label in enumerate(labels):
                values_ratio = []
                for d in undefined_data:
                    # d['ratio'] is the total undefined ratio relative to total molecules
                    # We need to calculate what portion of that comes from each type
                    total_undefined = d['1 undefined'] + d['2 undefined'] + d['3+ undefined']
                    if total_undefined > 0:
                        # Proportion of this type within all undefined chirality molecules
                        type_proportion = d[label] / total_undefined
                        # Multiply by total undefined ratio to get percentage relative to total database
                        ratio_pct = type_proportion * d['ratio']
                    else:
                        ratio_pct = 0
                    values_ratio.append(ratio_pct)

                bars = ax2.bar(x, values_ratio, width, bottom=bottoms_ratio,
                              color=[self.colors[j % len(self.colors)] for j in range(len(db_names))],
                              alpha=0.8, label=label)

                for bar in bars:
                    bar.set_hatch(patterns[i])

                # Remove individual value labels (we'll add total at the top instead)

                bottoms_ratio += values_ratio

            # Add total count and percentage labels at the top of each stacked bar for ratio plot
            for i, db_name in enumerate(db_names):
                total_count = undefined_data[i]['total']
                total_ratio = undefined_data[i]['ratio']
                if total_count > 0:
                    ax2.text(i, bottoms_ratio[i] + max(bottoms_ratio) * 0.02,
                            f'{total_ratio:.2f}%',
                            ha='center', va='bottom', fontsize=10)

            ax2.set_xlabel('Database', fontsize=21)
            ax2.set_ylabel('Undefined Chirality Ratio', fontsize=21)
            ax2.set_title('Undefined Chirality Ratio', fontsize=22)
            ax2.set_xticks(x)
            ax2.set_xticklabels(db_names, rotation=45, ha='right')
            ax2.tick_params(labelsize=18)
            # Set y-axis limit based on maximum total undefined ratio plus some margin
            max_ratio = max([d['ratio'] for d in undefined_data]) if undefined_data else 100
            ax2.set_ylim(0, max_ratio + 5)
            ax2.legend()
            
            plotting.plt.suptitle('Undefined Chirality Comparison Across Databases' + 
                        (' (Log Scale)' if use_log else '') + 
                        (' with Inset' if add_inset else ''), 
                        fontsize=23)
            
        plotting.plt.tight_layout()
        
        suffix = f"_{mode}"
        if use_log:
            suffix += "_log"
        if add_inset:
            suffix += "_inset"
        
        if save_path is None:
            save_path = os.path.join(self.comparison_dir, f'undefined_chirality{suffix}.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

