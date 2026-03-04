from .._common import *
from .. import plotting


class ScaffoldPlotMixin:
    def plot_top_scaffolds(self, db_name, top_n=10, figsize=(16, 10), save_path=None):
        """
        Plot top scaffolds with handling for empty scaffolds
        """
        if db_name not in self.scoring_results:
            print(f"Database {db_name} not found.")
            return
            
        snap = self.scoring_results[db_name]['snapshot']
        diversity_data = snap.analysis_results.get('Chemical Diversity', {})
        top_scaffolds = diversity_data.get('Top scaffolds', [])
        
        if not top_scaffolds:
            print(f"No scaffold data available for {db_name}")
            return
        
        top_scaffolds = top_scaffolds[:min(top_n, len(top_scaffolds))]
        
        n_scaffolds = len(top_scaffolds)
        n_cols = min(5, n_scaffolds)
        n_rows = (n_scaffolds + n_cols - 1) // n_cols
        
        fig = plotting.plt.figure(figsize=figsize)
        gs = plotting.GridSpec(n_rows, n_cols, figure=fig, hspace=0.4, wspace=0.1)
        
        for idx, scaffold_info in enumerate(top_scaffolds):
            row = idx // n_cols
            col = idx % n_cols
            ax = fig.add_subplot(gs[row, col])
            
            scaffold_smiles = scaffold_info['scaffold_smiles']
            count = scaffold_info['count']
            percentage = scaffold_info['percentage']
            
            # Handle empty scaffold
            if scaffold_smiles == '' or scaffold_smiles is None:
                ax.text(0.5, 0.5, 'No scaffold',
                       ha='center', va='center', transform=ax.transAxes,
                       fontsize=20)
                ax.axis('off')
                title = f"Rank {idx + 1}\n"
                title += f"Count: {count} ({percentage})"
                # Use explicit positioning to match image subplots
                ax.text(0.5, 0.91, title, fontsize=16, ha='center', va='bottom',
                       transform=ax.transAxes)
                continue
            
            try:
                mol = Chem.MolFromSmiles(scaffold_smiles)
                if mol:
                    from rdkit.Chem import AllChem
                    AllChem.Compute2DCoords(mol)
                    
                    from rdkit.Chem.Draw import rdMolDraw2D
                    from PIL import Image
                    import io
                    
                    drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
                    drawer.DrawMolecule(mol)
                    drawer.FinishDrawing()
                    
                    img_data = drawer.GetDrawingText()
                    img = Image.open(io.BytesIO(img_data))
                    
                    ax.imshow(img)
                    ax.axis('off')
                    
                    title = f"Rank {idx + 1}\n"
                    title += f"Count: {count} ({percentage})"
                    ax.set_title(title, fontsize=16, pad=15)
                else:
                    ax.text(0.5, 0.5, 'Failed to parse scaffold',
                           ha='right', va='center', transform=ax.transAxes)
                    ax.axis('off')
            except Exception as e:
                ax.text(0.5, 0.5, f'Error: {str(e)}',
                       ha='right', va='center', transform=ax.transAxes)
                ax.axis('off')
        
        fig.suptitle(f'Top {n_scaffolds} Scaffolds in {db_name}',
                    fontsize=23)
        
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.scoring_results[db_name]['output_dir'],
                                    f'top_{top_n}_scaffolds.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()

