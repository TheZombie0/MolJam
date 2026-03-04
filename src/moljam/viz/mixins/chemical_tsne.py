from .._common import *
from .. import plotting
from ..chem import calculate_fp_parallel


class ChemicalSpaceTsnePlotMixin:
    def plot_chemical_space_tsne(self, n_samples=None, figsize=(12, 10), save_path=None):
        """
        Plot t-SNE visualization of chemical space for each database separately
        With PCA preprocessing and visualization

        Args:
            n_samples: Number of molecules to sample per database. If None, use all molecules.
            figsize: Figure size for each plot
            save_path: Base path for saving plots (will be suffixed with database name)
        """
        try:
            from rdkit.Chem import AllChem
            from sklearn.manifold import TSNE
            from sklearn.decomposition import PCA
        except ImportError:
            print("scikit-learn required for t-SNE visualization")
            return

        if not self.scoring_results:
            print("No scoring results available.")
            return

        # Process each database separately
        for db_idx, (db_name, results) in enumerate(self.scoring_results.items()):
            snap = results['snapshot']

            if not hasattr(snap, 'valid_df') or snap.valid_df.empty:
                print(f"No valid data for {db_name}, skipping...")
                continue

            print(f"Processing t-SNE for {db_name}...")

            # Sample or use all molecules
            if n_samples is None:
                sampled_df = snap.valid_df
                sample_size = len(sampled_df)
            else:
                sample_size = min(n_samples, len(snap.valid_df))
                sampled_df = snap.valid_df.sample(n=sample_size, random_state=42)

            # Calculate fingerprints with parallel processing
            canonical_smiles = sampled_df['canonical_smiles'].tolist()

            # Use multiprocessing for fingerprint calculation
            print(f"  Calculating fingerprints using CPU cores...")
            if self.use_parallel:
                n_workers = min(cpu_count(), 100)
                with Pool(n_workers) as pool:
                    fps_results = pool.map(calculate_fp_parallel, canonical_smiles)
            else:
                # Serial processing
                fps_results = [calculate_fp_parallel(smiles) for smiles in canonical_smiles]

            # Filter out None values
            fps = [fp for fp in fps_results if fp is not None]

            if len(fps) < 10:
                print(f"Not enough valid molecules for t-SNE in {db_name} (only {len(fps)}), skipping...")
                continue

            # Convert to numpy array
            fp_array = np.array(fps, dtype=np.float32)

            # Perform PCA dimensionality reduction to 50 dimensions
            print(f"  Performing PCA reduction to 50 dimensions...")
            pca = PCA(n_components=30, random_state=42)
            fp_pca = pca.fit_transform(fp_array)
            explained_variance_total = pca.explained_variance_ratio_.sum()
            print(f"    PCA explained variance: {explained_variance_total:.2%}")

            # Plot 1: Explained Variance Ratio
            print(f"  Generating PCA explained variance plot...")
            fig, ax = plotting.plt.subplots(figsize=(12, 6))

            cumsum_var = np.cumsum(pca.explained_variance_ratio_)
            ax.bar(range(1, len(pca.explained_variance_ratio_) + 1),
                   pca.explained_variance_ratio_,
                   alpha=0.6, color=self.colors[db_idx % len(self.colors)],
                   label='Individual explained variance')
            ax.plot(range(1, len(cumsum_var) + 1), cumsum_var,
                   'o-', color='red', linewidth=2, markersize=4,
                   label='Cumulative explained variance')

            ax.set_xlabel('Principal Component', fontsize=14)
            ax.set_ylabel('Explained Variance Ratio', fontsize=14)
            ax.set_title(f'PCA Explained Variance - {db_name}\n(Total variance explained: {explained_variance_total:.2%})',
                        fontsize=16, fontweight='bold')
            ax.legend(fontsize=12)
            ax.grid(True, alpha=0.3)

            plotting.plt.tight_layout()

            pca_var_path = os.path.join(self.scoring_results[db_name]['output_dir'],
                                        f'pca_explained_variance_{db_name}.png')
            plotting.plt.savefig(pca_var_path, dpi=500, bbox_inches='tight')
            plotting.plt.close()
            print(f"    Saved to {pca_var_path}")

            # Plot 2: PCA 50-dimensional space (PC1 vs PC2)
            print(f"  Generating PCA space visualization...")
            fig, ax = plotting.plt.subplots(figsize=figsize)

            ax.scatter(fp_pca[:, 0], fp_pca[:, 1],
                      c=[self.colors[db_idx % len(self.colors)]],
                      alpha=0.6, s=5, edgecolors='white', linewidth=0.1)

            ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})', fontsize=21)
            ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})', fontsize=21)
            ax.set_title(f'PCA 30D Space Visualization - {db_name}\n(n={len(fps)} molecules)',
                        fontsize=23)
            ax.tick_params(labelsize=18)
            ax.grid(True, alpha=0.3)

            plotting.plt.tight_layout()

            pca_space_path = os.path.join(self.scoring_results[db_name]['output_dir'],
                                          f'pca_space_{db_name}.png')
            plotting.plt.savefig(pca_space_path, dpi=500, bbox_inches='tight')
            plotting.plt.close()
            print(f"    Saved to {pca_space_path}")

            # Plot 3: Perform t-SNE on PCA-reduced data
            print(f"  Performing t-SNE on PCA-reduced data ({len(fps)} molecules)...")
            perplexity = 30

            tsne = TSNE(n_components=2, random_state=42, init='random', perplexity=perplexity)
            embeddings = tsne.fit_transform(fp_pca)

            # Plot
            fig, ax = plotting.plt.subplots(figsize=figsize)

            ax.scatter(embeddings[:, 0], embeddings[:, 1],
                      c=[self.colors[db_idx % len(self.colors)]],
                      alpha=0.6, s=5, edgecolors='white', linewidth=0.1)

            ax.set_xlabel('t-SNE 1', fontsize=21)
            ax.set_ylabel('t-SNE 2', fontsize=21)
            ax.set_title(f'Chemical Space Visualization - {db_name}\n(n={len(fps)} molecules)',
                        fontsize=23)
            ax.tick_params(labelsize=18)
            ax.grid(True, alpha=0.3)

            plotting.plt.tight_layout()

            # Save plot to database-specific folder
            if save_path is None:
                # Save to the database-specific folder (same as top_scaffolds)
                save_path_db = os.path.join(self.scoring_results[db_name]['output_dir'],
                                           f'chemical_space_tsne.png')
            else:
                save_path_db = save_path.replace('.png', f'_{db_name}.png')

            plotting.plt.savefig(save_path_db, dpi=500, bbox_inches='tight')
            plotting.plt.close()
            print(f"  Saved to {save_path_db}")

            # Plot2
            fig, ax = plotting.plt.subplots(figsize=figsize)

            ax.scatter(embeddings[:, 0], embeddings[:, 1],
                      c=[self.colors[db_idx % len(self.colors)]],
                      alpha=0.6, s=5, edgecolors='white', linewidth=0.1)

            ax.set_xlabel('t-SNE 1', fontsize=21)
            ax.set_ylabel('t-SNE 2', fontsize=21)

            # Use axis ranges from combined t-SNE plot if available, otherwise use default
            if self.tsne_axis_ranges['xlim'] is not None and self.tsne_axis_ranges['ylim'] is not None:
                ax.set_xlim(self.tsne_axis_ranges['xlim'])
                ax.set_ylim(self.tsne_axis_ranges['ylim'])
                ax.tick_params(labelsize=18)
            else:
                # Fallback to original hardcoded values if ranges not available
                ax.set_xlim(-100, 100)
                ax.set_ylim(-100, 100)
                ax.tick_params(labelsize=18)

            ax.set_title(f'Chemical Space Visualization - {db_name}\n(n={len(fps)} molecules)',
                        fontsize=23)
            ax.grid(True, alpha=0.3)

            plotting.plt.tight_layout()

            # Save plot to database-specific folder
            if save_path is None:
                # Save to the database-specific folder (same as top_scaffolds)
                save_path_db = os.path.join(self.scoring_results[db_name]['output_dir'],
                                           f'chemical_space_tsne_allsize.png')
            else:
                save_path_db = save_path.replace('.png', f'_{db_name}.png')

            plotting.plt.savefig(save_path_db, dpi=500, bbox_inches='tight')
            plotting.plt.close()
            print(f"  Saved to {save_path_db}")

    def plot_all_databases_combined_tsne(self, n_samples_per_db=None, figsize=(14, 12), save_path=None,
                                         plot_individual=True):
        """
        Plot t-SNE visualization with all databases combined in one plot
        With PCA preprocessing and visualization
        Also plots individual t-SNE for each database showing its position in the combined space

        Args:
            n_samples_per_db: Number of molecules to sample per database. If None, use all molecules.
            figsize: Figure size
            save_path: Path for saving the plot
            plot_individual: If True, also generate individual plots for each database
        """
        try:
            from rdkit.Chem import AllChem
            from sklearn.manifold import TSNE
            from sklearn.decomposition import PCA
        except ImportError:
            print("scikit-learn required for t-SNE visualization")
            return

        if not self.scoring_results:
            print("No scoring results available.")
            return

        all_fps = []
        all_labels = []
        db_sample_counts = {}

        print("Collecting molecules from all databases...")
        for db_name, results in self.scoring_results.items():
            snap = results['snapshot']

            if hasattr(snap, 'valid_df') and not snap.valid_df.empty:
                # Sample or use all molecules
                if n_samples_per_db is None:
                    sampled_df = snap.valid_df
                else:
                    sample_size = min(n_samples_per_db, len(snap.valid_df))
                    sampled_df = snap.valid_df.sample(n=sample_size, random_state=42)

                print(f"  {db_name}: {len(sampled_df)} molecules")

                # Calculate fingerprints with parallel processing
                canonical_smiles = sampled_df['canonical_smiles'].tolist()

                # Use multiprocessing for fingerprint calculation
                if self.use_parallel:
                    n_workers = min(cpu_count(), 100)
                    with Pool(n_workers) as pool:
                        fps_results = pool.map(calculate_fp_parallel, canonical_smiles)
                else:
                    # Serial processing
                    fps_results = [calculate_fp_parallel(smiles) for smiles in canonical_smiles]

                # Filter out None values and add labels
                for fp in fps_results:
                    if fp is not None:
                        all_fps.append(fp)
                        all_labels.append(db_name)

                db_sample_counts[db_name] = sum(1 for fp in fps_results if fp is not None)

        if len(all_fps) < 10:
            print("Not enough valid molecules for t-SNE")
            return

        print(f"Total molecules collected: {len(all_fps)}")

        # Convert to numpy array
        fp_array = np.array(all_fps, dtype=np.float32)
        all_labels = np.array(all_labels)

        # Perform PCA dimensionality reduction to 50 dimensions
        print("Performing PCA reduction to 50 dimensions...")
        pca = PCA(n_components=30, random_state=42)
        fp_pca = pca.fit_transform(fp_array)
        explained_variance_total = pca.explained_variance_ratio_.sum()
        print(f"  PCA explained variance: {explained_variance_total:.2%}")

        # Plot 1: PCA Explained Variance Ratio
        print("Generating PCA explained variance plot...")
        fig, ax = plotting.plt.subplots(figsize=(15, 6))

        cumsum_var = np.cumsum(pca.explained_variance_ratio_)
        ax.bar(range(1, len(pca.explained_variance_ratio_) + 1),
               pca.explained_variance_ratio_,
               alpha=0.6, color='steelblue',
               label='Individual explained variance')
        ax.plot(range(1, len(cumsum_var) + 1), cumsum_var,
               'o-', color='red', linewidth=2, markersize=4,
               label='Cumulative explained variance')

        ax.set_xlabel('Principal Component', fontsize=14)
        ax.set_ylabel('Explained Variance Ratio', fontsize=14)
        ax.set_title(f'PCA Explained Variance (Combined Databases)\n(Total variance explained: {explained_variance_total:.2%})',
                    fontsize=16, fontweight='bold')
        ax.legend(fontsize=12)
        ax.grid(True, alpha=0.3)

        plotting.plt.tight_layout()

        pca_var_path = os.path.join(self.coverage_dir, 'pca_explained_variance_combined.png')
        plotting.plt.savefig(pca_var_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()
        print(f"  Saved to {pca_var_path}")

        # Plot 2: PCA 50-dimensional space (PC1 vs PC2) showing all databases
        print("Generating combined PCA space visualization...")
        fig, ax = plotting.plt.subplots(figsize=(14, 12))

        # Plot each database with different color
        for db_idx, db_name in enumerate(self.scoring_results.keys()):
            mask = all_labels == db_name
            if mask.any():
                ax.scatter(fp_pca[mask, 0], fp_pca[mask, 1],
                          c=[self.colors[db_idx % len(self.colors)]],
                          label=f'{db_name} (n={mask.sum()})',
                          alpha=0.6, s=5, edgecolors='white', linewidth=0.1)

        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})', fontsize=21)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})', fontsize=21)
        ax.set_title(f'Combined PCA 30D Space Visualization (All Databases)\n(n={len(all_fps)} molecules)',
                    fontsize=23)
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 0.8), markerscale=4, frameon=True, fancybox=True, fontsize=16)
        ax.tick_params(labelsize=18)
        ax.grid(True, alpha=0.3)

        plotting.plt.tight_layout()

        pca_space_path = os.path.join(self.coverage_dir, 'pca_space_combined.png')
        plotting.plt.savefig(pca_space_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()
        print(f"  Saved to {pca_space_path}")

        # Plot 3: Perform t-SNE on PCA-reduced data
        print("Performing t-SNE on PCA-reduced data...")
        perplexity = 30

        tsne = TSNE(n_components=2, random_state=42, init='random', perplexity=perplexity)
        embeddings = tsne.fit_transform(fp_pca)

        # Plot combined visualization
        fig, ax = plotting.plt.subplots(figsize=(16, 12))

        # Plot each database with different color
        for db_idx, db_name in enumerate(self.scoring_results.keys()):
            mask = all_labels == db_name
            if mask.any():
                ax.scatter(embeddings[mask, 0], embeddings[mask, 1],
                          c=[self.colors[db_idx % len(self.colors)]],
                          label=f'{db_name} (n={mask.sum()})',
                          alpha=0.6, s=5, edgecolors='white', linewidth=0.1)

        ax.set_xlabel('t-SNE 1', fontsize=21)
        ax.set_ylabel('t-SNE 2', fontsize=21)
        ax.set_title(f'Combined Chemical Space Visualization (t-SNE)\nTotal: {len(all_fps)} molecules',
                    fontsize=23)
        ax.legend(loc='upper left',bbox_to_anchor=(1.05, 0.8), markerscale=4, frameon=True, fancybox=True, fontsize=16)
        ax.tick_params(labelsize=18)
        ax.grid(True, alpha=0.3)

        plotting.plt.tight_layout()

        # Store the axis ranges for use by other t-SNE plots
        self.tsne_axis_ranges['xlim'] = ax.get_xlim()
        self.tsne_axis_ranges['ylim'] = ax.get_ylim()

        if save_path is None:
            save_path = os.path.join(self.coverage_dir, 'chemical_space_tsne_combined.png')
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()
        print(f"Saved combined t-SNE plot to {save_path}")

        # Plot individual database views if requested
        if plot_individual:
            print("\nGenerating individual database t-SNE views...")
            for db_idx, db_name in enumerate(self.scoring_results.keys()):
                mask = all_labels == db_name
                if not mask.any():
                    continue

                print(f"  Creating individual plot for {db_name}...")

                fig, ax = plotting.plt.subplots(figsize=figsize)

                # Plot all databases in gray as background
                background_label_added = False
                for other_idx, other_db in enumerate(self.scoring_results.keys()):
                    if other_db != db_name:
                        other_mask = all_labels == other_db
                        if other_mask.any():
                            label = 'background' if not background_label_added else ''
                            ax.scatter(embeddings[other_mask, 0], embeddings[other_mask, 1],
                                      c='lightgray', alpha=0.2, s=5, edgecolors='none',
                                      label=label)
                            if not background_label_added:
                                background_label_added = True
                            # ax.scatter(embeddings[other_mask, 0], embeddings[other_mask, 1],
                            #           c='lightgray', alpha=0.2, s=10, edgecolors='none',
                            #           label=f'{other_db} (background)' if other_idx == 0 else '')

                # Highlight current database
                ax.scatter(embeddings[mask, 0], embeddings[mask, 1],
                          c=[self.colors[db_idx % len(self.colors)]],
                          label=f'{db_name} (n={mask.sum()})',
                          alpha=0.7, s=30, edgecolors='white', linewidth=0.1)

                ax.set_xlabel('t-SNE 1', fontsize=21)
                ax.set_ylabel('t-SNE 2', fontsize=21)
                ax.set_title(f'Chemical Space Position - {db_name}\n(Shown in context of all databases)',
                            fontsize=23)
                ax.legend(loc='best', markerscale=2, frameon=True, fancybox=True)
                ax.tick_params(labelsize=18)
                ax.grid(True, alpha=0.3)

                plotting.plt.tight_layout()

                # Save to database-specific folder
                individual_save_path = os.path.join(self.scoring_results[db_name]['output_dir'],
                                                   f'chemical_space_tsne_in_context.png')
                plotting.plt.savefig(individual_save_path, dpi=500, bbox_inches='tight')
                plotting.plt.close()
                print(f"    Saved to {individual_save_path}")

        print("\nAll t-SNE visualizations completed!")

