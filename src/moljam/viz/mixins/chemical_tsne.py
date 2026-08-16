from .._common import *
from .. import plotting
from ..chem import calculate_fp_parallel


CHEMICAL_SPACE_FIGSIZE = (14, 8)
CHEMICAL_SPACE_DETAIL_FIGSIZE = (12, 10)
CHEMICAL_SPACE_AXIS_FONTSIZE = 21
CHEMICAL_SPACE_TITLE_FONTSIZE = 23
CHEMICAL_SPACE_TICK_FONTSIZE = 18
CHEMICAL_SPACE_LEGEND_FONTSIZE = 16
CHEMICAL_SPACE_BACKGROUND_COLOR = '#9e9e9e'
CHEMICAL_SPACE_BACKGROUND_ALPHA = 0.42


class ChemicalSpaceTsnePlotMixin:
    @staticmethod
    def _resolve_reduction_params(num_samples, feature_dim, max_pca_components=30, default_perplexity=30):
        """Pick PCA/t-SNE parameters that remain valid for small sample sets."""
        if num_samples < 3:
            return None, None

        pca_components = min(max_pca_components, num_samples, feature_dim)
        if pca_components < 2:
            return None, None

        perplexity = min(default_perplexity, num_samples - 1)
        if perplexity < 2:
            return None, None

        return pca_components, perplexity

    def plot_chemical_space_tsne(self, n_samples=5000, figsize=CHEMICAL_SPACE_FIGSIZE, save_path=None):
        """
        Plot t-SNE visualization of chemical space for each database separately
        With PCA preprocessing and visualization

        Args:
            n_samples: Number of molecules to sample per database. Defaults to 5000.
                       If None, use all molecules.
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
            scorer = results['scorer']

            if not hasattr(scorer, 'valid_df') or scorer.valid_df.empty:
                print(f"No valid data for {db_name}, skipping...")
                continue

            print(f"Processing t-SNE for {db_name}...")

            # Sample or use all molecules
            if n_samples is None:
                sampled_df = scorer.valid_df
                sample_size = len(sampled_df)
            else:
                sample_size = min(n_samples, len(scorer.valid_df))
                sampled_df = scorer.valid_df.sample(n=sample_size, random_state=42)

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
            if len(fps) < 3:
                print(f"Not enough valid molecules for t-SNE in {db_name} (only {len(fps)}), skipping...")
                continue

            # Convert to numpy array
            fp_array = np.array(fps, dtype=np.float32)
            pca_components, perplexity = self._resolve_reduction_params(len(fps), fp_array.shape[1])
            if pca_components is None:
                print(f"Not enough valid molecules for t-SNE in {db_name} (only {len(fps)}), skipping...")
                continue

            print(f"  Performing PCA reduction to {pca_components} dimensions...")
            pca = PCA(n_components=pca_components, random_state=42)
            fp_pca = pca.fit_transform(fp_array)
            explained_variance_total = pca.explained_variance_ratio_.sum()
            print(f"    PCA explained variance: {explained_variance_total:.2%}")

            # Plot 1: Explained Variance Ratio
            print(f"  Generating PCA explained variance plot...")
            fig, ax = plotting.plt.subplots(figsize=figsize)

            cumsum_var = np.cumsum(pca.explained_variance_ratio_)
            ax.bar(range(1, len(pca.explained_variance_ratio_) + 1),
                   pca.explained_variance_ratio_,
                   alpha=0.6, color=self.colors[db_idx % len(self.colors)],
                   label='Individual explained variance')
            ax.plot(range(1, len(cumsum_var) + 1), cumsum_var,
                   'o-', color='red', linewidth=2, markersize=4,
                   label='Cumulative explained variance')

            ax.set_xlabel('Principal Component', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
            ax.set_ylabel('Explained Variance Ratio', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
            ax.set_title(f'PCA Explained Variance - {db_name}\n(Total variance explained: {explained_variance_total:.2%})',
                        fontsize=CHEMICAL_SPACE_TITLE_FONTSIZE)
            ax.tick_params(labelsize=CHEMICAL_SPACE_TICK_FONTSIZE)
            ax.legend(fontsize=CHEMICAL_SPACE_LEGEND_FONTSIZE)
            ax.grid(True, alpha=0.3)

            plotting.plt.tight_layout()

            pca_var_path = os.path.join(self.scoring_results[db_name]['output_dir'],
                                        f'pca_explained_variance_{db_name}.png')
            plotting.plt.savefig(pca_var_path, dpi=500, bbox_inches='tight')
            plotting.plt.close()
            print(f"    Saved to {pca_var_path}")

            # Plot 2: PCA space (PC1 vs PC2)
            print(f"  Generating PCA space visualization...")
            fig, ax = plotting.plt.subplots(figsize=figsize)

            ax.scatter(fp_pca[:, 0], fp_pca[:, 1],
                      c=[self.colors[db_idx % len(self.colors)]],
                      alpha=0.6, s=5, edgecolors='white', linewidth=0.1)

            ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
            ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
            ax.set_title(f'PCA {pca_components}D Space Visualization - {db_name}\n(n={len(fps)} molecules)',
                        fontsize=CHEMICAL_SPACE_TITLE_FONTSIZE)
            ax.tick_params(labelsize=CHEMICAL_SPACE_TICK_FONTSIZE)
            ax.grid(True, alpha=0.3)

            plotting.plt.tight_layout()

            pca_space_path = os.path.join(self.scoring_results[db_name]['output_dir'],
                                          f'pca_space_{db_name}.png')
            plotting.plt.savefig(pca_space_path, dpi=500, bbox_inches='tight')
            plotting.plt.close()
            print(f"    Saved to {pca_space_path}")

            # Plot 3: Perform t-SNE on PCA-reduced data
            print(
                f"  Performing t-SNE on PCA-reduced data ({len(fps)} molecules, "
                f"perplexity={perplexity})..."
            )
            tsne = TSNE(n_components=2, random_state=42, init='random', perplexity=perplexity)
            embeddings = tsne.fit_transform(fp_pca)

            # Plot
            fig, ax = plotting.plt.subplots(figsize=figsize)

            ax.scatter(embeddings[:, 0], embeddings[:, 1],
                      c=[self.colors[db_idx % len(self.colors)]],
                      alpha=0.6, s=5, edgecolors='white', linewidth=0.1)

            ax.set_xlabel('t-SNE 1', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
            ax.set_ylabel('t-SNE 2', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
            ax.set_title(f'Chemical Space Visualization - {db_name}\n(n={len(fps)} molecules)',
                        fontsize=CHEMICAL_SPACE_TITLE_FONTSIZE)
            ax.tick_params(labelsize=CHEMICAL_SPACE_TICK_FONTSIZE)
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
            fig, ax = plotting.plt.subplots(figsize=CHEMICAL_SPACE_DETAIL_FIGSIZE)

            ax.scatter(embeddings[:, 0], embeddings[:, 1],
                      c=[self.colors[db_idx % len(self.colors)]],
                      alpha=0.6, s=5, edgecolors='white', linewidth=0.1)

            ax.set_xlabel('t-SNE 1', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
            ax.set_ylabel('t-SNE 2', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)

            # Use axis ranges from combined t-SNE plot if available, otherwise use default
            if self.tsne_axis_ranges['xlim'] is not None and self.tsne_axis_ranges['ylim'] is not None:
                ax.set_xlim(self.tsne_axis_ranges['xlim'])
                ax.set_ylim(self.tsne_axis_ranges['ylim'])
                ax.tick_params(labelsize=CHEMICAL_SPACE_TICK_FONTSIZE)
            else:
                # Fallback to original hardcoded values if ranges not available
                ax.set_xlim(-100, 100)
                ax.set_ylim(-100, 100)
                ax.tick_params(labelsize=CHEMICAL_SPACE_TICK_FONTSIZE)

            ax.set_title(f'Chemical Space Visualization - {db_name}\n(n={len(fps)} molecules)',
                        fontsize=CHEMICAL_SPACE_TITLE_FONTSIZE)
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

    def plot_all_databases_combined_tsne(self, n_samples_per_db=5000, figsize=CHEMICAL_SPACE_FIGSIZE, save_path=None,
                                         plot_individual=True):
        """
        Plot t-SNE visualization with all databases combined in one plot
        With PCA preprocessing and visualization
        Also plots individual t-SNE for each database showing its position in the combined space

        Args:
            n_samples_per_db: Number of molecules to sample per database. Defaults to 5000.
                              If None, use all molecules.
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
            scorer = results['scorer']

            if hasattr(scorer, 'valid_df') and not scorer.valid_df.empty:
                # Sample or use all molecules
                if n_samples_per_db is None:
                    sampled_df = scorer.valid_df
                else:
                    sample_size = min(n_samples_per_db, len(scorer.valid_df))
                    sampled_df = scorer.valid_df.sample(n=sample_size, random_state=42)

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

        if len(all_fps) < 3:
            print("Not enough valid molecules for t-SNE")
            return

        print(f"Total molecules collected: {len(all_fps)}")

        # Convert to numpy array
        fp_array = np.array(all_fps, dtype=np.float32)
        all_labels = np.array(all_labels)

        pca_components, perplexity = self._resolve_reduction_params(len(all_fps), fp_array.shape[1])
        if pca_components is None:
            print("Not enough valid molecules for PCA/t-SNE")
            return

        print(f"Performing PCA reduction to {pca_components} dimensions...")
        pca = PCA(n_components=pca_components, random_state=42)
        fp_pca = pca.fit_transform(fp_array)
        explained_variance_total = pca.explained_variance_ratio_.sum()
        print(f"  PCA explained variance: {explained_variance_total:.2%}")

        # Plot 1: PCA Explained Variance Ratio
        print("Generating PCA explained variance plot...")
        fig, ax = plotting.plt.subplots(figsize=figsize)

        cumsum_var = np.cumsum(pca.explained_variance_ratio_)
        ax.bar(range(1, len(pca.explained_variance_ratio_) + 1),
               pca.explained_variance_ratio_,
               alpha=0.6, color='steelblue',
               label='Individual explained variance')
        ax.plot(range(1, len(cumsum_var) + 1), cumsum_var,
               'o-', color='red', linewidth=2, markersize=4,
               label='Cumulative explained variance')

        ax.set_xlabel('Principal Component', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
        ax.set_ylabel('Explained Variance Ratio', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
        ax.set_title(f'PCA Explained Variance (Combined Databases)\n(Total variance explained: {explained_variance_total:.2%})',
                    fontsize=CHEMICAL_SPACE_TITLE_FONTSIZE)
        ax.tick_params(labelsize=CHEMICAL_SPACE_TICK_FONTSIZE)
        ax.legend(fontsize=CHEMICAL_SPACE_LEGEND_FONTSIZE)
        ax.grid(True, alpha=0.3)

        plotting.plt.tight_layout()

        pca_var_path = os.path.join(self.coverage_dir, 'pca_explained_variance_combined.png')
        plotting.plt.savefig(pca_var_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()
        print(f"  Saved to {pca_var_path}")

        # Plot 2: PCA space (PC1 vs PC2) showing all databases
        print("Generating combined PCA space visualization...")
        fig, ax = plotting.plt.subplots(figsize=figsize)

        # Plot each database with different color
        for db_idx, db_name in enumerate(self.scoring_results.keys()):
            mask = all_labels == db_name
            if mask.any():
                ax.scatter(fp_pca[mask, 0], fp_pca[mask, 1],
                          c=[self.colors[db_idx % len(self.colors)]],
                          label=f'{db_name} (n={mask.sum()})',
                          alpha=0.6, s=5, edgecolors='white', linewidth=0.1)

        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
        ax.set_title(
            f'Combined PCA {pca_components}D Space Visualization (All Databases)\n'
            f'(n={len(all_fps)} molecules)',
            fontsize=CHEMICAL_SPACE_TITLE_FONTSIZE,
        )
        ax.legend(
            loc='upper left',
            bbox_to_anchor=(1.05, 0.8),
            markerscale=4,
            frameon=True,
            fancybox=True,
            fontsize=CHEMICAL_SPACE_LEGEND_FONTSIZE,
        )
        ax.tick_params(labelsize=CHEMICAL_SPACE_TICK_FONTSIZE)
        ax.grid(True, alpha=0.3)

        plotting.plt.tight_layout()

        pca_space_path = os.path.join(self.coverage_dir, 'pca_space_combined.png')
        plotting.plt.savefig(pca_space_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()
        print(f"  Saved to {pca_space_path}")

        # Plot 3: Perform t-SNE on PCA-reduced data
        print(f"Performing t-SNE on PCA-reduced data (perplexity={perplexity})...")
        tsne = TSNE(n_components=2, random_state=42, init='random', perplexity=perplexity)
        embeddings = tsne.fit_transform(fp_pca)

        # Plot combined visualization
        fig, ax = plotting.plt.subplots(figsize=figsize)

        # Plot each database with different color
        for db_idx, db_name in enumerate(self.scoring_results.keys()):
            mask = all_labels == db_name
            if mask.any():
                ax.scatter(embeddings[mask, 0], embeddings[mask, 1],
                          c=[self.colors[db_idx % len(self.colors)]],
                          label=f'{db_name} (n={mask.sum()})',
                          alpha=0.6, s=5, edgecolors='white', linewidth=0.1)

        ax.set_xlabel('t-SNE 1', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
        ax.set_ylabel('t-SNE 2', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
        ax.set_title(f'Combined Chemical Space Visualization (t-SNE)\nTotal: {len(all_fps)} molecules',
                    fontsize=CHEMICAL_SPACE_TITLE_FONTSIZE)
        ax.legend(
            loc='upper left',
            bbox_to_anchor=(1.05, 0.8),
            markerscale=4,
            frameon=True,
            fancybox=True,
            fontsize=CHEMICAL_SPACE_LEGEND_FONTSIZE,
        )
        ax.tick_params(labelsize=CHEMICAL_SPACE_TICK_FONTSIZE)
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

                fig, ax = plotting.plt.subplots(figsize=CHEMICAL_SPACE_DETAIL_FIGSIZE)

                # Plot all databases in gray as background
                background_label_added = False
                for other_idx, other_db in enumerate(self.scoring_results.keys()):
                    if other_db != db_name:
                        other_mask = all_labels == other_db
                        if other_mask.any():
                            label = 'background' if not background_label_added else ''
                            ax.scatter(embeddings[other_mask, 0], embeddings[other_mask, 1],
                                      c=CHEMICAL_SPACE_BACKGROUND_COLOR,
                                      alpha=CHEMICAL_SPACE_BACKGROUND_ALPHA,
                                      s=5, edgecolors='none',
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

                ax.set_xlabel('t-SNE 1', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
                ax.set_ylabel('t-SNE 2', fontsize=CHEMICAL_SPACE_AXIS_FONTSIZE)
                ax.set_title(f'Chemical Space Position - {db_name}\n(Shown in context of all databases)',
                            fontsize=CHEMICAL_SPACE_TITLE_FONTSIZE)
                ax.legend(loc='best', markerscale=2, frameon=True, fancybox=True, fontsize=CHEMICAL_SPACE_LEGEND_FONTSIZE)
                ax.tick_params(labelsize=CHEMICAL_SPACE_TICK_FONTSIZE)
                ax.grid(True, alpha=0.3)

                plotting.plt.tight_layout()

                # Save to database-specific folder
                individual_save_path = os.path.join(self.scoring_results[db_name]['output_dir'],
                                                   f'chemical_space_tsne_in_context.png')
                plotting.plt.savefig(individual_save_path, dpi=500, bbox_inches='tight')
                plotting.plt.close()
                print(f"    Saved to {individual_save_path}")

        print("\nAll t-SNE visualizations completed!")
