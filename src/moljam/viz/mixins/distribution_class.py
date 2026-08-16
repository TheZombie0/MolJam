from .._common import *
from .. import plotting


class ClassDistributionPlotMixin:
    @staticmethod
    def _normalize_class_label(label):
        if isinstance(label, (float, np.floating)):
            numeric_label = float(label)
            if numeric_label in (0.0, 1.0) and numeric_label.is_integer():
                return str(int(numeric_label))

        label_str = str(label)
        if label_str in {'0.0', '1.0'}:
            return label_str[:-2]
        return label_str

    @classmethod
    def _format_class_labels(cls, labels):
        return [cls._normalize_class_label(label) for label in labels]

    @staticmethod
    def _shade_color(color, factor):
        from matplotlib import colors as mcolors

        rgb = np.array(mcolors.to_rgb(color))
        return tuple(np.clip(rgb * factor, 0, 1))

    @classmethod
    def _prepare_3d_pie_counts(cls, class_counts, max_slices=10):
        if len(class_counts) <= max_slices:
            return class_counts.copy()

        pie_counts = class_counts.head(max_slices - 1).copy()
        pie_counts.loc['Other'] = class_counts.iloc[max_slices - 1:].sum()
        return pie_counts

    @classmethod
    def _draw_pseudo_3d_pie(
        cls,
        ax,
        values,
        colors,
        *,
        start_angle=90,
        explode_distance=0.12,
        radius=1.0,
        depth=0.24,
        layer_count=28,
        vertical_scale=0.72,
    ):
        from matplotlib.patches import Wedge
        from matplotlib import transforms

        total = float(np.sum(values))
        current_angle = float(start_angle)
        pie_transform = transforms.Affine2D().scale(1.0, vertical_scale) + ax.transData

        for value, color in zip(values, colors):
            if value <= 0 or total <= 0:
                continue

            delta_angle = float(value / total * 360.0)
            theta_mid = np.deg2rad(current_angle + delta_angle / 2.0)
            offset_x = explode_distance * np.cos(theta_mid)
            offset_y = explode_distance * np.sin(theta_mid)

            for layer_offset in np.linspace(-depth, 0.0, layer_count):
                side_wedge = Wedge(
                    (offset_x, offset_y + layer_offset),
                    radius,
                    current_angle,
                    current_angle + delta_angle,
                )
                side_wedge.set_facecolor(cls._shade_color(color, 0.78))
                side_wedge.set_edgecolor('none')
                side_wedge.set_transform(pie_transform)
                ax.add_patch(side_wedge)

            top_wedge = Wedge(
                (offset_x, offset_y),
                radius,
                current_angle,
                current_angle + delta_angle,
            )
            top_wedge.set_facecolor(color)
            top_wedge.set_edgecolor('white')
            top_wedge.set_linewidth(1.2)
            top_wedge.set_transform(pie_transform)
            ax.add_patch(top_wedge)

            current_angle += delta_angle

    @staticmethod
    def _build_3d_pie_style(values):
        values = np.asarray(values, dtype=float)
        style = {
            'start_angle': 90.0,
            'explode_distance': 0.12,
            'radius': 1.0,
            'depth': 0.24,
            'layer_count': 28,
            'vertical_scale': 0.72,
        }

        # For binary pies, center the minority slice at the top so the gap reads symmetrically.
        if len(values) == 2 and values.sum() > 0 and np.count_nonzero(values) == 2:
            minority_angle = float(values.min() / values.sum() * 360.0)
            style.update(
                start_angle=90.0 + minority_angle / 2.0,
                explode_distance=0.18,
                depth=0.16,
            )

        return style

    @staticmethod
    def _build_3d_pie_legend_kwargs(num_slices):
        legend_kwargs = {
            'loc': 'upper center',
            'frameon': True,
            'fancybox': True,
            'handlelength': 1.4,
        }

        if num_slices == 2:
            legend_kwargs.update(
                bbox_to_anchor=(0.5, 0.055),
                ncol=1,
                fontsize=16,
                labelspacing=0.55,
            )
            return legend_kwargs

        legend_columns = min(3, max(2, int(np.ceil(num_slices / 4.0))))
        legend_kwargs.update(
            bbox_to_anchor=(0.5, 0.03),
            ncol=legend_columns,
            fontsize=16,
            columnspacing=1.1,
            labelspacing=0.8,
        )
        return legend_kwargs

    @staticmethod
    def _get_3d_pie_axis_limits(num_slices):
        if num_slices == 2:
            return (-1.6, 1.6), (-1.38, 1.02)
        return (-1.6, 1.6), (-1.32, 1.05)

    def _save_3d_pie_chart(self, db_name, col, class_counts, save_path):
        from matplotlib.patches import Patch

        pie_counts = self._prepare_3d_pie_counts(class_counts)
        labels = self._format_class_labels(pie_counts.index)
        values = pie_counts.values.astype(float)
        colors = [self.colors[i % len(self.colors)] for i in range(len(pie_counts))]
        total = values.sum()
        pie_style = self._build_3d_pie_style(values)
        legend_kwargs = self._build_3d_pie_legend_kwargs(len(values))
        xlim, ylim = self._get_3d_pie_axis_limits(len(values))

        fig, ax = plotting.plt.subplots(figsize=(14, 8))
        self._draw_pseudo_3d_pie(ax, values, colors, **pie_style)

        legend_handles = []
        for label, value, color in zip(labels, values, colors):
            percentage = (value / total * 100) if total > 0 else 0
            legend_handles.append(
                Patch(
                    facecolor=color,
                    edgecolor='white',
                    label=f'{label} ({int(value):,}, {percentage:.1f}%)',
                )
            )

        ax.set_title(f'Class Distribution - {db_name}\n({col})', fontsize=23, pad=18)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.legend(handles=legend_handles, **legend_kwargs)

        plotting.plt.tight_layout()
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close(fig)

    def plot_class_distribution(self, db_name, class_col=None, figsize=(14, 8), save_path=None):
        """Plot class distribution for a specific database"""
        if db_name not in self.scoring_results:
            print(f"Database {db_name} not found in scoring results.")
            return
        
        scorer = self.scoring_results[db_name]['scorer']
        class_cols = self.scoring_results[db_name]['class_cols']
        
        # Use specified class column or first class column
        col = class_col
        if col is None and class_cols:
            col = class_cols[0]
        
        if not col or col not in scorer.df.columns:
            print(f"Class column '{col}' not found in {db_name}.")
            return
        
        # Get class distribution
        class_counts = scorer.df[col].value_counts()
        formatted_labels = self._format_class_labels(class_counts.index)
        
        fig, ax = plotting.plt.subplots(figsize=figsize)
        
        # Create pie chart for binary/multi-class distribution
        if len(class_counts) <= 10:  # Reasonable number of classes for pie chart
            colors = [self.colors[i % len(self.colors)] for i in range(len(class_counts))]
            wedges, texts, autotexts = ax.pie(
                class_counts.values,
                labels=formatted_labels,
                autopct='%1.1f%%',
                startangle=90,
                colors=colors,
            )
            
            # Improve text formatting
            for text in texts:
                text.set_fontsize(18)
            for autotext in autotexts:
                autotext.set_fontsize(18)
                autotext.set_color('white')
                
            ax.set_title(f'Class Distribution - {db_name}\n({col})', fontsize=23)
            
        else:  # Too many classes, use bar chart
            # Limit to top 20 classes
            top_classes = class_counts.head(20)
            
            bars = ax.bar(range(len(top_classes)), top_classes.values,
                         color=self.colors[0], alpha=0.7)
            
            ax.set_xlabel('Class', fontsize=21)
            ax.set_ylabel('Count', fontsize=21)
            ax.set_title(f'Top 20 Class Distribution - {db_name}\n({col})', fontsize=23)
            ax.set_xticks(range(len(top_classes)))
            top_labels = [self._normalize_class_label(c)[:10] for c in top_classes.index]
            ax.set_xticklabels(top_labels, rotation=45, ha='right')
            ax.tick_params(labelsize=18)
        
        plotting.plt.tight_layout()
        
        if save_path is None:
            save_path = os.path.join(self.scoring_results[db_name]['output_dir'],
                                   f'class_distribution_{col}.png')
        save_path_3d = os.path.splitext(save_path)[0] + '_3d.png'
        plotting.plt.savefig(save_path, dpi=500, bbox_inches='tight')
        plotting.plt.close()
        self._save_3d_pie_chart(db_name, col, class_counts, save_path_3d)
