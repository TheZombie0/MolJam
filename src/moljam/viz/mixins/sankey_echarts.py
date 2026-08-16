from .._common import *
from .. import plotting
from ..sankey_backend import (
    EchartsSankey,
    PYECHARTS_AVAILABLE,
    SNAPSHOT_AVAILABLE,
    make_snapshot,
    opts,
    snapshot,
)


class SankeyEchartsPlotMixin:
    def plot_sankey_for_database(self, db_name, save_path=None):
        """
        Plot Sankey diagram for a single database (generates both versions by default)

        Args:
            db_name: Database name
            save_path: Base path for saving files (will generate both .html and .png if available)
        """
        # Generate both versions for comparison
        self.plot_sankey_both_versions(db_name, save_path)

    def plot_sankey_both_versions(self, db_name, save_path=None):
        """
        Generate both 2-layer and 3-layer Sankey diagrams for comparison

        Args:
            db_name: Database name
            save_path: Base path for saving files
        """
        print(f"\nGenerating both Sankey versions for {db_name}...")

        # Generate 2-layer version
        print("  - Creating 2-layer version (Category -> Targets with colored flows)...")
        self._plot_sankey_version(db_name, use_three_layer=False, save_path=save_path)

        # Generate 3-layer version
        print("  - Creating 3-layer version (Category -> Metrics -> Targets)...")
        self._plot_sankey_version(db_name, use_three_layer=True, save_path=save_path)

        print(f"  Both versions saved for {db_name}!")

    def _plot_sankey_version(self, db_name, use_three_layer=False, save_path=None):
        """
        Plot Sankey diagram for a single database showing score, deduction and penalty flows using pyecharts

        Args:
            db_name: Database name
            use_three_layer: If True, use 3-layer structure; if False, use 2-layer structure
            save_path: Base path for saving files (will generate both .html and .png if available)
        """
        if db_name not in self.scoring_results:
            print(f"Database {db_name} not found in scoring results.")
            return

        scorer = self.scoring_results[db_name]['scorer']

        if not PYECHARTS_AVAILABLE:
            print("PyEcharts not available, using matplotlib fallback.")
            self.plot_sankey_for_database_matplotlib_fallback(db_name, save_path=save_path)
            return

        # Get metric penalties
        metric_penalties = scorer.scores.get('Metric Penalties', {})

        # Prepare nodes and links
        nodes = []
        links = []

        # Define colors for different categories and their metrics
        category_colors = {
            'Structural Integrity': ['#FF6B6B', '#FF8E8E', '#FFB1B1'],  # Red variants
            'Data Quality': ['#4ECDC4', '#70D7D0'],  # Teal variants
            'Experimental Information Quality': ['#45B7D1', '#6BC5D8', '#8CD3E5'],  # Blue variants (3 metrics)
            'Chemical Space Coverage': ['#96CEB4', '#B8DBB8'],  # Green variants
            'Data Distribution': ['#FECA57', '#FED782']  # Yellow variants
        }

        if use_three_layer:
            # === APPROACH 3: Three-layer structure ===

            # Add left side category nodes (layer 1)
            for category, metrics in self.categories.items():
                if category not in scorer.scores:
                    continue
                colors = category_colors.get(category, ['#95A5A6'] * 10)
                nodes.append({
                    "name": category,
                    "itemStyle": {"color": colors[0]}  # Use primary color of category
                })

            # Add middle layer metric nodes (layer 2)
            node_mapping = {}
            for category, metrics in self.categories.items():
                if category not in scorer.scores:
                    continue
                colors = category_colors.get(category, ['#95A5A6'] * 10)

                for metric_idx, metric in enumerate(metrics):
                    color = colors[metric_idx % len(colors)]
                    nodes.append({
                        "name": metric,
                        "itemStyle": {"color": color}
                    })
                    source_key = f"{category}::{metric}"
                    node_mapping[source_key] = metric

            # Add target nodes (layer 3 - right side)
            actual_scores_node = "Actual Scores"
            score_deductions_node = "Score Deductions"
            penalty_deductions_node = "Penalty Deductions"

            nodes.extend([
                {"name": actual_scores_node, "itemStyle": {"color": "#27AE60"}},  # Green
                {"name": score_deductions_node, "itemStyle": {"color": "#E74C3C"}},  # Red
                {"name": penalty_deductions_node, "itemStyle": {"color": "#8E44AD"}}  # Purple
            ])

            # Create links: Category -> Metric (layer 1 -> 2)
            for category, metrics in self.categories.items():
                if category not in scorer.scores:
                    continue
                colors = category_colors.get(category, ['#95A5A6'] * 10)

                for metric_idx, metric in enumerate(metrics):
                    metric_color = colors[metric_idx % len(colors)]
                    source_key = f"{category}::{metric}"

                    # Calculate total flow for this metric (score + deduction + penalty)
                    actual_score = scorer.scores[category].get(metric, 0)
                    if actual_score is None:
                        actual_score = 0
                    score_deduction = 10 - actual_score
                    penalty = metric_penalties.get(source_key, 0)
                    total_flow = actual_score + score_deduction + penalty

                    if total_flow > 0:
                        links.append({
                            "source": category,
                            "target": metric,
                            "value": round(total_flow, 2),
                            "lineStyle": {"color": metric_color, "opacity": 0.6}
                        })

            # Create links: Metric -> Target (layer 2 -> 3)
            for category, metrics in self.categories.items():
                if category not in scorer.scores:
                    continue
                colors = category_colors.get(category, ['#95A5A6'] * 10)

                for metric_idx, metric in enumerate(metrics):
                    source_key = f"{category}::{metric}"
                    if source_key not in node_mapping:
                        continue

                    metric_color = colors[metric_idx % len(colors)]
                    actual_score = scorer.scores[category].get(metric, 0)
                    if actual_score is None:
                        actual_score = 0
                    score_deduction = 10 - actual_score
                    penalty = metric_penalties.get(source_key, 0)

                    # Links from metric to targets
                    if actual_score > 0:
                        links.append({
                            "source": metric,
                            "target": actual_scores_node,
                            "value": round(actual_score, 2),
                            "lineStyle": {"color": metric_color, "opacity": 0.7}
                        })

                    if score_deduction > 0:
                        links.append({
                            "source": metric,
                            "target": score_deductions_node,
                            "value": round(score_deduction, 2),
                            "lineStyle": {"color": metric_color, "opacity": 0.5}
                        })

                    if penalty > 0:
                        links.append({
                            "source": metric,
                            "target": penalty_deductions_node,
                            "value": round(penalty, 2),
                            "lineStyle": {"color": metric_color, "opacity": 0.8}
                        })

        else:
            # === APPROACH 2 IMPROVED: Two-layer with better visualization ===

            # Add left side category nodes with category colors
            for category, metrics in self.categories.items():
                if category not in scorer.scores:
                    continue
                colors = category_colors.get(category, ['#95A5A6'] * 10)
                # Use the primary color of the category with some transparency
                nodes.append({
                    "name": category,
                    "itemStyle": {"color": colors[0], "opacity": 0.8}
                })

            # Add target nodes (right side)
            actual_scores_node = "Actual Scores"
            score_deductions_node = "Score Deductions"
            penalty_deductions_node = "Penalty Deductions"

            nodes.extend([
                {"name": actual_scores_node, "itemStyle": {"color": "#27AE60"}},  # Green
                {"name": score_deductions_node, "itemStyle": {"color": "#E74C3C"}},  # Red
                {"name": penalty_deductions_node, "itemStyle": {"color": "#8E44AD"}}  # Purple
            ])

            # Create detailed links with metric names in tooltips
            for category, metrics in self.categories.items():
                if category not in scorer.scores:
                    continue

                colors = category_colors.get(category, ['#95A5A6'] * 10)

                for metric_idx, metric in enumerate(metrics):
                    source_key = f"{category}::{metric}"
                    metric_color = colors[metric_idx % len(colors)]

                    actual_score = scorer.scores[category].get(metric, 0)
                    if actual_score is None:
                        actual_score = 0
                    score_deduction = 10 - actual_score
                    penalty = metric_penalties.get(source_key, 0)

                    # Create links with metric names embedded
                    if actual_score > 0:
                        links.append({
                            "source": category,
                            "target": actual_scores_node,
                            "value": round(actual_score, 2),
                            "lineStyle": {"color": metric_color, "opacity": 0.7}
                        })

                    if score_deduction > 0:
                        links.append({
                            "source": category,
                            "target": score_deductions_node,
                            "value": round(score_deduction, 2),
                            "lineStyle": {"color": metric_color, "opacity": 0.5}
                        })

                    if penalty > 0:
                        links.append({
                            "source": category,
                            "target": penalty_deductions_node,
                            "value": round(penalty, 2),
                            "lineStyle": {"color": metric_color, "opacity": 0.8}
                        })

        # Create Sankey chart
        sankey = (
            EchartsSankey(init_opts=opts.InitOpts(
                width="1000px",
                height="900px" if use_three_layer else "800px"  # Taller canvas for 3-layer
            ))
            .add(
                "",  # Empty series name to avoid overlap
                nodes=nodes,
                links=links,
                linestyle_opt=opts.LineStyleOpts(opacity=0.6, curve=0.5),
                label_opts=opts.LabelOpts(position="right", font_size=11),
                node_gap=25 if not use_three_layer else 8,  # Much smaller gap for 3-layer to compress middle metrics
                node_width=15 if not use_three_layer else 8,  # Narrower nodes in 3-layer for better fit
            )
            .set_global_opts(
                title_opts=opts.TitleOpts(
                    title=f"Score Flow Analysis - {db_name}",
                    subtitle=f"Normalized: {scorer.scores.get('Normalized Score', 0):.1f}/100 | "
                             f"Final: {scorer.scores.get('Final Adjusted Score', 0):.1f}/100",
                    pos_left="center",
                    pos_top="0px" if use_three_layer else "20px"  # More space for 3-layer title
                ),
                tooltip_opts=opts.TooltipOpts(trigger="item", trigger_on="mousemove"),
            )
        )

        # Save both HTML and PNG (if available)
        suffix = "_3layer" if use_three_layer else "_2layer"
        base_path = save_path or os.path.join(self.scoring_results[db_name]['output_dir'], f'sankey_score_flow{suffix}')

        # Always save HTML
        html_path = f"{base_path}.html"
        try:
            sankey.render(html_path)
            print(f"    Sankey diagram HTML saved to: {html_path}")
        except Exception as e:
            print(f"    Failed to save HTML: {str(e)}")

        # Try to save PNG if selenium is available
        if SNAPSHOT_AVAILABLE:
            png_path = f"{base_path}.png"
            try:
                make_snapshot(snapshot, sankey.render(), png_path)
                print(f"    Sankey diagram PNG saved to: {png_path}")
            except Exception as e:
                print(f"    Failed to save PNG: {str(e)}")
        else:
            print("    selenium-snapshot not available, PNG output skipped.")
