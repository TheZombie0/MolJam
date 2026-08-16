from matplotlib.colors import to_rgb

from .._common import *
from .. import plotting
from .quality import QUALITY_HEATMAP_CATEGORY_BAR_COLOR_SCHEMES


RUNTIME_CATEGORY_COLORS = QUALITY_HEATMAP_CATEGORY_BAR_COLOR_SCHEMES["deep"]


class RuntimePlotMixin:
    def _runtime_results(self):
        plotting.ensure_plotting_imports()

        runtime_results = []
        for db_name, results in self.scoring_results.items():
            scorer = results["scorer"]
            profile = getattr(scorer, "runtime_profile", None)
            if not profile or not profile.get("metric_order"):
                continue
            runtime_results.append((db_name, results, scorer, profile))

        if not runtime_results:
            print("No runtime profile available. Please score databases with run_all_checks() first.")

        return runtime_results

    def _runtime_category_order(self, runtime_results):
        available = set()
        for _, _, _, profile in runtime_results:
            available.update(profile.get("category_metrics", {}))
        return [category for category in self.categories if category in available]

    def _runtime_metric_order(self, runtime_results, category_order):
        available = set()
        for _, _, _, profile in runtime_results:
            available.update(profile.get("metric_order", []))

        metric_order = []
        for category in category_order:
            for metric in self.categories.get(category, []):
                if metric in available:
                    metric_order.append(metric)
        return metric_order

    def _runtime_metric_label_map(self):
        label_map = {}
        for category, metrics in self.categories.items():
            labels = self.categories_labels.get(category, metrics)
            for metric, label in zip(metrics, labels):
                label_map[metric] = label
        return label_map

    def _runtime_category_label(self, category):
        return str(category).replace(" ", "\n")

    def _runtime_metric_labels(self, metric_order):
        label_map = self._runtime_metric_label_map()
        return [label_map.get(metric, metric) for metric in metric_order]

    def _runtime_default_save_path(self, filename):
        return os.path.join(self.runtime_dir, filename)

    def _runtime_prepare_save(self, save_path):
        save_dir = os.path.dirname(save_path)
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)

    def _runtime_save_figure(self, fig, save_path, tight_rect=None):
        self._runtime_prepare_save(save_path)
        if tight_rect is None:
            fig.tight_layout()
        else:
            fig.tight_layout(rect=tight_rect)
        fig.savefig(save_path, dpi=500, bbox_inches="tight")
        plotting.plt.close(fig)

    def _lighten_color(self, color, amount):
        rgb = np.array(to_rgb(color))
        return tuple(rgb + (1.0 - rgb) * amount)

    def _runtime_metric_colors(self, category_order):
        colors = {}
        for category in category_order:
            metrics = self.categories.get(category, [])
            if not metrics:
                continue
            lighten_values = np.linspace(0.10, 0.45, len(metrics))
            for metric, amount in zip(metrics, lighten_values):
                colors[metric] = self._lighten_color(RUNTIME_CATEGORY_COLORS[category], float(amount))
        return colors

    def _runtime_matrix(self, runtime_results, key, order):
        rows = []
        db_names = []
        for db_name, _, _, profile in runtime_results:
            values = profile.get(key, {})
            rows.append([float(values.get(item, 0.0)) for item in order])
            db_names.append(db_name)
        return db_names, np.array(rows, dtype=float)

    def plot_total_runtime_line(self, figsize=(14, 8), save_path=None):
        runtime_results = self._runtime_results()
        if not runtime_results:
            return

        db_names = [db_name for db_name, _, _, _ in runtime_results]
        total_seconds = [profile["total_seconds"] for _, _, _, profile in runtime_results]
        x_positions = np.arange(len(db_names))

        fig, ax = plotting.plt.subplots(figsize=figsize)
        ax.plot(x_positions, total_seconds, marker="o", linewidth=2.5, markersize=8, color="#2F5C8A")
        ax.set_xticks(x_positions)
        ax.set_xticklabels(db_names, rotation=45, ha="right", fontsize=18)
        ax.set_xlabel("Database", fontsize=21)
        ax.set_ylabel("Evaluation Runtime (s)", fontsize=21)
        ax.set_title("Evaluation Runtime by Database", fontsize=22)
        ax.tick_params(axis="y", labelsize=18)
        ax.grid(axis="y", linestyle="--", alpha=0.3)

        y_offset = max(total_seconds) * 0.02 if max(total_seconds) > 0 else 0.02
        for x_pos, seconds in zip(x_positions, total_seconds):
            ax.text(x_pos, seconds + y_offset, f"{seconds:.2f}s", ha="center", va="bottom", fontsize=9)

        fig.suptitle("Database Evaluation Runtime Comparison", fontsize=23)
        save_path = save_path or self._runtime_default_save_path("runtime_total_line.png")
        self._runtime_save_figure(fig, save_path, tight_rect=(0.0, 0.0, 1.0, 0.96))

    def plot_total_runtime_bar(self, figsize=(14, 8), save_path=None):
        runtime_results = self._runtime_results()
        if not runtime_results:
            return

        db_names = [db_name for db_name, _, _, _ in runtime_results]
        total_seconds = [profile["total_seconds"] for _, _, _, profile in runtime_results]

        fig, ax = plotting.plt.subplots(figsize=figsize)
        bars = ax.bar(
            db_names,
            total_seconds,
            color=[self.colors[index % len(self.colors)] for index in range(len(db_names))],
        )
        ax.set_xlabel("Database", fontsize=21)
        ax.set_ylabel("Evaluation Runtime (s)", fontsize=21)
        ax.set_title("Evaluation Runtime by Database", fontsize=22)
        ax.tick_params(axis="x", rotation=45, labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        plotting.plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")

        y_offset = max(total_seconds) * 0.02 if max(total_seconds) > 0 else 0.02
        for bar, seconds in zip(bars, total_seconds):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + y_offset,
                f"{seconds:.2f}s",
                ha="center",
                va="bottom",
                fontsize=9,
            )

        fig.suptitle("Database Evaluation Runtime Comparison", fontsize=23)
        save_path = save_path or self._runtime_default_save_path("runtime_total_bar.png")
        self._runtime_save_figure(fig, save_path, tight_rect=(0.0, 0.0, 1.0, 0.96))

    def plot_category_runtime_percentage_bars(self, figsize=(14, 8), save_path=None):
        runtime_results = self._runtime_results()
        if not runtime_results:
            return

        category_order = self._runtime_category_order(runtime_results)
        db_names, matrix = self._runtime_matrix(runtime_results, "category_percentages", category_order)
        x_positions = np.arange(len(category_order))
        width = min(0.78 / max(len(db_names), 1), 0.28)
        annotate = matrix.size <= 30

        fig, ax = plotting.plt.subplots(figsize=figsize)
        for index, db_name in enumerate(db_names):
            offsets = x_positions + (index - (len(db_names) - 1) / 2) * width
            bars = ax.bar(
                offsets,
                matrix[index],
                width=width,
                color=self.colors[index % len(self.colors)],
                label=db_name,
            )
            if annotate:
                for bar, percentage in zip(bars, matrix[index]):
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.6,
                        f"{percentage:.1f}%",
                        ha="center",
                        va="bottom",
                        fontsize=9,
                    )

        ax.set_xticks(x_positions)
        ax.set_xticklabels(
            [self._runtime_category_label(category) for category in category_order],
            fontsize=18,
        )
        ax.set_xlabel("Category", fontsize=21)
        ax.set_ylabel("Runtime Share (%)", fontsize=21)
        ax.set_title("Category Runtime Percentage by Database", fontsize=22)
        ax.tick_params(axis="y", labelsize=18)
        ax.legend(
            title="Database",
            fontsize=14,
            title_fontsize=14,
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
        )
        ax.set_ylim(0, max(100, float(matrix.max()) * 1.15 if matrix.size else 100))
        fig.suptitle("Runtime Share Across Major Scoring Categories", fontsize=23)

        save_path = save_path or self._runtime_default_save_path("runtime_category_percentages.png")
        self._runtime_save_figure(fig, save_path, tight_rect=(0.0, 0.0, 0.86, 0.96))

    def plot_metric_runtime_percentage_bars(self, figsize=(20, 8), save_path=None):
        runtime_results = self._runtime_results()
        if not runtime_results:
            return

        category_order = self._runtime_category_order(runtime_results)
        metric_order = self._runtime_metric_order(runtime_results, category_order)
        db_names, matrix = self._runtime_matrix(runtime_results, "metric_percentages", metric_order)
        metric_labels = self._runtime_metric_labels(metric_order)
        x_positions = np.arange(len(metric_order))
        width = min(0.78 / max(len(db_names), 1), 0.24)
        annotate = matrix.size <= 36

        fig, ax = plotting.plt.subplots(figsize=figsize)
        for index, db_name in enumerate(db_names):
            offsets = x_positions + (index - (len(db_names) - 1) / 2) * width
            bars = ax.bar(
                offsets,
                matrix[index],
                width=width,
                color=self.colors[index % len(self.colors)],
                label=db_name,
            )
            if annotate:
                for bar, percentage in zip(bars, matrix[index]):
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.4,
                        f"{percentage:.1f}%",
                        ha="center",
                        va="bottom",
                        fontsize=9,
                    )

        ax.set_xticks(x_positions)
        ax.set_xticklabels(metric_labels, rotation=45, ha="right", fontsize=18)
        ax.set_xlabel("Metric", fontsize=21)
        ax.set_ylabel("Runtime Share (%)", fontsize=21)
        ax.set_title("Metric Runtime Percentage by Database", fontsize=22)
        ax.tick_params(axis="y", labelsize=18)
        ax.legend(title="Database", fontsize=13, title_fontsize=13, frameon=False)
        ax.set_ylim(0, max(100, float(matrix.max()) * 1.20 if matrix.size else 100))
        fig.text(
            0.5,
            0.01,
            "Percentages are normalized within each database across the 12 scored metrics.",
            ha="center",
            fontsize=10,
        )
        fig.suptitle("Runtime Share Across All Scoring Metrics", fontsize=23)

        save_path = save_path or self._runtime_default_save_path("runtime_metric_percentages.png")
        self._runtime_save_figure(fig, save_path, tight_rect=(0.0, 0.03, 1.0, 0.96))

    def plot_category_runtime_stacked_bars(self, figsize=(14, 8), save_path=None):
        runtime_results = self._runtime_results()
        if not runtime_results:
            return

        category_order = self._runtime_category_order(runtime_results)
        db_names, matrix = self._runtime_matrix(runtime_results, "category_seconds", category_order)

        fig, ax = plotting.plt.subplots(figsize=figsize)
        bottom = np.zeros(len(db_names), dtype=float)
        for category_index, category in enumerate(category_order):
            values = matrix[:, category_index]
            ax.bar(
                db_names,
                values,
                bottom=bottom,
                color=RUNTIME_CATEGORY_COLORS[category],
                label=category,
            )
            bottom += values

        ax.set_xlabel("Database", fontsize=21)
        ax.set_ylabel("Runtime (s)", fontsize=21)
        ax.set_title("Category Runtime Composition", fontsize=22)
        ax.tick_params(axis="x", rotation=45, labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        plotting.plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
        ax.legend(
            fontsize=13,
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
        )

        y_offset = max(bottom) * 0.02 if np.max(bottom) > 0 else 0.02
        for index, total_seconds in enumerate(bottom):
            ax.text(index, total_seconds + y_offset, f"{total_seconds:.2f}s", ha="center", va="bottom", fontsize=9)

        fig.suptitle("Database Runtime Composition by Category", fontsize=23)
        save_path = save_path or self._runtime_default_save_path("runtime_category_stacked.png")
        self._runtime_save_figure(fig, save_path, tight_rect=(0.0, 0.0, 1.0, 0.96))

    def plot_category_runtime_stacked_percentage_bars(self, figsize=(14, 8), save_path=None):
        runtime_results = self._runtime_results()
        if not runtime_results:
            return

        category_order = self._runtime_category_order(runtime_results)
        db_names, matrix = self._runtime_matrix(runtime_results, "category_percentages", category_order)

        fig, ax = plotting.plt.subplots(figsize=figsize)
        bottom = np.zeros(len(db_names), dtype=float)
        for category_index, category in enumerate(category_order):
            values = matrix[:, category_index]
            bars = ax.bar(
                db_names,
                values,
                bottom=bottom,
                color=RUNTIME_CATEGORY_COLORS[category],
                label=category,
            )
            for bar, percentage, baseline in zip(bars, values, bottom):
                if percentage >= 6.0:
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        baseline + percentage / 2,
                        f"{percentage:.1f}%",
                        ha="center",
                        va="center",
                        fontsize=9,
                        color="white",
                        fontweight="bold",
                    )
            bottom += values

        ax.set_xlabel("Database", fontsize=21)
        ax.set_ylabel("Runtime Share (%)", fontsize=21)
        ax.set_title("Category Runtime Composition", fontsize=22)
        ax.tick_params(axis="x", rotation=45, labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        plotting.plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
        ax.legend(
            fontsize=13,
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
        )
        ax.set_ylim(0, 100)
        fig.suptitle("Database Runtime Composition by Category", fontsize=23)
        save_path = save_path or self._runtime_default_save_path("runtime_category_stacked_percent.png")
        self._runtime_save_figure(fig, save_path, tight_rect=(0.0, 0.0, 0.86, 0.96))

    def plot_metric_runtime_heatmap(self, figsize=(18, 10), save_path=None):
        runtime_results = self._runtime_results()
        if not runtime_results:
            return

        category_order = self._runtime_category_order(runtime_results)
        metric_order = self._runtime_metric_order(runtime_results, category_order)
        db_names, matrix = self._runtime_matrix(runtime_results, "metric_seconds", metric_order)
        metric_labels = self._runtime_metric_labels(metric_order)

        fig, ax = plotting.plt.subplots(figsize=figsize)
        heatmap = plotting.sns.heatmap(
            matrix,
            annot=True,
            fmt=".2f",
            cmap="YlOrBr",
            linewidths=0.5,
            linecolor="white",
            cbar_kws={"label": "Runtime (s)"},
            xticklabels=metric_labels,
            yticklabels=db_names,
            ax=ax,
        )
        ax.set_xlabel("Metric", fontsize=21)
        ax.set_ylabel("Database", fontsize=21)
        ax.set_title("Metric Runtime Heatmap", fontsize=23)
        ax.tick_params(axis="x", labelrotation=45, labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        plotting.plt.setp(ax.xaxis.get_majorticklabels(), ha="right")
        plotting.plt.setp(ax.yaxis.get_majorticklabels(), rotation=45, ha="right", va="center")
        colorbar = heatmap.collections[0].colorbar
        colorbar.ax.tick_params(labelsize=18)
        colorbar.set_label("Runtime (s)", fontsize=20)

        save_path = save_path or self._runtime_default_save_path("runtime_metric_heatmap.png")
        self._runtime_save_figure(fig, save_path, tight_rect=(0.08, 0.0, 1.0, 1.0))

    def _runtime_db_profile(self, db_name):
        runtime_results = self._runtime_results()
        if not runtime_results:
            return None

        for current_name, results, scorer, profile in runtime_results:
            if current_name == db_name:
                return results, scorer, profile

        print(f"No runtime profile found for database: {db_name}")
        return None

    def plot_runtime_gantt(self, db_name, figsize=(14, 8), save_path=None):
        payload = self._runtime_db_profile(db_name)
        if payload is None:
            return

        results, _, profile = payload
        metric_records = profile.get("metric_records", [])
        if not metric_records:
            print(f"No metric timeline available for database: {db_name}")
            return

        label_map = self._runtime_metric_label_map()
        fig, ax = plotting.plt.subplots(figsize=figsize)
        y_positions = np.arange(len(metric_records))

        for index, record in enumerate(metric_records):
            ax.barh(
                y_positions[index],
                record["seconds"],
                left=record["start_seconds"],
                color=RUNTIME_CATEGORY_COLORS[record["category"]],
                edgecolor="white",
                linewidth=1.0,
            )
            if record["seconds"] > profile["total_seconds"] * 0.04:
                ax.text(
                    record["start_seconds"] + record["seconds"] / 2,
                    y_positions[index],
                    f"{record['seconds']:.2f}s",
                    ha="center",
                    va="center",
                    fontsize=9,
                    color="white",
                    fontweight="bold",
                )

        ax.set_yticks(y_positions)
        ax.set_yticklabels(
            [label_map.get(record["metric"], record["metric"]) for record in metric_records],
            fontsize=14,
        )
        ax.invert_yaxis()
        ax.set_xlabel("Elapsed Time (s)", fontsize=21)
        ax.set_ylabel("Metric", fontsize=21)
        ax.set_title(f"Runtime Gantt Chart: {db_name}", fontsize=22)
        ax.tick_params(axis="x", labelsize=18)
        ax.grid(axis="x", linestyle="--", alpha=0.3)

        legend_handles = [
            plotting.mpatches.Patch(color=RUNTIME_CATEGORY_COLORS[category], label=category)
            for category in profile["category_order"]
        ]
        ax.legend(handles=legend_handles, fontsize=12, frameon=False, loc="lower left")

        save_path = save_path or os.path.join(results["output_dir"], "runtime_gantt.png")
        self._runtime_save_figure(fig, save_path)

    def plot_runtime_pareto(self, db_name, figsize=(14, 8), save_path=None):
        payload = self._runtime_db_profile(db_name)
        if payload is None:
            return

        results, _, profile = payload
        metric_records = sorted(
            profile.get("metric_records", []),
            key=lambda record: record["seconds"],
            reverse=True,
        )
        if not metric_records:
            print(f"No metric runtime data available for database: {db_name}")
            return

        label_map = self._runtime_metric_label_map()
        category_order = profile["category_order"]
        metric_colors = self._runtime_metric_colors(category_order)
        labels = [label_map.get(record["metric"], record["metric"]) for record in metric_records]
        values = np.array([record["seconds"] for record in metric_records], dtype=float)
        cumulative = np.cumsum(values) / values.sum() * 100.0 if values.sum() > 0 else np.zeros_like(values)

        fig, ax = plotting.plt.subplots(figsize=figsize)
        bars = ax.bar(
            labels,
            values,
            color=[metric_colors.get(record["metric"], "#888888") for record in metric_records],
        )
        ax.set_xlabel("Metric", fontsize=21)
        ax.set_ylabel("Runtime (s)", fontsize=21)
        ax.set_title(f"Runtime Pareto Analysis: {db_name}", fontsize=22)
        ax.tick_params(axis="x", rotation=45, labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        plotting.plt.setp(ax.xaxis.get_majorticklabels(), ha="right")

        y_offset = values.max() * 0.02 if values.size and values.max() > 0 else 0.02
        for bar, seconds in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + y_offset,
                f"{seconds:.2f}s",
                ha="center",
                va="bottom",
                fontsize=9,
            )

        ax2 = ax.twinx()
        ax2.plot(np.arange(len(labels)), cumulative, color="#7B1E3A", marker="o", linewidth=2.2)
        ax2.set_ylabel("Cumulative Share (%)", fontsize=21)
        ax2.tick_params(axis="y", labelsize=18)
        ax2.set_ylim(0, 105)

        fig.suptitle(f"Metric Runtime Pareto Chart: {db_name}", fontsize=23)
        save_path = save_path or os.path.join(results["output_dir"], "runtime_pareto.png")
        self._runtime_save_figure(fig, save_path, tight_rect=(0.0, 0.0, 1.0, 0.96))

    def plot_runtime_sunburst(self, db_name, figsize=(12, 12), save_path=None):
        payload = self._runtime_db_profile(db_name)
        if payload is None:
            return

        results, _, profile = payload
        category_order = profile["category_order"]
        metric_order = self._runtime_metric_order([(db_name, results, None, profile)], category_order)
        metric_label_map = self._runtime_metric_label_map()
        metric_colors = self._runtime_metric_colors(category_order)

        category_values = [profile["category_seconds"].get(category, 0.0) for category in category_order]
        outer_values = []
        outer_labels = []
        outer_colors = []
        for category in category_order:
            for metric in self.categories.get(category, []):
                if metric not in metric_order:
                    continue
                value = profile["metric_seconds"].get(metric, 0.0)
                if value <= 0:
                    continue
                outer_values.append(value)
                percentage = profile["metric_percentages"].get(metric, 0.0)
                outer_labels.append(metric_label_map.get(metric, metric) if percentage >= 4 else "")
                outer_colors.append(metric_colors.get(metric, RUNTIME_CATEGORY_COLORS[category]))

        if sum(category_values) <= 0 or not outer_values:
            print(f"No runtime hierarchy data available for database: {db_name}")
            return

        fig, ax = plotting.plt.subplots(figsize=figsize)
        ax.pie(
            category_values,
            radius=1.0,
            labels=[self._runtime_category_label(category) for category in category_order],
            colors=[RUNTIME_CATEGORY_COLORS[category] for category in category_order],
            labeldistance=0.68,
            textprops={"fontsize": 12, "fontweight": "bold"},
            wedgeprops={"width": 0.28, "edgecolor": "white"},
        )
        ax.pie(
            outer_values,
            radius=1.33,
            labels=outer_labels,
            colors=outer_colors,
            labeldistance=1.04,
            textprops={"fontsize": 10},
            wedgeprops={"width": 0.28, "edgecolor": "white"},
        )
        ax.text(
            0,
            0,
            f"{profile['total_seconds']:.2f}s\nTotal",
            ha="center",
            va="center",
            fontsize=18,
            fontweight="bold",
        )
        ax.set_title(f"Runtime Sunburst: {db_name}", fontsize=23, pad=18)
        ax.set_aspect("equal")

        save_path = save_path or os.path.join(results["output_dir"], "runtime_sunburst.png")
        self._runtime_save_figure(fig, save_path)

    def _runtime_slice_rectangles(self, values, x, y, width, height, horizontal):
        total = sum(values)
        rectangles = []
        if total <= 0:
            return rectangles

        offset = 0.0
        for value in values:
            fraction = value / total if total > 0 else 0.0
            if horizontal:
                rect_width = width * fraction
                rectangles.append((x + offset, y, rect_width, height))
                offset += rect_width
            else:
                rect_height = height * fraction
                rectangles.append((x, y + offset, width, rect_height))
                offset += rect_height
        return rectangles

    def plot_runtime_treemap(self, db_name, figsize=(14, 8), save_path=None):
        payload = self._runtime_db_profile(db_name)
        if payload is None:
            return

        results, _, profile = payload
        category_order = [
            category
            for category in profile["category_order"]
            if profile["category_seconds"].get(category, 0.0) > 0
        ]
        if not category_order:
            print(f"No runtime hierarchy data available for database: {db_name}")
            return

        metric_colors = self._runtime_metric_colors(category_order)
        fig, ax = plotting.plt.subplots(figsize=figsize)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

        category_values = [profile["category_seconds"][category] for category in category_order]
        category_rects = self._runtime_slice_rectangles(category_values, 0.0, 0.0, 1.0, 1.0, horizontal=True)

        for category, rect in zip(category_order, category_rects):
            x_pos, y_pos, rect_width, rect_height = rect
            category_patch = plotting.Rectangle(
                (x_pos, y_pos),
                rect_width,
                rect_height,
                facecolor=self._lighten_color(RUNTIME_CATEGORY_COLORS[category], 0.55),
                edgecolor="white",
                linewidth=2.0,
            )
            ax.add_patch(category_patch)

            metrics = [
                metric
                for metric in self.categories.get(category, [])
                if profile["metric_seconds"].get(metric, 0.0) > 0
            ]
            metric_values = [profile["metric_seconds"][metric] for metric in metrics]
            inner_horizontal = rect_width < rect_height
            metric_rects = self._runtime_slice_rectangles(
                metric_values,
                x_pos,
                y_pos,
                rect_width,
                rect_height,
                horizontal=inner_horizontal,
            )

            for metric, metric_rect in zip(metrics, metric_rects):
                mx, my, mwidth, mheight = metric_rect
                metric_patch = plotting.Rectangle(
                    (mx, my),
                    mwidth,
                    mheight,
                    facecolor=metric_colors.get(metric, RUNTIME_CATEGORY_COLORS[category]),
                    edgecolor="white",
                    linewidth=1.5,
                )
                ax.add_patch(metric_patch)
                area = mwidth * mheight
                if area >= 0.035:
                    ax.text(
                        mx + mwidth / 2,
                        my + mheight / 2,
                        f"{metric}\n{profile['metric_seconds'][metric]:.2f}s",
                        ha="center",
                        va="center",
                        fontsize=10,
                        color="white",
                        fontweight="bold",
                    )

            if rect_width * rect_height >= 0.09:
                ax.text(
                    x_pos + 0.01,
                    y_pos + rect_height - 0.02,
                    f"{category}\n{profile['category_seconds'][category]:.2f}s",
                    ha="left",
                    va="top",
                    fontsize=12,
                    color="black",
                    fontweight="bold",
                )

        ax.set_title(f"Runtime Treemap: {db_name}", fontsize=23, pad=12)
        save_path = save_path or os.path.join(results["output_dir"], "runtime_treemap.png")
        self._runtime_save_figure(fig, save_path)

    def generate_runtime_analysis_plots(self):
        runtime_results = self._runtime_results()
        if not runtime_results:
            return

        self.plot_total_runtime_line()
        self.plot_total_runtime_bar()
        self.plot_category_runtime_percentage_bars()
        self.plot_metric_runtime_percentage_bars()
        self.plot_category_runtime_stacked_bars()
        self.plot_category_runtime_stacked_percentage_bars()
        self.plot_metric_runtime_heatmap()

        for db_name, _, _, _ in runtime_results:
            self.plot_runtime_gantt(db_name)
            self.plot_runtime_pareto(db_name)
            self.plot_runtime_sunburst(db_name)
            self.plot_runtime_treemap(db_name)
