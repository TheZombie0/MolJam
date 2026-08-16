from .._common import *


class AnnotationQualityChecksMixin:
    ANNOTATION_SUPPORT_KEY = "Annotation Support Quality"
    CONTEXT_SUPPORT_THRESHOLD = 0.6
    CONTEXT_SUPPORT_ABUNDANCE_SCALE = 4.0
    CONTEXT_FAMILY_COVERAGE_SATURATION = 0.3
    CONTEXT_FAMILY_BREADTH_THRESHOLD = 0.4
    CORE_CONTEXT_FAMILIES = ("assay", "target", "biosystem", "protocol")
    OPTIONAL_CONTEXT_FAMILIES = ("variant",)

    def _resolve_annotation_column(self, column_name):
        if not column_name:
            return None
        if column_name in self.df.columns:
            return column_name

        lowercase_map = {str(col).lower(): col for col in self.df.columns}
        return lowercase_map.get(str(column_name).lower())

    def _build_explicit_annotation_roles(self):
        priorities = {}
        explicit_roles = {}

        def assign(column_name, role, priority):
            resolved_name = self._resolve_annotation_column(column_name)
            if resolved_name is None:
                return
            if priorities.get(resolved_name, -1) <= priority:
                priorities[resolved_name] = priority
                explicit_roles[resolved_name] = role

        assign(self.id_col, "excluded", 100)
        assign(self.name_col, "excluded", 100)
        assign(self.time_col, "time", 90)

        for column_name in self.activity_cols:
            assign(column_name, "activity", 80)

        for column_name in self.class_cols:
            assign(column_name, "label", 70)

        for column_name in self.experimental_method_cols:
            assign(column_name, "experimental_context", 60)

        return explicit_roles

    def _score_from_signal(self, quality_signal, max_score=10):
        bounded_signal = max(0.0, min(1.0, float(quality_signal)))
        error_rate = 100 * (1 - bounded_signal)
        return self.calculate_quality_score(error_rate, max_score=max_score), error_rate

    def _column_support_signal(self, column):
        return max(
            0.0,
            min(
                1.0,
                column["coverage"] * (0.8 * column["confidence"] + 0.2 * column["info_strength"]),
            ),
        )

    def _classify_context_family(self, column_name):
        column_name = str(column_name).lower()
        family_keywords = (
            ("variant", ("variant", "mutation", "accession", "genotype")),
            ("biosystem", ("cell", "tissue", "organism", "subcellular", "fraction", "species", "strain")),
            ("protocol", ("bao", "format", "parameter", "protocol", "method", "condition", "readout", "endpoint", "technology")),
            ("target", ("target",)),
            ("assay", ("assay",)),
        )

        for family_name, keywords in family_keywords:
            if any(keyword in column_name for keyword in keywords):
                return family_name

        return "other"

    def _context_family_presence_signal(self, columns):
        if not columns:
            return 0.0

        max_coverage = max(float(column["coverage"]) for column in columns)
        if max_coverage <= 0:
            return 0.0

        return float(
            max(
                0.0,
                min(
                    1.0,
                    np.sqrt(max_coverage / self.CONTEXT_FAMILY_COVERAGE_SATURATION),
                ),
            )
        )

    def _column_summary(self, column):
        return {
            "name": column["column"],
            "role": column["role"],
            "accepted": column["accepted"],
            "confidence": float(round(column["confidence"], 4)),
            "coverage": float(round(column["coverage"], 4)),
            "info_strength": float(round(column["info_strength"], 4)),
            "agreement": column["agreement"],
            "conflicted": column["conflicted"],
            "decision_type": column["decision_type"],
            "reason": column["reason"],
        }

    def _calculate_time_label_metrics(self):
        resolved_time_col = self._resolve_annotation_column(self.time_col)
        if resolved_time_col and resolved_time_col in self.df.columns:
            time_series = self.df[resolved_time_col]
            time_coverage = time_series.notna().sum() / self.num_molecules if self.num_molecules > 0 else 0.0
            unique_times = time_series.dropna().unique()
            num_time_points = len(unique_times)

            richness_component = min(1.0, np.log2(1 + num_time_points) / np.log2(1 + 4)) if num_time_points > 0 else 0.0
            quality_signal = 0.8 * time_coverage + 0.2 * richness_component
            score, error_rate = self._score_from_signal(quality_signal)

            info = {
                "Time label column": resolved_time_col,
                "Time coverage rate": f"{time_coverage * 100:.2f}%",
                "Number of unique time points": num_time_points,
                "Example time points": list(unique_times[:10]),
                "Coverage component": round(time_coverage, 4),
                "Richness component": round(richness_component, 4),
                "Quality signal": round(quality_signal, 4),
                "Error rate": f"{error_rate:.2f}%",
            }
            return score, info, quality_signal

        info = {
            "Time label column": self.time_col if self.time_col else "None",
            "Time coverage rate": "0.00%",
            "Coverage component": 0.0,
            "Richness component": 0.0,
            "Quality signal": 0.0,
            "Error rate": "100.00%",
            "Note": "No time label column provided or column does not exist",
        }
        return 0.0, info, 0.0

    def _calculate_annotation_support_from_context(self, context_columns):
        context_signal_details = [
            {
                "column": column["column"],
                "signal": round(self._column_support_signal(column), 4),
                "coverage": round(float(column["coverage"]), 4),
                "confidence": round(float(column["confidence"]), 4),
            }
            for column in context_columns
        ]
        qualifying_columns = [
            detail
            for detail in context_signal_details
            if detail["signal"] > self.CONTEXT_SUPPORT_THRESHOLD
        ]

        if not qualifying_columns:
            return {
                "quality_signal": 0.0,
                "strength_component": 0.0,
                "abundance_component": 0.0,
                "qualifying_count": 0,
                "qualified_columns": [],
                "all_context_columns": context_signal_details,
            }

        strength_component = float(np.mean([detail["signal"] for detail in qualifying_columns]))
        abundance_component = float(1 - np.exp(-len(qualifying_columns) / self.CONTEXT_SUPPORT_ABUNDANCE_SCALE))
        quality_signal = float(0.6 * strength_component + 0.4 * abundance_component)

        return {
            "quality_signal": quality_signal,
            "strength_component": strength_component,
            "abundance_component": abundance_component,
            "qualifying_count": len(qualifying_columns),
            "qualified_columns": qualifying_columns,
            "all_context_columns": context_signal_details,
        }

    def _calculate_context_breadth_signal(self, context_columns):
        family_columns = {
            family_name: []
            for family_name in self.CORE_CONTEXT_FAMILIES + self.OPTIONAL_CONTEXT_FAMILIES + ("other",)
        }
        for column in context_columns:
            family_columns[self._classify_context_family(column["column"])].append(column)

        family_presences = {
            family_name: round(self._context_family_presence_signal(columns), 4)
            for family_name, columns in family_columns.items()
        }
        core_presences = {
            family_name: family_presences[family_name]
            for family_name in self.CORE_CONTEXT_FAMILIES
        }
        present_core_count = sum(
            1
            for signal in core_presences.values()
            if signal >= self.CONTEXT_FAMILY_BREADTH_THRESHOLD
        )
        breadth_ratio = (
            present_core_count / len(self.CORE_CONTEXT_FAMILIES)
            if self.CORE_CONTEXT_FAMILIES
            else 0.0
        )
        present_core_signals = [signal for signal in core_presences.values() if signal > 0]
        coverage_depth = float(np.mean(present_core_signals)) if present_core_signals else 0.0
        variant_bonus = 0.05 * family_presences["variant"]
        quality_signal = float(
            min(1.0, 0.45 * breadth_ratio + 0.55 * coverage_depth + variant_bonus)
        )

        family_details = {
            family_name: {
                "columns": [column["column"] for column in columns],
                "max_coverage": round(
                    max((float(column["coverage"]) for column in columns), default=0.0),
                    4,
                ),
                "presence_signal": family_presences[family_name],
            }
            for family_name, columns in family_columns.items()
        }

        return {
            "quality_signal": quality_signal,
            "breadth_ratio": round(breadth_ratio, 4),
            "coverage_depth": round(coverage_depth, 4),
            "variant_bonus": round(variant_bonus, 4),
            "family_presences": family_presences,
            "present_core_count": present_core_count,
            "family_details": family_details,
        }

    def check_annotation_quality(self):
        """
        Score experimental information quality from a structured annotation profile.

        Time labels keep the original coverage/richness logic. The other 2 metrics
        focus on experimental context specifically:
        - Annotation Support Quality measures strong experimental-context support.
        - Type Diversity measures breadth across experimental-context families.
        """
        if not self.experimental_info:
            print("Annotation quality not included in scoring, skipping")
            self.completed_checks.add("check_annotation_quality")
            return 0

        total_score = 0.0
        time_label_elapsed = 0.0
        support_phase_elapsed = 0.0
        diversity_phase_elapsed = 0.0

        time_label_start = time.perf_counter()
        time_label_score, time_label_info, time_signal = self._calculate_time_label_metrics()
        time_label_elapsed = time.perf_counter() - time_label_start
        self.scores["Experimental Information Quality"]["Time Label Availability"] = time_label_score
        total_score += time_label_score

        try:
            profile_start = time.perf_counter()
            from ...classification import HybridColumnClassifier

            classifier = HybridColumnClassifier()
            explicit_roles = self._build_explicit_annotation_roles()
            profile = classifier.build_annotation_profile(
                self.df,
                smiles_col=self.smiles_col,
                explicit_roles=explicit_roles,
            )

            accepted_columns = profile["accepted_columns"]
            role_counts = profile["summary"]["role_counts"]
            accepted_annotation_columns = [
                column
                for column in accepted_columns
                if column["role"] in {"activity", "label", "experimental_context"}
            ]
            context_columns = [
                column
                for column in accepted_columns
                if column["role"] == "experimental_context"
            ]

            if accepted_annotation_columns:
                avg_confidence = sum(column["confidence"] for column in accepted_annotation_columns) / len(accepted_annotation_columns)
                agreement_rate = sum(1 for column in accepted_annotation_columns if column["agreement"]) / len(accepted_annotation_columns)
                conflict_rate = (
                    sum(
                        1
                        for column in accepted_annotation_columns
                        if column["conflicted"] or column["confidence"] < 0.6
                    )
                    / len(accepted_annotation_columns)
                )
            else:
                avg_confidence = 0.0
                agreement_rate = 0.0
                conflict_rate = 0.0

            shared_profile_elapsed = time.perf_counter() - profile_start

            support_calc_start = time.perf_counter()
            support_metrics = self._calculate_annotation_support_from_context(context_columns)
            support_signal = support_metrics["quality_signal"]
            annotation_support_quality, support_error_rate = self._score_from_signal(support_signal)
            support_calc_elapsed = time.perf_counter() - support_calc_start

            diversity_calc_start = time.perf_counter()
            diversity_metrics = self._calculate_context_breadth_signal(context_columns)
            diversity_signal = diversity_metrics["quality_signal"]
            diversity_score, diversity_error_rate = self._score_from_signal(diversity_signal)
            diversity_calc_elapsed = time.perf_counter() - diversity_calc_start

            self.scores["Experimental Information Quality"][self.ANNOTATION_SUPPORT_KEY] = float(annotation_support_quality)
            self.scores["Experimental Information Quality"]["Type Diversity"] = float(diversity_score)
            total_score += annotation_support_quality + diversity_score

            accepted_details = [
                self._column_summary(column)
                for column in sorted(
                    accepted_columns,
                    key=lambda item: (item["role"], -item["confidence"], -item["coverage"]),
                )
            ]
            context_details = [
                self._column_summary(column)
                for column in sorted(
                    context_columns,
                    key=lambda item: (-item["coverage"], -item["confidence"], item["column"]),
                )
            ]
            profile_preview = [
                self._column_summary(column)
                for column in profile["columns"][:20]
            ]

            types_found = [
                role_name
                for role_name in ["time", "activity", "label", "experimental_context"]
                if role_counts.get(role_name, 0) > 0
            ]

            self.analysis_results["Experimental Information Quality"] = {
                "Time Label": time_label_info,
                "Total columns analyzed": profile["summary"]["total_columns"],
                "Accepted columns": profile["summary"]["accepted_count"],
                "Support columns": len(context_columns),
                "Qualified support columns": support_metrics["qualifying_count"],
                "Excluded columns": role_counts.get("excluded", 0),
                "Derived/Predicted columns": role_counts.get("derived_or_predicted", 0),
                "Unknown columns": role_counts.get("unknown", 0),
                "Role counts": role_counts,
                "Average accepted coverage": f"{profile['summary']['avg_coverage'] * 100:.2f}%",
                "Average accepted confidence": f"{profile['summary']['avg_confidence']:.2f}",
                "Agreement rate": f"{agreement_rate * 100:.2f}%",
                "Conflict rate": f"{conflict_rate * 100:.2f}%",
                "Types found": types_found,
                "Metric signals": {
                    "Time Label Availability": {
                        "quality_signal": round(time_signal, 4),
                        "score": f"{time_label_score:.2f}/10",
                    },
                    self.ANNOTATION_SUPPORT_KEY: {
                        "quality_signal": round(support_signal, 4),
                        "signal_threshold": round(self.CONTEXT_SUPPORT_THRESHOLD, 4),
                        "strength_component": round(support_metrics["strength_component"], 4),
                        "abundance_component": round(support_metrics["abundance_component"], 4),
                        "qualifying_count": support_metrics["qualifying_count"],
                        "qualified_columns": support_metrics["qualified_columns"],
                        "all_context_columns": support_metrics["all_context_columns"],
                        "score": f"{annotation_support_quality:.2f}/10",
                        "error_rate": f"{support_error_rate:.2f}%",
                    },
                    "Type Diversity": {
                        "quality_signal": round(diversity_signal, 4),
                        "breadth_ratio": diversity_metrics["breadth_ratio"],
                        "coverage_depth": diversity_metrics["coverage_depth"],
                        "variant_bonus": diversity_metrics["variant_bonus"],
                        "family_presences": diversity_metrics["family_presences"],
                        "present_core_count": diversity_metrics["present_core_count"],
                        "family_details": diversity_metrics["family_details"],
                        "score": f"{diversity_score:.2f}/10",
                        "error_rate": f"{diversity_error_rate:.2f}%",
                    },
                },
                "Accepted column details": accepted_details[:20],
                "Support column details": context_details[:20],
                "Column profile preview": profile_preview,
            }

            print(
                "Annotation quality analysis: "
                f"time={time_label_score:.2f}/10, "
                f"support={annotation_support_quality:.2f}/10, "
                f"diversity={diversity_score:.2f}/10, "
                f"accepted={profile['summary']['accepted_count']}, "
                f"total={total_score:.2f}/30"
            )

            support_phase_elapsed = shared_profile_elapsed * 0.5 + support_calc_elapsed
            diversity_phase_elapsed = shared_profile_elapsed * 0.5 + diversity_calc_elapsed

        except ImportError as error:
            print(f"Warning: Could not import HybridColumnClassifier: {error}")
            self.scores["Experimental Information Quality"][self.ANNOTATION_SUPPORT_KEY] = 0
            self.scores["Experimental Information Quality"]["Type Diversity"] = 0
            self.analysis_results["Experimental Information Quality"] = {
                "Time Label": time_label_info,
                "Note": "Column classification skipped due to import error",
            }
            support_phase_elapsed = 0.0
            diversity_phase_elapsed = 0.0
        except Exception as error:
            print(f"Error during annotation quality check: {error}")
            self.scores["Experimental Information Quality"][self.ANNOTATION_SUPPORT_KEY] = 0
            self.scores["Experimental Information Quality"]["Type Diversity"] = 0
            self.analysis_results["Experimental Information Quality"] = {
                "Time Label": time_label_info,
                "Note": f"Annotation profile generation failed: {error}",
            }
            support_phase_elapsed = 0.0
            diversity_phase_elapsed = 0.0

        self._queue_runtime_category_allocations(
            "Experimental Information Quality",
            [
                ("Time Label Availability", time_label_elapsed),
                (self.ANNOTATION_SUPPORT_KEY, support_phase_elapsed),
                ("Type Diversity", diversity_phase_elapsed),
            ],
        )

        self.completed_checks.add("check_annotation_quality")
        return total_score
