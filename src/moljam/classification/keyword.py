"""
Column classifier based on column-name keywords.

This module now exposes a structured role vote in addition to the legacy
"useful/excluded/unknown" category mapping.
"""

import re

import pandas as pd


class KeywordColumnClassifier:
    """Column classifier driven by semantic signals in column names."""

    POSITIVE_ROLES = {"time", "activity", "label", "experimental_context"}

    def __init__(self, uniqueness_threshold=0.98):
        self.uniqueness_threshold = uniqueness_threshold
        self.min_categories = 2
        self.max_categories = 50
        self.cv_threshold = 0.1

        self._init_keywords()

    def _init_keywords(self):
        self.id_keywords = ["id", "key", "index", "number", "code"]
        self.id_exceptions = ["assay", "experiment", "target", "document", "cell", "bao", "format"]

        self.name_keywords = ["name", "title"]
        self.name_exceptions = ["assay", "target", "cell", "tissue"]

        self.time_keywords = [
            "time",
            "timepoint",
            "time point",
            "hour",
            "hours",
            "hr",
            "hrs",
            "day",
            "days",
            "week",
            "weeks",
            "month",
            "months",
            "duration",
        ]

        self.calculable_descriptors = [
            "molecular weight",
            "mw",
            "exact mass",
            "logp",
            "alogp",
            "xlogp",
            "mlogp",
            "tpsa",
            "polar surface area",
            "hbd",
            "h-bond donor",
            "hydrogen bond donor",
            "hba",
            "h-bond acceptor",
            "hydrogen bond acceptor",
            "rotatable bond",
            "rotation",
            "ring",
            "aromatic",
            "ro5",
            "lipinski",
            "violation",
            "degree",
            "atom count",
            "heavy atom",
            "ligand efficiency",
        ]

        self.relation_keywords = ["relation", "operator", "comparator"]
        self.metadata_keywords = ["journal", "year", "source description", "document", "comment", "properties"]

        self.experimental_keywords = [
            "assay",
            "experiment",
            "method",
            "protocol",
            "technique",
            "organism",
            "species",
            "strain",
            "tissue",
            "cell",
            "target",
            "enzyme",
            "condition",
            "temperature",
            "buffer",
            "bao",
            "bioassay",
            "screening",
            "variant",
            "mutation",
            "subcellular",
            "fraction",
            "medium",
        ]
        self.experimental_weak_keywords = ["parameter", "matrix"]

        self.activity_pattern_keywords = [
            (r"(?<![A-Za-z0-9])p?ic50(?![A-Za-z0-9])", "ic50"),
            (r"(?<![A-Za-z0-9])p?ec50(?![A-Za-z0-9])", "ec50"),
            (r"(?<![A-Za-z0-9])p?ed50(?![A-Za-z0-9])", "ed50"),
            (r"(?<![A-Za-z0-9])p?ki(?![A-Za-z0-9])", "ki"),
            (r"(?<![A-Za-z0-9])p?kd(?![A-Za-z0-9])", "kd"),
        ]

        self.activity_keywords = [
            "activity",
            "potency",
            "inhibition",
            "binding",
            "solubility",
            "measured",
            "observed",
            "expt",
            "concentration",
            "dose",
            "response",
            "pchembl",
            "affinity",
        ]

        self.prediction_keywords = [
            "calc",
            "calculated",
            "predicted",
            "pred",
            "estimate",
            "estimated",
            "computed",
            "simulation",
            "model",
        ]

        self.label_keywords = [
            "class",
            "type",
            "category",
            "group",
            "status",
            "label",
            "flag",
            "disorder",
            "disease",
            "issue",
            "outcome",
            "active",
            "inactive",
        ]

        self.ph_keywords = ["ph", "ph value", "ph level"]

    def _keyword_match(self, text, keywords, use_word_boundary=True, exceptions=None):
        text_lower = text.lower()

        for keyword in keywords:
            if use_word_boundary:
                pattern = r"(?<![A-Za-z0-9])" + re.escape(keyword) + r"(?![A-Za-z0-9])"

                if re.search(pattern, text_lower):
                    if exceptions and any(exc in text_lower for exc in exceptions):
                        continue
                    return True, keyword
            else:
                if keyword in text_lower:
                    if exceptions and any(exc in text_lower for exc in exceptions):
                        continue
                    return True, keyword

        return False, None

    def _match_activity_keyword(self, text):
        text_lower = text.lower()

        for pattern, keyword in self.activity_pattern_keywords:
            if re.search(pattern, text_lower):
                return True, keyword

        if "muv-" in text_lower:
            return True, "muv-"

        return self._keyword_match(text_lower, self.activity_keywords, use_word_boundary=True)

    def _match_experimental_keyword(self, text):
        text_lower = text.lower()

        matched, keyword = self._keyword_match(text_lower, self.experimental_keywords, use_word_boundary=True)
        if matched:
            return True, keyword

        has_strong_anchor = any(
            self._keyword_match(text_lower, [anchor], use_word_boundary=True)[0]
            for anchor in self.experimental_keywords
        )
        if not has_strong_anchor:
            return False, None

        return self._keyword_match(text_lower, self.experimental_weak_keywords, use_word_boundary=True)

    def looks_like_measurement_name(self, col_name):
        col_lower = str(col_name).lower().strip()

        if self._match_activity_keyword(col_lower)[0]:
            return True

        generic_measurement_names = {
            "value",
            "result",
            "measurement",
            "readout",
            "endpoint",
            "exp",
            "expt",
        }
        if col_lower in generic_measurement_names:
            return True

        measurement_phrases = [
            "standard value",
            "assay value",
            "activity value",
            "measured value",
            "observed value",
            "assay result",
            "activity result",
            "measurement value",
            "measurement result",
        ]
        return any(phrase in col_lower for phrase in measurement_phrases)

    def _is_float_measurement(self, col_data):
        if not pd.api.types.is_numeric_dtype(col_data):
            return False

        valid_data = col_data.dropna()
        if len(valid_data) == 0 or not pd.api.types.is_float_dtype(col_data):
            return False

        try:
            return (valid_data % 1 != 0).any()
        except Exception:
            return False

    def _role_to_legacy_category(self, role):
        if role in self.POSITIVE_ROLES:
            return "useful"
        if role in {"excluded", "derived_or_predicted"}:
            return "excluded"
        return "unknown"

    def _build_result(self, role, reason, confidence, col_name, col_data):
        return {
            "role": role,
            "reason": reason,
            "confidence": float(confidence),
            "details": self._get_column_details(col_name, col_data),
        }

    def is_excluded_column(self, col_name, col_data):
        valid_data = col_data.dropna()
        if len(valid_data) == 0:
            return True, "All values are missing", 1.0

        col_lower = col_name.lower()
        uniqueness_ratio = col_data.nunique(dropna=True) / len(col_data) if len(col_data) > 0 else 0

        if uniqueness_ratio > 0.95:
            if self._is_float_measurement(col_data):
                if uniqueness_ratio <= self.uniqueness_threshold:
                    return False, None, 0.0
            else:
                return True, f"Unique identifier (cardinality: {uniqueness_ratio:.1%})", 0.9

        matched, keyword = self._keyword_match(
            col_lower,
            self.id_keywords,
            use_word_boundary=True,
            exceptions=self.id_exceptions,
        )
        if matched:
            return True, f"Identifier column ('{keyword}')", 0.8

        matched, keyword = self._keyword_match(
            col_lower,
            self.name_keywords,
            use_word_boundary=True,
            exceptions=self.name_exceptions,
        )
        if matched:
            if uniqueness_ratio < 0.5:
                return False, None, 0.0
            return True, f"Name column ('{keyword}')", 0.7

        for descriptor in self.calculable_descriptors:
            if descriptor in col_lower:
                return True, f"Calculable descriptor ('{descriptor}')", 0.95

        matched, keyword = self._keyword_match(col_lower, self.relation_keywords, use_word_boundary=True)
        if matched:
            return True, f"Relation/Operator column ('{keyword}')", 0.9

        matched, keyword = self._keyword_match(col_lower, self.metadata_keywords, use_word_boundary=False)
        if matched:
            return True, f"Metadata column ('{keyword}')", 0.8

        if "duplicate" in col_lower:
            return True, "Duplicate flag column", 0.9

        return False, None, 0.0

    def is_useful_column(self, col_name, col_data):
        result = self.classify_role(col_name, col_data)
        is_useful = result["role"] in self.POSITIVE_ROLES
        if not is_useful:
            return False, None, 0.0
        return True, result["reason"], result["confidence"]

    def classify_role(self, col_name, col_data):
        col_lower = str(col_name).lower()

        matched, keyword = self._keyword_match(col_lower, self.prediction_keywords, use_word_boundary=True)
        is_excluded, reason, confidence = self.is_excluded_column(col_name, col_data)
        if matched and (not is_excluded or "Calculable descriptor" in (reason or "")):
            return self._build_result(
                "derived_or_predicted",
                f"Derived or predicted value ('{keyword}')",
                0.9,
                col_name,
                col_data,
            )
        if is_excluded:
            return self._build_result("excluded", reason, confidence, col_name, col_data)

        matched, keyword = self._keyword_match(col_lower, self.time_keywords, use_word_boundary=True)
        if matched:
            return self._build_result("time", f"Time annotation ('{keyword}')", 0.88, col_name, col_data)

        matched, keyword = self._keyword_match(col_lower, self.ph_keywords, use_word_boundary=True)
        if matched:
            return self._build_result(
                "experimental_context",
                f"Experimental condition ('{keyword}')",
                0.9,
                col_name,
                col_data,
            )

        matched, keyword = self._match_experimental_keyword(col_lower)
        if matched:
            return self._build_result(
                "experimental_context",
                f"Experimental context ('{keyword}')",
                0.9,
                col_name,
                col_data,
            )

        matched, keyword = self._match_activity_keyword(col_lower)
        if matched:
            return self._build_result("activity", f"Activity or measurement ('{keyword}')", 0.85, col_name, col_data)

        n_unique = col_data.nunique(dropna=True)
        matched, keyword = self._keyword_match(col_lower, self.label_keywords, use_word_boundary=True)
        if matched and self.min_categories <= n_unique <= self.max_categories:
            return self._build_result(
                "label",
                f"Classification label ('{keyword}', {n_unique} categories)",
                0.8,
                col_name,
                col_data,
            )
        if matched:
            return self._build_result(
                "label",
                f"Label-like column ('{keyword}') without stable category count",
                0.65,
                col_name,
                col_data,
            )

        return self._build_result("unknown", "No matching keyword rules", 0.0, col_name, col_data)

    def classify_column(self, col_name, col_data):
        role_result = self.classify_role(col_name, col_data)
        return {
            "category": self._role_to_legacy_category(role_result["role"]),
            "reason": role_result["reason"],
            "confidence": role_result["confidence"],
            "details": role_result["details"],
            "role": role_result["role"],
        }

    def _get_column_details(self, col_name, col_data):
        n_unique = col_data.nunique(dropna=True)
        n_total = len(col_data)
        if n_total == 0:
            return {
                "unique_values": 0,
                "uniqueness_ratio": 0.0,
                "missing_rate": 0.0,
                "coverage": 0.0,
                "dtype": str(col_data.dtype),
                "is_numeric": pd.api.types.is_numeric_dtype(col_data),
                "is_float": pd.api.types.is_float_dtype(col_data),
            }

        missing_rate = col_data.isnull().sum() / n_total
        coverage = 1 - missing_rate

        return {
            "unique_values": n_unique,
            "uniqueness_ratio": n_unique / n_total,
            "missing_rate": missing_rate,
            "coverage": coverage,
            "dtype": str(col_data.dtype),
            "is_numeric": pd.api.types.is_numeric_dtype(col_data),
            "is_float": pd.api.types.is_float_dtype(col_data),
        }

    def classify_columns(self, df, smiles_col="Smiles"):
        results = {"useful": [], "excluded": [], "unknown": []}

        smiles_col_actual = None
        for col in df.columns:
            if col.lower() == smiles_col.lower():
                smiles_col_actual = col
                break

        for col in df.columns:
            if smiles_col_actual and col == smiles_col_actual:
                continue

            result = self.classify_column(col, df[col])
            entry = (col, result["reason"], result["confidence"], result["details"])
            results[result["category"]].append(entry)

        return results
