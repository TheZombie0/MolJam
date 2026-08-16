"""
Pure statistical column classifier.

The classifier now emits structured role votes that can be merged with
keyword signals into a richer annotation profile.
"""

import re

import numpy as np
import pandas as pd


class StatisticalColumnClassifier:
    """Classifier based on statistical features only."""

    POSITIVE_ROLES = {"time", "activity", "label", "experimental_context"}

    def __init__(self):
        self.features_cache = {}

    def extract_features(self, col_name, col_data):
        n_total = len(col_data)
        n_valid = col_data.notna().sum()
        n_missing = n_total - n_valid

        features = {
            "n_total": n_total,
            "n_valid": n_valid,
            "n_missing": n_missing,
            "missing_rate": n_missing / n_total if n_total > 0 else 0,
            "n_unique": col_data.nunique(dropna=True),
            "uniqueness_ratio": col_data.nunique(dropna=True) / n_total if n_total > 0 else 0,
            "dtype": str(col_data.dtype),
            "is_numeric": pd.api.types.is_numeric_dtype(col_data),
            "is_integer": pd.api.types.is_integer_dtype(col_data),
            "is_float": pd.api.types.is_float_dtype(col_data),
            "is_object": pd.api.types.is_object_dtype(col_data),
        }

        if features["is_numeric"] and n_valid > 0:
            valid_data = col_data.dropna()

            features["mean"] = valid_data.mean()
            features["median"] = valid_data.median()
            features["std"] = valid_data.std()
            features["min"] = valid_data.min()
            features["max"] = valid_data.max()
            features["range"] = features["max"] - features["min"]
            features["cv"] = features["std"] / abs(features["mean"]) if features["mean"] not in (0, None) else 0

            try:
                features["skewness"] = valid_data.skew()
                features["kurtosis"] = valid_data.kurtosis()
            except Exception:
                features["skewness"] = 0
                features["kurtosis"] = 0

            if features["is_float"]:
                features["has_decimal"] = (valid_data % 1 != 0).any()
                if features["has_decimal"]:
                    decimal_counts = valid_data.apply(
                        lambda x: len(str(x).split(".")[-1]) if "." in str(x) else 0
                    )
                    features["avg_decimal_places"] = decimal_counts.mean()
                else:
                    features["avg_decimal_places"] = 0
            else:
                features["has_decimal"] = False
                features["avg_decimal_places"] = 0

            features["has_negative"] = (valid_data < 0).any()
            features["is_monotonic_increasing"] = bool(
                features["is_integer"] and n_valid > 1 and valid_data.is_monotonic_increasing
            )
        else:
            features.update(
                {
                    "mean": None,
                    "median": None,
                    "std": None,
                    "min": None,
                    "max": None,
                    "range": None,
                    "cv": None,
                    "skewness": None,
                    "kurtosis": None,
                    "has_decimal": False,
                    "avg_decimal_places": 0,
                    "has_negative": False,
                    "is_monotonic_increasing": False,
                }
            )

        if features["is_object"] and n_valid > 0:
            valid_data = col_data.dropna().astype(str)
            str_lengths = valid_data.str.len()
            features["avg_str_length"] = str_lengths.mean()
            features["std_str_length"] = str_lengths.std()
            features["min_str_length"] = str_lengths.min()
            features["max_str_length"] = str_lengths.max()
            features["is_fixed_length"] = features["std_str_length"] == 0

            prefix_number_pattern = r"^[A-Za-z]+[-_]?\d+$"
            features["has_id_pattern"] = valid_data.str.match(prefix_number_pattern).mean() > 0.8

            uuid_pattern = r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$"
            features["is_uuid"] = valid_data.str.match(uuid_pattern, flags=re.IGNORECASE).mean() > 0.8
        else:
            features.update(
                {
                    "avg_str_length": None,
                    "std_str_length": None,
                    "min_str_length": None,
                    "max_str_length": None,
                    "is_fixed_length": False,
                    "has_id_pattern": False,
                    "is_uuid": False,
                }
            )

        if features["n_unique"] <= 100:
            value_counts = col_data.value_counts(dropna=True)
            if len(value_counts) > 0:
                features["most_common_ratio"] = value_counts.iloc[0] / n_valid if n_valid > 0 else 0
                probs = value_counts / n_valid
                entropy = -np.sum(probs * np.log2(probs + 1e-10))
                max_entropy = np.log2(len(value_counts)) if len(value_counts) > 0 else 1
                features["category_balance"] = entropy / max_entropy if max_entropy > 0 else 0
            else:
                features["most_common_ratio"] = 0
                features["category_balance"] = 0
        else:
            features["most_common_ratio"] = None
            features["category_balance"] = None

        return features

    def _normalize_number(self, value):
        if value is None:
            return 0.0
        if isinstance(value, (float, np.floating)) and np.isnan(value):
            return 0.0
        return float(value)

    def compute_info_strength(self, features, role=None):
        if features.get("n_valid", 0) == 0:
            return 0.0

        role = role or "unknown"

        if role == "activity":
            cv = max(0.0, self._normalize_number(features.get("cv")))
            strength = min(1.0, cv / 0.5)
            if features.get("is_float") and features.get("has_decimal"):
                strength = min(1.0, strength + 0.1)
            return max(0.0, strength)

        if role in {"label", "experimental_context", "time"} or not features.get("is_numeric", False):
            entropy_norm = max(0.0, self._normalize_number(features.get("category_balance")))
            n_unique = max(int(features.get("n_unique", 0)), 0)
            cardinality_norm = 0.0
            if n_unique > 0:
                cardinality_norm = min(1.0, np.log2(1 + n_unique) / np.log2(1 + 10))
            return max(0.0, min(1.0, 0.6 * entropy_norm + 0.4 * cardinality_norm))

        cv = max(0.0, self._normalize_number(features.get("cv")))
        return max(0.0, min(1.0, cv / 0.75))

    def _role_to_legacy_category(self, role):
        if role in self.POSITIVE_ROLES:
            return "useful"
        if role in {"excluded", "derived_or_predicted"}:
            return "excluded"
        return "unknown"

    def classify_role_by_features(self, col_name, col_data, features=None):
        if features is None:
            features = self.extract_features(col_name, col_data)

        result = {"role": "unknown", "reason": "No clear pattern in statistical features", "confidence": 0.0}

        if features["missing_rate"] >= 1.0:
            result = {"role": "excluded", "reason": "All values missing", "confidence": 1.0}
        elif features["uniqueness_ratio"] > 0.95:
            if features["is_uuid"]:
                result = {"role": "excluded", "reason": "UUID identifier", "confidence": 0.95}
            elif features["has_id_pattern"]:
                result = {"role": "excluded", "reason": "ID pattern (PREFIX-NUMBER)", "confidence": 0.9}
            elif features["is_integer"] and features["is_monotonic_increasing"]:
                result = {"role": "excluded", "reason": "Sequential integer ID", "confidence": 0.85}
            elif features["is_integer"] and features["uniqueness_ratio"] > 0.98:
                result = {
                    "role": "excluded",
                    "reason": "High-cardinality integer (likely identifier)",
                    "confidence": 0.8,
                }
            elif not (features["is_float"] and features["has_decimal"]):
                if features["is_object"] and self._normalize_number(features.get("avg_str_length")) > 10:
                    result = {
                        "role": "excluded",
                        "reason": "High-cardinality text field (likely name or description)",
                        "confidence": 0.75,
                    }
        elif (
            features["is_object"]
            and 0.5 < features["uniqueness_ratio"] < 0.95
            and self._normalize_number(features.get("avg_str_length")) > 8
            and self._normalize_number(features.get("std_str_length")) > 3
        ):
            result = {"role": "excluded", "reason": "Free-text field (likely name or description)", "confidence": 0.7}
        elif features["n_unique"] == 1:
            result = {"role": "excluded", "reason": "Constant column", "confidence": 0.95}
        elif self._normalize_number(features.get("most_common_ratio")) > 0.99:
            result = {"role": "excluded", "reason": "Near-constant column", "confidence": 0.9}
        else:
            cv = self._normalize_number(features.get("cv"))

            if features["is_numeric"] and cv > 0.12:
                if features["is_float"] and features["has_decimal"]:
                    confidence = 0.85 if features["uniqueness_ratio"] > 0.7 else 0.75
                    result = {
                        "role": "activity",
                        "reason": f"Continuous measurement pattern (CV={cv:.2f})",
                        "confidence": confidence,
                    }
                elif features["is_integer"] and features["uniqueness_ratio"] > 0.25 and features["n_unique"] > 5:
                    result = {
                        "role": "activity",
                        "reason": f"Discrete measurement pattern (CV={cv:.2f})",
                        "confidence": 0.65,
                    }
            elif 2 <= features["n_unique"] <= 20:
                balance = self._normalize_number(features.get("category_balance"))
                if balance > 0.6 and features["uniqueness_ratio"] < 0.6:
                    result = {
                        "role": "label",
                        "reason": f"Balanced categorical pattern ({features['n_unique']} categories)",
                        "confidence": 0.65,
                    }
                elif balance > 0.3 and features["uniqueness_ratio"] < 0.8:
                    result = {
                        "role": "label",
                        "reason": f"Categorical pattern ({features['n_unique']} categories)",
                        "confidence": 0.55,
                    }
            elif 20 < features["n_unique"] <= 100 and features["uniqueness_ratio"] < 0.95:
                confidence = 0.55 if features["is_object"] else 0.45
                result = {
                    "role": "experimental_context",
                    "reason": f"Context-like categorical pattern ({features['n_unique']} categories)",
                    "confidence": confidence,
                }

        result["features"] = features
        result["info_strength"] = self.compute_info_strength(features, result["role"])
        return result

    def classify_by_features(self, col_name, col_data, features=None):
        result = self.classify_role_by_features(col_name, col_data, features)
        return self._role_to_legacy_category(result["role"]), result["reason"], result["confidence"]

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

            features = self.extract_features(col, df[col])
            role_result = self.classify_role_by_features(col, df[col], features)
            category = self._role_to_legacy_category(role_result["role"])
            entry = (col, role_result["reason"], role_result["confidence"], features)
            results[category].append(entry)

        return results
