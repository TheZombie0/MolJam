"""
Hybrid classifier combining keyword and statistical signals.

The hybrid output keeps legacy compatibility while also returning a structured
annotation profile used by the scoring pipeline.
"""

from .keyword import KeywordColumnClassifier
from .statistical import StatisticalColumnClassifier


class HybridColumnClassifier:
    """Merge keyword and statistical votes into a structured column profile."""

    ACCEPTED_ROLES = {"time", "activity", "label", "experimental_context"}

    def __init__(self, keyword_weight=0.6, statistical_weight=0.4):
        self.keyword_classifier = KeywordColumnClassifier()
        self.statistical_classifier = StatisticalColumnClassifier()

        self.keyword_weight = keyword_weight
        self.statistical_weight = statistical_weight

    def _clamp(self, value, lower=0.0, upper=1.0):
        return max(lower, min(upper, float(value)))

    def _statistical_activity_allowed(self, col_name):
        return self.keyword_classifier.looks_like_measurement_name(col_name)

    def _role_to_legacy_category(self, role, accepted):
        if accepted:
            return "useful"
        if role in {"excluded", "derived_or_predicted"}:
            return "excluded"
        return "unknown"

    def _merge_explicit_role(self, explicit_role, keyword_vote, statistical_vote):
        support_count = 0
        conflict_count = 0

        for vote in (keyword_vote, statistical_vote):
            if vote["role"] == explicit_role:
                support_count += 1
            elif vote["role"] not in {explicit_role, "unknown"}:
                conflict_count += 1

        agreement = support_count >= 1 or conflict_count == 0
        conflicted = conflict_count > 0
        confidence = 0.9 + 0.03 * support_count - 0.05 * conflict_count
        confidence = self._clamp(confidence)

        reason_parts = [f"Explicit prior role '{explicit_role}'"]
        if support_count > 0:
            reason_parts.append(f"supported by {support_count} classifier(s)")
        if conflicted:
            reason_parts.append("overrides conflicting learned signals")

        notes = ["explicit_prior"]
        if support_count > 0:
            notes.append("classifier_support")
        if conflicted:
            notes.append("explicit_conflict")

        return {
            "role": explicit_role,
            "confidence": confidence,
            "agreement": agreement,
            "conflicted": conflicted,
            "decision_type": "explicit_prior",
            "reason": " | ".join(reason_parts),
            "notes": notes,
        }

    def _merge_votes(self, col_name, keyword_vote, statistical_vote):
        keyword_role = keyword_vote["role"]
        keyword_conf = keyword_vote["confidence"]
        keyword_reason = keyword_vote["reason"]

        statistical_role = statistical_vote["role"]
        statistical_conf = statistical_vote["confidence"]
        statistical_reason = statistical_vote["reason"]

        if keyword_role == statistical_role and keyword_role != "unknown":
            final_conf = self._clamp(
                keyword_conf * self.keyword_weight + statistical_conf * self.statistical_weight + 0.05
            )
            return {
                "role": keyword_role,
                "confidence": final_conf,
                "agreement": True,
                "conflicted": False,
                "decision_type": "agreement",
                "reason": f"{keyword_reason} | Confirmed by statistics: {statistical_reason}",
                "notes": ["agreement"],
            }

        if keyword_conf >= 0.85 and keyword_role != "unknown":
            conflicted = statistical_role not in {keyword_role, "unknown"}
            final_conf = keyword_conf if not conflicted else keyword_conf * 0.9
            note = "keyword_override_conflict" if conflicted else "keyword_override"
            reason = keyword_reason
            if conflicted:
                reason = f"{keyword_reason} | Overrides statistical signal: {statistical_reason}"
            elif statistical_role == "unknown":
                reason = f"{keyword_reason} | No statistical support needed"
            return {
                "role": keyword_role,
                "confidence": self._clamp(final_conf),
                "agreement": False,
                "conflicted": conflicted,
                "decision_type": "keyword_override",
                "reason": reason,
                "notes": [note],
            }

        if statistical_conf >= 0.8 and keyword_role == "unknown" and statistical_role != "unknown":
            if statistical_role == "activity" and not self._statistical_activity_allowed(col_name):
                return {
                    "role": "unknown",
                    "confidence": 0.0,
                    "agreement": False,
                    "conflicted": False,
                    "decision_type": "statistical_rejected",
                    "reason": f"{statistical_reason} | Rejected because column name does not look like an activity measurement",
                    "notes": ["statistical_rejected"],
                }
            return {
                "role": statistical_role,
                "confidence": self._clamp(statistical_conf * 0.8),
                "agreement": False,
                "conflicted": False,
                "decision_type": "statistical_only",
                "reason": f"{statistical_reason} | No keyword signal",
                "notes": ["statistical_only"],
            }

        if keyword_role != "unknown" and statistical_role == "unknown":
            return {
                "role": keyword_role,
                "confidence": self._clamp(keyword_conf * 0.85),
                "agreement": False,
                "conflicted": False,
                "decision_type": "keyword_only",
                "reason": f"{keyword_reason} | Weak statistical signal",
                "notes": ["keyword_only"],
            }

        if keyword_role == "unknown" and statistical_role != "unknown":
            if statistical_role == "activity" and not self._statistical_activity_allowed(col_name):
                return {
                    "role": "unknown",
                    "confidence": 0.0,
                    "agreement": False,
                    "conflicted": False,
                    "decision_type": "statistical_rejected",
                    "reason": f"{statistical_reason} | Rejected because column name does not look like an activity measurement",
                    "notes": ["statistical_rejected"],
                }
            return {
                "role": statistical_role,
                "confidence": self._clamp(statistical_conf * 0.75),
                "agreement": False,
                "conflicted": False,
                "decision_type": "statistical_support",
                "reason": f"{statistical_reason} | No keyword match",
                "notes": ["statistical_support"],
            }

        if keyword_role != "unknown" and statistical_role != "unknown":
            keyword_score = keyword_conf * self.keyword_weight
            statistical_score = statistical_conf * self.statistical_weight
            if keyword_score >= statistical_score:
                final_role = keyword_role
                final_conf = keyword_score * 0.75
                final_reason = f"{keyword_reason} | Conflicts with statistical signal: {statistical_reason}"
            else:
                final_role = statistical_role
                final_conf = statistical_score * 0.75
                final_reason = f"{statistical_reason} | Conflicts with keyword signal: {keyword_reason}"

            return {
                "role": final_role,
                "confidence": self._clamp(final_conf),
                "agreement": False,
                "conflicted": True,
                "decision_type": "conflict",
                "reason": final_reason,
                "notes": ["conflict"],
            }

        return {
            "role": "unknown",
            "confidence": 0.0,
            "agreement": True,
            "conflicted": False,
            "decision_type": "unknown",
            "reason": "No clear signal from either method",
            "notes": ["unknown"],
        }

    def classify_column(self, col_name, col_data, explicit_role=None):
        keyword_vote = self.keyword_classifier.classify_role(col_name, col_data)
        statistical_features = self.statistical_classifier.extract_features(col_name, col_data)
        statistical_vote = self.statistical_classifier.classify_role_by_features(
            col_name,
            col_data,
            statistical_features,
        )

        if explicit_role is not None:
            merged = self._merge_explicit_role(explicit_role, keyword_vote, statistical_vote)
        else:
            merged = self._merge_votes(col_name, keyword_vote, statistical_vote)

        final_role = merged["role"]
        confidence = merged["confidence"]
        accepted = final_role in self.ACCEPTED_ROLES and (
            explicit_role is not None or confidence >= 0.55
        )
        coverage = statistical_features["n_valid"] / statistical_features["n_total"] if statistical_features["n_total"] > 0 else 0.0
        info_strength = self.statistical_classifier.compute_info_strength(statistical_features, final_role)
        legacy_category = self._role_to_legacy_category(final_role, accepted)

        return {
            "column": col_name,
            "category": legacy_category,
            "role": final_role,
            "accepted": accepted,
            "reason": merged["reason"],
            "confidence": confidence,
            "agreement": merged["agreement"],
            "conflicted": merged["conflicted"],
            "decision_type": merged["decision_type"],
            "coverage": coverage,
            "info_strength": info_strength,
            "keyword_vote": {
                "role": keyword_vote["role"],
                "reason": keyword_vote["reason"],
                "confidence": keyword_vote["confidence"],
            },
            "statistical_vote": {
                "role": statistical_vote["role"],
                "reason": statistical_vote["reason"],
                "confidence": statistical_vote["confidence"],
            },
            "details": keyword_vote["details"],
            "features": statistical_features,
            "notes": merged["notes"],
        }

    def build_annotation_profile(self, df, smiles_col="Smiles", explicit_roles=None):
        explicit_roles = explicit_roles or {}

        profile = {
            "columns": [],
            "accepted_columns": [],
            "by_role": {
                "time": [],
                "activity": [],
                "label": [],
                "experimental_context": [],
                "derived_or_predicted": [],
                "excluded": [],
                "unknown": [],
            },
            "summary": {},
        }

        smiles_col_actual = None
        for col in df.columns:
            if col.lower() == smiles_col.lower():
                smiles_col_actual = col
                break

        for col in df.columns:
            if smiles_col_actual and col == smiles_col_actual:
                continue

            result = self.classify_column(col, df[col], explicit_role=explicit_roles.get(col))
            profile["columns"].append(result)
            profile["by_role"].setdefault(result["role"], []).append(result)
            if result["accepted"]:
                profile["accepted_columns"].append(result)

        accepted_columns = profile["accepted_columns"]
        role_counts = {
            role: len(columns)
            for role, columns in profile["by_role"].items()
        }

        agreement_rate = (
            sum(1 for column in accepted_columns if column["agreement"]) / len(accepted_columns)
            if accepted_columns
            else 0.0
        )
        avg_confidence = (
            sum(column["confidence"] for column in accepted_columns) / len(accepted_columns)
            if accepted_columns
            else 0.0
        )
        avg_coverage = (
            sum(column["coverage"] for column in accepted_columns) / len(accepted_columns)
            if accepted_columns
            else 0.0
        )
        avg_info_strength = (
            sum(column["info_strength"] for column in accepted_columns) / len(accepted_columns)
            if accepted_columns
            else 0.0
        )

        profile["summary"] = {
            "total_columns": len(profile["columns"]),
            "accepted_count": len(accepted_columns),
            "agreement_rate": agreement_rate,
            "avg_confidence": avg_confidence,
            "avg_coverage": avg_coverage,
            "avg_info_strength": avg_info_strength,
            "role_counts": role_counts,
        }

        return profile

    def classify_columns(self, df, smiles_col="Smiles", explicit_roles=None):
        profile = self.build_annotation_profile(df, smiles_col=smiles_col, explicit_roles=explicit_roles)

        results = {"useful": [], "excluded": [], "unknown": []}
        for column in profile["columns"]:
            entry = (
                column["column"],
                column["reason"],
                column["confidence"],
                {
                    "role": column["role"],
                    "coverage": column["coverage"],
                    "info_strength": column["info_strength"],
                    "details": column["details"],
                },
                column["agreement"],
            )
            results[column["category"]].append(entry)

        results["profile"] = profile
        return results
