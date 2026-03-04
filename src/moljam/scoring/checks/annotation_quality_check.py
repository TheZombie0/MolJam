from ..._logging import get_logger

logger = get_logger(__name__)


class AnnotationQualityChecksMixin:
    def check_annotation_quality(self):
        """
        Check annotation quality using HybridColumnClassifier and time label availability.
        Evaluates: time label availability, useful column quality, classification confidence, and type diversity.
        """
        if not self.include_experimental_info:
            print("Annotation quality not included in scoring, skipping")
            self.completed_checks.add('check_annotation_quality')
            return 0

        total_score = 0

        # 1. Time Label Availability (10 points)
        time_label_score = 0
        time_label_info = {}
        if self.time_col and self.time_col in self.df.columns:
            time_count = self.df[self.time_col].notna().sum()
            time_coverage = (time_count / self.num_molecules) * 100 if self.num_molecules > 0 else 0
            unique_times = self.df[self.time_col].dropna().unique()
            num_time_points = len(unique_times)

            missing_time_rate = 100 - time_coverage
            time_label_score = self.calculate_quality_score(missing_time_rate, max_score=10)

            time_label_info = {
                'Time label column': self.time_col,
                'Time coverage rate': f"{time_coverage:.2f}%",
                'Number of unique time points': num_time_points,
                'Example time points': list(unique_times[:10])
            }
        else:
            time_label_info = {
                'Time label column': self.time_col if self.time_col else 'None',
                'Time coverage rate': "0.00%",
                'Note': "No time label column provided or column does not exist"
            }

        self.scores["Experimental Information Quality"]["Time Label Availability"] = time_label_score
        total_score += time_label_score

        # 2-4. Column classification using HybridColumnClassifier
        try:
            import sys
            import os
            origin_dir = os.path.dirname(os.path.abspath(__file__))
            if origin_dir not in sys.path:
                sys.path.insert(0, origin_dir)
            from ...classification import HybridColumnClassifier

            classifier = HybridColumnClassifier()
            results = classifier.classify_columns(self.df, smiles_col=self.smiles_col)

            useful_cols = results['useful']
            excluded_cols = results['excluded']
            unknown_cols = results['unknown']

            total_cols = len(useful_cols) + len(excluded_cols) + len(unknown_cols)

            if total_cols == 0:
                print("No columns to analyze for annotation quality")
                self.scores["Experimental Information Quality"]["Useful Column Quality"] = 0
                self.scores["Experimental Information Quality"]["Classification Confidence"] = 0
                self.scores["Experimental Information Quality"]["Type Diversity"] = 0
                self.completed_checks.add('check_annotation_quality')
                return time_label_score

            # 2. Useful Column Quality (10 points) = average of Coverage and Ratio
            # Coverage: average non-null ratio of useful columns
            if len(useful_cols) > 0:
                coverage_rates = []
                for col_info in useful_cols:
                    col_name = col_info[0]
                    if col_name in self.df.columns:
                        coverage = self.df[col_name].notna().sum() / len(self.df)
                        coverage_rates.append(coverage)
                avg_coverage = sum(coverage_rates) / len(coverage_rates) if coverage_rates else 0
            else:
                avg_coverage = 0

            # Ratio: useful columns / total columns
            useful_ratio = len(useful_cols) / total_cols if total_cols > 0 else 0

            # Combine coverage and ratio using average
            useful_column_quality = ((avg_coverage + useful_ratio) / 2) * 10  # 0-10 points
            self.scores["Experimental Information Quality"]["Useful Column Quality"] = useful_column_quality
            total_score += useful_column_quality

            # 3. Classification Confidence (10 points)
            if len(useful_cols) > 0:
                confidences = [col_info[2] for col_info in useful_cols]
                avg_confidence = sum(confidences) / len(confidences)
                confidence_score = avg_confidence * 10
            else:
                avg_confidence = 0
                confidence_score = 0

            self.scores["Experimental Information Quality"]["Classification Confidence"] = confidence_score
            total_score += confidence_score

            # 4. Type Diversity (10 points)
            type_keywords = {
                'activity': ['Activity', 'Measurement', 'Predicted', 'Calculated'],
                'label': ['label', 'categories', 'Classification'],
                'experimental': ['Experimental', 'condition']
            }

            found_types = set()
            for col_info in useful_cols:
                reason = col_info[1]
                for type_name, keywords in type_keywords.items():
                    if any(kw.lower() in reason.lower() for kw in keywords):
                        found_types.add(type_name)

            num_types = len(found_types)
            diversity_score = (num_types / 3) * 10
            self.scores["Experimental Information Quality"]["Type Diversity"] = diversity_score
            total_score += diversity_score

            # Store analysis results
            self.analysis_results['Experimental Information Quality'] = {
                'Time Label': time_label_info,
                'Total columns analyzed': total_cols,
                'Useful columns': len(useful_cols),
                'Excluded columns': len(excluded_cols),
                'Unknown columns': len(unknown_cols),
                'Average coverage of useful columns': f"{avg_coverage*100:.2f}%",
                'Useful column ratio': f"{useful_ratio*100:.2f}%",
                'Useful column quality score': f"{useful_column_quality:.2f}/5",
                'Average classification confidence': f"{avg_confidence:.2f}",
                'Types found': list(found_types),
                'Useful column details': [
                    {'name': col[0], 'reason': col[1], 'confidence': col[2]}
                    for col in useful_cols[:10]
                ]
            }

            print(f"Annotation quality analysis: time_label={time_label_score:.2f}/10, "
                  f"useful_cols={len(useful_cols)}, quality={useful_column_quality:.2f}/10, "
                  f"confidence={confidence_score:.2f}/10, types={num_types}, "
                  f"total: {total_score:.2f}/40")

        except ImportError as e:
            print(f"Warning: Could not import HybridColumnClassifier: {e}")
            self.scores["Experimental Information Quality"]["Useful Column Quality"] = 0
            self.scores["Experimental Information Quality"]["Classification Confidence"] = 0
            self.scores["Experimental Information Quality"]["Type Diversity"] = 0
            self.analysis_results['Experimental Information Quality'] = {
                'Time Label': time_label_info,
                'Note': 'Column classification skipped due to import error'
            }
        except Exception as e:
            print(f"Error during annotation quality check: {e}")
            self.scores["Experimental Information Quality"]["Useful Column Quality"] = 0
            self.scores["Experimental Information Quality"]["Classification Confidence"] = 0
            self.scores["Experimental Information Quality"]["Type Diversity"] = 0

        self.completed_checks.add('check_annotation_quality')
        return total_score

