from .._common import *


class ReportMixin:
    def _normalize_report_value(self, value):
        """Convert numpy scalars/arrays into native Python values for readable reports."""
        if isinstance(value, np.generic):
            return value.item()

        if isinstance(value, np.ndarray):
            return [self._normalize_report_value(item) for item in value.tolist()]

        if isinstance(value, dict):
            return {
                key: self._normalize_report_value(subvalue)
                for key, subvalue in value.items()
            }

        if isinstance(value, (list, tuple, set)):
            return [self._normalize_report_value(item) for item in value]

        return value

    def get_detailed_report(self):
        """Generate detailed report"""
        # Ensure all checks are completed
        self.run_missing_checks()

        report = "# Molecular Database Quality Assessment Report\n\n"

        report += f"## Dataset Information\n"
        report += f"- Total molecules: {self.num_molecules}\n"
        if hasattr(self, 'valid_mols'):
            report += f"- Valid molecules: {len(self.valid_mols)}\n"
        report += "\n"

        structure_max = 3 * 10  # 3 metrics, 10 points each
        quality_max = 2 * 10  # 2 metrics, 10 points each
        coverage_max = 2 * 10  # 2 metrics, 10 points each
        distribution_max = 2 * 10  # 2 metrics, 10 points each

        if self.experimental_info:
            annotation_max = 30  # 3 metrics × 10 points each
            total_max = structure_max + quality_max + annotation_max + coverage_max + distribution_max

            # Calculate normalized scores for 5-module system
            normalized_structure = (self.scores['Structural Integrity']['Total'] / structure_max * 20) if structure_max > 0 else 0
            normalized_quality = (self.scores['Data Quality']['Total'] / quality_max * 20) if quality_max > 0 else 0
            normalized_annotation = (self.scores['Experimental Information Quality']['Total'] / annotation_max * 20) if annotation_max > 0 else 0
            normalized_coverage = (self.scores['Chemical Space Coverage']['Total'] / coverage_max * 20) if coverage_max > 0 else 0
            normalized_distribution = (self.scores['Data Distribution']['Total'] / distribution_max * 20) if distribution_max > 0 else 0
        else:
            annotation_max = 0
            total_max = structure_max + quality_max + coverage_max + distribution_max

            # Calculate normalized scores for 4-module system
            normalized_structure = (self.scores['Structural Integrity']['Total'] / structure_max * 25) if structure_max > 0 else 0
            normalized_quality = (self.scores['Data Quality']['Total'] / quality_max * 25) if quality_max > 0 else 0
            normalized_annotation = 0
            normalized_coverage = (self.scores['Chemical Space Coverage']['Total'] / coverage_max * 25) if coverage_max > 0 else 0
            normalized_distribution = (self.scores['Data Distribution']['Total'] / distribution_max * 25) if distribution_max > 0 else 0

        report += "## Score Summary\n"
        if self.experimental_info:
            report += f"- Structural Integrity: {self.scores['Structural Integrity']['Total']:.2f}/{structure_max} (Normalized: {normalized_structure:.2f}/20)\n"
            report += f"- Data Quality: {self.scores['Data Quality']['Total']:.2f}/{quality_max} (Normalized: {normalized_quality:.2f}/20)\n"
            report += f"- Experimental Information Quality: {self.scores['Experimental Information Quality']['Total']:.2f}/{annotation_max} (Normalized: {normalized_annotation:.2f}/20)\n"
            report += f"- Chemical Space Coverage: {self.scores['Chemical Space Coverage']['Total']:.2f}/{coverage_max} (Normalized: {normalized_coverage:.2f}/20)\n"
            report += f"- Data Distribution: {self.scores['Data Distribution']['Total']:.2f}/{distribution_max} (Normalized: {normalized_distribution:.2f}/20)\n"
        else:
            report += f"- Structural Integrity: {self.scores['Structural Integrity']['Total']:.2f}/{structure_max} (Normalized: {normalized_structure:.2f}/25)\n"
            report += f"- Data Quality: {self.scores['Data Quality']['Total']:.2f}/{quality_max} (Normalized: {normalized_quality:.2f}/25)\n"
            report += f"- Chemical Space Coverage: {self.scores['Chemical Space Coverage']['Total']:.2f}/{coverage_max} (Normalized: {normalized_coverage:.2f}/25)\n"
            report += f"- Data Distribution: {self.scores['Data Distribution']['Total']:.2f}/{distribution_max} (Normalized: {normalized_distribution:.2f}/25)\n"
            report += f"  *Note: Experimental Information Quality not included in scoring (experimental_info=False)*\n"
        report += f"- **Total Score**: {self.scores['Total Score']:.2f}/{total_max}\n"
        report += f"- **Normalized Score**: {self.scores['Normalized Score']:.2f}/100\n"
        report += f"- **Final Adjusted Score**: {self.scores['Final Adjusted Score']:.2f}/100\n\n"

        # Add detailed analysis results
        report += "## Detailed Analysis Results\n\n"

        for category, results in self.analysis_results.items():
            results = self._normalize_report_value(results)
            report += f"### {category}\n"
            for key, value in results.items():
                if isinstance(value, dict):
                    report += f"- {key}:\n"
                    for subkey, subvalue in value.items():
                        if isinstance(subvalue, dict):
                            report += f"  - {subkey}:\n"
                            for subsubkey, subsubvalue in subvalue.items():
                                report += f"    - {subsubkey}: {subsubvalue}\n"
                        else:
                            report += f"  - {subkey}: {subvalue}\n"
                else:
                    report += f"- {key}: {value}\n"
            report += "\n"

        return report
