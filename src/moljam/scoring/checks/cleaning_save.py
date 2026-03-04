from ..._logging import get_logger

logger = get_logger(__name__)


class CleaningSaveMixin:
    def save_cleaned_database(self, output_path, remove_invalid_smiles=True,
                             remove_undefined_stereochemistry=True,
                             remove_conflicting_labels=True,
                             remove_consistent_duplicates=True,
                             verbose=True):
        """
        Clean the database and save to file.
        
        Parameters:
            output_path: str - Path to save the cleaned CSV file
            Other parameters same as clean_database()
            
        Returns:
            cleaned_df: pandas DataFrame - The cleaned dataset
            cleaning_report: dict - Details about cleaning
        """
        # Clean the database
        cleaned_df, cleaning_report = self.clean_database(
            remove_invalid_smiles=remove_invalid_smiles,
            remove_undefined_stereochemistry=remove_undefined_stereochemistry,
            remove_conflicting_labels=remove_conflicting_labels,
            remove_consistent_duplicates=remove_consistent_duplicates,
            verbose=verbose
        )
        
        # Save to file
        cleaned_df.to_csv(output_path, index=False)
        print(f"\nCleaned database saved to: {output_path}")
        
        # Save cleaning report
        report_path = output_path.replace('.csv', '_cleaning_report.txt')
        with open(report_path, 'w') as f:
            f.write("=== Molecular Database Cleaning Report ===\n\n")
            for category, details in cleaning_report.items():
                f.write(f"{category}:\n")
                for key, value in details.items():
                    f.write(f"  {key}: {value}\n")
                f.write("\n")
        print(f"Cleaning report saved to: {report_path}")
        
        return cleaned_df, cleaning_report

