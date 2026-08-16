from .._common import *
from ..protonation import derive_parent_form_references
import re


class RepresentationConsistencyMixin:
    _ACID_KEYWORDS = (
        "acetate",
        "acetic acid",
        "formate",
        "formic acid",
        "methanesulfonic acid",
        "mesylate",
        "trifluoroacetate",
        "trifluoroacetic acid",
        "tosylate",
        "p-toluenesulfonic acid",
        "besylate",
        "benzenesulfonic acid",
        "succinate",
        "succinic acid",
        "oxalate",
        "oxalic acid",
        "lactate",
        "lactic acid",
        "maleate",
        "maleic",
        "fumarate",
        "fumaric",
        "tartrate",
        "tartaric acid",
        "citrate",
        "citric acid",
    )

    @staticmethod
    def _representation_signature(row):
        return (
            row['canonical_smiles'],
            row['observed_parent_smiles'],
            tuple(row['removed_salts'] or []),
            tuple(row['removed_solvents'] or []),
            tuple(row['duplicate_parent_fragments'] or []),
        )

    @staticmethod
    def _extract_english_label(raw_name):
        match = re.search(r"\(([^()]*)\)", raw_name)
        label = match.group(1) if match else raw_name
        return label.replace("_", " ")

    @classmethod
    def _formal_charge(cls, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0
        return Chem.GetFormalCharge(mol)

    @staticmethod
    def _canonicalize_smiles_for_match(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

    @classmethod
    def _representation_tags(
        cls,
        canonical_smiles,
        observed_parent_smiles,
        parent_form_smiles,
        parent_form_candidates,
        removed_salts,
        removed_solvents,
        duplicate_parent_fragments,
    ):
        tags = []

        removed_salts_display = [cls._extract_english_label(item) for item in removed_salts]
        removed_solvents_display = [cls._extract_english_label(item) for item in removed_solvents]
        duplicate_parent_fragments_display = list(duplicate_parent_fragments or [])

        acid_hits = []
        salt_hits = []
        for item in removed_salts_display:
            lowered = item.lower()
            if any(keyword in lowered for keyword in cls._ACID_KEYWORDS):
                acid_hits.append(item)
            else:
                salt_hits.append(item)

        if salt_hits:
            tags.append('salt')
        if acid_hits:
            tags.append('acid adduct')
        if removed_solvents_display:
            tags.append('solvent stripping')

        if duplicate_parent_fragments_display:
            tags.append('duplicate-component')

        parent_form_candidates = {
            cls._canonicalize_smiles_for_match(candidate)
            for candidate in (parent_form_candidates or [parent_form_smiles])
        }
        canonical_smiles_for_match = cls._canonicalize_smiles_for_match(canonical_smiles)
        parent_form_matches = (
            canonical_smiles_for_match in parent_form_candidates
            and not removed_salts_display
            and not removed_solvents_display
            and not duplicate_parent_fragments_display
        )

        if not parent_form_matches:
            observed_parent_charge = cls._formal_charge(observed_parent_smiles)
            parent_form_charge = cls._formal_charge(parent_form_smiles)
            if observed_parent_charge > parent_form_charge:
                tags.append('protonated')
            elif observed_parent_charge < parent_form_charge:
                tags.append('deprotonated')

        if parent_form_matches:
            tags = ['parent form']
        elif not tags:
            tags.append('other non-parent form')

        return (
            tags,
            removed_salts_display,
            removed_solvents_display,
            duplicate_parent_fragments_display,
            parent_form_matches,
        )

    def check_representation_consistency(self):
        """
        Check for the same parent molecule appearing in different
        acid/base/salt/solvate representations.
        """
        if 'validate_smiles' not in self.completed_checks:
            self.validate_smiles()

        if len(self.valid_mols) == 0 or self.valid_df.empty:
            print("No valid molecules, skipping representation consistency check")
            score = 0
            self.representation_consistency_groups = []
            self.scores["Structural Integrity"]["Representation Consistency"] = score
            self.analysis_results['Representation Consistency'] = {
                'Molecules with multiple forms': 0,
                'Total redundant molecules': 0,
                'Redundancy rate': "0.00%",
                'Example groups': [],
                'Note': "No valid molecules"
            }
            self.completed_checks.add('check_representation_consistency')
            return score

        print("Starting representation consistency check...")
        start_time = time.time()

        grouped = self.valid_df.groupby('parent_smiles', sort=False)
        redundant_groups = []
        total_redundant_molecules = 0
        grouped_items = [
            (parent_smiles, group)
            for parent_smiles, group in grouped
            if len(group) > 1
        ]
        parent_form_references = derive_parent_form_references(
            [parent_smiles for parent_smiles, _ in grouped_items],
            backend=getattr(self, 'parent_form_backend', 'dimorphite_dl'),
            ph=getattr(self, 'parent_form_ph', 7.4),
            chemaxon_executable=getattr(self, 'chemaxon_executable', 'cxcalc'),
            dimorphite_python=getattr(self, 'dimorphite_python', None),
            dimorphite_conda_env=getattr(self, 'dimorphite_conda_env', 'dimorphite'),
        )

        for parent_smiles, group in grouped_items:
            signature_groups = defaultdict(list)
            for _, row in group.iterrows():
                signature = self._representation_signature(row)
                signature_groups[signature].append(
                        {
                            'smiles': row[self.smiles_col],
                            'canonical_smiles': row['canonical_smiles'],
                            'standardized_smiles': row['standardized_smiles'],
                            'observed_parent_smiles': row['observed_parent_smiles'],
                            'original_index': int(row['original_index']) + 1,
                            'removed_salts': list(row['removed_salts'] or []),
                            'removed_solvents': list(row['removed_solvents'] or []),
                            'duplicate_parent_fragments': list(row['duplicate_parent_fragments'] or []),
                        }
                    )

            if len(signature_groups) <= 1:
                continue

            parent_form_reference = parent_form_references[parent_smiles]

            variant_groups = []
            group_tags = set()
            for signature, examples in signature_groups.items():
                canonical_smiles, observed_parent_smiles, removed_salts, removed_solvents, duplicate_parent_fragments = signature
                (
                    representation_tags,
                    removed_salts_display,
                    removed_solvents_display,
                    duplicate_parent_fragments_display,
                    parent_form_matches,
                ) = self._representation_tags(
                    canonical_smiles,
                    observed_parent_smiles,
                    parent_form_reference.smiles,
                    parent_form_reference.candidates,
                    removed_salts,
                    removed_solvents,
                    duplicate_parent_fragments,
                )
                group_tags.update(representation_tags)
                variant_groups.append(
                    {
                        'canonical_smiles': canonical_smiles,
                        'observed_parent_smiles': observed_parent_smiles,
                        'removed_salts': list(removed_salts),
                        'removed_solvents': list(removed_solvents),
                        'duplicate_parent_fragments': list(duplicate_parent_fragments),
                        'removed_salts_display': removed_salts_display,
                        'removed_solvents_display': removed_solvents_display,
                        'duplicate_parent_fragments_display': duplicate_parent_fragments_display,
                        'representation_tags': representation_tags,
                        'parent_form_matches': parent_form_matches,
                        'examples': examples[:5],
                        'count': len(examples),
                    }
                )

            redundant_groups.append(
                {
                    'group_parent_smiles': parent_smiles,
                    'parent_form': parent_form_reference.smiles,
                    'parent_form_candidates': parent_form_reference.candidates,
                    'parent_form_backend_requested': parent_form_reference.backend_requested,
                    'parent_form_backend_used': parent_form_reference.backend_used,
                    'parent_form_comment': parent_form_reference.comment,
                    'parent_form_ph': getattr(self, 'parent_form_ph', 7.4),
                    'molecule_count': len(group),
                    'distinct_representations': len(signature_groups),
                    'group_tags': sorted(group_tags),
                    'variants': variant_groups,
                }
            )
            total_redundant_molecules += len(group) - 1

        # Use score_count_based_issues instead of score_low_count_issues
        redundant_count = len(redundant_groups)
        score = self.score_count_based_issues(
            redundant_count,
            len(self.valid_mols),
            max_score=10,
            severity='medium'  # Medium severity for representation consistency
        )

        # Calculate redundancy rate
        redundancy_rate = (total_redundant_molecules / len(self.valid_mols)) * 100 if self.valid_mols else 0

        # Prepare example groups
        example_groups = []
        for group in redundant_groups[:5]:
            example_groups.append({
                'parent_form': group['parent_form'],
                'distinct_representations': group['distinct_representations'],
                'group_tags': group['group_tags'],
                'variants': group['variants'],
                'count': group['molecule_count']
            })

        self.representation_consistency_groups = redundant_groups
        self.analysis_results['Representation Consistency'] = {
            'Molecules with multiple forms': len(redundant_groups),
            'Total redundant molecules': total_redundant_molecules,
            'Redundancy rate': f"{redundancy_rate:.2f}%",
            'Example groups': example_groups
        }

        self.scores["Structural Integrity"]["Representation Consistency"] = score

        elapsed_time = time.time() - start_time
        print(f"Representation consistency check completed in {elapsed_time:.2f} seconds")
        print(f"Representation consistency: {len(redundant_groups)} molecules with multiple forms")
        print(f"Redundancy rate: {redundancy_rate:.2f}%, score: {score:.2f}/10")

        self.completed_checks.add('check_representation_consistency')
        return score
