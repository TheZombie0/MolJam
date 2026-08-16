"""
RDKit-centric helpers used by the scoring pipeline.

These functions are kept at module scope so they can be used with
`multiprocessing.Pool`.
"""

from __future__ import annotations

from collections import Counter

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, QED
from rdkit.Chem import MolStandardize
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmilesFromSmiles

try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError:  # pragma: no cover - fallback for older RDKit layouts
    rdMolStandardize = None


_COMMON_SOLVENT_SMARTS = (
    ("水(water)", "O"),
    ("甲醇(methanol)", "CO"),
    ("乙醇(ethanol)", "CCO"),
    ("异丙醇(isopropanol)", "CC(C)O"),
    ("丙酮(acetone)", "CC(=O)C"),
    ("乙腈(acetonitrile)", "CC#N"),
    ("二甲基亚砜(DMSO)", "CS(=O)C"),
    ("N,N-二甲基甲酰胺(DMF)", "CN(C)C=O"),
    ("四氢呋喃(THF)", "C1CCOC1"),
    ("1,4-二氧六环(dioxane)", "O1CCOCC1"),
    ("乙酸乙酯(ethyl acetate)", "CCOC(=O)C"),
    ("正丁醇(n-butanol)", "CCCCO"),
    ("叔丁醇(tert-butanol)", "CC(C)(C)O"),
    ("正己烷(n-hexane)", "CCCCCC"),
    ("正庚烷(n-heptane)", "CCCCCCC"),
    ("环己烷(cyclohexane)", "C1CCCCC1"),
    ("苯(benzene)", "c1ccccc1"),
    ("甲苯(toluene)", "Cc1ccccc1"),
    ("二氯甲烷(dichloromethane)", "ClCCl"),
    ("氯仿(chloroform)", "ClC(Cl)Cl"),
    ("二乙醚(diethyl ether)", "CCOCC"),
    ("甲基叔丁基醚(MTBE)", "COC(C)(C)C"),
    ("吡啶(pyridine)", "c1ccncc1"),
    ("N-甲基吡咯烷酮(NMP)", "CN1CCCC1=O"),
    ("N,N-二甲基乙酰胺(DMAc)", "CC(=O)N(C)C"),
)

_COMMON_SALT_SMARTS = (
    ("氯化物(chloride)", "[Cl-]"),
    ("中性氯化物(chloride_neutral)", "Cl"),
    ("溴化物(bromide)", "[Br-]"),
    ("中性溴化物(bromide_neutral)", "Br"),
    ("碘化物(iodide)", "[I-]"),
    ("中性碘化物(iodide_neutral)", "I"),
    ("氟化物(fluoride)", "[F-]"),
    ("中性氟化物(fluoride_neutral)", "F"),
    ("钠盐(sodium)", "[Na+]"),
    ("中性钠(sodium_neutral)", "[Na]"),
    ("钾盐(potassium)", "[K+]"),
    ("中性钾(potassium_neutral)", "[K]"),
    ("锂盐(lithium)", "[Li+]"),
    ("中性锂(lithium_neutral)", "[Li]"),
    ("钙盐(calcium)", "[Ca+2]"),
    ("中性钙(calcium_neutral)", "[Ca]"),
    ("镁盐(magnesium)", "[Mg+2]"),
    ("中性镁(magnesium_neutral)", "[Mg]"),
    ("乙酸根(acetate)", "CC(=O)[O-]"),
    ("乙酸(acetic acid)", "CC(=O)O"),
    ("甲酸根(formate)", "C(=O)[O-]"),
    ("甲酸(formic acid)", "C(=O)O"),
    ("甲磺酸根(mesylate)", "CS(=O)(=O)[O-]"),
    ("甲磺酸(methanesulfonic acid)", "CS(=O)(=O)O"),
    ("三氟乙酸根(trifluoroacetate)", "O=C([O-])C(F)(F)F"),
    ("三氟乙酸(trifluoroacetic acid)", "O=C(O)C(F)(F)F"),
    ("硫酸根(sulfate)", "[O-]S(=O)(=O)[O-]"),
    ("磷酸根(phosphate)", "O=P([O-])([O-])[O-]"),
    ("对甲苯磺酸根(tosylate)", "Cc1ccc(S(=O)(=O)[O-])cc1"),
    ("对甲苯磺酸(p-toluenesulfonic acid)", "Cc1ccc(S(=O)(=O)O)cc1"),
    ("苯磺酸根(besylate)", "O=S(=O)([O-])c1ccccc1"),
    ("苯磺酸(benzenesulfonic acid)", "O=S(=O)(O)c1ccccc1"),
    ("琥珀酸根(succinate)", "O=C([O-])CCC(=O)[O-]"),
    ("琥珀酸(succinic acid)", "O=C(O)CCC(=O)O"),
    ("草酸根(oxalate)", "O=C([O-])C(=O)[O-]"),
    ("草酸(oxalic acid)", "O=C(O)C(=O)O"),
    ("乳酸根(lactate)", "CC(O)C(=O)[O-]"),
    ("乳酸(lactic acid)", "CC(O)C(=O)O"),
    ("马来酸/富马酸根(maleate_or_fumarate)", "O=C([O-])C=CC(=O)[O-]"),
    ("马来酸/富马酸(maleic_or_fumaric_acid)", "O=C(O)C=CC(=O)O"),
    ("酒石酸根(tartrate)", "O=C([O-])C(O)C(O)C(=O)[O-]"),
    ("酒石酸(tartaric acid)", "OC(=O)C(O)C(O)C(=O)O"),
    ("柠檬酸根(citrate)", "O=C([O-])CC(O)(CC(=O)[O-])C(=O)[O-]"),
    ("柠檬酸(citric acid)", "OC(=O)CC(O)(CC(=O)O)C(=O)O"),
)

_SOLVENT_PATTERNS = [(name, Chem.MolFromSmarts(smarts)) for name, smarts in _COMMON_SOLVENT_SMARTS]
_SALT_PATTERNS = [(name, Chem.MolFromSmarts(smarts)) for name, smarts in _COMMON_SALT_SMARTS]

if rdMolStandardize is not None:
    _NORMALIZER = rdMolStandardize.Normalizer()
    try:
        _UNCHARGER = rdMolStandardize.Uncharger(canonicalOrder=True)
    except TypeError:  # pragma: no cover - older RDKit API
        _UNCHARGER = rdMolStandardize.Uncharger()
else:  # pragma: no cover - fallback path
    _NORMALIZER = MolStandardize.normalize.Normalizer()
    _UNCHARGER = MolStandardize.charge.Uncharger()


def _canonical_smiles_from_mol(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _clear_substance_groups(mol: Chem.Mol) -> None:
    if hasattr(Chem, "ClearMolSubstanceGroups"):
        try:
            Chem.ClearMolSubstanceGroups(mol)
        except Exception:
            pass


def _kekulize_if_possible(mol: Chem.Mol) -> None:
    try:
        Chem.Kekulize(mol, clearAromaticFlags=False)
    except Exception:
        pass


def _standardize_molecule(mol: Chem.Mol) -> tuple[Chem.Mol | None, str]:
    """
    Build a conservative standardized molecule.

    This keeps stereochemistry and isotopes intact, and intentionally does not
    canonicalize tautomers at this stage.
    """
    try:
        standardized = Chem.Mol(mol)
        standardized.UpdatePropertyCache(strict=False)
        _clear_substance_groups(standardized)
        _kekulize_if_possible(standardized)
        standardized = _NORMALIZER.normalize(standardized)
        standardized.UpdatePropertyCache(strict=False)
        return standardized, "Standardized"
    except Exception as exc:
        return None, f"Standardization failed: {exc}"


def _neutralize_for_grouping(mol: Chem.Mol) -> tuple[Chem.Mol | None, str]:
    """Build a charge-neutralized copy used only for grouping related forms."""
    try:
        neutralized = Chem.Mol(mol)
        neutralized.UpdatePropertyCache(strict=False)
        neutralized = _UNCHARGER.uncharge(neutralized)
        neutralized.UpdatePropertyCache(strict=False)
        return neutralized, "Neutralized for grouping"
    except Exception as exc:
        return None, f"Grouping neutralization failed: {exc}"


def _prepare_fragment_for_matching(frag: Chem.Mol) -> Chem.Mol:
    prepared = Chem.RemoveHs(Chem.Mol(frag), sanitize=False)
    prepared.UpdatePropertyCache(strict=False)
    try:
        Chem.SetAromaticity(prepared)
    except Exception:
        pass
    return prepared


def _match_fragment(
    fragment: Chem.Mol,
    prepared_fragment: Chem.Mol,
    patterns: list[tuple[str, Chem.Mol]],
) -> str | None:
    for name, pattern in patterns:
        if pattern is None:
            continue
        if (
            prepared_fragment.GetNumAtoms() == pattern.GetNumAtoms()
            and prepared_fragment.GetNumBonds() == pattern.GetNumBonds()
            and prepared_fragment.HasSubstructMatch(pattern)
        ):
            return name
    return None


def _rebuild_parent_from_fragments(fragments: list[Chem.Mol]) -> Chem.Mol:
    if len(fragments) == 1:
        return Chem.Mol(fragments[0])

    rebuilt = Chem.Mol(fragments[0])
    for fragment in fragments[1:]:
        rebuilt = Chem.CombineMols(rebuilt, fragment)
    return rebuilt


def _extract_parent_molecule(
    mol: Chem.Mol,
) -> tuple[Chem.Mol, list[str], list[str], list[str], bool, str]:
    """
    Extract a parent structure by stripping solvent/salt fragments conservatively.

    The logic follows the S6-style fragment workflow more closely than the
    previous MolJam canonical-only pipeline:
    - split into graph fragments
    - remove known solvents first
    - remove known salts second
    - keep multiple meaningful fragments instead of taking the longest SMILES
    - fall back if stripping removes everything
    """
    fragments = list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False))
    if len(fragments) <= 1:
        return Chem.Mol(mol), [], [], [], False, "Single-fragment parent"

    fragment_info = []
    for fragment in fragments:
        prepared = _prepare_fragment_for_matching(fragment)
        solvent_name = _match_fragment(fragment, prepared, _SOLVENT_PATTERNS)
        fragment_info.append(
            {
                "fragment": Chem.Mol(fragment),
                "prepared": prepared,
                "solvent_name": solvent_name,
            }
        )

    kept_after_solvent = []
    removed_solvents = []
    for info in fragment_info:
        if info["solvent_name"] is not None:
            removed_solvents.append(info["solvent_name"])
        else:
            kept_after_solvent.append(info)

    if not kept_after_solvent:
        return (
            Chem.Mol(mol),
            [],
            removed_solvents,
            [],
            True,
            "Parent fallback: solvent stripping removed all fragments",
        )

    kept_after_salt = []
    removed_salts = []
    for info in kept_after_solvent:
        salt_name = _match_fragment(info["fragment"], info["prepared"], _SALT_PATTERNS)
        if salt_name is not None:
            removed_salts.append(salt_name)
        else:
            kept_after_salt.append(info)

    if not kept_after_salt:
        rebuilt = _rebuild_parent_from_fragments([info["fragment"] for info in kept_after_solvent])
        fragment_smiles = [_canonical_smiles_from_mol(info["fragment"]) for info in kept_after_solvent]
        duplicate_fragments = sorted(
            smiles for smiles, count in Counter(fragment_smiles).items() if count > 1
        )
        return (
            rebuilt,
            removed_salts,
            removed_solvents,
            duplicate_fragments,
            True,
            "Parent fallback: salt stripping removed all non-solvent fragments",
        )

    kept_fragments = [info["fragment"] for info in kept_after_salt]
    kept_fragment_smiles = [_canonical_smiles_from_mol(fragment) for fragment in kept_fragments]
    duplicate_fragments = sorted(
        smiles for smiles, count in Counter(kept_fragment_smiles).items() if count > 1
    )
    unique_by_smiles = {}
    for fragment in kept_fragments:
        unique_by_smiles.setdefault(_canonical_smiles_from_mol(fragment), fragment)

    parent = _rebuild_parent_from_fragments(list(unique_by_smiles.values()))

    if removed_salts or removed_solvents:
        comment = "Parent extracted after stripping salts/solvents"
    else:
        comment = "Parent extracted without fragment stripping"

    return parent, removed_salts, removed_solvents, duplicate_fragments, False, comment


def _build_pipeline_outputs(
    smiles: str,
) -> tuple[str | None, str | None, str | None, str, str, list[str], list[str], list[str], bool]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None, "Invalid SMILES", "Invalid SMILES", [], [], [], False

    standardized_mol, standardization_comment = _standardize_molecule(mol)
    if standardized_mol is None:
        canonical_smiles = _canonical_smiles_from_mol(mol)
        return (
            canonical_smiles,
            canonical_smiles,
            canonical_smiles,
            standardization_comment,
            "Parent fallback: standardization failed",
            [],
            [],
            [],
            True,
        )

    standardized_smiles = _canonical_smiles_from_mol(standardized_mol)
    observed_parent_mol, removed_salts, removed_solvents, duplicate_fragments, fallback_triggered, parent_comment = _extract_parent_molecule(standardized_mol)
    observed_parent_smiles = _canonical_smiles_from_mol(observed_parent_mol)

    grouping_mol, grouping_comment = _neutralize_for_grouping(standardized_mol)
    if grouping_mol is None:
        grouping_mol = Chem.Mol(standardized_mol)
        grouping_fallback = True
        parent_comment = f"{parent_comment}; {grouping_comment}"
    else:
        grouping_fallback = False
    grouping_parent_mol, _, _, _, grouping_parent_fallback, grouping_parent_comment = _extract_parent_molecule(grouping_mol)
    parent_smiles = _canonical_smiles_from_mol(grouping_parent_mol)
    if grouping_parent_comment != "Single-fragment parent":
        parent_comment = f"{parent_comment}; {grouping_parent_comment}"

    return (
        standardized_smiles,
        observed_parent_smiles,
        parent_smiles,
        standardization_comment,
        parent_comment,
        removed_salts,
        removed_solvents,
        duplicate_fragments,
        fallback_triggered or grouping_fallback or grouping_parent_fallback,
    )


def process_single_smiles(smiles_data):
    """Process a single SMILES string for validation"""
    idx, smiles = smiles_data
    if not isinstance(smiles, str) or pd.isna(smiles):
        return idx, None, None, True, False, None, None, None, "Invalid SMILES", "Invalid SMILES", [], [], [], False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return idx, None, None, True, False, None, None, None, "Invalid SMILES", "Invalid SMILES", [], [], [], False
    try:
        canonical_smiles = _canonical_smiles_from_mol(mol)
        is_non_canonical = canonical_smiles != smiles
        (
            standardized_smiles,
            observed_parent_smiles,
            parent_smiles,
            standardization_comment,
            parent_comment,
            removed_salts,
            removed_solvents,
            duplicate_fragments,
            parent_fallback,
        ) = _build_pipeline_outputs(smiles)
        return (
            idx,
            mol,
            canonical_smiles,
            False,
            is_non_canonical,
            standardized_smiles,
            observed_parent_smiles,
            parent_smiles,
            standardization_comment,
            parent_comment,
            removed_salts,
            removed_solvents,
            duplicate_fragments,
            parent_fallback,
        )
    except Exception as exc:
        return idx, mol, None, False, True, None, None, None, f"Processing failed: {exc}", "Parent not generated", [], [], [], True


def process_mol_for_consistency(mol_data):
    """Process a molecule for representation consistency check"""
    idx, smiles = mol_data
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return idx, smiles, smiles, 0
        normalizer = MolStandardize.normalize.Normalizer()
        uncharger = MolStandardize.charge.Uncharger()
        normalized_mol = normalizer.normalize(mol)
        neutral_mol = uncharger.uncharge(normalized_mol)
        neutral_smiles = Chem.MolToSmiles(neutral_mol)
        formal_charge = Chem.rdmolops.GetFormalCharge(mol)
        return idx, smiles, neutral_smiles, formal_charge
    except Exception:
        return idx, smiles, smiles, 0


def _summarize_molecule_stereochemistry(
    mol: Chem.Mol,
) -> tuple[int, int, int, int]:
    """Return chiral-center and stereogenic-double-bond counts for a molecule."""
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    num_chiral_centers = len(chiral_centers)
    undefined_chiral_centers = sum(1 for _, flag in chiral_centers if flag == '?')

    num_stereogenic_double_bonds = 0
    undefined_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue

        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetSymbol() != 'C' or atom2.GetSymbol() != 'C':
            continue

        if atom1.GetDegree() < 2 or atom2.GetDegree() < 2:
            continue

        num_stereogenic_double_bonds += 1
        stereo = bond.GetStereo()
        if stereo not in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
            undefined_double_bonds += 1

    return (
        num_chiral_centers,
        undefined_chiral_centers,
        num_stereogenic_double_bonds,
        undefined_double_bonds,
    )


def _empty_component_stereochemistry_summary(smiles: str | None = None) -> dict:
    return {
        "basis_smiles": smiles or "",
        "num_chiral_centers": 0,
        "undefined_chiral_centers": 0,
        "num_stereogenic_double_bonds": 0,
        "undefined_double_bonds": 0,
    }


def _is_better_component_candidate(
    candidate: dict,
    best: dict | None,
    undefined_key: str,
    total_key: str,
) -> bool:
    if best is None:
        return True
    if candidate[undefined_key] != best[undefined_key]:
        return candidate[undefined_key] > best[undefined_key]
    if candidate[total_key] != best[total_key]:
        return candidate[total_key] > best[total_key]
    return candidate["basis_smiles"] < best["basis_smiles"]


def _select_component_stereochemistry_summary(smiles: str | None) -> dict:
    """Choose one unique disconnected component as the detail-count basis."""
    if not smiles:
        return {
            "chirality": _empty_component_stereochemistry_summary(),
            "double_bond": _empty_component_stereochemistry_summary(),
        }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "chirality": _empty_component_stereochemistry_summary(smiles),
            "double_bond": _empty_component_stereochemistry_summary(smiles),
        }

    fragments = list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False))
    if not fragments:
        fragments = [mol]

    unique_component_smiles = sorted({_canonical_smiles_from_mol(fragment) for fragment in fragments})
    chirality_best = None
    double_bond_best = None

    for component_smiles in unique_component_smiles:
        component_mol = Chem.MolFromSmiles(component_smiles)
        if component_mol is None:
            continue

        (
            num_chiral_centers,
            undefined_chiral_centers,
            num_stereogenic_double_bonds,
            undefined_double_bonds,
        ) = _summarize_molecule_stereochemistry(component_mol)
        candidate = {
            "basis_smiles": component_smiles,
            "num_chiral_centers": num_chiral_centers,
            "undefined_chiral_centers": undefined_chiral_centers,
            "num_stereogenic_double_bonds": num_stereogenic_double_bonds,
            "undefined_double_bonds": undefined_double_bonds,
        }

        if _is_better_component_candidate(
            candidate,
            chirality_best,
            undefined_key="undefined_chiral_centers",
            total_key="num_chiral_centers",
        ):
            chirality_best = candidate
        if _is_better_component_candidate(
            candidate,
            double_bond_best,
            undefined_key="undefined_double_bonds",
            total_key="num_stereogenic_double_bonds",
        ):
            double_bond_best = candidate

    return {
        "chirality": chirality_best or _empty_component_stereochemistry_summary(smiles),
        "double_bond": double_bond_best or _empty_component_stereochemistry_summary(smiles),
    }


def check_mol_chirality(mol_data):
    """Check atom-centered chirality and double-bond stereochemistry for a molecule."""
    if len(mol_data) == 3:
        idx, smiles, detail_smiles = mol_data
    else:
        idx, smiles = mol_data
        detail_smiles = smiles

    default_result = {
        "idx": idx,
        "record_num_chiral_centers": 0,
        "record_undefined_chiral_centers": 0,
        "record_has_undefined_chiral": False,
        "record_num_stereogenic_double_bonds": 0,
        "record_undefined_double_bonds": 0,
        "record_has_undefined_double_bond": False,
        "detail_num_chiral_centers": 0,
        "detail_undefined_chiral_centers": 0,
        "detail_chirality_basis_smiles": "",
        "detail_chirality_basis_source": "observed_parent_smiles",
        "detail_num_stereogenic_double_bonds": 0,
        "detail_undefined_double_bonds": 0,
        "detail_double_bond_basis_smiles": "",
        "detail_double_bond_basis_source": "observed_parent_smiles",
    }
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return default_result

        (
            record_num_chiral_centers,
            record_undefined_chiral_centers,
            record_num_stereogenic_double_bonds,
            record_undefined_double_bonds,
        ) = _summarize_molecule_stereochemistry(mol)

        detail_component_summary = _select_component_stereochemistry_summary(detail_smiles)
        record_component_summary = (
            detail_component_summary
            if detail_smiles == smiles
            else _select_component_stereochemistry_summary(smiles)
        )

        detail_chirality = detail_component_summary["chirality"]
        detail_chirality_source = "observed_parent_smiles"
        if (
            record_undefined_chiral_centers > 0
            and detail_chirality["undefined_chiral_centers"] == 0
        ):
            detail_chirality = record_component_summary["chirality"]
            detail_chirality_source = "canonical_smiles"

        detail_double_bond = detail_component_summary["double_bond"]
        detail_double_bond_source = "observed_parent_smiles"
        if (
            record_undefined_double_bonds > 0
            and detail_double_bond["undefined_double_bonds"] == 0
        ):
            detail_double_bond = record_component_summary["double_bond"]
            detail_double_bond_source = "canonical_smiles"

        return {
            "idx": idx,
            "record_num_chiral_centers": record_num_chiral_centers,
            "record_undefined_chiral_centers": record_undefined_chiral_centers,
            "record_has_undefined_chiral": record_undefined_chiral_centers > 0,
            "record_num_stereogenic_double_bonds": record_num_stereogenic_double_bonds,
            "record_undefined_double_bonds": record_undefined_double_bonds,
            "record_has_undefined_double_bond": record_undefined_double_bonds > 0,
            "detail_num_chiral_centers": detail_chirality["num_chiral_centers"],
            "detail_undefined_chiral_centers": detail_chirality["undefined_chiral_centers"],
            "detail_chirality_basis_smiles": detail_chirality["basis_smiles"],
            "detail_chirality_basis_source": detail_chirality_source,
            "detail_num_stereogenic_double_bonds": detail_double_bond["num_stereogenic_double_bonds"],
            "detail_undefined_double_bonds": detail_double_bond["undefined_double_bonds"],
            "detail_double_bond_basis_smiles": detail_double_bond["basis_smiles"],
            "detail_double_bond_basis_source": detail_double_bond_source,
        }
    except Exception:
        return default_result


def calculate_mol_fingerprint(smiles):
    """Calculate Morgan fingerprint for a molecule"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        return None
    except Exception:
        return None


def calculate_mol_scaffold(smiles):
    """Calculate Murcko scaffold for a molecule"""
    try:
        return MurckoScaffoldSmilesFromSmiles(smiles)
    except Exception:
        return None


def calculate_qed_value(smiles):
    """Calculate QED value for a molecule"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return QED.qed(mol)
        return None
    except Exception:
        return None
