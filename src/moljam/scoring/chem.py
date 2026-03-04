"""
RDKit-centric helpers used by the scoring pipeline.

These functions are kept at module scope so they can be used with
`multiprocessing.Pool`.
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, QED
from rdkit.Chem import MolStandardize
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmilesFromSmiles


def process_single_smiles(smiles_data):
    """Process a single SMILES string for validation"""
    idx, smiles = smiles_data
    if not isinstance(smiles, str) or pd.isna(smiles):
        return idx, None, None, True, False  # idx, mol, canonical_smiles, is_invalid, is_non_canonical
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return idx, None, None, True, False
    try:
        canonical_smiles = Chem.MolToSmiles(mol)
        is_non_canonical = canonical_smiles != smiles
        return idx, mol, canonical_smiles, False, is_non_canonical
    except Exception:
        return idx, mol, None, False, True


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


def check_mol_chirality(mol_data):
    """Check chirality for a single molecule"""
    idx, smiles = mol_data
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return idx, 0, 0, False
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if len(chiral_centers) > 0:
            undefined = sum(1 for _, flag in chiral_centers if flag == '?')
            return idx, len(chiral_centers), undefined, undefined > 0
        return idx, 0, 0, False
    except Exception:
        return idx, 0, 0, False


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

