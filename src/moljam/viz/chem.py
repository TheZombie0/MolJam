from rdkit import Chem
from rdkit.Chem import Descriptors


def calculate_qed_parallel(smiles):
    """Calculate QED for a single SMILES"""
    try:
        from rdkit.Chem import QED

        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return QED.qed(mol)
    except Exception:
        pass
    return None


def calculate_properties_parallel(smiles):
    """Calculate molecular properties for a single SMILES"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return {
                'MW': Descriptors.MolWt(mol),
                'LogP': Descriptors.MolLogP(mol),
                'HBA': Descriptors.NumHAcceptors(mol),
                'HBD': Descriptors.NumHDonors(mol),
                'TPSA': Descriptors.TPSA(mol),
                'RotBonds': Descriptors.NumRotatableBonds(mol),
                'Lipinski': sum([
                    Descriptors.MolWt(mol) <= 500,
                    Descriptors.MolLogP(mol) <= 5,
                    Descriptors.NumHAcceptors(mol) <= 10,
                    Descriptors.NumHDonors(mol) <= 5
                ])
            }
    except Exception:
        pass
    return None


def calculate_fp_parallel(smiles):
    """Calculate Morgan fingerprint for a single SMILES (for t-SNE)"""
    try:
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            return list(fp)
    except Exception:
        pass
    return None

