from .._common import *


class FingerprintBatchMixin:
    def compute_fingerprints_batch(self, mols):
        """Compute fingerprints in batch for better performance"""
        fps = []
        for mol in mols:
            try:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fps.append(fp)
            except Exception:
                fps.append(None)
        return fps

