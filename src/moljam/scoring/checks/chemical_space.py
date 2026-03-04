from .chemical_space_diversity import ChemicalDiversityChecksMixin
from .chemical_space_druglikeness import DrugLikenessChecksMixin
from .chemical_space_fingerprints import FingerprintBatchMixin


class ChemicalSpaceChecksMixin(
    FingerprintBatchMixin,
    ChemicalDiversityChecksMixin,
    DrugLikenessChecksMixin,
):
    pass

