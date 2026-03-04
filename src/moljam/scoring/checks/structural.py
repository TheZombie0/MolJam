from .structural_representation import RepresentationConsistencyMixin
from .structural_stereochemistry import StereochemistryChecksMixin
from .structural_validation import SmilesValidationMixin


class StructuralChecksMixin(
    SmilesValidationMixin,
    RepresentationConsistencyMixin,
    StereochemistryChecksMixin,
):
    pass

