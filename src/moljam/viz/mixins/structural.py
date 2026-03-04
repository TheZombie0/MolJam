from .structural_chirality import ChiralityPlotMixin
from .structural_contradictions import ContradictoryLabelPlotMixin
from .structural_duplication import StructuralDuplicationPlotMixin
from .structural_smiles_quality import SmilesQualityPlotMixin


class StructuralPlotMixin(
    ChiralityPlotMixin,
    StructuralDuplicationPlotMixin,
    ContradictoryLabelPlotMixin,
    SmilesQualityPlotMixin,
):
    pass

