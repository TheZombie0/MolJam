from .chemical_activity import ActivityConsistencyPlotMixin
from .chemical_properties import MolecularPropertiesPlotMixin
from .chemical_scaffolds import ScaffoldPlotMixin
from .chemical_tsne import ChemicalSpaceTsnePlotMixin


class ChemicalPlotMixin(
    ScaffoldPlotMixin,
    ActivityConsistencyPlotMixin,
    MolecularPropertiesPlotMixin,
    ChemicalSpaceTsnePlotMixin,
):
    pass

