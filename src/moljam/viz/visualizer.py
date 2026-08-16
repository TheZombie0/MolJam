from .mixins import (
    AnnotationQualityPlotMixin,
    BarsPlotMixin,
    ChemicalPlotMixin,
    DistributionPlotMixin,
    QualityPlotMixin,
    QualityTableMixin,
    RadarPlotMixin,
    ReportingMixin,
    RuntimePlotMixin,
    SankeyPlotMixin,
    StructuralPlotMixin,
    UsefulColumnsPlotMixin,
    VisualizerCoreMixin,
)


class MoleculeDBVisualizer(
    VisualizerCoreMixin,
    BarsPlotMixin,
    RadarPlotMixin,
    SankeyPlotMixin,
    QualityPlotMixin,
    QualityTableMixin,
    StructuralPlotMixin,
    ChemicalPlotMixin,
    DistributionPlotMixin,
    RuntimePlotMixin,
    AnnotationQualityPlotMixin,
    UsefulColumnsPlotMixin,
    ReportingMixin,
):
    pass
