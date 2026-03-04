from .mixins import (
    AnnotationQualityPlotMixin,
    BarsPlotMixin,
    ChemicalPlotMixin,
    DistributionPlotMixin,
    QualityPlotMixin,
    RadarPlotMixin,
    ReportingMixin,
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
    StructuralPlotMixin,
    ChemicalPlotMixin,
    DistributionPlotMixin,
    AnnotationQualityPlotMixin,
    UsefulColumnsPlotMixin,
    ReportingMixin,
):
    pass
