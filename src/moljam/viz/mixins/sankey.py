from .sankey_echarts import SankeyEchartsPlotMixin
from .sankey_matplotlib import SankeyMatplotlibFallbackMixin


class SankeyPlotMixin(
    SankeyEchartsPlotMixin,
    SankeyMatplotlibFallbackMixin,
):
    pass

