from .distribution_activity import ActivityDistributionPlotMixin
from .distribution_balance import DistributionBalancePlotMixin
from .distribution_class import ClassDistributionPlotMixin
from .distribution_sizes import DatabaseSizePlotMixin


class DistributionPlotMixin(
    DistributionBalancePlotMixin,
    ActivityDistributionPlotMixin,
    DatabaseSizePlotMixin,
    ClassDistributionPlotMixin,
):
    pass

