from .distribution_balance import DataBalanceDistributionMixin
from .distribution_size import DataSizeChecksMixin


class DistributionChecksMixin(
    DataSizeChecksMixin,
    DataBalanceDistributionMixin,
):
    pass

