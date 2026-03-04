from .data_quality_consistency import DataConsistencyReliabilityMixin
from .data_quality_labels import LabelConsistencyChecksMixin


class DataQualityChecksMixin(
    DataConsistencyReliabilityMixin,
    LabelConsistencyChecksMixin,
):
    pass

