from .annotation_experimental import ExperimentalMethodChecksMixin
from .annotation_quality_check import AnnotationQualityChecksMixin
from .annotation_time_label import TimeLabelAvailabilityChecksMixin


class AnnotationChecksMixin(
    ExperimentalMethodChecksMixin,
    TimeLabelAvailabilityChecksMixin,
    AnnotationQualityChecksMixin,
):
    pass

