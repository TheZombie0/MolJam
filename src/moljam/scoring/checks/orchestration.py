from .orchestration_penalty import LowScorePenaltyMixin
from .orchestration_runner import OrchestrationRunnerMixin
from .orchestration_totals import TotalScoreCalculationMixin


class OrchestrationMixin(
    LowScorePenaltyMixin,
    TotalScoreCalculationMixin,
    OrchestrationRunnerMixin,
):
    pass

