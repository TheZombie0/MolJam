"""ScoringRunner: composition-based alternative to the mixin orchestration."""

from typing import List

from .check_protocol import ScoringCheck
from .config import DEFAULT_CONFIG, ScoringConfig
from .context import ScoringContext


class ScoringRunner:
    """Run a sequence of ScoringCheck instances against a ScoringContext.

    This is the composition-based alternative to MoleculeDBScorer's mixin chain.
    Existing mixin methods can be wrapped as ScoringCheck adapters.
    """

    def __init__(
        self,
        checks: List[ScoringCheck],
        config: ScoringConfig = DEFAULT_CONFIG,
    ):
        self.checks = checks
        self.config = config

    def run_all(self, ctx: ScoringContext) -> float:
        """Execute all checks in order and return the sum of scores."""
        total = 0.0
        for check in self.checks:
            score = check.run(ctx)
            ctx.completed_checks.add(check.name)
            total += score
        return total
