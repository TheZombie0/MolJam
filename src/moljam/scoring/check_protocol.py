"""Protocol for scoring checks (composition pattern)."""

from typing import Protocol, runtime_checkable

from .context import ScoringContext


@runtime_checkable
class ScoringCheck(Protocol):
    """Interface that all new-style scoring checks must implement."""

    name: str

    def run(self, ctx: ScoringContext) -> float:
        """Execute the check, mutating ctx.scores and ctx.analysis_results.

        Returns the score (0-10) for this check.
        """
        ...
