from dataclasses import dataclass, field


@dataclass(frozen=True)
class PenaltyConfig:
    very_low_threshold: float = 2.0
    low_threshold: float = 4.0
    medium_threshold: float = 6.0
    very_low_penalty: float = 1.5
    low_penalty: float = 0.5
    medium_penalty: float = 0.1
    ratio_cutoff: float = 0.3
    additional_multiplier: float = 1.2
    total_cap: float = 30.0


@dataclass(frozen=True)
class NormalizationConfig:
    points_per_module_4: float = 25.0
    points_per_module_5: float = 20.0
    max_score_per_metric: float = 10.0


@dataclass(frozen=True)
class ScoringConfig:
    penalty: PenaltyConfig = field(default_factory=PenaltyConfig)
    normalization: NormalizationConfig = field(default_factory=NormalizationConfig)
    max_parallel_workers: int = 100


DEFAULT_CONFIG = ScoringConfig()
