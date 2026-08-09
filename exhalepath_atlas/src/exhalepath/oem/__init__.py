"""OEM embed kit — score-sample + locked-split verify for LIMS / partner stacks."""

from .kit import (
    import_feature_table,
    score_sample,
    verify_locked_split,
)

__all__ = ["import_feature_table", "score_sample", "verify_locked_split"]
