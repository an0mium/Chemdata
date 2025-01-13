"""AlphaFold integration package."""

from .config import AlphaFoldConfig
from .predictor import AlphaFoldPredictor, AlphaFoldResult

__all__ = ["AlphaFoldConfig", "AlphaFoldPredictor", "AlphaFoldResult"]
