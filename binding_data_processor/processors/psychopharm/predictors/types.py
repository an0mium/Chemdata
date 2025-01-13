"""Type definitions for predictors.

This module contains type definitions and data classes used across predictors
to avoid circular dependencies.
"""

from dataclasses import dataclass
from typing import Any, Dict, Optional


@dataclass
class PredictionResult:
    """Result from a prediction operation."""

    value: Any
    confidence: float
    supporting_data: Optional[Dict] = None
