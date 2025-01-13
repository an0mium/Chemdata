"""Core data models for compound processing."""

from dataclasses import dataclass
from enum import Enum, auto
from typing import Any, Dict, Optional


class ErrorCode(Enum):
    """Error codes for MCP operations."""

    InvalidRequest = auto()
    MethodNotFound = auto()
    InvalidParams = auto()
    InternalError = auto()
    ServiceUnavailable = auto()
    RateLimitExceeded = auto()
    AuthenticationError = auto()
    AuthorizationError = auto()
    ResourceNotFound = auto()
    ResourceConflict = auto()
    ValidationError = auto()
    CircuitBreakerOpen = auto()


class McpError(Exception):
    """Base exception class for MCP errors."""

    def __init__(self, code: ErrorCode, message: str):
        """Initialize MCP error.

        Args:
            code: Error code
            message: Error message
        """
        self.code = code
        self.message = message
        super().__init__(f"{code.name}: {message}")


@dataclass
class PredictionResult:
    """Data class for storing prediction results."""

    value: Any
    confidence: float
    supporting_data: Dict[str, Any]


@dataclass
class CompoundData:
    """Data class for storing compound information."""

    name: str
    smiles: str
    molecular_weight: Optional[float] = None
    logp: Optional[float] = None
    tpsa: Optional[float] = None
    hbd: Optional[int] = None
    hba: Optional[int] = None
    rotatable_bonds: Optional[int] = None
