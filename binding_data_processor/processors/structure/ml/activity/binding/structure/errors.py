"""Error handling module for structure analysis system.

This module provides:
1. Custom exception classes
2. Error codes and messages
3. Error handling utilities
4. Debugging helpers
"""

import logging
from typing import Dict, List, Optional, Any, Type
from dataclasses import dataclass

logger = logging.getLogger(__name__)


class StructureError(Exception):
    """Base class for structure analysis errors."""

    def __init__(self, message: str, code: str = "STRUCTURE_ERROR"):
        self.message = message
        self.code = code
        super().__init__(self.message)


class ConfigurationError(StructureError):
    """Error in system configuration."""

    def __init__(self, message: str):
        super().__init__(message, code="CONFIG_ERROR")


class InitializationError(StructureError):
    """Error during system initialization."""

    def __init__(self, message: str):
        super().__init__(message, code="INIT_ERROR")


class ValidationError(StructureError):
    """Error in data validation."""

    def __init__(self, message: str):
        super().__init__(message, code="VALIDATION_ERROR")


class AnalysisError(StructureError):
    """Error during structure analysis."""

    def __init__(self, message: str):
        super().__init__(message, code="ANALYSIS_ERROR")


class PredictionError(StructureError):
    """Error during structure prediction."""

    def __init__(self, message: str):
        super().__init__(message, code="PREDICTION_ERROR")


class ComponentError(StructureError):
    """Error in component operation."""

    def __init__(self, component: str, message: str):
        self.component = component
        super().__init__(f"Error in {component}: {message}", code="COMPONENT_ERROR")


@dataclass
class ErrorContext:
    """Context information for errors."""

    error_type: Type[StructureError]
    message: str
    component: Optional[str] = None
    details: Optional[Dict[str, Any]] = None
    traceback: Optional[str] = None


class ErrorHandler:
    """Handler for structure analysis errors."""

    def __init__(self, debug: bool = False):
        """Initialize error handler.

        Args:
            debug: Whether to include debug information in errors
        """
        self.debug = debug
        self._error_contexts: List[ErrorContext] = []

    def handle_error(
        self,
        error: Exception,
        component: Optional[str] = None,
        details: Optional[Dict[str, Any]] = None,
    ) -> ErrorContext:
        """Handle an error.

        Args:
            error: Exception that occurred
            component: Component where error occurred
            details: Additional error details

        Returns:
            Error context
        """
        try:
            # Convert to structure error if needed
            if not isinstance(error, StructureError):
                error = self._convert_error(error)

            # Create error context
            context = ErrorContext(
                error_type=type(error),
                message=str(error),
                component=component,
                details=details,
                traceback=self._get_traceback() if self.debug else None,
            )

            # Log error
            self._log_error(context)

            # Store context
            self._error_contexts.append(context)

            return context

        except Exception as e:
            logger.error(f"Error in error handler: {str(e)}")
            return ErrorContext(
                error_type=StructureError,
                message="Error handling failed",
            )

    def get_last_error(self) -> Optional[ErrorContext]:
        """Get most recent error context.

        Returns:
            Most recent error context or None
        """
        return self._error_contexts[-1] if self._error_contexts else None

    def get_errors(self) -> List[ErrorContext]:
        """Get all error contexts.

        Returns:
            List of error contexts
        """
        return self._error_contexts.copy()

    def clear_errors(self):
        """Clear error history."""
        self._error_contexts.clear()

    def _convert_error(self, error: Exception) -> StructureError:
        """Convert generic exception to structure error.

        Args:
            error: Exception to convert

        Returns:
            Converted StructureError
        """
        if isinstance(error, ValueError):
            return ValidationError(str(error))
        if isinstance(error, RuntimeError):
            return AnalysisError(str(error))
        return StructureError(str(error))

    def _log_error(self, context: ErrorContext):
        """Log error context.

        Args:
            context: Error context to log
        """
        try:
            # Basic error info
            msg = f"Error [{context.error_type.__name__}]: {context.message}"
            if context.component:
                msg = f"{msg} (in {context.component})"

            # Add details if present
            if context.details:
                msg = f"{msg}\nDetails: {context.details}"

            # Add traceback in debug mode
            if self.debug and context.traceback:
                msg = f"{msg}\nTraceback:\n{context.traceback}"

            logger.error(msg)

        except Exception as e:
            logger.error(f"Error logging failed: {str(e)}")

    def _get_traceback(self) -> Optional[str]:
        """Get current exception traceback.

        Returns:
            Formatted traceback string or None
        """
        import traceback

        try:
            return "".join(traceback.format_exc())
        except Exception:
            return None


def handle_errors(debug: bool = False):
    """Decorator for error handling.

    Args:
        debug: Whether to include debug information

    Returns:
        Decorated function
    """

    def decorator(func):
        def wrapper(*args, **kwargs):
            handler = ErrorHandler(debug=debug)
            try:
                return func(*args, **kwargs)
            except Exception as e:
                context = handler.handle_error(e)
                if debug:
                    raise
                return None

        return wrapper

    return decorator


class ErrorCollection:
    """Collection of structure analysis errors."""

    def __init__(self):
        self.errors: List[ErrorContext] = []

    def add_error(
        self,
        message: str,
        error_type: Type[StructureError] = StructureError,
        component: Optional[str] = None,
        details: Optional[Dict[str, Any]] = None,
    ):
        """Add error to collection.

        Args:
            message: Error message
            error_type: Type of error
            component: Component where error occurred
            details: Additional error details
        """
        context = ErrorContext(
            error_type=error_type,
            message=message,
            component=component,
            details=details,
        )
        self.errors.append(context)

    def has_errors(self) -> bool:
        """Check if collection has errors.

        Returns:
            True if errors present
        """
        return len(self.errors) > 0

    def get_messages(self) -> List[str]:
        """Get all error messages.

        Returns:
            List of error messages
        """
        return [error.message for error in self.errors]

    def clear(self):
        """Clear all errors."""
        self.errors.clear()


def create_error_context(
    error: Exception,
    component: Optional[str] = None,
    details: Optional[Dict[str, Any]] = None,
    debug: bool = False,
) -> ErrorContext:
    """Create error context from exception.

    Args:
        error: Exception that occurred
        component: Component where error occurred
        details: Additional error details
        debug: Whether to include debug information

    Returns:
        Error context
    """
    handler = ErrorHandler(debug=debug)
    return handler.handle_error(error, component, details)
