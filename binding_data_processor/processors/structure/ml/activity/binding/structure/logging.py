"""Logging module for structure analysis system.

This module provides:
1. Logging configuration
2. Custom formatters
3. Context managers
4. Debugging utilities
"""

import logging
import sys
import time
import traceback
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Dict, List, Optional, Any, TextIO
from pathlib import Path

# Configure root logger
logger = logging.getLogger(__name__)


@dataclass
class LogContext:
    """Context information for log messages."""

    component: str
    operation: str
    details: Optional[Dict[str, Any]] = None
    start_time: Optional[float] = None
    end_time: Optional[float] = None

    @property
    def duration(self) -> Optional[float]:
        """Get operation duration in seconds."""
        if self.start_time and self.end_time:
            return self.end_time - self.start_time
        return None


class StructureFormatter(logging.Formatter):
    """Custom formatter for structure analysis logs."""

    def __init__(self, include_context: bool = True):
        """Initialize formatter.

        Args:
            include_context: Whether to include context information
        """
        super().__init__()
        self.include_context = include_context

    def format(self, record: logging.LogRecord) -> str:
        """Format log record.

        Args:
            record: Log record to format

        Returns:
            Formatted log message
        """
        # Basic message format
        msg = f"[{record.levelname}] {record.getMessage()}"

        # Add context if present
        if self.include_context and hasattr(record, "context"):
            context: LogContext = record.context
            msg = f"[{context.component}:{context.operation}] {msg}"
            if context.duration is not None:
                msg = f"{msg} (took {context.duration:.3f}s)"
            if context.details:
                msg = f"{msg}\nDetails: {context.details}"

        # Add exception info if present
        if record.exc_info:
            msg = f"{msg}\n{''.join(traceback.format_exception(*record.exc_info))}"

        return msg


class StructureLogger:
    """Logger for structure analysis system."""

    def __init__(
        self,
        name: str,
        level: int = logging.INFO,
        file: Optional[TextIO] = None,
        include_context: bool = True,
    ):
        """Initialize logger.

        Args:
            name: Logger name
            level: Logging level
            file: Optional file to write logs to
            include_context: Whether to include context in log messages
        """
        self.logger = logging.getLogger(name)
        self.logger.setLevel(level)

        # Create handlers
        handlers = []

        # Console handler
        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(StructureFormatter(include_context=include_context))
        handlers.append(console)

        # File handler if specified
        if file:
            file_handler = logging.FileHandler(file)
            file_handler.setFormatter(StructureFormatter(include_context=include_context))
            handlers.append(file_handler)

        # Add handlers
        for handler in handlers:
            self.logger.addHandler(handler)

        self._context_stack: List[LogContext] = []

    def debug(self, msg: str, *args, **kwargs):
        """Log debug message."""
        self._log(logging.DEBUG, msg, *args, **kwargs)

    def info(self, msg: str, *args, **kwargs):
        """Log info message."""
        self._log(logging.INFO, msg, *args, **kwargs)

    def warning(self, msg: str, *args, **kwargs):
        """Log warning message."""
        self._log(logging.WARNING, msg, *args, **kwargs)

    def error(self, msg: str, *args, **kwargs):
        """Log error message."""
        self._log(logging.ERROR, msg, *args, **kwargs)

    def critical(self, msg: str, *args, **kwargs):
        """Log critical message."""
        self._log(logging.CRITICAL, msg, *args, **kwargs)

    def _log(self, level: int, msg: str, *args, **kwargs):
        """Internal logging method.

        Args:
            level: Logging level
            msg: Message to log
            *args: Additional positional arguments
            **kwargs: Additional keyword arguments
        """
        # Add current context if available
        if self._context_stack:
            kwargs["extra"] = kwargs.get("extra", {})
            kwargs["extra"]["context"] = self._context_stack[-1]

        self.logger.log(level, msg, *args, **kwargs)

    @contextmanager
    def context(
        self,
        component: str,
        operation: str,
        details: Optional[Dict[str, Any]] = None,
    ):
        """Context manager for logging operations.

        Args:
            component: Component name
            operation: Operation name
            details: Optional operation details
        """
        context = LogContext(
            component=component,
            operation=operation,
            details=details,
            start_time=time.time(),
        )
        self._context_stack.append(context)

        try:
            yield context
        except Exception as e:
            self.error(f"Operation failed: {str(e)}", exc_info=True)
            raise
        finally:
            context.end_time = time.time()
            self._context_stack.pop()


def create_logger(
    name: str = "structure",
    level: int = logging.INFO,
    log_dir: Optional[Path] = None,
    include_context: bool = True,
) -> StructureLogger:
    """Create structure logger.

    Args:
        name: Logger name
        level: Logging level
        log_dir: Optional directory for log files
        include_context: Whether to include context in log messages

    Returns:
        Configured StructureLogger instance
    """
    # Create log file if directory specified
    file = None
    if log_dir:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        file = log_dir / f"{name}.log"

    return StructureLogger(
        name=name,
        level=level,
        file=file.open("a") if file else None,
        include_context=include_context,
    )


@contextmanager
def log_operation(
    logger: StructureLogger,
    component: str,
    operation: str,
    details: Optional[Dict[str, Any]] = None,
    level: int = logging.INFO,
):
    """Context manager for logging operations.

    Args:
        logger: Logger instance
        component: Component name
        operation: Operation name
        details: Optional operation details
        level: Logging level for start/end messages
    """
    with logger.context(component, operation, details) as context:
        logger._log(level, f"Starting {operation}")
        try:
            yield context
            logger._log(
                level,
                f"Completed {operation}",
                extra={"context": context},
            )
        except Exception as e:
            logger.error(
                f"Failed {operation}: {str(e)}",
                exc_info=True,
                extra={"context": context},
            )
            raise


def get_logger(name: str = "structure") -> StructureLogger:
    """Get or create structure logger.

    Args:
        name: Logger name

    Returns:
        StructureLogger instance
    """
    if name in logging.root.manager.loggerDict:
        logger = logging.getLogger(name)
        if isinstance(logger, StructureLogger):
            return logger
    return create_logger(name)
