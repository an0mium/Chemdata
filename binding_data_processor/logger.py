"""Logging configuration and management.

This module provides centralized logging configuration and management
through the LogManager singleton class with:
1. Consistent logging configuration across the package
2. Custom log formatting
3. Log file handling
4. Log level management
5. Module-specific logging
6. Singleton pattern to ensure consistent configuration
"""

import logging
import os
from pathlib import Path
from typing import Optional


class LogManager:
    """Manages logging configuration and provides logger instances."""

    _instance = None
    _initialized = False

    def __new__(cls, *args, **kwargs):
        """Ensure singleton instance."""
        if cls._instance is None:
            cls._instance = super(LogManager, cls).__new__(cls)
        return cls._instance

    def __init__(
        self,
        log_dir: Optional[str] = None,
        log_level: int = logging.INFO,
        console_format: Optional[str] = None,
        file_format: Optional[str] = None,
    ):
        """Initialize logging configuration if not already done.

        Args:
            log_dir: Optional directory for log files (default: "logs")
            log_level: Logging level (default: INFO)
            console_format: Optional custom console log format
            file_format: Optional custom file log format
        """
        # Only initialize once due to singleton pattern
        if LogManager._initialized:
            return

        self.log_dir = Path(log_dir) if log_dir else Path("logs")
        self.log_level = log_level
        self.console_format = console_format or "%(levelname)s - %(message)s"
        self.file_format = file_format or "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

        # Create log directory if it doesn't exist
        os.makedirs(self.log_dir, exist_ok=True)

        # Configure root logger
        self._configure_root_logger()

        # Create logger for this module
        self.logger = logging.getLogger(__name__)

        LogManager._initialized = True

    def _configure_root_logger(self) -> None:
        """Configure the root logger with file and console handlers."""
        root_logger = logging.getLogger()
        root_logger.setLevel(self.log_level)

        # Remove any existing handlers
        for handler in root_logger.handlers[:]:
            root_logger.removeHandler(handler)

        # Create formatters
        console_formatter = logging.Formatter(self.console_format)
        file_formatter = logging.Formatter(self.file_format)

        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(self.log_level)
        console_handler.setFormatter(console_formatter)
        root_logger.addHandler(console_handler)

        # File handler for main log file
        main_log_file = self.log_dir / "binding_data_processor.log"
        file_handler = logging.FileHandler(main_log_file)
        file_handler.setLevel(self.log_level)
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)

    def get_logger(self, name: str) -> logging.Logger:
        """Get a logger with the specified name.

        Args:
            name: Logger name (typically __name__ of the calling module)

        Returns:
            Logger instance configured with both file and console handlers
        """
        logger = logging.getLogger(name)

        # Add module-specific file handler if it doesn't exist
        has_file_handler = any(isinstance(h, logging.FileHandler) for h in logger.handlers)
        if not has_file_handler:
            module_log_file = self.log_dir / f"{name.replace('.', '_')}.log"
            file_handler = logging.FileHandler(module_log_file)
            file_handler.setLevel(self.log_level)
            file_handler.setFormatter(logging.Formatter(self.file_format))
            logger.addHandler(file_handler)

        return logger

    def set_level(self, level: int) -> None:
        """Set logging level for all handlers.

        Args:
            level: New logging level (e.g. logging.INFO, logging.DEBUG)
        """
        self.log_level = level
        root_logger = logging.getLogger()
        root_logger.setLevel(level)
        for handler in root_logger.handlers:
            handler.setLevel(level)

        # Update level for all module loggers
        for name in logging.root.manager.loggerDict:
            logger = logging.getLogger(name)
            logger.setLevel(level)
            for handler in logger.handlers:
                handler.setLevel(level)

    def add_file_handler(
        self,
        filename: str,
        level: Optional[int] = None,
        format_str: Optional[str] = None,
    ) -> None:
        """Add an additional file handler.

        Args:
            filename: Name of log file
            level: Optional logging level (defaults to manager's level)
            format_str: Optional log format (defaults to manager's file format)
        """
        handler = logging.FileHandler(self.log_dir / filename)
        handler.setLevel(level or self.log_level)
        handler.setFormatter(logging.Formatter(format_str or self.file_format))
        logging.getLogger().addHandler(handler)

    def remove_file_handler(self, filename: str) -> None:
        """Remove a specific file handler.

        Args:
            filename: Name of log file to remove handler for
        """
        root_logger = logging.getLogger()
        filepath = self.log_dir / filename

        for handler in root_logger.handlers[:]:
            if isinstance(handler, logging.FileHandler) and handler.baseFilename == str(filepath):
                root_logger.removeHandler(handler)
                handler.close()

    def update_formats(
        self,
        console_format: Optional[str] = None,
        file_format: Optional[str] = None,
    ) -> None:
        """Update log formats for all handlers.

        Args:
            console_format: New console log format
            file_format: New file log format
        """
        if console_format:
            self.console_format = console_format
        if file_format:
            self.file_format = file_format

        root_logger = logging.getLogger()
        for handler in root_logger.handlers:
            if isinstance(handler, logging.StreamHandler) and not isinstance(handler, logging.FileHandler):
                handler.setFormatter(logging.Formatter(self.console_format))
            elif isinstance(handler, logging.FileHandler):
                handler.setFormatter(logging.Formatter(self.file_format))
