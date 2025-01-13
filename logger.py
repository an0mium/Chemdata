"""Logging configuration and management.

This module provides centralized logging configuration and management
through the LogManager class.
"""

import logging
import os
from typing import Optional


class LogManager:
    """Manages logging configuration and provides logger instances."""

    _instance = None
    _initialized = False

    def __new__(cls):
        """Ensure singleton instance."""
        if cls._instance is None:
            cls._instance = super(LogManager, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        """Initialize logging configuration if not already done."""
        if not LogManager._initialized:
            # Configure basic logging
            logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO)
            LogManager._initialized = True

    def get_logger(self, name: str) -> logging.Logger:
        """Get a logger instance with the specified name.

        Args:
            name: Name for the logger, typically the module name

        Returns:
            Configured logger instance
        """
        logger = logging.getLogger(name)

        # Add file handler if log directory exists
        log_dir = os.path.join(os.getcwd(), "logs")
        if os.path.exists(log_dir):
            log_file = os.path.join(log_dir, f"{name}.log")
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
            logger.addHandler(file_handler)

        return logger
