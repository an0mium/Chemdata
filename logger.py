"""Logging configuration and utilities."""

import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

from config import LOG_DIR, LOG_FORMAT, LOG_LEVEL


class LogManager:
    """Manages application logging with file and console output."""
    
    _instance: Optional['LogManager'] = None
    
    def __new__(cls) -> 'LogManager':
        """Ensure singleton instance."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        """Initialize logging configuration if not already initialized."""
        if hasattr(self, '_initialized'):
            return
            
        self._initialized = True
        self.log_dir = LOG_DIR
        self.log_dir.mkdir(parents=True, exist_ok=True)
        
        # Create base logger
        self.logger = logging.getLogger('chemical_data_collector')
        self.logger.setLevel(LOG_LEVEL)
        
        # Prevent duplicate handlers
        if not self.logger.handlers:
            self._setup_handlers()

    def _setup_handlers(self):
        """Configure logging handlers for file and console output."""
        # Create formatters
        formatter = logging.Formatter(LOG_FORMAT)
        
        # File handler with daily rotation
        log_file = self.log_dir / f"chemical_data_{datetime.now():%Y%m%d}.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(LOG_LEVEL)
        
        # Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
        console_handler.setLevel(LOG_LEVEL)
        
        # Add handlers to logger
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)

    def get_logger(self, name: str = None) -> logging.Logger:
        """
        Get a logger instance.
        
        Args:
            name: Optional name for the logger (will be prefixed with base logger name)
            
        Returns:
            Configured logger instance
        """
        if name:
            return self.logger.getChild(name)
        return self.logger

    def archive_logs(self, days: int = 30) -> bool:
        """
        Archive logs older than specified days.
        
        Args:
            days: Number of days to keep logs (default: 30)
            
        Returns:
            True if successful, False on error
        """
        try:
            archive_dir = self.log_dir / "archive"
            archive_dir.mkdir(exist_ok=True)
            
            cutoff = datetime.now().timestamp() - (days * 24 * 60 * 60)
            
            for log_file in self.log_dir.glob("chemical_data_*.log"):
                if log_file.stat().st_mtime < cutoff:
                    archive_path = archive_dir / log_file.name
                    log_file.rename(archive_path)
                    
            return True
            
        except OSError as e:
            self.logger.error(f"Error archiving logs: {str(e)}")
            return False

    def get_log_stats(self) -> dict:
        """
        Get statistics about log files.
        
        Returns:
            Dictionary containing log statistics
        """
        log_files = list(self.log_dir.glob("chemical_data_*.log"))
        archive_files = list((self.log_dir / "archive").glob("chemical_data_*.log"))
        
        current_size = sum(f.stat().st_size for f in log_files)
        archive_size = sum(f.stat().st_size for f in archive_files)
        
        return {
            'current_logs': len(log_files),
            'archived_logs': len(archive_files),
            'current_size_bytes': current_size,
            'archive_size_bytes': archive_size,
            'total_size_bytes': current_size + archive_size
        }


# Global logger instance
logger = LogManager().get_logger()

# Convenience methods
def debug(msg: str, *args, **kwargs):
    """Log debug message."""
    logger.debug(msg, *args, **kwargs)

def info(msg: str, *args, **kwargs):
    """Log info message."""
    logger.info(msg, *args, **kwargs)

def warning(msg: str, *args, **kwargs):
    """Log warning message."""
    logger.warning(msg, *args, **kwargs)

def error(msg: str, *args, **kwargs):
    """Log error message."""
    logger.error(msg, *args, **kwargs)

def critical(msg: str, *args, **kwargs):
    """Log critical message."""
    logger.critical(msg, *args, **kwargs)

def exception(msg: str, *args, **kwargs):
    """Log exception with traceback."""
    logger.exception(msg, *args, **kwargs)
