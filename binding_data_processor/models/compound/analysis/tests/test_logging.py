"""Tests for compound analysis logging functionality."""

import pytest
import logging
import tempfile
import json
from pathlib import Path
from datetime import datetime

from ..logging import (
    AnalysisLogger,
    LogLevel,
    LogContext,
    LogFilter,
    format_log_message,
    setup_logging,
)


@pytest.fixture
def temp_log_file():
    """Create temporary log file."""
    with tempfile.NamedTemporaryFile(suffix=".log", delete=False) as f:
        yield Path(f.name)
        Path(f.name).unlink()


def test_log_levels():
    """Test log level handling."""
    # Test log level enum
    assert LogLevel.DEBUG.value == logging.DEBUG
    assert LogLevel.INFO.value == logging.INFO
    assert LogLevel.WARNING.value == logging.WARNING
    assert LogLevel.ERROR.value == logging.ERROR
    assert LogLevel.CRITICAL.value == logging.CRITICAL

    # Test level comparison
    assert LogLevel.DEBUG < LogLevel.INFO
    assert LogLevel.INFO < LogLevel.WARNING
    assert LogLevel.WARNING < LogLevel.ERROR
    assert LogLevel.ERROR < LogLevel.CRITICAL

    # Test level names
    assert LogLevel.DEBUG.name == "DEBUG"
    assert LogLevel.INFO.name == "INFO"
    assert LogLevel.WARNING.name == "WARNING"
    assert LogLevel.ERROR.name == "ERROR"
    assert LogLevel.CRITICAL.name == "CRITICAL"


def test_log_formatting():
    """Test log message formatting."""
    # Test basic message
    msg = format_log_message("Test message")
    assert isinstance(msg, str)
    assert "Test message" in msg
    assert datetime.now().strftime("%Y-%m-%d") in msg

    # Test with context
    context = {"value": 1.2, "unit": "nM"}
    msg = format_log_message("Test message", context)
    assert "value: 1.2" in msg
    assert "unit: nM" in msg

    # Test with level
    msg = format_log_message("Test message", level=LogLevel.ERROR)
    assert "[ERROR]" in msg
    assert "Test message" in msg


def test_log_context():
    """Test log context management."""
    # Test context creation
    context = LogContext()
    assert context.data == {}

    # Test context update
    context.update({"value": 1.2})
    assert context.data["value"] == 1.2

    # Test context manager
    with LogContext() as ctx:
        ctx.update({"test": True})
        assert ctx.data["test"] is True

    # Test nested context
    with LogContext() as outer:
        outer.update({"outer": True})
        with LogContext(parent=outer) as inner:
            inner.update({"inner": True})
            assert "outer" in inner.data
            assert "inner" in inner.data


def test_log_filter():
    """Test log filtering functionality."""
    # Create filter
    log_filter = LogFilter(
        min_level=LogLevel.INFO,
        include_patterns=["*ERROR*", "*WARNING*"],
        exclude_patterns=["*DEBUG*"]
    )

    # Test level filtering
    assert log_filter.should_log("Test message", LogLevel.INFO)
    assert not log_filter.should_log("Test message", LogLevel.DEBUG)

    # Test pattern filtering
    assert log_filter.should_log("ERROR: Test message", LogLevel.INFO)
    assert log_filter.should_log("WARNING: Test message", LogLevel.INFO)
    assert not log_filter.should_log("DEBUG: Test message", LogLevel.INFO)

    # Test combined filtering
    assert not log_filter.should_log("ERROR: Test message", LogLevel.DEBUG)
    assert log_filter.should_log("WARNING: Test message", LogLevel.WARNING)


def test_logger_setup(temp_log_file):
    """Test logger setup and configuration."""
    # Setup logger
    logger = setup_logging(
        log_file=temp_log_file,
        min_level=LogLevel.INFO,
        include_patterns=["*ERROR*", "*WARNING*"],
        exclude_patterns=["*DEBUG*"]
    )

    # Test logger configuration
    assert isinstance(logger, AnalysisLogger)
    assert logger.min_level == LogLevel.INFO
    assert "*ERROR*" in logger.include_patterns
    assert "*DEBUG*" in logger.exclude_patterns

    # Test log file creation
    assert temp_log_file.exists()
    assert temp_log_file.stat().st_size == 0


def test_logger_functionality(temp_log_file):
    """Test logger core functionality."""
    # Setup logger
    logger = setup_logging(log_file=temp_log_file)

    # Test logging at different levels
    logger.debug("Debug message")
    logger.info("Info message")
    logger.warning("Warning message")
    logger.error("Error message")
    logger.critical("Critical message")

    # Check log file content
    content = temp_log_file.read_text()
    assert "Debug message" in content
    assert "Info message" in content
    assert "Warning message" in content
    assert "Error message" in content
    assert "Critical message" in content


def test_structured_logging(temp_log_file):
    """Test structured logging functionality."""
    # Setup logger
    logger = setup_logging(
        log_file=temp_log_file,
        structured=True
    )

    # Log with context
    context = {"value": 1.2, "unit": "nM"}
    logger.info("Test message", context)

    # Check log file content
    content = temp_log_file.read_text()
    log_entry = json.loads(content.strip().split("\n")[-1])
    
    assert log_entry["message"] == "Test message"
    assert log_entry["context"]["value"] == 1.2
    assert log_entry["context"]["unit"] == "nM"
    assert "timestamp" in log_entry
    assert "level" in log_entry


def test_log_rotation(temp_log_file):
    """Test log file rotation."""
    # Setup logger with rotation
    logger = setup_logging(
        log_file=temp_log_file,
        max_size=1024,  # 1KB
        backup_count=3
    )

    # Generate enough logs to trigger rotation
    large_message = "X" * 512  # 512 bytes
    for _ in range(10):  # Should create multiple log files
        logger.info(large_message)

    # Check rotated files
    log_dir = temp_log_file.parent
    rotated_files = list(log_dir.glob(f"{temp_log_file.stem}.*"))
    assert len(rotated_files) > 0
    assert len(rotated_files) <= 3  # backup_count


def test_log_analysis():
    """Test log analysis functionality."""
    # Create in-memory log stream
    logger = setup_logging(use_memory=True)

    # Generate test logs
    logger.info("Test info")
    logger.warning("Test warning")
    logger.error("Test error")

    # Analyze logs
    analysis = logger.analyze_logs()
    
    assert analysis["total_entries"] == 3
    assert analysis["by_level"]["INFO"] == 1
    assert analysis["by_level"]["WARNING"] == 1
    assert analysis["by_level"]["ERROR"] == 1
    assert len(analysis["error_messages"]) == 1
    assert "Test error" in analysis["error_messages"][0]


def test_log_cleanup():
    """Test log cleanup functionality."""
    # Create temporary log directory
    with tempfile.TemporaryDirectory() as temp_dir:
        log_dir = Path(temp_dir)
        
        # Create multiple log files
        for i in range(5):
            log_file = log_dir / f"test_{i}.log"
            log_file.write_text(f"Test log {i}")

        # Clean old logs
        cleanup_result = AnalysisLogger.cleanup_logs(
            log_dir,
            max_age_days=0,
            dry_run=False
        )

        assert cleanup_result["scanned"] == 5
        assert cleanup_result["deleted"] == 5
        assert len(list(log_dir.glob("*.log"))) == 0
