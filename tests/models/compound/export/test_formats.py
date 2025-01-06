"""Tests for export formats module."""

import pytest
from binding_data_processor.models.compound.export import formats


@pytest.fixture
def mock_exporter():
    """Create a mock exporter for testing."""
    return formats.BaseExporter()


def test_exporter_initialization(mock_exporter):
    """Test exporter initialization."""
    assert isinstance(mock_exporter, formats.BaseExporter)
    assert hasattr(mock_exporter, "export")


def test_export_to_tsv(mock_exporter):
    """Test TSV export functionality."""
    with pytest.raises(NotImplementedError):
        mock_exporter.export_to_tsv(None)


def test_export_to_json(mock_exporter):
    """Test JSON export functionality."""
    with pytest.raises(NotImplementedError):
        mock_exporter.export_to_json(None)


def test_format_validation(mock_exporter):
    """Test format validation."""
    with pytest.raises(NotImplementedError):
        mock_exporter.validate_format("invalid_format")


def test_column_selection(mock_exporter):
    """Test column selection for export."""
    with pytest.raises(NotImplementedError):
        mock_exporter.select_columns(["col1", "col2"])
