"""Tests for data export functionality."""

import pytest
import pandas as pd
import json
from pathlib import Path
from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..data_export import (
    DataExporter,
    ExportConfig,
    ColumnSelector,
    FilterConfig,
    ExportResult,
    ExportFormat,
    TSVExporter,
    JSONExporter,
    BatchExporter,
)


@pytest.fixture
def test_compounds():
    """Create test compounds fixture."""
    compounds = []
    
    # Caffeine
    caffeine = PsychoactiveCompound(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    caffeine.psychoactive_class = PsychoactiveClass.STIMULANT
    caffeine.add_receptor_binding(
        "A2A",
        affinity=0.8,
        confidence=0.95,
        activity="antagonist"
    )
    caffeine.add_receptor_binding(
        "A1",
        affinity=0.7,
        confidence=0.9,
        activity="antagonist"
    )
    caffeine.effect_profile = {
        "stimulation": (0.8, 0.9),
        "focus": (0.7, 0.8),
        "wakefulness": (0.9, 0.95),
    }
    caffeine.safety_alerts = {
        "anxiety": RiskLevel.MODERATE,
        "insomnia": RiskLevel.HIGH,
        "tachycardia": RiskLevel.LOW,
    }
    compounds.append(caffeine)
    
    # Amphetamine
    amphetamine = PsychoactiveCompound(
        name="Amphetamine",
        smiles="CC(N)CC1=CC=CC=C1",
        cas_number="300-62-9",
    )
    amphetamine.psychoactive_class = PsychoactiveClass.STIMULANT
    amphetamine.add_receptor_binding(
        "DAT",
        affinity=0.05,
        confidence=0.95,
        activity="inhibitor"
    )
    amphetamine.add_receptor_binding(
        "NET",
        affinity=0.1,
        confidence=0.9,
        activity="inhibitor"
    )
    amphetamine.effect_profile = {
        "stimulation": (0.9, 0.95),
        "euphoria": (0.8, 0.9),
        "focus": (0.85, 0.9),
    }
    amphetamine.safety_alerts = {
        "addiction": RiskLevel.HIGH,
        "cardiovascular": RiskLevel.HIGH,
        "neurotoxicity": RiskLevel.MODERATE,
    }
    compounds.append(amphetamine)
    
    return compounds


@pytest.fixture
def export_config():
    """Create test export configuration fixture."""
    return ExportConfig(
        format=ExportFormat.TSV,
        columns=ColumnSelector(
            basic=["name", "cas_number", "smiles"],
            receptor_data=["binding_profiles"],
            activity_data=["effect_profiles"],
            safety_data=["safety_alerts"],
            predictions=["binding_predictions", "activity_predictions"],
            web_data=["experience_reports", "social_mentions"],
            patent_data=["patent_citations"],
        ),
        filters=FilterConfig(
            psychoactive_class=[PsychoactiveClass.STIMULANT],
            min_binding_affinity=0.7,
            min_risk_level=RiskLevel.MODERATE,
        ),
        batch_size=1000,
        include_metadata=True,
        output_dir=Path("exports"),
        compression=None,
        validation_enabled=True,
    )


@pytest.fixture
def data_exporter(export_config):
    """Create test data exporter fixture."""
    return DataExporter(config=export_config)


class TestBasicExport:
    """Tests for basic export functionality."""

    def test_initialization(self, data_exporter, export_config):
        """Test initialization of DataExporter."""
        assert isinstance(data_exporter.tsv_exporter, TSVExporter)
        assert isinstance(data_exporter.json_exporter, JSONExporter)
        assert isinstance(data_exporter.batch_exporter, BatchExporter)
        assert data_exporter.config == export_config
        assert data_exporter.stats == {}

    def test_basic_tsv_export(self, test_compounds, data_exporter, tmp_path):
        """Test basic TSV export."""
        # Export compounds
        output_file = tmp_path / "compounds.tsv"
        result = data_exporter.export_tsv(
            compounds=test_compounds,
            output_file=output_file,
            columns=["name", "cas_number", "smiles"]
        )
        
        # Check result
        assert isinstance(result, ExportResult)
        assert result.success
        assert output_file.exists()
        
        # Check content
        df = pd.read_csv(output_file, sep="\t")
        assert len(df) == len(test_compounds)
        assert list(df.columns) == ["name", "cas_number", "smiles"]
        assert df.iloc[0]["name"] == "Caffeine"
        assert df.iloc[0]["cas_number"] == "58-08-2"

    def test_basic_json_export(self, test_compounds, data_exporter, tmp_path):
        """Test basic JSON export."""
        # Export compounds
        output_file = tmp_path / "compounds.json"
        result = data_exporter.export_json(
            compounds=test_compounds,
            output_file=output_file
        )
        
        # Check result
        assert result.success
        assert output_file.exists()
        
        # Check content
        with open(output_file) as f:
            data = json.load(f)
            assert len(data["compounds"]) == len(test_compounds)
            assert data["compounds"][0]["name"] == "Caffeine"
            assert data["compounds"][0]["cas_number"] == "58-08-2"


class TestColumnSelection:
    """Tests for column selection functionality."""

    def test_basic_columns(self, test_compounds, data_exporter, tmp_path):
        """Test basic column selection."""
        # Configure columns
        config = ExportConfig(
            format=ExportFormat.TSV,
            columns=ColumnSelector(
                basic=["name", "cas_number"],
                receptor_data=[],
                activity_data=[],
                safety_data=[],
            ),
        )
        
        # Export with selected columns
        output_file = tmp_path / "basic_columns.tsv"
        result = data_exporter.export_tsv(
            compounds=test_compounds,
            output_file=output_file,
            config=config
        )
        
        # Check result
        assert result.success
        df = pd.read_csv(output_file, sep="\t")
        assert list(df.columns) == ["name", "cas_number"]

    def test_all_data_columns(self, test_compounds, data_exporter, tmp_path):
        """Test selection of all data columns."""
        # Configure columns
        config = ExportConfig(
            format=ExportFormat.TSV,
            columns=ColumnSelector(
                basic=["name", "cas_number"],
                receptor_data=["binding_profiles"],
                activity_data=["effect_profiles"],
                safety_data=["safety_alerts"],
                predictions=["binding_predictions"],
                web_data=["experience_reports"],
                patent_data=["patent_citations"],
            ),
        )
        
        # Export with all columns
        output_file = tmp_path / "all_columns.tsv"
        result = data_exporter.export_tsv(
            compounds=test_compounds,
            output_file=output_file,
            config=config
        )
        
        # Check result
        assert result.success
        df = pd.read_csv(output_file, sep="\t")
        assert "binding_profiles" in df.columns
        assert "effect_profiles" in df.columns
        assert "safety_alerts" in df.columns
        assert "binding_predictions" in df.columns
        assert "experience_reports" in df.columns
        assert "patent_citations" in df.columns


class TestFiltering:
    """Tests for filtering functionality."""

    def test_class_filter(self, test_compounds, data_exporter, tmp_path):
        """Test filtering by psychoactive class."""
        # Configure filter
        config = ExportConfig(
            format=ExportFormat.TSV,
            columns=ColumnSelector(basic=["name"]),
            filters=FilterConfig(
                psychoactive_class=[PsychoactiveClass.STIMULANT],
            ),
        )
        
        # Export with filter
        output_file = tmp_path / "class_filtered.tsv"
        result = data_exporter.export_tsv(
            compounds=test_compounds,
            output_file=output_file,
            config=config
        )
        
        # Check result
        assert result.success
        df = pd.read_csv(output_file, sep="\t")
        assert len(df) == 2  # Both compounds are stimulants

    def test_binding_filter(self, test_compounds, data_exporter, tmp_path):
        """Test filtering by binding affinity."""
        # Configure filter
        config = ExportConfig(
            format=ExportFormat.TSV,
            columns=ColumnSelector(basic=["name"]),
            filters=FilterConfig(
                min_binding_affinity=0.7,
            ),
        )
        
        # Export with filter
        output_file = tmp_path / "binding_filtered.tsv"
        result = data_exporter.export_tsv(
            compounds=test_compounds,
            output_file=output_file,
            config=config
        )
        
        # Check result
        assert result.success
        df = pd.read_csv(output_file, sep="\t")
        assert len(df) == 1  # Only caffeine has binding affinity >= 0.7

    def test_risk_filter(self, test_compounds, data_exporter, tmp_path):
        """Test filtering by risk level."""
        # Configure filter
        config = ExportConfig(
            format=ExportFormat.TSV,
            columns=ColumnSelector(basic=["name"]),
            filters=FilterConfig(
                min_risk_level=RiskLevel.HIGH,
            ),
        )
        
        # Export with filter
        output_file = tmp_path / "risk_filtered.tsv"
        result = data_exporter.export_tsv(
            compounds=test_compounds,
            output_file=output_file,
            config=config
        )
        
        # Check result
        assert result.success
        df = pd.read_csv(output_file, sep="\t")
        assert len(df) == 2  # Both compounds have high risks


class TestBatchProcessing:
    """Tests for batch processing functionality."""

    def test_batch_export(self, test_compounds, data_exporter, tmp_path):
        """Test batch export functionality."""
        # Create large compound list
        compounds = test_compounds * 50  # 100 compounds
        
        # Configure batch export
        config = ExportConfig(
            format=ExportFormat.TSV,
            columns=ColumnSelector(basic=["name"]),
            batch_size=10,
        )
        
        # Export in batches
        output_dir = tmp_path / "batch_export"
        output_dir.mkdir()
        result = data_exporter.export_batch(
            compounds=compounds,
            output_dir=output_dir,
            config=config
        )
        
        # Check result
        assert result.success
        assert len(list(output_dir.glob("*.tsv"))) == 10
        assert result.stats["total_batches"] == 10
        assert result.stats["compounds_exported"] == len(compounds)

    def test_parallel_export(self, test_compounds, data_exporter, tmp_path):
        """Test parallel batch export."""
        # Create large compound list
        compounds = test_compounds * 50  # 100 compounds
        
        # Configure parallel export
        config = ExportConfig(
            format=ExportFormat.TSV,
            columns=ColumnSelector(basic=["name"]),
            batch_size=10,
            parallel=True,
            num_workers=4,
        )
        
        # Export in parallel
        output_dir = tmp_path / "parallel_export"
        output_dir.mkdir()
        result = data_exporter.export_batch(
            compounds=compounds,
            output_dir=output_dir,
            config=config
        )
        
        # Check result
        assert result.success
        assert len(list(output_dir.glob("*.tsv"))) == 10
        assert result.stats["total_batches"] == 10
        assert result.stats["compounds_exported"] == len(compounds)
        assert result.stats["num_workers"] == 4


class TestValidation:
    """Tests for validation functionality."""

    def test_column_validation(self, test_compounds, data_exporter, tmp_path):
        """Test column validation."""
        # Try invalid columns
        with pytest.raises(ValueError):
            data_exporter.export_tsv(
                compounds=test_compounds,
                output_file=tmp_path / "invalid.tsv",
                columns=["invalid_column"]
            )

    def test_data_validation(self, test_compounds, data_exporter, tmp_path):
        """Test data validation."""
        # Create invalid compound
        invalid_compound = PsychoactiveCompound(
            name="Invalid",
            smiles="",  # Invalid SMILES
            cas_number="invalid-cas"  # Invalid CAS
        )
        compounds = test_compounds + [invalid_compound]
        
        # Configure validation
        config = ExportConfig(
            format=ExportFormat.TSV,
            columns=ColumnSelector(basic=["name", "smiles", "cas_number"]),
            validation_enabled=True,
        )
        
        # Export with validation
        result = data_exporter.export_tsv(
            compounds=compounds,
            output_file=tmp_path / "validated.tsv",
            config=config
        )
        
        # Check validation results
        assert not result.success
        assert len(result.validation_errors) == 2
        assert "SMILES" in str(result.validation_errors[0])
        assert "CAS" in str(result.validation_errors[1])

    def test_format_validation(self, test_compounds, data_exporter, tmp_path):
        """Test format validation."""
        # Try invalid format
        with pytest.raises(ValueError):
            data_exporter.export_batch(
                compounds=test_compounds,
                output_dir=tmp_path,
                format="invalid_format"
            )


class TestPerformance:
    """Tests for performance functionality."""

    def test_large_export(self, test_compounds, data_exporter, tmp_path):
        """Test large export performance."""
        # Create large compound list
        compounds = test_compounds * 1000  # 2000 compounds
        
        # Configure performance monitoring
        config = ExportConfig(
            format=ExportFormat.TSV,
            columns=ColumnSelector(basic=["name", "cas_number", "smiles"]),
            monitor_performance=True,
        )
        
        # Export with monitoring
        output_file = tmp_path / "large_export.tsv"
        result = data_exporter.export_tsv(
            compounds=compounds,
            output_file=output_file,
            config=config
        )
        
        # Check performance metrics
        assert result.success
        assert "export_time" in result.metrics
        assert "memory_usage" in result.metrics
        assert result.metrics["export_time"] > 0
        assert result.metrics["memory_usage"] > 0

    def test_compression(self, test_compounds, data_exporter, tmp_path):
        """Test compressed export."""
        # Configure compression
        config = ExportConfig(
            format=ExportFormat.TSV,
            columns=ColumnSelector(basic=["name", "cas_number", "smiles"]),
            compression="gzip",
        )
        
        # Export with compression
        output_file = tmp_path / "compounds.tsv.gz"
        result = data_exporter.export_tsv(
            compounds=test_compounds,
            output_file=output_file,
            config=config
        )
        
        # Check result
        assert result.success
        assert output_file.exists()
        
        # Check compressed size
        uncompressed_size = len(pd.DataFrame([
            {
                "name": c.name,
                "cas_number": c.cas_number,
                "smiles": c.smiles
            }
            for c in test_compounds
        ]).to_csv(sep="\t").encode())
        compressed_size = output_file.stat().st_size
        assert compressed_size < uncompressed_size


if __name__ == "__main__":
    pytest.main([__file__])
