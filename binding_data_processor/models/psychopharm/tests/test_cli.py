"""Tests for command-line interface functionality."""

import pytest
from unittest.mock import patch
from click.testing import CliRunner

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..cli import (
    CLI,
    CompoundCommand,
    SearchCommand,
    FilterCommand,
    ExportCommand,
    process_command,
    enrich_command,
    analyze_command,
    serve_command,
    validate_command,
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
    caffeine.effect_profile = {
        "stimulation": (0.8, 0.9),
        "focus": (0.7, 0.8),
    }
    caffeine.safety_alerts = {
        "anxiety": RiskLevel.MODERATE,
        "insomnia": RiskLevel.HIGH,
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
    amphetamine.effect_profile = {
        "stimulation": (0.9, 0.95),
        "euphoria": (0.8, 0.9),
    }
    amphetamine.safety_alerts = {
        "addiction": RiskLevel.HIGH,
        "cardiovascular": RiskLevel.HIGH,
    }
    compounds.append(amphetamine)
    
    return compounds


@pytest.fixture
def cli(test_compounds):
    """Create test CLI fixture."""
    cli = CLI()
    cli.load_compounds(test_compounds)
    return cli


@pytest.fixture
def runner():
    """Create test CLI runner fixture."""
    return CliRunner()


@pytest.fixture
def test_input(tmp_path, test_compounds):
    """Create test input file fixture."""
    input_file = tmp_path / "input.tsv"
    with input_file.open("w") as f:
        f.write("name\tcas_number\tsmiles\n")
        for compound in test_compounds:
            f.write(f"{compound.name}\t{compound.cas_number}\t{compound.smiles}\n")
    return input_file


class TestCLI:
    """Tests for CLI class."""

    def test_initialization(self, cli):
        """Test initialization of CLI."""
        assert isinstance(cli.compound_command, CompoundCommand)
        assert isinstance(cli.search_command, SearchCommand)
        assert isinstance(cli.filter_command, FilterCommand)
        assert isinstance(cli.export_command, ExportCommand)
        assert cli.stats == {}

    def test_version_command(self, runner, cli):
        """Test version command."""
        result = runner.invoke(cli.app, ["--version"])
        assert result.exit_code == 0
        assert cli.version in result.output

    def test_help_command(self, runner, cli):
        """Test help command."""
        result = runner.invoke(cli.app, ["--help"])
        assert result.exit_code == 0
        assert "Commands:" in result.output
        assert "compounds" in result.output
        assert "search" in result.output
        assert "filter" in result.output
        assert "export" in result.output
        assert "process" in result.output
        assert "enrich" in result.output
        assert "analyze" in result.output
        assert "serve" in result.output


class TestCompoundCommand:
    """Tests for CompoundCommand class."""

    def test_list_compounds(self, runner, cli):
        """Test listing compounds."""
        result = runner.invoke(cli.app, ["compounds", "list"])
        assert result.exit_code == 0
        assert "Caffeine" in result.output
        assert "58-08-2" in result.output
        assert "Amphetamine" in result.output
        assert "300-62-9" in result.output

    def test_show_compound(self, runner, cli):
        """Test showing compound details."""
        result = runner.invoke(cli.app, ["compounds", "show", "58-08-2"])  # Caffeine CAS
        assert result.exit_code == 0
        assert "Caffeine" in result.output
        assert "58-08-2" in result.output
        assert "STIMULANT" in result.output
        assert "A2A" in result.output
        assert "antagonist" in result.output

    def test_show_nonexistent_compound(self, runner, cli):
        """Test showing nonexistent compound."""
        result = runner.invoke(cli.app, ["compounds", "show", "invalid-cas"])
        assert result.exit_code == 1
        assert "not found" in result.output.lower()


class TestSearchCommand:
    """Tests for SearchCommand class."""

    def test_text_search(self, runner, cli):
        """Test text search."""
        # Search by name
        result = runner.invoke(cli.app, ["search", "text", "caffeine", "--type", "name"])
        assert result.exit_code == 0
        assert "Caffeine" in result.output
        assert "58-08-2" in result.output
        assert "Amphetamine" not in result.output
        
        # Search by receptor
        result = runner.invoke(cli.app, ["search", "text", "DAT", "--type", "receptor"])
        assert result.exit_code == 0
        assert "Amphetamine" in result.output
        assert "300-62-9" in result.output
        assert "Caffeine" not in result.output

    def test_structure_search(self, runner, cli):
        """Test structure search."""
        result = runner.invoke(cli.app, [
            "search", "structure",
            "--smiles", "CN1C=NC2=C1C(=O)N",
            "--threshold", "0.8"
        ])
        assert result.exit_code == 0
        assert "Caffeine" in result.output
        assert "58-08-2" in result.output
        assert "Amphetamine" not in result.output

    def test_invalid_search(self, runner, cli):
        """Test invalid search."""
        result = runner.invoke(cli.app, ["search", "text", "--type", "invalid"])
        assert result.exit_code == 2
        assert "invalid choice" in result.output.lower()


class TestFilterCommand:
    """Tests for FilterCommand class."""

    def test_property_filter(self, runner, cli):
        """Test property filter."""
        result = runner.invoke(cli.app, [
            "filter", "property",
            "--name", "binding_affinity",
            "--min", "0.7",
            "--max", "1.0"
        ])
        assert result.exit_code == 0
        assert "Caffeine" in result.output
        assert "58-08-2" in result.output
        assert "Amphetamine" not in result.output

    def test_class_filter(self, runner, cli):
        """Test class filter."""
        # Filter by class
        result = runner.invoke(cli.app, [
            "filter", "class",
            "--class", "STIMULANT"
        ])
        assert result.exit_code == 0
        assert "Caffeine" in result.output
        assert "Amphetamine" in result.output
        
        # Filter by risk level
        result = runner.invoke(cli.app, [
            "filter", "risk",
            "--min-level", "HIGH"
        ])
        assert result.exit_code == 0
        assert "Caffeine" in result.output  # Has high insomnia risk
        assert "Amphetamine" in result.output  # Has high addiction risk

    def test_invalid_filter(self, runner, cli):
        """Test invalid filter."""
        result = runner.invoke(cli.app, [
            "filter", "property",
            "--name", "invalid_property"
        ])
        assert result.exit_code == 2
        assert "invalid property" in result.output.lower()


class TestExportCommand:
    """Tests for ExportCommand class."""

    def test_export_tsv(self, runner, cli, tmp_path):
        """Test TSV export."""
        output_file = tmp_path / "compounds.tsv"
        result = runner.invoke(cli.app, [
            "export",
            "--format", "tsv",
            "--columns", "name,cas_number,smiles",
            "--output", str(output_file)
        ])
        assert result.exit_code == 0
        assert output_file.exists()
        content = output_file.read_text()
        assert "name\tcas_number\tsmiles" in content
        assert "Caffeine\t58-08-2" in content

    def test_export_json(self, runner, cli, tmp_path):
        """Test JSON export."""
        output_file = tmp_path / "compounds.json"
        result = runner.invoke(cli.app, [
            "export",
            "--format", "json",
            "--columns", "name,cas_number",
            "--output", str(output_file)
        ])
        assert result.exit_code == 0
        assert output_file.exists()
        content = output_file.read_text()
        assert "Caffeine" in content
        assert "58-08-2" in content

    def test_invalid_export(self, runner, cli):
        """Test invalid export."""
        result = runner.invoke(cli.app, [
            "export",
            "--format", "invalid_format"
        ])
        assert result.exit_code == 2
        assert "invalid choice" in result.output.lower()


class TestPipelineCommands:
    """Tests for pipeline commands."""

    def test_process_command(self, runner, test_input, tmp_path):
        """Test process command."""
        output_file = tmp_path / "output.tsv"
        result = runner.invoke(process_command, [
            "--input", str(test_input),
            "--output", str(output_file),
            "--batch-size", "10",
            "--parallel",
            "--n-jobs", "2"
        ])
        assert result.exit_code == 0
        assert output_file.exists()
        assert "Processing complete" in result.output

    def test_enrich_command(self, runner, test_input, tmp_path):
        """Test enrich command."""
        output_file = tmp_path / "enriched.tsv"
        result = runner.invoke(enrich_command, [
            "--input", str(test_input),
            "--output", str(output_file),
            "--sources", "chembl,pubchem",
            "--community",
            "--social"
        ])
        assert result.exit_code == 0
        assert output_file.exists()
        assert "Enrichment complete" in result.output

    def test_analyze_command(self, runner, test_input, tmp_path):
        """Test analyze command."""
        output_file = tmp_path / "analysis.tsv"
        result = runner.invoke(analyze_command, [
            "--input", str(test_input),
            "--output", str(output_file),
            "--binding",
            "--activity",
            "--safety",
            "--sar"
        ])
        assert result.exit_code == 0
        assert output_file.exists()
        assert "Analysis complete" in result.output

    def test_validate_command(self, runner, test_input):
        """Test validate command."""
        result = runner.invoke(validate_command, [
            "--input", str(test_input),
            "--check-structure",
            "--check-properties",
            "--check-data"
        ])
        assert result.exit_code == 0
        assert "Validation complete" in result.output

    @patch("uvicorn.run")
    def test_serve_command(self, mock_run, runner):
        """Test serve command."""
        result = runner.invoke(serve_command, [
            "--host", "localhost",
            "--port", "8000",
            "--reload"
        ])
        assert result.exit_code == 0
        mock_run.assert_called_once()
        args = mock_run.call_args[1]
        assert args["host"] == "localhost"
        assert args["port"] == 8000
        assert args["reload"] is True


class TestProgressMonitoring:
    """Tests for progress monitoring."""

    def test_verbose_output(self, runner, test_input, tmp_path):
        """Test verbose output mode."""
        result = runner.invoke(process_command, [
            "--input", str(test_input),
            "--output", str(tmp_path / "output.tsv"),
            "--verbose"
        ])
        assert result.exit_code == 0
        assert "DEBUG" in result.output
        assert "Processing compound:" in result.output
        assert "Caffeine" in result.output

    def test_progress_display(self, runner, test_input, tmp_path):
        """Test progress display."""
        result = runner.invoke(process_command, [
            "--input", str(test_input),
            "--output", str(tmp_path / "output.tsv"),
            "--show-progress"
        ])
        assert result.exit_code == 0
        assert "Progress:" in result.output
        assert "100%" in result.output


class TestErrorHandling:
    """Tests for error handling."""

    def test_missing_required_args(self, runner):
        """Test missing required arguments."""
        result = runner.invoke(process_command)
        assert result.exit_code != 0
        assert "Missing required argument" in result.output

    def test_invalid_input_file(self, runner, tmp_path):
        """Test invalid input file."""
        result = runner.invoke(process_command, [
            "--input", str(tmp_path / "nonexistent.tsv"),
            "--output", str(tmp_path / "output.tsv")
        ])
        assert result.exit_code != 0
        assert "Input file not found" in result.output

    def test_invalid_output_directory(self, runner, test_input):
        """Test invalid output directory."""
        result = runner.invoke(process_command, [
            "--input", str(test_input),
            "--output", "/invalid/path/output.tsv"
        ])
        assert result.exit_code != 0
        assert "Output directory not found" in result.output

    def test_invalid_batch_size(self, runner, test_input, tmp_path):
        """Test invalid batch size."""
        result = runner.invoke(process_command, [
            "--input", str(test_input),
            "--output", str(tmp_path / "output.tsv"),
            "--batch-size", "0"
        ])
        assert result.exit_code != 0
        assert "Invalid batch size" in result.output

    def test_invalid_n_jobs(self, runner, test_input, tmp_path):
        """Test invalid number of jobs."""
        result = runner.invoke(process_command, [
            "--input", str(test_input),
            "--output", str(tmp_path / "output.tsv"),
            "--parallel",
            "--n-jobs", "0"
        ])
        assert result.exit_code != 0
        assert "Invalid number of jobs" in result.output


if __name__ == "__main__":
    pytest.main([__file__])
