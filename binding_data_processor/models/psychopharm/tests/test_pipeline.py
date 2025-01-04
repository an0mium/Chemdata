"""Tests for pipeline functionality."""

import pytest
from datetime import datetime
from unittest.mock import patch, MagicMock
import pandas as pd

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..pipeline import (
    Pipeline,
    PsychopharmPipeline,
    ProcessingStage,
    EnrichmentStage,
    AnalysisStage,
    ExportStage,
    DataLoader,
    DataEnricher,
    DataAnalyzer,
    DataExporter,
    PipelineConfig,
    StageResult,
    PipelineResult,
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
def pipeline_config():
    """Create test pipeline configuration fixture."""
    return PipelineConfig(
        processing={"format": "tsv"},
        enrichment={"sources": ["chembl", "pubchem", "community"]},
        analysis={"types": ["binding", "activity", "safety"]},
        export={"format": "tsv", "columns": ["name", "cas_number", "receptor_profiles"]}
    )


@pytest.fixture
def pipeline(pipeline_config):
    """Create test pipeline fixture."""
    return Pipeline(config=pipeline_config)


@pytest.fixture
def psychopharm_pipeline():
    """Create test psychopharm pipeline fixture."""
    return PsychopharmPipeline()


@pytest.fixture
def test_data(tmp_path):
    """Create test input data fixture."""
    # Create test TSV file
    input_file = tmp_path / "input.tsv"
    data = pd.DataFrame({
        "name": ["Caffeine", "Amphetamine"],
        "smiles": [
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "CC(N)CC1=CC=CC=C1"
        ],
        "cas_number": ["58-08-2", "300-62-9"]
    })
    data.to_csv(input_file, sep="\t", index=False)
    return input_file


class TestPipeline:
    """Tests for Pipeline class."""

    def test_initialization(self, pipeline, pipeline_config):
        """Test initialization of Pipeline."""
        assert isinstance(pipeline.processing_stage, ProcessingStage)
        assert isinstance(pipeline.enrichment_stage, EnrichmentStage)
        assert isinstance(pipeline.analysis_stage, AnalysisStage)
        assert isinstance(pipeline.export_stage, ExportStage)
        assert isinstance(pipeline.loader, DataLoader)
        assert isinstance(pipeline.enricher, DataEnricher)
        assert isinstance(pipeline.analyzer, DataAnalyzer)
        assert isinstance(pipeline.exporter, DataExporter)
        assert pipeline.config == pipeline_config
        assert pipeline.stats == {}

    def test_run_pipeline(self, pipeline, test_data, tmp_path):
        """Test running complete pipeline."""
        # Configure pipeline
        config = PipelineConfig(
            input_file=test_data,
            output_file=tmp_path / "output.tsv",
            enrich_sources=["chembl", "pubchem", "community"],
            analyses=["binding", "activity", "safety"],
            export_format="tsv"
        )
        
        # Run pipeline
        result = pipeline.run(config)
        
        # Check result
        assert isinstance(result, PipelineResult)
        assert result.success
        assert len(result.compounds) == 2
        assert result.stats["processed_compounds"] == 2
        
        # Check output file
        output_file = tmp_path / "output.tsv"
        assert output_file.exists()
        df = pd.read_csv(output_file, sep="\t")
        assert len(df) == 2
        assert "name" in df.columns
        assert "cas_number" in df.columns
        assert "receptor_profiles" in df.columns
        assert "safety_alerts" in df.columns

    def test_pipeline_stages(self, pipeline, test_compounds):
        """Test individual pipeline stages."""
        # Test processing stage
        result = pipeline.run_processing(compounds=test_compounds)
        assert isinstance(result, StageResult)
        assert result.success
        assert len(result.processed_compounds) == 2
        assert result.stats["total_processed"] == 2
        
        # Test enrichment stage
        with patch("requests.get") as mock_get:
            mock_get.return_value.json.return_value = {
                "receptor_data": {"A2A": {"affinity": 0.8}},
                "safety_data": {"alerts": ["anxiety", "insomnia"]},
                "effect_data": {"stimulation": 0.8, "focus": 0.7}
            }
            result = pipeline.run_enrichment(compounds=test_compounds)
            assert result.success
            assert len(result.processed_compounds) == 2
            assert result.stats["total_enriched"] == 2
        
        # Test analysis stage
        result = pipeline.run_analysis(compounds=test_compounds)
        assert result.success
        assert "binding_analysis" in result.analysis_results
        assert "activity_analysis" in result.analysis_results
        assert "safety_analysis" in result.analysis_results
        
        # Test export stage
        result = pipeline.run_export(
            compounds=test_compounds,
            output_file=pipeline.config.output_file
        )
        assert result.success
        assert result.stats["total_exported"] == 2

    def test_pipeline_checkpointing(self, pipeline, test_data, tmp_path):
        """Test pipeline checkpointing functionality."""
        # Configure pipeline with checkpointing
        config = PipelineConfig(
            input_file=test_data,
            output_file=tmp_path / "output.tsv",
            checkpoint_dir=tmp_path / "checkpoints",
            checkpoint_interval=1
        )
        
        # Run pipeline
        result = pipeline.run(config)
        
        # Check checkpoints
        checkpoint_dir = tmp_path / "checkpoints"
        assert checkpoint_dir.exists()
        checkpoints = list(checkpoint_dir.glob("*.pkl"))
        assert len(checkpoints) > 0
        
        # Test checkpoint recovery
        pipeline = Pipeline()
        result = pipeline.recover_from_checkpoint(
            checkpoint_file=checkpoints[-1],
            config=config
        )
        assert result.success
        assert len(result.compounds) == 2

    def test_pipeline_error_handling(self, pipeline, test_data, tmp_path):
        """Test pipeline error handling."""
        # Configure pipeline with invalid settings
        config = PipelineConfig(
            input_file=test_data,
            output_file=tmp_path / "output.tsv",
            enrich_sources=["invalid_source"]
        )
        
        # Run pipeline
        result = pipeline.run(config)
        
        # Check error handling
        assert not result.success
        assert "invalid enrichment source" in str(result.error).lower()
        assert result.stats["failed_compounds"] > 0

    def test_parallel_processing(self, pipeline, test_compounds):
        """Test parallel processing functionality."""
        # Enable parallel processing
        pipeline.config.parallel = True
        pipeline.config.n_jobs = 2
        
        # Run processing stage
        result = pipeline.run_processing(compounds=test_compounds)
        
        # Check parallel execution
        assert result.success
        assert result.stats["n_jobs"] == 2
        assert "parallel_time" in result.stats

    def test_monitoring(self, pipeline, test_compounds):
        """Test pipeline monitoring."""
        # Add progress callback
        progress = []
        pipeline.config.progress_callback = lambda x: progress.append(x)
        
        # Run pipeline
        pipeline.run_processing(compounds=test_compounds)
        
        # Check progress monitoring
        assert len(progress) > 0
        assert progress[-1] == 100
        assert "processing_time" in pipeline.stats
        assert "memory_usage" in pipeline.stats


class TestPsychopharmPipeline:
    """Tests for PsychopharmPipeline class."""

    def test_initialization(self, psychopharm_pipeline):
        """Test initialization of PsychopharmPipeline."""
        assert psychopharm_pipeline.compounds == {}
        assert psychopharm_pipeline.last_run is None
        assert psychopharm_pipeline.stats == {}

    def test_add_compound(self, psychopharm_pipeline, test_compounds):
        """Test add_compound method."""
        # Add test compounds
        for compound in test_compounds:
            psychopharm_pipeline.add_compound(compound)
        
        # Check compounds were added
        assert len(psychopharm_pipeline.compounds) == 2
        assert "58-08-2" in psychopharm_pipeline.compounds  # Caffeine
        assert "300-62-9" in psychopharm_pipeline.compounds  # Amphetamine
        
        # Check compound data preserved
        caffeine = psychopharm_pipeline.compounds["58-08-2"]
        assert caffeine.name == "Caffeine"
        assert caffeine.psychoactive_class == PsychoactiveClass.STIMULANT
        assert "A2A" in caffeine.receptor_profiles
        assert "stimulation" in caffeine.effect_profile

    def test_process_compounds(self, psychopharm_pipeline, test_compounds):
        """Test process_compounds method."""
        # Add test compounds
        for compound in test_compounds:
            psychopharm_pipeline.add_compound(compound)
        
        # Process compounds
        psychopharm_pipeline.process_compounds()
        
        # Check processing stats
        assert psychopharm_pipeline.stats["total_compounds"] == 2
        assert psychopharm_pipeline.stats["processed_compounds"] == 2
        assert isinstance(psychopharm_pipeline.last_run, datetime)
        
        # Check compound enrichment
        caffeine = psychopharm_pipeline.compounds["58-08-2"]
        assert len(caffeine.receptor_profiles) >= 1
        assert caffeine.psychoactive_class == PsychoactiveClass.STIMULANT
        assert len(caffeine.safety_alerts) >= 2
        assert len(caffeine.effect_profile) >= 2
        
        # Check predictions
        assert "binding_predictions" in caffeine.predictions
        assert "activity_predictions" in caffeine.predictions
        assert "safety_predictions" in caffeine.predictions

    def test_analyze_compound_relationships(self, psychopharm_pipeline, test_compounds):
        """Test analyze_compound_relationships method."""
        # Add test compounds
        for compound in test_compounds:
            psychopharm_pipeline.add_compound(compound)
        
        # Analyze relationships
        relationships = psychopharm_pipeline.analyze_compound_relationships()
        
        # Check relationship data
        assert "similar_compounds" in relationships
        assert "shared_targets" in relationships
        assert "risk_patterns" in relationships
        assert "effect_patterns" in relationships
        
        # Check specific relationships
        similar = relationships["similar_compounds"]
        assert len(similar) > 0
        assert isinstance(similar[0], tuple)  # (compound1, compound2, similarity)
        
        targets = relationships["shared_targets"]
        assert isinstance(targets, dict)
        assert any(len(compounds) > 1 for compounds in targets.values())
        
        patterns = relationships["risk_patterns"]
        assert isinstance(patterns, list)
        assert len(patterns) > 0
        
        effects = relationships["effect_patterns"]
        assert isinstance(effects, dict)
        assert "stimulation" in effects  # Common effect between compounds


class TestDataLoader:
    """Tests for DataLoader class."""

    def test_load_tsv(self, pipeline, test_data):
        """Test loading TSV data."""
        result = pipeline.loader.load_tsv(test_data)
        assert result.success
        assert len(result.compounds) == 2
        assert all(isinstance(c, PsychoactiveCompound) for c in result.compounds)

    def test_load_json(self, pipeline, tmp_path):
        """Test loading JSON data."""
        # Create test JSON file
        input_file = tmp_path / "input.json"
        data = {
            "compounds": [
                {
                    "name": "Caffeine",
                    "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                    "cas_number": "58-08-2"
                }
            ]
        }
        input_file.write_text(pd.io.json.dumps(data))
        
        # Load data
        result = pipeline.loader.load_json(input_file)
        assert result.success
        assert len(result.compounds) == 1
        assert result.compounds[0].name == "Caffeine"

    def test_invalid_input(self, pipeline, tmp_path):
        """Test handling invalid input."""
        # Test nonexistent file
        result = pipeline.loader.load_tsv(tmp_path / "nonexistent.tsv")
        assert not result.success
        assert "file not found" in str(result.error).lower()
        
        # Test invalid format
        input_file = tmp_path / "invalid.txt"
        input_file.write_text("invalid data")
        result = pipeline.loader.load_tsv(input_file)
        assert not result.success
        assert "invalid format" in str(result.error).lower()


class TestDataEnricher:
    """Tests for DataEnricher class."""

    @patch("requests.get")
    def test_chembl_enrichment(self, mock_get, pipeline, test_compounds):
        """Test ChEMBL data enrichment."""
        # Mock ChEMBL API response
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "compounds": [{
                "molecule_chembl_id": "CHEMBL113",
                "molecule_properties": {
                    "alogp": "0.61",
                    "psa": "58.44"
                }
            }]
        }
        mock_get.return_value = mock_response
        
        # Enrich compounds
        result = pipeline.enricher.enrich_from_chembl(test_compounds[0])
        assert result.success
        assert "chembl_id" in result.data
        assert "properties" in result.data

    @patch("requests.get")
    def test_pubchem_enrichment(self, mock_get, pipeline, test_compounds):
        """Test PubChem data enrichment."""
        # Mock PubChem API response
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "PC_Compounds": [{
                "id": {"id": {"cid": 2519}},
                "props": [{
                    "urn": {"label": "Molecular Weight"},
                    "value": {"sval": "194.19"}
                }]
            }]
        }
        mock_get.return_value = mock_response
        
        # Enrich compounds
        result = pipeline.enricher.enrich_from_pubchem(test_compounds[0])
        assert result.success
        assert "pubchem_id" in result.data
        assert "properties" in result.data

    def test_community_enrichment(self, pipeline, test_compounds):
        """Test community data enrichment."""
        result = pipeline.enricher.enrich_from_community(test_compounds[0])
        assert result.success
        assert "community_data" in result.data
        assert "safety_alerts" in result.data


class TestDataAnalyzer:
    """Tests for DataAnalyzer class."""

    def test_binding_analysis(self, pipeline, test_compounds):
        """Test binding profile analysis."""
        result = pipeline.analyzer.analyze_binding(test_compounds[0])
        assert result.success
        assert "binding_profile" in result.data
        assert "target_selectivity" in result.data

    def test_activity_analysis(self, pipeline, test_compounds):
        """Test activity analysis."""
        result = pipeline.analyzer.analyze_activity(test_compounds[0])
        assert result.success
        assert "activity_profile" in result.data
        assert "mechanism_prediction" in result.data

    def test_safety_analysis(self, pipeline, test_compounds):
        """Test safety analysis."""
        result = pipeline.analyzer.analyze_safety(test_compounds[0])
        assert result.success
        assert "safety_profile" in result.data
        assert "risk_assessment" in result.data


class TestDataExporter:
    """Tests for DataExporter class."""

    def test_export_tsv(self, pipeline, test_compounds, tmp_path):
        """Test TSV export."""
        output_file = tmp_path / "export.tsv"
        result = pipeline.exporter.export_tsv(
            compounds=test_compounds,
            output_file=output_file,
            columns=["name", "cas_number", "receptor_profiles"]
        )
        assert result.success
        assert output_file.exists()
        df = pd.read_csv(output_file, sep="\t")
        assert len(df) == 2
        assert all(col in df.columns for col in ["name", "cas_number"])

    def test_export_json(self, pipeline, test_compounds, tmp_path):
        """Test JSON export."""
        output_file = tmp_path / "export.json"
        result = pipeline.exporter.export_json(
            compounds=test_compounds,
            output_file=output_file
        )
        assert result.success
        assert output_file.exists()
        data = pd.read_json(output_file)
        assert len(data) == 2
        assert "compounds" in data.columns

    def test_invalid_export(self, pipeline, test_compounds, tmp_path):
        """Test invalid export handling."""
        # Test invalid format
        result = pipeline.exporter.export(
            compounds=test_compounds,
            output_file=tmp_path / "export.txt",
            format="invalid"
        )
        assert not result.success
        assert "invalid format" in str(result.error).lower()


if __name__ == "__main__":
    pytest.main([__file__])
