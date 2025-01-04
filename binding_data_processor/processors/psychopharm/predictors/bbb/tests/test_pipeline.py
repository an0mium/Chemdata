"""Integration tests for BBB permeability prediction pipeline.

This script tests the entire BBB prediction pipeline by:
1. Loading example compounds
2. Running predictions with each predictor level
3. Validating web data enrichment
4. Verifying export functionality
"""

import logging
import tempfile
from pathlib import Path

import pytest

from ......models.core import CompoundData
from .. import (
    BBBPredictorBase,
    BBBPredictor,
    BBBPredictorEnhanced,
    BBBPredictorWebEnriched,
)


def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )


@pytest.fixture
def example_compounds():
    """Load example compounds from test data."""
    compounds = []
    example_file = (
        Path(__file__)
        .parent.parent.parent.parent.parent.parent.parent
        / "examples/data/example_compounds.tsv"
    )
    
    with open(example_file) as f:
        next(f)  # Skip header
        for line in f:
            if line.strip() and not line.startswith("#"):
                name, smiles, cas = line.strip().split("\t")
                compound = CompoundData(
                    name=name,
                    smiles=smiles,
                    cas_number=cas,
                )
                compounds.append(compound)
    
    return compounds


@pytest.fixture
def temp_dir():
    """Create temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


def test_base_predictor(example_compounds, temp_dir):
    """Test base BBB predictor functionality."""
    predictor = BBBPredictorBase(
        model_dir=str(temp_dir / "models"),
        cache_dir=str(temp_dir / "cache"),
    )
    
    for compound in example_compounds:
        result = predictor.predict(compound)
        
        # Validate result structure
        assert result.value is not None
        assert 0 <= result.confidence <= 1
        assert "permeability_score" in result.supporting_data
        
        # Log prediction
        logging.info(
            f"{compound.name}: {result.value} "
            f"(confidence: {result.confidence:.2f})"
        )


def test_transporter_predictor(example_compounds, temp_dir):
    """Test BBB predictor with transporter analysis."""
    predictor = BBBPredictor(
        model_dir=str(temp_dir / "models"),
        cache_dir=str(temp_dir / "cache"),
    )
    
    for compound in example_compounds:
        result = predictor.predict(compound)
        
        # Validate result structure
        assert result.value is not None
        assert 0 <= result.confidence <= 1
        assert "transporters" in result.supporting_data
        assert "receptor_transport" in result.supporting_data
        
        # Log prediction
        logging.info(
            f"{compound.name}: {result.value} "
            f"(confidence: {result.confidence:.2f})"
        )
        logging.debug(f"Transporters: {result.supporting_data['transporters']}")


def test_enhanced_predictor(example_compounds, temp_dir):
    """Test enhanced BBB predictor with ML models."""
    predictor = BBBPredictorEnhanced(
        model_dir=str(temp_dir / "models"),
        cache_dir=str(temp_dir / "cache"),
    )
    
    for compound in example_compounds:
        result = predictor.predict(compound)
        
        # Validate result structure
        assert result.value is not None
        assert 0 <= result.confidence <= 1
        assert "abuse" in result.supporting_data
        assert "toxicity" in result.supporting_data
        assert "receptors" in result.supporting_data
        
        # Log prediction
        logging.info(
            f"{compound.name}: {result.value} "
            f"(confidence: {result.confidence:.2f})"
        )
        logging.debug(f"ML predictions: {result.supporting_data}")


def test_web_enriched_predictor(example_compounds, temp_dir):
    """Test web-enriched BBB predictor."""
    predictor = BBBPredictorWebEnriched(
        model_dir=str(temp_dir / "models"),
        cache_dir=str(temp_dir / "cache"),
    )
    
    for compound in example_compounds:
        result = predictor.predict(compound)
        
        # Validate result structure
        assert result.value is not None
        assert 0 <= result.confidence <= 1
        assert "web_data" in result.supporting_data
        assert "bbb_data" in result.supporting_data["web_data"]
        
        # Log prediction
        logging.info(
            f"{compound.name}: {result.value} "
            f"(confidence: {result.confidence:.2f})"
        )
        logging.debug(f"Web data: {result.supporting_data['web_data']}")


def test_export_functionality(example_compounds, temp_dir):
    """Test prediction export functionality."""
    predictor = BBBPredictorWebEnriched(
        model_dir=str(temp_dir / "models"),
        cache_dir=str(temp_dir / "cache"),
    )
    
    # Make predictions
    for compound in example_compounds:
        predictor.predict(compound)
    
    # Test different export configurations
    export_configs = [
        {
            "output_file": temp_dir / "basic_export.tsv",
            "include_supporting_data": False,
            "include_web_data": False,
        },
        {
            "output_file": temp_dir / "full_export.tsv",
            "include_supporting_data": True,
            "include_web_data": True,
        },
        {
            "output_file": temp_dir / "custom_export.tsv",
            "columns": [
                "compound_name",
                "permeability_class",
                "confidence",
                "transporter_data",
            ],
            "include_supporting_data": True,
            "include_web_data": False,
        },
    ]
    
    for config in export_configs:
        predictor.export_predictions(**config)
        assert config["output_file"].exists()
        
        # Verify file content
        with open(config["output_file"]) as f:
            header = next(f).strip().split("\t")
            assert "compound_name" in header
            assert len(header) > 1
            
            # Check data rows
            data_rows = list(f)
            assert len(data_rows) == len(example_compounds)


def test_validation_compounds(example_compounds, temp_dir):
    """Test predictions on validation compounds."""
    predictor = BBBPredictorWebEnriched(
        model_dir=str(temp_dir / "models"),
        cache_dir=str(temp_dir / "cache"),
    )
    
    # Split compounds into CNS-active and peripherally selective
    cns_active = [c for c in example_compounds if not c.name.startswith("#")][:10]
    peripheral = [c for c in example_compounds if not c.name.startswith("#")][10:]
    
    # Test CNS-active compounds
    for compound in cns_active:
        result = predictor.predict(compound)
        logging.info(
            f"CNS-active - {compound.name}: {result.value} "
            f"(confidence: {result.confidence:.2f})"
        )
        # Should predict high BBB permeability
        assert result.value in ["HIGH", "MEDIUM-HIGH"]
    
    # Test peripherally selective compounds
    for compound in peripheral:
        result = predictor.predict(compound)
        logging.info(
            f"Peripheral - {compound.name}: {result.value} "
            f"(confidence: {result.confidence:.2f})"
        )
        # Should predict low BBB permeability
        assert result.value in ["LOW", "MEDIUM-LOW"]


if __name__ == "__main__":
    setup_logging()
    pytest.main([__file__, "-v"])
