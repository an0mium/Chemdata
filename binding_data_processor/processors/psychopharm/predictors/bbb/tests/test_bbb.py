"""Tests for BBB permeability prediction package."""

import pytest
from unittest.mock import Mock, patch

from ......models.core import CompoundData
from ......models.psychopharm import BBBPermeability
from .. import (
    BBBPredictorBase,
    BBBPredictor,
    BBBPredictorEnhanced,
    BBBPredictorWebEnriched,
)


@pytest.fixture
def test_compounds():
    """Create test compounds fixture."""
    compounds = []

    # Caffeine (BBB permeable)
    caffeine = CompoundData(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    compounds.append(caffeine)

    # L-DOPA (BBB permeable)
    ldopa = CompoundData(
        name="L-DOPA",
        smiles="NC(Cc1ccc(O)c(O)c1)C(=O)O",
        cas_number="59-92-7",
    )
    compounds.append(ldopa)

    # 5-HTP (BBB permeable)
    fivehtp = CompoundData(
        name="5-HTP",
        smiles="NC(Cc1c[nH]c2ccc(O)cc12)C(=O)O",
        cas_number="56-69-9",
    )
    compounds.append(fivehtp)

    return compounds


@pytest.fixture
def mock_web_clients():
    """Create mock web clients fixture."""
    return {
        "chembl": Mock(),
        "pubchem": Mock(),
        "swiss": Mock(),
        "community": Mock(),
        "social": Mock(),
        "web_search": Mock(),
    }


@pytest.fixture
def mock_llm_processor():
    """Create mock LLM processor fixture."""
    return Mock()


class TestBBBPredictorBase:
    """Tests for BBBPredictorBase class."""

    def test_initialization(self):
        """Test initialization of BBBPredictorBase."""
        predictor = BBBPredictorBase()
        assert predictor.model_dir is None
        assert predictor.cache_dir is None
        assert predictor.feature_types == ["fingerprints", "descriptors", "enhanced"]
        assert predictor.stats == {}

    def test_predict(self, test_compounds):
        """Test basic prediction functionality."""
        predictor = BBBPredictorBase()
        
        # Test caffeine (small molecule)
        result = predictor.predict(test_compounds[0])
        assert result.value in BBBPermeability
        assert 0 <= result.confidence <= 1
        assert "permeability_score" in result.supporting_data
        
        # Test L-DOPA (amino acid)
        result = predictor.predict(test_compounds[1])
        assert result.value in BBBPermeability
        assert 0 <= result.confidence <= 1
        assert "permeability_score" in result.supporting_data
        
        # Test 5-HTP (amino acid)
        result = predictor.predict(test_compounds[2])
        assert result.value in BBBPermeability
        assert 0 <= result.confidence <= 1
        assert "permeability_score" in result.supporting_data


class TestBBBPredictor:
    """Tests for BBBPredictor class."""

    def test_initialization(self):
        """Test initialization of BBBPredictor."""
        predictor = BBBPredictor()
        assert predictor.transporters
        assert predictor.receptor_transporters
        assert predictor.model_weights

    def test_predict_with_transporters(self, test_compounds):
        """Test prediction with transporter analysis."""
        predictor = BBBPredictor()
        
        # Test caffeine (passive diffusion)
        result = predictor.predict(test_compounds[0])
        assert result.value in BBBPermeability
        assert 0 <= result.confidence <= 1
        assert "transporters" in result.supporting_data
        assert "receptor_transport" in result.supporting_data
        
        # Test L-DOPA (LAT1 substrate)
        result = predictor.predict(test_compounds[1])
        assert result.value in BBBPermeability
        assert 0 <= result.confidence <= 1
        assert "transporters" in result.supporting_data
        assert "LAT1" in result.supporting_data["transporters"]
        
        # Test 5-HTP (LAT1 substrate)
        result = predictor.predict(test_compounds[2])
        assert result.value in BBBPermeability
        assert 0 <= result.confidence <= 1
        assert "transporters" in result.supporting_data
        assert "LAT1" in result.supporting_data["transporters"]


class TestBBBPredictorEnhanced:
    """Tests for BBBPredictorEnhanced class."""

    def test_initialization(self):
        """Test initialization of BBBPredictorEnhanced."""
        predictor = BBBPredictorEnhanced()
        assert predictor.predictors
        assert "abuse" in predictor.predictors
        assert "toxicity" in predictor.predictors
        assert "receptors" in predictor.predictors

    def test_predict_with_ml(self, test_compounds):
        """Test prediction with ML models."""
        predictor = BBBPredictorEnhanced()
        
        # Test all compounds
        for compound in test_compounds:
            result = predictor.predict(compound)
            assert result.value in BBBPermeability
            assert 0 <= result.confidence <= 1
            assert "abuse" in result.supporting_data
            assert "toxicity" in result.supporting_data
            assert "receptors" in result.supporting_data


class TestBBBPredictorWebEnriched:
    """Tests for BBBPredictorWebEnriched class."""

    @patch("web_enrichment.llm_utils.LLMProcessor")
    def test_initialization(self, mock_llm_class, mock_web_clients):
        """Test initialization of BBBPredictorWebEnriched."""
        predictor = BBBPredictorWebEnriched(web_clients=mock_web_clients)
        assert predictor.web_clients == mock_web_clients
        assert predictor.llm_processor is not None

    @patch("web_enrichment.llm_utils.LLMProcessor")
    def test_predict_with_web_data(
        self, mock_llm_class, test_compounds, mock_web_clients, mock_llm_processor
    ):
        """Test prediction with web data enrichment."""
        # Configure mocks
        mock_llm_class.return_value = mock_llm_processor
        mock_web_clients["chembl"].get_compound_data.return_value = {
            "bbb_data": {"permeability": "high"}
        }
        mock_llm_processor.extract_bbb_data.return_value = {
            "supports_prediction": True,
            "permeability": "high",
        }

        # Create predictor
        predictor = BBBPredictorWebEnriched(web_clients=mock_web_clients)

        # Test all compounds
        for compound in test_compounds:
            result = predictor.predict(compound)
            assert result.value in BBBPermeability
            assert 0 <= result.confidence <= 1
            assert "web_data" in result.supporting_data
            assert "bbb_data" in result.supporting_data["web_data"]

    def test_export_predictions(
        self, test_compounds, mock_web_clients, tmp_path
    ):
        """Test prediction export functionality."""
        # Create predictor
        predictor = BBBPredictorWebEnriched(web_clients=mock_web_clients)

        # Make predictions
        for compound in test_compounds:
            predictor.predict(compound)

        # Export predictions
        output_path = tmp_path / "predictions.tsv"
        predictor.export_predictions(
            output_path,
            include_supporting_data=True,
            include_web_data=True,
        )

        # Verify export
        assert output_path.exists()
        content = output_path.read_text()
        assert "compound_name" in content
        assert "Caffeine" in content
        assert "L-DOPA" in content
        assert "5-HTP" in content
        assert "chembl_data" in content
        assert "pubchem_data" in content


if __name__ == "__main__":
    pytest.main([__file__])
