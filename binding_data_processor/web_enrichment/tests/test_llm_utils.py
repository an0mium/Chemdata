"""Tests for LLM utilities."""

import pytest
from unittest.mock import AsyncMock, MagicMock, Mock, patch
from datetime import datetime
import torch

from ..llm_utils import (
    NootropicLLMExtractor,
    Effect,
    Mechanism,
    SafetyData,
    CommunityInsight,
)


def _has_required_fields(result, fields):
    """Check if result has all required fields."""
    return all(field in result for field in fields)


def _check_result_fields(results, required_fields):
    """Check if all results have required fields."""
    assert len(results) > 0
    for result in results:
        assert _has_required_fields(result, required_fields)


@pytest.fixture
def mock_models():
    """Mock ML models."""
    with patch("transformers.pipeline") as ner_mock, patch("sentence_transformers.SentenceTransformer") as st_mock:
        # Mock NER model
        ner = MagicMock()
        ner_mock.return_value = ner

        # Mock sentence transformers
        analyzer = MagicMock()
        analyzer.encode.return_value = torch.randn(768)  # Standard BERT embedding size
        st_mock.return_value = analyzer

        yield {
            "ner": ner,
            "analyzer": analyzer,
        }


@pytest.fixture
def extractor(mock_models):
    """Create test extractor."""
    return NootropicLLMExtractor(
        llm_provider="test_llm",
        api_token="test_token",
    )


def test_extract_chemical_info(extractor, mock_models):
    """Test chemical information extraction."""
    # Mock NER output
    mock_models["ner"].return_value = [
        {"word": "caffeine", "entity_group": "CHEMICAL"},
        {"word": "adenosine", "entity_group": "COMPOUND"},
    ]

    # Test text
    description = """
    Caffeine (SMILES: CN1C=NC2=C1C(=O)N(C(=O)N2C)C)
    MW: 194.19
    InChI: 1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3
    """
    claims = "The compound inhibits adenosine receptors."

    # Extract info
    result = extractor.extract_chemical_info(description, claims)

    # Verify results
    assert "caffeine" in result["chemicals"]
    assert "adenosine" in result["chemicals"]
    assert "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" in result["smiles"]
    assert "1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3" in result["inchi"]
    assert 194.19 in result["molecular_weights"]

    # Check required fields
    required_fields = ["chemicals", "smiles", "inchi", "molecular_weights"]
    assert _has_required_fields(result, required_fields)


def test_analyze_text(extractor, mock_models):
    """Test text analysis."""
    # Mock embeddings
    mock_models["analyzer"].encode.return_value = torch.randn(768)

    # Test text
    title = "Novel NMDA receptor antagonist"
    abstract = "IC50 = 50 nM for NMDA receptors"
    content = """
    The compound shows good BBB penetration.
    Safety studies show no significant toxicity.
    No signs of abuse potential in animal models.
    """

    # Analyze text
    result = extractor.analyze_text(title, abstract, content)

    # Verify results
    assert "relevance_scores" in result
    assert len(result["relevance_scores"]) > 0
    assert "Contains binding/activity data" in result["key_findings"]
    assert "Contains BBB permeability data" in result["key_findings"]
    assert "Contains safety/toxicity data" in result["key_findings"]
    assert "Contains abuse potential data" in result["key_findings"]

    # Check required fields
    required_fields = ["relevance_scores", "key_findings"]
    assert _has_required_fields(result, required_fields)


@pytest.mark.asyncio
async def test_extract_effects(extractor, mock_models):
    """Test effects extraction."""
    # Mock LLM response
    with patch.object(extractor, "_extract_with_llm") as mock_llm:
        mock_llm.return_value = [
            {
                "name": "Memory Enhancement",
                "description": "Improves working memory",
                "category": "cognitive",
                "evidence": ["Clinical study shows improvement"],
            }
        ]

        # Mock embeddings
        mock_models["analyzer"].encode.return_value = torch.randn(768)

        # Test text
        text = """
        The compound enhances memory and focus.
        Clinical studies show improved cognitive performance.
        """

        # Extract effects
        effects = await extractor.extract_effects(text)

        # Verify results
        assert len(effects) > 0
        effect = effects[0]
        assert isinstance(effect, Effect)
        assert effect.name == "Memory Enhancement"
        assert effect.description == "Improves working memory"
        assert effect.category == "cognitive"
        assert len(effect.evidence) > 0
        assert effect.confidence > 0.5
        assert effect.metadata["source"] == "llm"
        assert isinstance(effect.metadata["timestamp"], str)

        # Check all effects have required fields
        required_fields = ["name", "description", "category", "confidence", "evidence"]
        _check_result_fields([e.__dict__ for e in effects], required_fields)


@pytest.mark.asyncio
async def test_extract_mechanisms(extractor, mock_models):
    """Test mechanisms extraction."""
    # Mock LLM response
    with patch.object(extractor, "_extract_with_llm") as mock_llm:
        mock_llm.return_value = [
            {
                "name": "AMPA Modulation",
                "description": "Positive modulation of AMPA receptors",
                "category": "receptor",
                "evidence": ["Binding studies confirm interaction"],
            }
        ]

        # Mock embeddings
        mock_models["analyzer"].encode.return_value = torch.randn(768)

        # Test text
        text = """
        The compound acts as a positive modulator of AMPA receptors.
        Binding studies show high affinity.
        """

        # Extract mechanisms
        mechanisms = await extractor.extract_mechanisms(text)

        # Verify results
        assert len(mechanisms) > 0
        mechanism = mechanisms[0]
        assert isinstance(mechanism, Mechanism)
        assert mechanism.name == "AMPA Modulation"
        assert mechanism.description == "Positive modulation of AMPA receptors"
        assert mechanism.category == "receptor"
        assert len(mechanism.evidence) > 0
        assert mechanism.confidence > 0.5
        assert mechanism.metadata["source"] == "llm"
        assert isinstance(mechanism.metadata["timestamp"], str)

        # Check all mechanisms have required fields
        required_fields = ["name", "description", "category", "confidence", "evidence"]
        _check_result_fields([m.__dict__ for m in mechanisms], required_fields)


@pytest.mark.asyncio
async def test_extract_safety(extractor, mock_models):
    """Test safety data extraction."""
    # Mock LLM response
    with patch.object(extractor, "_extract_with_llm") as mock_llm:
        mock_llm.return_value = [
            {
                "warning": "May cause insomnia",
                "severity": "moderate",
                "category": "side_effects",
                "evidence": ["Reported in clinical trials"],
            }
        ]

        # Mock embeddings
        mock_models["analyzer"].encode.return_value = torch.randn(768)

        # Test text
        text = """
        Side effects include insomnia in some patients.
        No serious adverse events reported.
        """

        # Extract safety data
        safety_items = await extractor.extract_safety(text)

        # Verify results
        assert len(safety_items) > 0
        safety = safety_items[0]
        assert isinstance(safety, SafetyData)
        assert safety.warning == "May cause insomnia"
        assert safety.severity == "moderate"
        assert safety.category == "side_effects"
        assert len(safety.evidence) > 0
        assert safety.confidence > 0.5
        assert safety.metadata["source"] == "llm"
        assert isinstance(safety.metadata["timestamp"], str)

        # Check all safety items have required fields
        required_fields = ["warning", "severity", "category", "confidence", "evidence"]
        _check_result_fields([s.__dict__ for s in safety_items], required_fields)


@pytest.mark.asyncio
async def test_extract_community_insights(extractor):
    """Test community insights extraction."""
    # Mock LLM response
    with patch.object(extractor, "_extract_with_llm") as mock_llm:
        mock_llm.return_value = [
            {
                "insight": "Works best with choline",
                "sentiment": "positive",
                "sources": ["Forum discussion", "User reports"],
            }
        ]

        # Test text
        text = """
        Users report best results when combined with choline.
        Generally positive experiences reported.
        """

        # Extract insights
        insights = await extractor.extract_community_insights(text)

        # Verify results
        assert len(insights) > 0
        insight = insights[0]
        assert isinstance(insight, CommunityInsight)
        assert insight.insight == "Works best with choline"
        assert insight.sentiment == "positive"
        assert len(insight.sources) > 0
        assert insight.confidence > 0.5
        assert insight.metadata["source"] == "llm"
        assert isinstance(insight.metadata["timestamp"], str)

        # Check all insights have required fields
        required_fields = ["insight", "sentiment", "confidence", "sources"]
        _check_result_fields([i.__dict__ for i in insights], required_fields)


def test_extract_evidence(extractor):
    """Test evidence extraction."""
    text = """
    First sentence about something else.
    Memory enhancement was observed in clinical trials.
    Another unrelated sentence.
    """

    evidence = extractor._extract_evidence(text, "memory enhancement")
    assert evidence == "Memory enhancement was observed in clinical trials."

    # Test with no matching evidence
    evidence = extractor._extract_evidence(text, "nonexistent term")
    assert evidence is None


def test_calculate_confidence(extractor):
    """Test confidence calculation."""
    # Full data
    full_data = {
        "name": "Effect",
        "description": "Description",
        "evidence": ["Evidence 1", "Evidence 2"],
        "sources": ["Source 1", "Source 2"],
    }
    assert extractor._calculate_confidence(full_data) == 1.0

    # Partial data
    partial_data = {
        "name": "Effect",
        "description": "Description",
    }
    assert 0.2 < extractor._calculate_confidence(partial_data) < 0.8

    # Empty data
    empty_data = {}
    assert extractor._calculate_confidence(empty_data) == 0.0


def test_error_handling(extractor, mock_models):
    """Test error handling."""
    # Mock NER failure
    mock_models["ner"].side_effect = Exception("NER error")

    # Should handle error gracefully
    result = extractor.extract_chemical_info("test", "test")
    assert result == {}

    # Mock embedding failure
    mock_models["analyzer"].encode.side_effect = Exception("Embedding error")

    # Should handle error gracefully
    result = extractor.analyze_text("test", "test", "test")
    assert result == {}


def test_model_loading_error():
    """Test error handling when models fail to load."""
    with patch("transformers.pipeline", side_effect=Exception("Model error")):
        with pytest.raises(RuntimeError, match="Required models not loaded"):
            NootropicLLMExtractor()
