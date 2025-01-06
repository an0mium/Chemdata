"""Tests for enhanced validation with comprehensive coverage."""

import pytest
from unittest.mock import AsyncMock, MagicMock, patch
from datetime import datetime, timedelta

import torch
import numpy as np

from ..validation.enhanced import EnhancedValidator, ValidationConfig, ValidationResult
from ..validation.schema import BaseSchema


@pytest.fixture
def mock_llm_extractor():
    """Mock LLM extractor with detailed responses."""
    with patch("..llm_utils.NootropicLLMExtractor") as mock:
        extractor = AsyncMock()
        # Support both structured and simple responses
        extractor.extract_effects.return_value = [
            {"name": "improved memory", "confidence": 0.8},
            {"name": "increased focus", "confidence": 0.9},
        ]
        extractor.extract_mechanisms.return_value = [
            {"name": "adenosine antagonism", "confidence": 0.95},
            {"name": "phosphodiesterase inhibition", "confidence": 0.85},
        ]
        extractor.extract_safety.return_value = [
            {"warning": "mild side effects", "severity": "low", "confidence": 0.9},
            {"warning": "potential insomnia", "severity": "moderate", "confidence": 0.8},
        ]
        mock.return_value = extractor
        yield extractor


@pytest.fixture
def mock_semantic_model():
    """Mock semantic model with realistic embeddings."""
    with patch("sentence_transformers.SentenceTransformer") as mock:
        model = MagicMock()

        # Generate realistic embeddings
        def encode_text(text):
            # Hash text to generate deterministic but unique embeddings
            seed = sum(ord(c) for c in str(text))
            np.random.seed(seed)
            return torch.from_numpy(np.random.randn(768))  # Standard BERT size

        model.encode.side_effect = encode_text
        mock.return_value = model
        yield model


@pytest.fixture
def validator(mock_llm_extractor, mock_semantic_model):
    """Create test validator with comprehensive config."""
    return EnhancedValidator(
        config=ValidationConfig(
            min_source_agreement=0.7,
            min_temporal_consistency=0.8,
            min_semantic_validity=0.6,
            min_community_consensus=0.5,
            max_temporal_gap=timedelta(days=365),
            required_sources={"pubmed", "google_scholar"},
            semantic_threshold=0.8,
        ),
        llm_provider="test_llm",
        api_token="test_token",
    )


@pytest.fixture
def test_compound_data():
    """Create test compound data with comprehensive fields."""

    class TestCompoundSchema(BaseSchema):
        def __init__(self):
            self.title = "Test Compound Effects"
            self.abstract = "Study examining effects on memory and attention..."
            self.compounds = ["caffeine", "adenosine"]
            self.effects = ["improved memory", "increased focus"]
            self.mechanisms = ["adenosine antagonism"]
            self.safety = ["well tolerated", "mild side effects"]
            self.binding_data = [
                {"target": "5-HT2A", "value": 7.2, "type": "Ki"},
                {"target": "D2", "value": 6.8, "type": "Ki"},
            ]
            self.experimental_data = {
                "solubility": 21.6,
                "logP": -0.07,
                "pKa": 14.0,
            }
            self.references = [
                "PMID:12345678",
                "DOI:10.1234/test.1",
            ]

        def __str__(self):
            return (
                f"{self.title}\n"
                f"Abstract: {self.abstract}\n"
                f"Compounds: {', '.join(self.compounds)}\n"
                f"Effects: {', '.join(self.effects)}\n"
                f"Mechanisms: {', '.join(self.mechanisms)}\n"
                f"Safety: {', '.join(self.safety)}"
            )

    return TestCompoundSchema()


@pytest.mark.asyncio
async def test_validate_success(validator, test_compound_data):
    """Test successful validation with comprehensive data."""
    # Test data from different sources with detailed information
    source_data = {
        "pubmed": {
            "title": "Similar Study",
            "abstract": "Similar findings on memory and attention...",
            "compounds": ["caffeine"],
            "effects": ["memory enhancement"],
            "mechanisms": ["adenosine receptor binding"],
            "binding_data": [
                {"target": "5-HT2A", "value": 7.1, "type": "Ki"},
            ],
            "experimental_data": {
                "solubility": 21.4,
                "logP": -0.08,
            },
        },
        "google_scholar": {
            "title": "Related Research",
            "abstract": "Consistent results on cognitive effects...",
            "compounds": ["caffeine", "adenosine"],
            "effects": ["cognitive improvement"],
            "mechanisms": ["adenosine antagonism"],
            "binding_data": [
                {"target": "5-HT2A", "value": 7.3, "type": "Ki"},
                {"target": "D2", "value": 6.9, "type": "Ki"},
            ],
        },
    }

    # Historical data points with timestamps
    temporal_data = [
        {
            "timestamp": (datetime.now() - timedelta(days=30)).isoformat(),
            "data": {
                "title": "Previous Study",
                "abstract": "Earlier findings on cognitive effects...",
                "compounds": ["caffeine"],
                "effects": ["memory improvement"],
                "mechanisms": ["adenosine antagonism"],
                "binding_data": [
                    {"target": "5-HT2A", "value": 7.2, "type": "Ki"},
                ],
            },
            "confidence": 0.9,
        },
        {
            "timestamp": (datetime.now() - timedelta(days=60)).isoformat(),
            "data": {
                "title": "Initial Study",
                "abstract": "First investigation of effects...",
                "compounds": ["caffeine"],
                "effects": ["attention enhancement"],
                "mechanisms": ["adenosine binding"],
                "binding_data": [
                    {"target": "5-HT2A", "value": 7.0, "type": "Ki"},
                ],
            },
            "confidence": 0.8,
        },
    ]

    # Community data points with confidence scores
    community_data = [
        {
            "data": {
                "title": "Community Report",
                "abstract": "User experiences with cognitive enhancement...",
                "compounds": ["caffeine"],
                "effects": ["better memory", "more focus"],
                "mechanisms": ["adenosine antagonism"],
                "safety": ["well tolerated"],
            },
            "confidence": 0.8,
            "votes": 42,
            "verified": True,
        },
        {
            "data": {
                "title": "User Feedback",
                "abstract": "Reported effects on attention...",
                "compounds": ["caffeine"],
                "effects": ["improved focus"],
                "mechanisms": ["adenosine related"],
                "safety": ["mild side effects"],
            },
            "confidence": 0.6,
            "votes": 28,
            "verified": False,
        },
    ]

    # Validate data
    result = await validator.validate(
        data=test_compound_data,
        source_data=source_data,
        temporal_data=temporal_data,
        community_data=community_data,
    )

    # Verify comprehensive results
    assert isinstance(result, ValidationResult)
    assert result.is_valid
    assert result.confidence > 0.8
    assert len(result.source_agreements) == 2
    assert all(score > 0.7 for score in result.source_agreements.values())
    assert result.temporal_consistency > 0.8
    assert result.semantic_validity > 0.6
    assert result.community_consensus > 0.5

    # Verify detailed metadata
    assert "validation_date" in result.metadata
    assert "config" in result.metadata
    assert "validation_type" in result.metadata
    assert "data_sources" in result.metadata
    assert "temporal_points" in result.metadata
    assert "community_points" in result.metadata
    assert "binding_data_consistency" in result.metadata
    assert "experimental_data_consistency" in result.metadata


@pytest.mark.asyncio
async def test_validate_missing_sources(validator, test_compound_data):
    """Test validation with missing required sources."""
    # Test data with missing source
    source_data = {
        "pubmed": {
            "title": "Similar Study",
            "abstract": "Similar findings...",
            "compounds": ["caffeine"],
            "effects": ["memory enhancement"],
            "binding_data": [
                {"target": "5-HT2A", "value": 7.1, "type": "Ki"},
            ],
        },
    }

    # Validate data
    result = await validator.validate(
        data=test_compound_data,
        source_data=source_data,
        temporal_data=[],
        community_data=[],
    )

    # Verify results
    assert not result.is_valid
    assert result.confidence < 0.7
    assert len(result.source_agreements) == 1
    assert "missing_sources" in result.metadata
    assert "google_scholar" in result.metadata["missing_sources"]


@pytest.mark.asyncio
async def test_validate_outdated_data(validator, test_compound_data):
    """Test validation with outdated temporal data."""
    # Test data with old timestamp
    temporal_data = [
        {
            "timestamp": (datetime.now() - timedelta(days=400)).isoformat(),
            "data": {
                "title": "Old Study",
                "abstract": "Outdated findings...",
                "compounds": ["caffeine"],
                "effects": ["memory improvement"],
                "binding_data": [
                    {"target": "5-HT2A", "value": 7.0, "type": "Ki"},
                ],
            },
            "confidence": 0.9,
        },
    ]

    # Validate data
    result = await validator.validate(
        data=test_compound_data,
        source_data={},
        temporal_data=temporal_data,
        community_data=[],
    )

    # Verify results
    assert not result.is_valid
    assert result.temporal_consistency < 0.8
    assert "data_age" in result.metadata
    assert result.metadata["data_age"].days > 365


@pytest.mark.asyncio
async def test_validate_semantic_mismatch(validator, test_compound_data, mock_llm_extractor):
    """Test validation with semantic mismatches."""
    # Mock inconsistent semantic extraction
    mock_llm_extractor.extract_effects.return_value = [
        {"name": "drowsiness", "confidence": 0.8},
        {"name": "confusion", "confidence": 0.7},
    ]
    mock_llm_extractor.extract_mechanisms.return_value = [
        {"name": "GABA modulation", "confidence": 0.9},
    ]
    mock_llm_extractor.extract_safety.return_value = [
        {"warning": "severe side effects", "severity": "high", "confidence": 0.85},
    ]

    # Validate data
    result = await validator.validate(
        data=test_compound_data,
        source_data={},
        temporal_data=[],
        community_data=[],
    )

    # Verify results
    assert not result.is_valid
    assert result.semantic_validity < 0.6
    assert "semantic_conflicts" in result.metadata
    assert len(result.metadata["semantic_conflicts"]) > 0


@pytest.mark.asyncio
async def test_validate_low_community_consensus(validator, test_compound_data):
    """Test validation with low community consensus."""
    # Test data with low confidence and conflicting reports
    community_data = [
        {
            "data": {
                "title": "Unreliable Report",
                "abstract": "Questionable effects...",
                "compounds": ["caffeine"],
                "effects": ["different effects"],
                "safety": ["severe concerns"],
            },
            "confidence": 0.3,
            "votes": 5,
            "verified": False,
        },
        {
            "data": {
                "title": "Disputed Claims",
                "abstract": "Controversial findings...",
                "compounds": ["caffeine"],
                "effects": ["opposite effects"],
                "safety": ["major issues"],
            },
            "confidence": 0.4,
            "votes": 3,
            "verified": False,
        },
    ]

    # Validate data
    result = await validator.validate(
        data=test_compound_data,
        source_data={},
        temporal_data=[],
        community_data=community_data,
    )

    # Verify results
    assert not result.is_valid
    assert result.community_consensus < 0.5
    assert "community_disputes" in result.metadata
    assert "low_confidence_reports" in result.metadata


@pytest.mark.asyncio
async def test_validate_error_handling(validator, test_compound_data, mock_llm_extractor):
    """Test comprehensive error handling."""
    # Mock cascading errors
    mock_llm_extractor.extract_effects.side_effect = Exception("LLM error")
    validator.semantic_model.encode.side_effect = Exception("Model error")

    # Validate data
    result = await validator.validate(
        data=test_compound_data,
        source_data={},
        temporal_data=[],
        community_data=[],
    )

    # Verify results
    assert not result.is_valid
    assert result.confidence == 0.0
    assert "errors" in result.metadata
    assert len(result.metadata["errors"]) > 0
    assert "error_count" in result.metadata
    assert "error_types" in result.metadata


def test_validation_config():
    """Test comprehensive validation configuration."""
    # Test default config
    default_config = ValidationConfig()
    assert default_config.min_source_agreement == 0.7
    assert default_config.min_temporal_consistency == 0.8
    assert default_config.min_semantic_validity == 0.6
    assert default_config.min_community_consensus == 0.5
    assert default_config.max_temporal_gap == timedelta(days=365)
    assert default_config.required_sources is None
    assert default_config.semantic_threshold == 0.8

    # Test custom config with all parameters
    custom_config = ValidationConfig(
        min_source_agreement=0.8,
        min_temporal_consistency=0.9,
        min_semantic_validity=0.7,
        min_community_consensus=0.6,
        max_temporal_gap=timedelta(days=180),
        required_sources={"source1", "source2"},
        semantic_threshold=0.9,
    )

    # Verify all custom values
    assert custom_config.min_source_agreement == 0.8
    assert custom_config.min_temporal_consistency == 0.9
    assert custom_config.min_semantic_validity == 0.7
    assert custom_config.min_community_consensus == 0.6
    assert custom_config.max_temporal_gap == timedelta(days=180)
    assert custom_config.required_sources == {"source1", "source2"}
    assert custom_config.semantic_threshold == 0.9


def test_validation_result():
    """Test comprehensive validation result."""
    # Create result with detailed metadata
    result = ValidationResult(
        is_valid=True,
        confidence=0.9,
        source_agreements={"source1": 0.8, "source2": 0.9},
        temporal_consistency=0.85,
        semantic_validity=0.75,
        community_consensus=0.7,
        validation_date=datetime.now().isoformat(),
        metadata={
            "validation_type": "enhanced",
            "data_sources": ["source1", "source2"],
            "temporal_points": 3,
            "community_points": 2,
            "binding_data_consistency": 0.95,
            "experimental_data_consistency": 0.88,
            "semantic_analysis": {
                "effect_mechanism_correlation": 0.82,
                "mechanism_safety_correlation": 0.78,
            },
            "community_analysis": {
                "verified_reports": 12,
                "total_reports": 15,
                "average_confidence": 0.76,
            },
        },
    )

    # Verify all fields and metadata
    assert result.is_valid
    assert result.confidence == 0.9
    assert result.source_agreements == {"source1": 0.8, "source2": 0.9}
    assert result.temporal_consistency == 0.85
    assert result.semantic_validity == 0.75
    assert result.community_consensus == 0.7
    assert isinstance(result.validation_date, str)

    # Verify detailed metadata
    assert result.metadata["validation_type"] == "enhanced"
    assert len(result.metadata["data_sources"]) == 2
    assert result.metadata["temporal_points"] == 3
    assert result.metadata["community_points"] == 2
    assert result.metadata["binding_data_consistency"] == 0.95
    assert result.metadata["experimental_data_consistency"] == 0.88
    assert "semantic_analysis" in result.metadata
    assert "community_analysis" in result.metadata
