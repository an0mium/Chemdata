"""Tests for patent search functionality."""

import pytest
from unittest.mock import patch
from datetime import datetime

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..patent_search import (
    PatentSearcher,
    PatentExtractor,
    PatentAnalyzer,
    SearchResult,
    Patent,
)


@pytest.fixture
def test_compound():
    """Create test compound fixture."""
    compound = PsychoactiveCompound(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    compound.psychoactive_class = PsychoactiveClass.STIMULANT
    compound.add_receptor_binding(
        "A2A",
        affinity=0.8,
        confidence=0.95,
        activity="antagonist"
    )
    compound.safety_alerts = {
        "anxiety": RiskLevel.MODERATE,
        "insomnia": RiskLevel.HIGH,
    }
    return compound


@pytest.fixture
def test_patent():
    """Create test patent fixture."""
    return Patent(
        number="US20210123456",
        title="Novel xanthine derivatives as adenosine receptor antagonists",
        abstract=(
            "The present invention relates to novel xanthine derivatives "
            "that act as adenosine receptor antagonists, particularly "
            "selective for the A2A receptor subtype."
        ),
        filing_date=datetime(2020, 1, 1),
        publication_date=datetime(2021, 1, 1),
        assignee="PharmaCorp Inc.",
        inventors=["John Smith", "Jane Doe"],
        claims=[
            "1. A compound of Formula I: ...",
            "2. The compound of claim 1, wherein R1 is methyl.",
        ],
        compounds=[
            {
                "name": "Compound 1",
                "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                "activity": "A2A antagonist",
                "ic50": "0.8 nM",
            }
        ],
    )


@pytest.fixture
def patent_searcher():
    """Create test patent searcher fixture."""
    return PatentSearcher()


@pytest.fixture
def patent_extractor():
    """Create test patent extractor fixture."""
    return PatentExtractor()


@pytest.fixture
def patent_analyzer():
    """Create test patent analyzer fixture."""
    return PatentAnalyzer()


class TestPatentSearcher:
    """Tests for PatentSearcher class."""

    def test_initialization(self, patent_searcher):
        """Test initialization of PatentSearcher."""
        assert patent_searcher.base_url == "https://patents.google.com/api"
        assert patent_searcher.stats == {}

    def test_search_by_text(self, patent_searcher, test_compound):
        """Test searching patents by text query."""
        # Mock API response
        mock_response = {
            "patents": [{
                "publication_number": "US20210123456",
                "title": "Novel xanthine derivatives",
                "abstract": "The present invention relates to...",
                "filing_date": "2020-01-01",
                "publication_date": "2021-01-01",
                "assignee": "PharmaCorp Inc.",
                "inventors": ["John Smith", "Jane Doe"],
            }]
        }
        
        with patch("requests.get") as mock_get:
            mock_get.return_value.json.return_value = mock_response
            
            # Search patents
            result = patent_searcher.search_by_text(
                query="adenosine receptor antagonist",
                from_date="2020-01-01",
                to_date="2021-12-31"
            )
            
            # Check result
            assert isinstance(result, SearchResult)
            assert result.success
            assert len(result.patents) > 0
            assert result.stats["total_results"] > 0
            assert result.stats["api_calls"] == 1

    def test_search_by_structure(self, patent_searcher, test_compound):
        """Test searching patents by chemical structure."""
        # Mock API response
        mock_response = {
            "patents": [{
                "publication_number": "US20210123456",
                "title": "Novel xanthine derivatives",
                "abstract": "The present invention relates to...",
                "filing_date": "2020-01-01",
                "publication_date": "2021-01-01",
                "assignee": "PharmaCorp Inc.",
                "inventors": ["John Smith", "Jane Doe"],
                "similarity_score": 0.95,
            }]
        }
        
        with patch("requests.post") as mock_post:
            mock_post.return_value.json.return_value = mock_response
            
            # Search patents
            result = patent_searcher.search_by_structure(
                smiles=test_compound.smiles,
                min_similarity=0.8
            )
            
            # Check result
            assert isinstance(result, SearchResult)
            assert result.success
            assert len(result.patents) > 0
            assert all(p.similarity_score >= 0.8 for p in result.patents)
            assert result.stats["api_calls"] == 1

    def test_error_handling(self, patent_searcher):
        """Test error handling during search."""
        with patch("requests.get") as mock_get:
            mock_get.side_effect = Exception("API Error")
            
            # Attempt search
            result = patent_searcher.search_by_text(query="test")
            
            # Check error handling
            assert not result.success
            assert "API Error" in str(result.error)
            assert result.stats["failed_requests"] == 1


class TestPatentExtractor:
    """Tests for PatentExtractor class."""

    def test_initialization(self, patent_extractor):
        """Test initialization of PatentExtractor."""
        assert patent_extractor.stats == {}

    def test_extract_compounds(self, patent_extractor, test_patent):
        """Test extracting compounds from patent."""
        # Extract compounds
        result = patent_extractor.extract_compounds(test_patent)
        
        # Check result
        assert isinstance(result, SearchResult)
        assert result.success
        assert len(result.compounds) > 0
        assert all(c.smiles for c in result.compounds)
        assert all(c.activity for c in result.compounds)
        assert result.stats["compounds_found"] > 0

    def test_extract_from_claims(self, patent_extractor, test_patent):
        """Test extracting compounds from patent claims."""
        # Extract from claims
        result = patent_extractor.extract_from_claims(test_patent)
        
        # Check result
        assert isinstance(result, SearchResult)
        assert result.success
        assert len(result.compounds) > 0
        assert all(c.smiles for c in result.compounds)
        assert result.stats["claims_processed"] > 0

    def test_extract_from_examples(self, patent_extractor, test_patent):
        """Test extracting compounds from patent examples."""
        # Extract from examples
        result = patent_extractor.extract_from_examples(test_patent)
        
        # Check result
        assert isinstance(result, SearchResult)
        assert result.success
        assert len(result.compounds) > 0
        assert all(c.smiles for c in result.compounds)
        assert result.stats["examples_processed"] > 0


class TestPatentAnalyzer:
    """Tests for PatentAnalyzer class."""

    def test_initialization(self, patent_analyzer):
        """Test initialization of PatentAnalyzer."""
        assert patent_analyzer.stats == {}

    def test_analyze_patent_family(self, patent_analyzer, test_patent):
        """Test analyzing patent family relationships."""
        # Mock family data
        mock_family = [
            test_patent,
            Patent(
                number="EP20210123456",
                title="Novel xanthine derivatives",
                filing_date=datetime(2020, 2, 1),
                publication_date=datetime(2021, 2, 1),
            ),
        ]
        
        # Analyze family
        result = patent_analyzer.analyze_patent_family(mock_family)
        
        # Check result
        assert isinstance(result, SearchResult)
        assert result.success
        assert "priority_date" in result.data
        assert "family_members" in result.data
        assert len(result.data["family_members"]) > 1

    def test_analyze_citation_network(self, patent_analyzer, test_patent):
        """Test analyzing patent citation network."""
        # Mock citation data
        mock_citations = {
            "forward": [
                Patent(
                    number="US20220123456",
                    title="Improved xanthine derivatives",
                    filing_date=datetime(2021, 1, 1),
                    publication_date=datetime(2022, 1, 1),
                ),
            ],
            "backward": [
                Patent(
                    number="US20200123456",
                    title="Original xanthine derivatives",
                    filing_date=datetime(2019, 1, 1),
                    publication_date=datetime(2020, 1, 1),
                ),
            ],
        }
        
        # Analyze citations
        result = patent_analyzer.analyze_citation_network(
            test_patent,
            citations=mock_citations
        )
        
        # Check result
        assert isinstance(result, SearchResult)
        assert result.success
        assert "forward_citations" in result.data
        assert "backward_citations" in result.data
        assert len(result.data["forward_citations"]) > 0
        assert len(result.data["backward_citations"]) > 0

    def test_analyze_assignee_portfolio(self, patent_analyzer, test_patent):
        """Test analyzing assignee patent portfolio."""
        # Mock portfolio data
        mock_portfolio = [
            test_patent,
            Patent(
                number="US20210654321",
                title="Other xanthine derivatives",
                filing_date=datetime(2020, 6, 1),
                publication_date=datetime(2021, 6, 1),
                assignee="PharmaCorp Inc.",
            ),
        ]
        
        # Analyze portfolio
        result = patent_analyzer.analyze_assignee_portfolio(
            assignee="PharmaCorp Inc.",
            patents=mock_portfolio
        )
        
        # Check result
        assert isinstance(result, SearchResult)
        assert result.success
        assert "total_patents" in result.data
        assert "technology_areas" in result.data
        assert "filing_trends" in result.data
        assert result.data["total_patents"] > 1


if __name__ == "__main__":
    pytest.main([__file__])
