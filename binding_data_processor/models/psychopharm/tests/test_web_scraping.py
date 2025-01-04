"""Tests for web scraping functionality."""

import pytest
from unittest.mock import patch

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..web_scraping import (
    WebScrapingPipeline,
    PubChemScraper,
    ChEMBLScraper,
    CommunityDataScraper,
    SocialMediaScraper,
    ScrapingResult,
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
def scraping_pipeline():
    """Create test web scraping pipeline fixture."""
    return WebScrapingPipeline()


@pytest.fixture
def pubchem_scraper():
    """Create test PubChem scraper fixture."""
    return PubChemScraper()


@pytest.fixture
def chembl_scraper():
    """Create test ChEMBL scraper fixture."""
    return ChEMBLScraper()


@pytest.fixture
def community_scraper():
    """Create test community data scraper fixture."""
    return CommunityDataScraper(
        sources=["psychonaut", "erowid", "tripsit"]
    )


@pytest.fixture
def social_scraper():
    """Create test social media scraper fixture."""
    return SocialMediaScraper(
        platforms=["reddit", "twitter", "bluesky"]
    )


class TestWebScrapingPipeline:
    """Tests for WebScrapingPipeline class."""

    def test_initialization(self, scraping_pipeline):
        """Test initialization of WebScrapingPipeline."""
        assert scraping_pipeline.pubchem_scraper is not None
        assert scraping_pipeline.chembl_scraper is not None
        assert scraping_pipeline.community_scraper is not None
        assert scraping_pipeline.social_scraper is not None
        assert scraping_pipeline.stats == {}

    def test_data_integration(self, test_compound, scraping_pipeline):
        """Test integration of data from multiple sources."""
        # Mock data from different sources
        pubchem_data = {
            "iupac_name": "1,3,7-trimethylxanthine",
            "molecular_weight": 194.19,
            "3d_coords": {"x": [0.0], "y": [0.0], "z": [0.0]},
        }
        chembl_data = {
            "activities": [{
                "target_chembl_id": "CHEMBL2056",
                "standard_type": "Ki",
                "standard_value": 0.8,
                "standard_units": "nM",
            }],
        }
        community_data = {
            "psychonaut": {
                "class": "Stimulant",
                "roas": {"oral": {"common_dose": 200}},
                "effects": {
                    "positive": ["Stimulation"],
                    "negative": ["Anxiety"],
                },
            },
            "erowid": {
                "reports": [{
                    "dose": "200mg oral",
                    "effects": "Strong stimulation",
                    "duration": "4 hours",
                }],
            },
            "tripsit": {
                "dose": {"common": "200-400"},
                "effects": {
                    "positive": ["Focus"],
                    "negative": ["Insomnia"],
                },
            },
        }
        social_data = {
            "reddit": {
                "posts": [
                    {"title": "Caffeine Experience", "content": "Great for focus"},
                    {"title": "Side Effects", "content": "Some anxiety noted"},
                ],
            },
            "twitter": {
                "tweets": [
                    {"text": "Caffeine helps productivity"},
                    {"text": "Watch out for insomnia"},
                ],
            },
        }
        
        # Mock scraping methods
        with patch.multiple(
                scraping_pipeline,
                scrape_pubchem=lambda x: pubchem_data,
                scrape_chembl=lambda x: chembl_data,
                scrape_community=lambda x: community_data,
                scrape_social=lambda x: social_data):
            
            # Integrate data
            enriched_compound = scraping_pipeline.enrich_compound(test_compound)
            
            # Check integrated data
            assert enriched_compound.psychoactive_class == PsychoactiveClass.STIMULANT
            assert enriched_compound.iupac_name == "1,3,7-trimethylxanthine"
            assert "oral" in enriched_compound.dosage_data
            assert len(enriched_compound.effect_profile) > 0
            assert len(enriched_compound.experience_reports) > 0
            assert len(enriched_compound.social_mentions) > 0
            assert enriched_compound.last_updated is not None

    def test_error_handling(self, test_compound, scraping_pipeline):
        """Test error handling during scraping."""
        # Mock failed request
        with patch("requests.get") as mock_get:
            mock_get.side_effect = Exception("Connection error")
            
            # Attempt scraping
            enriched_compound = scraping_pipeline.enrich_compound(test_compound)
            
            # Check error handling
            assert enriched_compound.enrichment_errors is not None
            assert "connection_error" in scraping_pipeline.stats
            assert scraping_pipeline.stats["failed_scrapes"] > 0

    def test_data_validation(self, test_compound, scraping_pipeline):
        """Test validation of scraped data."""
        # Test data with invalid values
        invalid_data = {
            "pubchem": {
                "iupac_name": None,
                "molecular_weight": "invalid",
            },
            "community": {
                "class": "Invalid",
                "roas": {"oral": {"common_dose": "invalid"}},
                "effects": None,
            },
        }
        
        # Mock scraping with invalid data
        with patch.multiple(
                scraping_pipeline,
                scrape_pubchem=lambda x: invalid_data["pubchem"],
                scrape_community=lambda x: invalid_data["community"]):
            
            # Attempt data integration
            enriched_compound = scraping_pipeline.enrich_compound(test_compound)
            
            # Check data validation
            assert enriched_compound.iupac_name is None
            assert enriched_compound.molecular_weight is None
            assert enriched_compound.psychoactive_class == PsychoactiveClass.UNKNOWN
            assert enriched_compound.dosage_data == {}
            assert "validation_errors" in scraping_pipeline.stats


class TestPubChemScraper:
    """Tests for PubChemScraper class."""

    def test_initialization(self, pubchem_scraper):
        """Test initialization of PubChemScraper."""
        assert pubchem_scraper.base_url == "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        assert pubchem_scraper.stats == {}

    def test_fetch_compound_data(self, pubchem_scraper, test_compound):
        """Test fetching compound data from PubChem."""
        # Mock PubChem API response
        mock_response = {
            "PC_Compounds": [{
                "props": [{
                    "urn": {"label": "IUPAC Name"},
                    "value": {"sval": "1,3,7-trimethylxanthine"}
                }],
                "coords": [{
                    "type": "3d",
                    "aid": 1,
                    "conformers": [{"x": [0.0], "y": [0.0], "z": [0.0]}]
                }]
            }]
        }
        
        with patch("requests.get") as mock_get:
            mock_get.return_value.json.return_value = mock_response
            
            # Fetch data
            result = pubchem_scraper.fetch_compound_data(test_compound)
            
            # Check result
            assert isinstance(result, ScrapingResult)
            assert result.success
            assert "iupac_name" in result.data
            assert "3d_coords" in result.data
            assert result.stats["api_calls"] == 1


class TestChEMBLScraper:
    """Tests for ChEMBLScraper class."""

    def test_initialization(self, chembl_scraper):
        """Test initialization of ChEMBLScraper."""
        assert chembl_scraper.base_url == "https://www.ebi.ac.uk/chembl/api/data"
        assert chembl_scraper.stats == {}

    def test_fetch_bioactivity_data(self, chembl_scraper, test_compound):
        """Test fetching bioactivity data from ChEMBL."""
        # Mock ChEMBL API response
        mock_response = {
            "activities": [{
                "target_chembl_id": "CHEMBL2056",
                "standard_type": "Ki",
                "standard_value": 0.8,
                "standard_units": "nM",
                "target_organism": "Homo sapiens"
            }]
        }
        
        with patch("requests.get") as mock_get:
            mock_get.return_value.json.return_value = mock_response
            
            # Fetch data
            result = chembl_scraper.fetch_bioactivity_data(test_compound)
            
            # Check result
            assert isinstance(result, ScrapingResult)
            assert result.success
            assert "activities" in result.data
            assert len(result.data["activities"]) > 0
            assert result.stats["api_calls"] == 1


class TestCommunityDataScraper:
    """Tests for CommunityDataScraper class."""

    def test_initialization(self, community_scraper):
        """Test initialization of CommunityDataScraper."""
        assert "psychonaut" in community_scraper.sources
        assert "erowid" in community_scraper.sources
        assert "tripsit" in community_scraper.sources
        assert community_scraper.stats == {}

    def test_fetch_experience_reports(self, community_scraper, test_compound):
        """Test fetching experience reports."""
        # Mock community site responses
        mock_responses = {
            "psychonaut": {"reports": [{"title": "Test Report", "content": "Test content"}]},
            "erowid": {"experiences": [{"title": "Test Experience", "text": "Test text"}]},
            "tripsit": {"reports": [{"title": "Test Trip", "description": "Test desc"}]}
        }
        
        with patch.object(community_scraper, "_fetch_from_source") as mock_fetch:
            mock_fetch.side_effect = lambda src, *args: mock_responses[src]
            
            # Fetch data
            result = community_scraper.fetch_experience_reports(test_compound)
            
            # Check result
            assert isinstance(result, ScrapingResult)
            assert result.success
            assert all(src in result.data for src in community_scraper.sources)
            assert result.stats["total_reports"] > 0

    def test_fetch_safety_info(self, community_scraper, test_compound):
        """Test fetching safety information."""
        # Mock safety data responses
        mock_responses = {
            "psychonaut": {"warnings": ["Test warning"], "interactions": ["Test interaction"]},
            "tripsit": {"alerts": ["Test alert"], "combos": ["Test combo"]}
        }
        
        with patch.object(community_scraper, "_fetch_from_source") as mock_fetch:
            mock_fetch.side_effect = lambda src, *args: mock_responses[src]
            
            # Fetch data
            result = community_scraper.fetch_safety_info(test_compound)
            
            # Check result
            assert isinstance(result, ScrapingResult)
            assert result.success
            assert "warnings" in result.data
            assert "interactions" in result.data
            assert result.stats["sources_checked"] == len(community_scraper.sources)


class TestSocialMediaScraper:
    """Tests for SocialMediaScraper class."""

    def test_initialization(self, social_scraper):
        """Test initialization of SocialMediaScraper."""
        assert "reddit" in social_scraper.platforms
        assert "twitter" in social_scraper.platforms
        assert "bluesky" in social_scraper.platforms
        assert social_scraper.stats == {}

    def test_fetch_reddit_data(self, social_scraper, test_compound):
        """Test fetching data from Reddit."""
        # Mock Reddit API response
        mock_response = {
            "data": {
                "children": [
                    {"data": {"title": "Test Post", "selftext": "Test content"}},
                    {"data": {"title": "Another Post", "selftext": "More content"}}
                ]
            }
        }
        
        with patch("requests.get") as mock_get:
            mock_get.return_value.json.return_value = mock_response
            
            # Fetch data
            result = social_scraper.fetch_reddit_data(
                test_compound,
                subreddits=["researchchemicals", "nootropics"]
            )
            
            # Check result
            assert isinstance(result, ScrapingResult)
            assert result.success
            assert "posts" in result.data
            assert len(result.data["posts"]) > 0
            assert result.stats["subreddits_checked"] == 2

    def test_fetch_twitter_data(self, social_scraper, test_compound):
        """Test fetching data from Twitter."""
        # Mock Twitter API response
        mock_response = {
            "data": [
                {"text": "Test tweet about compound"},
                {"text": "Another relevant tweet"}
            ]
        }
        
        with patch("requests.get") as mock_get:
            mock_get.return_value.json.return_value = mock_response
            
            # Fetch data
            result = social_scraper.fetch_twitter_data(
                test_compound,
                days_back=7
            )
            
            # Check result
            assert isinstance(result, ScrapingResult)
            assert result.success
            assert "tweets" in result.data
            assert len(result.data["tweets"]) > 0
            assert result.stats["days_searched"] == 7

    def test_rate_limiting(self, social_scraper, test_compound):
        """Test rate limiting functionality."""
        # Make multiple requests
        for _ in range(5):
            social_scraper.fetch_reddit_data(
                test_compound,
                subreddits=["researchchemicals"]
            )
        
        # Check rate limiting stats
        assert "rate_limit_remaining" in social_scraper.stats
        assert "rate_limit_reset" in social_scraper.stats
        assert social_scraper.stats["requests_made"] == 5


if __name__ == "__main__":
    pytest.main([__file__])
