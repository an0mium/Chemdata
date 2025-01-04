"""Tests for data enrichment functionality."""

import pytest
from unittest.mock import patch
from pathlib import Path

from ..base import PsychoactiveClass, RiskLevel
from ..compound import PsychoactiveCompound
from ..enrichment import (
    DataEnrichmentPipeline,
    DataEnricher,
    EnrichmentConfig,
    DataSource,
    ValidationConfig,
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
def enrichment_config():
    """Create test enrichment configuration fixture."""
    return EnrichmentConfig(
        data_sources=[
            DataSource.CHEMBL,
            DataSource.PUBCHEM,
            DataSource.PSYCHONAUT_WIKI,
            DataSource.EROWID,
            DataSource.REDDIT,
            DataSource.TWITTER,
            DataSource.PATENTS,
            DataSource.PUBMED,
        ],
        validation=ValidationConfig(
            validate_web_data=True,
            validate_patent_data=True,
            validate_literature_data=True,
            min_confidence=0.8,
        ),
        cache_dir=Path("cache"),
        max_retries=3,
        timeout=30,
    )


@pytest.fixture
def enrichment_pipeline():
    """Create test enrichment pipeline fixture."""
    return DataEnrichmentPipeline()


@pytest.fixture
def data_enricher(enrichment_config):
    """Create test data enricher fixture."""
    return DataEnricher(config=enrichment_config)


class TestDataEnrichmentPipeline:
    """Tests for DataEnrichmentPipeline class."""

    def test_initialization(self, enrichment_pipeline):
        """Test initialization of DataEnrichmentPipeline."""
        assert enrichment_pipeline.chembl_client is not None
        assert enrichment_pipeline.pubchem_client is not None
        assert enrichment_pipeline.community_client is not None
        assert enrichment_pipeline.social_client is not None

    @patch('chembl_webresource_client.new_client.new_client')
    def test_chembl_enrichment(self, mock_chembl, test_compounds, enrichment_pipeline):
        """Test ChEMBL data enrichment."""
        # Mock ChEMBL response
        mock_chembl.molecule.filter.return_value = [{
            'molecule_chembl_id': 'CHEMBL113',
            'pref_name': 'CAFFEINE',
            'molecule_properties': {
                'alogp': -0.07,
                'psa': 58.44,
            },
            'molecule_type': 'Small molecule',
        }]
        mock_chembl.activity.filter.return_value = [{
            'target_chembl_id': 'CHEMBL2096',
            'standard_type': 'Ki',
            'standard_value': 13000.0,
            'standard_units': 'nM',
            'target_pref_name': 'Adenosine A2a receptor',
        }]
        
        # Enrich compound
        enrichment_pipeline.enrich_from_chembl(test_compounds[0])
        
        # Check enriched data
        assert test_compounds[0].external_ids["chembl_id"] == "CHEMBL113"
        assert test_compounds[0].properties["alogp"] == -0.07
        assert test_compounds[0].properties["psa"] == 58.44
        assert "A2A" in test_compounds[0].receptor_profiles

    @patch('pubchempy.get_compounds')
    def test_pubchem_enrichment(self, mock_pubchem, test_compounds, enrichment_pipeline):
        """Test PubChem data enrichment."""
        # Mock PubChem response
        mock_compound = type('PubChemCompound', (), {
            'cid': 2519,
            'iupac_name': '1,3,7-trimethylpurine-2,6-dione',
            'molecular_weight': 194.19,
            'xlogp': -0.07,
            'rotatable_bond_count': 0,
            'h_bond_donor_count': 0,
            'h_bond_acceptor_count': 6,
        })
        mock_pubchem.return_value = [mock_compound]
        
        # Enrich compound
        enrichment_pipeline.enrich_from_pubchem(test_compounds[0])
        
        # Check enriched data
        assert test_compounds[0].external_ids["pubchem_cid"] == "2519"
        assert test_compounds[0].properties["molecular_weight"] == 194.19
        assert test_compounds[0].properties["xlogp"] == -0.07
        assert test_compounds[0].properties["rotatable_bonds"] == 0
        assert test_compounds[0].properties["h_bond_donors"] == 0
        assert test_compounds[0].properties["h_bond_acceptors"] == 6


class TestWebEnrichment:
    """Tests for web data enrichment."""

    def test_psychonaut_wiki(self, test_compounds, data_enricher):
        """Test PsychonautWiki data enrichment."""
        # Mock PsychonautWiki API response
        mock_response = {
            "data": {
                "substances": [{
                    "name": "Caffeine",
                    "effects": ["Stimulation", "Focus", "Anxiety"],
                    "toxicity": ["Insomnia", "Tachycardia"],
                    "interactions": ["Alcohol", "MAOIs"],
                }]
            }
        }
        
        with patch("requests.get") as mock_get:
            mock_get.return_value.json.return_value = mock_response
            mock_get.return_value.status_code = 200
            
            # Enrich compounds
            result = data_enricher.web_enricher.enrich_psychonaut_wiki(
                compounds=test_compounds
            )
            
            # Check result
            assert result.success
            assert len(result.enriched_compounds) == 1
            compound = result.enriched_compounds[0]
            assert "Stimulation" in compound.web_data["psychonaut_wiki"]["effects"]
            assert "Insomnia" in compound.web_data["psychonaut_wiki"]["toxicity"]

    def test_erowid(self, test_compounds, data_enricher):
        """Test Erowid data enrichment."""
        # Mock Erowid data
        mock_reports = [
            {
                "substance": "Caffeine",
                "effects": ["Increased energy", "Mental clarity"],
                "side_effects": ["Anxiety", "Sleep issues"],
                "dosage": "200-400mg",
            }
        ]
        
        with patch.object(data_enricher.web_enricher, "_scrape_erowid") as mock_scrape:
            mock_scrape.return_value = mock_reports
            
            # Enrich compounds
            result = data_enricher.web_enricher.enrich_erowid(
                compounds=test_compounds
            )
            
            # Check result
            assert result.success
            assert len(result.enriched_compounds) == 1
            compound = result.enriched_compounds[0]
            assert "Increased energy" in compound.web_data["erowid"]["effects"]
            assert "200-400mg" == compound.web_data["erowid"]["dosage"]

    def test_reddit(self, test_compounds, data_enricher):
        """Test Reddit data enrichment."""
        # Mock Reddit API response
        mock_posts = [
            {
                "title": "Caffeine experience",
                "selftext": "Great for focus and productivity",
                "subreddit": "Nootropics",
                "score": 100,
            }
        ]
        
        with patch("praw.Reddit") as mock_reddit:
            mock_reddit.return_value.subreddit.return_value.search.return_value = (
                mock_posts
            )
            
            # Enrich compounds
            result = data_enricher.web_enricher.enrich_reddit(
                compounds=test_compounds,
                subreddits=["Nootropics", "Drugs"]
            )
            
            # Check result
            assert result.success
            assert len(result.enriched_compounds) == 1
            compound = result.enriched_compounds[0]
            assert "focus" in compound.web_data["reddit"]["mentions"][0]["text"].lower()

    def test_twitter(self, test_compounds, data_enricher):
        """Test Twitter data enrichment."""
        # Mock Twitter API response
        mock_tweets = [
            {
                "text": "Caffeine is great for studying #nootropics",
                "created_at": "2023-01-01",
                "retweet_count": 10,
                "favorite_count": 20,
            }
        ]
        
        with patch("tweepy.Client") as mock_client:
            mock_client.return_value.search_recent_tweets.return_value.data = mock_tweets
            
            # Enrich compounds
            result = data_enricher.web_enricher.enrich_twitter(
                compounds=test_compounds
            )
            
            # Check result
            assert result.success
            assert len(result.enriched_compounds) == 1
            compound = result.enriched_compounds[0]
            assert "studying" in compound.web_data["twitter"]["mentions"][0]["text"]


class TestPatentEnrichment:
    """Tests for patent data enrichment."""

    def test_patent_search(self, test_compounds, data_enricher):
        """Test patent search enrichment."""
        # Mock patent search response
        mock_patents = [
            {
                "patent_number": "US1234567",
                "title": "Caffeine formulation",
                "abstract": "Novel caffeine formulation for improved focus",
                "claims": ["A composition comprising caffeine"],
                "filing_date": "2023-01-01",
            }
        ]
        
        with patch.object(data_enricher.patent_enricher, "search_patents") as mock_search:
            mock_search.return_value = mock_patents
            
            # Enrich compounds
            result = data_enricher.patent_enricher.enrich_patents(
                compounds=test_compounds
            )
            
            # Check result
            assert result.success
            assert len(result.enriched_compounds) == 1
            compound = result.enriched_compounds[0]
            assert "US1234567" in compound.patent_data["patent_numbers"]
            assert "focus" in compound.patent_data["applications"].lower()

    def test_structure_search(self, test_compounds, data_enricher):
        """Test structure-based patent search."""
        # Mock structure search response
        mock_similar = [
            {
                "patent_number": "US7654321",
                "similarity": 0.85,
                "structure_type": "scaffold",
            }
        ]
        
        with patch.object(
            data_enricher.patent_enricher, "search_similar_structures"
        ) as mock_search:
            mock_search.return_value = mock_similar
            
            # Enrich compounds
            result = data_enricher.patent_enricher.enrich_structure_patents(
                compounds=test_compounds,
                min_similarity=0.8
            )
            
            # Check result
            assert result.success
            assert len(result.enriched_compounds) == 1
            compound = result.enriched_compounds[0]
            assert "US7654321" in compound.patent_data["similar_structures"]


class TestLiteratureEnrichment:
    """Tests for literature data enrichment."""

    def test_pubmed_search(self, test_compounds, data_enricher):
        """Test PubMed search enrichment."""
        # Mock PubMed response
        mock_articles = [
            {
                "pmid": "12345678",
                "title": "Caffeine effects on cognition",
                "abstract": "Study of caffeine's impact on focus",
                "journal": "Journal of Psychopharmacology",
                "year": 2023,
            }
        ]
        
        with patch.object(data_enricher.literature_enricher, "search_pubmed") as mock_search:
            mock_search.return_value = mock_articles
            
            # Enrich compounds
            result = data_enricher.literature_enricher.enrich_pubmed(
                compounds=test_compounds
            )
            
            # Check result
            assert result.success
            assert len(result.enriched_compounds) == 1
            compound = result.enriched_compounds[0]
            assert "12345678" in compound.literature_data["pubmed_ids"]
            assert "cognition" in compound.literature_data["topics"]

    def test_citation_analysis(self, test_compounds, data_enricher):
        """Test citation analysis."""
        # Mock citation data
        mock_citations = {
            "12345678": {
                "citing_articles": ["87654321", "98765432"],
                "citation_count": 2,
                "impact_factor": 3.5,
            }
        }
        
        with patch.object(data_enricher.literature_enricher, "analyze_citations") as mock_analyze:
            mock_analyze.return_value = mock_citations
            
            # Enrich compounds
            result = data_enricher.literature_enricher.enrich_citations(
                compounds=test_compounds
            )
            
            # Check result
            assert result.success
            assert len(result.enriched_compounds) == 1
            compound = result.enriched_compounds[0]
            assert compound.literature_data["citation_counts"]["12345678"] == 2


class TestPropertyPrediction:
    """Tests for property prediction enrichment."""

    def test_binding_prediction(self, test_compounds, data_enricher):
        """Test binding affinity prediction."""
        # Mock prediction model
        mock_predictions = {
            "A2A": (0.85, 0.9),  # (affinity, confidence)
            "D2": (0.3, 0.85),
        }
        
        with patch.object(data_enricher.property_predictor, "predict_binding") as mock_predict:
            mock_predict.return_value = mock_predictions
            
            # Enrich compounds
            result = data_enricher.property_predictor.enrich_binding_predictions(
                compounds=test_compounds
            )
            
            # Check result
            assert result.success
            assert len(result.enriched_compounds) == 1
            compound = result.enriched_compounds[0]
            assert "A2A" in compound.predicted_properties["binding_affinities"]
            assert compound.predicted_properties["binding_affinities"]["A2A"][0] == 0.85

    def test_activity_prediction(self, test_compounds, data_enricher):
        """Test activity prediction."""
        # Mock activity model
        mock_predictions = {
            "stimulant": 0.9,
            "nootropic": 0.7,
            "psychedelic": 0.1,
        }
        
        with patch.object(data_enricher.property_predictor, "predict_activity") as mock_predict:
            mock_predict.return_value = mock_predictions
            
            # Enrich compounds
            result = data_enricher.property_predictor.enrich_activity_predictions(
                compounds=test_compounds
            )
            
            # Check result
            assert result.success
            assert len(result.enriched_compounds) == 1
            compound = result.enriched_compounds[0]
            assert "stimulant" in compound.predicted_properties["activities"]
            assert compound.predicted_properties["activities"]["stimulant"] == 0.9


class TestValidation:
    """Tests for enrichment validation."""

    def test_web_data_validation(self, test_compounds, data_enricher):
        """Test web data validation."""
        # Create invalid web data
        invalid_data = {
            "effects": None,  # Should be list
            "dosage": 123,  # Should be string
            "confidence": "high",  # Should be float
        }
        
        with patch.object(data_enricher.web_enricher, "_get_web_data") as mock_get:
            mock_get.return_value = invalid_data
            
            # Try enrichment
            result = data_enricher.web_enricher.enrich_psychonaut_wiki(
                compounds=test_compounds
            )
            
            # Check validation
            assert not result.success
            assert len(result.validation_errors) == 3
            assert "effects" in str(result.validation_errors[0])
            assert "dosage" in str(result.validation_errors[1])

    def test_confidence_threshold(self, test_compounds, data_enricher):
        """Test confidence threshold validation."""
        # Mock low-confidence predictions
        mock_predictions = {
            "A2A": (0.8, 0.7),  # Below confidence threshold
            "D2": (0.3, 0.6),
        }
        
        with patch.object(data_enricher.property_predictor, "predict_binding") as mock_predict:
            mock_predict.return_value = mock_predictions
            
            # Try enrichment
            result = data_enricher.property_predictor.enrich_binding_predictions(
                compounds=test_compounds
            )
            
            # Check validation
            assert not result.success
            assert len(result.validation_errors) == 2
            assert "confidence threshold" in str(result.validation_errors[0])


class TestDataStandardization:
    """Tests for data standardization."""

    def test_standardize_web_data(self, test_compounds, data_enricher):
        """Test web data standardization."""
        # Add non-standardized data
        compound = test_compounds[0]
        compound.web_data = {
            "psychonaut": {
                "effects": ["STIMULATION", "Focus"],
                "duration": "4-6 hours"
            },
            "reddit": {
                "reported_effects": ["Energy", "FOCUS"]
            }
        }
        
        # Standardize data
        result = data_enricher.standardizer.standardize_web_data(compound)
        
        # Check result
        assert result.success
        assert all(e.islower() for e in result.data["psychonaut"]["effects"])
        assert all(e.islower() for e in result.data["reddit"]["reported_effects"])

    def test_standardize_literature_data(self, test_compounds, data_enricher):
        """Test literature data standardization."""
        # Add non-standardized data
        compound = test_compounds[0]
        compound.literature_data = {
            "topics": ["COGNITION", "Focus", "MEMORY"],
            "keywords": ["STIMULANT", "Nootropic"]
        }
        
        # Standardize data
        result = data_enricher.standardizer.standardize_literature_data(compound)
        
        # Check result
        assert result.success
        assert all(t.islower() for t in result.data["topics"])
        assert all(k.islower() for k in result.data["keywords"])


if __name__ == "__main__":
    pytest.main([__file__])
