"""Tests for data enrichment functionality."""

import pytest
from datetime import datetime

from ..enrichment import EnrichmentMixin


class TestEnrichmentMixin:
    """Tests for EnrichmentMixin class."""

    def setup_method(self):
        """Set up test instance."""
        self.enrichment = EnrichmentMixin()
        
        # Set up test data
        self.enrichment.data_sources = {
            "chembl": {
                "last_update": datetime(2023, 1, 1),
                "version": "31",
                "compounds": ["CHEMBL1234"],
            },
            "pubchem": {
                "last_update": datetime(2023, 1, 1),
                "compounds": ["CID123456"],
            },
        }
        
        self.enrichment.external_ids = {
            "chembl_id": "CHEMBL1234",
            "pubchem_cid": "CID123456",
            "drugbank_id": "DB00123",
        }
        
        self.enrichment.literature_data = {
            "pubmed_ids": ["12345678", "23456789"],
            "patent_numbers": ["US1234567A"],
            "citations": {
                "12345678": {
                    "title": "Test Study",
                    "year": 2023,
                    "journal": "Test Journal",
                    "findings": {
                        "binding_data": {
                            "5-HT2A": {"ki": 1.2, "confidence": 0.9},
                        },
                        "activity_data": {
                            "psychedelic": {"score": 0.8, "confidence": 0.9},
                        },
                    },
                },
            },
        }
        
        self.enrichment.regulatory_data = {
            "status": {
                "FDA": "Investigational",
                "EMA": "Not Approved",
            },
            "scheduling": {
                "US": "Schedule I",
                "UK": "Class A",
            },
        }

    def test_initialization(self):
        """Test initialization of EnrichmentMixin."""
        enrichment = EnrichmentMixin()
        assert enrichment.data_sources == {}
        assert enrichment.external_ids == {}
        assert enrichment.literature_data == {}
        assert enrichment.regulatory_data == {}

    def test_get_enrichment_dict(self):
        """Test get_enrichment_dict method."""
        data = self.enrichment.get_enrichment_dict()
        
        assert isinstance(data, dict)
        assert "sources" in data
        assert "identifiers" in data
        assert "literature" in data
        assert "regulatory" in data
        
        # Check source data
        sources = data["sources"]
        assert "chembl" in sources
        assert sources["chembl"]["version"] == "31"
        assert "CHEMBL1234" in sources["chembl"]["compounds"]
        
        # Check identifier data
        identifiers = data["identifiers"]
        assert identifiers["chembl_id"] == "CHEMBL1234"
        assert identifiers["pubchem_cid"] == "CID123456"
        
        # Check literature data
        literature = data["literature"]
        assert "12345678" in literature["pubmed_ids"]
        assert "US1234567A" in literature["patent_numbers"]
        
        # Check regulatory data
        regulatory = data["regulatory"]
        assert regulatory["status"]["FDA"] == "Investigational"
        assert regulatory["scheduling"]["US"] == "Schedule I"

    def test_add_literature_data(self):
        """Test add_literature_data method."""
        citation_data = {
            "pubmed_id": "34567890",
            "title": "New Study",
            "year": 2023,
            "journal": "Another Journal",
            "findings": {
                "binding_data": {
                    "D2": {"ki": 5.6, "confidence": 0.8},
                },
                "safety_data": {
                    "cardiotoxicity": {"risk": "HIGH", "confidence": 0.8},
                },
            },
        }
        
        # Add citation
        self.enrichment.add_literature_data(citation_data)
        
        # Check citation tracking
        assert "34567890" in self.enrichment.literature_data["pubmed_ids"]
        assert "34567890" in self.enrichment.literature_data["citations"]
        
        # Check citation data
        citation = self.enrichment.literature_data["citations"]["34567890"]
        assert citation["title"] == "New Study"
        assert citation["findings"]["binding_data"]["D2"]["ki"] == 5.6

    @pytest.mark.parametrize("region,status,scheduling,restrictions", [
        ("EU", "Clinical Trials", "Schedule II", ["Research Only"]),
        ("CA", "Approved", "Schedule III", ["Prescription Only"]),
        ("JP", "Under Review", "Controlled", ["Hospital Use Only"]),
    ])
    def test_add_regulatory_data(self, region, status, scheduling, restrictions):
        """Test add_regulatory_data method with different scenarios."""
        regulatory_data = {
            "region": region,
            "status": status,
            "scheduling": scheduling,
            "restrictions": restrictions,
            "last_update": datetime(2023, 2, 1),
        }
        
        # Add regulatory data
        self.enrichment.add_regulatory_data(regulatory_data)
        
        # Check regulatory data updates
        assert self.enrichment.regulatory_data["status"][region] == status
        assert self.enrichment.regulatory_data["scheduling"][region] == scheduling
        assert restrictions[0] in self.enrichment.regulatory_data.get("restrictions", [])
        
        # Check existing data preserved
        assert self.enrichment.regulatory_data["status"]["FDA"] == "Investigational"
        assert self.enrichment.regulatory_data["scheduling"]["US"] == "Schedule I"
