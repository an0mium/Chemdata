"""Tests for web interface functionality."""

import pytest

from ..base import (
    PsychoactiveClass,
    RiskLevel,
)
from ..compound import PsychopharmCompound
from ..pipeline import PsychopharmPipeline
from ..web import WebInterface


class TestWebInterface:
    """Tests for WebInterface class."""

    def setup_method(self):
        """Set up test instance."""
        # Set up pipeline with test data
        self.pipeline = PsychopharmPipeline()
        self.test_compounds = {
            "caffeine": PsychopharmCompound(
                name="Caffeine",
                smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                cas_number="58-08-2",
            ),
            "amphetamine": PsychopharmCompound(
                name="Amphetamine",
                smiles="CC(N)CC1=CC=CC=C1",
                cas_number="300-62-9",
            ),
        }

        # Add test data
        for compound in self.test_compounds.values():
            self.pipeline.add_compound(compound)
        self.pipeline.process_compounds()

        # Initialize web interface
        self.web = WebInterface(self.pipeline)

    def test_initialization(self):
        """Test initialization of WebInterface."""
        assert self.web.pipeline is self.pipeline
        assert self.web.current_view == "list"
        assert self.web.current_filters == {}
        assert self.web.current_sort == ("name", "asc")

    def test_list_view(self):
        """Test list view functionality."""
        # Get default list view
        compounds = self.web.get_compound_list()
        assert len(compounds) == 2
        assert all(isinstance(c, dict) for c in compounds)
        assert all("name" in c for c in compounds)

        # Test filtering
        stimulants = self.web.get_compound_list(filters={"psychoactive_class": PsychoactiveClass.STIMULANT})
        assert len(stimulants) == 2  # Both are stimulants

        high_risk = self.web.get_compound_list(filters={"min_risk_level": RiskLevel.HIGH})
        assert len(high_risk) == 2  # Both have high risks

        # Test sorting
        by_name = self.web.get_compound_list(sort=("name", "asc"))
        assert by_name[0]["name"] == "Amphetamine"
        assert by_name[1]["name"] == "Caffeine"

        by_name_desc = self.web.get_compound_list(sort=("name", "desc"))
        assert by_name_desc[0]["name"] == "Caffeine"
        assert by_name_desc[1]["name"] == "Amphetamine"

    def test_detail_view(self):
        """Test detail view functionality."""
        # Get compound details
        caffeine = self.web.get_compound_details("58-08-2")  # Caffeine
        assert caffeine["name"] == "Caffeine"
        assert caffeine["cas_number"] == "58-08-2"
        assert caffeine["psychoactive_class"] == PsychoactiveClass.STIMULANT.value

        # Check detailed data sections
        assert "receptor_profiles" in caffeine
        assert "effect_profile" in caffeine
        assert "safety_alerts" in caffeine
        assert "predictions" in caffeine

        # Check visualization data
        assert "structure_svg" in caffeine
        assert "binding_plot" in caffeine
        assert "effect_plot" in caffeine
        assert "risk_plot" in caffeine

    def test_export_functionality(self, tmp_path):
        """Test export functionality."""
        # Export all compounds
        output_file = tmp_path / "all_compounds.tsv"
        self.web.export_compounds(output_file)

        assert output_file.exists()
        content = output_file.read_text()
        assert "cas_number" in content
        assert "58-08-2" in content  # Caffeine
        assert "300-62-9" in content  # Amphetamine

        # Export filtered compounds
        filtered_file = tmp_path / "filtered.tsv"
        self.web.export_compounds(filtered_file, filters={"min_risk_level": RiskLevel.HIGH})

        assert filtered_file.exists()
        filtered_content = filtered_file.read_text()
        assert "cas_number" in filtered_content
        assert "58-08-2" in filtered_content  # Has high insomnia risk
        assert "300-62-9" in filtered_content  # Has high addiction risk

    def test_relationship_visualization(self):
        """Test relationship visualization functionality."""
        # Get relationship graph data
        graph = self.web.get_relationship_graph()

        # Check graph structure
        assert "nodes" in graph
        assert "edges" in graph
        assert len(graph["nodes"]) == 2  # Two compounds
        assert len(graph["edges"]) > 0  # Should have some relationships

        # Check node data
        nodes = {n["id"]: n for n in graph["nodes"]}
        assert "58-08-2" in nodes  # Caffeine
        assert "300-62-9" in nodes  # Amphetamine
        assert all("label" in n for n in nodes.values())
        assert all("class" in n for n in nodes.values())

        # Check edge data
        assert all("source" in e for e in graph["edges"])
        assert all("target" in e for e in graph["edges"])
        assert all("type" in e for e in graph["edges"])
        assert all("weight" in e for e in graph["edges"])

    def test_search_functionality(self):
        """Test search functionality."""
        # Search by name
        results = self.web.search_compounds("caffeine")
        assert len(results) == 1
        assert results[0]["name"] == "Caffeine"

        # Search by CAS
        results = self.web.search_compounds("58-08-2")
        assert len(results) == 1
        assert results[0]["cas_number"] == "58-08-2"

        # Search by SMILES substructure
        results = self.web.search_compounds("CC1=CC=CC=C1")  # Benzene ring
        assert len(results) == 1
        assert results[0]["name"] == "Amphetamine"

        # Search by effect
        results = self.web.search_compounds("stimulation")
        assert len(results) == 2  # Both have stimulation effects


if __name__ == "__main__":
    pytest.main([__file__])
