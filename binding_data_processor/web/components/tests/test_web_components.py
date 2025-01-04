"""Tests for web interface components."""

import pytest
from unittest.mock import Mock
from pathlib import Path

from ....models.compound.enhanced import EnhancedCompound
from ..compound_list import CompoundList, CompoundListConfig
from ..compound_details import CompoundDetails, CompoundDetailsConfig
from ..compound_search import CompoundSearch, CompoundSearchConfig


@pytest.fixture
def test_compounds():
    """Create test compounds fixture."""
    compounds = []
    
    # Caffeine
    caffeine = EnhancedCompound(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    caffeine.molecular_weight = 194.19
    caffeine.logp = -0.07
    caffeine.hbd = 0
    caffeine.hba = 6
    caffeine.tpsa = 58.4
    compounds.append(caffeine)
    
    # Amphetamine
    amphetamine = EnhancedCompound(
        name="Amphetamine",
        smiles="CC(N)CC1=CC=CC=C1",
        cas_number="300-62-9",
    )
    amphetamine.molecular_weight = 135.21
    amphetamine.logp = 1.76
    amphetamine.hbd = 2
    amphetamine.hba = 1
    amphetamine.tpsa = 26.02
    compounds.append(amphetamine)
    
    return compounds


@pytest.fixture
def compound_list():
    """Create compound list component fixture."""
    return CompoundList(
        config=CompoundListConfig(
            page_size=10,
            show_structures=True,
            show_predictions=True,
        )
    )


@pytest.fixture
def compound_details():
    """Create compound details component fixture."""
    return CompoundDetails(
        config=CompoundDetailsConfig(
            show_structure=True,
            show_predictions=True,
            structure_width=400,
            structure_height=400,
        )
    )


@pytest.fixture
def compound_search():
    """Create compound search component fixture."""
    return CompoundSearch(
        config=CompoundSearchConfig(
            enable_text_search=True,
            enable_structure_search=True,
            enable_property_search=True,
        )
    )


class TestCompoundList:
    """Tests for CompoundList component."""

    def test_initialization(self, compound_list):
        """Test initialization of compound list."""
        assert compound_list.compounds == []
        assert compound_list.filtered_compounds == []
        assert compound_list.current_page == 1
        assert compound_list.sort_by is None

    def test_update_compounds(self, compound_list, test_compounds):
        """Test updating compound list."""
        # Update compounds
        compound_list.update_compounds(test_compounds)
        
        # Check state
        assert len(compound_list.compounds) == 2
        assert len(compound_list.filtered_compounds) == 2
        assert compound_list.current_page == 1

    def test_filtering(self, compound_list, test_compounds):
        """Test compound filtering."""
        # Update compounds
        compound_list.update_compounds(test_compounds)
        
        # Filter by source
        compound_list.set_filters({"source": "bindingdb"})
        assert len(compound_list.filtered_compounds) == 0
        
        # Filter by predictions
        compound_list.set_filters({"has_predictions": True})
        assert len(compound_list.filtered_compounds) == 0
        
        # Clear filters
        compound_list.set_filters({})
        assert len(compound_list.filtered_compounds) == 2

    def test_sorting(self, compound_list, test_compounds):
        """Test compound sorting."""
        # Update compounds
        compound_list.update_compounds(test_compounds)
        
        # Sort by name ascending
        compound_list.set_sort("name", ascending=True)
        assert compound_list.filtered_compounds[0].name == "Amphetamine"
        assert compound_list.filtered_compounds[1].name == "Caffeine"
        
        # Sort by name descending
        compound_list.set_sort("name", ascending=False)
        assert compound_list.filtered_compounds[0].name == "Caffeine"
        assert compound_list.filtered_compounds[1].name == "Amphetamine"

    def test_pagination(self, compound_list, test_compounds):
        """Test compound pagination."""
        # Update compounds
        compound_list.update_compounds(test_compounds)
        
        # Set page size
        compound_list.config.page_size = 1
        
        # Get pages
        page1 = compound_list.get_page(1)
        assert len(page1) == 1
        
        page2 = compound_list.get_page(2)
        assert len(page2) == 1
        
        # Check total pages
        assert compound_list.get_total_pages() == 2


class TestCompoundDetails:
    """Tests for CompoundDetails component."""

    def test_initialization(self, compound_details):
        """Test initialization of compound details."""
        assert compound_details.compound is None
        assert compound_details.active_tab == "info"
        assert compound_details.plot_data == {}

    def test_set_compound(self, compound_details, test_compounds):
        """Test setting compound."""
        # Set compound
        compound_details.set_compound(test_compounds[0])
        
        # Check state
        assert compound_details.compound == test_compounds[0]
        assert compound_details.plot_data != {}

    def test_basic_info(self, compound_details, test_compounds):
        """Test getting basic info."""
        # Set compound
        compound_details.set_compound(test_compounds[0])
        
        # Get info
        info = compound_details.get_basic_info()
        assert info["name"] == "Caffeine"
        assert info["smiles"] == "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        assert info["cas_number"] == "58-08-2"

    def test_export(self, compound_details, test_compounds, tmp_path):
        """Test exporting details."""
        # Set compound
        compound_details.set_compound(test_compounds[0])
        
        # Export JSON
        json_file = tmp_path / "details.json"
        compound_details.export_details(json_file, format="json")
        assert json_file.exists()
        
        # Export TSV
        tsv_file = tmp_path / "details.tsv"
        compound_details.export_details(tsv_file, format="tsv")
        assert tsv_file.exists()


class TestCompoundSearch:
    """Tests for CompoundSearch component."""

    def test_initialization(self, compound_search):
        """Test initialization of compound search."""
        assert compound_search.compounds == []
        assert compound_search.search_results == []
        assert compound_search.current_query is None
        assert compound_search.current_filters == {}

    def test_text_search(self, compound_search, test_compounds):
        """Test text search."""
        # Update compounds
        compound_search.update_compounds(test_compounds)
        
        # Search by name
        results = compound_search.text_search("caffeine")
        assert len(results) == 1
        assert results[0].name == "Caffeine"
        
        # Search by CAS
        results = compound_search.text_search("58-08-2")
        assert len(results) == 1
        assert results[0].name == "Caffeine"
        
        # Search with no matches
        results = compound_search.text_search("xyz")
        assert len(results) == 0

    def test_structure_search(self, compound_search, test_compounds):
        """Test structure search."""
        # Update compounds
        compound_search.update_compounds(test_compounds)
        
        # Search by exact structure
        results = compound_search.structure_search(
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        )
        assert len(results) == 1
        assert results[0].name == "Caffeine"
        
        # Search with no matches
        results = compound_search.structure_search("C1=CC=CC=C1")
        assert len(results) == 0

    def test_property_search(self, compound_search, test_compounds):
        """Test property search."""
        # Update compounds
        compound_search.update_compounds(test_compounds)
        
        # Search by molecular weight
        results = compound_search.property_search({
            "molecular_weight": {"min": 190, "max": 200}
        })
        assert len(results) == 1
        assert results[0].name == "Caffeine"
        
        # Search by multiple properties
        results = compound_search.property_search({
            "molecular_weight": {"min": 100, "max": 200},
            "logp": {"min": -1, "max": 0}
        })
        assert len(results) == 1
        assert results[0].name == "Caffeine"

    def test_combined_search(self, compound_search, test_compounds):
        """Test combined search."""
        # Update compounds
        compound_search.update_compounds(test_compounds)
        
        # Search by text and properties
        results = compound_search.combined_search(
            text_query="caffeine",
            property_filters={
                "molecular_weight": {"min": 190, "max": 200}
            }
        )
        assert len(results) == 1
        assert results[0].name == "Caffeine"
        
        # Search with no matches
        results = compound_search.combined_search(
            text_query="xyz",
            property_filters={
                "molecular_weight": {"min": 0, "max": 100}
            }
        )
        assert len(results) == 0


class TestComponentIntegration:
    """Tests for component integration."""

    def test_search_to_list(
        self,
        compound_search,
        compound_list,
        test_compounds,
    ):
        """Test search results to list."""
        # Update compounds
        compound_search.update_compounds(test_compounds)
        compound_list.update_compounds(test_compounds)
        
        # Perform search
        results = compound_search.text_search("caffeine")
        
        # Update list with results
        compound_list.update_compounds(results)
        
        # Check list
        assert len(compound_list.compounds) == 1
        assert compound_list.compounds[0].name == "Caffeine"

    def test_list_to_details(
        self,
        compound_list,
        compound_details,
        test_compounds,
    ):
        """Test list selection to details."""
        # Update list
        compound_list.update_compounds(test_compounds)
        
        # Select compound
        compound = compound_list.compounds[0]
        compound_details.set_compound(compound)
        
        # Check details
        assert compound_details.compound == compound
        assert compound_details.get_basic_info()["name"] == "Caffeine"

    def test_search_to_details(
        self,
        compound_search,
        compound_details,
        test_compounds,
    ):
        """Test search results to details."""
        # Update compounds
        compound_search.update_compounds(test_compounds)
        
        # Perform search
        results = compound_search.text_search("caffeine")
        
        # Show details for first result
        compound_details.set_compound(results[0])
        
        # Check details
        assert compound_details.compound == results[0]
        assert compound_details.get_basic_info()["name"] == "Caffeine"
