"""Tests for enrich_compounds.py example script."""

import json
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from ..scripts.enrich_compounds import (
    Category,
    create_compound_from_json,
    create_compound_from_tsv,
    load_compounds,
    load_compounds_from_json,
    load_compounds_from_tsv,
    parse_category_from_comment,
    save_enriched_data,
)
from binding_data_processor.models.compound import Compound


@pytest.fixture
def example_json_path(tmp_path):
    """Create example JSON file."""
    data = {
        "compounds": [
            {
                "category": "CNS-Active (BBB Permeable)",
                "compounds": [
                    {
                        "name": "Caffeine",
                        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                        "cas_number": "58-08-2",
                    },
                ],
            },
        ],
    }
    path = tmp_path / "compounds.json"
    with open(path, "w") as f:
        json.dump(data, f)
    return path


@pytest.fixture
def example_tsv_path(tmp_path):
    """Create example TSV file."""
    content = (
        "name\tsmiles\tcas_number\n"
        "# CNS-Active (BBB Permeable)\n"
        "Caffeine\tCN1C=NC2=C1C(=O)N(C(=O)N2C)C\t58-08-2\n"
    )
    path = tmp_path / "compounds.tsv"
    with open(path, "w") as f:
        f.write(content)
    return path


def test_parse_category_from_comment():
    """Test parsing category from comment."""
    assert parse_category_from_comment("# Category") == "Category"
    assert parse_category_from_comment("#Category") == "Category"
    assert parse_category_from_comment("# ") is None
    assert parse_category_from_comment("#") is None
    assert parse_category_from_comment("") is None


def test_create_compound_from_json():
    """Test creating compound from JSON data."""
    data = {
        "name": "Caffeine",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "cas_number": "58-08-2",
    }
    compound = create_compound_from_json(data)
    assert isinstance(compound, Compound)
    assert compound.name == "Caffeine"
    assert compound.smiles == "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    assert compound.cas_number == "58-08-2"


def test_create_compound_from_tsv():
    """Test creating compound from TSV fields."""
    fields = ["Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "58-08-2"]
    compound = create_compound_from_tsv(
        fields=fields,
        name_idx=0,
        smiles_idx=1,
        cas_idx=2,
    )
    assert isinstance(compound, Compound)
    assert compound.name == "Caffeine"
    assert compound.smiles == "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    assert compound.cas_number == "58-08-2"


def test_load_compounds_from_json(example_json_path):
    """Test loading compounds from JSON file."""
    categories = load_compounds_from_json(example_json_path)
    assert len(categories) == 1
    assert isinstance(categories[0], Category)
    assert categories[0].name == "CNS-Active (BBB Permeable)"
    assert len(categories[0].compounds) == 1
    assert categories[0].compounds[0].name == "Caffeine"


def test_load_compounds_from_tsv(example_tsv_path):
    """Test loading compounds from TSV file."""
    categories = load_compounds_from_tsv(example_tsv_path)
    assert len(categories) == 1
    assert isinstance(categories[0], Category)
    assert categories[0].name == "CNS-Active (BBB Permeable)"
    assert len(categories[0].compounds) == 1
    assert categories[0].compounds[0].name == "Caffeine"


def test_load_compounds(example_json_path, example_tsv_path):
    """Test loading compounds from different file formats."""
    # Test JSON
    categories = load_compounds(example_json_path)
    assert len(categories) == 1
    assert categories[0].name == "CNS-Active (BBB Permeable)"

    # Test TSV
    categories = load_compounds(example_tsv_path)
    assert len(categories) == 1
    assert categories[0].name == "CNS-Active (BBB Permeable)"

    # Test unsupported format
    with pytest.raises(ValueError) as exc:
        load_compounds(Path("test.xyz"))
    assert "Unsupported file format" in str(exc.value)


def test_save_enriched_data(tmp_path):
    """Test saving enriched data."""
    # Create test data
    compound = Compound(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
    )
    compound.swiss_data = {
        "targets": [
            {
                "target": "Adenosine A2a receptor",
                "probability": 0.95,
            },
        ],
    }
    compound.community_data = {
        "reports": [
            {
                "source": "PsychonautWiki",
                "text": "Test report",
            },
        ],
    }
    compound.social_data = {
        "posts": [
            {
                "platform": "Reddit",
                "text": "Test post",
            },
        ],
    }
    compound.enrichment_metadata = {
        "timestamp": "2024-01-01T12:00:00",
        "category": "CNS-Active (BBB Permeable)",
        "sources": ["swiss", "community", "social"],
    }

    category = Category(
        name="CNS-Active (BBB Permeable)",
        compounds=[compound],
        description="Test category",
    )

    # Save data
    output_path = tmp_path / "enriched.json"
    save_enriched_data([category], output_path)

    # Load and verify
    with open(output_path) as f:
        data = json.load(f)

    assert "compounds" in data
    assert len(data["compounds"]) == 1
    assert data["compounds"][0]["category"] == "CNS-Active (BBB Permeable)"
    assert data["compounds"][0]["description"] == "Test category"
    assert len(data["compounds"][0]["compounds"]) == 1
    assert data["compounds"][0]["compounds"][0]["name"] == "Caffeine"
    assert "swiss_data" in data["compounds"][0]["compounds"][0]
    assert "community_data" in data["compounds"][0]["compounds"][0]
    assert "social_data" in data["compounds"][0]["compounds"][0]
    assert "enrichment_metadata" in data["compounds"][0]["compounds"][0]


def test_main_with_mock_manager(example_json_path, tmp_path):
    """Test main function with mock manager."""
    output_path = tmp_path / "enriched.json"

    # Mock manager and its methods
    mock_manager = Mock()
    mock_manager.get_metrics.return_value = {
        "processed_compounds": 1,
        "failed_compounds": 0,
        "clients": {
            "swiss": {"requests": 1},
            "community": {"requests": 1},
            "social": {"requests": 1},
        },
    }

    with patch(
        "binding_data_processor.web_enrichment.manager.WebEnrichmentManager",
        return_value=mock_manager,
    ), patch(
        "sys.argv",
        [
            "enrich_compounds.py",
            str(example_json_path),
            str(output_path),
            "--verbose",
        ],
    ):
        from ..scripts.enrich_compounds import main
        main()

    # Verify manager was used correctly
    assert mock_manager.enrich_compounds.called
    assert mock_manager.get_metrics.called
    assert mock_manager.close.called

    # Verify output file was created
    assert output_path.exists()
