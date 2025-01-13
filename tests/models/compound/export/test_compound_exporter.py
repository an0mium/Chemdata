"""Integration tests for CompoundExporter."""

import pytest
import pandas as pd
import json
import os
from pathlib import Path
from rdkit import Chem
from datetime import datetime

from binding_data_processor.models.compound import CompoundType, CompoundData
from scripts.export_compounds import CompoundExporter


@pytest.fixture
def test_output_dir(tmp_path):
    """Create temporary output directory."""
    return tmp_path / "test_output"


@pytest.fixture
def test_compounds():
    """Create test compounds with real data."""
    compounds = [
        CompoundData(
            name="Rapamycin",
            cas_number="53123-88-9",
            smiles="CC1CCC2CC(=O)C(=C(O)C2(O1)C(=O)C1=C(O)CC(OC)C(=O)C=C1OC)C(=O)OC",
            compound_type=CompoundType.LONGEVITY,
            molecular_weight=914.172,
            logp=4.85,
            tpsa=195.5,
            hbd=3,
            hba=14,
            rotatable_bonds=6,
            literature_refs=["10.1038/nrd3439"],
            patent_refs=["US7456224B2"],
            safety_reports={"clinical_trials": 245, "adverse_events": 12},
            validation_status="validated",
        ),
        CompoundData(
            name="2C-B",
            cas_number="66142-81-2",
            smiles="CC(NC)CC1=CC(=C(C=C1)BR)OC",
            compound_type=CompoundType.PSYCHOACTIVE,
            molecular_weight=260.13,
            logp=2.3,
            tpsa=38.7,
            hbd=1,
            hba=2,
            rotatable_bonds=4,
            literature_refs=["10.1124/jpet.116.233312"],
            community_refs=["erowid_2cb_report_1"],
            safety_reports={"adverse_events": 156},
            validation_status="validated",
        ),
    ]
    return compounds


@pytest.fixture
def exporter(test_output_dir):
    """Create CompoundExporter instance."""
    exporter = CompoundExporter(cache_dir=test_output_dir / "cache")
    exporter.setup_directories()
    return exporter


@pytest.mark.asyncio
async def test_collect_compounds(exporter):
    """Test compound collection from all sources."""
    compounds = await exporter.collect_compounds()
    assert len(compounds) > 0
    assert any(c.compound_type == CompoundType.LONGEVITY for c in compounds)
    assert any(c.compound_type == CompoundType.PSYCHOACTIVE for c in compounds)


def test_standardize_structure(exporter):
    """Test structure standardization."""
    # Test with valid SMILES
    mol = Chem.MolFromSmiles("CC(NC)CC1=CC(=C(C=C1)BR)OC")
    std_mol = exporter.standardize_structure(mol)
    assert std_mol is not None
    assert Chem.MolToSmiles(std_mol) != ""

    # Test with invalid SMILES
    mol = Chem.MolFromSmiles("INVALID")
    std_mol = exporter.standardize_structure(mol)
    assert std_mol is None


@pytest.mark.asyncio
async def test_enrich_compounds(exporter, test_compounds):
    """Test compound enrichment with external data."""
    enriched = await exporter.enrich_compounds(test_compounds)
    assert len(enriched) == len(test_compounds)

    for compound in enriched:
        assert hasattr(compound, "patent_refs")
        assert hasattr(compound, "literature_refs")
        assert hasattr(compound, "community_refs")
        assert hasattr(compound, "safety_reports")


def test_merge_duplicates(exporter, test_compounds):
    """Test duplicate compound merging."""
    # Create duplicate with different data
    duplicate = test_compounds[0]
    duplicate.patent_refs = ["US9876543B2"]

    merged = exporter.merge_duplicates(test_compounds + [duplicate])
    assert len(merged) == len(test_compounds)

    # Check merged compound has combined references
    merged_compound = next(c for c in merged if c.name == "Rapamycin")
    assert len(merged_compound.patent_refs) == 2


def test_export_tsv(exporter, test_compounds, test_output_dir):
    """Test TSV export format."""
    output_path = test_output_dir / "test.tsv"
    exporter.export_tsv(test_compounds, output_path)

    assert output_path.exists()
    df = pd.read_csv(output_path, sep="\t")
    assert len(df) == len(test_compounds)
    assert all(col in df.columns for col in ["name", "cas_number", "smiles", "compound_type", "molecular_weight", "logp", "tpsa", "validation_status"])


def test_export_json(exporter, test_compounds, test_output_dir):
    """Test JSON export format."""
    output_path = test_output_dir / "test.json"
    exporter.export_json(test_compounds, output_path)

    assert output_path.exists()
    with open(output_path) as f:
        data = json.load(f)
    assert len(data) == len(test_compounds)
    assert all(isinstance(item, dict) for item in data)


def test_export_excel(exporter, test_compounds, test_output_dir):
    """Test Excel export format."""
    output_path = test_output_dir / "test.xlsx"
    exporter.export_excel(test_compounds, output_path)

    assert output_path.exists()
    df = pd.read_excel(output_path)
    assert len(df) == len(test_compounds)


def test_export_sdf(exporter, test_compounds, test_output_dir):
    """Test SDF export format."""
    output_path = test_output_dir / "test.sdf"
    exporter.export_sdf(test_compounds, output_path)

    assert output_path.exists()
    with Chem.SDMolSupplier(str(output_path)) as supplier:
        mols = [mol for mol in supplier if mol]
    assert len(mols) == len(test_compounds)


def test_export_mol_files(exporter, test_compounds, test_output_dir):
    """Test individual MOL file export."""
    output_dir = test_output_dir / "mol"
    exporter.export_mol_files(test_compounds, output_dir)

    assert output_dir.exists()
    mol_files = list(output_dir.glob("*.mol"))
    assert len(mol_files) == len(test_compounds)

    # Verify each MOL file
    for mol_file in mol_files:
        mol = Chem.MolFromMolFile(str(mol_file))
        assert mol is not None


def test_export_data_all_formats(exporter, test_compounds, test_output_dir):
    """Test export in all formats with timestamp."""
    formats = ["tsv", "json", "excel", "sdf", "mol"]
    exporter.export_data(test_compounds, test_output_dir, formats)

    timestamp = datetime.now().strftime("%Y%m%d")

    # Check each format was exported
    assert any(f.name.startswith("compounds_") and f.name.endswith(".tsv") for f in test_output_dir.glob("*.tsv"))
    assert any(f.name.startswith("compounds_") and f.name.endswith(".json") for f in test_output_dir.glob("*.json"))
    assert any(f.name.startswith("compounds_") and f.name.endswith(".xlsx") for f in test_output_dir.glob("*.xlsx"))
    assert any(f.name.startswith("compounds_") and f.name.endswith(".sdf") for f in test_output_dir.glob("*.sdf"))
    assert (test_output_dir / "mol").exists()


def test_export_performance(exporter, test_compounds, test_output_dir):
    """Test export performance with large dataset."""
    # Create large test dataset
    large_compounds = test_compounds * 1000  # 2000 compounds

    start_time = datetime.now()
    exporter.export_data(large_compounds, test_output_dir)
    duration = (datetime.now() - start_time).total_seconds()

    # Export should complete in reasonable time
    assert duration < 30  # seconds


@pytest.mark.asyncio
async def test_custom_compounds_file(exporter, test_output_dir):
    """Test loading compounds from custom file."""
    # Create test TSV file
    tsv_path = test_output_dir / "custom.tsv"
    pd.DataFrame({"name": ["Test Compound"], "cas_number": ["123-45-6"], "smiles": ["CC(=O)OC1=CC=CC=C1C(=O)O"], "compound_type": ["NOOTROPIC"]}).to_csv(
        tsv_path, sep="\t", index=False
    )

    compounds = await exporter.load_custom_compounds_file(tsv_path)
    assert len(compounds) == 1
    assert compounds[0].name == "Test Compound"


@pytest.mark.asyncio
async def test_full_export_pipeline(exporter, test_output_dir):
    """Test complete export pipeline end-to-end."""
    # Collect compounds
    compounds = await exporter.collect_compounds()
    assert len(compounds) > 0

    # Enrich compounds
    enriched = await exporter.enrich_compounds(compounds)
    assert len(enriched) == len(compounds)

    # Export in all formats
    formats = ["tsv", "json", "excel", "sdf", "mol"]
    exporter.export_data(enriched, test_output_dir, formats)

    # Verify exports
    assert len(list(test_output_dir.glob("*"))) >= len(formats)

    # Check data quality
    tsv_path = next(test_output_dir.glob("*.tsv"))
    df = pd.read_csv(tsv_path, sep="\t")
    assert not df.empty
    assert df["validation_status"].notna().all()
    assert df["last_updated"].notna().all()


def test_error_handling(exporter, test_output_dir):
    """Test error handling in export process."""
    # Test with invalid compound
    invalid_compound = CompoundData(name="Invalid", cas_number="invalid", smiles="INVALID", compound_type=CompoundType.OTHER)

    # Should handle invalid structure gracefully
    exporter.export_data([invalid_compound], test_output_dir)

    # Check error was logged but export completed
    assert any(test_output_dir.glob("*.tsv"))


def test_concurrent_exports(exporter, test_compounds, test_output_dir):
    """Test concurrent export operations."""
    import asyncio
    import concurrent.futures

    async def export_task():
        exporter.export_data(test_compounds, test_output_dir / str(id(asyncio.current_task())))

    # Run multiple exports concurrently
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        tasks = [export_task() for _ in range(4)]
        asyncio.gather(*tasks)

    # Verify all exports completed
    assert len(list(test_output_dir.glob("*/*.tsv"))) == 4
