"""Compound exporter functionality.

This module provides the CompoundExporter class for exporting compound data
in various formats with enhanced capabilities:
- TSV export with comprehensive data
- JSON export with full metadata
- Excel export for easy viewing
- SDF export for chemical structure software
- MOL files for individual structures
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Any, Union
import json
import csv
import time
import os
import asyncio
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from tqdm.asyncio import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdDeprotect, SDWriter
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Crippen

from ..base.types import CompoundType
from ..base.core import CompoundData
from ..analysis import BindingAnalyzer, ActivityAnalyzer, SafetyAnalyzer, PropertyAnalyzer, SARAnalyzer


class CompoundExporter:
    """Handles export of compound data with enrichment from multiple sources."""

    def __init__(self, cache_dir: Optional[Path] = None, n_workers: int = 4):
        self.cache_dir = cache_dir or Path("cache")
        self.n_workers = n_workers
        self.setup_standardizer()

    def setup_standardizer(self):
        """Initialize structure standardization parameters."""
        self.standardizer = rdMolStandardize.Standardizer()
        self.uncharger = rdMolStandardize.Uncharger()
        self.normalizer = rdMolStandardize.Normalizer()
        self.reionizer = rdMolStandardize.Reionizer()

    def setup_directories(self) -> Dict[str, Path]:
        """Create necessary directories."""
        base_dir = Path.cwd()
        data_dir = base_dir / "data"
        cache_dir = base_dir / "cache"
        output_dir = base_dir / "output"
        log_dir = base_dir / "logs"

        for d in [data_dir, cache_dir, output_dir, log_dir]:
            d.mkdir(parents=True, exist_ok=True)

        return {"data": data_dir, "cache": cache_dir, "output": output_dir, "logs": log_dir}

    def standardize_structure(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Enhanced structure standardization with retry logic."""
        try:
            if not mol:
                return None

            # Remove salt
            mol = self.uncharger.uncharge(mol)

            # Normalize structure
            mol = self.normalizer.normalize(mol)

            # Reionize
            mol = self.reionizer.reionize(mol)

            # Standardize
            mol = self.standardizer.standardize(mol)

            # Generate 3D conformation if needed
            if not mol.GetConformer().Is3D():
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol)

            return mol

        except Exception as e:
            print(f"Error standardizing structure: {str(e)}")
            return None

    def export_data(self, compounds: List[CompoundData], output_dir: Path, formats: Optional[List[str]] = None) -> None:
        """Export compounds in multiple formats with progress tracking."""
        formats = formats or ["tsv", "json", "excel", "sdf", "mol"]
        output_dir.mkdir(parents=True, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        with tqdm(total=len(formats), desc="Exporting formats") as pbar:
            for fmt in formats:
                try:
                    if fmt == "tsv":
                        self.export_tsv(compounds, output_dir / f"compounds_{timestamp}.tsv")
                    elif fmt == "json":
                        self.export_json(compounds, output_dir / f"compounds_{timestamp}.json")
                    elif fmt == "excel":
                        self.export_excel(compounds, output_dir / f"compounds_{timestamp}.xlsx")
                    elif fmt == "sdf":
                        self.export_sdf(compounds, output_dir / f"compounds_{timestamp}.sdf")
                    elif fmt == "mol":
                        self.export_mol_files(compounds, output_dir / "mol")
                    pbar.update(1)
                except Exception as e:
                    print(f"Error exporting {fmt} format: {str(e)}")

    def export_tsv(self, compounds: List[CompoundData], output_path: Path) -> None:
        """Export compounds to TSV format."""
        fieldnames = [
            "name",
            "cas_number",
            "smiles",
            "inchi",
            "inchi_key",
            "compound_type",
            "molecular_weight",
            "logp",
            "tpsa",
            "hbd",
            "hba",
            "rotatable_bonds",
            "target",
            "activity_type",
            "activity_value",
            "patent_count",
            "community_mentions",
            "safety_score",
            "data_quality",
            "doi_refs",
            "pmid_refs",
            "patent_refs",
            "source",
            "prediction_confidence",
            "activity_predictions",
            "toxicity_predictions",
            "abuse_predictions",
            "community_safety_score",
            "literature_confidence",
            "patent_confidence",
            "social_confidence",
            "total_references",
            "validation_status",
            "last_updated",
        ]

        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for compound in compounds:
                row = compound.to_dict()

                # Add reference counts
                row["doi_refs"] = ";".join(compound.literature_refs or [])
                row["pmid_refs"] = ";".join(compound.pubmed_refs or [])
                row["patent_refs"] = ";".join(compound.patent_refs or [])
                row["total_references"] = len((compound.literature_refs or []) + (compound.pubmed_refs or []) + (compound.patent_refs or []))

                # Add prediction data
                if hasattr(compound, "activity_predictions"):
                    row["activity_predictions"] = json.dumps(compound.activity_predictions)
                if hasattr(compound, "toxicity_predictions"):
                    row["toxicity_predictions"] = json.dumps(compound.toxicity_predictions)
                if hasattr(compound, "abuse_predictions"):
                    row["abuse_predictions"] = json.dumps(compound.abuse_predictions)

                # Add confidence scores
                row["prediction_confidence"] = getattr(compound, "prediction_confidence", None)
                row["literature_confidence"] = getattr(compound, "literature_confidence", None)
                row["patent_confidence"] = getattr(compound, "patent_confidence", None)
                row["social_confidence"] = getattr(compound, "social_confidence", None)

                # Add metadata
                row["validation_status"] = getattr(compound, "validation_status", "unvalidated")
                row["last_updated"] = datetime.now().isoformat()

                writer.writerow(row)

    def export_json(self, compounds: List[CompoundData], output_path: Path) -> None:
        """Export compounds to JSON format with full data."""
        with open(output_path, "w") as f:
            json.dump([c.to_dict(include_all=True) for c in compounds], f, indent=2)

    def export_excel(self, compounds: List[CompoundData], output_path: Path) -> None:
        """Export compounds to Excel format for easier viewing."""
        df = pd.DataFrame([c.to_dict() for c in compounds])
        df.to_excel(output_path, index=False, engine="openpyxl")

    def export_sdf(self, compounds: List[CompoundData], output_path: Path) -> None:
        """Export compounds to SDF format."""
        writer = SDWriter(str(output_path))
        for compound in compounds:
            if compound.smiles:
                mol = Chem.MolFromSmiles(compound.smiles)
                if mol:
                    mol = self.standardize_structure(mol)
                    if mol:
                        for key, value in compound.to_dict().items():
                            if value is not None:
                                mol.SetProp(key, str(value))
                        writer.write(mol)
        writer.close()

    def export_mol_files(self, compounds: List[CompoundData], output_dir: Path) -> None:
        """Export individual MOL files for each compound."""
        output_dir.mkdir(parents=True, exist_ok=True)
        for compound in compounds:
            if compound.smiles:
                mol = Chem.MolFromSmiles(compound.smiles)
                if mol:
                    mol = self.standardize_structure(mol)
                    if mol:
                        # Create a safe filename from compound name
                        safe_name = "".join(c for c in compound.name if c.isalnum() or c in (" ", "-", "_")).rstrip()
                        mol_path = output_dir / f"{safe_name}.mol"
                        try:
                            Chem.MolToMolFile(mol, str(mol_path))
                        except Exception as e:
                            print(f"Error saving MOL file for {compound.name}: {str(e)}")
