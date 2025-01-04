"""Data processing and enrichment for the dashboard.

This module provides:
1. Data loading and validation
2. Web source enrichment
3. Data cleaning and standardization
4. Format conversion and export
5. Cache management
"""

import logging
from typing import Dict, List, Optional, Set, Tuple
import pandas as pd
import numpy as np
from rdkit import Chem

from ...models.compound import CompoundData
from ...web_enrichment import WebEnrichment
from ...web_enrichment.data_sources.swiss import SwissClient
from ...web_enrichment.data_sources.chembl import ChEMBLClient
from ...web_enrichment.data_sources.pubchem import PubChemClient
from ...processors.structure import StructureProcessor


class DataProcessor:
    """Handles data processing and enrichment."""

    # Required columns for different data types
    REQUIRED_COLUMNS = {
        "structure": ["smiles", "name"],
        "activity": ["target", "activity_type", "activity_value"],
        "binding": ["target", "affinity_type", "affinity_value"],
    }

    # Column type mappings
    COLUMN_TYPES = {
        "smiles": str,
        "name": str,
        "target": str,
        "activity_type": str,
        "activity_value": float,
        "affinity_type": str,
        "affinity_value": float,
    }

    def __init__(
        self,
        web_enrichment: Optional[WebEnrichment] = None,
        swiss_client: Optional[SwissClient] = None,
        chembl_client: Optional[ChEMBLClient] = None,
        pubchem_client: Optional[PubChemClient] = None,
        structure_processor: Optional[StructureProcessor] = None,
    ):
        """Initialize data processor.

        Args:
            web_enrichment: Web enrichment client
            swiss_client: Swiss tools client
            chembl_client: ChEMBL client
            pubchem_client: PubChem client
            structure_processor: Structure processor
        """
        self.logger = logging.getLogger(__name__)
        self.web_enrichment = web_enrichment or WebEnrichment()
        self.swiss_client = swiss_client or SwissClient()
        self.chembl_client = chembl_client or ChEMBLClient()
        self.pubchem_client = pubchem_client or PubChemClient()
        self.structure_processor = structure_processor or StructureProcessor()

    def load_data(self, file_path: str, data_type: str = "structure") -> pd.DataFrame:
        """Load and validate input data.

        Args:
            file_path: Path to input file
            data_type: Type of data being loaded

        Returns:
            Validated DataFrame
        """
        try:
            # Read file based on extension
            if file_path.endswith(".csv"):
                df = pd.read_csv(file_path)
            elif file_path.endswith(".tsv"):
                df = pd.read_csv(file_path, sep="\t")
            elif file_path.endswith(".xlsx"):
                df = pd.read_excel(file_path)
            else:
                raise ValueError(f"Unsupported file format: {file_path}")

            # Validate required columns
            missing_cols = set(self.REQUIRED_COLUMNS[data_type]) - set(df.columns)
            if missing_cols:
                raise ValueError(f"Missing required columns: {missing_cols}")

            # Convert column types
            for col, dtype in self.COLUMN_TYPES.items():
                if col in df.columns:
                    df[col] = df[col].astype(dtype)

            # Validate SMILES if present
            if "smiles" in df.columns:
                df = self._validate_structures(df)

            return df

        except Exception as e:
            self.logger.error(f"Error loading data: {str(e)}")
            raise

    def enrich_data(
        self,
        df: pd.DataFrame,
        enrichment_types: List[str],
        progress_callback: Optional[callable] = None,
    ) -> pd.DataFrame:
        """Enrich data with additional information.

        Args:
            df: Input DataFrame
            enrichment_types: Types of enrichment to perform
            progress_callback: Optional callback for progress updates

        Returns:
            Enriched DataFrame
        """
        try:
            enriched_df = df.copy()

            total_compounds = len(df)
            processed = 0

            for _, row in df.iterrows():
                try:
                    # Create compound data object
                    compound = CompoundData(
                        smiles=row.get("smiles", ""),
                        name=row.get("name", ""),
                        inchi=row.get("inchi", ""),
                        inchi_key=row.get("inchi_key", ""),
                    )

                    # Perform requested enrichments
                    enriched_data = {}
                    for enrich_type in enrichment_types:
                        if enrich_type == "web":
                            web_data = self._get_web_data(compound)
                            enriched_data.update(web_data)
                        elif enrich_type == "swiss":
                            swiss_data = self._get_swiss_data(compound)
                            enriched_data.update(swiss_data)
                        elif enrich_type == "chembl":
                            chembl_data = self._get_chembl_data(compound)
                            enriched_data.update(chembl_data)
                        elif enrich_type == "pubchem":
                            pubchem_data = self._get_pubchem_data(compound)
                            enriched_data.update(pubchem_data)

                    # Update row with enriched data
                    for key, value in enriched_data.items():
                        enriched_df.at[processed, key] = value

                except Exception as e:
                    self.logger.error(
                        f"Error enriching compound {row.get('name', '')}: {str(e)}"
                    )

                processed += 1
                if progress_callback:
                    progress_callback(processed / total_compounds)

            return enriched_df

        except Exception as e:
            self.logger.error(f"Error enriching data: {str(e)}")
            return df

    def _validate_structures(self, df: pd.DataFrame) -> pd.DataFrame:
        """Validate and standardize chemical structures.

        Args:
            df: Input DataFrame

        Returns:
            DataFrame with validated structures
        """
        valid_indices = []
        for idx, row in df.iterrows():
            try:
                smiles = str(row["smiles"])
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    # Standardize SMILES
                    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
                    df.at[idx, "smiles"] = canonical_smiles
                    valid_indices.append(idx)
                else:
                    self.logger.warning(f"Invalid SMILES: {smiles}")
            except Exception as e:
                self.logger.error(f"Error validating structure: {str(e)}")

        return df.loc[valid_indices]

    def _get_web_data(self, compound: CompoundData) -> Dict:
        """Get data from web sources.

        Args:
            compound: Compound to enrich

        Returns:
            Dictionary of web data
        """
        try:
            web_data = {}

            # Get common names
            names = self.web_enrichment.get_common_names(
                compound.name,
                smiles=compound.smiles,
                inchi=compound.inchi,
            )
            if names:
                web_data["common_names"] = [n["name"] for n in names]

            # Get legal status
            legal = self.web_enrichment.get_legal_status(compound.name)
            if legal:
                web_data["legal_status"] = legal

            # Get pharmacology
            pharm = self.web_enrichment.get_pharmacology(compound.name)
            if pharm:
                web_data["pharmacology"] = pharm

            # Get reference URLs
            urls = self.web_enrichment.get_reference_urls(compound.name)
            if urls:
                web_data["reference_urls"] = urls

            return web_data

        except Exception as e:
            self.logger.error(f"Error getting web data: {str(e)}")
            return {}

    def _get_swiss_data(self, compound: CompoundData) -> Dict:
        """Get data from Swiss tools.

        Args:
            compound: Compound to enrich

        Returns:
            Dictionary of Swiss tools data
        """
        try:
            swiss_data = {}

            # Get target predictions
            targets = self.swiss_client.get_target_predictions(compound.smiles)
            if targets:
                swiss_data["target_predictions"] = targets

            # Get ADME properties
            adme = self.swiss_client.get_adme_properties(compound.smiles)
            if adme:
                swiss_data["adme_properties"] = adme

            return swiss_data

        except Exception as e:
            self.logger.error(f"Error getting Swiss data: {str(e)}")
            return {}

    def _get_chembl_data(self, compound: CompoundData) -> Dict:
        """Get data from ChEMBL.

        Args:
            compound: Compound to enrich

        Returns:
            Dictionary of ChEMBL data
        """
        try:
            chembl_data = {}

            # Get compound info
            info = self.chembl_client.get_compound_info(compound.smiles, compound.name)
            if info:
                chembl_data["chembl_info"] = info

            # Get bioactivity data
            bioactivity = self.chembl_client.get_bioactivity_data(compound.smiles)
            if bioactivity:
                chembl_data["bioactivity"] = bioactivity

            return chembl_data

        except Exception as e:
            self.logger.error(f"Error getting ChEMBL data: {str(e)}")
            return {}

    def _get_pubchem_data(self, compound: CompoundData) -> Dict:
        """Get data from PubChem.

        Args:
            compound: Compound to enrich

        Returns:
            Dictionary of PubChem data
        """
        try:
            pubchem_data = {}

            # Get compound info
            info = self.pubchem_client.get_compound_info(compound.smiles, compound.name)
            if info:
                pubchem_data["pubchem_info"] = info

            # Get bioassay data
            bioassay = self.pubchem_client.get_bioassay_data(compound.smiles)
            if bioassay:
                pubchem_data["bioassay"] = bioassay

            return pubchem_data

        except Exception as e:
            self.logger.error(f"Error getting PubChem data: {str(e)}")
            return {}

    def export_data(
        self,
        df: pd.DataFrame,
        output_path: str,
        format: str = "tsv",
        include_cols: Optional[List[str]] = None,
    ) -> bool:
        """Export processed data.

        Args:
            df: DataFrame to export
            output_path: Path to save exported file
            format: Export format (tsv, csv, xlsx)
            include_cols: Optional list of columns to include

        Returns:
            True if export successful
        """
        try:
            # Filter columns if specified
            if include_cols:
                df = df[include_cols]

            # Export based on format
            if format == "tsv":
                df.to_csv(output_path, sep="\t", index=False)
            elif format == "csv":
                df.to_csv(output_path, index=False)
            elif format == "xlsx":
                df.to_excel(output_path, index=False)
            else:
                raise ValueError(f"Unsupported export format: {format}")

            return True

        except Exception as e:
            self.logger.error(f"Error exporting data: {str(e)}")
            return False
