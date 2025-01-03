"""Binding data processing and BindingDB integration.

This module handles:
1. Loading and processing binding data from BindingDB
2. Gathering 5-HT2 receptor ligands from multiple sources
3. Enriching compound data with web sources
4. Patent searching and analysis
5. Structure validation and standardization
"""

import os
import re
import csv
import json
import requests
from typing import Any, Dict, List, Optional, Set, Tuple
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

from logger import LogManager
from models import CompoundData
from structure_utils import StructureUtils
from web_enrichment import WebEnrichment

logger = LogManager().get_logger("binding_data_processor")


class BindingDataProcessor:
    """Handles binding data processing and BindingDB integration."""

    # Target patterns for receptor ligands
    TARGET_PATTERNS = [
        # 5-HT2 receptor patterns
        "5-HT2",
        "serotonin 2",
        "HTR2",
        "5-hydroxytryptamine receptor 2",
        "serotonin receptor type 2",
        "serotonin receptor subtype 2",
        "5-HT2A",
        "5-HT2B",
        "5-HT2C",
        "HTR2A",
        "HTR2B",
        "HTR2C",
        "serotonin 2A",
        "serotonin 2B",
        "serotonin 2C",
        "HT2A_HUMAN",
        "HT2B_HUMAN",
        "HT2C_HUMAN",
        "5-hydroxytryptamine2A",
        "5-hydroxytryptamine2B",
        "5-hydroxytryptamine2C",
        "serotonin receptor 2A",
        "serotonin receptor 2B",
        "serotonin receptor 2C",
        "5-hydroxytryptamine receptor 2A",
        "5-hydroxytryptamine receptor 2B",
        "5-hydroxytryptamine receptor 2C",
        "5HT2A",
        "5HT2B",
        "5HT2C",
        "5-HT-2A",
        "5-HT-2B",
        "5-HT-2C",
        # NMDA receptor patterns
        "NMDA",
        "N-methyl-D-aspartate",
        "GRIN1",
        "GRIN2A",
        "GRIN2B",
        "GRIN2C",
        "GRIN2D",
        "GluN1",
        "GluN2A",
        "GluN2B",
        "GluN2C",
        "GluN2D",
        "NR1",
        "NR2A",
        "NR2B",
        "NR2C",
        "NR2D",
        "NMDAR1",
        "NMDAR2A",
        "NMDAR2B",
        "NMDAR2C",
        "NMDAR2D",
        "glutamate receptor ionotropic NMDA",
        "glutamate [NMDA] receptor",
        # NMDA antagonist patterns
        "3-HO-PCP",
        "3-MeO-PCP",
        "3-Methyl-PCPy",
        "4-MeO-PCP",
        "ACE mixture",
        "Agmatine",
        "Alaproclate",
        "Alazocine",
        "Amantadine",
        "AP-7",
        "AP5",
        "Apigenin",
        "Aptiganel",
        "Arketamine",
        "Atomoxetine",
        "Besonprodil",
        "Budipine",
        "Bumetanide",
        "Buphenine",
        "Carisoprodol",
        "Caroverine",
        "CGP-37849",
        "CGP-39551",
        "4-Chlorokynurenine",
        "CNQX",
        "Conantokin",
        "Coronaridine",
        "Crocetin",
        "Cyclopropane",
        "Delucemine",
        "Deschloroketamine",
        "Dextrallorphan",
        "Dextromethorphan",
        "Dextropropoxyphene",
        "Dextrorphan",
        "1,3-Diaminopropane",
        "5,7-Dichlorokynurenic acid",
        "Diethyl ether",
        "Diethylenetriamine",
        "Dieticyclidine",
        "Diphenidine",
        "Dizocilpine",
        "DNQX",
        "Eliprodil",
        "Α-Endopsychosin",
        "Enflurane",
        "Ephenidine",
        "Esketamine",
        "Esmethadone",
        "NEFA",
        "Eticyclidine",
        "EVT-101",
        "EVT-103",
        "Felbamate",
        "Flufenamic acid",
        "2-Fluorodeschloroketamine",
        "Fluorolintane",
        "Flupirtine",
        "Fourphit",
        "Furosemide",
        "Gacyclidine",
        "Gavestinel",
        "HA-966",
        "Haloperidol",
        "Halothane",
        "Hemantane",
        "Hodgkinsine",
        "Huperzine A",
        "Hydroxynorketamine",
        "Ibogaine",
        "Ibogamine",
        "Ifenprodil",
        "Indantadol",
        "Indeloxazine",
        "Isoflurane",
        "Isoxsuprine",
        "Kaitocephalin",
        "Ketamine",
        "Ketobemidone",
        "Ketofol",
        "Kynurenic acid",
        "Kynurenine",
        "L-701324",
        "Lanicemine",
        "Levomethadone",
        "Levomethorphan",
        "Levomilnacipran",
        "Levorphanol",
        "Licostinel",
        "Lubeluzole",
        "LY-235959",
        "Memantine",
        "Meprobamate",
        "Metaphit",
        "Methoxetamine",
        "Methoxphenidine",
        "18-Methoxycoronaridine",
        "Methoxyflurane",
        "Midafotel",
        "Milnacipran",
        "Minocycline",
        "Nelonemdaz",
        "Neramexane",
        "Niflumic acid",
        "Nitromemantine",
        "Nitrous oxide",
        "Noribogaine",
        "Norketamine",
        "Nortilidine",
        "NPDPA",
        "Onfasprodil",
        "Orphenadrine",
        "PCPr",
        "PD-137889",
        "PEAQX",
        "Pentamidine",
        "Perzinfotel",
        "Pethidine",
        "Phencyclidine",
        "8A-PDHQ",
        "Piretanide",
        "Promethazine",
        "Psychotridine",
        "Putrescine",
        "Racemorphan",
        "Ralfinamide",
        "Remacemide",
        "Rhynchophylline",
        "Rislenemdaz",
        "Rolicyclidine",
        "Sabeluzole",
        "Selfotel",
        "Sevoflurane",
        "SN 35210",
        "Spasmolytic A29",
        "Tabernanthine",
        "Tenocyclidine",
        "Tiletamine",
        "Tramadol",
        "Traxoprodil",
        "2,2,2-Trichloroethanol",
        "Trichloroethylene",
        "Xenon",
        "XW10508",
        "ZD-9379",
        # Nootropic targets
        "acetylcholinesterase",
        "AChE",
        "nicotinic acetylcholine receptor",
        "nAChR",
        "muscarinic acetylcholine receptor",
        "mAChR",
        "AMPA receptor",
        "GRIA",
        "GluA",
        "dopamine transporter",
        "DAT",
        "SLC6A3",
        "norepinephrine transporter",
        "NET",
        "SLC6A2",
        # Analgesic targets
        "mu opioid receptor",
        "OPRM1",
        "MOR",
        "delta opioid receptor",
        "OPRD1",
        "DOR",
        "kappa opioid receptor",
        "OPRK1",
        "KOR",
        "cannabinoid receptor",
        "CNR1",
        "CNR2",
        "CB1",
        "CB2",
        "cyclooxygenase",
        "COX-1",
        "COX-2",
        "PTGS1",
        "PTGS2",
        # Antidepressant targets
        "serotonin transporter",
        "SERT",
        "SLC6A4",
        "noradrenaline transporter",
        "NET",
        "SLC6A2",
        "monoamine oxidase",
        "MAO-A",
        "MAO-B",
        # Anxiolytic targets
        "GABA receptor",
        "GABRA",
        "GABRB",
        "GABRG",
        "benzodiazepine receptor",
        "TSPO",
        "5-HT1A",
        "HTR1A",
        "serotonin 1A",
        # Stimulant targets
        "dopamine transporter",
        "DAT",
        "SLC6A3",
        "norepinephrine transporter",
        "NET",
        "SLC6A2",
        "trace amine receptor",
        "TAAR1",
        # Psychedelic targets
        "5-HT2A",
        "HTR2A",
        "serotonin 2A",
        "5-HT1A",
        "HTR1A",
        "serotonin 1A",
        "sigma receptor",
        "SIGMAR1",
        "sigma-1",
        # Dissociative targets
        "NMDA receptor",
        "glutamate [NMDA] receptor",
        "sigma receptor",
        "SIGMAR1",
        "sigma-1",
        "kappa opioid receptor",
        "OPRK1",
        "KOR",
        # Entactogen targets
        "serotonin transporter",
        "SERT",
        "SLC6A4",
        "vesicular monoamine transporter",
        "VMAT2",
        "SLC18A2",
        "5-HT2A",
        "HTR2A",
        "serotonin 2A",
        "5-HT1A",
        "HTR1A",
        "serotonin 1A",
    ]

    # Activity type patterns with more detail
    ACTIVITY_PATTERNS = {
        "superagonist": [
            r"super.?agonist",
            r"high.?efficacy.?agonist",
            r"full.?agonist.+high.?efficacy",
            r"efficacy\s*>\s*100%",
            r"super.?potent.?agonist",
        ],
        "full_agonist": [
            r"full.?agonist",
            r"complete.?agonist",
            r"full.?receptor.?activation",
            r"efficacy\s*[~≈≃]\s*100%",
            r"maximal.?response",
        ],
        "partial_agonist": [
            r"partial.?agonist",
            r"submaximal.?activation",
            r"partial.?receptor.?activation",
            r"efficacy\s*[<≈]\s*\d{1,2}%",
            r"partial.?response",
        ],
        "weak_partial_agonist": [
            r"weak.?partial.?agonist",
            r"low.?efficacy.?partial",
            r"weak.?partial.?activation",
            r"efficacy\s*<\s*20%",
            r"minimal.?agonist",
        ],
        "mixed_agonist_antagonist": [
            r"mixed.?agonist.?antagonist",
            r"partial.?agonist.?antagonist",
            r"dual.?activity",
            r"context.?dependent",
            r"tissue.?dependent",
        ],
        "antagonist": [
            r"antagonist",
            r"blocker",
            r"inhibitor",
            r"neutral.?antagonist",
            r"competitive.?antagonist",
        ],
        "inverse_agonist": [
            r"inverse.?agonist",
            r"negative.?agonist",
            r"inverse.?activity",
            r"negative.?efficacy",
            r"constitutive.?inhibitor",
        ],
        "positive_allosteric_modulator": [
            r"positive.?allosteric",
            r"PAM",
            r"positive.?modulator",
            r"allosteric.?potentiator",
            r"positive.?cooperativity",
        ],
        "negative_allosteric_modulator": [
            r"negative.?allosteric",
            r"NAM",
            r"negative.?modulator",
            r"allosteric.?inhibitor",
            r"negative.?cooperativity",
        ],
        "complex_modulator": [
            r"complex.?modulator",
            r"mixed.?modulator",
            r"complex.?pharmacology",
            r"bitopic",
            r"dual.?mechanism",
        ],
        "enzyme_inhibitor": [
            r"enzyme.?inhibitor",
            r"inhibits?.?\w+.?enzyme",
            r"inhibits?.?\w+.?activity",
            r"reduces?.?enzyme.?activity",
            r"blocks?.?enzyme.?function",
        ],
        "enzyme_inducer": [
            r"enzyme.?inducer",
            r"induces?.?\w+.?enzyme",
            r"increases?.?enzyme.?activity",
            r"enhances?.?enzyme.?function",
            r"upregulates?.?enzyme",
        ],
    }

    def __init__(self, pubmed_client, web_client=None, http_client=None):
        """
        Initialize binding data processor.

        Args:
            pubmed_client: PubMed client for relevance scoring
            web_client: Optional WebEnrichment client for web scraping
            http_client: Optional HTTP client for Swiss tools
        """
        self.logger = LogManager().get_logger("binding_data_processor")
        self.pubmed_client = pubmed_client
        self.web_client = web_client or WebEnrichment()
        self.structure_utils = StructureUtils()

        # Initialize Swiss tools client if HTTP client provided
        from web_enrichment.data_sources.swiss import SwissClient

        self.swiss_client = SwissClient(http_client) if http_client else None

        # Path to downloaded BindingDB TSV file
        self.bindingdb_path = os.path.join(
            os.path.dirname(__file__), "BindingDB_All.tsv"
        )

        # Initialize checkpoint manager
        from checkpoint_manager import CheckpointManager

        self.checkpoint_manager = CheckpointManager()

    def _determine_activity_type(self, text: str) -> str:
        """
        Determine activity type from text description.

        Args:
            text: Text to analyze

        Returns:
            Activity type or 'unknown' if not determined
        """
        text = text.lower()
        for activity_type, patterns in self.ACTIVITY_PATTERNS.items():
            for pattern in patterns:
                if re.search(pattern, text, re.I):
                    return activity_type
        return "unknown"

    def _ensure_bindingdb_file(self) -> None:
        """Download and extract BindingDB TSV file if it doesn't exist."""
        if not os.path.exists(self.bindingdb_path):
            self.logger.info("BindingDB TSV file not found. Downloading...")
            url = "https://bindingdb.org/bind/downloads/BindingDB_All_202501_tsv.zip"
            zip_path = os.path.join(
                os.path.dirname(__file__), "BindingDB_All_202501_tsv.zip"
            )

            try:
                # Download zip file with progress bar
                response = requests.get(url, stream=True)
                response.raise_for_status()
                total_size = int(response.headers.get("content-length", 0))

                progress = tqdm(
                    total=total_size,
                    unit="iB",
                    unit_scale=True,
                    desc="Downloading BindingDB data",
                )

                with open(zip_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        size = f.write(chunk)
                        progress.update(size)
                progress.close()

                # Extract TSV file with progress
                self.logger.info("Extracting BindingDB data...")
                import zipfile

                with zipfile.ZipFile(zip_path, "r") as zip_ref:
                    # Get list of files to extract
                    files = zip_ref.namelist()
                    extract_progress = tqdm(
                        files, desc="Extracting files", unit="files"
                    )

                    for file in extract_progress:
                        zip_ref.extract(file, os.path.dirname(__file__))
                        extract_progress.set_description(f"Extracted {file}")

                # Remove zip file
                os.remove(zip_path)
                self.logger.info(
                    "BindingDB TSV file downloaded and extracted successfully."
                )

            except requests.exceptions.RequestException as e:
                self.logger.error(f"Error downloading BindingDB data: {str(e)}")
                if os.path.exists(zip_path):
                    os.remove(zip_path)
                raise
            except zipfile.BadZipFile as e:
                self.logger.error(f"Error extracting BindingDB data: {str(e)}")
                if os.path.exists(zip_path):
                    os.remove(zip_path)
                raise
            except Exception as e:
                self.logger.error(f"Unexpected error: {str(e)}")
                if os.path.exists(zip_path):
                    os.remove(zip_path)
                raise

    def load_bindingdb_data(
        self, compound_name: str, smiles: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """
        Load binding data from BindingDB TSV file for a specific compound.

        Args:
            compound_name: Name of compound to search for
            smiles: Optional SMILES string for additional matching

        Returns:
            List of binding data dictionaries containing:
            - target_common_name: Common name of target protein
            - target_protein_name: UniProt recommended name
            - target_gene_name: Gene name if available
            - affinity_value: Binding affinity value
            - affinity_unit: Units for affinity (typically nM)
            - affinity_type: Type of measurement (Ki, IC50, Kd, EC50)
            - activity_type: Type of activity (agonist, antagonist, etc.)
            - source: Data source (BindingDB)
            - doi: Article DOI if available
            - pmid: PubMed ID if available
        """
        try:
            # Read only needed columns
            needed_columns = [
                "BindingDB Ligand Name",
                "Ligand SMILES",
                "Ligand InChI",
                "Ligand InChI Key",
                "Target Name",
                "Target Source Organism According to Curator or DataSource",
                "Ki (nM)",
                "IC50 (nM)",
                "Kd (nM)",
                "EC50 (nM)",
                "kon (M-1-s-1)",
                "koff (s-1)",
                "pH",
                "Temp (C)",
                "Curation/DataSource",
                "Article DOI",
                "BindingDB Entry DOI",
                "PMID",
                "PubChem AID",
                "Patent Number",
                "Authors",
                "Institution",
                "Link to Ligand in BindingDB",
                "Link to Target in BindingDB",
                "Link to Ligand-Target Pair in BindingDB",
                "Ligand HET ID in PDB",
                "PDB ID(s) for Ligand-Target Complex",
                "PubChem CID",
                "PubChem SID",
                "ChEBI ID of Ligand",
                "ChEMBL ID of Ligand",
                "DrugBank ID of Ligand",
                "IUPHAR_GRAC ID of Ligand",
                "KEGG ID of Ligand",
                "ZINC ID of Ligand",
                "Number of Protein Chains in Target (>1 implies a multichain complex)",
                "UniProt (SwissProt) Recommended Name of Target Chain",
                "UniProt (SwissProt) Entry Name of Target Chain",
                "UniProt (SwissProt) Primary ID of Target Chain",
                "UniProt (SwissProt) Secondary ID(s) of Target Chain",
                "UniProt (SwissProt) Alternative ID(s) of Target Chain",
                "UniProt (TrEMBL) Submitted Name of Target Chain",
                "UniProt (TrEMBL) Entry Name of Target Chain",
                "UniProt (TrEMBL) Primary ID of Target Chain",
                "UniProt (TrEMBL) Secondary ID(s) of Target Chain",
                "UniProt (TrEMBL) Alternative ID(s) of Target Chain",
            ]

            # Download BindingDB TSV file if needed
            self._ensure_bindingdb_file()

            self.logger.info(
                "Reading BindingDB TSV file (this may take a few minutes)..."
            )
            total_lines = sum(1 for _ in open(self.bindingdb_path))
            self.logger.info(f"Total lines in TSV file: {total_lines:,}")
            # Read in chunks to handle large file
            chunk_size = 100000
            chunks = []
            total_chunks = 0

            # Define dtypes for problematic columns
            dtype_dict = {
                "Ki (nM)": str,
                "IC50 (nM)": str,
                "Kd (nM)": str,
                "EC50 (nM)": str,
                "kon (M-1-s-1)": str,
                "koff (s-1)": str,
                "pH": str,
                "Temp (C)": str,
                "PubChem CID": str,
                "PubChem SID": str,
                "PubChem AID": str,
                "ChEBI ID of Ligand": str,
                "ChEMBL ID of Ligand": str,
                "DrugBank ID of Ligand": str,
                "IUPHAR_GRAC ID of Ligand": str,
                "KEGG ID of Ligand": str,
                "ZINC ID of Ligand": str,
            }

            # Process chunks with progress bar
            from tqdm import tqdm

            total_chunks = (total_lines + chunk_size - 1) // chunk_size
            chunk_progress = tqdm(
                pd.read_csv(
                    self.bindingdb_path,
                    sep="\t",
                    usecols=needed_columns,
                    dtype=dtype_dict,
                    on_bad_lines="skip",
                    chunksize=chunk_size,
                    low_memory=False,
                ),
                desc="Reading BindingDB chunks",
                total=total_chunks,
                unit="chunks",
            )

            for chunk in chunk_progress:
                chunk_progress.set_description(
                    f"Reading chunk {len(chunks) + 1}/{total_chunks}"
                )
                chunks.append(chunk)

            # Combine chunks with progress bar
            self.logger.info(f"Combining {len(chunks)} chunks...")
            from tqdm import tqdm

            concat_progress = tqdm(
                total=len(chunks), desc="Combining chunks", unit="chunks"
            )

            df = pd.concat(chunks, ignore_index=True, copy=False)  # Reduce memory usage
            concat_progress.update(len(chunks))
            concat_progress.close()

            self.logger.info("BindingDB TSV file loaded successfully.")

            # Search for compound by name and SMILES
            mask = df["BindingDB Ligand Name"].str.contains(
                str(compound_name), case=False, na=False
            )
            if smiles:
                mask |= df["Ligand SMILES"] == str(smiles)

            matches = df[mask]

            binding_data = []
            for _, row in matches.iterrows():
                # Skip non-mammalian targets
                organism = str(
                    row.get(
                        "Target Source Organism According to Curator or DataSource", ""
                    )
                ).lower()
                if not any(x in organism for x in ["human", "mouse", "rat", "mammal"]):
                    self.logger.debug(f"Skipping non-mammalian target: {organism}")
                    continue

                data = {
                    "target_common_name": str(row["Target Name"]),
                    "target_protein_name": str(
                        row.get(
                            "UniProt (SwissProt) Recommended Name of Target Chain",
                            "N/A",
                        )
                    ),
                    "target_gene_name": str(
                        row.get("UniProt (SwissProt) Primary ID of Target Chain", "N/A")
                    ),
                    "affinity_value": None,
                    "affinity_unit": "nM",
                    "affinity_type": None,
                    "activity_type": self._determine_activity_type(
                        str(row.get("Assay Description", ""))
                    ),
                    "source": "BindingDB",
                    "doi": str(row.get("Article DOI", "N/A")),
                    "pmid": str(row.get("PMID", "N/A")),
                }

                # Get strongest affinity value
                for aff_type in ["Ki", "IC50", "Kd", "EC50"]:
                    col = f"{aff_type} (nM)"
                    if pd.notna(row[col]):
                        value = str(row[col])
                        # Handle '>' and '<' in affinity values
                        if value.startswith(">"):
                            # For '>' values, use the number * 10 to sort them after exact values
                            data["affinity_value"] = float(value[1:]) * 10
                            data["affinity_modifier"] = ">"
                        elif value.startswith("<"):
                            # For '<' values, use the number / 10 to sort them before exact values
                            data["affinity_value"] = float(value[1:]) / 10
                            data["affinity_modifier"] = "<"
                        else:
                            data["affinity_value"] = float(value)
                            data["affinity_modifier"] = ""
                        data["affinity_type"] = aff_type
                        break

                if data["affinity_value"] is not None:
                    binding_data.append(data)

            return binding_data

        except Exception as e:
            self.logger.error(f"Error loading BindingDB data: {str(e)}")
            return []

    def sort_binding_data(self, compound: CompoundData) -> None:
        """
        Sort binding data by PubMed relevance and organize into columns.

        Args:
            compound: Compound to process

        The method:
        1. Loads binding data from BindingDB
        2. Combines with existing binding data
        3. Scores targets by PubMed relevance
        4. Updates compound with sorted data (up to 12 targets)
        """
        # First get BindingDB data
        binding_data = self.load_bindingdb_data(compound.name, compound.smiles)

        # Add existing binding data from other sources
        for i in range(1, 13):
            target_name = getattr(compound, f"target_{i}_common_name", "N/A")
            if target_name != "N/A":
                data = {
                    "target_common_name": target_name,
                    "target_protein_name": getattr(
                        compound, f"target_{i}_protein_name", "N/A"
                    ),
                    "target_gene_name": getattr(
                        compound, f"target_{i}_gene_name", "N/A"
                    ),
                    "affinity_value": getattr(compound, f"target_{i}_affinity", "N/A"),
                    "affinity_unit": getattr(
                        compound, f"target_{i}_affinity_unit", "N/A"
                    ),
                    "affinity_type": getattr(
                        compound, f"target_{i}_affinity_type", "N/A"
                    ),
                    "activity_type": getattr(
                        compound, f"target_{i}_activity_type", "unknown"
                    ),
                    "source": getattr(compound, f"target_{i}_source", "N/A"),
                    "doi": getattr(compound, f"target_{i}_doi", "N/A"),
                    "pmid": getattr(compound, f"target_{i}_pmid", "N/A"),
                }
                binding_data.append(data)

        # Get PubMed relevance scores for each target
        scored_data = []
        for data in binding_data:
            pubmed_results = self.pubmed_client.get_binding_relevance(
                compound.name, data["target_common_name"]
            )
            scored_data.append((data, pubmed_results))

        # Sort by PubMed results
        scored_data.sort(key=lambda x: x[1], reverse=True)

        # Update compound with sorted binding data (up to 12 targets)
        for i, (data, _) in enumerate(scored_data[:12], 1):
            setattr(compound, f"target_{i}_common_name", data["target_common_name"])
            setattr(compound, f"target_{i}_protein_name", data["target_protein_name"])
            setattr(compound, f"target_{i}_gene_name", data["target_gene_name"])
            setattr(compound, f"target_{i}_affinity", str(data["affinity_value"]))
            setattr(compound, f"target_{i}_affinity_unit", data["affinity_unit"])
            setattr(compound, f"target_{i}_affinity_type", data["affinity_type"])
            setattr(compound, f"target_{i}_activity_type", data["activity_type"])
            setattr(compound, f"target_{i}_source", data["source"])
            setattr(compound, f"target_{i}_doi", data["doi"])
            setattr(compound, f"target_{i}_pmid", data["pmid"])

    def gather_ligands(self, llm_api_key: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Gather receptor ligands from multiple sources.

        Args:
            llm_api_key: Optional API key for LLM-powered analysis

        Returns:
            List of dictionaries containing compound information:
            - name: Compound name
            - smiles: SMILES structure
            - inchi: InChI identifier
            - inchi_key: InChI key
            - chembl_id: ChEMBL ID if available
            - cas_number: CAS number if found
            - binding_data: Dictionary of binding data
            - reference_urls: Dictionary of reference URLs
            - patents: List of relevant patents (if llm_api_key provided)
        """
        try:
            # Clear any existing checkpoints
            self.checkpoint_manager.clear_checkpoints()

            self.logger.info("Reading BindingDB data with enhanced columns...")
            needed_columns = [
                "BindingDB Ligand Name",
                "Ligand SMILES",
                "Ligand InChI",
                "Ligand InChI Key",
                "Target Name",
                "Target Source Organism According to Curator or DataSource",
                "Ki (nM)",
                "IC50 (nM)",
                "Kd (nM)",
                "EC50 (nM)",
                "kon (M-1-s-1)",
                "koff (s-1)",
                "pH",
                "Temp (C)",
                "Curation/DataSource",
                "Article DOI",
                "BindingDB Entry DOI",
                "PMID",
                "PubChem AID",
                "Patent Number",
                "Authors",
                "Institution",
                "Link to Ligand in BindingDB",
                "Link to Target in BindingDB",
                "Link to Ligand-Target Pair in BindingDB",
                "Ligand HET ID in PDB",
                "PDB ID(s) for Ligand-Target Complex",
                "PubChem CID",
                "PubChem SID",
                "ChEBI ID of Ligand",
                "ChEMBL ID of Ligand",
                "DrugBank ID of Ligand",
                "IUPHAR_GRAC ID of Ligand",
                "KEGG ID of Ligand",
                "ZINC ID of Ligand",
                "Number of Protein Chains in Target (>1 implies a multichain complex)",
                "UniProt (SwissProt) Recommended Name of Target Chain",
                "UniProt (SwissProt) Entry Name of Target Chain",
                "UniProt (SwissProt) Primary ID of Target Chain",
                "UniProt (SwissProt) Secondary ID(s) of Target Chain",
                "UniProt (SwissProt) Alternative ID(s) of Target Chain",
                "UniProt (TrEMBL) Submitted Name of Target Chain",
                "UniProt (TrEMBL) Entry Name of Target Chain",
                "UniProt (TrEMBL) Primary ID of Target Chain",
                "UniProt (TrEMBL) Secondary ID(s) of Target Chain",
                "UniProt (TrEMBL) Alternative ID(s) of Target Chain",
            ]

            # Download BindingDB TSV file if needed
            self._ensure_bindingdb_file()

            # Check for BindingDB parsing checkpoint
            df = None
            if self.checkpoint_manager.is_step_completed("parse_bindingdb"):
                self.logger.info("Loading parsed BindingDB data from checkpoint...")
                df = self.checkpoint_manager.load_step_data("parse_bindingdb")
            else:
                self.logger.info(
                    "Reading BindingDB TSV file (this may take a few minutes)..."
                )
                total_lines = sum(1 for _ in open(self.bindingdb_path))
                self.logger.info(f"Total lines in TSV file: {total_lines:,}")

                # Read in chunks to handle large file
                chunk_size = 100000
                chunks = []
                total_chunks = 0

                # Define dtypes for problematic columns
                dtype_dict = {
                    "Ki (nM)": str,
                    "IC50 (nM)": str,
                    "Kd (nM)": str,
                    "EC50 (nM)": str,
                    "kon (M-1-s-1)": str,
                    "koff (s-1)": str,
                    "pH": str,
                    "Temp (C)": str,
                    "PubChem CID": str,
                    "PubChem SID": str,
                    "PubChem AID": str,
                    "ChEBI ID of Ligand": str,
                    "ChEMBL ID of Ligand": str,
                    "DrugBank ID of Ligand": str,
                    "IUPHAR_GRAC ID of Ligand": str,
                    "KEGG ID of Ligand": str,
                    "ZINC ID of Ligand": str,
                }

                for chunk in pd.read_csv(
                    self.bindingdb_path,
                    sep="\t",
                    usecols=needed_columns,
                    dtype=dtype_dict,
                    on_bad_lines="skip",
                    chunksize=chunk_size,
                    low_memory=False,
                ):
                    total_chunks += 1
                    self.logger.info(
                        f"Processing chunk {total_chunks} ({chunk_size:,} rows per chunk)"
                    )
                    chunks.append(chunk)

                self.logger.info(f"Combining {len(chunks)} chunks...")
                df = pd.concat(chunks, ignore_index=True)
                self.logger.info("BindingDB TSV file loaded successfully.")

                # Save parsing checkpoint
                self.checkpoint_manager.save_checkpoint("parse_bindingdb", df)

            # Find 5-HT2 receptor targets
            matches = None
            if self.checkpoint_manager.is_step_completed("find_targets"):
                self.logger.info("Loading target matches from checkpoint...")
                matches = self.checkpoint_manager.load_step_data("find_targets")
            else:
                self.logger.info(
                    "Searching for 5-HT2 receptor targets in BindingDB data..."
                )
                pattern = "|".join(self.TARGET_PATTERNS)

                # Convert columns to string and handle NaN values
                self.logger.info("Converting columns to string format...")
                df["Target Name"] = df["Target Name"].fillna("").astype(str)
                df["UniProt (SwissProt) Recommended Name of Target Chain"] = (
                    df["UniProt (SwissProt) Recommended Name of Target Chain"]
                    .fillna("")
                    .astype(str)
                )
                df["UniProt (SwissProt) Primary ID of Target Chain"] = (
                    df["UniProt (SwissProt) Primary ID of Target Chain"]
                    .fillna("")
                    .astype(str)
                )

                # Compile regex pattern once
                import re

                compiled_pattern = re.compile(pattern, re.IGNORECASE)

                # Search in all relevant columns with enhanced progress tracking
                from tqdm import tqdm

                self.logger.info("Searching for 5-HT2 patterns in target columns...")
                total_rows = len(df)
                matches_mask = pd.Series(False, index=df.index)

                columns = [
                    "Target Name",
                    "UniProt (SwissProt) Recommended Name of Target Chain",
                    "UniProt (SwissProt) Primary ID of Target Chain",
                ]
                for col in tqdm(columns, desc="Searching columns", unit="col"):
                    self.logger.info(f"Searching column: {col}")
                    # Process in chunks with progress bar
                    chunk_size = 10000
                    chunks = range(0, total_rows, chunk_size)
                    with tqdm(
                        total=total_rows, desc=f"Processing {col}", unit="rows"
                    ) as pbar:
                        for i in chunks:
                            chunk = df[col].iloc[i : i + chunk_size]
                            # Use vectorized operations where possible
                            chunk_matches = (
                                chunk.fillna("")
                                .astype(str)
                                .apply(lambda x: bool(compiled_pattern.search(x)))
                            )
                            matches_mask.iloc[i : i + chunk_size] |= chunk_matches
                            pbar.update(len(chunk))

                matches = df[matches_mask]
                self.logger.info(f"Found {len(matches):,} initial matches in BindingDB")

                # Save target matches checkpoint
                self.checkpoint_manager.save_checkpoint("find_targets", matches)

            self.logger.info("Filtering matches by organism and structure...")

            # Process compounds with enhanced progress tracking
            from tqdm import tqdm

            compounds = []
            seen_inchikeys = set()
            processed_count = 0
            compound_groups = list(matches.groupby("BindingDB Ligand Name"))
            total_compounds = len(compound_groups)

            self.logger.info(f"Processing {total_compounds:,} unique compounds...")
            progress_bar = tqdm(
                total=total_compounds, desc="Processing compounds", unit="compounds"
            )

            # Check for compound processing checkpoint
            if self.checkpoint_manager.is_step_completed("process_compounds"):
                self.logger.info("Loading processed compounds from checkpoint...")
                compounds = self.checkpoint_manager.load_step_data("process_compounds")
                processed_count = len(compounds)

            # Process compounds with progress bar
            from tqdm import tqdm

            progress_bar = tqdm(
                compound_groups,
                desc="Processing compounds",
                unit="compounds",
                total=total_compounds,
            )

            for name, group in progress_bar:
                processed_count += 1
                progress_bar.set_description(f"Processing {name[:30]}...")

                # Save intermediate checkpoint every 100 compounds
                if processed_count % 100 == 0:
                    self.checkpoint_manager.save_checkpoint(
                        "process_compounds",
                        compounds,
                        metadata={"processed_count": processed_count},
                    )

                # Skip if we've seen this compound
                inchikey = str(group.iloc[0]["Ligand InChI Key"])
                if pd.notna(inchikey) and inchikey in seen_inchikeys:
                    continue
                seen_inchikeys.add(inchikey)

                # Get all affinity values
                affinities = []
                for _, row in group.iterrows():
                    for aff_type in ["Ki", "IC50", "Kd", "EC50"]:
                        col = f"{aff_type} (nM)"
                        if pd.notna(row[col]):
                            value = str(row[col])
                            # Handle '>' and '<' in affinity values
                            if value.startswith(">"):
                                value = float(value[1:]) * 10
                                modifier = ">"
                            elif value.startswith("<"):
                                value = float(value[1:]) / 10
                                modifier = "<"
                            else:
                                value = float(value)
                                modifier = ""

                            # Include assay conditions
                            conditions = {
                                "pH": row.get("pH", "N/A"),
                                "temperature": row.get("Temp (C)", "N/A"),
                                "source": row.get("Curation/DataSource", "N/A"),
                            }

                            affinities.append(
                                {
                                    "type": aff_type,
                                    "value": value,
                                    "modifier": modifier,
                                    "conditions": conditions,
                                }
                            )

                if affinities:
                    # Sort by value (lower is better)
                    affinities.sort(key=lambda x: x["value"])
                    best_affinity = affinities[0]

                    # Get first row for compound info
                    row = group.iloc[0]

                    # Convert SMILES to mol for structure analysis
                    mol = None
                    if pd.notna(row["Ligand SMILES"]):
                        mol = Chem.MolFromSmiles(str(row["Ligand SMILES"]))

                        # Check if molecule has 5-HT2 ligand-like features
                        if mol is not None:
                            is_potential, patterns = (
                                self.structure_utils.is_potential_ligand(mol)
                            )
                            if not is_potential:
                                continue

                    # Clean compound name for web searches
                    clean_name = self._clean_name(str(name))

                    # Get identifiers using name utils
                    from web_enrichment.name_utils import (
                        extract_identifiers,
                        clean_name as clean_name_func,
                    )

                    # Try multiple name variations
                    name_variations = [
                        str(name),  # Original name
                        clean_name,  # Cleaned name
                        str(row.get("BindingDB Ligand Name", "")),  # BindingDB name
                        str(row.get("ChEMBL ID of Ligand", "")),  # ChEMBL ID
                        str(row.get("PubChem CID", "")),  # PubChem CID
                    ]

                    # Try each variation until we get identifiers
                    identifiers = None
                    for variation in name_variations:
                        if variation and variation.strip():
                            try:
                                clean_name, chembl_id, cas_number, patent_id = (
                                    extract_identifiers(variation)
                                )
                                if any([chembl_id, cas_number, patent_id]):
                                    identifiers = (
                                        clean_name,
                                        chembl_id,
                                        cas_number,
                                        patent_id,
                                    )
                                    break
                            except Exception as e:
                                self.logger.debug(
                                    f"Error extracting identifiers from {variation}: {str(e)}"
                                )
                                continue

                    # If no identifiers found, use cleaned original name
                    if not identifiers:
                        clean_name = clean_name_func(clean_name)  # Additional cleaning
                        identifiers = (clean_name, None, None, None)
                    else:
                        clean_name, chembl_id, cas_number, patent_id = identifiers

                    # If patent example, extract identifiers using LLM
                    if patent_id and llm_api_key:
                        example_match = re.search(
                            r"(?:Example|Compound)\s*([A-Za-z0-9-]+)", str(name)
                        )
                        if example_match:
                            example_id = f"Example {example_match.group(1)}"
                            from web_enrichment.llm_utils import extract_patent_compound

                            patent_data = extract_patent_compound(
                                example_id, llm_api_key
                            )
                            if patent_data:
                                cas_number = patent_data.get("cas_number", cas_number)
                                smiles = patent_data.get(
                                    "smiles", str(row["Ligand SMILES"])
                                )
                                inchi = patent_data.get(
                                    "inchi", str(row["Ligand InChI"])
                                )
                                clean_name = patent_data.get("name", clean_name)

                    self.logger.info(f"Processing compound: {clean_name}")

                    try:
                        # Get PubChem data
                        self.logger.info("Fetching PubChem data...")
                        pubchem_data = self.web_client._get_pubchem_data(
                            cas=cas_number,
                            smiles=str(row["Ligand SMILES"]),
                            inchi=str(row["Ligand InChI"]),
                        )

                        self.logger.info("Creating compound data structure...")

                        # Get Swiss tools data if available
                        swiss_data = {}
                        if self.swiss_client and pd.notna(row["Ligand SMILES"]):
                            self.logger.info("Getting Swiss tools predictions...")

                            # Get target predictions
                            target_data = self.swiss_client.get_target_predictions(
                                str(row["Ligand SMILES"])
                            )
                            if target_data:
                                swiss_data["target_predictions"] = target_data[
                                    "predictions"
                                ]
                                swiss_data["target_predictions_url"] = target_data[
                                    "url"
                                ]

                            # Get ADME properties
                            adme_data = self.swiss_client.get_adme_properties(
                                str(row["Ligand SMILES"])
                            )
                            if adme_data:
                                swiss_data["adme_properties"] = adme_data

                            # Get similar compounds
                            similar_data = self.swiss_client.search_similar_compounds(
                                str(row["Ligand SMILES"])
                            )
                            if similar_data:
                                swiss_data["similar_compounds"] = similar_data[
                                    "similar_compounds"
                                ]
                                swiss_data["similar_compounds_url"] = similar_data[
                                    "url"
                                ]

                        compound_data = {
                            "name": clean_name,
                            "original_name": str(name),
                            "smiles": str(row["Ligand SMILES"]),
                            "inchi": str(row["Ligand InChI"]),
                            "inchi_key": str(row["Ligand InChI Key"]),
                            "affinity_type": best_affinity["type"],
                            "affinity_value": best_affinity["value"],
                            "affinity_unit": "nM",
                            "affinity_modifier": best_affinity["modifier"],
                            "assay_conditions": best_affinity["conditions"],
                            "activity_type": self._determine_activity_type(
                                str(row.get("Assay Description", ""))
                            ),
                            "target": str(row["Target Name"]),
                            "target_gene": str(
                                row.get(
                                    "UniProt (SwissProt) Primary ID of Target Chain",
                                    "N/A",
                                )
                            ),
                            "kinetics": {
                                "kon": str(row.get("kon (M-1-s-1)", "N/A")),
                                "koff": str(row.get("koff (s-1)", "N/A")),
                            },
                            "identifiers": {
                                "pubchem_cid": str(row["PubChem CID"]),
                                "pubchem_sid": str(row["PubChem SID"]),
                                "pubchem_aid": str(row.get("PubChem AID", "N/A")),
                                "chembl_id": str(row["ChEMBL ID of Ligand"]),
                                "chebi_id": str(row.get("ChEBI ID of Ligand", "N/A")),
                                "drugbank_id": str(
                                    row.get("DrugBank ID of Ligand", "N/A")
                                ),
                                "iuphar_id": str(
                                    row.get("IUPHAR_GRAC ID of Ligand", "N/A")
                                ),
                                "kegg_id": str(row.get("KEGG ID of Ligand", "N/A")),
                                "zinc_id": str(row.get("ZINC ID of Ligand", "N/A")),
                                "pdb_het": str(row.get("Ligand HET ID in PDB", "N/A")),
                                "pdb_complexes": str(
                                    row.get(
                                        "PDB ID(s) for Ligand-Target Complex", "N/A"
                                    )
                                ).split(","),
                            },
                            "references": {
                                "doi": str(row.get("Article DOI", "N/A")),
                                "pmid": str(row.get("PMID", "N/A")),
                                "authors": str(row.get("Authors", "N/A")),
                                "institution": str(row.get("Institution", "N/A")),
                                "bindingdb_ligand": str(
                                    row.get("Link to Ligand in BindingDB", "N/A")
                                ),
                                "bindingdb_target": str(
                                    row.get("Link to Target in BindingDB", "N/A")
                                ),
                                "bindingdb_pair": str(
                                    row.get(
                                        "Link to Ligand-Target Pair in BindingDB", "N/A"
                                    )
                                ),
                            },
                            "cas_number": cas_number,
                            "patent_number": str(row.get("Patent Number", "N/A")),
                            "reference_urls": self.web_client.get_reference_urls(
                                clean_name
                            ),
                            "matching_patterns": patterns if mol is not None else [],
                            "species": str(
                                row.get(
                                    "Target Source Organism According to Curator", "N/A"
                                )
                            ),
                            "data_source": str(row.get("Curation/DataSource", "N/A")),
                        }

                        # Get data from web sources with progress bar
                        from tqdm import tqdm

                        web_sources = [
                            (
                                "common_names",
                                lambda: self.web_client.get_common_names(
                                    clean_name,
                                    chembl_id=chembl_id,
                                    smiles=str(row["Ligand SMILES"]),
                                    inchi=str(row["Ligand InChI"]),
                                ),
                            ),
                            (
                                "legal_status",
                                lambda: self.web_client.get_legal_status(clean_name),
                            ),
                            (
                                "pharmacology",
                                lambda: self.web_client.get_pharmacology(clean_name),
                            ),
                        ]

                        web_progress = tqdm(
                            web_sources, desc="Fetching web data", unit="sources"
                        )

                        for source_name, source_func in web_progress:
                            web_progress.set_description(f"Fetching {source_name}")
                            if source_name == "common_names":
                                data = source_func()
                                compound_data["common_names"] = [
                                    name["name"] for name in data
                                ]
                            else:
                                compound_data[source_name] = source_func()

                        # Add patent information if API key provided
                        if llm_api_key:
                            self.logger.info("Searching patents...")
                            patent_info = self._search_patents(
                                clean_name, str(row["Ligand SMILES"]), llm_api_key
                            )
                            compound_data.update(patent_info)

                        # Add if valid structure
                        if self._validate_structure(compound_data):
                            compounds.append(compound_data)
                            self.logger.info(
                                f"Successfully processed compound: {clean_name}"
                            )
                        else:
                            self.logger.warning(
                                f"Invalid structure for compound: {clean_name}"
                            )

                    except Exception as e:
                        self.logger.error(
                            f"Error processing compound {clean_name}: {str(e)}"
                        )
                        continue

            # Add compounds from patents if API key provided
            if llm_api_key:
                self.logger.info("Searching patents for additional compounds...")
                # Check for patent compounds checkpoint
                patent_compounds = []
                if self.checkpoint_manager.is_step_completed("patent_compounds"):
                    self.logger.info("Loading patent compounds from checkpoint...")
                    patent_compounds = self.checkpoint_manager.load_step_data(
                        "patent_compounds"
                    )
                else:
                    patent_compounds = self._extract_compounds_from_patents(
                        "5-HT2 receptor ligand", llm_api_key
                    )
                    # Save patent compounds checkpoint
                    self.checkpoint_manager.save_checkpoint(
                        "patent_compounds", patent_compounds
                    )
                compounds.extend(patent_compounds)

            # Deduplicate by InChI Key
            self.logger.info("Deduplicating compounds by InChI Key...")
            unique_compounds = []
            seen_keys = set()
            for comp in compounds:
                if comp["inchi_key"] not in seen_keys:
                    seen_keys.add(comp["inchi_key"])
                    unique_compounds.append(comp)

            self.logger.info(f"Found {len(unique_compounds)} unique compounds")

            # Save final checkpoint
            self.checkpoint_manager.save_checkpoint("gather_ligands", unique_compounds)

            return unique_compounds

        except Exception as e:
            self.logger.error(f"Error gathering 5-HT2 ligands: {str(e)}")
            return []

    def save_5ht2_ligands(
        self, output_path: str, llm_api_key: Optional[str] = None
    ) -> None:
        """
        Save receptor ligands to TSV file.

        Args:
            output_path: Path to save TSV file
            llm_api_key: Optional API key for LLM-powered patent analysis
        """
        self.logger.info("Starting ligand collection process...")
        self.logger.info("Step 1/5: Gathering ligands from BindingDB...")
        compounds = self.gather_ligands(llm_api_key)
        if compounds:
            # Define TSV columns with enhanced fields
            columns = [
                # Basic identifiers
                "name",
                "chembl_id",
                "cas_number",
                "smiles",
                "inchi",
                "inchi_key",
                "common_name_1",
                "common_name_2",
                "common_name_3",
                # Target information
                "target",
                "target_gene",
                "target_protein",
                # Binding data
                "affinity_type",
                "affinity_value",
                "affinity_unit",
                "affinity_modifier",
                "activity_type",
                "kon",
                "koff",
                "pH",
                "temperature",
                # Bioassay data (up to 12)
            ]

            # Add bioassay columns
            for i in range(1, 13):
                prefix = f"bioassay_{i}_"
                columns.extend(
                    [
                        f"{prefix}aid",
                        f"{prefix}name",
                        f"{prefix}type",
                        f"{prefix}activity_type",
                        f"{prefix}target",
                        f"{prefix}value",
                        f"{prefix}unit",
                    ]
                )

            # Add Swiss tools columns
            columns.extend(
                [
                    # Target predictions (top 5)
                    "predicted_target_1",
                    "predicted_target_1_probability",
                    "predicted_target_2",
                    "predicted_target_2_probability",
                    "predicted_target_3",
                    "predicted_target_3_probability",
                    "predicted_target_4",
                    "predicted_target_4_probability",
                    "predicted_target_5",
                    "predicted_target_5_probability",
                    "target_predictions_url",
                    # ADME properties
                    "molecular_weight",
                    "logp",
                    "hbd",
                    "hba",
                    "tpsa",
                    "gi_absorption",
                    "bbb_permeant",
                    "pgp_substrate",
                    "cyp1a2_inhibitor",
                    "cyp2c19_inhibitor",
                    "cyp2c9_inhibitor",
                    "cyp2d6_inhibitor",
                    "cyp3a4_inhibitor",
                    "lipinski_violations",
                    "druglikeness_score",
                    # Similar compounds (top 3)
                    "similar_compound_1",
                    "similar_compound_1_similarity",
                    "similar_compound_2",
                    "similar_compound_2_similarity",
                    "similar_compound_3",
                    "similar_compound_3_similarity",
                    "similar_compounds_url",
                    # Add remaining columns
                    # Species and structure info
                    "species",
                    "pdb_het",
                    "pdb_complexes",
                    # Database IDs
                    "pubchem_cid",
                    "pubchem_sid",
                    "pubchem_aid",
                    "chebi_id",
                    "drugbank_id",
                    "iuphar_id",
                    "kegg_id",
                    "zinc_id",
                    # References
                    "doi",
                    "pmid",
                    "authors",
                    "institution",
                    # Additional data
                    "legal_status",
                    "mechanism_of_action",
                    "metabolism",
                    "toxicity",
                    "pharmacokinetics",
                    "drug_interactions",
                    "contraindications",
                    # URLs
                    "chembl_url",
                    "pubchem_url",
                    "wikipedia_url",
                    "psychonaut_url",
                    "erowid_url",
                    "bindingdb_ligand_url",
                    "bindingdb_target_url",
                    "bindingdb_pair_url",
                ]
            )

            self.logger.info("Step 2/5: Converting data to DataFrame...")
            # Convert to DataFrame for easier handling
            df = pd.DataFrame(compounds)
            total_compounds = len(df)
            self.logger.info(f"Found {total_compounds:,} unique compounds")

            self.logger.info("Step 3/5: Processing common names...")
            # Add common names columns with progress
            for i in range(3):
                self.logger.info(f"Processing common name column {i + 1}/3...")
                df[f"common_name_{i + 1}"] = df["common_names"].apply(
                    lambda x: x[i] if i < len(x) else ""
                )

            self.logger.info("Step 4/6: Processing additional data columns...")
            # Add legal status and pharmacology columns with progress bar
            from tqdm import tqdm

            progress_bar = tqdm(
                df.iterrows(),
                total=total_compounds,
                desc="Processing compound data",
                unit="compounds",
            )

            for idx, row in progress_bar:
                progress_bar.set_description(
                    f"Processing {row.get('name', '')[:30]}..."
                )

                # Legal status
                df.at[idx, "legal_status"] = (
                    "; ".join(
                        f"{s['jurisdiction']}: {s['schedule']}"
                        for s in row["legal_status"]["scheduling"]
                    )
                    if isinstance(row.get("legal_status"), dict)
                    else ""
                )

                # Pharmacology
                if isinstance(row.get("pharmacology"), dict):
                    df.at[idx, "mechanism_of_action"] = "; ".join(
                        row["pharmacology"].get("mechanism_of_action", [])
                    )
                    df.at[idx, "metabolism"] = "; ".join(
                        row["pharmacology"].get("metabolism", [])
                    )
                    df.at[idx, "toxicity"] = "; ".join(
                        row["pharmacology"].get("toxicity", [])
                    )

            self.logger.info("Step 5/6: Adding Swiss tools data...")
            # Add Swiss tools data with progress bar
            swiss_progress = tqdm(
                df.iterrows(),
                total=total_compounds,
                desc="Processing Swiss tools data",
                unit="compounds",
            )

            for idx, row in swiss_progress:
                swiss_progress.set_description(
                    f"Processing {row.get('name', '')[:30]}..."
                )

                # Add target predictions
                predictions = row.get("swiss_data", {}).get("target_predictions", [])
                for i, pred in enumerate(predictions[:5], 1):
                    df.at[idx, f"predicted_target_{i}"] = pred.get("target", "")
                    df.at[idx, f"predicted_target_{i}_probability"] = pred.get(
                        "probability", 0.0
                    )
                df.at[idx, "target_predictions_url"] = row.get("swiss_data", {}).get(
                    "target_predictions_url", ""
                )

                # Add ADME properties
                adme = row.get("swiss_data", {}).get("adme_properties", {})
                physchem = adme.get("physicochemical", {})
                df.at[idx, "molecular_weight"] = physchem.get("mw", "")
                df.at[idx, "logp"] = physchem.get("logp", "")
                df.at[idx, "hbd"] = physchem.get("hbd", "")
                df.at[idx, "hba"] = physchem.get("hba", "")
                df.at[idx, "tpsa"] = physchem.get("tpsa", "")

                absorption = adme.get("absorption", {})
                df.at[idx, "gi_absorption"] = absorption.get("gi_absorption", "")
                df.at[idx, "bbb_permeant"] = absorption.get("bbb_permeant", "")
                df.at[idx, "pgp_substrate"] = absorption.get("pgp_substrate", "")

                metabolism = adme.get("metabolism", {}).get("cyp_inhibition", {})
                df.at[idx, "cyp1a2_inhibitor"] = metabolism.get("cyp1a2", "")
                df.at[idx, "cyp2c19_inhibitor"] = metabolism.get("cyp2c19", "")
                df.at[idx, "cyp2c9_inhibitor"] = metabolism.get("cyp2c9", "")
                df.at[idx, "cyp2d6_inhibitor"] = metabolism.get("cyp2d6", "")
                df.at[idx, "cyp3a4_inhibitor"] = metabolism.get("cyp3a4", "")

                druglike = adme.get("druglikeness", {})
                df.at[idx, "lipinski_violations"] = druglike.get(
                    "Lipinski Violations", ""
                )
                df.at[idx, "druglikeness_score"] = druglike.get(
                    "Druglikeness Score", ""
                )

                # Add similar compounds
                similars = row.get("swiss_data", {}).get("similar_compounds", [])
                for i, sim in enumerate(similars[:3], 1):
                    df.at[idx, f"similar_compound_{i}"] = sim.get("name", "")
                    df.at[idx, f"similar_compound_{i}_similarity"] = sim.get(
                        "similarity", 0.0
                    )
                df.at[idx, "similar_compounds_url"] = row.get("swiss_data", {}).get(
                    "similar_compounds_url", ""
                )

            self.logger.info("Step 6/6: Adding reference URLs...")
            # Add reference URLs with progress bar
            from tqdm import tqdm

            url_types = ["chembl", "pubchem", "wikipedia", "psychonaut", "erowid"]
            url_progress = tqdm(url_types, desc="Processing URL types", unit="type")

            for url_type in url_progress:
                url_progress.set_description(f"Processing {url_type} URLs")
                df[f"{url_type}_url"] = df.apply(
                    lambda row: (
                        row["reference_urls"].get(f"{url_type}_url", "")
                        if isinstance(row["reference_urls"], dict)
                        else ""
                    ),
                    axis=1,
                )

            # Save as TSV
            df[columns].to_csv(output_path, sep="\t", index=False)

            # Log statistics
            self.logger.info(
                f"Saved {len(compounds)} 5-HT2 receptor ligands to {output_path}"
            )

            # Collect statistics with progress bar
            from tqdm import tqdm

            sources = set()
            with_cas = 0
            with_patents = 0
            activity_types = {}
            species_counts = {}

            stats_progress = tqdm(
                compounds, desc="Collecting statistics", unit="compounds"
            )

            for comp in stats_progress:
                stats_progress.set_description(
                    f"Analyzing {comp.get('name', '')[:30]}..."
                )

                if "source" in comp:
                    sources.add(comp["source"])
                if comp.get("cas_number", "N/A") != "N/A":
                    with_cas += 1
                if comp.get("patents", []):
                    with_patents += 1

                activity = comp.get("activity_type", "unknown")
                activity_types[activity] = activity_types.get(activity, 0) + 1

                species = comp.get("species", "N/A")
                species_counts[species] = species_counts.get(species, 0) + 1

            self.logger.info(
                f"Data sources used: {', '.join(sources)}\n"
                f"Compounds with CAS numbers: {with_cas}\n"
                f"Compounds with patent data: {with_patents}\n"
                f"Activity type distribution:\n"
                + "\n".join(f"  {k}: {v}" for k, v in activity_types.items())
                + "\nSpecies distribution:\n"
                + "\n".join(f"  {k}: {v}" for k, v in species_counts.items())
            )
        else:
            self.logger.warning("No 5-HT2 receptor ligands found")

    def _validate_structure(self, compound: Dict[str, Any]) -> bool:
        """
        Validate chemical structure using RDKit.

        Args:
            compound: Compound dictionary containing SMILES

        Returns:
            True if structure is valid, False otherwise
        """
        if not compound.get("smiles"):
            return False

        try:
            mol = Chem.MolFromSmiles(str(compound["smiles"]))
            if mol is None:
                return False

            # Generate 3D conformation
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)

            return True
        except Exception as e:
            self.logger.error(f"Error validating structure: {str(e)}")
            return False

    def _clean_name(self, name: str) -> str:
        """
        Clean compound name for web searches.

        Args:
            name: Raw compound name

        Returns:
            Cleaned compound name
        """
        # Remove ChEMBL ID
        name = re.sub(r"::CHEMBL\d+$", "", name)

        # Remove stereochemistry
        name = re.sub(r"\([^)]*\)", "", name)
        name = re.sub(r"[RS]-", "", name)
        name = re.sub(r"\b[RS]\b", "", name)

        # Remove special characters but keep some important ones
        name = re.sub(r"[^\w\s\-\+]", "", name)

        # Normalize whitespace
        name = " ".join(name.split())

        # Remove common prefixes/suffixes
        name = re.sub(r"^(?:\+\/\-|\+|\-)\s*", "", name)
        name = re.sub(
            r"\s*(?:hydrochloride|HCl|salt|hydrate|solvate)$", "", name, flags=re.I
        )

        return name

    def _get_cas_number(self, name: str, smiles: str, mol: Optional[Chem.Mol]) -> str:
        """
        Get CAS number from various sources.

        Args:
            name: Compound name
            smiles: SMILES string
            mol: RDKit molecule object

        Returns:
            CAS number or 'N/A' if not found
        """
        try:
            # Try PubChem first
            if pd.notna(smiles) and mol is not None:
                # Use canonical SMILES
                canonical_smiles = self.structure_utils.standardize_smiles(smiles)
                if canonical_smiles:
                    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{canonical_smiles}/synonyms/JSON"
                    response = self.web_client._make_request(url)
                    if response:
                        data = response.json()
                        if "InformationList" in data:
                            synonyms = data["InformationList"]["Information"][0][
                                "Synonym"
                            ]
                            for syn in synonyms:
                                # Look for CAS pattern
                                if re.match(r"^\d{1,7}-\d{2}-\d$", syn):
                                    return syn

            # Try web search
            search_results = self.web_client.get_common_names(name)
            for result in search_results:
                if re.match(r"^\d{1,7}-\d{2}-\d$", result["name"]):
                    return result["name"]

        except Exception as e:
            self.logger.error(f"Error getting CAS number: {str(e)}")

        return "N/A"

    def _search_patents(self, name: str, smiles: str, api_key: str) -> Dict[str, Any]:
        """
        Search patents for compound information.

        Args:
            name: Compound name
            smiles: SMILES string
            api_key: API key for LLM service

        Returns:
            Dictionary containing:
            - patents: List of relevant patents
            - patent_count: Total number of patents found
        """
        try:
            # Use Google Patents API
            url = f"https://patents.google.com/api/query?q={name}"
            headers = {"Authorization": f"Bearer {api_key}"}
            response = requests.get(url, headers=headers)

            if response.status_code == 200:
                data = response.json()
                patents = []

                for result in data.get("results", []):
                    patent = {
                        "number": result.get("patent_number"),
                        "title": result.get("title"),
                        "abstract": result.get("abstract"),
                        "url": f"https://patents.google.com/patent/{result['patent_number']}",
                    }
                    patents.append(patent)

                return {
                    "patents": patents[:5],  # Top 5 most relevant patents
                    "patent_count": len(data.get("results", [])),
                }

        except Exception as e:
            self.logger.error(f"Error searching patents: {str(e)}")

        return {"patents": [], "patent_count": 0}

    def _extract_compounds_from_patents(
        self, search_term: str, api_key: str
    ) -> List[Dict[str, Any]]:
        """
        Extract novel compounds from patents.

        Args:
            search_term: Search term for patents
            api_key: API key for LLM service

        Returns:
            List of compound dictionaries
        """
        compounds = []
        try:
            # Search patents
            url = f"https://patents.google.com/api/query?q={search_term}"
            headers = {"Authorization": f"Bearer {api_key}"}
            response = requests.get(url, headers=headers)

            if response.status_code == 200:
                data = response.json()

                for result in data.get("results", []):
                    # Extract text and chemical structures
                    text = f"{result.get('title', '')} {result.get('abstract', '')}"

                    # Use LLM to extract chemical information
                    extracted = self._extract_chemical_info(text, api_key)

                    for compound in extracted:
                        if "smiles" in compound:
                            # Validate structure
                            mol = Chem.MolFromSmiles(compound["smiles"])
                            if mol is not None:
                                compounds.append(
                                    {
                                        "name": compound.get("name", "Unknown"),
                                        "smiles": compound["smiles"],
                                        "inchi": Chem.MolToInchi(mol),
                                        "inchi_key": Chem.MolToInchiKey(mol),
                                        "source": "Patent",
                                        "patent_number": result.get("patent_number"),
                                        "patent_url": f"https://patents.google.com/patent/{result['patent_number']}",
                                    }
                                )

        except Exception as e:
            self.logger.error(f"Error extracting compounds from patents: {str(e)}")

        return compounds

    def _extract_chemical_info(self, text: str, api_key: str) -> List[Dict[str, Any]]:
        """
        Use LLM to extract chemical information from text.

        Args:
            text: Text to analyze
            api_key: API key for LLM service

        Returns:
            List of dictionaries containing extracted chemical information
        """
        try:
            # Prepare prompt with enhanced instructions
            prompt = f"""Extract chemical compounds and their properties from the following text. For each compound, provide:

1. Chemical names (IUPAC, common names, and any synonyms)
2. Chemical structure information:
   - SMILES string
   - InChI string (if available)
   - Molecular formula
   - Structural features (rings, functional groups)
3. Biological activity:
   - Target receptors/proteins
   - Activity type (agonist, antagonist, etc.)
   - Binding affinity values (Ki, IC50, etc.)
4. Physical properties:
   - Molecular weight
   - Melting/boiling points
   - Solubility
5. Patent-specific information:
   - Example numbers
   - Synthesis methods
   - Claims coverage

Text: {text}

Format the response as a list of JSON objects with these fields. Include any numerical values with their units. For chemical structures, prioritize standardized identifiers (SMILES, InChI) over descriptive text."""

            # Make API request
            response = requests.post(
                "https://api.openai.com/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {api_key}",
                    "Content-Type": "application/json",
                },
                json={
                    "model": "gpt-4",
                    "messages": [{"role": "user", "content": prompt}],
                    "temperature": 0.7,
                },
            )

            if response.status_code == 200:
                result = response.json()
                try:
                    # Parse the response as JSON
                    compounds = json.loads(result["choices"][0]["message"]["content"])
                    return compounds
                except json.JSONDecodeError:
                    return []

        except Exception as e:
            self.logger.error(f"Error extracting chemical info: {str(e)}")

        return []
