"""Binding data processing and BindingDB integration.

This module handles:
1. Loading and processing binding data from BindingDB
2. Gathering 5-HT2 receptor ligands from multiple sources
3. Enriching compound data with web sources
4. Patent searching and analysis
5. Structure validation and standardization
6. Swiss tools integration for target prediction and ADME properties
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
        
        # Initialize Swiss tools client
        from web_enrichment.data_sources.swiss import SwissClient
        self.swiss_client = SwissClient(http_client) if http_client else None
        
        # Path to downloaded BindingDB TSV file
        self.bindingdb_path = os.path.join(
            os.path.dirname(__file__),
            "BindingDB_All.tsv"
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
            zip_path = os.path.join(os.path.dirname(__file__), "BindingDB_All_202501_tsv.zip")
            
            try:
                # Download zip file with progress bar
                response = requests.get(url, stream=True)
                response.raise_for_status()
                total_size = int(response.headers.get("content-length", 0))
                
                progress = tqdm(
                    total=total_size,
                    unit="iB",
                    unit_scale=True,
                    desc="Downloading BindingDB data"
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
                        files,
                        desc="Extracting files",
                        unit="files"
                    )
                    
                    for file in extract_progress:
                        zip_ref.extract(file, os.path.dirname(__file__))
                        extract_progress.set_description(f"Extracted {file}")
                
                # Remove zip file
                os.remove(zip_path)
                self.logger.info("BindingDB TSV file downloaded and extracted successfully.")
                
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
        self,
        compound_name: str,
        smiles: Optional[str] = None
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
