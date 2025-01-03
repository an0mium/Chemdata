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

from logger import LogManager
from models import CompoundData
from structure_utils import StructureUtils
from web_enrichment import WebEnrichment


class BindingDataProcessor:
    """Handles binding data processing and BindingDB integration."""
    
    # Target patterns for receptor ligands
    TARGET_PATTERNS = [
        # 5-HT2 receptor patterns
        '5-HT2', 'serotonin 2', 'HTR2',
        '5-hydroxytryptamine receptor 2',
        'serotonin receptor type 2',
        'serotonin receptor subtype 2',
        '5-HT2A', '5-HT2B', '5-HT2C',
        'HTR2A', 'HTR2B', 'HTR2C',
        'serotonin 2A', 'serotonin 2B', 'serotonin 2C',
        'HT2A_HUMAN', 'HT2B_HUMAN', 'HT2C_HUMAN',
        '5-hydroxytryptamine2A', '5-hydroxytryptamine2B', '5-hydroxytryptamine2C',
        'serotonin receptor 2A', 'serotonin receptor 2B', 'serotonin receptor 2C',
        '5-hydroxytryptamine receptor 2A', '5-hydroxytryptamine receptor 2B', '5-hydroxytryptamine receptor 2C',
        '5HT2A', '5HT2B', '5HT2C', '5-HT-2A', '5-HT-2B', '5-HT-2C',
        
        # NMDA receptor patterns
        'NMDA', 'N-methyl-D-aspartate',
        'GRIN1', 'GRIN2A', 'GRIN2B', 'GRIN2C', 'GRIN2D',
        'GluN1', 'GluN2A', 'GluN2B', 'GluN2C', 'GluN2D',
        'NR1', 'NR2A', 'NR2B', 'NR2C', 'NR2D',
        'NMDAR1', 'NMDAR2A', 'NMDAR2B', 'NMDAR2C', 'NMDAR2D',
        'glutamate receptor ionotropic NMDA',
        'glutamate [NMDA] receptor',
        
        # NMDA antagonist patterns
        '3-HO-PCP', '3-MeO-PCP', '3-Methyl-PCPy', '4-MeO-PCP',
        'ACE mixture', 'Agmatine', 'Alaproclate', 'Alazocine',
        'Amantadine', 'AP-7', 'AP5', 'Apigenin', 'Aptiganel',
        'Arketamine', 'Atomoxetine', 'Besonprodil', 'Budipine',
        'Bumetanide', 'Buphenine', 'Carisoprodol', 'Caroverine',
        'CGP-37849', 'CGP-39551', '4-Chlorokynurenine', 'CNQX',
        'Conantokin', 'Coronaridine', 'Crocetin', 'Cyclopropane',
        'Delucemine', 'Deschloroketamine', 'Dextrallorphan',
        'Dextromethorphan', 'Dextropropoxyphene', 'Dextrorphan',
        '1,3-Diaminopropane', '5,7-Dichlorokynurenic acid',
        'Diethyl ether', 'Diethylenetriamine', 'Dieticyclidine',
        'Diphenidine', 'Dizocilpine', 'DNQX', 'Eliprodil',
        'Α-Endopsychosin', 'Enflurane', 'Ephenidine', 'Esketamine',
        'Esmethadone', 'NEFA', 'Eticyclidine', 'EVT-101', 'EVT-103',
        'Felbamate', 'Flufenamic acid', '2-Fluorodeschloroketamine',
        'Fluorolintane', 'Flupirtine', 'Fourphit', 'Furosemide',
        'Gacyclidine', 'Gavestinel', 'HA-966', 'Haloperidol',
        'Halothane', 'Hemantane', 'Hodgkinsine', 'Huperzine A',
        'Hydroxynorketamine', 'Ibogaine', 'Ibogamine', 'Ifenprodil',
        'Indantadol', 'Indeloxazine', 'Isoflurane', 'Isoxsuprine',
        'Kaitocephalin', 'Ketamine', 'Ketobemidone', 'Ketofol',
        'Kynurenic acid', 'Kynurenine', 'L-701324', 'Lanicemine',
        'Levomethadone', 'Levomethorphan', 'Levomilnacipran',
        'Levorphanol', 'Licostinel', 'Lubeluzole', 'LY-235959',
        'Memantine', 'Meprobamate', 'Metaphit', 'Methoxetamine',
        'Methoxphenidine', '18-Methoxycoronaridine', 'Methoxyflurane',
        'Midafotel', 'Milnacipran', 'Minocycline', 'Nelonemdaz',
        'Neramexane', 'Niflumic acid', 'Nitromemantine', 'Nitrous oxide',
        'Noribogaine', 'Norketamine', 'Nortilidine', 'NPDPA',
        'Onfasprodil', 'Orphenadrine', 'PCPr', 'PD-137889', 'PEAQX',
        'Pentamidine', 'Perzinfotel', 'Pethidine', 'Phencyclidine',
        '8A-PDHQ', 'Piretanide', 'Promethazine', 'Psychotridine',
        'Putrescine', 'Racemorphan', 'Ralfinamide', 'Remacemide',
        'Rhynchophylline', 'Rislenemdaz', 'Rolicyclidine', 'Sabeluzole',
        'Selfotel', 'Sevoflurane', 'SN 35210', 'Spasmolytic A29',
        'Tabernanthine', 'Tenocyclidine', 'Tiletamine', 'Tramadol',
        'Traxoprodil', '2,2,2-Trichloroethanol', 'Trichloroethylene',
        'Xenon', 'XW10508', 'ZD-9379',
        
        # Nootropic targets
        'acetylcholinesterase', 'AChE',
        'nicotinic acetylcholine receptor', 'nAChR',
        'muscarinic acetylcholine receptor', 'mAChR',
        'AMPA receptor', 'GRIA', 'GluA',
        'dopamine transporter', 'DAT', 'SLC6A3',
        'norepinephrine transporter', 'NET', 'SLC6A2',
        
        # Analgesic targets
        'mu opioid receptor', 'OPRM1', 'MOR',
        'delta opioid receptor', 'OPRD1', 'DOR',
        'kappa opioid receptor', 'OPRK1', 'KOR',
        'cannabinoid receptor', 'CNR1', 'CNR2', 'CB1', 'CB2',
        'cyclooxygenase', 'COX-1', 'COX-2', 'PTGS1', 'PTGS2',
        
        # Antidepressant targets
        'serotonin transporter', 'SERT', 'SLC6A4',
        'noradrenaline transporter', 'NET', 'SLC6A2',
        'monoamine oxidase', 'MAO-A', 'MAO-B',
        
        # Anxiolytic targets
        'GABA receptor', 'GABRA', 'GABRB', 'GABRG',
        'benzodiazepine receptor', 'TSPO',
        '5-HT1A', 'HTR1A', 'serotonin 1A',
        
        # Stimulant targets
        'dopamine transporter', 'DAT', 'SLC6A3',
        'norepinephrine transporter', 'NET', 'SLC6A2',
        'trace amine receptor', 'TAAR1',
        
        # Psychedelic targets
        '5-HT2A', 'HTR2A', 'serotonin 2A',
        '5-HT1A', 'HTR1A', 'serotonin 1A',
        'sigma receptor', 'SIGMAR1', 'sigma-1',
        
        # Dissociative targets
        'NMDA receptor', 'glutamate [NMDA] receptor',
        'sigma receptor', 'SIGMAR1', 'sigma-1',
        'kappa opioid receptor', 'OPRK1', 'KOR',
        
        # Entactogen targets
        'serotonin transporter', 'SERT', 'SLC6A4',
        'vesicular monoamine transporter', 'VMAT2', 'SLC18A2',
        '5-HT2A', 'HTR2A', 'serotonin 2A',
        '5-HT1A', 'HTR1A', 'serotonin 1A'
    ]
    
    # Activity type patterns with more detail
    ACTIVITY_PATTERNS = {
        'superagonist': [
            r'super.?agonist',
            r'high.?efficacy.?agonist',
            r'full.?agonist.+high.?efficacy',
            r'efficacy\s*>\s*100%',
            r'super.?potent.?agonist'
        ],
        'full_agonist': [
            r'full.?agonist',
            r'complete.?agonist',
            r'full.?receptor.?activation',
            r'efficacy\s*[~≈≃]\s*100%',
            r'maximal.?response'
        ],
        'partial_agonist': [
            r'partial.?agonist',
            r'submaximal.?activation',
            r'partial.?receptor.?activation',
            r'efficacy\s*[<≈]\s*\d{1,2}%',
            r'partial.?response'
        ],
        'weak_partial_agonist': [
            r'weak.?partial.?agonist',
            r'low.?efficacy.?partial',
            r'weak.?partial.?activation',
            r'efficacy\s*<\s*20%',
            r'minimal.?agonist'
        ],
        'mixed_agonist_antagonist': [
            r'mixed.?agonist.?antagonist',
            r'partial.?agonist.?antagonist',
            r'dual.?activity',
            r'context.?dependent',
            r'tissue.?dependent'
        ],
        'antagonist': [
            r'antagonist',
            r'blocker',
            r'inhibitor',
            r'neutral.?antagonist',
            r'competitive.?antagonist'
        ],
        'inverse_agonist': [
            r'inverse.?agonist',
            r'negative.?agonist',
            r'inverse.?activity',
            r'negative.?efficacy',
            r'constitutive.?inhibitor'
        ],
        'positive_allosteric_modulator': [
            r'positive.?allosteric',
            r'PAM',
            r'positive.?modulator',
            r'allosteric.?potentiator',
            r'positive.?cooperativity'
        ],
        'negative_allosteric_modulator': [
            r'negative.?allosteric',
            r'NAM',
            r'negative.?modulator',
            r'allosteric.?inhibitor',
            r'negative.?cooperativity'
        ],
        'complex_modulator': [
            r'complex.?modulator',
            r'mixed.?modulator',
            r'complex.?pharmacology',
            r'bitopic',
            r'dual.?mechanism'
        ]
    }
    
    def __init__(self, pubmed_client, web_client=None):
        """
        Initialize binding data processor.
        
        Args:
            pubmed_client: PubMed client for relevance scoring
            web_client: Optional WebEnrichment client for web scraping
        """
        self.logger = LogManager().get_logger("binding_data_processor")
        self.pubmed_client = pubmed_client
        self.web_client = web_client or WebEnrichment()
        self.structure_utils = StructureUtils()
        
        # Path to downloaded BindingDB TSV file
        self.bindingdb_path = os.path.join(
            os.path.dirname(__file__),
            'BindingDB_All.tsv'
        )

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
        return 'unknown'

    def _ensure_bindingdb_file(self) -> None:
        """Download and extract BindingDB TSV file if it doesn't exist."""
        if not os.path.exists(self.bindingdb_path):
            self.logger.info("BindingDB TSV file not found. Downloading...")
            url = "https://bindingdb.org/bind/downloads/BindingDB_All_202501_tsv.zip"
            zip_path = os.path.join(os.path.dirname(__file__), 'BindingDB_All_202501_tsv.zip')
            
            # Download zip file
            response = requests.get(url, stream=True)
            response.raise_for_status()
            with open(zip_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            # Extract TSV file
            import zipfile
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(os.path.dirname(__file__))
            
            # Remove zip file
            os.remove(zip_path)
            self.logger.info("BindingDB TSV file downloaded and extracted.")

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
                'BindingDB Ligand Name',
                'Ligand SMILES',
                'Ligand InChI',
                'Ligand InChI Key',
                'Target Name',
                'Target Source Organism According to Curator or DataSource',
                'Ki (nM)',
                'IC50 (nM)',
                'Kd (nM)',
                'EC50 (nM)',
                'kon (M-1-s-1)',
                'koff (s-1)',
                'pH',
                'Temp (C)',
                'Curation/DataSource',
                'Article DOI',
                'BindingDB Entry DOI',
                'PMID',
                'PubChem AID',
                'Patent Number',
                'Authors',
                'Institution',
                'Link to Ligand in BindingDB',
                'Link to Target in BindingDB',
                'Link to Ligand-Target Pair in BindingDB',
                'Ligand HET ID in PDB',
                'PDB ID(s) for Ligand-Target Complex',
                'PubChem CID',
                'PubChem SID',
                'ChEBI ID of Ligand',
                'ChEMBL ID of Ligand',
                'DrugBank ID of Ligand',
                'IUPHAR_GRAC ID of Ligand',
                'KEGG ID of Ligand',
                'ZINC ID of Ligand',
                'Number of Protein Chains in Target (>1 implies a multichain complex)',
                'UniProt (SwissProt) Recommended Name of Target Chain',
                'UniProt (SwissProt) Entry Name of Target Chain',
                'UniProt (SwissProt) Primary ID of Target Chain',
                'UniProt (SwissProt) Secondary ID(s) of Target Chain',
                'UniProt (SwissProt) Alternative ID(s) of Target Chain',
                'UniProt (TrEMBL) Submitted Name of Target Chain',
                'UniProt (TrEMBL) Entry Name of Target Chain',
                'UniProt (TrEMBL) Primary ID of Target Chain',
                'UniProt (TrEMBL) Secondary ID(s) of Target Chain',
                'UniProt (TrEMBL) Alternative ID(s) of Target Chain'
            ]

            
            # Download BindingDB TSV file if needed
            self._ensure_bindingdb_file()
            
            self.logger.info("Reading BindingDB TSV file (this may take a few minutes)...")
            total_lines = sum(1 for _ in open(self.bindingdb_path))
            self.logger.info(f"Total lines in TSV file: {total_lines:,}")
            df = pd.read_csv(
                self.bindingdb_path,
                sep='\t',
                usecols=needed_columns,
                on_bad_lines='skip',
                low_memory=False
            )
            self.logger.info("BindingDB TSV file loaded successfully.")
            
            # Search for compound by name and SMILES
            mask = df['BindingDB Ligand Name'].str.contains(
                str(compound_name), case=False, na=False
            )
            if smiles:
                mask |= df['Ligand SMILES'] == str(smiles)
                
            matches = df[mask]
            
            binding_data = []
            for _, row in matches.iterrows():
                # Skip non-mammalian targets
                organism = str(row.get('Target Source Organism According to Curator or DataSource', '')).lower()
                if not any(x in organism for x in ['human', 'mouse', 'rat', 'mammal']):
                    self.logger.debug(f"Skipping non-mammalian target: {organism}")
                    continue
                    
                data = {
                    'target_common_name': str(row['Target Name']),
                    'target_protein_name': str(row.get(
                        'UniProt (SwissProt) Recommended Name of Target Chain',
                        'N/A'
                    )),
                    'target_gene_name': str(row.get('UniProt (SwissProt) Primary ID of Target Chain', 'N/A')),
                    'affinity_value': None,
                    'affinity_unit': 'nM',
                    'affinity_type': None,
                    'activity_type': self._determine_activity_type(
                        str(row.get('Assay Description', ''))
                    ),
                    'source': 'BindingDB',
                    'doi': str(row.get('Article DOI', 'N/A')),
                    'pmid': str(row.get('PMID', 'N/A'))
                }
                
                # Get strongest affinity value
                for aff_type in ['Ki', 'IC50', 'Kd', 'EC50']:
                    col = f'{aff_type} (nM)'
                    if pd.notna(row[col]):
                        value = str(row[col])
                        # Handle '>' and '<' in affinity values
                        if value.startswith('>'):
                            # For '>' values, use the number * 10 to sort them after exact values
                            data['affinity_value'] = float(value[1:]) * 10
                            data['affinity_modifier'] = '>'
                        elif value.startswith('<'):
                            # For '<' values, use the number / 10 to sort them before exact values
                            data['affinity_value'] = float(value[1:]) / 10
                            data['affinity_modifier'] = '<'
                        else:
                            data['affinity_value'] = float(value)
                            data['affinity_modifier'] = ''
                        data['affinity_type'] = aff_type
                        break
                    
                if data['affinity_value'] is not None:
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
            target_name = getattr(compound, f'target_{i}_common_name', 'N/A')
            if target_name != "N/A":
                data = {
                    'target_common_name': target_name,
                    'target_protein_name': getattr(
                        compound, f'target_{i}_protein_name', 'N/A'
                    ),
                    'target_gene_name': getattr(
                        compound, f'target_{i}_gene_name', 'N/A'
                    ),
                    'affinity_value': getattr(
                        compound, f'target_{i}_affinity', 'N/A'
                    ),
                    'affinity_unit': getattr(
                        compound, f'target_{i}_affinity_unit', 'N/A'
                    ),
                    'affinity_type': getattr(
                        compound, f'target_{i}_affinity_type', 'N/A'
                    ),
                    'activity_type': getattr(
                        compound, f'target_{i}_activity_type', 'unknown'
                    ),
                    'source': getattr(compound, f'target_{i}_source', 'N/A'),
                    'doi': getattr(compound, f'target_{i}_doi', 'N/A'),
                    'pmid': getattr(compound, f'target_{i}_pmid', 'N/A')
                }
                binding_data.append(data)
                
        # Get PubMed relevance scores for each target
        scored_data = []
        for data in binding_data:
            pubmed_results = self.pubmed_client.get_binding_relevance(
                compound.name,
                data['target_common_name']
            )
            scored_data.append((data, pubmed_results))
            
        # Sort by PubMed results
        scored_data.sort(key=lambda x: x[1], reverse=True)
        
        # Update compound with sorted binding data (up to 12 targets)
        for i, (data, _) in enumerate(scored_data[:12], 1):
            setattr(compound, f'target_{i}_common_name', data['target_common_name'])
            setattr(compound, f'target_{i}_protein_name', data['target_protein_name'])
            setattr(compound, f'target_{i}_gene_name', data['target_gene_name'])
            setattr(compound, f'target_{i}_affinity', str(data['affinity_value']))
            setattr(compound, f'target_{i}_affinity_unit', data['affinity_unit'])
            setattr(compound, f'target_{i}_affinity_type', data['affinity_type'])
            setattr(compound, f'target_{i}_activity_type', data['activity_type'])
            setattr(compound, f'target_{i}_source', data['source'])
            setattr(compound, f'target_{i}_doi', data['doi'])
            setattr(compound, f'target_{i}_pmid', data['pmid'])



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
            # Read BindingDB data with enhanced columns
            needed_columns = [
                'BindingDB Ligand Name',
                'Ligand SMILES',
                'Ligand InChI',
                'Ligand InChI Key',
                'Target Name',
                'Target Source Organism According to Curator or DataSource',
                'Ki (nM)',
                'IC50 (nM)',
                'Kd (nM)',
                'EC50 (nM)',
                'kon (M-1-s-1)',
                'koff (s-1)',
                'pH',
                'Temp (C)',
                'Curation/DataSource',
                'Article DOI',
                'BindingDB Entry DOI',
                'PMID',
                'PubChem AID',
                'Patent Number',
                'Authors',
                'Institution',
                'Link to Ligand in BindingDB',
                'Link to Target in BindingDB',
                'Link to Ligand-Target Pair in BindingDB',
                'Ligand HET ID in PDB',
                'PDB ID(s) for Ligand-Target Complex',
                'PubChem CID',
                'PubChem SID',
                'ChEBI ID of Ligand',
                'ChEMBL ID of Ligand',
                'DrugBank ID of Ligand',
                'IUPHAR_GRAC ID of Ligand',
                'KEGG ID of Ligand',
                'ZINC ID of Ligand',
                'Number of Protein Chains in Target (>1 implies a multichain complex)',
                'UniProt (SwissProt) Recommended Name of Target Chain',
                'UniProt (SwissProt) Entry Name of Target Chain',
                'UniProt (SwissProt) Primary ID of Target Chain',
                'UniProt (SwissProt) Secondary ID(s) of Target Chain',
                'UniProt (SwissProt) Alternative ID(s) of Target Chain',
                'UniProt (TrEMBL) Submitted Name of Target Chain',
                'UniProt (TrEMBL) Entry Name of Target Chain',
                'UniProt (TrEMBL) Primary ID of Target Chain',
                'UniProt (TrEMBL) Secondary ID(s) of Target Chain',
                'UniProt (TrEMBL) Alternative ID(s) of Target Chain'
            ]

            
            # Download BindingDB TSV file if needed
            self._ensure_bindingdb_file()
            
            self.logger.info("Reading BindingDB TSV file (this may take a few minutes)...")
            total_lines = sum(1 for _ in open(self.bindingdb_path))
            self.logger.info(f"Total lines in TSV file: {total_lines:,}")
            df = pd.read_csv(
                self.bindingdb_path,
                sep='\t',
                usecols=needed_columns,
                on_bad_lines='skip',
                low_memory=False
            )
            self.logger.info("BindingDB TSV file loaded successfully.")
            
            # Find 5-HT2 receptor targets
            self.logger.info("Searching for 5-HT2 receptor targets in BindingDB data...")
            pattern = '|'.join(self.TARGET_PATTERNS)
            
            # Convert columns to string and handle NaN values
            self.logger.info("Converting columns to string format...")
            df['Target Name'] = df['Target Name'].fillna('').astype(str)
            df['UniProt (SwissProt) Recommended Name of Target Chain'] = df['UniProt (SwissProt) Recommended Name of Target Chain'].fillna('').astype(str)
            df['UniProt (SwissProt) Primary ID of Target Chain'] = df['UniProt (SwissProt) Primary ID of Target Chain'].fillna('').astype(str)
            
            # Search in all relevant columns
            self.logger.info("Searching for 5-HT2 patterns in target columns...")
            mask = df['Target Name'].str.contains(pattern, case=False, na=False)
            mask |= df['UniProt (SwissProt) Recommended Name of Target Chain'].str.contains(pattern, case=False, na=False)
            mask |= df['UniProt (SwissProt) Primary ID of Target Chain'].str.contains(pattern, case=False, na=False)
            matches = df[mask]
            
            self.logger.info(f"Found {len(matches):,} initial matches in BindingDB")
            self.logger.info("Filtering matches by organism and structure...")
            
            # Group by compound and combine data
            compounds = []
            seen_inchikeys = set()
            processed_count = 0
            total_compounds = len(matches.groupby('BindingDB Ligand Name'))
            
            for name, group in matches.groupby('BindingDB Ligand Name'):
                processed_count += 1
                if processed_count % 10 == 0:
                    self.logger.info(f"Processing compound {processed_count}/{total_compounds}")
                
                # Skip if we've seen this compound
                inchikey = str(group.iloc[0]['Ligand InChI Key'])
                if pd.notna(inchikey) and inchikey in seen_inchikeys:
                    continue
                seen_inchikeys.add(inchikey)
                
                # Get all affinity values
                affinities = []
                for _, row in group.iterrows():
                    for aff_type in ['Ki', 'IC50', 'Kd', 'EC50']:
                        col = f'{aff_type} (nM)'
                        if pd.notna(row[col]):
                            value = str(row[col])
                            # Handle '>' and '<' in affinity values
                            if value.startswith('>'):
                                value = float(value[1:]) * 10
                                modifier = '>'
                            elif value.startswith('<'):
                                value = float(value[1:]) / 10
                                modifier = '<'
                            else:
                                value = float(value)
                                modifier = ''
                            
                            # Include assay conditions
                            conditions = {
                                'pH': row.get('pH', 'N/A'),
                                'temperature': row.get('Temp (C)', 'N/A'),
                                'source': row.get('Curation/DataSource', 'N/A')
                            }
                            
                            affinities.append({
                                'type': aff_type,
                                'value': value,
                                'modifier': modifier,
                                'conditions': conditions
                            })
                
                if affinities:
                    # Sort by value (lower is better)
                    affinities.sort(key=lambda x: x['value'])
                    best_affinity = affinities[0]
                    
                    # Get first row for compound info
                    row = group.iloc[0]
                    
                    # Convert SMILES to mol for structure analysis
                    mol = None
                    if pd.notna(row['Ligand SMILES']):
                        mol = Chem.MolFromSmiles(str(row['Ligand SMILES']))
                        
                        # Check if molecule has 5-HT2 ligand-like features
                        if mol is not None:
                            is_potential, patterns = self.structure_utils.is_potential_ligand(mol)
                            if not is_potential:
                                continue
                    
                    # Clean compound name for web searches
                    clean_name = self._clean_name(str(name))
                    
                    # Get identifiers
                    identifiers = self.web_client._extract_identifiers(str(name))
                    clean_name = identifiers[0]
                    chembl_id = identifiers[1]
                    cas_number = identifiers[2]
                    patent_id = identifiers[3]
                    
                    # If patent example, extract identifiers using LLM
                    if patent_id and llm_api_key:
                        example_match = re.search(r'(?:Example|Compound)\s*([A-Za-z0-9-]+)', str(name))
                        if example_match:
                            example_id = f"Example {example_match.group(1)}"
                            patent_data = self.web_client._extract_patent_compound(
                                patent_id,
                                example_id,
                                llm_api_key
                            )
                            if patent_data:
                                cas_number = patent_data.get('cas_number', cas_number)
                                smiles = patent_data.get('smiles', str(row['Ligand SMILES']))
                                inchi = patent_data.get('inchi', str(row['Ligand InChI']))
                                clean_name = patent_data.get('name', clean_name)
                    
                    self.logger.info(f"Processing compound: {clean_name}")
                    
                    try:
                        # Get PubChem data
                        self.logger.info("Fetching PubChem data...")
                        pubchem_data = self.web_client._get_pubchem_data(
                            cas=cas_number,
                            smiles=str(row['Ligand SMILES']),
                            inchi=str(row['Ligand InChI'])
                        )
                        
                        self.logger.info("Creating compound data structure...")

                        compound_data = {
                            'name': clean_name,
                            'original_name': str(name),
                            'smiles': str(row['Ligand SMILES']),
                            'inchi': str(row['Ligand InChI']),
                            'inchi_key': str(row['Ligand InChI Key']),
                            'affinity_type': best_affinity['type'],
                            'affinity_value': best_affinity['value'],
                            'affinity_unit': 'nM',
                            'affinity_modifier': best_affinity['modifier'],
                            'assay_conditions': best_affinity['conditions'],
                            'activity_type': self._determine_activity_type(
                                str(row.get('Assay Description', ''))
                            ),
                            'target': str(row['Target Name']),
                            'target_gene': str(row.get('UniProt (SwissProt) Primary ID of Target Chain', 'N/A')),
                            'kinetics': {
                                'kon': str(row.get('kon (M-1-s-1)', 'N/A')),
                                'koff': str(row.get('koff (s-1)', 'N/A'))
                            },
                            'identifiers': {
                                'pubchem_cid': str(row['PubChem CID']),
                                'pubchem_sid': str(row['PubChem SID']),
                                'pubchem_aid': str(row.get('PubChem AID', 'N/A')),
                                'chembl_id': str(row['ChEMBL ID of Ligand']),
                                'chebi_id': str(row.get('ChEBI ID of Ligand', 'N/A')),
                                'drugbank_id': str(row.get('DrugBank ID of Ligand', 'N/A')),
                                'iuphar_id': str(row.get('IUPHAR_GRAC ID of Ligand', 'N/A')),
                                'kegg_id': str(row.get('KEGG ID of Ligand', 'N/A')),
                                'zinc_id': str(row.get('ZINC ID of Ligand', 'N/A')),
                                'pdb_het': str(row.get('Ligand HET ID in PDB', 'N/A')),
                                'pdb_complexes': str(row.get('PDB ID(s) for Ligand-Target Complex', 'N/A')).split(',')
                            },
                            'references': {
                                'doi': str(row.get('Article DOI', 'N/A')),
                                'pmid': str(row.get('PMID', 'N/A')),
                                'authors': str(row.get('Authors', 'N/A')),
                                'institution': str(row.get('Institution', 'N/A')),
                                'bindingdb_ligand': str(row.get('Link to Ligand in BindingDB', 'N/A')),
                                'bindingdb_target': str(row.get('Link to Target in BindingDB', 'N/A')),
                                'bindingdb_pair': str(row.get('Link to Ligand-Target Pair in BindingDB', 'N/A'))
                            },
                            'cas_number': cas_number,
                            'patent_number': str(row.get('Patent Number', 'N/A')),
                            'reference_urls': self.web_client.get_reference_urls(clean_name),
                            'matching_patterns': patterns if mol is not None else [],
                            'species': str(row.get('Target Source Organism According to Curator', 'N/A')),
                            'data_source': str(row.get('Curation/DataSource', 'N/A'))
                        }

                        
                        # Get common names
                        self.logger.info("Fetching common names...")
                        common_names = self.web_client.get_common_names(
                            clean_name,
                            chembl_id=chembl_id,
                            smiles=str(row['Ligand SMILES']),
                            inchi=str(row['Ligand InChI'])
                        )
                        compound_data['common_names'] = [name['name'] for name in common_names]
                        
                        # Get legal status
                        self.logger.info("Fetching legal status...")
                        compound_data['legal_status'] = self.web_client.get_legal_status(clean_name)
                        
                        # Get pharmacology
                        self.logger.info("Fetching pharmacology data...")
                        compound_data['pharmacology'] = self.web_client.get_pharmacology(clean_name)
                        
                        # Add patent information if API key provided
                        if llm_api_key:
                            self.logger.info("Searching patents...")
                            patent_info = self._search_patents(
                                clean_name,
                                str(row['Ligand SMILES']),
                                llm_api_key
                            )
                            compound_data.update(patent_info)
                        
                        # Add if valid structure
                        if self._validate_structure(compound_data):
                            compounds.append(compound_data)
                            self.logger.info(f"Successfully processed compound: {clean_name}")
                        else:
                            self.logger.warning(f"Invalid structure for compound: {clean_name}")
                            
                    except Exception as e:
                        self.logger.error(f"Error processing compound {clean_name}: {str(e)}")
                        continue
            
            # Add compounds from patents if API key provided
            if llm_api_key:
                self.logger.info("Searching patents for additional compounds...")
                patent_compounds = self._extract_compounds_from_patents(
                    '5-HT2 receptor ligand',
                    llm_api_key
                )
                compounds.extend(patent_compounds)
            
            # Deduplicate by InChI Key
            self.logger.info("Deduplicating compounds by InChI Key...")
            unique_compounds = []
            seen_keys = set()
            for comp in compounds:
                if comp['inchi_key'] not in seen_keys:
                    seen_keys.add(comp['inchi_key'])
                    unique_compounds.append(comp)
            
            self.logger.info(f"Found {len(unique_compounds)} unique compounds")
            return unique_compounds
            
        except Exception as e:
            self.logger.error(f"Error gathering 5-HT2 ligands: {str(e)}")
            return []


    def save_5ht2_ligands(self, output_path: str, llm_api_key: Optional[str] = None) -> None:
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
                'name',
                'chembl_id',
                'cas_number',
                'smiles',
                'inchi',
                'inchi_key',
                'common_name_1',
                'common_name_2',
                'common_name_3',
                'target',
                'target_gene',
                'affinity_type',
                'affinity_value',
                'affinity_unit',
                'affinity_modifier',
                'activity_type',
                'kon',
                'koff',
                'pH',
                'temperature',
                'species',
                'pdb_het',
                'pdb_complexes',
                'pubchem_cid',
                'pubchem_sid',
                'pubchem_aid',
                'chebi_id',
                'drugbank_id',
                'iuphar_id',
                'kegg_id',
                'zinc_id',
                'doi',
                'pmid',
                'authors',
                'institution',
                'legal_status',
                'mechanism_of_action',
                'metabolism',
                'toxicity',
                'chembl_url',
                'pubchem_url',
                'wikipedia_url',
                'psychonaut_url',
                'erowid_url',
                'bindingdb_ligand_url',
                'bindingdb_target_url',
                'bindingdb_pair_url'
            ]

            
            self.logger.info("Step 2/5: Converting data to DataFrame...")
            # Convert to DataFrame for easier handling
            df = pd.DataFrame(compounds)
            total_compounds = len(df)
            self.logger.info(f"Found {total_compounds:,} unique compounds")
            
            self.logger.info("Step 3/5: Processing common names...")
            # Add common names columns with progress
            for i in range(3):
                self.logger.info(f"Processing common name column {i + 1}/3...")
                df[f'common_name_{i + 1}'] = df['common_names'].apply(
                    lambda x: x[i] if i < len(x) else ''
                )
            
            self.logger.info("Step 4/5: Processing additional data columns...")
            # Add legal status and pharmacology columns with progress
            processed = 0
            for idx, row in df.iterrows():
                processed += 1
                if processed % 100 == 0:
                    self.logger.info(f"Processed {processed}/{total_compounds} compounds...")
                
                # Legal status
                df.at[idx, 'legal_status'] = '; '.join(
                    f"{s['jurisdiction']}: {s['schedule']}"
                    for s in row['legal_status']['scheduling']
                ) if isinstance(row.get('legal_status'), dict) else ''
                
                # Pharmacology
                if isinstance(row.get('pharmacology'), dict):
                    df.at[idx, 'mechanism_of_action'] = '; '.join(row['pharmacology'].get('mechanism_of_action', []))
                    df.at[idx, 'metabolism'] = '; '.join(row['pharmacology'].get('metabolism', []))
                    df.at[idx, 'toxicity'] = '; '.join(row['pharmacology'].get('toxicity', []))
            
            self.logger.info("Step 5/5: Adding reference URLs...")
            # Add reference URLs with progress
            url_types = ['chembl', 'pubchem', 'wikipedia', 'psychonaut', 'erowid']
            for i, url_type in enumerate(url_types, 1):
                self.logger.info(f"Processing URL type {i} / {len(url_types)}: {url_type}")
                df[f'{url_type}_url'] = df.apply(
                    lambda row: row['reference_urls'].get(f'{url_type}_url', '')
                    if isinstance(row['reference_urls'], dict) else '',
                    axis=1
                )
            
            # Save as TSV
            df[columns].to_csv(output_path, sep='\t', index=False)
            
            # Log statistics
            self.logger.info(
                f"Saved {len(compounds)} 5-HT2 receptor ligands to {output_path}"
            )
            
            sources = set()
            with_cas = 0
            with_patents = 0
            activity_types = {}
            species_counts = {}
            
            for comp in compounds:
                if 'source' in comp:
                    sources.add(comp['source'])
                if comp.get('cas_number', 'N/A') != 'N/A':
                    with_cas += 1
                if comp.get('patents', []):
                    with_patents += 1
                    
                activity = comp.get('activity_type', 'unknown')
                activity_types[activity] = activity_types.get(activity, 0) + 1
                
                species = comp.get('species', 'N/A')
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
        if not compound.get('smiles'):
            return False
            
        try:
            mol = Chem.MolFromSmiles(str(compound['smiles']))
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
        name = re.sub(r'::CHEMBL\d+$', '', name)
        
        # Remove stereochemistry
        name = re.sub(r'\([^)]*\)', '', name)
        name = re.sub(r'[RS]-', '', name)
        name = re.sub(r'\b[RS]\b', '', name)
        
        # Remove special characters but keep some important ones
        name = re.sub(r'[^\w\s\-\+]', '', name)
        
        # Normalize whitespace
        name = ' '.join(name.split())
        
        # Remove common prefixes/suffixes
        name = re.sub(r'^(?:\+\/\-|\+|\-)\s*', '', name)
        name = re.sub(r'\s*(?:hydrochloride|HCl|salt|hydrate|solvate)$', '', name, flags=re.I)
        
        return name

    def _get_cas_number(
        self,
        name: str,
        smiles: str,
        mol: Optional[Chem.Mol]
    ) -> str:
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
                        if 'InformationList' in data:
                            synonyms = data['InformationList']['Information'][0]['Synonym']
                            for syn in synonyms:
                                # Look for CAS pattern
                                if re.match(r'^\d{1,7}-\d{2}-\d$', syn):
                                    return syn
            
            # Try web search
            search_results = self.web_client.get_common_names(name)
            for result in search_results:
                if re.match(r'^\d{1,7}-\d{2}-\d$', result['name']):
                    return result['name']
                    
        except Exception as e:
            self.logger.error(f"Error getting CAS number: {str(e)}")
            
        return 'N/A'

    def _search_patents(
        self,
        name: str,
        smiles: str,
        api_key: str
    ) -> Dict[str, Any]:
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
            headers = {'Authorization': f'Bearer {api_key}'}
            response = requests.get(url, headers=headers)
            
            if response.status_code == 200:
                data = response.json()
                patents = []
                
                for result in data.get('results', []):
                    patent = {
                        'number': result.get('patent_number'),
                        'title': result.get('title'),
                        'abstract': result.get('abstract'),
                        'url': f"https://patents.google.com/patent/{result['patent_number']}"
                    }
                    patents.append(patent)
                    
                return {
                    'patents': patents[:5],  # Top 5 most relevant patents
                    'patent_count': len(data.get('results', []))
                }
                
        except Exception as e:
            self.logger.error(f"Error searching patents: {str(e)}")
            
        return {'patents': [], 'patent_count': 0}

    def _extract_compounds_from_patents(
        self,
        search_term: str,
        api_key: str
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
            headers = {'Authorization': f'Bearer {api_key}'}
            response = requests.get(url, headers=headers)
            
            if response.status_code == 200:
                data = response.json()
                
                for result in data.get('results', []):
                    # Extract text and chemical structures
                    text = f"{result.get('title', '')} {result.get('abstract', '')}"
                    
                    # Use LLM to extract chemical information
                    extracted = self._extract_chemical_info(text, api_key)
                    
                    for compound in extracted:
                        if 'smiles' in compound:
                            # Validate structure
                            mol = Chem.MolFromSmiles(compound['smiles'])
                            if mol is not None:
                                compounds.append({
                                    'name': compound.get('name', 'Unknown'),
                                    'smiles': compound['smiles'],
                                    'inchi': Chem.MolToInchi(mol),
                                    'inchi_key': Chem.MolToInchiKey(mol),
                                    'source': 'Patent',
                                    'patent_number': result.get('patent_number'),
                                    'patent_url': f"https://patents.google.com/patent/{result['patent_number']}"
                                })
                                
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
            # Prepare prompt
            prompt = f"""Extract chemical compounds from the following text. For each compound, provide:
            1. Chemical name
            2. SMILES string (if available)
            3. Any mentioned activity or properties

            Text: {text}

            Format the response as a list of JSON objects."""

            # Make API request
            response = requests.post(
                "https://api.openai.com/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {api_key}",
                    "Content-Type": "application/json"
                },
                json={
                    "model": "gpt-4",
                    "messages": [{"role": "user", "content": prompt}],
                    "temperature": 0.7
                }
            )
            
            if response.status_code == 200:
                result = response.json()
                try:
                    # Parse the response as JSON
                    compounds = json.loads(
                        result['choices'][0]['message']['content']
                    )
                    return compounds
                except json.JSONDecodeError:
                    return []
                    
        except Exception as e:
            self.logger.error(f"Error extracting chemical info: {str(e)}")
            
        return []
