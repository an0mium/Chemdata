"""PubChem data source functionality.

This module handles:
1. Compound data retrieval
2. Structure search
3. Name lookup
4. Pharmacological data extraction
5. Bioassay data retrieval
6. Activity data analysis
"""

import time
from typing import Dict, Any, List, Optional, Tuple
from urllib.parse import quote
import re
from tqdm import tqdm

from bs4 import BeautifulSoup
from logger import LogManager
from ..http_client import HttpClient

logger = LogManager().get_logger("web_enrichment.data_sources.pubchem")


class PubChemClient:
    """Client for interacting with PubChem API."""
    
    # Bioassay types of interest
    BIOASSAY_TYPES = {
        'binding': [
            'binding assay',
            'radioligand binding',
            'competition binding',
            'saturation binding',
            'displacement assay',
            'affinity determination',
            'receptor binding',
            'ligand binding'
        ],
        'functional': [
            'functional assay',
            'cell-based assay',
            'reporter gene',
            'calcium flux',
            'gtpγs binding',
            'camp assay',
            'β-arrestin',
            'electrophysiology',
            'patch clamp',
            'membrane potential',
            'fluorescence imaging'
        ],
        'enzymatic': [
            'enzyme assay',
            'kinase assay',
            'phosphorylation',
            'catalytic activity',
            'substrate turnover',
            'enzyme inhibition'
        ],
        'pharmacological': [
            'agonist activity',
            'antagonist activity',
            'inverse agonist',
            'allosteric modulation',
            'potentiation assay',
            'inhibition assay'
        ]
    }
    
    # Activity type patterns
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

    
    def __init__(self, http_client: HttpClient):
        """Initialize PubChem client."""
        self.http = http_client

    def get_compound_data(
        self,
        cas: Optional[str] = None,
        smiles: Optional[str] = None,
        inchi: Optional[str] = None,
        include_bioassays: bool = True
    ) -> Optional[Dict[str, Any]]:
        """
        Get compound data from PubChem using identifiers.
        
        Args:
            cas: CAS number
            smiles: SMILES string
            inchi: InChI string
            include_bioassays: Whether to include bioassay data
            
        Returns:
            Dictionary containing PubChem data or None
        """
        try:
            # Set up progress bar
            progress = tqdm(total=4, desc="Fetching PubChem data", unit="steps")
            
            cid = None
            progress.set_description("Searching by structure")
            
            # Try structure search first (more reliable)
            if smiles:
                cid = self._search_by_structure('smiles', smiles)
            
            # Try InChI
            if not cid and inchi:
                cid = self._search_by_structure('inchi', inchi)
            
            # Try CAS as a fallback
            if not cid and cas:
                cid = self._search_by_identifier(cas)
            
            progress.update(1)
                
            if cid:
                # Get basic compound data
                progress.set_description("Getting compound data")
                compound_data = self.get_compound(cid)
                progress.update(1)
                
                if compound_data:
                    # Add bioassay data if requested
                    if include_bioassays:
                        progress.set_description("Getting bioassay data")
                        bioassays = self.get_bioassay_data(cid)
                        if bioassays:
                            compound_data['bioassays'] = bioassays
                        progress.update(1)
                        
                    # Get pharmacology data
                    progress.set_description("Getting pharmacology data")
                    pharmacology = self.get_pharmacology(str(cid))
                    if pharmacology:
                        compound_data['pharmacology'] = pharmacology
                    progress.update(1)
                    
                    progress.close()
                    return compound_data
                    
        except Exception as e:
            logger.error(f"Error getting PubChem data: {str(e)}")
            if 'progress' in locals():
                progress.close()
            
        return None

    def get_bioassay_data(self, cid: str) -> Optional[Dict[str, Any]]:
        """
        Get bioassay data for compound.
        
        Args:
            cid: PubChem Compound ID
            
        Returns:
            Dictionary containing bioassay data or None
        """
        try:
            # Get active assays for compound
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
            response = self.http.make_request(url)
            if not response:
                return None
                
            data = response.json()
            if 'AssaySummaries' not in data:
                return None
                
            assays = []
            activities = {
                'binding': [],
                'functional': [],
                'enzymatic': [],
                'pharmacological': []
            }
            
            # Process assays with progress bar
            assay_progress = tqdm(
                data['AssaySummaries'],
                desc="Processing bioassays",
                unit="assays"
            )
            
            for summary in assay_progress:
                aid = summary.get('AID')
                if not aid:
                    continue
                    
                assay_progress.set_description(f"Processing assay {aid}")
                
                # Get detailed assay data
                assay_data = self._get_assay_details(aid)
                if assay_data:
                    assay_type = self._classify_assay(assay_data)
                    if assay_type:  # Only include relevant assays
                        activity_type = self._determine_activity_type(
                            assay_data.get('description', '')
                        )
                        
                        assay_info = {
                            'aid': aid,
                            'name': assay_data.get('name', ''),
                            'description': assay_data.get('description', ''),
                            'type': assay_type,
                            'activity_type': activity_type,
                            'target': assay_data.get('target', ''),
                            'target_gene': assay_data.get('target_gene', ''),
                            'target_protein': assay_data.get('target_protein', ''),
                            'activity': summary.get('ActivityOutcome', ''),
                            'value': summary.get('ActivityValue'),
                            'unit': summary.get('ActivityUnit'),
                            'conditions': assay_data.get('conditions', {}),
                            'url': f"https://pubchem.ncbi.nlm.nih.gov/bioassay/{aid}"
                        }
                        
                        assays.append(assay_info)
                        activities[assay_type].append(assay_info)
            
            assay_progress.close()
            return {
                'assay_count': len(assays),
                'assays': assays,
                'activities': activities
            }
            
        except Exception as e:
            logger.error(f"Error getting bioassay data: {str(e)}")
            if 'assay_progress' in locals():
                assay_progress.close()
            
        return None


    def _search_by_structure(
        self,
        structure_type: str,
        structure: str
    ) -> Optional[str]:
        """Search PubChem by chemical structure."""
        try:
            # Clean and encode structure
            structure = structure.strip()
            if structure_type == 'smiles':
                # Remove any whitespace and escape special characters
                structure = re.sub(r'\s+', '', structure)
                # URL encode the structure
                structure = quote(structure)
            elif structure_type == 'inchi':
                # Ensure proper InChI format
                if not structure.startswith('InChI='):
                    structure = f"InChI={structure}"
                # URL encode the structure
                structure = quote(structure)
                    
            # Use POST request for long structures
            if len(structure) > 1000:
                url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/structure/cids/JSON"
                response = self.http.make_request(
                    url,
                    method='POST',
                    data={structure_type: structure}
                )
            else:
                url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/structure/cids/JSON"
                response = self.http.make_request(
                    url,
                    params={structure_type: structure}
                )
                
            if response and 'IdentifierList' in response.json():
                return response.json()['IdentifierList']['CID'][0]
                
        except Exception as e:
            logger.error(f"Error searching by structure: {str(e)}")
            
        return None

    def _search_by_identifier(self, identifier: str) -> Optional[str]:
        """Search PubChem by identifier with async support."""
        try:
            url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/search/JSON"
            response = self.http.make_request(url, params={'name': identifier})
            
            if response:
                data = response.json()
                if 'Waiting' in data:
                    # Handle asynchronous search
                    listkey = data['Waiting']['ListKey']
                    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{listkey}/cids/JSON"
                    
                    # Poll for results
                    for _ in range(5):  # Try up to 5 times
                        time.sleep(1)
                        response = self.http.make_request(url)
                        if response:
                            data = response.json()
                            if 'IdentifierList' in data:
                                return data['IdentifierList']['CID'][0]
                            elif 'Waiting' not in data:
                                break
                                
        except Exception as e:
            logger.error(f"Error searching by identifier: {str(e)}")
        return None

    def get_compound(self, cid: str) -> Optional[Dict[str, Any]]:
        """
        Get compound data from PubChem CID.
        
        Args:
            cid: PubChem Compound ID
            
        Returns:
            Dictionary containing compound data or None
        """
        try:
            # Get basic compound data
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
            response = self.http.make_request(url)
            if not response:
                return None
                
            data = response.json()
            if 'PC_Compounds' not in data:
                return None
                
            compound = data['PC_Compounds'][0]
            
            # Extract properties
            properties = self._extract_properties(compound)
            
            # Get computed properties
            computed = self._get_computed_properties(cid)
            if computed:
                properties.update(computed)
            
            # Get synonyms
            synonyms = self.get_compound_names(cid)
            
            # Get references
            references = self._get_references(cid)
            
            return {
                'cid': cid,
                'url': f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
                'properties': properties,
                'synonyms': synonyms,
                'references': references,
                'raw_data': compound  # Include raw data for further processing
            }
            
        except Exception as e:
            logger.error(f"Error getting PubChem compound: {str(e)}")
        return None

    def _extract_properties(self, compound: Dict) -> Dict[str, Any]:
        """Extract chemical properties from compound data."""
        properties = {}
        
        try:
            if 'props' in compound:
                for prop in compound['props']:
                    if 'urn' in prop:
                        # Extract property name from URN
                        name = prop['urn'].get('label', '').lower()
                        if 'value' in prop:
                            if 'sval' in prop['value']:
                                properties[name] = prop['value']['sval']
                            elif 'fval' in prop['value']:
                                properties[name] = prop['value']['fval']
                            elif 'ival' in prop['value']:
                                properties[name] = prop['value']['ival']
                                
        except Exception as e:
            logger.error(f"Error extracting properties: {str(e)}")
            
        return properties

    def _get_computed_properties(self, cid: str) -> Optional[Dict[str, Any]]:
        """Get computed properties from PubChem."""
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight,XLogP,TPSA,RotatableBondCount,HBondDonorCount,HBondAcceptorCount,Complexity/JSON"
            response = self.http.make_request(url)
            if response and 'PropertyTable' in response.json():
                props = response.json()['PropertyTable']['Properties'][0]
                return {
                    'molecular_weight': props.get('MolecularWeight'),
                    'xlogp': props.get('XLogP'),
                    'tpsa': props.get('TPSA'),
                    'rotatable_bonds': props.get('RotatableBondCount'),
                    'hbond_donors': props.get('HBondDonorCount'),
                    'hbond_acceptors': props.get('HBondAcceptorCount'),
                    'complexity': props.get('Complexity')
                }
        except Exception as e:
            logger.error(f"Error getting computed properties: {str(e)}")
        return None

    def get_compound_names(self, cid: str) -> List[Dict[str, Any]]:
        """
        Get compound names from PubChem CID.
        
        Args:
            cid: PubChem Compound ID
            
        Returns:
            List of dictionaries containing name information
        """
        names = []
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
            response = self.http.make_request(url)
            if response and 'InformationList' in response.json():
                synonyms = response.json()['InformationList']['Information'][0]['Synonym']
                
                for name in synonyms:
                    # Skip very long names
                    if len(name) > 250:
                        continue
                        
                    # Determine name type and relevance
                    name_type, relevance = self._classify_name(name)
                    
                    names.append({
                        'name': name,
                        'type': name_type,
                        'source': 'PubChem',
                        'relevance': relevance
                    })
                    
        except Exception as e:
            logger.error(f"Error getting PubChem names for CID {cid}: {str(e)}")
            
        return sorted(names, key=lambda x: x['relevance'], reverse=True)

    def _classify_name(self, name: str) -> Tuple[str, int]:
        """Classify chemical name and assign relevance score."""
        name_lower = name.lower()
        
        # Check for CAS number
        if re.match(r'^\d{1,7}-\d{2}-\d$', name):
            return 'cas', 100
            
        # Check for systematic name
        if any(x in name_lower for x in ['iupac', 'systematic']):
            return 'systematic', 90
            
        # Check for registry numbers
        if re.match(r'^[A-Z]{1,3}-\d+$', name):
            return 'registry', 85
            
        # Check for common name indicators
        if len(name) < 30:
            return 'common', 80
            
        # Check for semi-systematic name
        if re.search(r'[0-9]|acid|amine|phenyl', name_lower):
            return 'semi-systematic', 70
            
        return 'other', 60

    def get_pharmacology(
        self,
        name: str,
        cas_number: Optional[str] = None
    ) -> Optional[Dict[str, List[str]]]:
        """
        Get pharmacological information from PubChem.
        
        Args:
            name: Compound name
            cas_number: Optional CAS number
            
        Returns:
            Dictionary containing pharmacological data or None
        """
        try:
            search_terms = [name]
            if cas_number:
                search_terms.append(cas_number)
                
            for term in search_terms:
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(term)}/JSON"
                response = self.http.make_request(url)
                if response:
                    data = response.json()
                    if 'PC_Compounds' in data:
                        compound = data['PC_Compounds'][0]
                        info = {
                            'mechanism': [],
                            'metabolism': [],
                            'toxicity': [],
                            'pharmacokinetics': [],
                            'drug_interactions': []
                        }
                        
                        # Extract pharmacology sections
                        if 'Section' in compound:
                            for section in compound['Section']:
                                if section.get('TOCHeading') == 'Pharmacology':
                                    self._extract_pharmacology_data(section, info)
                                elif section.get('TOCHeading') == 'Toxicity':
                                    info['toxicity'].extend(
                                        self._extract_text_from_section(section)
                                    )
                                    
                        return info
                        
        except Exception as e:
            logger.error(f"Error getting PubChem pharmacology: {str(e)}")
            
        return None

    def _extract_pharmacology_data(
        self,
        section: Dict,
        info: Dict[str, List[str]]
    ) -> None:
        """Extract detailed pharmacology data from section."""
        for subsection in section.get('Section', []):
            heading = subsection.get('TOCHeading', '').lower()
            
            if 'mechanism' in heading or 'action' in heading:
                info['mechanism'].extend(
                    self._extract_text_from_section(subsection)
                )
            elif 'metabolism' in heading:
                info['metabolism'].extend(
                    self._extract_text_from_section(subsection)
                )
            elif 'pharmacokinetics' in heading or 'absorption' in heading:
                info['pharmacokinetics'].extend(
                    self._extract_text_from_section(subsection)
                )
            elif 'interaction' in heading:
                info['drug_interactions'].extend(
                    self._extract_text_from_section(subsection)
                )


    def get_bioassay_data(self, cid: str) -> Optional[Dict[str, Any]]:
        """
        Get bioassay data for compound.
        
        Args:
            cid: PubChem Compound ID
            
        Returns:
            Dictionary containing bioassay data or None
        """
        try:
            # Get active assays for compound
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
            response = self.http.make_request(url)
            if not response:
                return None
                
            data = response.json()
            if 'AssaySummaries' not in data:
                return None
                
            assays = []
            activities = {
                'binding': [],
                'functional': [],
                'enzymatic': [],
                'pharmacological': []
            }
            
            for summary in data['AssaySummaries']:
                aid = summary.get('AID')
                if not aid:
                    continue
                    
                # Get detailed assay data
                assay_data = self._get_assay_details(aid)
                if assay_data:
                    assay_type = self._classify_assay(assay_data)
                    if assay_type:  # Only include relevant assays
                        activity_type = self._determine_activity_type(
                            assay_data.get('description', '')
                        )
                        
                        assay_info = {
                            'aid': aid,
                            'name': assay_data.get('name', ''),
                            'description': assay_data.get('description', ''),
                            'type': assay_type,
                            'activity_type': activity_type,
                            'target': assay_data.get('target', ''),
                            'target_gene': assay_data.get('target_gene', ''),
                            'target_protein': assay_data.get('target_protein', ''),
                            'activity': summary.get('ActivityOutcome', ''),
                            'value': summary.get('ActivityValue'),
                            'unit': summary.get('ActivityUnit'),
                            'conditions': assay_data.get('conditions', {}),
                            'url': f"https://pubchem.ncbi.nlm.nih.gov/bioassay/{aid}"
                        }
                        
                        assays.append(assay_info)
                        activities[assay_type].append(assay_info)
                        
            return {
                'assay_count': len(assays),
                'assays': assays,
                'activities': activities
            }
            
        except Exception as e:
            logger.error(f"Error getting bioassay data: {str(e)}")
            
        return None

    def _get_assay_details(self, aid: str) -> Optional[Dict[str, Any]]:
        """Get detailed assay information."""
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/description/JSON"
            response = self.http.make_request(url)
            if response and 'PC_AssayContainer' in response.json():
                assay = response.json()['PC_AssayContainer'][0]
                
                # Extract assay details
                description = self._extract_text_from_section(
                    assay.get('Description', {})
                )
                
                # Extract target information
                target_info = self._extract_target_info(assay)
                
                # Extract assay conditions
                conditions = self._extract_assay_conditions(assay)
                
                return {
                    'name': assay.get('Name', ''),
                    'description': ' '.join(description),
                    'target': target_info.get('name', ''),
                    'target_gene': target_info.get('gene', ''),
                    'target_protein': target_info.get('protein', ''),
                    'conditions': conditions,
                    'raw_data': assay  # Include raw data for further processing
                }
                
        except Exception as e:
            logger.error(f"Error getting assay details: {str(e)}")
            
        return None

    def _extract_target_info(self, assay: Dict) -> Dict[str, str]:
        """Extract detailed target information."""
        info = {
            'name': '',
            'gene': '',
            'protein': ''
        }
        
        try:
            if 'Target' in assay:
                target = assay['Target']
                info['name'] = target.get('name', '')
                
                # Extract gene information
                if 'mol_id' in target:
                    mol_id = target['mol_id']
                    gene_url = f"https://www.ncbi.nlm.nih.gov/gene/{mol_id}"
                    response = self.http.make_request(gene_url)
                    if response:
                        # Parse gene page for symbol and name
                        soup = BeautifulSoup(response.text, 'html.parser')
                        gene_symbol = soup.find('dd', {'id': 'gene-symbol'})
                        if gene_symbol:
                            info['gene'] = gene_symbol.text.strip()
                
                # Extract protein information
                if 'protein_id' in target:
                    protein_id = target['protein_id']
                    protein_url = f"https://www.ncbi.nlm.nih.gov/protein/{protein_id}"
                    response = self.http.make_request(protein_url)
                    if response:
                        # Parse protein page for name
                        soup = BeautifulSoup(response.text, 'html.parser')
                        protein_name = soup.find('dt', text='Protein name')
                        if protein_name and protein_name.find_next('dd'):
                            info['protein'] = protein_name.find_next('dd').text.strip()
                            
        except Exception as e:
            logger.error(f"Error extracting target info: {str(e)}")
            
        return info

    def _extract_assay_conditions(self, assay: Dict) -> Dict[str, Any]:
        """Extract assay conditions and parameters."""
        conditions = {
            'temperature': None,
            'ph': None,
            'buffer': None,
            'incubation_time': None,
            'substrate': None,
            'cofactors': [],
            'detection_method': None
        }
        
        try:
            if 'Description' in assay:
                description = ' '.join(
                    self._extract_text_from_section(assay['Description'])
                )
                
                # Extract temperature
                temp_match = re.search(
                    r'(?:at|temperature)\s*(?:of)?\s*(\d+)[°\s]*[Cc]',
                    description
                )
                if temp_match:
                    conditions['temperature'] = float(temp_match.group(1))
                
                # Extract pH
                ph_match = re.search(
                    r'pH\s*(?:of)?\s*(\d+\.?\d*)',
                    description
                )
                if ph_match:
                    conditions['ph'] = float(ph_match.group(1))
                
                # Extract buffer
                buffer_match = re.search(
                    r'(?:in|using)\s+([^.]+?buffer)',
                    description,
                    re.I
                )
                if buffer_match:
                    conditions['buffer'] = buffer_match.group(1).strip()
                
                # Extract incubation time
                time_match = re.search(
                    r'(?:for|incubated?\s+for)\s+(\d+)\s*(min|hour|h|hrs?)',
                    description,
                    re.I
                )
                if time_match:
                    value = int(time_match.group(1))
                    unit = time_match.group(2).lower()
                    if unit.startswith('h'):
                        value *= 60  # Convert to minutes
                    conditions['incubation_time'] = value
                
                # Extract detection method
                for method in ['fluorescence', 'luminescence', 'absorbance', 'radioactivity']:
                    if method in description.lower():
                        conditions['detection_method'] = method
                        break
                        
        except Exception as e:
            logger.error(f"Error extracting assay conditions: {str(e)}")
            
        return conditions

    def _determine_activity_type(self, description: str) -> Optional[str]:
        """Determine activity type from assay description."""
        description = description.lower()
        
        for activity_type, patterns in self.ACTIVITY_PATTERNS.items():
            for pattern in patterns:
                if re.search(pattern, description, re.I):
                    return activity_type
                    
        return None


    def _classify_assay(self, assay_data: Dict[str, Any]) -> Optional[str]:
        """Classify assay type based on description."""
        description = assay_data.get('description', '').lower()
        
        # Check each assay type
        for assay_type, patterns in self.BIOASSAY_TYPES.items():
            if any(pattern in description for pattern in patterns):
                return assay_type
                
        return None

    def _get_references(self, cid: str) -> Dict[str, List[Dict[str, Any]]]:
        """Get literature references for compound."""
        references = {
            'patents': [],
            'articles': []
        }
        
        try:
            # Get patent references
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/PatentID/JSON"
            response = self.http.make_request(url)
            if response and 'InformationList' in response.json():
                for ref in response.json()['InformationList']['Information']:
                    if 'PatentID' in ref:
                        references['patents'].append({
                            'id': ref['PatentID'],
                            'source': 'PubChem'
                        })
                        
            # Get article references
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/PubMedID/JSON"
            response = self.http.make_request(url)
            if response and 'InformationList' in response.json():
                for ref in response.json()['InformationList']['Information']:
                    if 'PubMedID' in ref:
                        references['articles'].append({
                            'pmid': ref['PubMedID'],
                            'source': 'PubChem'
                        })
                        
        except Exception as e:
            logger.error(f"Error getting references: {str(e)}")
            
        return references

    def _extract_text_from_section(self, section: Dict) -> List[str]:
        """Extract text content from PubChem section."""
        texts = []
        if 'Information' in section:
            for info in section['Information']:
                if 'Value' in info and 'StringWithMarkup' in info['Value']:
                    for markup in info['Value']['StringWithMarkup']:
                        if 'String' in markup:
                            texts.append(markup['String'])
        return texts
