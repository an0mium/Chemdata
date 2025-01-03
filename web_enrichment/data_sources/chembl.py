"""ChEMBL data source functionality.

This module handles:
1. Compound name retrieval
2. Pharmacological data extraction
3. Target information lookup
4. Bioactivity data analysis
5. Structure-activity relationships
"""

from typing import Dict, Any, List, Optional, Tuple
import re

from logger import LogManager
from ..http_client import HttpClient

logger = LogManager().get_logger("web_enrichment.data_sources.chembl")


class ChEMBLClient:
    """Client for interacting with ChEMBL API."""
    
    # Target patterns for 5-HT2 receptors
    TARGET_PATTERNS = {
        '5HT2A': [
            r'5-HT2A',
            r'HTR2A',
            r'serotonin 2A',
            r'5-hydroxytryptamine 2A',
            r'5HT2A',
            r'5-HT-2A',
            r'HT2A_HUMAN',
            r'serotonin receptor 2A',
            r'5-hydroxytryptamine receptor 2A'
        ],
        '5HT2B': [
            r'5-HT2B',
            r'HTR2B',
            r'serotonin 2B',
            r'5-hydroxytryptamine 2B',
            r'5HT2B',
            r'5-HT-2B',
            r'HT2B_HUMAN',
            r'serotonin receptor 2B',
            r'5-hydroxytryptamine receptor 2B'
        ],
        '5HT2C': [
            r'5-HT2C',
            r'HTR2C',
            r'serotonin 2C',
            r'5-hydroxytryptamine 2C',
            r'5HT2C',
            r'5-HT-2C',
            r'HT2C_HUMAN',
            r'serotonin receptor 2C',
            r'5-hydroxytryptamine receptor 2C'
        ]
    }
    
    # Activity type patterns with expanded categories
    ACTIVITY_PATTERNS = {
        'Ki': [
            r'(?:inhibition constant|Ki value)',
            r'binding affinity',
            r'equilibrium dissociation constant'
        ],
        'IC50': [
            r'(?:IC50|half maximal inhibitory concentration)',
            r'inhibitory concentration 50',
            r'concentration for 50% inhibition'
        ],
        'EC50': [
            r'(?:EC50|half maximal effective concentration)',
            r'effective concentration 50',
            r'concentration for 50% effect'
        ],
        'Kd': [
            r'(?:dissociation constant|Kd value)',
            r'binding constant',
            r'equilibrium constant'
        ],
        'Potency': [
            r'(?:potency|activity)',
            r'binding potency',
            r'functional potency'
        ],
        'Efficacy': [
            r'(?:efficacy|maximal response|Emax)',
            r'intrinsic activity',
            r'maximal effect'
        ],
        'Selectivity': [
            r'selectivity',
            r'binding selectivity',
            r'receptor selectivity'
        ]
    }
    
    # Activity mechanism patterns
    MECHANISM_PATTERNS = {
        'superagonist': [
            r'super.?agonist',
            r'high.?efficacy.?agonist',
            r'full.?agonist.+high.?efficacy',
            r'efficacy\s*>\s*100%'
        ],
        'full_agonist': [
            r'full.?agonist',
            r'complete.?agonist',
            r'full.?receptor.?activation',
            r'efficacy\s*[~≈≃]\s*100%'
        ],
        'partial_agonist': [
            r'partial.?agonist',
            r'submaximal.?activation',
            r'partial.?receptor.?activation',
            r'efficacy\s*[<≈]\s*\d{1,2}%'
        ],
        'weak_partial_agonist': [
            r'weak.?partial.?agonist',
            r'low.?efficacy.?partial',
            r'weak.?partial.?activation',
            r'efficacy\s*<\s*20%'
        ],
        'mixed_agonist_antagonist': [
            r'mixed.?agonist.?antagonist',
            r'partial.?agonist.?antagonist',
            r'dual.?activity',
            r'context.?dependent'
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
            r'negative.?efficacy'
        ],
        'positive_allosteric_modulator': [
            r'positive.?allosteric',
            r'PAM',
            r'positive.?modulator',
            r'allosteric.?potentiator'
        ],
        'negative_allosteric_modulator': [
            r'negative.?allosteric',
            r'NAM',
            r'negative.?modulator',
            r'allosteric.?inhibitor'
        ],
        'complex_modulator': [
            r'complex.?modulator',
            r'mixed.?modulator',
            r'complex.?pharmacology',
            r'bitopic'
        ]
    }

    
    def __init__(self, http_client: HttpClient):
        """Initialize ChEMBL client."""
        self.http = http_client

    def get_compound_data(
        self,
        chembl_id: str,
        include_bioactivities: bool = True
    ) -> Optional[Dict[str, Any]]:
        """
        Get comprehensive compound data from ChEMBL.
        
        Args:
            chembl_id: ChEMBL ID
            include_bioactivities: Whether to include bioactivity data
            
        Returns:
            Dictionary containing compound data or None
        """
        try:
            # Get basic compound data
            url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}"
            response = self.http.make_request(url)
            if not response:
                return None
                
            data = response.json()
            
            # Extract compound information
            compound_data = {
                'chembl_id': chembl_id,
                'url': self.get_compound_url(chembl_id),
                'names': self.get_compound_names(chembl_id),
                'properties': self._extract_properties(data),
                'cross_references': self._extract_cross_references(data)
            }
            
            # Get mechanism of action data
            mechanism = self.get_pharmacology(chembl_id)
            if mechanism:
                compound_data['pharmacology'] = mechanism
            
            # Get bioactivity data if requested
            if include_bioactivities:
                bioactivities = self.get_bioactivity_data(chembl_id)
                if bioactivities:
                    compound_data['bioactivities'] = bioactivities
                    
            return compound_data
            
        except Exception as e:
            logger.error(f"Error getting ChEMBL compound data: {str(e)}")
            return None

    def get_compound_names(
        self,
        chembl_id: str
    ) -> List[Dict[str, Any]]:
        """
        Get compound names from ChEMBL.
        
        Args:
            chembl_id: ChEMBL ID
            
        Returns:
            List of dictionaries containing name information
        """
        names = []
        try:
            url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}"
            response = self.http.make_request(url)
            if response:
                data = response.json()
                
                # Add preferred name
                if 'pref_name' in data:
                    names.append({
                        'name': data['pref_name'],
                        'type': 'preferred',
                        'source': 'ChEMBL',
                        'relevance': 100
                    })
                    
                # Add systematic name
                if 'molecule_properties' in data:
                    if 'full_molformula' in data['molecule_properties']:
                        names.append({
                            'name': data['molecule_properties']['full_molformula'],
                            'type': 'systematic',
                            'source': 'ChEMBL',
                            'relevance': 90
                        })
                        
                # Get synonyms with types
                if 'molecule_synonyms' in data:
                    for syn in data['molecule_synonyms']:
                        if 'synonym' in syn:
                            syn_type = syn.get('syn_type', 'other').lower()
                            relevance = self._get_name_relevance(syn_type)
                            
                            names.append({
                                'name': syn['synonym'],
                                'type': syn_type,
                                'source': 'ChEMBL',
                                'relevance': relevance
                            })
                            
                # Get cross-references
                if 'cross_references' in data:
                    for ref in data['cross_references']:
                        if ref.get('xref_name'):
                            names.append({
                                'name': ref['xref_name'],
                                'type': 'cross_reference',
                                'source': ref.get('xref_src', 'ChEMBL'),
                                'relevance': 85
                            })
                            
        except Exception as e:
            logger.error(f"Error getting ChEMBL names: {str(e)}")
            
        return sorted(names, key=lambda x: x['relevance'], reverse=True)

    def _get_name_relevance(self, name_type: str) -> int:
        """Determine name relevance score based on type."""
        relevance_scores = {
            'iupac': 95,
            'inn': 90,
            'research_code': 85,
            'trade_name': 80,
            'common': 75,
            'other': 70
        }
        return relevance_scores.get(name_type, 70)

    def get_pharmacology(self, chembl_id: str) -> Optional[Dict[str, Any]]:
        """
        Get pharmacological information from ChEMBL.
        
        Args:
            chembl_id: ChEMBL ID
            
        Returns:
            Dictionary containing pharmacological data or None
        """
        try:
            # Get mechanism of action data
            url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism/{chembl_id}"
            response = self.http.make_request(url)
            if not response:
                return None
                
            data = response.json()
            info = {
                'mechanisms': [],
                'targets': [],
                'actions': [],
                'binding_sites': []
            }
            
            if 'mechanisms' in data:
                for mech in data['mechanisms']:
                    mechanism = {
                        'mechanism': mech.get('mechanism_of_action'),
                        'target': mech.get('target_name'),
                        'action_type': mech.get('action_type'),
                        'binding_site': mech.get('binding_site_name'),
                        'references': []
                    }
                    
                    # Add mechanism references
                    if 'mechanism_refs' in mech:
                        for ref in mech['mechanism_refs']:
                            if 'ref_type' in ref and 'ref_id' in ref:
                                mechanism['references'].append({
                                    'type': ref['ref_type'],
                                    'id': ref['ref_id']
                                })
                                
                    info['mechanisms'].append(mechanism)
                    
                    # Add to individual lists
                    if mech.get('target_name'):
                        info['targets'].append(mech['target_name'])
                    if mech.get('action_type'):
                        info['actions'].append(mech['action_type'])
                    if mech.get('binding_site_name'):
                        info['binding_sites'].append(mech['binding_site_name'])
                        
            # Remove duplicates
            info['targets'] = list(set(info['targets']))
            info['actions'] = list(set(info['actions']))
            info['binding_sites'] = list(set(info['binding_sites']))
            
            return info
            
        except Exception as e:
            logger.error(f"Error getting ChEMBL pharmacology: {str(e)}")
            return None


    def get_bioactivity_data(
        self,
        chembl_id: str,
        target_type: Optional[str] = None,
        include_raw: bool = False
    ) -> Optional[Dict[str, Any]]:
        """
        Get bioactivity data from ChEMBL.
        
        Args:
            chembl_id: ChEMBL ID
            target_type: Optional target type to filter (e.g., '5HT2A')
            include_raw: Whether to include raw activity data
            
        Returns:
            Dictionary containing bioactivity data or None
        """
        try:
            url = f"https://www.ebi.ac.uk/chembl/api/data/activity/{chembl_id}"
            response = self.http.make_request(url)
            if not response:
                return None
                
            data = response.json()
            activities = []
            
            if 'activities' in data:
                for activity in data['activities']:
                    # Extract target information
                    target_info = self._extract_target_info(activity)
                    if target_type and not self._matches_target_type(
                        target_info.get('name', ''),
                        target_type
                    ):
                        continue
                        
                    # Extract activity values
                    activity_data = {
                        'target': target_info,
                        'type': activity.get('standard_type'),
                        'relation': activity.get('standard_relation'),
                        'value': activity.get('standard_value'),
                        'units': activity.get('standard_units'),
                        'activity_comment': activity.get('activity_comment'),
                        'assay': {
                            'type': activity.get('assay_type'),
                            'description': activity.get('assay_description'),
                            'organism': activity.get('assay_organism'),
                            'cell_line': activity.get('assay_cell_type'),
                            'subcellular_fraction': activity.get('assay_subcellular_fraction'),
                            'parameters': activity.get('assay_parameters')
                        },
                        'source': {
                            'type': activity.get('src_id'),
                            'description': activity.get('src_description')
                        }
                    }
                    
                    # Add reference if available
                    if 'document_chembl_id' in activity:
                        activity_data['reference'] = {
                            'chembl_id': activity['document_chembl_id'],
                            'year': activity.get('document_year'),
                            'journal': activity.get('journal'),
                            'volume': activity.get('volume'),
                            'issue': activity.get('issue'),
                            'first_page': activity.get('first_page')
                        }
                    
                    # Determine activity mechanism
                    if activity.get('assay_description'):
                        mechanism = self._determine_mechanism(
                            activity['assay_description']
                        )
                        if mechanism:
                            activity_data['mechanism'] = mechanism
                            
                    # Add raw data if requested
                    if include_raw:
                        activity_data['raw_data'] = activity
                        
                    activities.append(activity_data)
                    
            # Group activities by target and type
            grouped = self._group_activities(activities)
            
            # Calculate summary statistics
            summary = self._calculate_activity_summary(activities)
            
            return {
                'activity_count': len(activities),
                'activities': activities,
                'grouped': grouped,
                'summary': summary
            }
            
        except Exception as e:
            logger.error(f"Error getting ChEMBL bioactivity data: {str(e)}")
            return None

    def _determine_mechanism(self, description: str) -> Optional[str]:
        """Determine activity mechanism from description."""
        description = description.lower()
        
        for mechanism, patterns in self.MECHANISM_PATTERNS.items():
            for pattern in patterns:
                if re.search(pattern, description, re.I):
                    return mechanism
                    
        return None

    def _calculate_activity_summary(
        self,
        activities: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """Calculate summary statistics for activities."""
        summary = {
            'activity_types': {},
            'mechanisms': {},
            'target_organisms': {},
            'assay_types': {},
            'value_ranges': {}
        }
        
        for activity in activities:
            # Count activity types
            act_type = activity['type']
            summary['activity_types'][act_type] = summary['activity_types'].get(act_type, 0) + 1
            
            # Count mechanisms
            if 'mechanism' in activity:
                mech = activity['mechanism']
                summary['mechanisms'][mech] = summary['mechanisms'].get(mech, 0) + 1
            
            # Count target organisms
            org = activity['target']['organism']
            if org:
                summary['target_organisms'][org] = summary['target_organisms'].get(org, 0) + 1
            
            # Count assay types
            assay_type = activity['assay']['type']
            if assay_type:
                summary['assay_types'][assay_type] = summary['assay_types'].get(assay_type, 0) + 1
            
            # Track value ranges
            if activity['value'] is not None:
                if act_type not in summary['value_ranges']:
                    summary['value_ranges'][act_type] = {
                        'min': float('inf'),
                        'max': float('-inf'),
                        'count': 0,
                        'sum': 0
                    }
                
                value = float(activity['value'])
                stats = summary['value_ranges'][act_type]
                stats['min'] = min(stats['min'], value)
                stats['max'] = max(stats['max'], value)
                stats['count'] += 1
                stats['sum'] += value
                
        # Calculate averages
        for act_type, stats in summary['value_ranges'].items():
            if stats['count'] > 0:
                stats['average'] = stats['sum'] / stats['count']
                del stats['sum']
                
        return summary


    def _extract_target_info(self, activity: Dict[str, Any]) -> Dict[str, Any]:
        """Extract target information from activity data."""
        return {
            'name': activity.get('target_pref_name'),
            'type': activity.get('target_type'),
            'organism': activity.get('target_organism'),
            'chembl_id': activity.get('target_chembl_id')
        }

    def _matches_target_type(self, target_name: str, target_type: str) -> bool:
        """Check if target name matches target type patterns."""
        if target_type in self.TARGET_PATTERNS:
            patterns = self.TARGET_PATTERNS[target_type]
            return any(
                re.search(pattern, target_name, re.I)
                for pattern in patterns
            )
        return False

    def _group_activities(
        self,
        activities: List[Dict[str, Any]]
    ) -> Dict[str, Dict[str, List[Dict[str, Any]]]]:
        """Group activities by target and activity type."""
        grouped = {}
        
        for activity in activities:
            target = activity['target']['name']
            if target not in grouped:
                grouped[target] = {}
                
            act_type = activity['type']
            if act_type not in grouped[target]:
                grouped[target][act_type] = []
                
            grouped[target][act_type].append(activity)
            
        return grouped

    def _extract_properties(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Extract compound properties from ChEMBL data."""
        properties = {}
        
        if 'molecule_properties' in data:
            props = data['molecule_properties']
            
            # Map property names
            property_map = {
                'full_molformula': 'molecular_formula',
                'full_mwt': 'molecular_weight',
                'alogp': 'alogp',
                'psa': 'polar_surface_area',
                'rtb': 'rotatable_bonds',
                'ro3_pass': 'rule_of_three_compliant',
                'num_ro5_violations': 'rule_of_five_violations',
                'cx_logp': 'logp',
                'cx_logd': 'logd',
                'aromatic_rings': 'aromatic_rings',
                'hba': 'hbond_acceptors',
                'hbd': 'hbond_donors'
            }
            
            for chembl_name, std_name in property_map.items():
                if chembl_name in props:
                    properties[std_name] = props[chembl_name]
                    
        return properties

    def _extract_cross_references(self, data: Dict[str, Any]) -> Dict[str, str]:
        """Extract cross-references from ChEMBL data."""
        references = {}
        
        if 'cross_references' in data:
            for ref in data['cross_references']:
                if 'xref_id' in ref and 'xref_src' in ref:
                    references[ref['xref_src'].lower()] = ref['xref_id']
                    
        return references

    def get_compound_url(self, chembl_id: str) -> str:
        """
        Get ChEMBL compound report card URL.
        
        Args:
            chembl_id: ChEMBL ID
            
        Returns:
            URL to ChEMBL compound report card
        """
        return f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}"
