"""Community data source functionality.

This module handles:
1. PsychonautWiki data extraction
2. Erowid data extraction
3. Pharmacological and toxicity information
4. Experience reports analysis
5. Dosage information
6. Route of administration data
"""

from typing import Dict, List, Optional, Any, Set, Tuple
from urllib.parse import quote, urljoin
from bs4 import BeautifulSoup
import re
import json

from logger import LogManager
from ..http_client import HttpClient
from ..llm_utils import analyze_content_with_llm

logger = LogManager().get_logger("web_enrichment.data_sources.community")


class CommunityClient:
    """Client for interacting with community data sources."""
    
    # Sections of interest with expanded coverage
    SECTIONS = {
        'chemistry': [
            'chemistry',
            'chemical structure',
            'molecular structure',
            'synthesis',
            'chemical composition',
            'physical properties',
            'solubility',
            'stability'
        ],
        'pharmacology': [
            'pharmacology',
            'mechanism of action',
            'pharmacodynamics',
            'pharmacokinetics',
            'metabolism',
            'receptor binding',
            'neurotransmitters',
            'enzyme interactions',
            'serotonin',
            '5-HT',
            'binding affinity'
        ],
        'effects': [
            'effects',
            'subjective effects',
            'physical effects',
            'cognitive effects',
            'psychological effects',
            'sensory effects',
            'auditory effects',
            'visual effects',
            'tactile effects',
            'emotional effects',
            'hallucinogenic',
            'psychedelic'
        ],
        'toxicity': [
            'toxicity',
            'safety',
            'health concerns',
            'adverse effects',
            'side effects',
            'overdose',
            'contraindications',
            'interactions',
            'harm reduction'
        ],
        'dosage': [
            'dosage',
            'dose',
            'dosing',
            'threshold',
            'common doses',
            'strong doses',
            'heavy doses',
            'recommended doses',
            'dosage forms'
        ],
        'roa': [
            'route of administration',
            'administration',
            'roa',
            'routes',
            'oral',
            'sublingual',
            'insufflation',
            'inhalation',
            'injection',
            'rectal'
        ],
        'duration': [
            'duration',
            'onset',
            'comeup',
            'peak',
            'offset',
            'after effects',
            'hangover',
            'timeline'
        ],
        'tolerance': [
            'tolerance',
            'dependence',
            'addiction',
            'withdrawal',
            'cross-tolerance',
            'tolerance reduction'
        ]
    }
    
    # Effect categories with 5-HT2 specific effects
    EFFECT_CATEGORIES = {
        'physical': [
            'stimulation',
            'sedation',
            'analgesia',
            'nausea',
            'muscle tension',
            'respiratory',
            'cardiovascular',
            'pupil dilation',
            'body temperature',
            'blood pressure'
        ],
        'cognitive': [
            'focus',
            'memory',
            'thought acceleration',
            'thought deceleration',
            'creativity',
            'analysis enhancement',
            'confusion',
            'conceptual thinking',
            'introspection'
        ],
        'sensory': [
            'visual',
            'auditory',
            'tactile',
            'olfactory',
            'gustatory',
            'proprioception',
            'synesthesia',
            'geometric patterns',
            'color enhancement'
        ],
        'emotional': [
            'euphoria',
            'anxiety',
            'empathy',
            'mood lift',
            'mood suppression',
            'emotional enhancement',
            'emotional suppression',
            'spiritual experiences',
            'ego dissolution'
        ]
    }
    
    # Receptor activity patterns with 5-HT2 specifics
    RECEPTOR_PATTERNS = {
        'agonist': [
            r'(?:full|partial|super)?\s*agonist',
            r'activates?',
            r'stimulates?',
            r'binds?\s+to\s+and\s+activates?',
            r'5-HT2\s*(?:A|B|C)?\s*agonist'
        ],
        'antagonist': [
            r'antagonist',
            r'blocks?',
            r'inhibits?',
            r'prevents?\s+activation',
            r'5-HT2\s*(?:A|B|C)?\s*antagonist'
        ],
        'modulator': [
            r'(?:positive|negative)?\s*allosteric\s*modulator',
            r'(?:PAM|NAM)',
            r'modulates?',
            r'enhances?|reduces?\s+activity',
            r'5-HT2\s*(?:A|B|C)?\s*modulator'
        ],
        'reuptake': [
            r'reuptake\s*(?:inhibitor|enhancer)',
            r'(?:blocks?|enhances?)\s*reuptake',
            r'(?:SRI|NRI|DRI)',
            r'monoamine(?:\s*reuptake)?\s*(?:inhibitor|enhancer)',
            r'serotonin\s*reuptake'
        ],
        'releaser': [
            r'releasing\s*agent',
            r'releases?',
            r'increases?\s*release',
            r'promotes?\s*release',
            r'serotonin\s*release'
        ]
    }
    
    # 5-HT2 receptor patterns
    RECEPTOR_SUBTYPES = {
        '5HT2A': [
            r'5-HT2A',
            r'5-HT-2A',
            r'5HT2A',
            r'HTR2A',
            r'serotonin\s*2A',
            r'serotonin\s*receptor\s*2A'
        ],
        '5HT2B': [
            r'5-HT2B',
            r'5-HT-2B',
            r'5HT2B',
            r'HTR2B',
            r'serotonin\s*2B',
            r'serotonin\s*receptor\s*2B'
        ],
        '5HT2C': [
            r'5-HT2C',
            r'5-HT-2C',
            r'5HT2C',
            r'HTR2C',
            r'serotonin\s*2C',
            r'serotonin\s*receptor\s*2C'
        ]
    }

    
    def __init__(self, http_client: HttpClient):
        """Initialize community client."""
        self.http = http_client

    def get_compound_data(
        self,
        name: str,
        cas_number: Optional[str] = None,
        llm_api_key: Optional[str] = None
    ) -> Optional[Dict[str, Any]]:
        """
        Get comprehensive compound data from community sources.
        
        Args:
            name: Compound name
            cas_number: Optional CAS number
            llm_api_key: Optional API key for LLM analysis
            
        Returns:
            Dictionary containing compound data or None
        """
        data = {
            'chemistry': [],
            'pharmacology': [],
            'effects': [],
            'toxicity': [],
            'dosage': {},
            'roa': {},
            'experience_reports': [],
            'sources': [],
            'urls': {}
        }
        
        # Get PsychonautWiki data
        psychonaut_data = self._get_psychonaut_data(name, cas_number, llm_api_key)
        if psychonaut_data:
            self._merge_data(data, psychonaut_data)
            data['sources'].append('PsychonautWiki')
            
        # Get Erowid data
        erowid_data = self._get_erowid_data(name, cas_number, llm_api_key)
        if erowid_data:
            self._merge_data(data, erowid_data)
            data['sources'].append('Erowid')
            
        # Get URLs
        urls = self.get_urls(name, cas_number)
        if urls:
            data['urls'].update(urls)
            
        return data if data['sources'] else None

    def _get_psychonaut_data(
        self,
        name: str,
        cas_number: Optional[str],
        llm_api_key: Optional[str]
    ) -> Optional[Dict[str, Any]]:
        """Get comprehensive data from PsychonautWiki."""
        data = {
            'chemistry': [],
            'pharmacology': [],
            'effects': [],
            'toxicity': [],
            'dosage': {},
            'roa': {},
            'experience_reports': []
        }
        
        try:
            search_terms = [name]
            if cas_number:
                search_terms.append(cas_number)
                
            for term in search_terms:
                # Try direct page
                url = f"https://psychonautwiki.org/wiki/{quote(term)}"
                response = self.http.make_request(url)
                if not response:
                    continue
                    
                soup = BeautifulSoup(response.text, 'html.parser')
                
                # Extract each section
                for section_type, patterns in self.SECTIONS.items():
                    section_data = self._extract_psychonaut_section(
                        soup,
                        patterns,
                        section_type
                    )
                    if section_data:
                        if isinstance(section_data, list):
                            data[section_type].extend(section_data)
                        elif isinstance(section_data, dict):
                            data[section_type].update(section_data)
                            
                # Extract dosage information
                dosage_data = self._extract_psychonaut_dosage(soup)
                if dosage_data:
                    data['dosage'].update(dosage_data)
                    
                # Extract ROA information
                roa_data = self._extract_psychonaut_roa(soup)
                if roa_data:
                    data['roa'].update(roa_data)
                    
                # Use LLM for additional analysis if available
                if llm_api_key:
                    llm_data = analyze_content_with_llm(
                        str(soup),
                        {'name': name},
                        llm_api_key
                    )
                    if llm_data:
                        self._merge_llm_data(data, llm_data)
                        
                if any(v for v in data.values()):
                    return data
                    
        except Exception as e:
            logger.error(f"Error getting PsychonautWiki data: {str(e)}")
            
        return None

    def _extract_psychonaut_section(
        self,
        soup: BeautifulSoup,
        patterns: List[str],
        section_type: str
    ) -> Optional[Any]:
        """Extract section data from PsychonautWiki."""
        try:
            for pattern in patterns:
                section = soup.find('span', {'id': re.compile(pattern, re.I)})
                if section:
                    content_div = section.find_next('div')
                    if not content_div:
                        continue
                        
                    if section_type in ['dosage', 'roa']:
                        # Extract structured data
                        return self._extract_structured_data(content_div, section_type)
                    elif section_type == 'effects':
                        # Extract categorized effects
                        return self._extract_effects(content_div)
                    elif section_type == 'pharmacology':
                        # Extract receptor activity
                        return self._extract_receptor_activity(content_div)
                    elif section_type == 'duration':
                        # Extract duration information
                        return self._extract_duration(content_div)
                    elif section_type == 'tolerance':
                        # Extract tolerance information
                        return self._extract_tolerance(content_div)
                    else:
                        # Extract text content
                        paragraphs = []
                        for p in content_div.find_all('p'):
                            text = p.text.strip()
                            if text and len(text) > 20:  # Skip short lines
                                paragraphs.append(text)
                        return paragraphs
                        
        except Exception as e:
            logger.error(f"Error extracting PsychonautWiki section: {str(e)}")
            
        return None

    def _extract_effects(self, element: BeautifulSoup) -> Dict[str, List[Dict[str, Any]]]:
        """Extract categorized effects."""
        effects = {category: [] for category in self.EFFECT_CATEGORIES}
        
        try:
            # First try to find effects in lists/tables
            for category, patterns in self.EFFECT_CATEGORIES.items():
                # Look for category headers
                for pattern in patterns:
                    headers = element.find_all(
                        ['h3', 'h4', 'strong'],
                        string=re.compile(pattern, re.I)
                    )
                    
                    for header in headers:
                        # Get the list that follows
                        effect_list = header.find_next('ul')
                        if effect_list:
                            for item in effect_list.find_all('li'):
                                effect_text = item.text.strip()
                                if effect_text:
                                    # Look for intensity/duration in parentheses
                                    intensity = None
                                    duration = None
                                    notes = None
                                    
                                    # Extract intensity
                                    intensity_match = re.search(
                                        r'\((mild|moderate|strong|intense|extreme)\)',
                                        effect_text,
                                        re.I
                                    )
                                    if intensity_match:
                                        intensity = intensity_match.group(1).lower()
                                        effect_text = re.sub(
                                            r'\([^)]*\)',
                                            '',
                                            effect_text
                                        ).strip()
                                    
                                    # Extract duration
                                    duration_match = re.search(
                                        r'\((\d+(?:-\d+)?\s*(?:min|hour|day)s?)\)',
                                        effect_text,
                                        re.I
                                    )
                                    if duration_match:
                                        duration = duration_match.group(1).lower()
                                        effect_text = re.sub(
                                            r'\([^)]*\)',
                                            '',
                                            effect_text
                                        ).strip()
                                    
                                    # Extract notes
                                    notes_match = re.search(
                                        r'\((.*?)\)',
                                        effect_text
                                    )
                                    if notes_match:
                                        notes = notes_match.group(1).strip()
                                        effect_text = re.sub(
                                            r'\([^)]*\)',
                                            '',
                                            effect_text
                                        ).strip()
                                    
                                    effects[category].append({
                                        'effect': effect_text,
                                        'intensity': intensity,
                                        'duration': duration,
                                        'notes': notes
                                    })
                                    
            # Then look for effects in paragraphs
            for p in element.find_all('p'):
                text = p.text.strip().lower()
                if len(text) > 20:  # Skip short lines
                    # Try to categorize the effect
                    for category, patterns in self.EFFECT_CATEGORIES.items():
                        for pattern in patterns:
                            if re.search(pattern, text, re.I):
                                effects[category].append({
                                    'effect': text,
                                    'intensity': None,
                                    'duration': None,
                                    'notes': None
                                })
                                break
                                
        except Exception as e:
            logger.error(f"Error extracting effects: {str(e)}")
            
        return effects


    def _extract_duration(self, element: BeautifulSoup) -> Dict[str, Dict[str, str]]:
        """Extract duration information."""
        duration_data = {}
        try:
            # Look for duration tables
            tables = element.find_all('table')
            for table in tables:
                roa = None
                # Try to find ROA from table header or previous heading
                header = table.find_previous(['h3', 'h4'])
                if header:
                    roa = header.text.strip().lower()
                
                if not roa:
                    # Try to find ROA from first row
                    first_row = table.find('tr')
                    if first_row:
                        first_cell = first_row.find('td')
                        if first_cell:
                            roa = first_cell.text.strip().lower()
                
                if roa:
                    duration_data[roa] = {
                        'onset': None,
                        'comeup': None,
                        'peak': None,
                        'offset': None,
                        'after_effects': None,
                        'total': None
                    }
                    
                    # Extract duration values
                    for row in table.find_all('tr'):
                        cols = row.find_all(['td', 'th'])
                        if len(cols) >= 2:
                            phase = cols[0].text.strip().lower()
                            value = cols[1].text.strip()
                            
                            # Map phase names
                            phase_map = {
                                'onset': 'onset',
                                'come up': 'comeup',
                                'peak': 'peak',
                                'offset': 'offset',
                                'after effects': 'after_effects',
                                'total': 'total'
                            }
                            
                            for key, mapped_key in phase_map.items():
                                if key in phase:
                                    duration_data[roa][mapped_key] = value
                                    break
                                    
        except Exception as e:
            logger.error(f"Error extracting duration data: {str(e)}")
            
        return duration_data

    def _extract_tolerance(self, element: BeautifulSoup) -> List[Dict[str, Any]]:
        """Extract tolerance information."""
        tolerance_data = []
        try:
            # Look for tolerance patterns
            patterns = {
                'onset': [
                    r'tolerance\s+(?:onset|develops?|builds?)',
                    r'(?:develops?|builds?)\s+tolerance'
                ],
                'duration': [
                    r'tolerance\s+(?:duration|lasts?|persists?)',
                    r'tolerance\s+for\s+(\d+[\s-]*(?:day|week|month))'
                ],
                'reduction': [
                    r'tolerance\s+(?:reduction|decrease|diminish)',
                    r'(?:reduce|decrease|reset)\s+tolerance'
                ],
                'cross_tolerance': [
                    r'cross[- ]tolerance',
                    r'tolerance\s+(?:with|to)\s+other'
                ]
            }
            
            text = element.text.lower()
            
            # Extract tolerance information
            for tolerance_type, type_patterns in patterns.items():
                for pattern in type_patterns:
                    matches = re.finditer(pattern, text, re.I)
                    for match in matches:
                        # Get surrounding context
                        sentence = self._get_surrounding_sentence(text, match.start())
                        if sentence:
                            tolerance_data.append({
                                'type': tolerance_type,
                                'description': sentence.strip()
                            })
                            
            # Look for specific time periods
            time_patterns = [
                r'(\d+[\s-]*(?:day|week|month))s?\s+(?:of|to|for|until)',
                r'(?:after|within|in)\s+(\d+[\s-]*(?:day|week|month))s?'
            ]
            
            for pattern in time_patterns:
                matches = re.finditer(pattern, text, re.I)
                for match in matches:
                    sentence = self._get_surrounding_sentence(text, match.start())
                    if sentence:
                        tolerance_data.append({
                            'type': 'timeline',
                            'period': match.group(1),
                            'description': sentence.strip()
                        })
                        
        except Exception as e:
            logger.error(f"Error extracting tolerance data: {str(e)}")
            
        return tolerance_data


    def _extract_psychonaut_dosage(
        self,
        soup: BeautifulSoup
    ) -> Optional[Dict[str, Any]]:
        """Extract dosage information from PsychonautWiki."""
        dosage_data = {}
        try:
            dosage_table = soup.find('table', {'class': 'dosage-table'})
            if dosage_table:
                for row in dosage_table.find_all('tr')[1:]:  # Skip header
                    cols = row.find_all('td')
                    if len(cols) >= 3:
                        roa = cols[0].text.strip().lower()
                        if roa not in dosage_data:
                            dosage_data[roa] = {}
                            
                        dose_type = cols[1].text.strip().lower()
                        value = cols[2].text.strip()
                        dosage_data[roa][dose_type] = value
                        
        except Exception as e:
            logger.error(f"Error extracting PsychonautWiki dosage: {str(e)}")
            
        return dosage_data

    def _extract_psychonaut_roa(
        self,
        soup: BeautifulSoup
    ) -> Optional[Dict[str, Any]]:
        """Extract route of administration data from PsychonautWiki."""
        roa_data = {}
        try:
            roa_section = soup.find('span', {'id': re.compile('routes?.*administration', re.I)})
            if roa_section:
                content_div = roa_section.find_next('div')
                if content_div:
                    for roa_div in content_div.find_all('div', recursive=False):
                        roa_title = roa_div.find('h4')
                        if roa_title:
                            roa = roa_title.text.strip().lower()
                            roa_data[roa] = {
                                'bioavailability': None,
                                'onset': None,
                                'duration': None,
                                'after_effects': None
                            }
                            
                            # Extract details
                            details = roa_div.find_all('p')
                            for detail in details:
                                text = detail.text.lower()
                                if 'bioavailability' in text:
                                    roa_data[roa]['bioavailability'] = self._extract_value(text)
                                elif 'onset' in text:
                                    roa_data[roa]['onset'] = self._extract_value(text)
                                elif 'duration' in text:
                                    roa_data[roa]['duration'] = self._extract_value(text)
                                elif 'after' in text:
                                    roa_data[roa]['after_effects'] = self._extract_value(text)
                                    
        except Exception as e:
            logger.error(f"Error extracting PsychonautWiki ROA: {str(e)}")
            
        return roa_data

    def _get_erowid_data(
        self,
        name: str,
        cas_number: Optional[str],
        llm_api_key: Optional[str]
    ) -> Optional[Dict[str, Any]]:
        """Get comprehensive data from Erowid."""
        data = {
            'chemistry': [],
            'pharmacology': [],
            'effects': [],
            'toxicity': [],
            'dosage': {},
            'roa': {},
            'experience_reports': []
        }
        
        try:
            search_terms = [name]
            if cas_number:
                search_terms.append(cas_number)
                
            for term in search_terms:
                # Try search page first
                search_url = f"https://erowid.org/search.php?q={quote(term)}"
                response = self.http.make_request(search_url, verify=False)
                if not response:
                    continue
                    
                soup = BeautifulSoup(response.text, 'html.parser')
                
                # Find main compound page
                compound_link = soup.find('a', href=re.compile(r'/chemicals/[^/]+/[^/]+\.shtml'))
                if not compound_link:
                    continue
                    
                # Get compound page
                compound_url = urljoin('https://erowid.org', compound_link['href'])
                response = self.http.make_request(compound_url, verify=False)
                if not response:
                    continue
                    
                soup = BeautifulSoup(response.text, 'html.parser')
                
                # Extract each section
                for section_type, patterns in self.SECTIONS.items():
                    section_data = self._extract_erowid_section(
                        soup,
                        patterns,
                        section_type
                    )
                    if section_data:
                        if isinstance(section_data, list):
                            data[section_type].extend(section_data)
                        elif isinstance(section_data, dict):
                            data[section_type].update(section_data)
                            
                # Get experience reports
                reports = self._get_erowid_reports(compound_url)
                if reports:
                    data['experience_reports'].extend(reports)
                    
                # Use LLM for additional analysis if available
                if llm_api_key:
                    llm_data = analyze_content_with_llm(
                        str(soup),
                        {'name': name},
                        llm_api_key
                    )
                    if llm_data:
                        self._merge_llm_data(data, llm_data)
                        
                if any(v for v in data.values()):
                    return data
                    
        except Exception as e:
            logger.error(f"Error getting Erowid data: {str(e)}")
            
        return None

    def _extract_erowid_section(
        self,
        soup: BeautifulSoup,
        patterns: List[str],
        section_type: str
    ) -> Optional[Any]:
        """Extract section data from Erowid."""
        try:
            for pattern in patterns:
                section = soup.find('div', {'class': re.compile(f'.*{pattern}.*', re.I)})
                if not section:
                    section = soup.find('h2', string=re.compile(pattern, re.I))
                    if section:
                        section = section.find_next('div')
                        
                if section:
                    if section_type in ['dosage', 'roa']:
                        # Extract structured data
                        return self._extract_structured_data(section, section_type)
                    else:
                        # Extract text content
                        paragraphs = []
                        for p in section.find_all('p'):
                            text = p.text.strip()
                            if text and len(text) > 20:  # Skip short lines
                                paragraphs.append(text)
                        return paragraphs
                        
        except Exception as e:
            logger.error(f"Error extracting Erowid section: {str(e)}")
            
        return None

    def _get_erowid_reports(self, base_url: str) -> List[Dict[str, Any]]:
        """Get experience reports from Erowid."""
        reports = []
        try:
            # Get reports index page
            reports_url = base_url.replace('index.shtml', 'experiences/')
            response = self.http.make_request(reports_url, verify=False)
            if not response:
                return []
                
            soup = BeautifulSoup(response.text, 'html.parser')
            
            # Get report links
            for report_link in soup.find_all('a', href=re.compile(r'/exp/\d+\.shtml')):
                try:
                    report_url = urljoin('https://erowid.org', report_link['href'])
                    response = self.http.make_request(report_url, verify=False)
                    if not response:
                        continue
                        
                    report_soup = BeautifulSoup(response.text, 'html.parser')
                    
                    # Extract report data
                    report = {
                        'title': report_soup.find('h1').text.strip(),
                        'date': None,
                        'author': None,
                        'substance': None,
                        'dose': None,
                        'route': None,
                        'content': [],
                        'url': report_url
                    }
                    
                    # Get metadata
                    metadata = report_soup.find('div', {'class': 'report-metadata'})
                    if metadata:
                        for item in metadata.find_all('div'):
                            text = item.text.lower()
                            if 'date' in text:
                                report['date'] = item.text.split(':')[1].strip()
                            elif 'author' in text:
                                report['author'] = item.text.split(':')[1].strip()
                            elif 'substance' in text:
                                report['substance'] = item.text.split(':')[1].strip()
                            elif 'dose' in text:
                                report['dose'] = item.text.split(':')[1].strip()
                            elif 'route' in text:
                                report['route'] = item.text.split(':')[1].strip()
                                
                    # Get report content
                    content_div = report_soup.find('div', {'class': 'report-text'})
                    if content_div:
                        for p in content_div.find_all('p'):
                            text = p.text.strip()
                            if text:
                                report['content'].append(text)
                                
                    reports.append(report)
                    
                    if len(reports) >= 5:  # Limit to 5 reports
                        break
                        
                except Exception as e:
                    logger.error(f"Error processing report {report_url}: {str(e)}")
                    continue
                    
        except Exception as e:
            logger.error(f"Error getting Erowid reports: {str(e)}")
            
        return reports

    def _extract_structured_data(
        self,
        element: BeautifulSoup,
        data_type: str
    ) -> Dict[str, Any]:
        """Extract structured data from HTML element."""
        data = {}
        try:
            if data_type == 'dosage':
                # Look for dosage tables or lists
                tables = element.find_all('table')
                for table in tables:
                    for row in table.find_all('tr'):
                        cols = row.find_all(['td', 'th'])
                        if len(cols) >= 2:
                            key = cols[0].text.strip().lower()
                            value = cols[1].text.strip()
                            data[key] = value
                            
            elif data_type == 'roa':
                # Look for ROA information
                for p in element.find_all('p'):
                    text = p.text.lower()
                    if ':' in text:
                        key, value = text.split(':', 1)
                        key = key.strip()
                        value = value.strip()
                        if key and value:
                            data[key] = value
                            
        except Exception as e:
            logger.error(f"Error extracting structured data: {str(e)}")
            
        return data

    def _extract_value(self, text: str) -> Optional[str]:
        """Extract value from text."""
        try:
            value = text.split(':', 1)[1].strip()
            return value if value else None
        except (IndexError, AttributeError) as e:
            logger.error(f"Error extracting value from text: {str(e)}")
            return None

    def _merge_data(
        self,
        target: Dict[str, Any],
        source: Dict[str, Any]
    ) -> None:
        """Merge source data into target data."""
        for key, value in source.items():
            if isinstance(value, list):
                if key not in target:
                    target[key] = []
                target[key].extend(value)
            elif isinstance(value, dict):
                if key not in target:
                    target[key] = {}
                target[key].update(value)
            else:
                target[key] = value

    def _merge_llm_data(
        self,
        target: Dict[str, Any],
        llm_data: Dict[str, Any]
    ) -> None:
        """Merge LLM-extracted data into target data."""
        # Map LLM fields to our fields
        field_map = {
            'chemical_properties': 'chemistry',
            'pharmacological_properties': 'pharmacology',
            'effects': 'effects',
            'safety_data': 'toxicity'
        }
        
        for llm_field, our_field in field_map.items():
            if llm_field in llm_data:
                if isinstance(llm_data[llm_field], list):
                    if our_field not in target:
                        target[our_field] = []
                    target[our_field].extend(llm_data[llm_field])
                elif isinstance(llm_data[llm_field], dict):
                    if our_field not in target:
                        target[our_field] = {}
                    target[our_field].update(llm_data[llm_field])


    def _extract_receptor_activity(self, element: BeautifulSoup) -> Dict[str, Dict[str, List[str]]]:
        """Extract receptor activity information."""
        activities = {
            subtype: {
                'agonist': [],
                'antagonist': [],
                'modulator': [],
                'reuptake': [],
                'releaser': []
            }
            for subtype in self.RECEPTOR_SUBTYPES
        }
        
        try:
            text = element.text.lower()
            
            # First identify receptor subtype context
            for subtype, patterns in self.RECEPTOR_SUBTYPES.items():
                for pattern in patterns:
                    matches = re.finditer(pattern, text, re.I)
                    for match in matches:
                        # Get surrounding text
                        sentence = self._get_surrounding_sentence(text, match.start())
                        if not sentence:
                            continue
                            
                        # Look for activity patterns in this context
                        for activity_type, act_patterns in self.RECEPTOR_PATTERNS.items():
                            for act_pattern in act_patterns:
                                if re.search(act_pattern, sentence, re.I):
                                    activities[subtype][activity_type].append(sentence.strip())
                                    
        except Exception as e:
            logger.error(f"Error extracting receptor activity: {str(e)}")
            
        return activities

    def _get_surrounding_sentence(self, text: str, pos: int) -> Optional[str]:
        """Get the sentence containing the given position."""
        try:
            # Find sentence boundaries
            start = text.rfind('.', 0, pos) + 1
            end = text.find('.', pos)
            
            if start == -1:
                start = 0
            if end == -1:
                end = len(text)
                
            return text[start:end].strip()
            
        except Exception as e:
            logger.error(f"Error getting surrounding sentence: {str(e)}")
            return None

    def get_urls(
        self,
        name: str,
        cas_number: Optional[str]
    ) -> Dict[str, str]:
        """
        Get URLs for community sources.
        
        Args:
            name: Compound name
            cas_number: Optional CAS number
            
        Returns:
            Dictionary of source URLs
        """
        urls = {}
        
        # Try PsychonautWiki
        try:
            url = f"https://psychonautwiki.org/wiki/{quote(name)}"
            response = self.http.make_request(url)
            if response and 'search' not in response.url:
                urls['psychonaut_url'] = response.url
        except Exception as e:
            logger.error(f"Error getting PsychonautWiki URL: {str(e)}")
            
        # Try Erowid
        try:
            # Try search first
            search_url = f"https://erowid.org/search.php?q={quote(name)}"
            response = self.http.make_request(search_url, verify=False)
            if response:
                soup = BeautifulSoup(response.text, 'html.parser')
                compound_link = soup.find('a', href=re.compile(r'/chemicals/[^/]+/[^/]+\.shtml'))
                if compound_link:
                    urls['erowid_url'] = urljoin('https://erowid.org', compound_link['href'])
        except Exception as e:
            logger.error(f"Error getting Erowid URL: {str(e)}")
            
        return urls
