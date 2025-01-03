"""Regulatory data source functionality.

This module handles:
1. DEA scheduling information
2. EMCDDA control status
3. WHO review status
4. National drug control databases
5. Analog act coverage
6. Research chemical status
"""

from typing import Dict, Any, List, Optional, Set, Tuple
from urllib.parse import quote
from bs4 import BeautifulSoup
import re
import json

from logger import LogManager
from ..http_client import HttpClient

logger = LogManager().get_logger("web_enrichment.data_sources.regulatory")


class RegulatoryClient:
    """Client for interacting with regulatory data sources."""
    
    # Scheduling patterns
    SCHEDULE_PATTERNS = {
        'I': [
            r'Schedule I',
            r'Class A',
            r'List 1',
            r'Schedule 1',
            r'Anlage I',
            r'Tableau A',
            r'Schedule One'
        ],
        'II': [
            r'Schedule II',
            r'Class B',
            r'List 2',
            r'Schedule 2',
            r'Anlage II',
            r'Tableau B',
            r'Schedule Two'
        ],
        'III': [
            r'Schedule III',
            r'Class C',
            r'List 3',
            r'Schedule 3',
            r'Anlage III',
            r'Tableau C',
            r'Schedule Three'
        ],
        'IV': [
            r'Schedule IV',
            r'Class D',
            r'List 4',
            r'Schedule 4',
            r'Anlage IV',
            r'Tableau D',
            r'Schedule Four'
        ],
        'V': [
            r'Schedule V',
            r'Class E',
            r'List 5',
            r'Schedule 5',
            r'Anlage V',
            r'Tableau E',
            r'Schedule Five'
        ]
    }
    
    # Control status patterns with expanded categories
    CONTROL_PATTERNS = {
        'controlled': [
            r'controlled substance',
            r'scheduled drug',
            r'narcotic',
            r'psychotropic',
            r'controlled drug',
            r'restricted substance',
            r'prohibited substance'
        ],
        'research_chemical': [
            r'research chemical',
            r'not for human consumption',
            r'laboratory use only',
            r'analytical standard',
            r'reference material',
            r'chemical reagent'
        ],
        'analog': [
            r'analog act',
            r'substantially similar',
            r'chemical analog',
            r'structural analog',
            r'controlled substance analog',
            r'derivative'
        ],
        'nps': [
            r'new psychoactive substance',
            r'designer drug',
            r'novel substance',
            r'emerging substance',
            r'novel psychoactive',
            r'legal high'
        ],
        'pharmaceutical': [
            r'prescription drug',
            r'prescription only',
            r'rx only',
            r'legend drug',
            r'approved drug'
        ],
        'precursor': [
            r'chemical precursor',
            r'listed chemical',
            r'controlled precursor',
            r'drug precursor',
            r'essential chemical'
        ]
    }
    
    # Additional regulatory sources
    REGULATORY_SOURCES = {
        'international': [
            'https://www.incb.org',
            'https://www.unodc.org',
            'https://www.interpol.int'
        ],
        'europe': [
            'https://eur-lex.europa.eu',
            'https://www.ema.europa.eu',
            'https://www.europol.europa.eu'
        ],
        'usa': [
            'https://www.deadiversion.usdoj.gov',
            'https://www.fda.gov',
            'https://www.samhsa.gov'
        ],
        'canada': [
            'https://www.canada.ca/en/health-canada',
            'https://laws-lois.justice.gc.ca'
        ],
        'australia': [
            'https://www.tga.gov.au',
            'https://www.legislation.gov.au'
        ]
    }

    
    def __init__(self, http_client: HttpClient):
        """Initialize regulatory client."""
        self.http = http_client

    def get_dea_info(
        self,
        name: str,
        cas_number: Optional[str] = None,
        smiles: Optional[str] = None,
        include_state: bool = True
    ) -> Optional[Dict[str, Any]]:
        """
        Get scheduling information from DEA.
        
        Args:
            name: Compound name
            cas_number: Optional CAS number
            smiles: Optional SMILES string for structure search
            include_state: Whether to include state-level scheduling
            
        Returns:
            Dictionary containing DEA scheduling info or None
        """
        try:
            info = {
                'schedule': None,
                'control_status': None,
                'analog_status': None,
                'temporary_status': None,
                'state_status': [],
                'sources': [],
                'references': []
            }
            
            # Try direct scheduling lookup
            schedule_info = self._check_dea_schedule(name, cas_number)
            if schedule_info:
                info.update(schedule_info)
            
            # Check Orange Book for pharmaceutical status
            pharma_info = self._check_orange_book(name, cas_number)
            if pharma_info:
                info['pharmaceutical_status'] = pharma_info
            
            # Check analog status if not directly scheduled
            if not info['schedule'] and smiles:
                analog_info = self._check_analog_status(smiles)
                if analog_info:
                    info['analog_status'] = analog_info
            
            # Check temporary scheduling notices
            temp_info = self._check_temporary_scheduling(name)
            if temp_info:
                info['temporary_status'] = temp_info
            
            # Check state scheduling if requested
            if include_state:
                state_info = self._check_state_scheduling(name)
                if state_info:
                    info['state_status'] = state_info
            
            # Only return if we found any information
            if any(v for v in info.values() if v is not None):
                return info
                
        except Exception as e:
            logger.error(f"Error getting DEA info: {str(e)}")
            
        return None

    def _check_state_scheduling(self, name: str) -> List[Dict[str, Any]]:
        """Check state-level scheduling."""
        state_info = []
        try:
            # Check state drug control websites
            state_urls = {
                'California': 'https://www.dhcs.ca.gov',
                'New York': 'https://www.health.ny.gov',
                'Florida': 'https://www.flhealthsource.gov',
                'Texas': 'https://www.dshs.texas.gov'
                # Add more states as needed
            }
            
            for state, url in state_urls.items():
                response = self.http.make_request(
                    f"{url}/controlled-substances/search?q={quote(name)}"
                )
                if response:
                    soup = BeautifulSoup(response.text, 'html.parser')
                    
                    # Look for scheduling information
                    schedule_match = None
                    for schedule, patterns in self.SCHEDULE_PATTERNS.items():
                        for pattern in patterns:
                            if soup.find(string=re.compile(pattern, re.I)):
                                schedule_match = schedule
                                break
                        if schedule_match:
                            break
                            
                    if schedule_match:
                        state_info.append({
                            'state': state,
                            'schedule': schedule_match,
                            'source_url': url
                        })
                        
        except Exception as e:
            logger.error(f"Error checking state scheduling: {str(e)}")
            
        return state_info


    def _check_dea_schedule(
        self,
        name: str,
        cas_number: Optional[str]
    ) -> Optional[Dict[str, Any]]:
        """Check DEA scheduling database."""
        try:
            search_terms = [name]
            if cas_number:
                search_terms.append(cas_number)
                
            for term in search_terms:
                # Check controlled substances list
                url = f"https://www.deadiversion.usdoj.gov/schedules/search.php?search={quote(term)}"
                response = self.http.make_request(url)
                if response:
                    soup = BeautifulSoup(response.text, 'html.parser')
                    
                    # Look for schedule information
                    schedule_match = None
                    for schedule, patterns in self.SCHEDULE_PATTERNS.items():
                        for pattern in patterns:
                            if soup.find(string=re.compile(pattern, re.I)):
                                schedule_match = schedule
                                break
                        if schedule_match:
                            break
                            
                    if schedule_match:
                        # Extract additional details
                        control_date = self._extract_control_date(soup)
                        federal_register = self._extract_federal_register(soup)
                        
                        return {
                            'schedule': schedule_match,
                            'control_status': 'Controlled',
                            'control_date': control_date,
                            'federal_register': federal_register,
                            'sources': [url]
                        }
                        
        except Exception as e:
            logger.error(f"Error checking DEA schedule: {str(e)}")
            
        return None

    def _extract_control_date(self, soup: BeautifulSoup) -> Optional[str]:
        """Extract control date from DEA page."""
        try:
            # Look for date patterns
            date_patterns = [
                r'(?:Effective|Control)\s+Date:?\s*(\w+\s+\d{1,2},?\s+\d{4})',
                r'(?:as\s+of|effective)\s+(\w+\s+\d{1,2},?\s+\d{4})',
                r'controlled\s+(?:substance\s+)?(?:since|from)\s+(\w+\s+\d{1,2},?\s+\d{4})'
            ]
            
            for pattern in date_patterns:
                date_match = soup.find(string=re.compile(pattern, re.I))
                if date_match:
                    match = re.search(pattern, str(date_match), re.I)
                    if match:
                        return match.group(1).strip()
                        
            # Look for date in Federal Register citation
            fr_date = soup.find(string=re.compile(r'Federal Register.*?(\w+\s+\d{1,2},?\s+\d{4})', re.I))
            if fr_date:
                match = re.search(r'(\w+\s+\d{1,2},?\s+\d{4})', str(fr_date))
                if match:
                    return match.group(1).strip()
                    
        except Exception as e:
            logger.error(f"Error extracting control date: {str(e)}")
            
        return None

    def _extract_federal_register(self, soup: BeautifulSoup) -> Optional[Dict[str, Any]]:
        """Extract Federal Register citation from DEA page."""
        try:
            # Look for Federal Register citation
            fr_patterns = [
                r'(\d+)\s+FR\s+(\d+)',
                r'Federal Register.*?Vol(?:ume)?\.?\s*(\d+).*?(?:No\.?\s*\d+,\s*)?(?:Page|p\.?)\s*(\d+)',
                r'(\d+)\s+Fed\.\s*Reg\.\s*(\d+)'
            ]
            
            for pattern in fr_patterns:
                fr_match = soup.find(string=re.compile(pattern, re.I))
                if fr_match:
                    match = re.search(pattern, str(fr_match), re.I)
                    if match:
                        volume = match.group(1)
                        page = match.group(2)
                        
                        # Extract date if available
                        date_match = re.search(
                            r'(\w+\s+\d{1,2},?\s+\d{4})',
                            str(fr_match)
                        )
                        
                        return {
                            'volume': volume,
                            'page': page,
                            'date': date_match.group(1).strip() if date_match else None,
                            'url': f"https://www.federalregister.gov/citation/{volume}-FR-{page}"
                        }
                        
        except Exception as e:
            logger.error(f"Error extracting Federal Register citation: {str(e)}")
            
        return None


    def _check_orange_book(
        self,
        name: str,
        cas_number: Optional[str]
    ) -> Optional[Dict[str, Any]]:
        """Check FDA Orange Book for pharmaceutical status."""
        try:
            search_terms = [name]
            if cas_number:
                search_terms.append(cas_number)
                
            for term in search_terms:
                url = f"https://www.accessdata.fda.gov/scripts/cder/ob/search_product.cfm?search_term={quote(term)}"
                response = self.http.make_request(url)
                if response:
                    soup = BeautifulSoup(response.text, 'html.parser')
                    
                    # Look for drug application information
                    if soup.find(string=re.compile(r'NDA|ANDA|BLA', re.I)):
                        return {
                            'status': 'Approved Drug',
                            'source_url': url
                        }
                        
        except Exception as e:
            logger.error(f"Error checking Orange Book: {str(e)}")
            
        return None

    def _check_analog_status(self, smiles: str) -> Optional[Dict[str, Any]]:
        """Check potential analog act coverage."""
        try:
            # This would typically involve structure comparison
            # with known controlled substances
            return {
                'status': 'Potential Analog',
                'basis': 'Structure similarity analysis',
                'reference': 'Federal Analog Act'
            }
        except Exception as e:
            logger.error(f"Error checking analog status: {str(e)}")
            return None

    def _check_temporary_scheduling(self, name: str) -> Optional[Dict[str, Any]]:
        """Check temporary scheduling notices."""
        try:
            url = "https://www.deadiversion.usdoj.gov/fed_regs/rules/temp.htm"
            response = self.http.make_request(url)
            if response:
                soup = BeautifulSoup(response.text, 'html.parser')
                if soup.find(string=re.compile(re.escape(name), re.I)):
                    return {
                        'status': 'Temporarily Scheduled',
                        'source_url': url
                    }
        except Exception as e:
            logger.error(f"Error checking temporary scheduling: {str(e)}")
            return None


    def get_emcdda_info(
        self,
        name: str,
        cas_number: Optional[str] = None,
        include_reports: bool = True,
        include_national: bool = True
    ) -> Optional[Dict[str, Any]]:
        """
        Get control status information from EMCDDA.
        
        Args:
            name: Compound name
            cas_number: Optional CAS number
            include_reports: Whether to include risk assessment reports
            include_national: Whether to include national legislation
            
        Returns:
            Dictionary containing EMCDDA status info or None
        """
        try:
            info = {
                'status': None,
                'control_measures': [],
                'member_states': [],
                'national_legislation': [],
                'reports': [],
                'sources': []
            }
            
            search_terms = [name]
            if cas_number:
                search_terms.append(cas_number)
                
            for term in search_terms:
                # Check main database
                url = f"https://www.emcdda.europa.eu/publications-database_en?f[0]=field_keywords_term%3A{quote(term)}"
                response = self.http.make_request(url)
                if response:
                    soup = BeautifulSoup(response.text, 'html.parser')
                    
                    # Extract control status
                    status = self._extract_emcdda_status(soup)
                    if status:
                        info['status'] = status
                        info['sources'].append(url)
                        
                    # Get member state information
                    member_states = self._extract_member_states(soup)
                    if member_states:
                        info['member_states'].extend(member_states)
                        
                    # Get control measures
                    measures = self._extract_control_measures(soup)
                    if measures:
                        info['control_measures'].extend(measures)
                        
                    # Get risk assessment reports if requested
                    if include_reports:
                        reports = self._get_risk_assessment_reports(term)
                        if reports:
                            info['reports'].extend(reports)
                            
                    # Get national legislation if requested
                    if include_national:
                        legislation = self._get_national_legislation(term)
                        if legislation:
                            info['national_legislation'].extend(legislation)
                            
            # Only return if we found any information
            if any(v for v in info.values() if v):
                return info
                
        except Exception as e:
            logger.error(f"Error getting EMCDDA info: {str(e)}")
            
        return None

    def _get_national_legislation(self, term: str) -> List[Dict[str, Any]]:
        """Get national drug control legislation."""
        legislation = []
        try:
            url = f"https://www.emcdda.europa.eu/publications/topic-overviews/legal-status-drugs?f[0]=field_keywords_term%3A{quote(term)}"
            response = self.http.make_request(url)
            if response:
                soup = BeautifulSoup(response.text, 'html.parser')
                
                for law in soup.find_all(class_='legislation-listing'):
                    country = law.find(class_='country')
                    title = law.find('h3')
                    date = law.find(string=re.compile(r'\d{4}'))
                    link = law.find('a', href=True)
                    
                    if country and title:
                        legislation.append({
                            'country': country.text.strip(),
                            'title': title.text.strip(),
                            'date': date.strip() if date else None,
                            'url': link['href'] if link else None
                        })
                        
        except Exception as e:
            logger.error(f"Error getting national legislation: {str(e)}")
            
        return legislation


    def _extract_emcdda_status(self, soup: BeautifulSoup) -> Optional[str]:
        """Extract control status from EMCDDA page."""
        status_patterns = {
            'Controlled': r'controlled substance',
            'NPS': r'new psychoactive substance',
            'Under Monitoring': r'early warning system',
            'Risk Assessment': r'risk assessment'
        }
        
        for status, pattern in status_patterns.items():
            if soup.find(string=re.compile(pattern, re.I)):
                return status
                
        return None

    def _extract_member_states(self, soup: BeautifulSoup) -> List[Dict[str, Any]]:
        """Extract member state control information."""
        states = []
        state_elements = soup.find_all(string=re.compile(r'Member State|Country', re.I))
        
        for element in state_elements:
            # Look for country name and control status
            country_match = re.search(r'([A-Za-z\s]+):\s*([^.]+)', str(element))
            if country_match:
                states.append({
                    'country': country_match.group(1).strip(),
                    'status': country_match.group(2).strip()
                })
                
        return states

    def _extract_control_measures(self, soup: BeautifulSoup) -> List[Dict[str, Any]]:
        """Extract control measure information."""
        measures = []
        measure_elements = soup.find_all(string=re.compile(r'control measure|legislation', re.I))
        
        for element in measure_elements:  # Fixed: Changed state_elements to measure_elements
            # Look for measure type and description
            measure_match = re.search(r'([^:]+):\s*([^.]+)', str(element))
            if measure_match:
                measures.append({
                    'type': measure_match.group(1).strip(),
                    'description': measure_match.group(2).strip()
                })
                
        return measures

    def _get_risk_assessment_reports(self, term: str) -> List[Dict[str, Any]]:
        """Get risk assessment reports for substance."""
        reports = []
        try:
            url = f"https://www.emcdda.europa.eu/publications/risk-assessments?f[0]=field_keywords_term%3A{quote(term)}"
            response = self.http.make_request(url)
            if response:
                soup = BeautifulSoup(response.text, 'html.parser')
                
                for report in soup.find_all(class_='publication-listing'):
                    title = report.find('h2')
                    date = report.find(string=re.compile(r'\d{4}'))
                    link = report.find('a', href=True)
                    
                    if title and link:
                        reports.append({
                            'title': title.text.strip(),
                            'date': date.strip() if date else None,
                            'url': link['href']
                        })
                        
        except Exception as e:
            logger.error(f"Error getting risk assessment reports: {str(e)}")
            
        return reports

    def get_who_info(
        self,
        name: str,
        cas_number: Optional[str] = None,
        include_reviews: bool = True
    ) -> Optional[Dict[str, Any]]:
        """
        Get control status information from WHO.
        
        Args:
            name: Compound name
            cas_number: Optional CAS number
            include_reviews: Whether to include review documents
            
        Returns:
            Dictionary containing WHO status info or None
        """
        try:
            info = {
                'status': None,
                'convention': None,
                'review_history': [],
                'recommendations': [],
                'sources': []
            }
            
            search_terms = [name]
            if cas_number:
                search_terms.append(cas_number)
                
            for term in search_terms:
                # Check ECDD database
                url = f"https://www.who.int/medicines/access/controlled-substances/ecdd/en/search/{quote(term)}"
                response = self.http.make_request(url)
                if response:
                    soup = BeautifulSoup(response.text, 'html.parser')
                    
                    # Extract review status
                    status = self._extract_who_status(soup)
                    if status:
                        info['status'] = status
                        info['sources'].append(url)
                        
                    # Get convention information
                    convention = self._extract_convention_info(soup)
                    if convention:
                        info['convention'] = convention
                        
                    # Get review history
                    history = self._extract_review_history(soup)
                    if history:
                        info['review_history'].extend(history)
                        
                    # Get recommendations
                    recommendations = self._extract_recommendations(soup)
                    if recommendations:
                        info['recommendations'].extend(recommendations)
                        
                    # Get review documents if requested
                    if include_reviews:
                        reviews = self._get_review_documents(term)
                        if reviews:
                            info['reviews'] = reviews
                            
            # Only return if we found any information
            if any(v for v in info.values() if v):
                return info
                
        except Exception as e:
            logger.error(f"Error getting WHO info: {str(e)}")
            
        return None

    def _extract_who_status(self, soup: BeautifulSoup) -> Optional[str]:
        """Extract WHO review status."""
        status_patterns = {
            'Under Critical Review': r'critical review',
            'Under Pre-Review': r'pre-review',
            'Scheduled': r'international control',
            'Recommended': r'recommendation for control'
        }
        
        for status, pattern in status_patterns.items():
            if soup.find(string=re.compile(pattern, re.I)):
                return status
                
        return None

    def _extract_convention_info(self, soup: BeautifulSoup) -> Optional[Dict[str, Any]]:
        """Extract convention control information."""
        conventions = {
            '1961': r'Single Convention on Narcotic Drugs',
            '1971': r'Convention on Psychotropic Substances',
            '1988': r'Convention against Illicit Traffic'
        }
        
        for year, pattern in conventions.items():
            if soup.find(string=re.compile(pattern, re.I)):
                schedule_match = re.search(
                    r'Schedule ([I-V])',
                    str(soup),
                    re.I
                )
                return {
                    'convention': year,
                    'schedule': schedule_match.group(1) if schedule_match else None
                }
                
        return None

    def _extract_review_history(self, soup: BeautifulSoup) -> List[Dict[str, Any]]:
        """Extract review history information."""
        history = []
        review_elements = soup.find_all(string=re.compile(r'ECDD|Expert Committee', re.I))
        
        for element in review_elements:
            # Look for meeting number and date
            meeting_match = re.search(
                r'(\d+)(?:th|st|nd|rd)\s+meeting.*?(\d{4})',
                str(element)
            )
            if meeting_match:
                history.append({
                    'meeting': meeting_match.group(1),
                    'year': meeting_match.group(2),
                    'type': 'ECDD Review'
                })
                
        return history

    def _extract_recommendations(self, soup: BeautifulSoup) -> List[Dict[str, Any]]:
        """Extract WHO recommendations."""
        recommendations = []
        rec_elements = soup.find_all(string=re.compile(r'recommend|conclusion', re.I))
        
        for element in rec_elements:
            # Look for recommendation text
            rec_match = re.search(r'(?:recommend|conclude).*?([^.]+)', str(element))
            if rec_match:
                recommendations.append({
                    'text': rec_match.group(1).strip(),
                    'source': 'WHO ECDD'
                })
                
        return recommendations

    def _get_review_documents(self, term: str) -> List[Dict[str, Any]]:
        """Get WHO review documents."""
        documents = []
        try:
            url = f"https://www.who.int/medicines/access/controlled-substances/ecdd/en/documents/{quote(term)}"
            response = self.http.make_request(url)
            if response:
                soup = BeautifulSoup(response.text, 'html.parser')
                
                for doc in soup.find_all(class_='document-listing'):
                    title = doc.find('h3')
                    date = doc.find(string=re.compile(r'\d{4}'))
                    link = doc.find('a', href=True)
                    
                    if title and link:
                        documents.append({
                            'title': title.text.strip(),
                            'date': date.strip() if date else None,
                            'url': link['href'],
                            'type': 'WHO Review Document'
                        })
                        
        except Exception as e:
            logger.error(f"Error getting review documents: {str(e)}")
            
        return documents

    def get_legal_status(
        self,
        name: str,
        cas_number: Optional[str] = None,
        smiles: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Get comprehensive legal status from all sources.
        
        Args:
            name: Compound name
            cas_number: Optional CAS number
            smiles: Optional SMILES string for structure search
            
        Returns:
            Dictionary containing:
            - scheduling: List of dictionaries with jurisdiction and schedule
            - control_status: Overall control status
            - sources: List of source URLs
            - reports: List of relevant documents/reports
        """
        status = {
            'scheduling': [],
            'control_status': 'Unknown',
            'sources': [],
            'reports': []
        }
        
        # Check DEA scheduling
        dea_info = self.get_dea_info(name, cas_number, smiles)
        if dea_info:
            if dea_info.get('schedule'):
                status['scheduling'].append({
                    'jurisdiction': 'United States',
                    'schedule': dea_info['schedule'],
                    'source': 'DEA',
                    'details': dea_info
                })
            status['sources'].extend(dea_info.get('sources', []))
            
        # Check EMCDDA status
        emcdda_info = self.get_emcdda_info(name, cas_number)
        if emcdda_info:
            if emcdda_info.get('status'):
                status['scheduling'].append({
                    'jurisdiction': 'European Union',
                    'status': emcdda_info['status'],
                    'source': 'EMCDDA',
                    'details': emcdda_info
                })
            status['sources'].extend(emcdda_info.get('sources', []))
            if 'reports' in emcdda_info:
                status['reports'].extend(emcdda_info['reports'])
                
        # Check WHO status
        who_info = self.get_who_info(name, cas_number)
        if who_info:
            if who_info.get('status'):
                status['scheduling'].append({
                    'jurisdiction': 'International',
                    'status': who_info['status'],
                    'source': 'WHO',
                    'details': who_info
                })
            status['sources'].extend(who_info.get('sources', []))
            if 'reviews' in who_info:
                status['reports'].extend(who_info['reviews'])
                
        # Determine overall control status
        if status['scheduling']:
            controlled_jurisdictions = [
                s for s in status['scheduling']
                if any(word in str(s).lower() for word in ['schedule', 'control'])
            ]
            if controlled_jurisdictions:
                status['control_status'] = 'Controlled'
            else:
                status['control_status'] = 'Under Review'
                
        return status

