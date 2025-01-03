"""LLM utilities for web enrichment.

This module handles:
1. Content analysis using LLM
2. Patent information extraction
3. Structured data extraction
4. PDF document processing
5. Chemical structure recognition
"""

import json
import re
from typing import Dict, Any, List, Optional, Tuple
import requests
import PyPDF2
from io import BytesIO

from logger import LogManager
from rdkit import Chem
from rdkit.Chem import AllChem

logger = LogManager().get_logger("web_enrichment.llm_utils")


def analyze_content_with_llm(
    content: str,
    compound_data: Dict[str, Any],
    llm_api_key: str
) -> Dict[str, Any]:
    """
    Use LLM to analyze webpage content and extract structured data.
    
    Args:
        content: Webpage text content
        compound_data: Known compound identifiers
        llm_api_key: API key for LLM service
        
    Returns:
        Dictionary of extracted data
    """
    try:
        # Pre-process content to highlight chemical information
        processed_content = _preprocess_chemical_content(content)
        
        # Construct enhanced prompt for LLM
        prompt = f"""Analyze this chemical compound webpage content and extract relevant information:

Known identifiers:
Name: {compound_data.get('name', 'Unknown')}
CAS: {compound_data.get('cas', 'Unknown')}
SMILES: {compound_data.get('smiles', 'Unknown')}
InChI: {compound_data.get('inchi', 'Unknown')}

Webpage content:
{processed_content[:4000]}  # Increased context window

Please extract the following information if present:
1. Alternative names/synonyms (including systematic names, trade names, and common names)
2. Chemical identifiers:
   - CAS numbers
   - SMILES strings
   - InChI strings
   - Registry numbers
3. Chemical properties:
   - Molecular weight
   - Melting point
   - Boiling point
   - LogP
   - pKa
4. Pharmacological properties:
   - Mechanism of action
   - Target receptors
   - Binding affinities
   - Activity type (agonist/antagonist)
   - Efficacy data
5. Safety data:
   - Toxicity (LD50, etc.)
   - Side effects
   - Handling precautions
   - Safety classifications
6. Regulatory status:
   - Scheduling status
   - Legal status by country
   - Control measures
7. References:
   - Scientific papers
   - Patents
   - Database entries

Format the response as a JSON object with these fields."""

        # Make LLM API request with increased max tokens
        response = requests.post(
            "https://api.openai.com/v1/chat/completions",
            headers={
                "Authorization": f"Bearer {llm_api_key}",
                "Content-Type": "application/json"
            },
            json={
                "model": "gpt-4",
                "messages": [{"role": "user", "content": prompt}],
                "temperature": 0.7,
                "max_tokens": 2000
            }
        )
        
        if response.status_code == 200:
            result = response.json()
            try:
                extracted_data = json.loads(
                    result['choices'][0]['message']['content']
                )
                # Post-process and validate extracted data
                validated_data = _validate_chemical_data(extracted_data)
                return validated_data
            except json.JSONDecodeError:
                logger.error("Failed to parse LLM response as JSON")
                
    except Exception as e:
        logger.error(f"Error analyzing content with LLM: {str(e)}")
        
    return {}


def extract_patent_compound(
    example_text: str,
    llm_api_key: str,
    include_synthesis: bool = True
) -> Dict[str, Any]:
    """
    Extract compound information from patent using LLM.
    
    Args:
        example_text: Patent example text
        llm_api_key: API key for LLM service
        include_synthesis: Whether to extract synthesis details
        
    Returns:
        Dictionary containing extracted identifiers and data
    """
    try:
        # Pre-process patent text
        processed_text = _preprocess_patent_text(example_text)
        
        # Enhanced prompt for chemical information
        prompt = f"""Extract detailed chemical information from this patent example:

Example text: {processed_text}

Please provide:
1. Chemical identifiers:
   - CAS number
   - SMILES string
   - InChI string
   - Registry numbers
2. Chemical names:
   - Systematic name (IUPAC)
   - Common names
   - Trade names
3. Structure details:
   - Molecular formula
   - Molecular weight
   - Stereochemistry
   - Structural features
4. Physical properties:
   - Melting point
   - Boiling point
   - Solubility
   - Physical form
5. Preparation method:
   - Synthetic route
   - Reaction conditions
   - Yields
   - Purification
6. Analytical data:
   - NMR data
   - Mass spec data
   - Elemental analysis
7. Biological activity:
   - Target receptors
   - Activity values
   - Assay conditions

Format the response as a JSON object with these fields."""

        # Make LLM API request with increased context
        response = requests.post(
            "https://api.openai.com/v1/chat/completions",
            headers={
                "Authorization": f"Bearer {llm_api_key}",
                "Content-Type": "application/json"
            },
            json={
                "model": "gpt-4",
                "messages": [{"role": "user", "content": prompt}],
                "temperature": 0.7,
                "max_tokens": 2000
            }
        )
        
        if response.status_code == 200:
            result = response.json()
            try:
                extracted_data = json.loads(
                    result['choices'][0]['message']['content']
                )
                # Validate and standardize extracted data
                validated_data = _validate_patent_data(extracted_data)
                return validated_data
            except json.JSONDecodeError:
                logger.error("Failed to parse LLM response as JSON")
                
    except Exception as e:
        logger.error(f"Error extracting patent compound: {str(e)}")
        
    return {}


def process_pdf_content(
    pdf_content: bytes,
    llm_api_key: str
) -> List[Dict[str, Any]]:
    """
    Process PDF content to extract chemical information.
    
    Args:
        pdf_content: PDF file content as bytes
        llm_api_key: API key for LLM service
        
    Returns:
        List of dictionaries containing extracted compound data
    """
    try:
        # Read PDF content
        pdf_file = BytesIO(pdf_content)
        pdf_reader = PyPDF2.PdfReader(pdf_file)
        
        # Extract text from each page
        full_text = ""
        for page in pdf_reader.pages:
            full_text += page.extract_text() + "\n"
            
        # Pre-process text to highlight chemical information
        processed_text = _preprocess_chemical_content(full_text)
        
        # Split into sections for processing
        sections = _split_into_sections(processed_text)
        
        compounds = []
        for section in sections:
            # Extract compounds from each section
            prompt = f"""Extract chemical compound information from this text section:

Text:
{section}

Please identify and extract information about any chemical compounds mentioned, including:
1. Chemical names and identifiers
2. Structural information
3. Properties and characteristics
4. Biological activity
5. Synthesis details if present

Format each compound as a JSON object."""

            # Process with LLM
            response = requests.post(
                "https://api.openai.com/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {llm_api_key}",
                    "Content-Type": "application/json"
                },
                json={
                    "model": "gpt-4",
                    "messages": [{"role": "user", "content": prompt}],
                    "temperature": 0.7,
                    "max_tokens": 2000
                }
            )
            
            if response.status_code == 200:
                result = response.json()
                try:
                    section_compounds = json.loads(
                        result['choices'][0]['message']['content']
                    )
                    if isinstance(section_compounds, list):
                        compounds.extend(section_compounds)
                    elif isinstance(section_compounds, dict):
                        compounds.append(section_compounds)
                except json.JSONDecodeError:
                    continue
                    
        # Deduplicate and validate compounds
        unique_compounds = _deduplicate_compounds(compounds)
        validated_compounds = [
            _validate_chemical_data(comp)
            for comp in unique_compounds
        ]
        
        return validated_compounds
        
    except Exception as e:
        logger.error(f"Error processing PDF content: {str(e)}")
        return []


def _preprocess_chemical_content(text: str) -> str:
    """
    Pre-process text to highlight chemical information.
    
    Args:
        text: Raw text content
        
    Returns:
        Processed text with highlighted chemical information
    """
    # Highlight potential CAS numbers
    text = re.sub(
        r'(\d{1,7}-\d{2}-\d)',
        r'CAS_NUMBER: \1',
        text
    )
    
    # Highlight potential SMILES strings
    text = re.sub(
        r'([CN]\[C@H\].+?(?=\s|$))',
        r'SMILES: \1',
        text
    )
    
    # Highlight potential InChI strings
    text = re.sub(
        r'(InChI=1S?.+?(?=\s|$))',
        r'INCHI: \1',
        text
    )
    
    # Highlight chemical formulas
    text = re.sub(
        r'([A-Z][a-z]?\d*)+',
        r'FORMULA: \1',
        text
    )
    
    return text


def _preprocess_patent_text(text: str) -> str:
    """
    Pre-process patent text for better chemical extraction.
    
    Args:
        text: Raw patent text
        
    Returns:
        Processed text
    """
    # Remove patent boilerplate
    text = re.sub(r'^\s*\[?\d+\]?\s*', '', text)
    
    # Highlight example sections
    text = re.sub(
        r'(?i)(example\s+\d+[.:]\s*)',
        r'\nEXAMPLE_START\1',
        text
    )
    
    # Highlight preparation sections
    text = re.sub(
        r'(?i)(preparation of .+?:)',
        r'\nPREPARATION_START\1',
        text
    )
    
    # Highlight analytical data
    text = re.sub(
        r'(?i)((?:^|\n)(?:1H |13C )?NMR[^:]*:.+)',
        r'\nANALYTICAL_DATA\1',
        text
    )
    
    return text


def _split_into_sections(text: str, max_length: int = 4000) -> List[str]:
    """
    Split text into processable sections.
    
    Args:
        text: Text to split
        max_length: Maximum section length
        
    Returns:
        List of text sections
    """
    sections = []
    
    # Try to split on natural boundaries
    boundaries = [
        r'\n\s*(?:example|preparation|compound)\s+\d+',
        r'\n\s*\d+\.\s+',
        r'\n\s*[A-Z][^a-z]+:',
        r'\n\s*\n',
        r'\.\s+'
    ]
    
    current_section = ""
    for line in text.split('\n'):
        if len(current_section) + len(line) > max_length:
            # Try each boundary pattern
            for pattern in boundaries:
                if re.search(pattern, current_section, re.I):
                    split_sections = re.split(
                        pattern,
                        current_section,
                        maxsplit=1
                    )
                    if len(split_sections) > 1:
                        sections.append(split_sections[0].strip())
                        current_section = split_sections[1].strip()
                        break
            else:
                # If no natural boundary found, split at max length
                sections.append(current_section.strip())
                current_section = ""
        current_section += line + "\n"
    
    if current_section.strip():
        sections.append(current_section.strip())
    
    return sections


def _validate_chemical_data(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate and standardize extracted chemical data.
    
    Args:
        data: Dictionary of extracted data
        
    Returns:
        Validated and standardized data
    """
    validated = {}
    
    # Validate CAS number
    if 'cas_number' in data:
        cas = str(data['cas_number'])
        if re.match(r'^\d{1,7}-\d{2}-\d$', cas):
            validated['cas_number'] = cas
            
    # Validate SMILES
    if 'smiles' in data:
        smiles = str(data['smiles'])
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            validated['smiles'] = Chem.MolToSmiles(mol, isomericSmiles=True)
            validated['inchi'] = Chem.MolToInchi(mol)
            validated['inchi_key'] = Chem.MolToInchiKey(mol)
            
    # Validate chemical names
    if 'names' in data and isinstance(data['names'], list):
        validated['names'] = [
            str(name) for name in data['names']
            if isinstance(name, (str, int, float))
        ]
        
    # Validate properties
    if 'properties' in data and isinstance(data['properties'], dict):
        validated['properties'] = {}
        for key, value in data['properties'].items():
            if isinstance(value, (int, float, str)):
                validated['properties'][str(key)] = value
                
    # Copy other fields
    for key in ['pharmacology', 'safety', 'regulatory', 'references']:
        if key in data:
            validated[key] = data[key]
            
    return validated


def _validate_patent_data(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate and standardize patent extracted data.
    
    Args:
        data: Dictionary of extracted data
        
    Returns:
        Validated and standardized data
    """
    validated = _validate_chemical_data(data)
    
    # Validate synthesis data if present
    if 'synthesis' in data and isinstance(data['synthesis'], dict):
        validated['synthesis'] = {
            'reagents': data['synthesis'].get('reagents', []),
            'conditions': data['synthesis'].get('conditions', {}),
            'yield': data['synthesis'].get('yield'),
            'procedure': data['synthesis'].get('procedure', '')
        }
        
    # Validate analytical data
    if 'analytical' in data and isinstance(data['analytical'], dict):
        validated['analytical'] = {
            'nmr': data['analytical'].get('nmr', {}),
            'ms': data['analytical'].get('ms', {}),
            'elemental': data['analytical'].get('elemental', {})
        }
        
    return validated


def _deduplicate_compounds(
    compounds: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """
    Deduplicate compounds based on structure.
    
    Args:
        compounds: List of compound dictionaries
        
    Returns:
        Deduplicated list of compounds
    """
    unique_compounds = {}
    
    for compound in compounds:
        # Try to get a unique identifier
        identifier = None
        
        # First try InChI Key
        if 'inchi_key' in compound:
            identifier = compound['inchi_key']
        # Then try CAS number
        elif 'cas_number' in compound:
            identifier = compound['cas_number']
        # Finally try SMILES
        elif 'smiles' in compound:
            mol = Chem.MolFromSmiles(compound['smiles'])
            if mol is not None:
                identifier = Chem.MolToInchiKey(mol)
                
        if identifier:
            if identifier not in unique_compounds:
                unique_compounds[identifier] = compound
            else:
                # Merge with existing compound
                existing = unique_compounds[identifier]
                for key, value in compound.items():
                    if key not in existing:
                        existing[key] = value
                    elif isinstance(value, list):
                        existing[key].extend(value)
                        existing[key] = list(set(existing[key]))
                    elif isinstance(value, dict):
                        existing[key].update(value)
                        
    return list(unique_compounds.values())
