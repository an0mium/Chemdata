"""Utilities for handling chemical names and identifiers.

This module handles:
1. Chemical name cleaning
2. Identifier extraction
3. Name normalization
4. IUPAC name standardization
5. Common name generation
6. Systematic name parsing
"""

import re
from typing import Tuple, Optional, List, Dict, Set, Any


class ChemicalNameNormalizer:
    """Handles chemical name normalization and standardization."""
    
    # Common chemical group patterns
    MULTIPLIERS = {
        'di': 2, 'tri': 3, 'tetra': 4, 'penta': 5,
        'hexa': 6, 'hepta': 7, 'octa': 8, 'nona': 9, 'deca': 10
    }
    
    POSITION_PREFIXES = {
        'ortho': '1,2-', 'meta': '1,3-', 'para': '1,4-',
        'o-': '1,2-', 'm-': '1,3-', 'p-': '1,4-'
    }
    
    GREEK_LETTERS = {
        'alpha': 'α', 'beta': 'β', 'gamma': 'γ', 'delta': 'δ',
        'epsilon': 'ε', 'zeta': 'ζ', 'eta': 'η', 'theta': 'θ',
        'iota': 'ι', 'kappa': 'κ', 'lambda': 'λ', 'mu': 'μ',
        'nu': 'ν', 'xi': 'ξ', 'omicron': 'ο', 'pi': 'π',
        'rho': 'ρ', 'sigma': 'σ', 'tau': 'τ', 'upsilon': 'υ',
        'phi': 'φ', 'chi': 'χ', 'psi': 'ψ', 'omega': 'ω'
    }
    
    SALT_FORMS = [
        'hydrochloride', 'HCl', 'hydrobromide', 'HBr',
        'hydroiodide', 'HI', 'sulfate', 'sulphate',
        'phosphate', 'acetate', 'citrate', 'maleate',
        'tartrate', 'fumarate', 'oxalate', 'besylate',
        'mesylate', 'tosylate', 'salt', 'hydrate',
        'solvate', 'dihydrate', 'trihydrate', 'tetrahydrate',
        'pentahydrate', 'hexahydrate', 'heptahydrate',
        'octahydrate', 'nonahydrate', 'decahydrate'
    ]
    
    COMMON_SUBSTITUENTS = {
        'methyl': 'Me', 'ethyl': 'Et', 'propyl': 'Pr',
        'butyl': 'Bu', 'phenyl': 'Ph', 'benzyl': 'Bn',
        'acetyl': 'Ac', 'hydroxyl': 'OH', 'amino': 'NH2',
        'carboxyl': 'COOH', 'methoxy': 'OMe', 'ethoxy': 'OEt'
    }
    
    def __init__(self):
        """Initialize name normalizer with compiled patterns."""
        # Compile common regex patterns
        self.cas_pattern = re.compile(r'\d{1,7}-\d{2}-\d')
        self.chembl_pattern = re.compile(r'CHEMBL\d+')
        self.patent_pattern = re.compile(r'[A-Z]{2}\d+[A-Z]\d?')
        self.position_pattern = re.compile(r'(\d+)[,-](\d+)')
        self.stereochem_pattern = re.compile(r'\([^)]*\)|[RS]-|\b[RS]\b')
        
        # Build salt form pattern
        salt_pattern = '|'.join(map(re.escape, self.SALT_FORMS))
        self.salt_pattern = re.compile(rf'\s*(?:{salt_pattern})$', re.I)
        
        # Build substituent pattern
        subst_pattern = '|'.join(map(re.escape, self.COMMON_SUBSTITUENTS.keys()))
        self.subst_pattern = re.compile(rf'\b({subst_pattern})\b', re.I)


def extract_identifiers(name: str) -> Tuple[str, Optional[str], Optional[str], Optional[str]]:
    """
    Extract identifiers from compound name.
    
    Args:
        name: Raw compound name
        
    Returns:
        Tuple of (clean_name, chembl_id, cas_number, patent_id)
    """
    normalizer = ChemicalNameNormalizer()
    
    # Extract ChEMBL ID
    chembl_match = re.search(r'::(?:CHEMBL\d+)', name)
    chembl_id = chembl_match.group(0)[2:] if chembl_match else None
    
    # Extract patent info
    patent_match = re.search(r'::(US\d+)', name)
    patent_id = patent_match.group(1) if patent_match else None
    
    # Try to extract CAS number
    cas_match = normalizer.cas_pattern.search(name)
    cas_number = cas_match.group(0) if cas_match else None
    
    # Clean name
    clean_name = name
    
    # Remove patent and example numbers
    clean_name = re.sub(r'::[A-Z0-9]+(?:,?\s*(?:Example|Compound)\s*[A-Za-z0-9-]+)*$', '', clean_name)
    clean_name = re.sub(r'US\d+\s*(?:Example|Compound)\s*[A-Za-z0-9-]+', '', clean_name)
    
    # Remove stereochemistry markers
    clean_name = normalizer.stereochem_pattern.sub('', clean_name)
    
    # Remove salt forms
    clean_name = normalizer.salt_pattern.sub('', clean_name)
    
    # Remove common prefixes/suffixes
    clean_name = re.sub(r'^(?:\+\/\-|\+|\-)\s*', '', clean_name)
    
    # Remove special characters but keep some important ones
    clean_name = re.sub(r'[^\w\s\-\+]', '', clean_name)
    
    # Normalize whitespace
    clean_name = ' '.join(clean_name.split())
    
    return clean_name, chembl_id, cas_number, patent_id


def clean_name(name: str) -> str:
    """
    Clean chemical name for comparison.
    
    Args:
        name: Chemical name to clean
        
    Returns:
        Cleaned and normalized name
    """
    normalizer = ChemicalNameNormalizer()
    
    # Remove stereochemistry
    name = normalizer.stereochem_pattern.sub('', name)
    
    # Remove salt forms
    name = normalizer.salt_pattern.sub('', name)
    
    # Remove special characters but keep some important ones
    name = re.sub(r'[^\w\s\-\+]', '', name)
    
    # Normalize whitespace
    return ' '.join(name.split()).lower()


def standardize_chemical_name(name: str) -> str:
    """
    Standardize chemical name to IUPAC-like format.
    
    Args:
        name: Chemical name to standardize
        
    Returns:
        Standardized name
    """
    normalizer = ChemicalNameNormalizer()
    
    # First clean the name
    name = clean_name(name)
    
    # Replace position prefixes
    for prefix, replacement in normalizer.POSITION_PREFIXES.items():
        name = re.sub(rf'\b{prefix}\b', replacement, name, flags=re.I)
    
    # Replace Greek letters
    for letter, symbol in normalizer.GREEK_LETTERS.items():
        name = re.sub(rf'\b{letter}\b', symbol, name, flags=re.I)
    
    # Handle common chemical group formatting
    for multiplier in normalizer.MULTIPLIERS.keys():
        name = re.sub(rf'\b{multiplier}(\w)', rf'{multiplier}-\1', name, flags=re.I)
    
    # Handle position numbers
    name = normalizer.position_pattern.sub(r'\1,\2-', name)
    
    # Handle substituent abbreviations
    for full, abbrev in normalizer.COMMON_SUBSTITUENTS.items():
        name = re.sub(rf'\b{full}\b', abbrev, name, flags=re.I)
    
    # Capitalize first letter if it's not a number
    if name and name[0].isalpha():
        name = name[0].upper() + name[1:]
    
    return name


def get_name_variants(name: str, include_abbreviations: bool = True) -> List[str]:
    """
    Generate variants of a chemical name for searching.
    
    Args:
        name: Base chemical name
        include_abbreviations: Whether to include abbreviated forms
        
    Returns:
        List of name variants
    """
    normalizer = ChemicalNameNormalizer()
    variants = set()
    
    # Add original name
    variants.add(name)
    
    # Add cleaned name
    clean = clean_name(name)
    variants.add(clean)
    
    # Add standardized name
    std = standardize_chemical_name(name)
    variants.add(std)
    
    # Replace Greek letters
    for letter, symbol in normalizer.GREEK_LETTERS.items():
        if letter in name.lower():
            variants.add(name.lower().replace(letter, symbol))
    
    # Handle position prefixes
    for prefix, replacement in normalizer.POSITION_PREFIXES.items():
        if prefix in name.lower():
            variants.add(name.lower().replace(prefix, replacement))
    
    # Handle position numbers
    numbered = normalizer.position_pattern.sub(r'\1,\2-', name)
    if numbered != name:
        variants.add(numbered)
    
    # Add abbreviated forms if requested
    if include_abbreviations:
        name_lower = name.lower()
        for full, abbrev in normalizer.COMMON_SUBSTITUENTS.items():
            if full in name_lower:
                variants.add(name_lower.replace(full, abbrev))
    
    # Remove variants that are None or empty
    variants = {v for v in variants if v and len(v.strip()) > 0}
    
    return list(variants)


def validate_cas_number(cas: str) -> bool:
    """
    Validate CAS Registry Number format and checksum.
    
    Args:
        cas: CAS number to validate
        
    Returns:
        True if valid, False otherwise
    """
    if not cas or not isinstance(cas, str):
        return False
    
    # Check format
    if not re.match(r'^\d{1,7}-\d{2}-\d$', cas):
        return False
    
    # Remove hyphens
    digits = cas.replace('-', '')
    
    # Calculate checksum
    total = 0
    for i in range(len(digits) - 1):
        total += int(digits[i]) * (len(digits) - 1 - i)
    
    check_digit = total % 10
    
    return check_digit == int(digits[-1])


def parse_systematic_name(name: str) -> Dict[str, Any]:
    """
    Parse systematic chemical name into components.
    
    Args:
        name: Systematic chemical name
        
    Returns:
        Dictionary containing parsed components
    """
    normalizer = ChemicalNameNormalizer()
    components = {
        'base': None,
        'substituents': [],
        'positions': [],
        'multipliers': [],
        'stereochemistry': [],
        'locants': []
    }
    
    # First clean and standardize
    name = standardize_chemical_name(name)
    
    # Extract stereochemistry
    stereo_matches = normalizer.stereochem_pattern.finditer(name)
    for match in stereo_matches:
        components['stereochemistry'].append(match.group())
    name = normalizer.stereochem_pattern.sub('', name)
    
    # Extract position numbers
    pos_matches = normalizer.position_pattern.finditer(name)
    for match in pos_matches:
        components['positions'].append(match.group())
    
    # Extract multipliers
    for multiplier in normalizer.MULTIPLIERS:
        if re.search(rf'\b{multiplier}\b', name, re.I):
            components['multipliers'].append(multiplier)
    
    # Extract substituents
    for subst in normalizer.COMMON_SUBSTITUENTS:
        if re.search(rf'\b{subst}\b', name, re.I):
            components['substituents'].append(subst)
    
    # What's left might be the base
    # Remove extracted parts and clean up
    for key in ['stereochemistry', 'positions', 'multipliers', 'substituents']:
        for item in components[key]:
            name = re.sub(rf'\b{re.escape(item)}\b', '', name)
    
    # Clean up and set base
    base = ' '.join(name.split())
    if base:
        components['base'] = base
    
    return components
