"""Configuration settings for chemical data collection."""

import os
from pathlib import Path

# API Configuration
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
BINDINGDB_BASE_URL = "https://bindingdb.org/bind/downloads"

# Rate Limiting
API_RATE_LIMIT = 0.2  # seconds between requests
MAX_RETRIES = 3
RETRY_DELAY = 1.0  # seconds

# Batch Processing
BATCH_SIZE = 50
MAX_WORKERS = 4  # for parallel processing

# Cache Configuration
CACHE_DIR = Path(os.path.expanduser("~")) / ".chemical_data_collector" / "cache"
CACHE_EXPIRY = 24 * 60 * 60  # 24 hours in seconds

# Data Export
OUTPUT_FORMATS = ["csv", "json", "excel"]
DEFAULT_OUTPUT_FORMAT = "csv"

# Logging
LOG_DIR = Path(os.path.expanduser("~")) / ".chemical_data_collector" / "logs"
LOG_LEVEL = "INFO"
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

# Data Sources
DATA_SOURCES = {
    'bindingdb': {
        'name': 'BindingDB',
        'description': 'Binding data and chemical properties',
        'enabled': True,
        'priority': 1
    },
    'pubchem': {
        'name': 'PubChem',
        'description': 'Additional chemical data and pharmacology',
        'enabled': True,
        'priority': 2
    },
    'wikipedia': {
        'name': 'Wikipedia',
        'description': 'General information and references',
        'enabled': True,
        'priority': 3
    },
    'psychonaut': {
        'name': 'PsychonautWiki',
        'description': 'User reports and effects',
        'enabled': True,
        'priority': 4
    },
    'erowid': {
        'name': 'Erowid',
        'description': 'Experience reports and safety data',
        'enabled': True,
        'priority': 5
    }
}

# Validation Settings
REQUIRED_FIELDS = [
    'cas',
    'name',
    'smiles',
    'inchi_key',
    'molecular_weight'
]

# Circuit Breaker Settings
CIRCUIT_BREAKER = {
    'failure_threshold': 5,
    'reset_timeout': 60,  # seconds
    'half_open_timeout': 30  # seconds
}
