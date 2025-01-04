"""Binding data processing and BindingDB integration.

This module provides a comprehensive framework for:
1. Loading and processing binding data from BindingDB
2. Gathering receptor ligands from multiple sources
3. Enriching compound data with web sources
4. Patent searching and analysis
5. Structure validation and standardization
6. Swiss tools integration for target prediction and ADME properties
7. Machine learning predictions for binding and activity
8. Community data integration and analysis
"""

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
import logging
import pandas as pd
from tqdm import tqdm
from rich.logging import RichHandler

__version__ = "0.2.0"

# Import core components
from .models.compound import Compound, CompoundData, CompoundType, LegalStatus
from .models.psychopharm import PsychoactiveCompound
from .pipeline.base import PipelineManager
from .web_enrichment.manager import WebEnrichmentManager
from .processors.patent import PatentProcessor

# Import processors
from .processors.bindingdb import BindingDBProcessor
from .processors.activity import ActivityProcessor
from .processors.structure import StructureProcessor
from .processors.pubmed import PubMedProcessor
from .processors.swiss import SwissToolsProcessor

# Import utilities
from .utils.names import NameCleaner
from .utils.identifiers import IdentifierLookup
from .utils.checkpoints import CheckpointManager
from logger import LogManager

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)]
)

# Suppress noisy third-party loggers
logging.getLogger("urllib3").setLevel(logging.WARNING)
logging.getLogger("selenium").setLevel(logging.WARNING)
logging.getLogger("pytorch_lightning").setLevel(logging.WARNING)
logging.getLogger("rdkit").setLevel(logging.WARNING)

# Create logger for this package
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Log version on import
logger.info(f"Binding Data Processor v{__version__}")

# Import optional dependencies
try:
    import rdkit
    logger.debug("RDKit available for structure processing")
except ImportError:
    logger.warning("RDKit not available - structure processing will be limited")

try:
    import torch
    logger.debug("PyTorch available for ML models")
except ImportError:
    logger.warning("PyTorch not available - ML functionality will be limited")

try:
    import dash
    logger.debug("Dash available for web interface")
except ImportError:
    logger.warning("Dash not available - web interface will be disabled")

# Configure default settings
DEFAULT_CACHE_DIR = Path(os.getenv("CHEMDATA_CACHE_DIR", "~/.cache/chemdata")).expanduser()
DEFAULT_CONFIG_DIR = Path(os.getenv("CHEMDATA_CONFIG_DIR", "~/.config/chemdata")).expanduser()
DEFAULT_DATA_DIR = Path(os.getenv("CHEMDATA_DATA_DIR", "~/.local/share/chemdata")).expanduser()

# Create directories if they don't exist
for directory in [DEFAULT_CACHE_DIR, DEFAULT_CONFIG_DIR, DEFAULT_DATA_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# Export default settings
__default_settings__ = {
    "cache_dir": DEFAULT_CACHE_DIR,
    "config_dir": DEFAULT_CONFIG_DIR,
    "data_dir": DEFAULT_DATA_DIR,
    "log_level": "INFO",
    "n_workers": os.cpu_count() or 1,
    "batch_size": 100,
    "rate_limit": 2,
    "timeout": 30,
    "max_retries": 3,
}


class BindingDataProcessor:
    """Main class for handling binding data processing and integration."""

    def __init__(self, pubmed_client, web_client=None, http_client=None):
        """
        Initialize binding data processor.

        Args:
            pubmed_client: PubMed client for relevance scoring
            web_client: Optional WebEnrichment client for web scraping
            http_client: Optional HTTP client for Swiss tools
        """
        self.logger = LogManager().get_logger("binding_data_processor")

        # Initialize sub-processors
        self.bindingdb = BindingDBProcessor()
        self.activity = ActivityProcessor()
        self.structure = StructureProcessor()
        self.patent = PatentProcessor(http_client) if http_client else None
        self.pubmed = PubMedProcessor(pubmed_client)
        self.swiss = SwissToolsProcessor(http_client) if http_client else None

        # Initialize utilities
        self.name_cleaner = NameCleaner()
        self.id_lookup = IdentifierLookup(web_client) if web_client else None
        self.checkpoint_manager = CheckpointManager()

    def process_compound(self, name: str, smiles: Optional[str] = None) -> CompoundData:
        """Process a single compound to gather all available data."""
        try:
            # Clean compound name
            clean_name = self.name_cleaner.clean_name(name)

            # Get binding data from BindingDB
            binding_data = self.bindingdb.load_data(clean_name, smiles)

            # Create compound object
            compound = CompoundData(
                name=clean_name,
                smiles=smiles,
                compound_type=CompoundType.OTHER,
                binding_data=binding_data
            )

            # Validate structure if SMILES provided
            if smiles:
                if not self.structure.validate_structure(smiles):
                    self.logger.warning(f"Invalid structure for compound: {clean_name}")
                    return None

            # Enrich with additional data if available
            if self.swiss:
                self.swiss.enrich_compound(compound)
            if self.id_lookup:
                self.id_lookup.enrich_identifiers(compound)
            if self.pubmed:
                self.pubmed.enrich_references(compound)
            if self.patent:
                self.patent.enrich_patent_data(compound)

            return compound

        except Exception as e:
            self.logger.error(f"Error processing compound {name}: {str(e)}")
            return None

    def gather_ligands(
        self,
        target_pattern: Optional[str] = None,
        llm_api_key: Optional[str] = None,
        use_checkpoints: bool = True,
    ) -> List[CompoundData]:
        """Gather receptor ligands from multiple sources."""
        try:
            compounds = []

            # Clear checkpoints if using them
            if use_checkpoints:
                self.checkpoint_manager.clear_checkpoints()

            # Get compounds from BindingDB
            self.logger.info("Gathering compounds from BindingDB...")
            bindingdb_compounds = self.bindingdb.gather_ligands(
                target_pattern, use_checkpoints=use_checkpoints
            )
            compounds.extend(bindingdb_compounds)

            # Search patents if API key provided
            if llm_api_key and self.patent:
                self.logger.info("Searching patents for compounds...")
                patent_compounds = self.patent.search_compounds(
                    target_pattern, llm_api_key, use_checkpoints=use_checkpoints
                )
                compounds.extend(patent_compounds)

            # Search PubMed
            if self.pubmed:
                self.logger.info("Searching PubMed for compounds...")
                pubmed_compounds = self.pubmed.search_compounds(
                    target_pattern, use_checkpoints=use_checkpoints
                )
                compounds.extend(pubmed_compounds)

            # Deduplicate compounds
            self.logger.info("Deduplicating compounds...")
            unique_compounds = self._deduplicate_compounds(compounds)

            # Enrich with Swiss tools data
            if self.swiss:
                self.logger.info("Enriching compounds with Swiss tools data...")
                for compound in tqdm(unique_compounds, desc="Swiss tools enrichment"):
                    try:
                        self.swiss.enrich_compound(compound)
                    except Exception as e:
                        self.logger.error(
                            f"Error enriching compound {compound.name}: {str(e)}"
                        )

            return unique_compounds

        except Exception as e:
            self.logger.error(f"Error gathering ligands: {str(e)}")
            return []

    def _deduplicate_compounds(
        self, compounds: List[CompoundData]
    ) -> List[CompoundData]:
        """Deduplicate compounds by InChI Key."""
        unique_compounds = []
        seen_keys = set()

        for compound in compounds:
            if not compound.inchi_key:
                continue

            if compound.inchi_key not in seen_keys:
                seen_keys.add(compound.inchi_key)
                unique_compounds.append(compound)

        return unique_compounds

    def save_compounds(
        self,
        compounds: List[CompoundData],
        output_path: str,
        include_swiss_data: bool = True,
    ) -> None:
        """Save compounds to TSV file."""
        try:
            # Ensure output directory exists
            os.makedirs(os.path.dirname(output_path), exist_ok=True)

            # Convert compounds to rows with progress bar
            self.logger.info("Converting compounds to rows...")
            rows = []
            for compound in tqdm(compounds, desc="Processing compounds"):
                row = compound.to_dict(include_predictions=include_swiss_data)
                rows.append(row)

            # Write TSV file
            df = pd.DataFrame(rows)
            df.to_csv(output_path, sep="\t", index=False)

            # Log statistics
            self.logger.info(f"Saved {len(compounds)} compounds to {output_path}")

            activity_types = {}
            sources = set()
            with_cas = 0
            with_patents = 0

            for compound in compounds:
                if compound.primary_activity:
                    activity_types[compound.primary_activity] = (
                        activity_types.get(compound.primary_activity, 0) + 1
                    )
                if hasattr(compound, 'data_sources'):
                    sources.update(compound.data_sources)
                if compound.cas_number:
                    with_cas += 1
                if compound.patent_data:
                    with_patents += 1

            self.logger.info(
                f"Data sources used: {', '.join(sources)}\n"
                f"Compounds with CAS numbers: {with_cas}\n"
                f"Compounds with patent data: {with_patents}\n"
                f"Activity type distribution:\n"
                + "\n".join(f" {k}: {v}" for k, v in activity_types.items())
            )

        except Exception as e:
            self.logger.error(f"Error saving compounds: {str(e)}")


# Expose key classes
__all__ = [
    "BindingDataProcessor",
    "Compound",
    "CompoundData",
    "CompoundType",
    "LegalStatus",
    "PsychoactiveCompound",
    "PipelineManager",
    "WebEnrichmentManager",
    "PatentProcessor",
    "ActivityProcessor",
    "StructureProcessor",
    "PubMedProcessor",
    "SwissToolsProcessor",
]
