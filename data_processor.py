"""Data processing with batch operations and parallel execution."""

import concurrent.futures
from pathlib import Path
from dataclasses import asdict
from typing import Any, Dict, List, Optional
import pandas as pd
from tqdm import tqdm

from binding_data_processor import BindingDataProcessor
from chemical_properties import ChemicalProperties
from web_enrichment import WebEnrichment
from api_client import PubChemClient, PubMedClient
from cache_manager import CacheManager
from config import BATCH_SIZE, DATA_SOURCES, IDENTIFIER_FIELDS, MAX_WORKERS
from logger import LogManager
from models import CompoundData, ValidationError


class DataProcessor:
    """Handles batch processing and parallel execution of data collection."""

    def __init__(self):
        """Initialize data processor."""
        self.logger = LogManager().get_logger("data_processor")
        self.cache = CacheManager()
        self.pubchem = PubChemClient()
        self.pubmed = PubMedClient()
        
        # Initialize specialized processors
        self.binding_processor = BindingDataProcessor(self.pubmed)
        self.chemical_properties = ChemicalProperties()
        self.web_enrichment = WebEnrichment()

    def validate_compound(self, compound: CompoundData) -> List[str]:
        """
        Validate compound data.
        
        Args:
            compound: Compound data to validate
            
        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []

        # Check that at least one identifier is present
        has_identifier = False
        for field in IDENTIFIER_FIELDS:
            value = getattr(compound, field)
            if value is not None and value != "N/A" and value != "":
                has_identifier = True
                break

        if not has_identifier:
            errors.append("At least one identifier (CAS, name, or SMILES) is required")

        # Validate structure if SMILES present
        if compound.smiles != "N/A":
            is_valid, error = self.chemical_properties.validate_structure(compound.smiles)
            if not is_valid:
                errors.append(f"Invalid structure: {error}")

        return errors

    def process_batch(
        self,
        compounds: List[CompoundData],
        sources: Optional[List[str]] = None
    ) -> List[CompoundData]:
        """
        Process a batch of compounds.
        
        Args:
            compounds: List of compounds to process
            sources: Optional list of data sources to use
            
        Returns:
            List of processed compounds
        """
        if sources is None:
            sources = [
                source for source, config in DATA_SOURCES.items()
                if config['enabled']
            ]

        processed = []
        errors = []

        # Process compounds in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            future_to_compound = {
                executor.submit(
                    self._process_single_compound, compound, sources
                ): compound
                for compound in compounds
            }

            for future in concurrent.futures.as_completed(future_to_compound):
                compound = future_to_compound[future]
                try:
                    result = future.result()
                    if result:
                        processed.append(result)
                except Exception as e:
                    self.logger.error(
                        f"Error processing compound {compound.name}: {str(e)}"
                    )
                    errors.append((compound, str(e)))

        # Log errors summary
        if errors:
            self.logger.warning(
                f"Failed to process {len(errors)} compounds:"
                f"\n" + "\n".join(
                    f"- {c.name}: {e}" for c, e in errors
                )
            )

        return processed

    def _process_single_compound(
        self,
        compound: CompoundData,
        sources: List[str]
    ) -> Optional[CompoundData]:
        """
        Process a single compound using specified data sources.
        
        Args:
            compound: Compound to process
            sources: List of data sources to use
            
        Returns:
            Processed compound or None if processing failed
        """
        try:
            # Validate input
            errors = self.validate_compound(compound)
            if errors:
                raise ValidationError(
                    f"Validation failed for {compound.name}: {', '.join(errors)}"
                )

            # Calculate chemical properties if SMILES available
            if compound.smiles != "N/A":
                props = self.chemical_properties.calculate_properties(compound.smiles)
                if props:
                    for key, value in props.items():
                        setattr(compound, key, value)

            # Rank common names by search results
            common_names = self.web_enrichment.get_common_names(compound.name)
            for i, name_data in enumerate(common_names, 1):
                setattr(compound, f'common_name_{i}', name_data['name'])
                setattr(compound, f'common_name_{i}_source', name_data['source'])
                setattr(compound, f'common_name_{i}_relevance', name_data['relevance'])

            # Sort and enrich binding data
            self.binding_processor.sort_binding_data(compound)

            # Get legal status
            legal_status = self.web_enrichment.get_legal_status(compound.name)
            if legal_status:
                compound.scheduling.update(
                    {s['jurisdiction']: s['schedule'] for s in legal_status['scheduling']}
                )
                compound.data_sources.update(legal_status['sources'])

            # Get pharmacology
            pharm_info = self.web_enrichment.get_pharmacology(compound.name)
            if pharm_info:
                if pharm_info['mechanism_of_action']:
                    compound.mechanism_of_action = '; '.join(pharm_info['mechanism_of_action'])
                if pharm_info['toxicity']:
                    compound.toxicity = '; '.join(pharm_info['toxicity'])
                compound.data_sources.update(pharm_info['sources'])

            # Get reference URLs
            urls = self.web_enrichment.get_reference_urls(compound.name)
            for key, value in urls.items():
                setattr(compound, key, value)

            return compound

        except Exception as e:
            self.logger.error(
                f"Error processing compound {compound.name}: {str(e)}"
            )
            return None

    def process_file(
        self,
        input_file: str,
        output_file: str,
        sources: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """
        Process compounds from input file.
        
        Args:
            input_file: Input TSV/CSV file path
            output_file: Output TSV/CSV file path
            sources: Optional list of data sources to use
            
        Returns:
            Processing statistics
        """
        # Determine file format from extension
        input_ext = Path(input_file).suffix.lower()
        output_ext = Path(output_file).suffix.lower()
        
        # Read input file
        if input_ext == '.tsv':
            df = pd.read_csv(input_file, sep='\t')
        else:
            df = pd.read_csv(input_file)
            
        compounds = [
            CompoundData(**row)
            for _, row in df.iterrows()
        ]

        # Process in batches
        processed = []
        for i in tqdm(
            range(0, len(compounds), BATCH_SIZE),
            desc="Processing compounds"
        ):
            batch = compounds[i:i + BATCH_SIZE]
            processed.extend(self.process_batch(batch, sources))

        # Save results
        result_df = pd.DataFrame([
            asdict(compound) for compound in processed
        ])
        
        # Save in appropriate format
        if output_ext == '.tsv':
            result_df.to_csv(output_file, sep='\t', index=False)
        else:
            result_df.to_csv(output_file, index=False)

        # Calculate statistics
        stats = {
            'total_compounds': len(compounds),
            'processed_compounds': len(processed),
            'success_rate': len(processed) / len(compounds) * 100,
            'sources_used': set().union(*(
                c.data_sources for c in processed
            )),
            'missing_identifiers': sum(
                1 for c in processed
                if any(
                    getattr(c, f) == "N/A"
                    for f in IDENTIFIER_FIELDS
                )
            )
        }

        self.logger.info(
            f"Processing completed:\n"
            f"- Total compounds: {stats['total_compounds']}\n"
            f"- Successfully processed: {stats['processed_compounds']}\n"
            f"- Success rate: {stats['success_rate']:.1f}%\n"
            f"- Sources used: {', '.join(stats['sources_used'])}\n"
            f"- Missing identifiers: {stats['missing_identifiers']}"
        )

        return stats
