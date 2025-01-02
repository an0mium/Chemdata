"""Data processing with batch operations and parallel execution."""

import concurrent.futures
from dataclasses import asdict
from typing import Dict, List, Optional, Set

import pandas as pd
from tqdm import tqdm

from api_client import McpError, PubChemClient
from cache_manager import CacheManager
from config import BATCH_SIZE, DATA_SOURCES, MAX_WORKERS, REQUIRED_FIELDS
from logger import LogManager
from models import CompoundData, ValidationError


class DataProcessor:
    """Handles batch processing and parallel execution of data collection."""
    
    def __init__(self):
        """Initialize data processor."""
        self.logger = LogManager().get_logger("data_processor")
        self.cache = CacheManager()
        self.pubchem = PubChemClient()
        
    def validate_compound(self, compound: CompoundData) -> List[str]:
        """
        Validate compound data.
        
        Args:
            compound: Compound data to validate
            
        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []
        
        # Check required fields
        for field in REQUIRED_FIELDS:
            value = getattr(compound, field)
            if value is None or value == "N/A" or value == "":
                errors.append(f"Missing required field: {field}")
                
        # Validate CAS number format if present
        if compound.cas != "N/A":
            if not self._validate_cas_format(compound.cas):
                errors.append(f"Invalid CAS number format: {compound.cas}")
                
        # Validate molecular weight is positive
        if compound.molecular_weight <= 0:
            errors.append(
                f"Invalid molecular weight: {compound.molecular_weight}"
            )
            
        return errors

    def _validate_cas_format(self, cas: str) -> bool:
        """
        Validate CAS number format.
        
        Args:
            cas: CAS number to validate
            
        Returns:
            True if valid, False otherwise
        """
        import re
        pattern = r'^\d{1,7}-\d{2}-\d$'
        return bool(re.match(pattern, cas))

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
                
            # Process each source in priority order
            for source in sorted(
                sources,
                key=lambda s: DATA_SOURCES[s]['priority']
            ):
                try:
                    if source == 'pubchem':
                        self._enrich_from_pubchem(compound)
                    # Add other data sources here
                    
                except McpError as e:
                    self.logger.warning(
                        f"Error enriching {compound.name} from {source}: {str(e)}"
                    )
                    
            return compound
            
        except Exception as e:
            self.logger.error(
                f"Error processing compound {compound.name}: {str(e)}"
            )
            return None

    def _enrich_from_pubchem(self, compound: CompoundData):
        """
        Enrich compound data from PubChem.
        
        Args:
            compound: Compound to enrich
        """
        # Try to get PubChem data
        if compound.pubchem_cid != "N/A":
            data = self.pubchem.get_compound_by_cid(compound.pubchem_cid)
        else:
            data = self.pubchem.get_compound_by_name(compound.name)
            
        if not data:
            return
            
        # Update compound with PubChem data
        compound.pubchem_cid = str(data.get('id', compound.pubchem_cid))
        compound.iupac_name = data.get('iupac_name', compound.iupac_name)
        compound.molecular_weight = float(
            data.get('molecular_weight', compound.molecular_weight)
        )
        compound.data_sources.add('PubChem')

    def process_file(
        self,
        input_file: str,
        output_file: str,
        sources: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """
        Process compounds from input file.
        
        Args:
            input_file: Input CSV file path
            output_file: Output CSV file path
            sources: Optional list of data sources to use
            
        Returns:
            Processing statistics
        """
        # Read input file
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
        result_df.to_csv(output_file, index=False)
        
        # Calculate statistics
        stats = {
            'total_compounds': len(compounds),
            'processed_compounds': len(processed),
            'success_rate': len(processed) / len(compounds) * 100,
            'sources_used': set().union(*(
                c.data_sources for c in processed
            )),
            'missing_required_fields': sum(
                1 for c in processed
                if any(
                    getattr(c, f) == "N/A"
                    for f in REQUIRED_FIELDS
                )
            )
        }
        
        self.logger.info(
            f"Processing completed:\n"
            f"- Total compounds: {stats['total_compounds']}\n"
            f"- Successfully processed: {stats['processed_compounds']}\n"
            f"- Success rate: {stats['success_rate']:.1f}%\n"
            f"- Sources used: {', '.join(stats['sources_used'])}\n"
            f"- Missing required fields: {stats['missing_required_fields']}"
        )
        
        return stats
