"""Web enrichment pipeline.

This module provides the WebPipeline class that:
1. Coordinates web data enrichment
2. Manages multiple data sources
3. Handles data validation
4. Tracks enrichment statistics
5. Integrates with ML predictions
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from datetime import datetime

from ...models.compound.enhanced import EnhancedCompound
from ...web_enrichment.manager import WebEnrichmentManager, EnrichmentConfig


@dataclass
class WebPipelineConfig:
    """Web pipeline configuration."""
    
    # Enrichment configuration
    enrichment_config: EnrichmentConfig = field(default_factory=EnrichmentConfig)
    
    # Pipeline settings
    batch_size: int = 32
    n_workers: int = 4
    use_cache: bool = True
    
    # Validation settings
    validate_data: bool = True
    min_sources: int = 2
    min_confidence: float = 0.5
    
    # Integration settings
    update_predictions: bool = True
    track_novel_compounds: bool = True


@dataclass
class EnrichmentStats:
    """Web enrichment statistics."""
    
    # Processing stats
    total_compounds: int = 0
    processed_compounds: int = 0
    failed_compounds: int = 0
    
    # Source stats
    swiss_compounds: int = 0
    community_compounds: int = 0
    social_compounds: int = 0
    
    # Data stats
    total_references: int = 0
    novel_compounds: int = 0
    validation_errors: int = 0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "processing": {
                "total": self.total_compounds,
                "processed": self.processed_compounds,
                "failed": self.failed_compounds,
                "success_rate": self._get_success_rate(),
            },
            "sources": {
                "swiss": self.swiss_compounds,
                "community": self.community_compounds,
                "social": self.social_compounds,
            },
            "data": {
                "references": self.total_references,
                "novel_compounds": self.novel_compounds,
                "validation_errors": self.validation_errors,
            },
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get processing success rate."""
        if not self.total_compounds:
            return None
        return self.processed_compounds / self.total_compounds


class WebPipeline:
    """Pipeline for web enrichment of compounds."""

    def __init__(
        self,
        config: Optional[WebPipelineConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize web pipeline.
        
        Args:
            config: Optional pipeline configuration
            logger: Optional logger instance
        """
        self.config = config or WebPipelineConfig()
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        
        # Initialize enrichment manager
        self.enrichment = WebEnrichmentManager(
            config=self.config.enrichment_config,
            logger=self.logger,
        )
        
        # Initialize stats
        self.stats = EnrichmentStats()
        
        # Track novel compounds
        self.novel_compounds: Dict[str, Dict] = {}

    def process_compounds(
        self,
        compounds: List[EnhancedCompound],
        batch_size: Optional[int] = None,
    ) -> List[EnhancedCompound]:
        """Process multiple compounds.
        
        Args:
            compounds: List of compounds to process
            batch_size: Optional batch size override
            
        Returns:
            List of processed compounds
        """
        self.stats.total_compounds += len(compounds)
        batch_size = batch_size or self.config.batch_size
        
        # Process in batches
        for i in range(0, len(compounds), batch_size):
            batch = compounds[i:i + batch_size]
            try:
                self._process_batch(batch)
                self.stats.processed_compounds += len(batch)
            except Exception as e:
                self.logger.error(f"Failed to process batch: {str(e)}")
                self.stats.failed_compounds += len(batch)
        
        return compounds

    def process_compound(
        self,
        compound: EnhancedCompound,
    ) -> EnhancedCompound:
        """Process a single compound.
        
        Args:
            compound: Compound to process
            
        Returns:
            Processed compound
        """
        self.stats.total_compounds += 1
        
        try:
            # Enrich with web data
            self.enrichment.enrich_compounds(
                [compound],
                use_cache=self.config.use_cache,
            )
            
            # Update stats
            self._update_source_stats(compound)
            
            # Validate enriched data
            if self.config.validate_data:
                self._validate_enriched_data(compound)
            
            # Update predictions if enabled
            if self.config.update_predictions:
                self._update_predictions(compound)
            
            # Track novel compounds if enabled
            if self.config.track_novel_compounds:
                self._track_novel_compounds(compound)
            
            self.stats.processed_compounds += 1
            return compound
            
        except Exception as e:
            self.logger.error(
                f"Failed to process compound {compound.name}: {str(e)}"
            )
            self.stats.failed_compounds += 1
            raise

    def _process_batch(
        self,
        compounds: List[EnhancedCompound],
    ) -> None:
        """Process a batch of compounds."""
        # Enrich compounds
        self.enrichment.enrich_compounds(
            compounds,
            use_cache=self.config.use_cache,
        )
        
        # Process enriched data
        for compound in compounds:
            try:
                # Update stats
                self._update_source_stats(compound)
                
                # Validate data
                if self.config.validate_data:
                    self._validate_enriched_data(compound)
                
                # Update predictions
                if self.config.update_predictions:
                    self._update_predictions(compound)
                
                # Track novel compounds
                if self.config.track_novel_compounds:
                    self._track_novel_compounds(compound)
                    
            except Exception as e:
                self.logger.error(
                    f"Error processing enriched data for {compound.name}: {str(e)}"
                )

    def _update_source_stats(
        self,
        compound: EnhancedCompound,
    ) -> None:
        """Update source statistics."""
        if compound.swiss_data:
            self.stats.swiss_compounds += 1
        if compound.community_data:
            self.stats.community_compounds += 1
        if compound.social_data:
            self.stats.social_compounds += 1
            
        # Update reference count
        self.stats.total_references += len(compound.reference_dois)
        self.stats.total_references += len(compound.reference_pmids)
        self.stats.total_references += len(compound.reference_urls)

    def _validate_enriched_data(
        self,
        compound: EnhancedCompound,
    ) -> None:
        """Validate enriched compound data."""
        try:
            # Check source count
            sources = []
            if compound.swiss_data:
                sources.append("swiss")
            if compound.community_data:
                sources.append("community")
            if compound.social_data:
                sources.append("social")
                
            if len(sources) < self.config.min_sources:
                self.logger.warning(
                    f"Limited sources for {compound.name}: {len(sources)}"
                )
                self.stats.validation_errors += 1
            
            # Check confidence scores
            for source in sources:
                data = getattr(compound, f"{source}_data")
                if not data:
                    continue
                    
                for key, value in data.items():
                    if isinstance(value, dict) and "confidence" in value:
                        if value["confidence"] < self.config.min_confidence:
                            self.logger.warning(
                                f"Low confidence {source} data for {compound.name}: "
                                f"{key} = {value['confidence']:.3f}"
                            )
                            self.stats.validation_errors += 1
            
        except Exception as e:
            self.logger.error(
                f"Error validating data for {compound.name}: {str(e)}"
            )
            self.stats.validation_errors += 1

    def _update_predictions(
        self,
        compound: EnhancedCompound,
    ) -> None:
        """Update predictions with web data."""
        predictions = compound.get_predictions()
        if not predictions:
            return
            
        # Update each prediction
        for pred_type, pred in predictions.items():
            # Get supporting web data
            support = self._get_web_support(compound, pred_type)
            if support:
                # Update confidence
                pred["confidence"] = max(
                    pred["confidence"],
                    support["confidence"]
                )
                # Add supporting data
                pred["web_support"] = support

    def _get_web_support(
        self,
        compound: EnhancedCompound,
        pred_type: str,
    ) -> Optional[Dict]:
        """Get supporting web data for prediction."""
        support = {
            "sources": [],
            "confidence": 0.0,
            "evidence": [],
        }
        
        # Check each data source
        for source in ["swiss", "community", "social"]:
            data = getattr(compound, f"{source}_data")
            if not data or pred_type not in data:
                continue
                
            source_data = data[pred_type]
            if not isinstance(source_data, dict):
                continue
                
            # Add source
            support["sources"].append(source)
            
            # Update confidence
            if "confidence" in source_data:
                support["confidence"] = max(
                    support["confidence"],
                    source_data["confidence"]
                )
            
            # Add evidence
            if "evidence" in source_data:
                support["evidence"].extend(source_data["evidence"])
        
        return support if support["sources"] else None

    def _track_novel_compounds(
        self,
        compound: EnhancedCompound,
    ) -> None:
        """Track novel compounds from social data."""
        if not compound.social_data:
            return
            
        # Check for novel compounds
        if "novel_mentions" in compound.social_data:
            for mention in compound.social_data["novel_mentions"]:
                if mention not in self.novel_compounds:
                    self.novel_compounds[mention] = {
                        "first_seen": datetime.now().isoformat(),
                        "sources": [],
                        "references": [],
                    }
                
                # Update tracking data
                tracking = self.novel_compounds[mention]
                if "sources" in compound.social_data:
                    tracking["sources"].extend(
                        source for source in compound.social_data["sources"]
                        if source not in tracking["sources"]
                    )
                if "references" in compound.social_data:
                    tracking["references"].extend(
                        ref for ref in compound.social_data["references"]
                        if ref not in tracking["references"]
                    )
            
            self.stats.novel_compounds = len(self.novel_compounds)

    def get_novel_compounds(self) -> Dict[str, Dict]:
        """Get tracked novel compounds."""
        return self.novel_compounds

    def close(self) -> None:
        """Close pipeline and cleanup."""
        self.enrichment.close()
