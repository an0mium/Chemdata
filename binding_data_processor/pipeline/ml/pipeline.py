"""ML pipeline for compound processing.

This module provides the MLPipeline class that:
1. Coordinates ML predictions
2. Manages model ensembles
3. Handles prediction caching
4. Tracks prediction statistics
5. Integrates with web enrichment
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field

from ...models.compound.enhanced import EnhancedCompound
from .manager import MLManager, ModelConfig, PredictionStats


@dataclass
class MLPipelineConfig:
    """ML pipeline configuration."""
    
    # Model configuration
    model_config: ModelConfig = field(default_factory=ModelConfig)
    
    # Pipeline settings
    batch_size: int = 32
    n_workers: int = 4
    use_cache: bool = True
    
    # Prediction thresholds
    min_confidence: float = 0.5
    min_support: int = 3
    
    # Web enrichment integration
    enrich_predictions: bool = True
    validate_web_data: bool = True


@dataclass
class PipelineStats:
    """Pipeline statistics."""
    
    # Processing stats
    total_compounds: int = 0
    processed_compounds: int = 0
    failed_compounds: int = 0
    
    # Prediction stats
    prediction_stats: PredictionStats = field(default_factory=PredictionStats)
    
    # Web enrichment stats
    enriched_compounds: int = 0
    enrichment_errors: int = 0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "processing": {
                "total": self.total_compounds,
                "processed": self.processed_compounds,
                "failed": self.failed_compounds,
                "success_rate": self._get_success_rate(),
            },
            "predictions": self.prediction_stats.to_dict(),
            "enrichment": {
                "enriched": self.enriched_compounds,
                "errors": self.enrichment_errors,
            },
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get processing success rate."""
        if not self.total_compounds:
            return None
        return self.processed_compounds / self.total_compounds


class MLPipeline:
    """Pipeline for ML processing of compounds."""

    def __init__(
        self,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        config: Optional[MLPipelineConfig] = None,
    ):
        """Initialize ML pipeline.
        
        Args:
            model_dir: Optional directory containing ML models
            cache_dir: Optional directory for caching
            config: Optional pipeline configuration
        """
        self.logger = logging.getLogger(self.__class__.__name__)
        self.model_dir = Path(model_dir) if model_dir else Path("models")
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.config = config or MLPipelineConfig()
        
        # Initialize ML manager
        self.ml_manager = MLManager(
            model_dir=self.model_dir,
            cache_dir=self.cache_dir,
            config=self.config.model_config,
        )
        
        # Initialize stats
        self.stats = PipelineStats()

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
            # Generate predictions
            compound = self.ml_manager.predict_compound(
                compound,
                use_cache=self.config.use_cache,
            )
            self.stats.prediction_stats = self.ml_manager.stats
            
            # Enrich with web data if enabled
            if self.config.enrich_predictions:
                compound = self._enrich_predictions(compound)
            
            # Validate predictions
            self._validate_predictions(compound)
            
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
        # Generate predictions
        for compound in compounds:
            try:
                self.ml_manager.predict_compound(
                    compound,
                    use_cache=self.config.use_cache,
                )
            except Exception as e:
                self.logger.error(
                    f"Failed predictions for {compound.name}: {str(e)}"
                )
        
        # Enrich predictions if enabled
        if self.config.enrich_predictions:
            for compound in compounds:
                try:
                    self._enrich_predictions(compound)
                except Exception as e:
                    self.logger.error(
                        f"Failed enrichment for {compound.name}: {str(e)}"
                    )
        
        # Validate predictions
        for compound in compounds:
            try:
                self._validate_predictions(compound)
            except Exception as e:
                self.logger.error(
                    f"Failed validation for {compound.name}: {str(e)}"
                )

    def _enrich_predictions(
        self,
        compound: EnhancedCompound,
    ) -> EnhancedCompound:
        """Enrich predictions with web data."""
        try:
            # Get web data
            web_data = compound.get_web_data()
            if not web_data:
                return compound
            
            # Update predictions
            for pred_type, pred in compound.get_predictions().items():
                if not pred:
                    continue
                
                # Get supporting data
                support = self._get_web_support(web_data, pred_type)
                if support:
                    pred["web_support"] = support
                    pred["confidence"] = max(
                        pred["confidence"],
                        support["confidence"]
                    )
            
            self.stats.enriched_compounds += 1
            return compound
            
        except Exception as e:
            self.logger.error(
                f"Failed to enrich predictions for {compound.name}: {str(e)}"
            )
            self.stats.enrichment_errors += 1
            return compound

    def _get_web_support(
        self,
        web_data: Dict,
        pred_type: str,
    ) -> Optional[Dict]:
        """Get supporting web data for prediction."""
        support = {
            "sources": [],
            "confidence": 0.0,
            "evidence": [],
        }
        
        # Check community data
        if "community" in web_data:
            community = web_data["community"]
            if pred_type in community:
                support["sources"].append("community")
                support["confidence"] = max(
                    support["confidence"],
                    community[pred_type].get("confidence", 0)
                )
                if "evidence" in community[pred_type]:
                    support["evidence"].extend(
                        community[pred_type]["evidence"]
                    )
        
        # Check literature data
        if "literature" in web_data:
            literature = web_data["literature"]
            if pred_type in literature:
                support["sources"].append("literature")
                support["confidence"] = max(
                    support["confidence"],
                    literature[pred_type].get("confidence", 0)
                )
                if "evidence" in literature[pred_type]:
                    support["evidence"].extend(
                        literature[pred_type]["evidence"]
                    )
        
        return support if support["sources"] else None

    def _validate_predictions(
        self,
        compound: EnhancedCompound,
    ) -> None:
        """Validate compound predictions."""
        predictions = compound.get_predictions()
        if not predictions:
            return
        
        # Check confidence threshold
        for pred_type, pred in predictions.items():
            if pred["confidence"] < self.config.min_confidence:
                self.logger.warning(
                    f"Low confidence prediction for {compound.name}: "
                    f"{pred_type} = {pred['confidence']:.3f}"
                )
        
        # Check web support
        if self.config.validate_web_data:
            for pred_type, pred in predictions.items():
                if "web_support" not in pred:
                    self.logger.warning(
                        f"No web support for {compound.name}: {pred_type}"
                    )
                elif len(pred["web_support"]["sources"]) < self.config.min_support:
                    self.logger.warning(
                        f"Limited web support for {compound.name}: "
                        f"{pred_type} = {len(pred['web_support']['sources'])} sources"
                    )

    def get_model_versions(self) -> Dict[str, str]:
        """Get model version information."""
        return self.ml_manager.get_model_versions()

    def clear_cache(self) -> None:
        """Clear prediction cache."""
        self.ml_manager.clear_cache()
