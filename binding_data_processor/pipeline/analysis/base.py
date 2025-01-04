"""Base analysis functionality.

This module provides:
1. Core analysis configuration
2. Analysis statistics tracking
3. Component coordination
4. Error handling
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field
from datetime import datetime

from ...models import CompoundData
from .binding import BindingAnalyzer
from .activity import ActivityAnalyzer
from .safety import SafetyAnalyzer
from .sar import SARAnalyzer
from .properties import PropertyAnalyzer


@dataclass
class AnalysisConfig:
    """Analysis configuration."""
    
    # Component configs
    binding: Dict[str, Any] = field(default_factory=lambda: {
        "enabled": True,
        "confidence_threshold": 0.7,
        "affinity_threshold": 1000,  # nM
    })
    
    activity: Dict[str, Any] = field(default_factory=lambda: {
        "enabled": True,
        "confidence_threshold": 0.7,
        "score_threshold": 0.5,
    })
    
    safety: Dict[str, Any] = field(default_factory=lambda: {
        "enabled": True,
        "confidence_threshold": 0.7,
        "toxicity_threshold": 0.5,
    })
    
    sar: Dict[str, Any] = field(default_factory=lambda: {
        "enabled": True,
        "similarity_threshold": 0.7,
        "pharmacophore_confidence": 0.7,
    })
    
    properties: Dict[str, Any] = field(default_factory=lambda: {
        "enabled": True,
        "ranges": {
            "molecular_weight": {"min": 0, "max": 2000},
            "logp": {"min": -10, "max": 10},
            "tpsa": {"min": 0, "max": 500},
            "hbd": {"min": 0, "max": 20},
            "hba": {"min": 0, "max": 20},
            "rotatable_bonds": {"min": 0, "max": 50},
        },
    })


@dataclass
class AnalysisStats:
    """Analysis statistics."""
    
    # Analysis counts
    total_analyses: int = 0
    successful_analyses: int = 0
    failed_analyses: int = 0
    
    # Component stats
    binding_stats: Dict[str, Any] = field(default_factory=dict)
    activity_stats: Dict[str, Any] = field(default_factory=dict)
    safety_stats: Dict[str, Any] = field(default_factory=dict)
    sar_stats: Dict[str, Any] = field(default_factory=dict)
    property_stats: Dict[str, Any] = field(default_factory=dict)
    
    # Error tracking
    errors: List[Dict[str, Any]] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary format."""
        return {
            "analyses": {
                "total": self.total_analyses,
                "successful": self.successful_analyses,
                "failed": self.failed_analyses,
                "success_rate": self._get_success_rate(),
            },
            "components": {
                "binding": self.binding_stats,
                "activity": self.activity_stats,
                "safety": self.safety_stats,
                "sar": self.sar_stats,
                "properties": self.property_stats,
            },
            "errors": self.errors,
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get analysis success rate."""
        if not self.total_analyses:
            return None
        return self.successful_analyses / self.total_analyses


class AnalysisManager:
    """Manager for compound analysis."""

    def __init__(
        self,
        cache_dir: Optional[Union[str, Path]] = None,
        config: Optional[AnalysisConfig] = None,
    ):
        """Initialize analysis manager.
        
        Args:
            cache_dir: Optional directory for caching
            config: Optional analysis configuration
        """
        self.logger = logging.getLogger(self.__class__.__name__)
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.config = config or AnalysisConfig()
        
        # Initialize analyzers
        self._init_analyzers()
        
        # Initialize stats
        self.stats = AnalysisStats()

    def _init_analyzers(self) -> None:
        """Initialize component analyzers."""
        try:
            # Initialize analyzers with component configs
            self.binding = BindingAnalyzer(
                cache_dir=self.cache_dir,
                config=self.config.binding,
            )
            self.activity = ActivityAnalyzer(
                cache_dir=self.cache_dir,
                config=self.config.activity,
            )
            self.safety = SafetyAnalyzer(
                cache_dir=self.cache_dir,
                config=self.config.safety,
            )
            self.sar = SARAnalyzer(
                cache_dir=self.cache_dir,
                config=self.config.sar,
            )
            self.properties = PropertyAnalyzer(
                cache_dir=self.cache_dir,
                config=self.config.properties,
            )
            
            self.logger.info("Successfully initialized analyzers")
            
        except Exception as e:
            self.logger.error(f"Failed to initialize analyzers: {str(e)}")
            raise

    def analyze_compound(
        self,
        compound: CompoundData,
    ) -> CompoundData:
        """Analyze compound data.
        
        Args:
            compound: CompoundData instance to analyze
            
        Returns:
            Analyzed CompoundData instance
        """
        try:
            self.stats.total_analyses += 1
            
            # Run component analyzers
            if self.config.binding["enabled"]:
                compound = self.binding.analyze(compound)
                self.stats.binding_stats.update(self.binding.stats)
            
            if self.config.activity["enabled"]:
                compound = self.activity.analyze(compound)
                self.stats.activity_stats.update(self.activity.stats)
            
            if self.config.safety["enabled"]:
                compound = self.safety.analyze(compound)
                self.stats.safety_stats.update(self.safety.stats)
            
            if self.config.sar["enabled"]:
                compound = self.sar.analyze(compound)
                self.stats.sar_stats.update(self.sar.stats)
            
            if self.config.properties["enabled"]:
                compound = self.properties.analyze(compound)
                self.stats.property_stats.update(self.properties.stats)
            
            self.stats.successful_analyses += 1
            return compound
            
        except Exception as e:
            self.logger.error(
                f"Failed to analyze compound {compound.name}: {str(e)}"
            )
            self.stats.failed_analyses += 1
            self.stats.errors.append({
                "type": "analysis_error",
                "compound": compound.name,
                "error": str(e),
                "timestamp": datetime.now().isoformat(),
            })
            raise

    def get_component_info(self) -> Dict[str, Any]:
        """Get component information."""
        return {
            "binding": self.binding.get_info(),
            "activity": self.activity.get_info(),
            "safety": self.safety.get_info(),
            "sar": self.sar.get_info(),
            "properties": self.properties.get_info(),
        }
