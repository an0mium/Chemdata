"""Analysis manager for pipeline.

This module provides the AnalysisManager class that handles:
1. Binding analysis
2. Activity analysis
3. Safety analysis
4. SAR analysis
5. Property analysis

The manager supports multiple analysis types:
- Binding profile analysis
- Activity pattern detection
- Safety assessment
- Structure-activity relationships
- Property relationships
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Set, Tuple
from dataclasses import dataclass, field
from datetime import datetime

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from rdkit.Chem.MolStandardize import rdMolStandardize

from ..models import CompoundData
from ..processors.structure import (
    StructureValidator,
    PharmacophoreDetector,
    SimilaritySearcher,
)


@dataclass
class AnalysisConfig:
    """Analysis configuration."""
    
    # Binding analysis
    analyze_binding: bool = True
    binding_confidence_threshold: float = 0.7
    binding_affinity_threshold: float = 1000  # nM
    
    # Activity analysis
    analyze_activity: bool = True
    activity_confidence_threshold: float = 0.7
    activity_score_threshold: float = 0.5
    
    # Safety analysis
    analyze_safety: bool = True
    safety_confidence_threshold: float = 0.7
    toxicity_threshold: float = 0.5
    
    # SAR analysis
    analyze_sar: bool = True
    similarity_threshold: float = 0.7
    pharmacophore_confidence: float = 0.7
    
    # Property analysis
    analyze_properties: bool = True
    property_ranges: Dict[str, Dict[str, float]] = field(default_factory=lambda: {
        "molecular_weight": {"min": 0, "max": 2000},
        "logp": {"min": -10, "max": 10},
        "tpsa": {"min": 0, "max": 500},
        "hbd": {"min": 0, "max": 20},
        "hba": {"min": 0, "max": 20},
        "rotatable_bonds": {"min": 0, "max": 50},
    })


@dataclass
class AnalysisStats:
    """Analysis statistics."""
    
    # Analysis counts
    total_analyses: int = 0
    successful_analyses: int = 0
    failed_analyses: int = 0
    
    # Component stats
    binding_stats: Dict[str, int] = field(default_factory=dict)
    activity_stats: Dict[str, int] = field(default_factory=dict)
    safety_stats: Dict[str, int] = field(default_factory=dict)
    sar_stats: Dict[str, int] = field(default_factory=dict)
    property_stats: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
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
        """Initialize analyzers."""
        try:
            # Structure analyzers
            self.structure_validator = StructureValidator()
            self.pharmacophore_detector = PharmacophoreDetector()
            self.similarity_searcher = SimilaritySearcher()
            
            # Standardization
            self.standardizer = rdMolStandardize.Standardizer()
            
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
            
            # Binding analysis
            if self.config.analyze_binding:
                self._analyze_binding(compound)
            
            # Activity analysis
            if self.config.analyze_activity:
                self._analyze_activity(compound)
            
            # Safety analysis
            if self.config.analyze_safety:
                self._analyze_safety(compound)
            
            # SAR analysis
            if self.config.analyze_sar:
                self._analyze_sar(compound)
            
            # Property analysis
            if self.config.analyze_properties:
                self._analyze_properties(compound)
            
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

    def _analyze_binding(self, compound: CompoundData) -> None:
        """Analyze binding data."""
        try:
            # Get binding predictions
            binding_data = compound.get_prediction("receptor")
            if not binding_data:
                return
            
            # Analyze strongest binding
            strongest = self._find_strongest_binding(binding_data)
            if strongest:
                self._update_binding_stats("strongest", strongest)
            
            # Analyze binding patterns
            patterns = self._find_binding_patterns(binding_data)
            for pattern in patterns:
                self._update_binding_stats("pattern", pattern)
            
            # Analyze target selectivity
            selectivity = self._analyze_selectivity(binding_data)
            self._update_binding_stats("selectivity", selectivity)
            
        except Exception as e:
            self.logger.error(f"Binding analysis error: {str(e)}")
            raise

    def _analyze_activity(self, compound: CompoundData) -> None:
        """Analyze activity data."""
        try:
            # Get activity predictions
            activity_data = compound.get_prediction("activity_ensemble")
            if not activity_data:
                return
            
            # Analyze primary activity
            primary = self._find_primary_activity(activity_data)
            if primary:
                self._update_activity_stats("primary", primary)
            
            # Analyze activity patterns
            patterns = self._find_activity_patterns(activity_data)
            for pattern in patterns:
                self._update_activity_stats("pattern", pattern)
            
            # Analyze mechanisms
            mechanisms = self._analyze_mechanisms(activity_data)
            for mechanism in mechanisms:
                self._update_activity_stats("mechanism", mechanism)
            
        except Exception as e:
            self.logger.error(f"Activity analysis error: {str(e)}")
            raise

    def _analyze_safety(self, compound: CompoundData) -> None:
        """Analyze safety data."""
        try:
            # Get safety predictions
            safety_data = compound.get_prediction("safety_ensemble")
            if not safety_data:
                return
            
            # Analyze toxicity risks
            risks = self._find_toxicity_risks(safety_data)
            for risk in risks:
                self._update_safety_stats("risk", risk)
            
            # Analyze safety patterns
            patterns = self._find_safety_patterns(safety_data)
            for pattern in patterns:
                self._update_safety_stats("pattern", pattern)
            
            # Analyze interactions
            interactions = self._analyze_interactions(safety_data)
            for interaction in interactions:
                self._update_safety_stats("interaction", interaction)
            
        except Exception as e:
            self.logger.error(f"Safety analysis error: {str(e)}")
            raise

    def _analyze_sar(self, compound: CompoundData) -> None:
        """Analyze structure-activity relationships."""
        try:
            # Get structure
            mol = Chem.MolFromSmiles(compound.smiles)
            if not mol:
                return
            
            # Standardize structure
            mol = self.standardizer.standardize(mol)
            
            # Detect pharmacophores
            pharmacophores = self.pharmacophore_detector.detect(mol)
            for pharm in pharmacophores:
                self._update_sar_stats("pharmacophore", pharm)
            
            # Find similar compounds
            similar = self.similarity_searcher.search(
                mol,
                threshold=self.config.similarity_threshold,
            )
            for sim in similar:
                self._update_sar_stats("similar", sim)
            
            # Analyze activity cliffs
            cliffs = self._find_activity_cliffs(compound, similar)
            for cliff in cliffs:
                self._update_sar_stats("cliff", cliff)
            
        except Exception as e:
            self.logger.error(f"SAR analysis error: {str(e)}")
            raise

    def _analyze_properties(self, compound: CompoundData) -> None:
        """Analyze chemical properties."""
        try:
            # Get structure
            mol = Chem.MolFromSmiles(compound.smiles)
            if not mol:
                return
            
            # Calculate properties
            for name, ranges in self.config.property_ranges.items():
                try:
                    value = getattr(Descriptors, f"Calc{name}")(mol)
                    if ranges["min"] <= value <= ranges["max"]:
                        if name not in self.stats.property_stats:
                            self.stats.property_stats[name] = {
                                "count": 0,
                                "sum": 0,
                                "min": float("inf"),
                                "max": float("-inf"),
                            }
                        stats = self.stats.property_stats[name]
                        stats["count"] += 1
                        stats["sum"] += value
                        stats["min"] = min(stats["min"], value)
                        stats["max"] = max(stats["max"], value)
                except Exception as e:
                    self.logger.debug(f"Property calculation error: {str(e)}")
            
        except Exception as e:
            self.logger.error(f"Property analysis error: {str(e)}")
            raise

    def _find_strongest_binding(
        self,
        binding_data: Dict[str, Any],
    ) -> Optional[Dict[str, Any]]:
        """Find strongest binding target."""
        if not binding_data.get("targets"):
            return None
        
        strongest = None
        strongest_affinity = float("inf")
        
        for target in binding_data["targets"]:
            affinity = target.get("affinity", float("inf"))
            confidence = target.get("confidence", 0)
            
            if (
                affinity < strongest_affinity
                and confidence >= self.config.binding_confidence_threshold
                and affinity <= self.config.binding_affinity_threshold
            ):
                strongest = target
                strongest_affinity = affinity
        
        return strongest

    def _find_binding_patterns(
        self,
        binding_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Find binding patterns."""
        patterns = []
        if not binding_data.get("targets"):
            return patterns
        
        # Group by target family
        families = {}
        for target in binding_data["targets"]:
            family = target.get("family")
            if family:
                if family not in families:
                    families[family] = []
                families[family].append(target)
        
        # Find patterns within families
        for family, targets in families.items():
            if len(targets) >= 2:
                patterns.append({
                    "type": "family_pattern",
                    "family": family,
                    "targets": targets,
                })
        
        return patterns

    def _analyze_selectivity(
        self,
        binding_data: Dict[str, Any],
    ) -> Optional[Dict[str, Any]]:
        """Analyze target selectivity."""
        if not binding_data.get("targets"):
            return None
        
        # Calculate selectivity ratios
        ratios = []
        primary = self._find_strongest_binding(binding_data)
        if not primary:
            return None
        
        primary_affinity = primary["affinity"]
        for target in binding_data["targets"]:
            if target != primary:
                affinity = target.get("affinity", float("inf"))
                if affinity != float("inf"):
                    ratio = affinity / primary_affinity
                    ratios.append({
                        "target": target["name"],
                        "ratio": ratio,
                    })
        
        return {
            "primary": primary["name"],
            "ratios": sorted(ratios, key=lambda x: x["ratio"]),
        }

    def _find_primary_activity(
        self,
        activity_data: Dict[str, Any],
    ) -> Optional[Dict[str, Any]]:
        """Find primary activity type."""
        if not activity_data.get("activities"):
            return None
        
        primary = None
        highest_score = 0
        
        for activity in activity_data["activities"]:
            score = activity.get("score", 0)
            confidence = activity.get("confidence", 0)
            
            if (
                score > highest_score
                and confidence >= self.config.activity_confidence_threshold
                and score >= self.config.activity_score_threshold
            ):
                primary = activity
                highest_score = score
        
        return primary

    def _find_activity_patterns(
        self,
        activity_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Find activity patterns."""
        patterns = []
        if not activity_data.get("activities"):
            return patterns
        
        # Group by activity class
        classes = {}
        for activity in activity_data["activities"]:
            cls = activity.get("class")
            if cls:
                if cls not in classes:
                    classes[cls] = []
                classes[cls].append(activity)
        
        # Find patterns within classes
        for cls, activities in classes.items():
            if len(activities) >= 2:
                patterns.append({
                    "type": "class_pattern",
                    "class": cls,
                    "activities": activities,
                })
        
        return patterns

    def _analyze_mechanisms(
        self,
        activity_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Analyze activity mechanisms."""
        mechanisms = []
        if not activity_data.get("activities"):
            return mechanisms
        
        # Extract mechanisms
        for activity in activity_data["activities"]:
            mechanism = activity.get("mechanism")
            if mechanism:
                mechanisms.append({
                    "type": "mechanism",
                    "activity": activity["name"],
                    "mechanism": mechanism,
                })
        
        return mechanisms

    def _find_toxicity_risks(
        self,
        safety_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Find toxicity risks."""
        risks = []
        if not safety_data.get("risks"):
            return risks
        
        for risk in safety_data["risks"]:
            score = risk.get("score", 0)
            confidence = risk.get("confidence", 0)
            
            if (
                score >= self.config.toxicity_threshold
                and confidence >= self.config.safety_confidence_threshold
            ):
                risks.append(risk)
        
        return risks

    def _find_safety_patterns(
        self,
        safety_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Find safety patterns."""
        patterns = []
        if not safety_data.get("risks"):
            return patterns
        
        # Group by risk type
        types = {}
        for risk in safety_data["risks"]:
            type_ = risk.get("type")
            if type_:
                if type_ not in types:
                    types[type_] = []
                types[type_].append(risk)
        
        # Find patterns within types
        for type_, risks in types.items():
            if len(risks) >= 2:
                patterns.append({
                    "type": "risk_pattern",
                    "risk_type": type_,
                    "risks": risks,
                })
        
        return patterns

    def _analyze_interactions(
        self,
        safety_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Analyze drug interactions."""
        interactions = []
        if not safety_data.get("interactions"):
            return interactions
        
        # Extract interactions
        for interaction in safety_data["interactions"]:
            severity = interaction.get("severity", 0)
            confidence = interaction.get("confidence", 0)
            
            if (
                severity >= self.config.toxicity_threshold
                and confidence >= self.config.safety_confidence_threshold
            ):
                interactions.append(interaction)
        
        return interactions

    def _find_activity_cliffs(
        self,
        compound: CompoundData,
        similar: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """Find activity cliffs."""
        cliffs = []
        
        # Get compound activity
        activity = compound.get_prediction("activity_ensemble")
        if not activity:
            return cliffs
        
        # Compare with similar compounds
        for sim in similar:
            sim_activity = sim.get("activity")
            if sim_activity:
                # Calculate activity difference
                diff = abs(
                    activity.get("score", 0) -
                    sim_activity.get("score", 0)
                )
                
                if diff >= 0.5:  # Significant activity difference
                    cliffs.append({
                        "compound": sim["name"],
                        "similarity": sim["similarity"],
                        "activity_diff": diff,
                    })
        
        return cliffs

    def _update_binding_stats(
        self,
        type_: str,
        data: Dict[str, Any],
    ) -> None:
        """Update binding statistics."""
        if type_ not in self.stats.binding_stats:
            self.stats.binding_stats[type_] = 0
        self.stats.binding_stats[type_] += 1

    def _update_activity_stats(
        self,
        type_: str,
        data: Dict[str, Any],
    ) -> None:
        """Update activity statistics."""
        if type_ not in self.stats.activity_stats:
            self.stats.activity_stats[type_] = 0
        self.stats.activity_stats[type_] += 1

    def _update_safety_stats(
        self,
        type_: str,
        data: Dict[str, Any],
    ) -> None:
        """Update safety statistics."""
        if type_ not in self.stats.safety_stats:
            self.stats.safety_stats[type_] = 0
        self.stats.safety_stats[type_] += 1

    def _update_sar_stats(
        self,
        type_: str,
        data: Dict[str, Any],
    ) -> None:
        """Update SAR statistics."""
        if type_ not in self.stats.sar_stats:
            self.stats.sar_stats[type_] = 0
        self.stats.sar_stats[type_] += 1
