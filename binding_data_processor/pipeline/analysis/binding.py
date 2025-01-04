"""Binding analysis functionality.

This module provides the BindingAnalyzer class that handles:
1. Binding profile analysis
2. Target selectivity analysis
3. Binding pattern detection
4. Affinity analysis
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field
from datetime import datetime

from ...models import CompoundData


@dataclass
class BindingStats:
    """Binding analysis statistics."""
    
    # Analysis counts
    total_analyses: int = 0
    successful_analyses: int = 0
    failed_analyses: int = 0
    
    # Binding stats
    total_targets: int = 0
    primary_targets: int = 0
    target_families: Dict[str, int] = field(default_factory=dict)
    binding_types: Dict[str, int] = field(default_factory=dict)
    affinity_ranges: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    # Pattern stats
    total_patterns: int = 0
    pattern_types: Dict[str, int] = field(default_factory=dict)
    
    # Selectivity stats
    selectivity_ratios: List[float] = field(default_factory=list)
    
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
            "binding": {
                "total_targets": self.total_targets,
                "primary_targets": self.primary_targets,
                "target_families": self.target_families,
                "binding_types": self.binding_types,
                "affinity_ranges": self.affinity_ranges,
            },
            "patterns": {
                "total": self.total_patterns,
                "types": self.pattern_types,
            },
            "selectivity": {
                "ratios": self._get_selectivity_stats(),
            },
            "errors": self.errors,
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get analysis success rate."""
        if not self.total_analyses:
            return None
        return self.successful_analyses / self.total_analyses
    
    def _get_selectivity_stats(self) -> Dict[str, float]:
        """Get selectivity statistics."""
        if not self.selectivity_ratios:
            return {}
        
        ratios = sorted(self.selectivity_ratios)
        return {
            "min": ratios[0],
            "max": ratios[-1],
            "median": ratios[len(ratios) // 2],
            "mean": sum(ratios) / len(ratios),
        }


class BindingAnalyzer:
    """Analyzer for binding data."""

    def __init__(
        self,
        cache_dir: Optional[Union[str, Path]] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        """Initialize binding analyzer.
        
        Args:
            cache_dir: Optional directory for caching
            config: Optional analyzer configuration
        """
        self.logger = logging.getLogger(self.__class__.__name__)
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.config = config or {}
        
        # Get config values
        self.confidence_threshold = self.config.get(
            "confidence_threshold", 0.7
        )
        self.affinity_threshold = self.config.get(
            "affinity_threshold", 1000
        )
        
        # Initialize stats
        self.stats = BindingStats()

    def analyze(
        self,
        compound: CompoundData,
    ) -> CompoundData:
        """Analyze binding data.
        
        Args:
            compound: CompoundData instance to analyze
            
        Returns:
            Analyzed CompoundData instance
        """
        try:
            self.stats.total_analyses += 1
            
            # Get binding predictions
            binding_data = compound.get_prediction("receptor")
            if not binding_data:
                return compound
            
            # Analyze binding profile
            profile = self._analyze_binding_profile(binding_data)
            compound.binding_profile = profile
            
            # Update binding stats
            self._update_binding_stats(profile)
            
            # Analyze binding patterns
            patterns = self._analyze_binding_patterns(binding_data)
            compound.binding_patterns = patterns
            
            # Update pattern stats
            self._update_pattern_stats(patterns)
            
            # Analyze target selectivity
            selectivity = self._analyze_target_selectivity(binding_data)
            compound.target_selectivity = selectivity
            
            # Update selectivity stats
            self._update_selectivity_stats(selectivity)
            
            self.stats.successful_analyses += 1
            return compound
            
        except Exception as e:
            self.logger.error(
                f"Failed to analyze binding data for {compound.name}: {str(e)}"
            )
            self.stats.failed_analyses += 1
            self.stats.errors.append({
                "type": "binding_analysis_error",
                "compound": compound.name,
                "error": str(e),
                "timestamp": datetime.now().isoformat(),
            })
            raise

    def _analyze_binding_profile(
        self,
        binding_data: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Analyze binding profile.
        
        Args:
            binding_data: Binding prediction data
            
        Returns:
            Binding profile analysis
        """
        profile = {
            "targets": [],
            "primary_targets": [],
            "target_families": set(),
            "binding_types": {},
            "affinity_ranges": {},
        }
        
        # Process each target
        for target in binding_data.get("targets", []):
            # Check confidence
            confidence = target.get("confidence", 0)
            if confidence < self.confidence_threshold:
                continue
            
            # Get target info
            target_info = {
                "name": target["name"],
                "family": target.get("family", "unknown"),
                "type": target.get("type", "unknown"),
                "affinity": target.get("affinity", float("inf")),
                "confidence": confidence,
            }
            
            # Add to profile
            profile["targets"].append(target_info)
            
            # Check if primary target
            if target_info["affinity"] <= self.affinity_threshold:
                profile["primary_targets"].append(target_info)
            
            # Update families
            profile["target_families"].add(target_info["family"])
            
            # Update binding types
            binding_type = target_info["type"]
            if binding_type not in profile["binding_types"]:
                profile["binding_types"][binding_type] = 0
            profile["binding_types"][binding_type] += 1
            
            # Update affinity ranges
            family = target_info["family"]
            if family not in profile["affinity_ranges"]:
                profile["affinity_ranges"][family] = {
                    "min": float("inf"),
                    "max": float("-inf"),
                }
            ranges = profile["affinity_ranges"][family]
            affinity = target_info["affinity"]
            ranges["min"] = min(ranges["min"], affinity)
            ranges["max"] = max(ranges["max"], affinity)
        
        # Convert sets to lists
        profile["target_families"] = list(profile["target_families"])
        
        return profile

    def _analyze_binding_patterns(
        self,
        binding_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Analyze binding patterns.
        
        Args:
            binding_data: Binding prediction data
            
        Returns:
            List of binding patterns
        """
        patterns = []
        
        # Group targets by family
        families = {}
        for target in binding_data.get("targets", []):
            family = target.get("family")
            if not family:
                continue
            
            if family not in families:
                families[family] = []
            families[family].append(target)
        
        # Find patterns within families
        for family, targets in families.items():
            if len(targets) < 2:
                continue
            
            # Sort by affinity
            targets = sorted(
                targets,
                key=lambda x: x.get("affinity", float("inf")),
            )
            
            # Check for patterns
            patterns.extend(self._find_family_patterns(family, targets))
        
        return patterns

    def _find_family_patterns(
        self,
        family: str,
        targets: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """Find binding patterns within a target family.
        
        Args:
            family: Target family name
            targets: List of targets in family
            
        Returns:
            List of binding patterns
        """
        patterns = []
        
        # Check affinity range pattern
        min_affinity = targets[0].get("affinity", float("inf"))
        max_affinity = targets[-1].get("affinity", float("inf"))
        if max_affinity / min_affinity >= 10:
            patterns.append({
                "type": "affinity_range",
                "family": family,
                "min_affinity": min_affinity,
                "max_affinity": max_affinity,
                "ratio": max_affinity / min_affinity,
            })
        
        # Check binding type pattern
        types = {}
        for target in targets:
            binding_type = target.get("type", "unknown")
            if binding_type not in types:
                types[binding_type] = 0
            types[binding_type] += 1
        
        if len(types) > 1:
            patterns.append({
                "type": "binding_types",
                "family": family,
                "types": types,
            })
        
        return patterns

    def _analyze_target_selectivity(
        self,
        binding_data: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Analyze target selectivity.
        
        Args:
            binding_data: Binding prediction data
            
        Returns:
            Target selectivity analysis
        """
        selectivity = {
            "primary_target": None,
            "selectivity_ratios": {},
        }
        
        # Find primary target
        primary = None
        primary_affinity = float("inf")
        
        for target in binding_data.get("targets", []):
            affinity = target.get("affinity", float("inf"))
            confidence = target.get("confidence", 0)
            
            if (
                affinity < primary_affinity
                and confidence >= self.confidence_threshold
                and affinity <= self.affinity_threshold
            ):
                primary = target
                primary_affinity = affinity
        
        if not primary:
            return selectivity
        
        # Set primary target
        selectivity["primary_target"] = {
            "name": primary["name"],
            "family": primary.get("family", "unknown"),
            "type": primary.get("type", "unknown"),
            "affinity": primary_affinity,
            "confidence": primary.get("confidence", 0),
        }
        
        # Calculate selectivity ratios
        for target in binding_data.get("targets", []):
            if target == primary:
                continue
            
            affinity = target.get("affinity", float("inf"))
            if affinity == float("inf"):
                continue
            
            ratio = affinity / primary_affinity
            selectivity["selectivity_ratios"][target["name"]] = {
                "ratio": ratio,
                "affinity": affinity,
            }
        
        return selectivity

    def _update_binding_stats(
        self,
        profile: Dict[str, Any],
    ) -> None:
        """Update binding statistics."""
        # Update target counts
        self.stats.total_targets += len(profile["targets"])
        self.stats.primary_targets += len(profile["primary_targets"])
        
        # Update family stats
        for family in profile["target_families"]:
            if family not in self.stats.target_families:
                self.stats.target_families[family] = 0
            self.stats.target_families[family] += 1
        
        # Update binding type stats
        for type_, count in profile["binding_types"].items():
            if type_ not in self.stats.binding_types:
                self.stats.binding_types[type_] = 0
            self.stats.binding_types[type_] += count
        
        # Update affinity ranges
        for family, ranges in profile["affinity_ranges"].items():
            if family not in self.stats.affinity_ranges:
                self.stats.affinity_ranges[family] = {
                    "min": float("inf"),
                    "max": float("-inf"),
                }
            family_ranges = self.stats.affinity_ranges[family]
            family_ranges["min"] = min(
                family_ranges["min"],
                ranges["min"],
            )
            family_ranges["max"] = max(
                family_ranges["max"],
                ranges["max"],
            )

    def _update_pattern_stats(
        self,
        patterns: List[Dict[str, Any]],
    ) -> None:
        """Update pattern statistics."""
        self.stats.total_patterns += len(patterns)
        
        for pattern in patterns:
            pattern_type = pattern["type"]
            if pattern_type not in self.stats.pattern_types:
                self.stats.pattern_types[pattern_type] = 0
            self.stats.pattern_types[pattern_type] += 1

    def _update_selectivity_stats(
        self,
        selectivity: Dict[str, Any],
    ) -> None:
        """Update selectivity statistics."""
        for target_data in selectivity.get("selectivity_ratios", {}).values():
            self.stats.selectivity_ratios.append(target_data["ratio"])

    def get_info(self) -> Dict[str, Any]:
        """Get analyzer information."""
        return {
            "config": {
                "confidence_threshold": self.confidence_threshold,
                "affinity_threshold": self.affinity_threshold,
            },
            "stats": self.stats.to_dict(),
        }
