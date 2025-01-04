"""Activity analysis functionality.

This module provides the ActivityAnalyzer class that handles:
1. Activity type analysis
2. Activity pattern detection
3. Mechanism analysis
4. Effect analysis
5. Integration with property and binding analysis
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field
from datetime import datetime

from ...models import CompoundData


@dataclass
class ActivityStats:
    """Activity analysis statistics."""
    
    # Analysis counts
    total_analyses: int = 0
    successful_analyses: int = 0
    failed_analyses: int = 0
    
    # Activity stats
    total_activities: int = 0
    primary_activities: int = 0
    activity_types: Dict[str, int] = field(default_factory=dict)
    activity_scores: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    # Pattern stats
    total_patterns: int = 0
    pattern_types: Dict[str, int] = field(default_factory=dict)
    
    # Mechanism stats
    total_mechanisms: int = 0
    mechanism_types: Dict[str, int] = field(default_factory=dict)
    
    # Effect stats
    total_effects: int = 0
    effect_types: Dict[str, int] = field(default_factory=dict)
    
    # Property-activity stats
    property_correlations: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    # Binding-activity stats
    binding_correlations: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
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
            "activity": {
                "total": self.total_activities,
                "primary": self.primary_activities,
                "types": self.activity_types,
                "scores": self.activity_scores,
            },
            "patterns": {
                "total": self.total_patterns,
                "types": self.pattern_types,
            },
            "mechanisms": {
                "total": self.total_mechanisms,
                "types": self.mechanism_types,
            },
            "effects": {
                "total": self.total_effects,
                "types": self.effect_types,
            },
            "correlations": {
                "property": self.property_correlations,
                "binding": self.binding_correlations,
            },
            "errors": self.errors,
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get analysis success rate."""
        if not self.total_analyses:
            return None
        return self.successful_analyses / self.total_analyses


class ActivityAnalyzer:
    """Analyzer for activity data."""

    def __init__(
        self,
        cache_dir: Optional[Union[str, Path]] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        """Initialize activity analyzer.
        
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
        self.score_threshold = self.config.get(
            "score_threshold", 0.5
        )
        self.correlation_threshold = self.config.get(
            "correlation_threshold", 0.3
        )
        
        # Initialize stats
        self.stats = ActivityStats()

    def analyze(
        self,
        compound: CompoundData,
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
    ) -> CompoundData:
        """Analyze activity data.
        
        Args:
            compound: CompoundData instance to analyze
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            
        Returns:
            Analyzed CompoundData instance
        """
        try:
            self.stats.total_analyses += 1
            
            # Get activity predictions
            activity_data = compound.get_prediction("activity_ensemble")
            if not activity_data:
                return compound
            
            # Analyze activity profile
            profile = self._analyze_activity_profile(
                activity_data,
                property_data,
                binding_data,
            )
            compound.activity_profile = profile
            
            # Update activity stats
            self._update_activity_stats(profile)
            
            # Analyze activity patterns
            patterns = self._analyze_activity_patterns(
                activity_data,
                property_data,
                binding_data,
            )
            compound.activity_patterns = patterns
            
            # Update pattern stats
            self._update_pattern_stats(patterns)
            
            # Analyze mechanisms
            mechanisms = self._analyze_mechanisms(
                activity_data,
                property_data,
                binding_data,
            )
            compound.activity_mechanisms = mechanisms
            
            # Update mechanism stats
            self._update_mechanism_stats(mechanisms)
            
            # Analyze effects
            effects = self._analyze_effects(
                activity_data,
                property_data,
                binding_data,
            )
            compound.activity_effects = effects
            
            # Update effect stats
            self._update_effect_stats(effects)
            
            # Analyze property correlations
            if property_data:
                correlations = self._analyze_property_correlations(
                    activity_data,
                    property_data,
                )
                compound.property_activity_correlations = correlations
                
                # Update correlation stats
                self._update_property_correlation_stats(correlations)
            
            # Analyze binding correlations
            if binding_data:
                correlations = self._analyze_binding_correlations(
                    activity_data,
                    binding_data,
                )
                compound.binding_activity_correlations = correlations
                
                # Update correlation stats
                self._update_binding_correlation_stats(correlations)
            
            self.stats.successful_analyses += 1
            return compound
            
        except Exception as e:
            self.logger.error(
                f"Failed to analyze activity data for {compound.name}: {str(e)}"
            )
            self.stats.failed_analyses += 1
            self.stats.errors.append({
                "type": "activity_analysis_error",
                "compound": compound.name,
                "error": str(e),
                "timestamp": datetime.now().isoformat(),
            })
            raise

    def _analyze_activity_profile(
        self,
        activity_data: Dict[str, Any],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Analyze activity profile.
        
        Args:
            activity_data: Activity prediction data
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            
        Returns:
            Activity profile analysis
        """
        profile = {
            "activities": [],
            "primary_activities": [],
            "activity_types": set(),
            "activity_scores": {},
        }
        
        # Process each activity
        for activity in activity_data.get("activities", []):
            # Check confidence
            confidence = activity.get("confidence", 0)
            if confidence < self.confidence_threshold:
                continue
            
            # Get activity info
            activity_info = {
                "name": activity["name"],
                "type": activity.get("type", "unknown"),
                "score": activity.get("score", 0),
                "confidence": confidence,
                "mechanism": activity.get("mechanism"),
                "effects": activity.get("effects", []),
            }
            
            # Enhance with property data
            if property_data:
                activity_info["property_effects"] = self._get_property_effects(
                    activity_info,
                    property_data,
                )
            
            # Enhance with binding data
            if binding_data:
                activity_info["binding_effects"] = self._get_binding_effects(
                    activity_info,
                    binding_data,
                )
            
            # Add to profile
            profile["activities"].append(activity_info)
            
            # Check if primary activity
            if activity_info["score"] >= self.score_threshold:
                profile["primary_activities"].append(activity_info)
            
            # Update types
            profile["activity_types"].add(activity_info["type"])
            
            # Update scores
            activity_type = activity_info["type"]
            if activity_type not in profile["activity_scores"]:
                profile["activity_scores"][activity_type] = {
                    "min": float("inf"),
                    "max": float("-inf"),
                    "sum": 0,
                    "count": 0,
                }
            scores = profile["activity_scores"][activity_type]
            score = activity_info["score"]
            scores["min"] = min(scores["min"], score)
            scores["max"] = max(scores["max"], score)
            scores["sum"] += score
            scores["count"] += 1
        
        # Convert sets to lists
        profile["activity_types"] = list(profile["activity_types"])
        
        return profile

    def _get_property_effects(
        self,
        activity_info: Dict[str, Any],
        property_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Get property-based effects on activity.
        
        Args:
            activity_info: Activity information
            property_data: Property analysis data
            
        Returns:
            List of property effects
        """
        effects = []
        
        # Check size effects
        if property_data.get("molecular_weight", 0) > 500:
            effects.append({
                "type": "size",
                "effect": "reduced_activity",
                "confidence": 0.8,
            })
        
        # Check polarity effects
        logp = property_data.get("logp", 0)
        if not (-2 <= logp <= 5):
            effects.append({
                "type": "polarity",
                "effect": "reduced_activity",
                "confidence": 0.7,
            })
        
        # Check BBB effects
        if property_data.get("bbb", {}).get("pass"):
            effects.append({
                "type": "bbb",
                "effect": "cns_activity",
                "confidence": 0.9,
            })
        
        return effects

    def _get_binding_effects(
        self,
        activity_info: Dict[str, Any],
        binding_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Get binding-based effects on activity.
        
        Args:
            activity_info: Activity information
            binding_data: Binding analysis data
            
        Returns:
            List of binding effects
        """
        effects = []
        
        # Check target effects
        if binding_data.get("primary_targets"):
            effects.append({
                "type": "target",
                "effect": "direct_binding",
                "confidence": 0.9,
                "targets": binding_data["primary_targets"],
            })
        
        # Check selectivity effects
        if binding_data.get("selectivity_ratios", {}):
            effects.append({
                "type": "selectivity",
                "effect": "target_selective",
                "confidence": 0.8,
                "ratios": binding_data["selectivity_ratios"],
            })
        
        return effects

    def _analyze_activity_patterns(
        self,
        activity_data: Dict[str, Any],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Analyze activity patterns.
        
        Args:
            activity_data: Activity prediction data
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            
        Returns:
            List of activity patterns
        """
        patterns = []
        
        # Group activities by type
        types = {}
        for activity in activity_data.get("activities", []):
            activity_type = activity.get("type")
            if not activity_type:
                continue
            
            if activity_type not in types:
                types[activity_type] = []
            types[activity_type].append(activity)
        
        # Find patterns within types
        for activity_type, activities in types.items():
            if len(activities) < 2:
                continue
            
            # Sort by score
            activities = sorted(
                activities,
                key=lambda x: x.get("score", 0),
                reverse=True,
            )
            
            # Check for patterns
            patterns.extend(
                self._find_type_patterns(
                    activity_type,
                    activities,
                    property_data,
                    binding_data,
                )
            )
        
        return patterns

    def _find_type_patterns(
        self,
        activity_type: str,
        activities: List[Dict[str, Any]],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Find activity patterns within an activity type.
        
        Args:
            activity_type: Activity type name
            activities: List of activities of this type
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            
        Returns:
            List of activity patterns
        """
        patterns = []
        
        # Check score range pattern
        max_score = activities[0].get("score", 0)
        min_score = activities[-1].get("score", 0)
        if max_score - min_score >= 0.5:  # Significant score range
            patterns.append({
                "type": "score_range",
                "activity_type": activity_type,
                "min_score": min_score,
                "max_score": max_score,
                "range": max_score - min_score,
            })
        
        # Check mechanism pattern
        mechanisms = {}
        for activity in activities:
            mechanism = activity.get("mechanism", "unknown")
            if mechanism not in mechanisms:
                mechanisms[mechanism] = 0
            mechanisms[mechanism] += 1
        
        if len(mechanisms) > 1:
            patterns.append({
                "type": "mechanisms",
                "activity_type": activity_type,
                "mechanisms": mechanisms,
            })
        
        # Check effect pattern
        effects = {}
        for activity in activities:
            for effect in activity.get("effects", []):
                if effect not in effects:
                    effects[effect] = 0
                effects[effect] += 1
        
        if len(effects) > 1:
            patterns.append({
                "type": "effects",
                "activity_type": activity_type,
                "effects": effects,
            })
        
        # Check property-based patterns
        if property_data:
            property_patterns = self._find_property_patterns(
                activity_type,
                activities,
                property_data,
            )
            patterns.extend(property_patterns)
        
        # Check binding-based patterns
        if binding_data:
            binding_patterns = self._find_binding_patterns(
                activity_type,
                activities,
                binding_data,
            )
            patterns.extend(binding_patterns)
        
        return patterns

    def _find_property_patterns(
        self,
        activity_type: str,
        activities: List[Dict[str, Any]],
        property_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Find property-based activity patterns.
        
        Args:
            activity_type: Activity type name
            activities: List of activities
            property_data: Property analysis data
            
        Returns:
            List of property-based patterns
        """
        patterns = []
        
        # Check size-activity pattern
        if property_data.get("molecular_weight", 0) > 500:
            high_scores = [
                a for a in activities
                if a.get("score", 0) >= self.score_threshold
            ]
            if high_scores:
                patterns.append({
                    "type": "size_activity",
                    "activity_type": activity_type,
                    "molecular_weight": property_data["molecular_weight"],
                    "high_activity_count": len(high_scores),
                })
        
        # Check polarity-activity pattern
        logp = property_data.get("logp", 0)
        if not (-2 <= logp <= 5):
            high_scores = [
                a for a in activities
                if a.get("score", 0) >= self.score_threshold
            ]
            if high_scores:
                patterns.append({
                    "type": "polarity_activity",
                    "activity_type": activity_type,
                    "logp": logp,
                    "high_activity_count": len(high_scores),
                })
        
        return patterns

    def _find_binding_patterns(
        self,
        activity_type: str,
        activities: List[Dict[str, Any]],
        binding_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Find binding-based activity patterns.
        
        Args:
            activity_type: Activity type name
            activities: List of activities
            binding_data: Binding analysis data
            
        Returns:
            List of binding-based patterns
        """
        patterns = []
        
        # Check target-activity pattern
        if binding_data.get("primary_targets"):
            high_scores = [
                a for a in activities
                if a.get("score", 0) >= self.score_threshold
            ]
            if high_scores:
                patterns.append({
                    "type": "target_activity",
                    "activity_type": activity_type,
                    "targets": binding_data["primary_targets"],
                    "high_activity_count": len(high_scores),
                })
        
        # Check selectivity-activity pattern
        if binding_data.get("selectivity_ratios", {}):
            high_scores = [
                a for a in activities
                if a.get("score", 0) >= self.score_threshold
            ]
            if high_scores:
                patterns.append({
                    "type": "selectivity_activity",
                    "activity_type": activity_type,
                    "selectivity": binding_data["selectivity_ratios"],
                    "high_activity_count": len(high_scores),
                })
        
        return patterns

    def _analyze_mechanisms(
        self,
        activity_data: Dict[str, Any],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Analyze activity mechanisms.
        
        Args:
            activity_data: Activity prediction data
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            
        Returns:
            List of activity mechanisms
        """
        mechanisms = []
        
        # Extract mechanisms
        for activity in activity_data.get("activities", []):
            mechanism = activity.get("mechanism")
            if not mechanism:
                continue
            
            confidence = activity.get("confidence", 0)
            if confidence < self.confidence_threshold:
                continue
            
            mechanism_info = {
                "name": mechanism,
                "type": activity.get("type", "unknown"),
                "score": activity.get("score", 0),
                "confidence": confidence,
                "effects": activity.get("effects", []),
            }
            
            # Enhance with property data
            if property_data:
                mechanism_info["property_effects"] = self._get_property_effects(
                    activity,
                    property_data,
                )
            
            # Enhance with binding data
            if binding_data:
                mechanism_info["binding_effects"] = self._get_binding_effects(
                    activity,
                    binding_data,
                )
            
            mechanisms.append(mechanism_info)
        
        return mechanisms

    def _analyze_effects(
        self,
        activity_data: Dict[str, Any],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Analyze activity effects.
        
        Args:
            activity_data: Activity prediction data
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            
        Returns:
            List of activity effects
        """
        effects = []
        
        # Extract effects
        for activity in activity_data.get("activities", []):
            activity_effects = activity.get("effects", [])
            if not activity_effects:
                continue
            
            confidence = activity.get("confidence", 0)
            if confidence < self.confidence_threshold:
                continue
            
            for effect in activity_effects:
                effect_info = {
                    "name": effect,
                    "type": activity.get("type", "unknown"),
                    "score": activity.get("score", 0),
                    "confidence": confidence,
                    "mechanism": activity.get("mechanism"),
                }
                
                # Enhance with property data
                if property_data:
                    effect_info["property_effects"] = self._get_property_effects(
                        activity,
                        property_data,
                    )
                
                # Enhance with binding data
                if binding_data:
                    effect_info["binding_effects"] = self._get_binding_effects(
                        activity,
                        binding_data,
                    )
                
                effects.append(effect_info)
        
        return effects

    def _analyze_property_correlations(
        self,
        activity_data: Dict[str, Any],
        property_data: Dict[str, Any],
    ) -> Dict[str, Dict[str, float]]:
        """Analyze property-activity correlations.
        
        Args:
            activity_data: Activity prediction data
            property_data: Property analysis data
            
        Returns:
            Dictionary of property-activity correlations
        """
        correlations = {}
        
        # Get activity scores by type
        activity_scores = {}
        for activity in activity_data.get("activities", []):
            activity_type = activity.get("type", "unknown")
            if activity_type not in activity_scores:
                activity_scores[activity_type] = []
            activity_scores[activity_type].append(
                activity.get("score", 0)
            )
        
        # Calculate correlations
        for prop_name, prop_value in property_data.items():
            if not isinstance(prop_value, (int, float)):
                continue
            
            correlations[prop_name] = {}
            for activity_type, scores in activity_scores.items():
                correlation = self._calculate_correlation(
                    [prop_value] * len(scores),
                    scores,
                )
                if abs(correlation) >= self.correlation_threshold:
                    correlations[prop_name][activity_type] = correlation
        
        return correlations

    def _analyze_binding_correlations(
        self,
        activity_data: Dict[str, Any],
        binding_data: Dict[str, Any],
    ) -> Dict[str, Dict[str, float]]:
        """Analyze binding-activity correlations.
        
        Args:
            activity_data: Activity prediction data
            binding_data: Binding analysis data
            
        Returns:
            Dictionary of binding-activity correlations
        """
        correlations = {}
        
        # Get activity scores by type
        activity_scores = {}
        for activity in activity_data.get("activities", []):
            activity_type = activity.get("type", "unknown")
            if activity_type not in activity_scores:
                activity_scores[activity_type] = []
            activity_scores[activity_type].append(
                activity.get("score", 0)
            )
        
        # Calculate correlations for each target
        for target in binding_data.get("targets", []):
            target_name = target["name"]
            target_affinity = target.get("affinity", float("inf"))
            
            correlations[target_name] = {}
            for activity_type, scores in activity_scores.items():
                correlation = self._calculate_correlation(
                    [target_affinity] * len(scores),
                    scores,
                )
                if abs(correlation) >= self.correlation_threshold:
                    correlations[target_name][activity_type] = correlation
        
        return correlations

    def _calculate_correlation(
        self,
        x: List[float],
        y: List[float],
    ) -> float:
        """Calculate correlation coefficient.
        
        Args:
            x: First list of values
            y: Second list of values
            
        Returns:
            Correlation coefficient
        """
        if len(x) != len(y) or len(x) < 2:
            return 0.0
        
        n = len(x)
        sum_x = sum(x)
        sum_y = sum(y)
        sum_xy = sum(xi * yi for xi, yi in zip(x, y))
        sum_x2 = sum(xi * xi for xi in x)
        sum_y2 = sum(yi * yi for yi in y)
        
        numerator = n * sum_xy - sum_x * sum_y
        denominator = ((n * sum_x2 - sum_x * sum_x) *
                      (n * sum_y2 - sum_y * sum_y)) ** 0.5
        
        if denominator == 0:
            return 0.0
        
        return numerator / denominator

    def _update_activity_stats(
        self,
        profile: Dict[str, Any],
    ) -> None:
        """Update activity statistics."""
        # Update activity counts
        self.stats.total_activities += len(profile["activities"])
        self.stats.primary_activities += len(profile["primary_activities"])
        
        # Update type stats
        for activity_type in profile["activity_types"]:
            if activity_type not in self.stats.activity_types:
                self.stats.activity_types[activity_type] = 0
            self.stats.activity_types[activity_type] += 1
        
        # Update score stats
        for type_, scores in profile["activity_scores"].items():
            if type_ not in self.stats.activity_scores:
                self.stats.activity_scores[type_] = {
                    "min": float("inf"),
                    "max": float("-inf"),
                    "sum": 0,
                    "count": 0,
                }
            type_scores = self.stats.activity_scores[type_]
            type_scores["min"] = min(
                type_scores["min"],
                scores["min"],
            )
            type_scores["max"] = max(
                type_scores["max"],
                scores["max"],
            )
            type_scores["sum"] += scores["sum"]
            type_scores["count"] += scores["count"]

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

    def _update_mechanism_stats(
        self,
        mechanisms: List[Dict[str, Any]],
    ) -> None:
        """Update mechanism statistics."""
        self.stats.total_mechanisms += len(mechanisms)
        
        for mechanism in mechanisms:
            mechanism_type = mechanism["type"]
            if mechanism_type not in self.stats.mechanism_types:
                self.stats.mechanism_types[mechanism_type] = 0
            self.stats.mechanism_types[mechanism_type] += 1

    def _update_effect_stats(
        self,
        effects: List[Dict[str, Any]],
    ) -> None:
        """Update effect statistics."""
        self.stats.total_effects += len(effects)
        
        for effect in effects:
            effect_type = effect["type"]
            if effect_type not in self.stats.effect_types:
                self.stats.effect_types[effect_type] = 0
            self.stats.effect_types[effect_type] += 1

    def _update_property_correlation_stats(
        self,
        correlations: Dict[str, Dict[str, float]],
    ) -> None:
        """Update property correlation statistics."""
        for prop_name, prop_corrs in correlations.items():
            if prop_name not in self.stats.property_correlations:
                self.stats.property_correlations[prop_name] = {}
            
            for activity_type, correlation in prop_corrs.items():
                if activity_type not in self.stats.property_correlations[prop_name]:
                    self.stats.property_correlations[prop_name][activity_type] = correlation

    def _update_binding_correlation_stats(
        self,
        correlations: Dict[str, Dict[str, float]],
    ) -> None:
        """Update binding correlation statistics."""
        for target_name, target_corrs in correlations.items():
            if target_name not in self.stats.binding_correlations:
                self.stats.binding_correlations[target_name] = {}
            
            for activity_type, correlation in target_corrs.items():
                if activity_type not in self.stats.binding_correlations[target_name]:
                    self.stats.binding_correlations[target_name][activity_type] = correlation

    def get_info(self) -> Dict[str, Any]:
        """Get analyzer information."""
        return {
            "config": {
                "confidence_threshold": self.confidence_threshold,
                "score_threshold": self.score_threshold,
                "correlation_threshold": self.correlation_threshold,
            },
            "stats": self.stats.to_dict(),
        }
