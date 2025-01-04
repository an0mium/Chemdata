"""Safety analysis functionality.

This module provides the SafetyAnalyzer class that handles:
1. Toxicity analysis
2. Risk assessment
3. Safety pattern detection
4. Interaction analysis
5. Integration with property, binding, and activity analysis
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field
from datetime import datetime

from ...models import CompoundData


@dataclass
class SafetyStats:
    """Safety analysis statistics."""
    
    # Analysis counts
    total_analyses: int = 0
    successful_analyses: int = 0
    failed_analyses: int = 0
    
    # Toxicity stats
    total_risks: int = 0
    high_risks: int = 0
    risk_types: Dict[str, int] = field(default_factory=dict)
    risk_scores: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    # Pattern stats
    total_patterns: int = 0
    pattern_types: Dict[str, int] = field(default_factory=dict)
    
    # Interaction stats
    total_interactions: int = 0
    interaction_types: Dict[str, int] = field(default_factory=dict)
    
    # Warning stats
    total_warnings: int = 0
    warning_types: Dict[str, int] = field(default_factory=dict)
    
    # Property-safety stats
    property_correlations: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    # Binding-safety stats
    binding_correlations: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    # Activity-safety stats
    activity_correlations: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
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
            "toxicity": {
                "total_risks": self.total_risks,
                "high_risks": self.high_risks,
                "risk_types": self.risk_types,
                "risk_scores": self.risk_scores,
            },
            "patterns": {
                "total": self.total_patterns,
                "types": self.pattern_types,
            },
            "interactions": {
                "total": self.total_interactions,
                "types": self.interaction_types,
            },
            "warnings": {
                "total": self.total_warnings,
                "types": self.warning_types,
            },
            "correlations": {
                "property": self.property_correlations,
                "binding": self.binding_correlations,
                "activity": self.activity_correlations,
            },
            "errors": self.errors,
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get analysis success rate."""
        if not self.total_analyses:
            return None
        return self.successful_analyses / self.total_analyses


class SafetyAnalyzer:
    """Analyzer for safety data."""

    def __init__(
        self,
        cache_dir: Optional[Union[str, Path]] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        """Initialize safety analyzer.
        
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
        self.toxicity_threshold = self.config.get(
            "toxicity_threshold", 0.5
        )
        self.correlation_threshold = self.config.get(
            "correlation_threshold", 0.3
        )
        
        # Initialize stats
        self.stats = SafetyStats()

    def analyze(
        self,
        compound: CompoundData,
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
    ) -> CompoundData:
        """Analyze safety data.
        
        Args:
            compound: CompoundData instance to analyze
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            activity_data: Optional activity analysis data
            
        Returns:
            Analyzed CompoundData instance
        """
        try:
            self.stats.total_analyses += 1
            
            # Get safety predictions
            safety_data = compound.get_prediction("safety_ensemble")
            if not safety_data:
                return compound
            
            # Analyze toxicity risks
            risks = self._analyze_toxicity_risks(
                safety_data,
                property_data,
                binding_data,
                activity_data,
            )
            compound.toxicity_risks = risks
            
            # Update risk stats
            self._update_risk_stats(risks)
            
            # Analyze safety patterns
            patterns = self._analyze_safety_patterns(
                safety_data,
                property_data,
                binding_data,
                activity_data,
            )
            compound.safety_patterns = patterns
            
            # Update pattern stats
            self._update_pattern_stats(patterns)
            
            # Analyze interactions
            interactions = self._analyze_interactions(
                safety_data,
                property_data,
                binding_data,
                activity_data,
            )
            compound.safety_interactions = interactions
            
            # Update interaction stats
            self._update_interaction_stats(interactions)
            
            # Analyze warnings
            warnings = self._analyze_warnings(
                safety_data,
                property_data,
                binding_data,
                activity_data,
            )
            compound.safety_warnings = warnings
            
            # Update warning stats
            self._update_warning_stats(warnings)
            
            # Analyze property correlations
            if property_data:
                correlations = self._analyze_property_correlations(
                    safety_data,
                    property_data,
                )
                compound.property_safety_correlations = correlations
                
                # Update correlation stats
                self._update_property_correlation_stats(correlations)
            
            # Analyze binding correlations
            if binding_data:
                correlations = self._analyze_binding_correlations(
                    safety_data,
                    binding_data,
                )
                compound.binding_safety_correlations = correlations
                
                # Update correlation stats
                self._update_binding_correlation_stats(correlations)
            
            # Analyze activity correlations
            if activity_data:
                correlations = self._analyze_activity_correlations(
                    safety_data,
                    activity_data,
                )
                compound.activity_safety_correlations = correlations
                
                # Update correlation stats
                self._update_activity_correlation_stats(correlations)
            
            self.stats.successful_analyses += 1
            return compound
            
        except Exception as e:
            self.logger.error(
                f"Failed to analyze safety data for {compound.name}: {str(e)}"
            )
            self.stats.failed_analyses += 1
            self.stats.errors.append({
                "type": "safety_analysis_error",
                "compound": compound.name,
                "error": str(e),
                "timestamp": datetime.now().isoformat(),
            })
            raise

    def _analyze_toxicity_risks(
        self,
        safety_data: Dict[str, Any],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Analyze toxicity risks.
        
        Args:
            safety_data: Safety prediction data
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            activity_data: Optional activity analysis data
            
        Returns:
            List of toxicity risks
        """
        risks = []
        
        # Process each risk
        for risk in safety_data.get("risks", []):
            # Check confidence
            confidence = risk.get("confidence", 0)
            if confidence < self.confidence_threshold:
                continue
            
            # Get risk info
            risk_info = {
                "type": risk.get("type", "unknown"),
                "score": risk.get("score", 0),
                "confidence": confidence,
                "severity": risk.get("severity", "unknown"),
                "mechanism": risk.get("mechanism"),
                "organs": risk.get("organs", []),
                "warnings": risk.get("warnings", []),
            }
            
            # Enhance with property data
            if property_data:
                risk_info["property_effects"] = self._get_property_effects(
                    risk_info,
                    property_data,
                )
            
            # Enhance with binding data
            if binding_data:
                risk_info["binding_effects"] = self._get_binding_effects(
                    risk_info,
                    binding_data,
                )
            
            # Enhance with activity data
            if activity_data:
                risk_info["activity_effects"] = self._get_activity_effects(
                    risk_info,
                    activity_data,
                )
            
            # Add to risks
            risks.append(risk_info)
        
        return risks

    def _get_property_effects(
        self,
        risk_info: Dict[str, Any],
        property_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Get property-based effects on safety.
        
        Args:
            risk_info: Risk information
            property_data: Property analysis data
            
        Returns:
            List of property effects
        """
        effects = []
        
        # Check size effects
        if property_data.get("molecular_weight", 0) > 500:
            effects.append({
                "type": "size",
                "effect": "increased_toxicity",
                "confidence": 0.8,
            })
        
        # Check polarity effects
        logp = property_data.get("logp", 0)
        if not (-2 <= logp <= 5):
            effects.append({
                "type": "polarity",
                "effect": "altered_distribution",
                "confidence": 0.7,
            })
        
        # Check BBB effects
        if property_data.get("bbb", {}).get("pass"):
            effects.append({
                "type": "bbb",
                "effect": "cns_toxicity",
                "confidence": 0.9,
            })
        
        return effects

    def _get_binding_effects(
        self,
        risk_info: Dict[str, Any],
        binding_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Get binding-based effects on safety.
        
        Args:
            risk_info: Risk information
            binding_data: Binding analysis data
            
        Returns:
            List of binding effects
        """
        effects = []
        
        # Check target effects
        if binding_data.get("primary_targets"):
            effects.append({
                "type": "target",
                "effect": "off_target_toxicity",
                "confidence": 0.9,
                "targets": binding_data["primary_targets"],
            })
        
        # Check selectivity effects
        if binding_data.get("selectivity_ratios", {}):
            effects.append({
                "type": "selectivity",
                "effect": "non_selective_toxicity",
                "confidence": 0.8,
                "ratios": binding_data["selectivity_ratios"],
            })
        
        return effects

    def _get_activity_effects(
        self,
        risk_info: Dict[str, Any],
        activity_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Get activity-based effects on safety.
        
        Args:
            risk_info: Risk information
            activity_data: Activity analysis data
            
        Returns:
            List of activity effects
        """
        effects = []
        
        # Check mechanism effects
        if activity_data.get("mechanisms"):
            effects.append({
                "type": "mechanism",
                "effect": "mechanism_toxicity",
                "confidence": 0.9,
                "mechanisms": activity_data["mechanisms"],
            })
        
        # Check effect overlap
        if activity_data.get("effects"):
            effects.append({
                "type": "effects",
                "effect": "effect_toxicity",
                "confidence": 0.8,
                "effects": activity_data["effects"],
            })
        
        return effects

    def _analyze_safety_patterns(
        self,
        safety_data: Dict[str, Any],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Analyze safety patterns.
        
        Args:
            safety_data: Safety prediction data
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            activity_data: Optional activity analysis data
            
        Returns:
            List of safety patterns
        """
        patterns = []
        
        # Group risks by type
        types = {}
        for risk in safety_data.get("risks", []):
            risk_type = risk.get("type")
            if not risk_type:
                continue
            
            if risk_type not in types:
                types[risk_type] = []
            types[risk_type].append(risk)
        
        # Find patterns within types
        for risk_type, risks in types.items():
            if len(risks) < 2:
                continue
            
            # Sort by score
            risks = sorted(
                risks,
                key=lambda x: x.get("score", 0),
                reverse=True,
            )
            
            # Check for patterns
            patterns.extend(
                self._find_type_patterns(
                    risk_type,
                    risks,
                    property_data,
                    binding_data,
                    activity_data,
                )
            )
        
        return patterns

    def _find_type_patterns(
        self,
        risk_type: str,
        risks: List[Dict[str, Any]],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Find safety patterns within a risk type.
        
        Args:
            risk_type: Risk type name
            risks: List of risks of this type
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            activity_data: Optional activity analysis data
            
        Returns:
            List of safety patterns
        """
        patterns = []
        
        # Check score range pattern
        max_score = risks[0].get("score", 0)
        min_score = risks[-1].get("score", 0)
        if max_score - min_score >= 0.5:  # Significant score range
            patterns.append({
                "type": "score_range",
                "risk_type": risk_type,
                "min_score": min_score,
                "max_score": max_score,
                "range": max_score - min_score,
            })
        
        # Check severity pattern
        severities = {}
        for risk in risks:
            severity = risk.get("severity", "unknown")
            if severity not in severities:
                severities[severity] = 0
            severities[severity] += 1
        
        if len(severities) > 1:
            patterns.append({
                "type": "severities",
                "risk_type": risk_type,
                "severities": severities,
            })
        
        # Check organ pattern
        organs = {}
        for risk in risks:
            for organ in risk.get("organs", []):
                if organ not in organs:
                    organs[organ] = 0
                organs[organ] += 1
        
        if len(organs) > 1:
            patterns.append({
                "type": "organs",
                "risk_type": risk_type,
                "organs": organs,
            })
        
        # Check property-based patterns
        if property_data:
            property_patterns = self._find_property_patterns(
                risk_type,
                risks,
                property_data,
            )
            patterns.extend(property_patterns)
        
        # Check binding-based patterns
        if binding_data:
            binding_patterns = self._find_binding_patterns(
                risk_type,
                risks,
                binding_data,
            )
            patterns.extend(binding_patterns)
        
        # Check activity-based patterns
        if activity_data:
            activity_patterns = self._find_activity_patterns(
                risk_type,
                risks,
                activity_data,
            )
            patterns.extend(activity_patterns)
        
        return patterns

    def _find_property_patterns(
        self,
        risk_type: str,
        risks: List[Dict[str, Any]],
        property_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Find property-based safety patterns.
        
        Args:
            risk_type: Risk type name
            risks: List of risks
            property_data: Property analysis data
            
        Returns:
            List of property-based patterns
        """
        patterns = []
        
        # Check size-toxicity pattern
        if property_data.get("molecular_weight", 0) > 500:
            high_scores = [
                r for r in risks
                if r.get("score", 0) >= self.toxicity_threshold
            ]
            if high_scores:
                patterns.append({
                    "type": "size_toxicity",
                    "risk_type": risk_type,
                    "molecular_weight": property_data["molecular_weight"],
                    "high_risk_count": len(high_scores),
                })
        
        # Check polarity-toxicity pattern
        logp = property_data.get("logp", 0)
        if not (-2 <= logp <= 5):
            high_scores = [
                r for r in risks
                if r.get("score", 0) >= self.toxicity_threshold
            ]
            if high_scores:
                patterns.append({
                    "type": "polarity_toxicity",
                    "risk_type": risk_type,
                    "logp": logp,
                    "high_risk_count": len(high_scores),
                })
        
        return patterns

    def _find_binding_patterns(
        self,
        risk_type: str,
        risks: List[Dict[str, Any]],
        binding_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Find binding-based safety patterns.
        
        Args:
            risk_type: Risk type name
            risks: List of risks
            binding_data: Binding analysis data
            
        Returns:
            List of binding-based patterns
        """
        patterns = []
        
        # Check target-toxicity pattern
        if binding_data.get("primary_targets"):
            high_scores = [
                r for r in risks
                if r.get("score", 0) >= self.toxicity_threshold
            ]
            if high_scores:
                patterns.append({
                    "type": "target_toxicity",
                    "risk_type": risk_type,
                    "targets": binding_data["primary_targets"],
                    "high_risk_count": len(high_scores),
                })
        
        # Check selectivity-toxicity pattern
        if binding_data.get("selectivity_ratios", {}):
            high_scores = [
                r for r in risks
                if r.get("score", 0) >= self.toxicity_threshold
            ]
            if high_scores:
                patterns.append({
                    "type": "selectivity_toxicity",
                    "risk_type": risk_type,
                    "selectivity": binding_data["selectivity_ratios"],
                    "high_risk_count": len(high_scores),
                })
        
        return patterns

    def _find_activity_patterns(
        self,
        risk_type: str,
        risks: List[Dict[str, Any]],
        activity_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Find activity-based safety patterns.
        
        Args:
            risk_type: Risk type name
            risks: List of risks
            activity_data: Activity analysis data
            
        Returns:
            List of activity-based patterns
        """
        patterns = []
        
        # Check mechanism-toxicity pattern
        if activity_data.get("mechanisms"):
            high_scores = [
                r for r in risks
                if r.get("score", 0) >= self.toxicity_threshold
            ]
            if high_scores:
                patterns.append({
                    "type": "mechanism_toxicity",
                    "risk_type": risk_type,
                    "mechanisms": activity_data["mechanisms"],
                    "high_risk_count": len(high_scores),
                })
        
        # Check effect-toxicity pattern
        if activity_data.get("effects"):
            high_scores = [
                r for r in risks
                if r.get("score", 0) >= self.toxicity_threshold
            ]
            if high_scores:
                patterns.append({
                    "type": "effect_toxicity",
                    "risk_type": risk_type,
                    "effects": activity_data["effects"],
                    "high_risk_count": len(high_scores),
                })
        
        return patterns

    def _analyze_interactions(
        self,
        safety_data: Dict[str, Any],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Analyze drug interactions.
        
        Args:
            safety_data: Safety prediction data
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            activity_data: Optional activity analysis data
            
        Returns:
            List of drug interactions
        """
        interactions = []
        
        # Process each interaction
        for interaction in safety_data.get("interactions", []):
            # Check confidence
            confidence = interaction.get("confidence", 0)
            if confidence < self.confidence_threshold:
                continue
            
            # Get interaction info
            interaction_info = {
                "type": interaction.get("type", "unknown"),
                "severity": interaction.get("severity", "unknown"),
                "confidence": confidence,
                "mechanism": interaction.get("mechanism"),
                "effects": interaction.get("effects", []),
                "warnings": interaction.get("warnings", []),
            }
            
            # Enhance with property data
            if property_data:
                interaction_info["property_effects"] = self._get_property_effects(
                    interaction_info,
                    property_data,
                )
            
            # Enhance with binding data
            if binding_data:
                interaction_info["binding_effects"] = self._get_binding_effects(
                    interaction_info,
                    binding_data,
                )
            
            # Enhance with activity data
            if activity_data:
                interaction_info["activity_effects"] = self._get_activity_effects(
                    interaction_info,
                    activity_data,
                )
            
            # Add to interactions
            interactions.append(interaction_info)
        
        return interactions

    def _analyze_warnings(
        self,
        safety_data: Dict[str, Any],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Analyze safety warnings.
        
        Args:
            safety_data: Safety prediction data
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            activity_data: Optional activity analysis data
            
        Returns:
            List of safety warnings
        """
        warnings = []
        
        # Process each warning
        for warning in safety_data.get("warnings", []):
            # Check confidence
            confidence = warning.get("confidence", 0)
            if confidence < self.confidence_threshold:
                continue
            
            # Get warning info
            warning_info = {
                "type": warning.get("type", "unknown"),
                "severity": warning.get("severity", "unknown"),
                "confidence": confidence,
                "description": warning.get("description"),
                "recommendations": warning.get("recommendations", []),
            }
            
            # Enhance with property data
            if property_data:
                warning_info["property_effects"] = self._get_property_effects(
                    warning_info,
                    property_data,
                )
            
            # Enhance with binding data
            if binding_data:
                warning_info["binding_effects"] = self._get_binding_effects(
                    warning_info,
                    binding_data,
                )
            
            # Enhance with activity data
            if activity_data:
                warning_info["activity_effects"] = self._get_activity_effects(
                    warning_info,
                    activity_data,
                )
            
            # Add to warnings
            warnings.append(warning_info)
        
        return warnings

    def _analyze_property_correlations(
        self,
        safety_data: Dict[str, Any],
        property_data: Dict[str, Any],
    ) -> Dict[str, Dict[str, float]]:
        """Analyze property-safety correlations.
        
        Args:
            safety_data: Safety prediction data
            property_data: Property analysis data
            
        Returns:
            Dictionary of property-safety correlations
        """
        correlations = {}
        
        # Get risk scores by type
        risk_scores = {}
        for risk in safety_data.get("risks", []):
            risk_type = risk.get("type", "unknown")
            if risk_type not in risk_scores:
                risk_scores[risk_type] = []
            risk_scores[risk_type].append(
                risk.get("score", 0)
            )
        
        # Calculate correlations
        for prop_name, prop_value in property_data.items():
            if not isinstance(prop_value, (int, float)):
                continue
            
            correlations[prop_name] = {}
            for risk_type, scores in risk_scores.items():
                correlation = self._calculate_correlation(
                    [prop_value] * len(scores),
                    scores,
                )
                if abs(correlation) >= self.correlation_threshold:
                    correlations[prop_name][risk_type] = correlation
        
        return correlations

    def _analyze_binding_correlations(
        self,
        safety_data: Dict[str, Any],
        binding_data: Dict[str, Any],
    ) -> Dict[str, Dict[str, float]]:
        """Analyze binding-safety correlations.
        
        Args:
            safety_data: Safety prediction data
            binding_data: Binding analysis data
            
        Returns:
            Dictionary of binding-safety correlations
        """
        correlations = {}
        
        # Get risk scores by type
        risk_scores = {}
        for risk in safety_data.get("risks", []):
            risk_type = risk.get("type", "unknown")
            if risk_type not in risk_scores:
                risk_scores[risk_type] = []
            risk_scores[risk_type].append(
                risk.get("score", 0)
            )
        
        # Calculate correlations for each target
        for target in binding_data.get("targets", []):
            target_name = target["name"]
            target_affinity = target.get("affinity", float("inf"))
            
            correlations[target_name] = {}
            for risk_type, scores in risk_scores.items():
                correlation = self._calculate_correlation(
                    [target_affinity] * len(scores),
                    scores,
                )
                if abs(correlation) >= self.correlation_threshold:
                    correlations[target_name][risk_type] = correlation
        
        return correlations

    def _analyze_activity_correlations(
        self,
        safety_data: Dict[str, Any],
        activity_data: Dict[str, Any],
    ) -> Dict[str, Dict[str, float]]:
        """Analyze activity-safety correlations.
        
        Args:
            safety_data: Safety prediction data
            activity_data: Activity analysis data
            
        Returns:
            Dictionary of activity-safety correlations
        """
        correlations = {}
        
        # Get risk scores by type
        risk_scores = {}
        for risk in safety_data.get("risks", []):
            risk_type = risk.get("type", "unknown")
            if risk_type not in risk_scores:
                risk_scores[risk_type] = []
            risk_scores[risk_type].append(
                risk.get("score", 0)
            )
        
        # Calculate correlations for each activity
        for activity in activity_data.get("activities", []):
            activity_name = activity["name"]
            activity_score = activity.get("score", 0)
            
            correlations[activity_name] = {}
            for risk_type, scores in risk_scores.items():
                correlation = self._calculate_correlation(
                    [activity_score] * len(scores),
                    scores,
                )
                if abs(correlation) >= self.correlation_threshold:
                    correlations[activity_name][risk_type] = correlation
        
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

    def _update_risk_stats(
        self,
        risks: List[Dict[str, Any]],
    ) -> None:
        """Update risk statistics."""
        self.stats.total_risks += len(risks)
        
        # Update high risk count
        self.stats.high_risks += sum(
            1 for risk in risks
            if risk["score"] >= self.toxicity_threshold
        )
        
        # Update type stats
        for risk in risks:
            risk_type = risk["type"]
            if risk_type not in self.stats.risk_types:
                self.stats.risk_types[risk_type] = 0
            self.stats.risk_types[risk_type] += 1
            
            # Update score stats
            if risk_type not in self.stats.risk_scores:
                self.stats.risk_scores[risk_type] = {
                    "min": float("inf"),
                    "max": float("-inf"),
                    "sum": 0,
                    "count": 0,
                }
            scores = self.stats.risk_scores[risk_type]
            score = risk["score"]
            scores["min"] = min(scores["min"], score)
            scores["max"] = max(scores["max"], score)
            scores["sum"] += score
            scores["count"] += 1

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

    def _update_interaction_stats(
        self,
        interactions: List[Dict[str, Any]],
    ) -> None:
        """Update interaction statistics."""
        self.stats.total_interactions += len(interactions)
        
        for interaction in interactions:
            interaction_type = interaction["type"]
            if interaction_type not in self.stats.interaction_types:
                self.stats.interaction_types[interaction_type] = 0
            self.stats.interaction_types[interaction_type] += 1

    def _update_warning_stats(
        self,
        warnings: List[Dict[str, Any]],
    ) -> None:
        """Update warning statistics."""
        self.stats.total_warnings += len(warnings)
        
        for warning in warnings:
            warning_type = warning["type"]
            if warning_type not in self.stats.warning_types:
                self.stats.warning_types[warning_type] = 0
            self.stats.warning_types[warning_type] += 1

    def get_info(self) -> Dict[str, Any]:
        """Get analyzer information."""
        return {
            "config": {
                "confidence_threshold": self.confidence_threshold,
                "toxicity_threshold": self.toxicity_threshold,
            },
            "stats": self.stats.to_dict(),
        }
