"""SAR analysis functionality.

This module provides the SARAnalyzer class that handles:
1. Pharmacophore analysis
2. Activity cliff detection
3. Structure similarity analysis
4. Property relationship analysis
5. Integration with property, binding, activity, and safety analysis
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field
from datetime import datetime

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from rdkit.Chem.MolStandardize import rdMolStandardize

from ...models import CompoundData
from ..processors.structure import (
    StructureValidator,
    PharmacophoreDetector,
    SimilaritySearcher,
)


@dataclass
class SARStats:
    """SAR analysis statistics."""
    
    # Analysis counts
    total_analyses: int = 0
    successful_analyses: int = 0
    failed_analyses: int = 0
    
    # Pharmacophore stats
    total_pharmacophores: int = 0
    pharmacophore_types: Dict[str, int] = field(default_factory=dict)
    pharmacophore_scores: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    # Activity cliff stats
    total_cliffs: int = 0
    cliff_types: Dict[str, int] = field(default_factory=dict)
    cliff_scores: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    # Similarity stats
    total_similarities: int = 0
    similarity_ranges: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    # Property stats
    total_relationships: int = 0
    relationship_types: Dict[str, int] = field(default_factory=dict)
    
    # Property-SAR stats
    property_correlations: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    # Binding-SAR stats
    binding_correlations: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    # Activity-SAR stats
    activity_correlations: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
    # Safety-SAR stats
    safety_correlations: Dict[str, Dict[str, float]] = field(default_factory=dict)
    
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
            "pharmacophores": {
                "total": self.total_pharmacophores,
                "types": self.pharmacophore_types,
                "scores": self.pharmacophore_scores,
            },
            "cliffs": {
                "total": self.total_cliffs,
                "types": self.cliff_types,
                "scores": self.cliff_scores,
            },
            "similarities": {
                "total": self.total_similarities,
                "ranges": self.similarity_ranges,
            },
            "relationships": {
                "total": self.total_relationships,
                "types": self.relationship_types,
            },
            "correlations": {
                "property": self.property_correlations,
                "binding": self.binding_correlations,
                "activity": self.activity_correlations,
                "safety": self.safety_correlations,
            },
            "errors": self.errors,
        }
    
    def _get_success_rate(self) -> Optional[float]:
        """Get analysis success rate."""
        if not self.total_analyses:
            return None
        return self.successful_analyses / self.total_analyses


class SARAnalyzer:
    """Analyzer for SAR data."""

    def __init__(
        self,
        cache_dir: Optional[Union[str, Path]] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        """Initialize SAR analyzer.
        
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
        self.similarity_threshold = self.config.get(
            "similarity_threshold", 0.7
        )
        self.pharmacophore_confidence = self.config.get(
            "pharmacophore_confidence", 0.7
        )
        self.correlation_threshold = self.config.get(
            "correlation_threshold", 0.3
        )
        
        # Initialize analyzers
        self._init_analyzers()
        
        # Initialize stats
        self.stats = SARStats()

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

    def analyze(
        self,
        compound: CompoundData,
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
        safety_data: Optional[Dict[str, Any]] = None,
    ) -> CompoundData:
        """Analyze SAR data.
        
        Args:
            compound: CompoundData instance to analyze
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            activity_data: Optional activity analysis data
            safety_data: Optional safety analysis data
            
        Returns:
            Analyzed CompoundData instance
        """
        try:
            self.stats.total_analyses += 1
            
            # Get structure
            mol = Chem.MolFromSmiles(compound.smiles)
            if not mol:
                return compound
            
            # Standardize structure
            mol = self.standardizer.standardize(mol)
            
            # Analyze pharmacophores
            pharmacophores = self._analyze_pharmacophores(
                mol,
                property_data,
                binding_data,
                activity_data,
                safety_data,
            )
            compound.pharmacophores = pharmacophores
            
            # Update pharmacophore stats
            self._update_pharmacophore_stats(pharmacophores)
            
            # Analyze activity cliffs
            cliffs = self._analyze_activity_cliffs(
                compound,
                property_data,
                binding_data,
                activity_data,
                safety_data,
            )
            compound.activity_cliffs = cliffs
            
            # Update cliff stats
            self._update_cliff_stats(cliffs)
            
            # Analyze similarities
            similarities = self._analyze_similarities(
                mol,
                property_data,
                binding_data,
                activity_data,
                safety_data,
            )
            compound.structure_similarities = similarities
            
            # Update similarity stats
            self._update_similarity_stats(similarities)
            
            # Analyze property relationships
            relationships = self._analyze_relationships(
                compound,
                property_data,
                binding_data,
                activity_data,
                safety_data,
            )
            compound.property_relationships = relationships
            
            # Update relationship stats
            self._update_relationship_stats(relationships)
            
            # Analyze property correlations
            if property_data:
                correlations = self._analyze_property_correlations(
                    compound,
                    property_data,
                )
                compound.property_sar_correlations = correlations
                
                # Update correlation stats
                self._update_property_correlation_stats(correlations)
            
            # Analyze binding correlations
            if binding_data:
                correlations = self._analyze_binding_correlations(
                    compound,
                    binding_data,
                )
                compound.binding_sar_correlations = correlations
                
                # Update correlation stats
                self._update_binding_correlation_stats(correlations)
            
            # Analyze activity correlations
            if activity_data:
                correlations = self._analyze_activity_correlations(
                    compound,
                    activity_data,
                )
                compound.activity_sar_correlations = correlations
                
                # Update correlation stats
                self._update_activity_correlation_stats(correlations)
            
            # Analyze safety correlations
            if safety_data:
                correlations = self._analyze_safety_correlations(
                    compound,
                    safety_data,
                )
                compound.safety_sar_correlations = correlations
                
                # Update correlation stats
                self._update_safety_correlation_stats(correlations)
            
            self.stats.successful_analyses += 1
            return compound
            
        except Exception as e:
            self.logger.error(
                f"Failed to analyze SAR data for {compound.name}: {str(e)}"
            )
            self.stats.failed_analyses += 1
            self.stats.errors.append({
                "type": "sar_analysis_error",
                "compound": compound.name,
                "error": str(e),
                "timestamp": datetime.now().isoformat(),
            })
            raise

    def _analyze_pharmacophores(
        self,
        mol: Chem.Mol,
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
        safety_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Analyze pharmacophores.
        
        Args:
            mol: RDKit molecule
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            activity_data: Optional activity analysis data
            safety_data: Optional safety analysis data
            
        Returns:
            List of pharmacophores
        """
        pharmacophores = []
        
        # Detect pharmacophores
        detected = self.pharmacophore_detector.detect(mol)
        
        # Process each pharmacophore
        for pharm in detected:
            # Check confidence
            confidence = pharm.get("confidence", 0)
            if confidence < self.pharmacophore_confidence:
                continue
            
            # Get pharmacophore info
            pharm_info = {
                "type": pharm.get("type", "unknown"),
                "score": pharm.get("score", 0),
                "confidence": confidence,
                "atoms": pharm.get("atoms", []),
                "features": pharm.get("features", []),
            }
            
            # Enhance with property data
            if property_data:
                pharm_info["property_effects"] = self._get_property_effects(
                    pharm_info,
                    property_data,
                )
            
            # Enhance with binding data
            if binding_data:
                pharm_info["binding_effects"] = self._get_binding_effects(
                    pharm_info,
                    binding_data,
                )
            
            # Enhance with activity data
            if activity_data:
                pharm_info["activity_effects"] = self._get_activity_effects(
                    pharm_info,
                    activity_data,
                )
            
            # Enhance with safety data
            if safety_data:
                pharm_info["safety_effects"] = self._get_safety_effects(
                    pharm_info,
                    safety_data,
                )
            
            # Add to pharmacophores
            pharmacophores.append(pharm_info)
        
        return pharmacophores

    def _get_property_effects(
        self,
        pharm_info: Dict[str, Any],
        property_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Get property-based effects on pharmacophores.
        
        Args:
            pharm_info: Pharmacophore information
            property_data: Property analysis data
            
        Returns:
            List of property effects
        """
        effects = []
        
        # Check size effects
        if property_data.get("molecular_weight", 0) > 500:
            effects.append({
                "type": "size",
                "effect": "reduced_binding",
                "confidence": 0.8,
            })
        
        # Check polarity effects
        logp = property_data.get("logp", 0)
        if not (-2 <= logp <= 5):
            effects.append({
                "type": "polarity",
                "effect": "altered_binding",
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
        pharm_info: Dict[str, Any],
        binding_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Get binding-based effects on pharmacophores.
        
        Args:
            pharm_info: Pharmacophore information
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

    def _get_activity_effects(
        self,
        pharm_info: Dict[str, Any],
        activity_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Get activity-based effects on pharmacophores.
        
        Args:
            pharm_info: Pharmacophore information
            activity_data: Activity analysis data
            
        Returns:
            List of activity effects
        """
        effects = []
        
        # Check mechanism effects
        if activity_data.get("mechanisms"):
            effects.append({
                "type": "mechanism",
                "effect": "mechanism_activity",
                "confidence": 0.9,
                "mechanisms": activity_data["mechanisms"],
            })
        
        # Check effect overlap
        if activity_data.get("effects"):
            effects.append({
                "type": "effects",
                "effect": "effect_activity",
                "confidence": 0.8,
                "effects": activity_data["effects"],
            })
        
        return effects

    def _get_safety_effects(
        self,
        pharm_info: Dict[str, Any],
        safety_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Get safety-based effects on pharmacophores.
        
        Args:
            pharm_info: Pharmacophore information
            safety_data: Safety analysis data
            
        Returns:
            List of safety effects
        """
        effects = []
        
        # Check toxicity effects
        if safety_data.get("risks"):
            effects.append({
                "type": "toxicity",
                "effect": "toxicity_risk",
                "confidence": 0.9,
                "risks": safety_data["risks"],
            })
        
        # Check interaction effects
        if safety_data.get("interactions"):
            effects.append({
                "type": "interactions",
                "effect": "interaction_risk",
                "confidence": 0.8,
                "interactions": safety_data["interactions"],
            })
        
        return effects

    def _analyze_activity_cliffs(
        self,
        compound: CompoundData,
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
        safety_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Analyze activity cliffs.
        
        Args:
            compound: CompoundData instance
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            activity_data: Optional activity analysis data
            safety_data: Optional safety analysis data
            
        Returns:
            List of activity cliffs
        """
        cliffs = []
        
        # Get activity data
        activity_data = activity_data or compound.get_prediction("activity_ensemble")
        if not activity_data:
            return cliffs
        
        # Get similar compounds
        mol = Chem.MolFromSmiles(compound.smiles)
        if not mol:
            return cliffs
        
        similar = self.similarity_searcher.search(
            mol,
            threshold=self.similarity_threshold,
        )
        
        # Compare activities
        for sim in similar:
            sim_activity = sim.get("activity")
            if not sim_activity:
                continue
            
            # Calculate activity difference
            diff = abs(
                activity_data.get("score", 0) -
                sim_activity.get("score", 0)
            )
            
            if diff >= 0.5:  # Significant activity difference
                cliff_info = {
                    "compound": sim["name"],
                    "similarity": sim["similarity"],
                    "activity_diff": diff,
                    "type": "activity_cliff",
                }
                
                # Enhance with property data
                if property_data:
                    cliff_info["property_effects"] = self._get_property_effects(
                        cliff_info,
                        property_data,
                    )
                
                # Enhance with binding data
                if binding_data:
                    cliff_info["binding_effects"] = self._get_binding_effects(
                        cliff_info,
                        binding_data,
                    )
                
                # Enhance with activity data
                if activity_data:
                    cliff_info["activity_effects"] = self._get_activity_effects(
                        cliff_info,
                        activity_data,
                    )
                
                # Enhance with safety data
                if safety_data:
                    cliff_info["safety_effects"] = self._get_safety_effects(
                        cliff_info,
                        safety_data,
                    )
                
                cliffs.append(cliff_info)
        
        return cliffs

    def _analyze_similarities(
        self,
        mol: Chem.Mol,
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
        safety_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Analyze structure similarities.
        
        Args:
            mol: RDKit molecule
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            activity_data: Optional activity analysis data
            safety_data: Optional safety analysis data
            
        Returns:
            List of structure similarities
        """
        similarities = []
        
        # Search similar compounds
        similar = self.similarity_searcher.search(
            mol,
            threshold=self.similarity_threshold,
        )
        
        # Process each similar compound
        for sim in similar:
            # Get similarity info
            sim_info = {
                "compound": sim["name"],
                "similarity": sim["similarity"],
                "type": sim.get("type", "unknown"),
                "features": sim.get("features", []),
            }
            
            # Enhance with property data
            if property_data:
                sim_info["property_effects"] = self._get_property_effects(
                    sim_info,
                    property_data,
                )
            
            # Enhance with binding data
            if binding_data:
                sim_info["binding_effects"] = self._get_binding_effects(
                    sim_info,
                    binding_data,
                )
            
            # Enhance with activity data
            if activity_data:
                sim_info["activity_effects"] = self._get_activity_effects(
                    sim_info,
                    activity_data,
                )
            
            # Enhance with safety data
            if safety_data:
                sim_info["safety_effects"] = self._get_safety_effects(
                    sim_info,
                    safety_data,
                )
            
            # Add to similarities
            similarities.append(sim_info)
        
        return similarities

    def _analyze_relationships(
        self,
        compound: CompoundData,
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
        safety_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Analyze property relationships.
        
        Args:
            compound: CompoundData instance
            property_data: Optional property analysis data
            binding_data: Optional binding analysis data
            activity_data: Optional activity analysis data
            safety_data: Optional safety analysis data
            
        Returns:
            List of property relationships
        """
        relationships = []
        
        # Get structure
        mol = Chem.MolFromSmiles(compound.smiles)
        if not mol:
            return relationships
        
        # Calculate basic properties
        properties = self._calculate_basic_properties(mol)
        
        # Analyze different relationship types
        self._analyze_lipinski_relationship(
            properties, relationships, property_data, binding_data, 
            activity_data, safety_data
        )
        self._analyze_bbb_relationship(
            properties, relationships, property_data, binding_data,
            activity_data, safety_data
        )
        self._analyze_oral_relationship(
            properties, relationships, property_data, binding_data,
            activity_data, safety_data
        )
        
        return relationships

    def _calculate_basic_properties(
        self,
        mol: Chem.Mol,
    ) -> Dict[str, float]:
        """Calculate basic molecular properties.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Dictionary of property values
        """
        return {
            "mw": Descriptors.ExactMolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "tpsa": Descriptors.TPSA(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "rotatable": Descriptors.NumRotatableBonds(mol),
        }

    def _analyze_lipinski_relationship(
        self,
        properties: Dict[str, float],
        relationships: List[Dict[str, Any]],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
        safety_data: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Analyze Lipinski's Rule of Five relationship."""
        mw = properties["mw"]
        logp = properties["logp"]
        hbd = properties["hbd"]
        hba = properties["hba"]
        
        if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
            rel_info = {
                "type": "lipinski",
                "score": 1.0,
                "properties": {
                    "mw": mw,
                    "logp": logp,
                    "hbd": hbd,
                    "hba": hba,
                },
            }
            self._enhance_relationship_info(
                rel_info, property_data, binding_data,
                activity_data, safety_data
            )
            relationships.append(rel_info)

    def _analyze_bbb_relationship(
        self,
        properties: Dict[str, float],
        relationships: List[Dict[str, Any]],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
        safety_data: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Analyze blood-brain barrier relationship."""
        logp = properties["logp"]
        mw = properties["mw"]
        tpsa = properties["tpsa"]
        
        if logp >= 0 and logp <= 6 and mw <= 400 and tpsa <= 90:
            rel_info = {
                "type": "bbb",
                "score": 1.0,
                "properties": {
                    "logp": logp,
                    "mw": mw,
                    "tpsa": tpsa,
                },
            }
            self._enhance_relationship_info(
                rel_info, property_data, binding_data,
                activity_data, safety_data
            )
            relationships.append(rel_info)

    def _analyze_oral_relationship(
        self,
        properties: Dict[str, float],
        relationships: List[Dict[str, Any]],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
        safety_data: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Analyze oral bioavailability relationship."""
        mw = properties["mw"]
        logp = properties["logp"]
        rotatable = properties["rotatable"]
        tpsa = properties["tpsa"]
        hbd = properties["hbd"]
        hba = properties["hba"]
        
        if (mw <= 500 and logp <= 5 and rotatable <= 10
                and tpsa <= 140 and hbd <= 5 and hba <= 10):
            rel_info = {
                "type": "oral",
                "score": 1.0,
                "properties": {
                    "mw": mw,
                    "logp": logp,
                    "rotatable": rotatable,
                    "tpsa": tpsa,
                    "hbd": hbd,
                    "hba": hba,
                },
            }
            self._enhance_relationship_info(
                rel_info, property_data, binding_data,
                activity_data, safety_data
            )
            relationships.append(rel_info)

    def _enhance_relationship_info(
        self,
        rel_info: Dict[str, Any],
        property_data: Optional[Dict[str, Any]] = None,
        binding_data: Optional[Dict[str, Any]] = None,
        activity_data: Optional[Dict[str, Any]] = None,
        safety_data: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Enhance relationship info with additional data."""
        if property_data:
            rel_info["property_effects"] = self._get_property_effects(
                rel_info,
                property_data,
            )
        
        if binding_data:
            rel_info["binding_effects"] = self._get_binding_effects(
                rel_info,
                binding_data,
            )
        
        if activity_data:
            rel_info["activity_effects"] = self._get_activity_effects(
                rel_info,
                activity_data,
            )
        
        if safety_data:
            rel_info["safety_effects"] = self._get_safety_effects(
                rel_info,
                safety_data,
            )


    def _analyze_property_correlations(
        self,
        compound: CompoundData,
        property_data: Dict[str, Any],
    ) -> Dict[str, Dict[str, float]]:
        """Analyze property-SAR correlations.
        
        Args:
            compound: CompoundData instance
            property_data: Property analysis data
            
        Returns:
            Dictionary of property-SAR correlations
        """
        correlations = {}
        
        # Get SAR scores by type
        sar_scores = {}
        for pharm in compound.pharmacophores:
            pharm_type = pharm.get("type", "unknown")
            if pharm_type not in sar_scores:
                sar_scores[pharm_type] = []
            sar_scores[pharm_type].append(
                pharm.get("score", 0)
            )
        
        # Calculate correlations
        for prop_name, prop_value in property_data.items():
            if not isinstance(prop_value, (int, float)):
                continue
            
            correlations[prop_name] = {}
            for sar_type, scores in sar_scores.items():
                correlation = self._calculate_correlation(
                    [prop_value] * len(scores),
                    scores,
                )
                if abs(correlation) >= self.correlation_threshold:
                    correlations[prop_name][sar_type] = correlation
        
        return correlations

    def _analyze_binding_correlations(
        self,
        compound: CompoundData,
        binding_data: Dict[str, Any],
    ) -> Dict[str, Dict[str, float]]:
        """Analyze binding-SAR correlations.
        
        Args:
            compound: CompoundData instance
            binding_data: Binding analysis data
            
        Returns:
            Dictionary of binding-SAR correlations
        """
        correlations = {}
        
        # Get SAR scores by type
        sar_scores = {}
        for pharm in compound.pharmacophores:
            pharm_type = pharm.get("type", "unknown")
            if pharm_type not in sar_scores:
                sar_scores[pharm_type] = []
            sar_scores[pharm_type].append(
                pharm.get("score", 0)
            )
        
        # Calculate correlations for each target
        for target in binding_data.get("targets", []):
            target_name = target["name"]
            target_affinity = target.get("affinity", float("inf"))
            
            correlations[target_name] = {}
            for sar_type, scores in sar_scores.items():
                correlation = self._calculate_correlation(
                    [target_affinity] * len(scores),
                    scores,
                )
                if abs(correlation) >= self.correlation_threshold:
                    correlations[target_name][sar_type] = correlation
        
        return correlations

    def _analyze_activity_correlations(
        self,
        compound: CompoundData,
        activity_data: Dict[str, Any],
    ) -> Dict[str, Dict[str, float]]:
        """Analyze activity-SAR correlations.
        
        Args:
            compound: CompoundData instance
            activity_data: Activity analysis data
            
        Returns:
            Dictionary of activity-SAR correlations
        """
        correlations = {}
        
        # Get SAR scores by type
        sar_scores = {}
        for pharm in compound.pharmacophores:
            pharm_type = pharm.get("type", "unknown")
            if pharm_type not in sar_scores:
                sar_scores[pharm_type] = []
            sar_scores[pharm_type].append(
                pharm.get("score", 0)
            )
        
        # Calculate correlations for each activity
        for activity in activity_data.get("activities", []):
            activity_name = activity["name"]
            activity_score = activity.get("score", 0)
            
            correlations[activity_name] = {}
            for sar_type, scores in sar_scores.items():
                correlation = self._calculate_correlation(
                    [activity_score] * len(scores),
                    scores,
                )
                if abs(correlation) >= self.correlation_threshold:
                    correlations[activity_name][sar_type] = correlation
        
        return correlations

    def _analyze_safety_correlations(
        self,
        compound: CompoundData,
        safety_data: Dict[str, Any],
    ) -> Dict[str, Dict[str, float]]:
        """Analyze safety-SAR correlations.
        
        Args:
            compound: CompoundData instance
            safety_data: Safety analysis data
            
        Returns:
            Dictionary of safety-SAR correlations
        """
        correlations = {}
        
        # Get SAR scores by type
        sar_scores = {}
        for pharm in compound.pharmacophores:
            pharm_type = pharm.get("type", "unknown")
            if pharm_type not in sar_scores:
                sar_scores[pharm_type] = []
            sar_scores[pharm_type].append(
                pharm.get("score", 0)
            )
        
        # Calculate correlations for each risk
        for risk in safety_data.get("risks", []):
            risk_type = risk.get("type", "unknown")
            risk_score = risk.get("score", 0)
            
            correlations[risk_type] = {}
            for sar_type, scores in sar_scores.items():
                correlation = self._calculate_correlation(
                    [risk_score] * len(scores),
                    scores,
                )
                if abs(correlation) >= self.correlation_threshold:
                    correlations[risk_type][sar_type] = correlation
        
        return correlations

    def _update_property_correlation_stats(
        self,
        correlations: Dict[str, Dict[str, float]],
    ) -> None:
        """Update property correlation statistics."""
        for prop_name, prop_corrs in correlations.items():
            if prop_name not in self.stats.property_correlations:
                self.stats.property_correlations[prop_name] = {}
            
            for sar_type, correlation in prop_corrs.items():
                if sar_type not in self.stats.property_correlations[prop_name]:
                    self.stats.property_correlations[prop_name][sar_type] = correlation

    def _update_binding_correlation_stats(
        self,
        correlations: Dict[str, Dict[str, float]],
    ) -> None:
        """Update binding correlation statistics."""
        for target_name, target_corrs in correlations.items():
            if target_name not in self.stats.binding_correlations:
                self.stats.binding_correlations[target_name] = {}
            
            for sar_type, correlation in target_corrs.items():
                if sar_type not in self.stats.binding_correlations[target_name]:
                    self.stats.binding_correlations[target_name][sar_type] = correlation

    def _update_activity_correlation_stats(
        self,
        correlations: Dict[str, Dict[str, float]],
    ) -> None:
        """Update activity correlation statistics."""
        for activity_name, activity_corrs in correlations.items():
            if activity_name not in self.stats.activity_correlations:
                self.stats.activity_correlations[activity_name] = {}
            
            for sar_type, correlation in activity_corrs.items():
                if sar_type not in self.stats.activity_correlations[activity_name]:
                    self.stats.activity_correlations[activity_name][sar_type] = correlation

    def _update_safety_correlation_stats(
        self,
        correlations: Dict[str, Dict[str, float]],
    ) -> None:
        """Update safety correlation statistics."""
        for risk_type, risk_corrs in correlations.items():
            if risk_type not in self.stats.safety_correlations:
                self.stats.safety_correlations[risk_type] = {}
            
            for sar_type, correlation in risk_corrs.items():
                if sar_type not in self.stats.safety_correlations[risk_type]:
                    self.stats.safety_correlations[risk_type][sar_type] = correlation


    def _update_pharmacophore_stats(
        self,
        pharmacophores: List[Dict[str, Any]],
    ) -> None:
        """Update pharmacophore statistics."""
        self.stats.total_pharmacophores += len(pharmacophores)
        
        for pharm in pharmacophores:
            # Update type stats
            pharm_type = pharm["type"]
            if pharm_type not in self.stats.pharmacophore_types:
                self.stats.pharmacophore_types[pharm_type] = 0
            self.stats.pharmacophore_types[pharm_type] += 1
            
            # Update score stats
            if pharm_type not in self.stats.pharmacophore_scores:
                self.stats.pharmacophore_scores[pharm_type] = {
                    "min": float("inf"),
                    "max": float("-inf"),
                    "sum": 0,
                    "count": 0,
                }
            scores = self.stats.pharmacophore_scores[pharm_type]
            score = pharm["score"]
            scores["min"] = min(scores["min"], score)
            scores["max"] = max(scores["max"], score)
            scores["sum"] += score
            scores["count"] += 1

    def _update_cliff_stats(
        self,
        cliffs: List[Dict[str, Any]],
    ) -> None:
        """Update cliff statistics."""
        self.stats.total_cliffs += len(cliffs)
        
        for cliff in cliffs:
            # Update type stats
            cliff_type = cliff["type"]
            if cliff_type not in self.stats.cliff_types:
                self.stats.cliff_types[cliff_type] = 0
            self.stats.cliff_types[cliff_type] += 1
            
            # Update score stats
            if cliff_type not in self.stats.cliff_scores:
                self.stats.cliff_scores[cliff_type] = {
                    "min": float("inf"),
                    "max": float("-inf"),
                    "sum": 0,
                    "count": 0,
                }
            scores = self.stats.cliff_scores[cliff_type]
            score = cliff["activity_diff"]
            scores["min"] = min(scores["min"], score)
            scores["max"] = max(scores["max"], score)
            scores["sum"] += score
            scores["count"] += 1

    def _update_similarity_stats(
        self,
        similarities: List[Dict[str, Any]],
    ) -> None:
        """Update similarity statistics."""
        self.stats.total_similarities += len(similarities)
        
        for sim in similarities:
            # Update range stats
            sim_type = sim["type"]
            if sim_type not in self.stats.similarity_ranges:
                self.stats.similarity_ranges[sim_type] = {
                    "min": float("inf"),
                    "max": float("-inf"),
                    "sum": 0,
                    "count": 0,
                }
            ranges = self.stats.similarity_ranges[sim_type]
            similarity = sim["similarity"]
            ranges["min"] = min(ranges["min"], similarity)
            ranges["max"] = max(ranges["max"], similarity)
            ranges["sum"] += similarity
            ranges["count"] += 1

    def _update_relationship_stats(
        self,
        relationships: List[Dict[str, Any]],
    ) -> None:
        """Update relationship statistics."""
        self.stats.total_relationships += len(relationships)
        
        for rel in relationships:
            rel_type = rel["type"]
            if rel_type not in self.stats.relationship_types:
                self.stats.relationship_types[rel_type] = 0
            self.stats.relationship_types[rel_type] += 1

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

    def get_info(self) -> Dict[str, Any]:
        """Get analyzer information."""
        return {
            "config": {
                "confidence_threshold": self.confidence_threshold,
                "similarity_threshold": self.similarity_threshold,
                "pharmacophore_confidence": self.pharmacophore_confidence,
                "correlation_threshold": self.correlation_threshold,
            },
            "stats": self.stats.to_dict(),
        }

