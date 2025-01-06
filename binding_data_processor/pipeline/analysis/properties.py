"""Property analysis functionality.

This module provides the PropertyAnalyzer class that handles:
1. Chemical property calculation and prediction
2. Property relationship analysis
3. Property pattern detection
4. Property statistics
5. ML-based predictions
6. Web data enrichment
7. Integration with binding analysis
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Tuple
from dataclasses import dataclass, field
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor
import functools
import json

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, AllChem, rdDecomposition, rdMolTransforms
import numpy as np
from sklearn.ensemble import RandomForestRegressor
import joblib

from ...models import CompoundData
from ..ml import MLPredictor
from ..web import WebEnricher
from ...processors.structure.pharmacophore import PharmacophoreDetector
from ...processors.psychopharm.predictors.bbb.base import BBBPredictor


@dataclass
class PropertyStats:
    """Property analysis statistics."""

    # Analysis counts
    total_analyses: int = 0
    successful_analyses: int = 0
    failed_analyses: int = 0

    # Property stats
    property_ranges: Dict[str, Dict[str, float]] = field(default_factory=dict)
    property_distributions: Dict[str, Dict[str, int]] = field(default_factory=dict)

    # Pattern stats
    total_patterns: int = 0
    pattern_types: Dict[str, int] = field(default_factory=dict)

    # Relationship stats
    total_relationships: int = 0
    relationship_types: Dict[str, int] = field(default_factory=dict)

    # ML stats
    ml_predictions: Dict[str, Dict[str, float]] = field(default_factory=dict)
    prediction_accuracies: Dict[str, float] = field(default_factory=dict)

    # Web enrichment stats
    web_enrichments: Dict[str, int] = field(default_factory=dict)

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
            "properties": {
                "ranges": self.property_ranges,
                "distributions": self.property_distributions,
            },
            "patterns": {
                "total": self.total_patterns,
                "types": self.pattern_types,
            },
            "relationships": {
                "total": self.total_relationships,
                "types": self.relationship_types,
            },
            "ml": {
                "predictions": self.ml_predictions,
                "accuracies": self.prediction_accuracies,
            },
            "web": {
                "enrichments": self.web_enrichments,
            },
            "errors": self.errors,
        }

    def _get_success_rate(self) -> Optional[float]:
        """Get analysis success rate."""
        if not self.total_analyses:
            return None
        return self.successful_analyses / self.total_analyses


class PropertyAnalyzer:
    """Analyzer for chemical properties."""

    def __init__(
        self,
        cache_dir: Optional[Union[str, Path]] = None,
        config: Optional[Dict[str, Any]] = None,
        n_jobs: int = 4,
    ):
        """Initialize property analyzer.

        Args:
            cache_dir: Optional directory for caching
            config: Optional analyzer configuration
            n_jobs: Number of parallel jobs
        """
        self.logger = logging.getLogger(self.__class__.__name__)
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.config = config or {}
        self.n_jobs = n_jobs

        # Get config values
        self.property_ranges = self.config.get(
            "property_ranges",
            {
                "molecular_weight": {"min": 0, "max": 2000},
                "logp": {"min": -10, "max": 10},
                "tpsa": {"min": 0, "max": 500},
                "hbd": {"min": 0, "max": 20},
                "hba": {"min": 0, "max": 20},
                "rotatable_bonds": {"min": 0, "max": 50},
            },
        )

        # Initialize components
        self._init_ml_components()
        self._init_web_components()
        self._init_structure_components()

        # Initialize cache
        self._init_cache()

        # Initialize stats
        self.stats = PropertyStats()

    def _init_ml_components(self) -> None:
        """Initialize ML components."""
        self.ml_predictor = MLPredictor(
            cache_dir=self.cache_dir,
            config=self.config.get("ml_config"),
        )
        self.bbb_predictor = BBBPredictor(
            cache_dir=self.cache_dir,
            config=self.config.get("bbb_config"),
        )

    def _init_web_components(self) -> None:
        """Initialize web components."""
        self.web_enricher = WebEnricher(
            cache_dir=self.cache_dir,
            config=self.config.get("web_config"),
        )

    def _init_structure_components(self) -> None:
        """Initialize structure components."""
        self.pharmacophore_detector = PharmacophoreDetector(
            cache_dir=self.cache_dir,
            config=self.config.get("pharmacophore_config"),
        )

    def _init_cache(self) -> None:
        """Initialize cache."""
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            self.property_cache = self.cache_dir / "properties.json"
            if self.property_cache.exists():
                with open(self.property_cache) as f:
                    self._cache = json.load(f)
            else:
                self._cache = {}
        else:
            self._cache = {}

    def _cache_result(func):
        """Cache function results."""

        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            # Generate cache key
            key = f"{func.__name__}:{hash(str(args))}"

            # Check cache
            if key in self._cache:
                return self._cache[key]

            # Calculate result
            result = func(self, *args, **kwargs)

            # Update cache
            self._cache[key] = result
            if self.cache_dir:
                with open(self.property_cache, "w") as f:
                    json.dump(self._cache, f)

            return result

        return wrapper

    def analyze(
        self,
        compound: CompoundData,
        binding_data: Optional[Dict[str, Any]] = None,
        enrich_web: bool = True,
    ) -> CompoundData:
        """Analyze chemical properties.

        Args:
            compound: CompoundData instance to analyze
            binding_data: Optional binding analysis data
            enrich_web: Whether to enrich with web data

        Returns:
            Analyzed CompoundData instance
        """
        try:
            self.stats.total_analyses += 1

            # Get structure
            mol = Chem.MolFromSmiles(compound.smiles)
            if not mol:
                return compound

            # Calculate properties
            properties = self._calculate_properties(mol)
            compound.properties = properties

            # Update property stats
            self._update_property_stats(properties)

            # Analyze property patterns
            patterns = self._analyze_property_patterns(properties, binding_data)
            compound.property_patterns = patterns

            # Update pattern stats
            self._update_pattern_stats(patterns)

            # Analyze property relationships
            relationships = self._analyze_property_relationships(properties, binding_data)
            compound.property_relationships = relationships

            # Update relationship stats
            self._update_relationship_stats(relationships)

            # ML predictions
            predictions = self._make_ml_predictions(mol, properties, binding_data)
            compound.property_predictions = predictions

            # Update ML stats
            self._update_ml_stats(predictions)

            # Web enrichment
            if enrich_web:
                enrichments = self._enrich_web_data(compound)
                compound.property_enrichments = enrichments

                # Update web stats
                self._update_web_stats(enrichments)

            self.stats.successful_analyses += 1
            return compound

        except Exception as e:
            self.logger.error(f"Failed to analyze properties for {compound.name}: {str(e)}")
            self.stats.failed_analyses += 1
            self.stats.errors.append(
                {
                    "type": "property_analysis_error",
                    "compound": compound.name,
                    "error": str(e),
                    "timestamp": datetime.now().isoformat(),
                }
            )
            raise

    @_cache_result
    def _calculate_properties(
        self,
        mol: Chem.Mol,
    ) -> Dict[str, Any]:
        """Calculate chemical properties.

        Args:
            mol: RDKit molecule

        Returns:
            Dictionary of calculated properties
        """
        properties = {}

        # Basic properties
        properties["molecular_weight"] = Descriptors.ExactMolWt(mol)
        properties["logp"] = Crippen.MolLogP(mol)
        properties["tpsa"] = Descriptors.TPSA(mol)
        properties["hbd"] = rdMolDescriptors.CalcNumHBD(mol)
        properties["hba"] = rdMolDescriptors.CalcNumHBA(mol)
        properties["rotatable_bonds"] = Descriptors.NumRotatableBonds(mol)

        # Extended properties
        properties["aromatic_rings"] = Descriptors.NumAromaticRings(mol)
        properties["heavy_atoms"] = mol.GetNumHeavyAtoms()
        properties["fraction_csp3"] = Descriptors.FractionCSP3(mol)
        properties["qed"] = Descriptors.qed(mol)
        properties["sas"] = Descriptors.SAS(mol)
        properties["charge"] = Descriptors.NumRadicalElectrons(mol)

        # 3D properties
        try:
            mol_3d = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol_3d, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol_3d)
            properties["inertial_shape"] = rdMolTransforms.ComputePrincipalAxesAndMoments(mol_3d)
            properties["radius_gyration"] = rdMolTransforms.ComputeRadiusOfGyration(mol_3d)
        except (ValueError, RuntimeError) as e:
            self.logger.warning(f"Failed to calculate 3D properties: {str(e)}")

        # Structural decomposition
        try:
            properties["scaffolds"] = rdDecomposition.GetScaffolds(mol)
            properties["fragments"] = rdDecomposition.GetFragments(mol)
        except (ValueError, RuntimeError) as e:
            self.logger.warning(f"Failed to calculate structural decomposition: {str(e)}")

        # Drug-likeness properties
        properties["lipinski"] = self._check_lipinski(properties)
        properties["veber"] = self._check_veber(properties)
        properties["ghose"] = self._check_ghose(properties)
        properties["muegge"] = self._check_muegge(properties)
        properties["bbb"] = self._check_bbb(properties)
        properties["pfizer"] = self._check_pfizer(properties)
        properties["gsk"] = self._check_gsk(properties)
        properties["golden_triangle"] = self._check_golden_triangle(properties)

        return properties

    def _check_lipinski(
        self,
        properties: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Check Lipinski's Rule of Five.

        Args:
            properties: Dictionary of properties

        Returns:
            Lipinski analysis results
        """
        violations = []

        if properties["molecular_weight"] > 500:
            violations.append("molecular_weight > 500")
        if properties["logp"] > 5:
            violations.append("logp > 5")
        if properties["hbd"] > 5:
            violations.append("hbd > 5")
        if properties["hba"] > 10:
            violations.append("hba > 10")

        return {
            "pass": len(violations) <= 1,
            "violations": violations,
            "score": 1 - (len(violations) / 4),
        }

    def _check_veber(
        self,
        properties: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Check Veber rules.

        Args:
            properties: Dictionary of properties

        Returns:
            Veber analysis results
        """
        violations = []

        if properties["rotatable_bonds"] > 10:
            violations.append("rotatable_bonds > 10")
        if properties["tpsa"] > 140:
            violations.append("tpsa > 140")

        return {
            "pass": len(violations) == 0,
            "violations": violations,
            "score": 1 - (len(violations) / 2),
        }

    def _check_ghose(
        self,
        properties: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Check Ghose filter.

        Args:
            properties: Dictionary of properties

        Returns:
            Ghose analysis results
        """
        violations = []

        if not (160 <= properties["molecular_weight"] <= 480):
            violations.append("molecular_weight outside [160, 480]")
        if not (-0.4 <= properties["logp"] <= 5.6):
            violations.append("logp outside [-0.4, 5.6]")
        if not (20 <= properties["heavy_atoms"] <= 70):
            violations.append("heavy_atoms outside [20, 70]")

        return {
            "pass": len(violations) == 0,
            "violations": violations,
            "score": 1 - (len(violations) / 3),
        }

    def _check_muegge(
        self,
        properties: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Check Muegge filter.

        Args:
            properties: Dictionary of properties

        Returns:
            Muegge analysis results
        """
        violations = []

        if not (200 <= properties["molecular_weight"] <= 600):
            violations.append("molecular_weight outside [200, 600]")
        if not (-2 <= properties["logp"] <= 5):
            violations.append("logp outside [-2, 5]")
        if properties["tpsa"] < 75:
            violations.append("tpsa < 75")
        if properties["rotatable_bonds"] > 15:
            violations.append("rotatable_bonds > 15")

        return {
            "pass": len(violations) == 0,
            "violations": violations,
            "score": 1 - (len(violations) / 4),
        }

    def _check_bbb(
        self,
        properties: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Check blood-brain barrier rules.

        Args:
            properties: Dictionary of properties

        Returns:
            BBB analysis results
        """
        violations = []

        if properties["molecular_weight"] > 400:
            violations.append("molecular_weight > 400")
        if not (0 <= properties["logp"] <= 6):
            violations.append("logp outside [0, 6]")
        if properties["tpsa"] > 90:
            violations.append("tpsa > 90")
        if properties["hbd"] > 3:
            violations.append("hbd > 3")

        return {
            "pass": len(violations) == 0,
            "violations": violations,
            "score": 1 - (len(violations) / 4),
        }

    def _check_pfizer(
        self,
        properties: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Check Pfizer 3/75 rules.

        Args:
            properties: Dictionary of properties

        Returns:
            Pfizer analysis results
        """
        violations = []

        if properties["logp"] > 3:
            violations.append("logp > 3")
        if properties["tpsa"] < 75:
            violations.append("tpsa < 75")

        return {
            "pass": len(violations) == 0,
            "violations": violations,
            "score": 1 - (len(violations) / 2),
        }

    def _check_gsk(
        self,
        properties: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Check GSK 4/400 rules.

        Args:
            properties: Dictionary of properties

        Returns:
            GSK analysis results
        """
        violations = []

        if properties["molecular_weight"] > 400:
            violations.append("molecular_weight > 400")
        if properties["logp"] > 4:
            violations.append("logp > 4")

        return {
            "pass": len(violations) == 0,
            "violations": violations,
            "score": 1 - (len(violations) / 2),
        }

    def _check_golden_triangle(
        self,
        properties: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Check Golden Triangle rules.

        Args:
            properties: Dictionary of properties

        Returns:
            Golden Triangle analysis results
        """
        violations = []

        if not (200 <= properties["molecular_weight"] <= 500):
            violations.append("molecular_weight outside [200, 500]")
        if not (-2 <= properties["logp"] <= 5):
            violations.append("logp outside [-2, 5]")
        if properties["fraction_csp3"] < 0.4:
            violations.append("fraction_csp3 < 0.4")

        return {
            "pass": len(violations) == 0,
            "violations": violations,
            "score": 1 - (len(violations) / 3),
        }

    def _analyze_property_patterns(
        self,
        properties: Dict[str, Any],
        binding_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Analyze property patterns.

        Args:
            properties: Dictionary of properties
            binding_data: Optional binding analysis data

        Returns:
            List of property patterns
        """
        patterns = []

        # Basic property patterns
        patterns.extend(self._analyze_basic_patterns(properties))

        # Drug-likeness patterns
        patterns.extend(self._analyze_druglike_patterns(properties))

        # Structural patterns
        patterns.extend(self._analyze_structural_patterns(properties))

        # Binding-related patterns
        if binding_data:
            patterns.extend(self._analyze_binding_patterns(properties, binding_data))

        return patterns

    def _analyze_basic_patterns(
        self,
        properties: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Analyze basic property patterns.

        Args:
            properties: Dictionary of properties

        Returns:
            List of basic property patterns
        """
        patterns = []

        # Check property ranges
        for name, value in properties.items():
            if name not in self.property_ranges:
                continue

            ranges = self.property_ranges[name]
            if ranges["min"] <= value <= ranges["max"]:
                patterns.append(
                    {
                        "type": "range",
                        "property": name,
                        "value": value,
                        "range": ranges,
                    }
                )

        return patterns

    def _analyze_druglike_patterns(
        self,
        properties: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Analyze drug-likeness patterns.

        Args:
            properties: Dictionary of properties

        Returns:
            List of drug-likeness patterns
        """
        patterns = []

        # Check drug-likeness rules
        for rule in ["lipinski", "veber", "ghose", "muegge", "bbb", "pfizer", "gsk", "golden_triangle"]:
            if properties[rule]["pass"]:
                patterns.append(
                    {
                        "type": "drug_likeness",
                        "rule": rule,
                        "score": properties[rule]["score"],
                    }
                )

        return patterns

    def _analyze_structural_patterns(
        self,
        properties: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Analyze structural patterns.

        Args:
            properties: Dictionary of properties

        Returns:
            List of structural patterns
        """
        patterns = []

        # Check scaffolds
        if "scaffolds" in properties:
            patterns.append(
                {
                    "type": "structural",
                    "scaffolds": properties["scaffolds"],
                }
            )

        # Check fragments
        if "fragments" in properties:
            patterns.append(
                {
                    "type": "structural",
                    "fragments": properties["fragments"],
                }
            )

        return patterns

    def _analyze_binding_patterns(
        self,
        properties: Dict[str, Any],
        binding_data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Analyze binding-related patterns.

        Args:
            properties: Dictionary of properties
            binding_data: Binding analysis data

        Returns:
            List of binding-related patterns
        """
        patterns = []

        # Size-binding relationship
        if (
            properties["molecular_weight"] < 500
            and properties["rotatable_bonds"] <= 10
            and binding_data.get("primary_targets")
        ):
            patterns.append(
                {
                    "type": "size_binding",
                    "score": 1.0,
                    "properties": {
                        "molecular_weight": properties["molecular_weight"],
                        "rotatable_bonds": properties["rotatable_bonds"],
                    },
                    "binding": {
                        "targets": len(binding_data["primary_targets"]),
                    },
                }
            )

        # Polarity-binding relationship
        if 0 <= properties["logp"] <= 5 and properties["tpsa"] < 140 and binding_data.get("target_families"):
            patterns.append(
                {
                    "type": "polarity_binding",
                    "score": 1.0,
                    "properties": {
                        "logp": properties["logp"],
                        "tpsa": properties["tpsa"],
                    },
                    "binding": {
                        "families": binding_data["target_families"],
                    },
                }
            )

        return patterns

    def _analyze_property_relationships(
        self,
        properties: Dict[str, Any],
        binding_data: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Analyze property relationships.

        Args:
            properties: Dictionary of properties
            binding_data: Optional binding analysis data

        Returns:
            List of property relationships
        """
        relationships = []

        # Check size-complexity relationship
        if (
            properties["molecular_weight"] > 300
            and properties["rotatable_bonds"] > 5
            and properties["aromatic_rings"] > 2
        ):
            relationships.append(
                {
                    "type": "size_complexity",
                    "score": 1.0,
                    "properties": {
                        "molecular_weight": properties["molecular_weight"],
                        "rotatable_bonds": properties["rotatable_bonds"],
                        "aromatic_rings": properties["aromatic_rings"],
                    },
                }
            )

        # Check polarity-permeability relationship
        if (
            properties["tpsa"] < 140
            and 0 <= properties["logp"] <= 5
            and properties["hbd"] <= 5
            and properties["hba"] <= 10
        ):
            relationships.append(
                {
                    "type": "polarity_permeability",
                    "score": 1.0,
                    "properties": {
                        "tpsa": properties["tpsa"],
                        "logp": properties["logp"],
                        "hbd": properties["hbd"],
                        "hba": properties["hba"],
                    },
                }
            )

        # Check drug-likeness relationship
        if properties["lipinski"]["pass"] and properties["veber"]["pass"] and properties["qed"] > 0.5:
            relationships.append(
                {
                    "type": "drug_likeness",
                    "score": properties["qed"],
                    "properties": {
                        "lipinski": properties["lipinski"]["score"],
                        "veber": properties["veber"]["score"],
                        "qed": properties["qed"],
                    },
                }
            )

        # Check 3D shape relationship
        if "inertial_shape" in properties:
            relationships.append(
                {
                    "type": "3d_shape",
                    "score": 1.0,
                    "properties": {
                        "inertial_shape": properties["inertial_shape"],
                        "radius_gyration": properties.get("radius_gyration"),
                    },
                }
            )

        # Check binding-related relationships
        if binding_data:
            # Structure-activity relationship
            if (
                "scaffolds" in properties
                and binding_data.get("primary_targets")
                and binding_data.get("selectivity_ratios")
            ):
                relationships.append(
                    {
                        "type": "structure_activity",
                        "score": 1.0,
                        "properties": {
                            "scaffolds": properties["scaffolds"],
                            "fragments": properties.get("fragments", []),
                        },
                        "binding": {
                            "primary_targets": binding_data["primary_targets"],
                            "selectivity": binding_data["selectivity_ratios"],
                        },
                    }
                )

            # 3D-binding relationship
            if "inertial_shape" in properties and binding_data.get("binding_types"):
                relationships.append(
                    {
                        "type": "3d_binding",
                        "score": 1.0,
                        "properties": {
                            "inertial_shape": properties["inertial_shape"],
                            "radius_gyration": properties.get("radius_gyration"),
                        },
                        "binding": {
                            "types": binding_data["binding_types"],
                        },
                    }
                )

        return relationships

    def _make_ml_predictions(
        self,
        mol: Chem.Mol,
        properties: Dict[str, Any],
        binding_data: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Make ML predictions.

        Args:
            mol: RDKit molecule
            properties: Dictionary of properties
            binding_data: Optional binding analysis data

        Returns:
            Dictionary of predictions
        """
        predictions = {}

        # BBB predictions
        bbb_pred = self.bbb_predictor.predict(
            mol,
            properties=properties,
            binding_data=binding_data,
        )
        predictions["bbb"] = {
            "probability": bbb_pred.probability,
            "confidence": bbb_pred.confidence,
            "mechanisms": bbb_pred.mechanisms,
        }

        # Property predictions
        prop_preds = self.ml_predictor.predict_properties(
            mol,
            binding_data=binding_data,
        )
        predictions["properties"] = prop_preds

        return predictions

    def _enrich_web_data(
        self,
        compound: CompoundData,
    ) -> Dict[str, Any]:
        """Enrich with web data.

        Args:
            compound: CompoundData instance

        Returns:
            Dictionary of web enrichments
        """
        enrichments = {}

        # Get web data
        web_data = self.web_enricher.enrich(compound)

        # Extract property-related data
        if "properties" in web_data:
            enrichments["properties"] = web_data["properties"]

        if "drug_likeness" in web_data:
            enrichments["drug_likeness"] = web_data["drug_likeness"]

        if "safety" in web_data:
            enrichments["safety"] = web_data["safety"]

        return enrichments

    def _update_property_stats(
        self,
        properties: Dict[str, Any],
    ) -> None:
        """Update property statistics."""
        for name, value in properties.items():
            if isinstance(value, (int, float)):
                # Update range stats
                if name not in self.stats.property_ranges:
                    self.stats.property_ranges[name] = {
                        "min": float("inf"),
                        "max": float("-inf"),
                        "sum": 0,
                        "count": 0,
                    }
                ranges = self.stats.property_ranges[name]
                ranges["min"] = min(ranges["min"], value)
                ranges["max"] = max(ranges["max"], value)
                ranges["sum"] += value
                ranges["count"] += 1

                # Update distribution stats
                if name not in self.stats.property_distributions:
                    self.stats.property_distributions[name] = {}
                dist = self.stats.property_distributions[name]
                bin_value = round(value, 1)
                if bin_value not in dist:
                    dist[bin_value] = 0
                dist[bin_value] += 1

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

    def _update_ml_stats(
        self,
        predictions: Dict[str, Any],
    ) -> None:
        """Update ML statistics."""
        for pred_type, pred_data in predictions.items():
            if pred_type not in self.stats.ml_predictions:
                self.stats.ml_predictions[pred_type] = {
                    "total": 0,
                    "positive": 0,
                    "negative": 0,
                }
            stats = self.stats.ml_predictions[pred_type]
            stats["total"] += 1
            if isinstance(pred_data, dict) and "probability" in pred_data:
                if pred_data["probability"] > 0.5:
                    stats["positive"] += 1
                else:
                    stats["negative"] += 1

    def _update_web_stats(
        self,
        enrichments: Dict[str, Any],
    ) -> None:
        """Update web statistics."""
        for enrich_type in enrichments:
            if enrich_type not in self.stats.web_enrichments:
                self.stats.web_enrichments[enrich_type] = 0
            self.stats.web_enrichments[enrich_type] += 1

    def get_info(self) -> Dict[str, Any]:
        """Get analyzer information."""
        return {
            "config": {
                "property_ranges": self.property_ranges,
                "cache_dir": str(self.cache_dir) if self.cache_dir else None,
                "n_jobs": self.n_jobs,
            },
            "components": {
                "ml_predictor": self.ml_predictor.get_info(),
                "bbb_predictor": self.bbb_predictor.get_info(),
                "web_enricher": self.web_enricher.get_info(),
                "pharmacophore_detector": self.pharmacophore_detector.get_info(),
            },
            "stats": self.stats.to_dict(),
        }
