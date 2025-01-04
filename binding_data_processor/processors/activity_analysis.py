"""Activity data analysis and SAR processing.

This module provides comprehensive functionality for:
1. Structure-activity relationship (SAR) analysis
2. Activity pattern analysis and statistics
3. Activity cliff and pharmacophore detection
4. Activity comparison and trend analysis
5. Statistical analysis and visualization
"""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union
import numpy as np
from scipy import stats
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, ChemicalFeatures
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from .activity_types import ActivityTypes
from .structure import StructureProcessor

logger = logging.getLogger(__name__)


class ActivityAnalysis:
    """Handles comprehensive analysis of activity data."""

    # Unit conversion factors to nM
    UNIT_CONVERSIONS = {
        "M": 1e9,
        "mM": 1e6,
        "μM": 1e3,
        "uM": 1e3,
        "nM": 1.0,
        "pM": 1e-3,
        "fM": 1e-6,
    }

    def __init__(
        self,
        structure_processor: Optional[StructureProcessor] = None,
        activity_types: Optional[ActivityTypes] = None,
    ):
        """
        Initialize activity analyzer.

        Args:
            structure_processor: Optional StructureProcessor instance
            activity_types: Optional ActivityTypes instance
        """
        self.logger = logging.getLogger(__name__)
        self.structure_processor = structure_processor or StructureProcessor()
        self.activity_types = activity_types or ActivityTypes()
        self.feature_factory = ChemicalFeatures.BuildFeatureFactory()

    def analyze_activity_patterns(
        self, compounds: List[Dict]
    ) -> Dict[str, Dict[str, float]]:
        """
        Analyze activity patterns across compounds.

        Args:
            compounds: List of compound dictionaries

        Returns:
            Dictionary of activity pattern statistics
        """
        try:
            # Get activity descriptions
            descriptions = [
                c.get("activity_description", "")
                for c in compounds
                if "activity_description" in c
            ]

            # Analyze distribution
            distribution = self.activity_types.analyze_activity_distribution(
                descriptions
            )

            # Calculate additional statistics
            analysis_stats = {
                "distribution": distribution,
                "summary": self._calculate_summary_stats(distribution),
                "trends": self._analyze_activity_trends(compounds),
            }

            return analysis_stats

        except Exception as e:
            self.logger.error(f"Error analyzing activity patterns: {str(e)}")
            return {}

    def process_activity_value(
        self, value: str, unit: str = "nM"
    ) -> Tuple[Optional[float], str, Optional[str]]:
        """
        Process activity value with modifiers and unit conversion.

        Args:
            value: Activity value string (e.g. ">100", "<0.1")
            unit: Unit of measurement

        Returns:
            Tuple of (numeric_value, standardized_unit, modifier)
        """
        try:
            # Handle modifiers
            modifier = None
            numeric_value = value.strip()

            if numeric_value.startswith(">"):
                modifier = ">"
                numeric_value = numeric_value[1:]
            elif numeric_value.startswith("<"):
                modifier = "<"
                numeric_value = numeric_value[1:]
            elif numeric_value.startswith("~"):
                modifier = "~"
                numeric_value = numeric_value[1:]

            # Convert to float
            numeric_value = float(numeric_value)

            # Standardize unit to nM
            unit = unit.strip()
            if unit in self.UNIT_CONVERSIONS:
                numeric_value *= self.UNIT_CONVERSIONS[unit]
                unit = "nM"
            else:
                self.logger.warning(f"Unknown unit: {unit}")

            return numeric_value, unit, modifier

        except ValueError as e:
            self.logger.error(f"Error converting value '{value}': {str(e)}")
            return None, unit, None
        except Exception as e:
            self.logger.error(f"Error processing activity value: {str(e)}")
            return None, unit, None

    def normalize_activities(
        self,
        values: List[Union[float, str]],
        reference_value: float,
        log_transform: bool = False,
    ) -> List[Optional[float]]:
        """
        Normalize activity values relative to reference.

        Args:
            values: List of activity values
            reference_value: Reference activity value
            log_transform: Whether to log transform values

        Returns:
            List of normalized values
        """
        try:
            normalized = []
            for value in values:
                # Handle string values with modifiers
                if isinstance(value, str):
                    processed_value, _, modifier = self.process_activity_value(value)
                    if processed_value is None:
                        normalized.append(None)
                        continue
                    value = processed_value

                    # Adjust for modifiers
                    if modifier == ">":
                        value *= 10  # Less active
                    elif modifier == "<":
                        value *= 0.1  # More active

                # Calculate normalized value
                if log_transform:
                    norm_value = np.log10(value / reference_value)
                else:
                    norm_value = value / reference_value

                normalized.append(norm_value)

            return normalized

        except Exception as e:
            self.logger.error(f"Error normalizing activities: {str(e)}")
            return [None] * len(values)

    def analyze_sar(
        self,
        compounds: List[Dict],
        activity_threshold: float = 100.0,
        similarity_threshold: float = 0.7,
    ) -> Dict:
        """
        Analyze structure-activity relationships.

        Args:
            compounds: List of compound dictionaries
            activity_threshold: Activity threshold in nM
            similarity_threshold: Structural similarity threshold

        Returns:
            Dictionary containing:
            - active_features: Frequency of structural features in active compounds
            - inactive_features: Frequency of structural features in inactive compounds
            - activity_ranges: Statistical ranges for each structural feature
            - structure_clusters: Structural similarity clusters
            - feature_contributions: Estimated contribution of each feature to activity
            - activity_cliffs: Identified activity cliffs between similar structures
            - pharmacophores: Common pharmacophore patterns
        """
        try:
            results = {
                "active_features": {},
                "inactive_features": {},
                "activity_ranges": {},
                "structure_clusters": [],
                "feature_contributions": {},
                "activity_cliffs": [],
                "pharmacophores": [],
            }

            # Process compounds and collect data
            processed_compounds = []
            all_features = set()
            feature_activities = {}

            for compound in compounds:
                if not all(k in compound for k in ["smiles", "activity_value"]):
                    continue

                # Convert SMILES to mol
                try:
                    mol = Chem.MolFromSmiles(compound["smiles"])
                    if mol is None:
                        self.logger.warning(f"Invalid SMILES: {compound['smiles']}")
                        continue
                except Exception as e:
                    self.logger.error(f"Error processing SMILES: {str(e)}")
                    continue

                # Get structural features
                features = self.structure_processor.analyze_structure(mol)
                all_features.update(features.keys())

                # Track activities per feature
                value = float(compound["activity_value"])
                for feature, count in features.items():
                    if feature not in feature_activities:
                        feature_activities[feature] = []
                    feature_activities[feature].append((value, count))

                # Store processed data
                processed_compounds.append(
                    {
                        "compound": compound,
                        "features": features,
                        "mol": mol,
                        "value": value,
                        "is_active": value <= activity_threshold,
                    }
                )

            if not processed_compounds:
                return results

            # Split into active/inactive
            active = [c for c in processed_compounds if c["is_active"]]
            inactive = [c for c in processed_compounds if not c["is_active"]]

            # Calculate feature frequencies
            if active:
                for compound in active:
                    for feature, count in compound["features"].items():
                        results["active_features"][feature] = (
                            results["active_features"].get(feature, 0) + count
                        )

            if inactive:
                for compound in inactive:
                    for feature, count in compound["features"].items():
                        results["inactive_features"][feature] = (
                            results["inactive_features"].get(feature, 0) + count
                        )

            # Calculate activity statistics per feature
            for feature in all_features:
                if feature in feature_activities:
                    activities = [a[0] for a in feature_activities[feature]]
                    if activities:
                        results["activity_ranges"][feature] = {
                            "min": float(np.min(activities)),
                            "max": float(np.max(activities)),
                            "mean": float(np.mean(activities)),
                            "median": float(np.median(activities)),
                            "std": float(np.std(activities)),
                            "count": len(activities),
                            "enrichment_factor": self._calculate_enrichment(
                                feature,
                                results["active_features"],
                                results["inactive_features"],
                                len(active),
                                len(inactive),
                            ),
                        }

            # Calculate feature contributions
            results["feature_contributions"] = self._calculate_contributions(
                all_features, active, inactive
            )

            # Detect activity cliffs
            results["activity_cliffs"] = self._detect_activity_cliffs(
                processed_compounds, similarity_threshold
            )

            # Generate pharmacophore patterns
            results["pharmacophores"] = self._analyze_pharmacophores(
                [c for c in processed_compounds if c["is_active"]]
            )

            # Perform clustering
            results["structure_clusters"] = self._cluster_compounds(
                processed_compounds, similarity_threshold
            )

            return results

        except Exception as e:
            self.logger.error(f"Error analyzing SAR: {str(e)}")
            return {}

    def calculate_activity_stats(
        self, values: List[Union[float, str]], by_type: bool = False
    ) -> Dict[str, Dict[str, float]]:
        """
        Calculate comprehensive activity statistics.

        Args:
            values: List of activity values
            by_type: Whether to group by activity type

        Returns:
            Dictionary of activity statistics
        """
        try:
            stats_dict = {}

            # Process values and handle modifiers
            processed_values = []
            for value in values:
                if isinstance(value, str):
                    processed_value, _, modifier = self.process_activity_value(value)
                    if processed_value is not None:
                        if modifier == ">":
                            processed_value *= 10
                        elif modifier == "<":
                            processed_value *= 0.1
                        processed_values.append(processed_value)
                else:
                    processed_values.append(float(value))

            if not processed_values:
                return stats_dict

            # Calculate basic statistics
            stats_dict["basic"] = {
                "count": len(processed_values),
                "min": float(np.min(processed_values)),
                "max": float(np.max(processed_values)),
                "mean": float(np.mean(processed_values)),
                "median": float(np.median(processed_values)),
                "std": float(np.std(processed_values)),
                "variance": float(np.var(processed_values)),
            }

            # Calculate percentiles
            percentiles = [10, 25, 50, 75, 90]
            stats_dict["percentiles"] = {
                f"p{p}": float(np.percentile(processed_values, p)) for p in percentiles
            }

            # Calculate distribution statistics
            stats_dict["distribution"] = {
                "skewness": float(stats.skew(processed_values)),
                "kurtosis": float(stats.kurtosis(processed_values)),
                "shapiro_stat": float(stats.shapiro(processed_values)[0]),
                "shapiro_p": float(stats.shapiro(processed_values)[1]),
            }

            return stats_dict

        except Exception as e:
            self.logger.error(f"Error calculating activity stats: {str(e)}")
            return {}

    def _calculate_summary_stats(
        self, distribution: Dict[str, Dict[str, float]]
    ) -> Dict[str, float]:
        """Calculate summary statistics from activity distribution."""
        try:
            total = sum(d["count"] for d in distribution.values())
            if total == 0:
                return {}

            # Calculate category frequencies
            category_counts = {}
            for activity, activity_stats in distribution.items():
                category = activity_stats.get("category")
                if category:
                    category_counts[category] = (
                        category_counts.get(category, 0) + activity_stats["count"]
                    )

            # Calculate diversity metrics
            activity_freqs = [
                activity_stats["frequency"] for activity_stats in distribution.values()
            ]
            entropy = stats.entropy(activity_freqs)

            return {
                "total_compounds": total,
                "unique_activities": len(distribution),
                "entropy": float(entropy),
                "category_distribution": {
                    cat: count / total for cat, count in category_counts.items()
                },
            }

        except Exception as e:
            self.logger.error(f"Error calculating summary stats: {str(e)}")
            return {}

    def _analyze_activity_trends(self, compounds: List[Dict]) -> Dict[str, List[Dict]]:
        """Analyze trends in activity data."""
        try:
            trends = {
                "activity_cliffs": [],
                "activity_switches": [],
                "activity_patterns": [],
            }

            # Sort compounds by activity value
            compounds = sorted(
                [c for c in compounds if "activity_value" in c],
                key=lambda x: float(x["activity_value"]),
            )

            if not compounds:
                return trends

            # Detect activity cliffs (large changes in activity)
            for i in range(len(compounds) - 1):
                curr_val = float(compounds[i]["activity_value"])
                next_val = float(compounds[i + 1]["activity_value"])
                if next_val / curr_val >= 10:  # 10-fold difference
                    trends["activity_cliffs"].append(
                        {
                            "compound1": compounds[i]["name"],
                            "compound2": compounds[i + 1]["name"],
                            "activity_ratio": next_val / curr_val,
                        }
                    )

            # Detect activity type switches
            for i in range(len(compounds) - 1):
                curr_type = compounds[i].get("activity_type")
                next_type = compounds[i + 1].get("activity_type")
                if curr_type and next_type and curr_type != next_type:
                    trends["activity_switches"].append(
                        {
                            "compound1": compounds[i]["name"],
                            "compound2": compounds[i + 1]["name"],
                            "type1": curr_type,
                            "type2": next_type,
                        }
                    )

            # Analyze activity patterns
            activity_values = [float(c["activity_value"]) for c in compounds]
            if len(activity_values) > 2:
                # Calculate basic statistics
                mean = np.mean(activity_values)
                std = np.std(activity_values)
                median = np.median(activity_values)

                # Detect patterns
                trends["activity_patterns"] = [
                    {
                        "pattern": "range",
                        "min": float(min(activity_values)),
                        "max": float(max(activity_values)),
                        "span": float(max(activity_values) - min(activity_values)),
                    },
                    {
                        "pattern": "distribution",
                        "mean": float(mean),
                        "median": float(median),
                        "std": float(std),
                        "skewness": float(stats.skew(activity_values)),
                        "kurtosis": float(stats.kurtosis(activity_values)),
                    },
                ]

            return trends

        except Exception as e:
            self.logger.error(f"Error analyzing activity trends: {str(e)}")
            return {
                "activity_cliffs": [],
                "activity_switches": [],
                "activity_patterns": [],
            }

    def _calculate_enrichment(
        self,
        feature: str,
        active_features: Dict[str, int],
        inactive_features: Dict[str, int],
        n_active: int,
        n_inactive: int,
    ) -> float:
        """Calculate enrichment factor for a feature."""
        try:
            if feature not in active_features or feature not in inactive_features:
                return 0.0

            active_freq = active_features[feature] / max(n_active, 1)
            inactive_freq = inactive_features[feature] / max(n_inactive, 1)

            return active_freq / max(inactive_freq, 0.001)

        except Exception as e:
            self.logger.error(f"Error calculating enrichment: {str(e)}")
            return 0.0

    def _calculate_contributions(
        self,
        features: Set[str],
        active: List[Dict],
        inactive: List[Dict],
    ) -> Dict[str, Dict]:
        """Calculate feature contributions to activity."""
        try:
            contributions = {}
            n_active = max(len(active), 1)
            n_inactive = max(len(inactive), 1)

            for feature in features:
                active_count = sum(1 for c in active if feature in c["features"])
                inactive_count = sum(1 for c in inactive if feature in c["features"])

                if active_count + inactive_count > 0:
                    contribution = (active_count / n_active) - (
                        inactive_count / n_inactive
                    )
                    contributions[feature] = {
                        "contribution_score": contribution,
                        "active_frequency": active_count / n_active,
                        "inactive_frequency": inactive_count / n_inactive,
                    }

            return contributions

        except Exception as e:
            self.logger.error(f"Error calculating contributions: {str(e)}")
            return {}

    def _detect_activity_cliffs(
        self,
        compounds: List[Dict],
        similarity_threshold: float = 0.7,
    ) -> List[Dict]:
        """Detect activity cliffs between compounds."""
        try:
            cliffs = []

            for i, comp1 in enumerate(compounds):
                fp1 = AllChem.GetMorganFingerprintAsBitVect(comp1["mol"], 2)

                for comp2 in compounds[i + 1 :]:
                    fp2 = AllChem.GetMorganFingerprintAsBitVect(comp2["mol"], 2)
                    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

                    if similarity >= similarity_threshold:
                        activity_ratio = abs(comp1["value"] / comp2["value"])
                        if activity_ratio >= 10:  # Significant activity difference
                            cliffs.append(
                                {
                                    "compound1": comp1["compound"]["name"],
                                    "compound2": comp2["compound"]["name"],
                                    "similarity": similarity,
                                    "activity_ratio": activity_ratio,
                                }
                            )

            return cliffs

        except Exception as e:
            self.logger.error(f"Error detecting activity cliffs: {str(e)}")
            return []

    def _analyze_pharmacophores(
        self,
        active_compounds: List[Dict],
        min_support: float = 0.25,
    ) -> List[Dict]:
        """Analyze pharmacophore patterns in active compounds."""
        try:
            if not active_compounds:
                return []

            common_features = {}
            n_compounds = len(active_compounds)

            for compound in active_compounds:
                features = self.feature_factory.GetFeaturesForMol(compound["mol"])
                pattern = tuple(sorted(f.GetFamily() for f in features))
                common_features[pattern] = common_features.get(pattern, 0) + 1

            # Filter by minimum support
            min_count = int(n_compounds * min_support)
            return [
                {
                    "pattern": list(pattern),
                    "frequency": count,
                    "support": count / n_compounds,
                }
                for pattern, count in common_features.items()
                if count >= min_count
            ]

        except Exception as e:
            self.logger.error(f"Error analyzing pharmacophores: {str(e)}")
            return []

    def _cluster_compounds(
        self,
        compounds: List[Dict],
        similarity_threshold: float = 0.7,
    ) -> List[Dict]:
        """Cluster compounds by structural similarity."""
        try:
            n_compounds = len(compounds)
            similarity_matrix = np.zeros((n_compounds, n_compounds))

            # Calculate similarity matrix
            for i in range(n_compounds):
                fp1 = AllChem.GetMorganFingerprintAsBitVect(compounds[i]["mol"], 2)
                for j in range(i + 1, n_compounds):
                    fp2 = AllChem.GetMorganFingerprintAsBitVect(compounds[j]["mol"], 2)
                    sim = DataStructs.TanimotoSimilarity(fp1, fp2)
                    similarity_matrix[i, j] = similarity_matrix[j, i] = sim

            # Perform clustering
            distances = 1 - similarity_matrix
            linkage = hierarchy.linkage(squareform(distances), method="complete")
            clusters = hierarchy.fcluster(
                linkage, 1 - similarity_threshold, criterion="distance"
            )

            # Organize clusters
            cluster_dict = {}
            for i, cluster_id in enumerate(clusters):
                if cluster_id not in cluster_dict:
                    cluster_dict[cluster_id] = []
                cluster_dict[cluster_id].append(compounds[i]["compound"])

            return [
                {
                    "cluster_id": cluster_id,
                    "compounds": compounds,
                    "mean_activity": float(
                        np.mean([float(c["activity_value"]) for c in compounds])
                    ),
                }
                for cluster_id, compounds in cluster_dict.items()
            ]

        except Exception as e:
            self.logger.error(f"Error clustering compounds: {str(e)}")
            return []
