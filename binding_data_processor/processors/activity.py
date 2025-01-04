"""Binding activity data processing and analysis.

This module provides comprehensive functionality for:
1. Processing binding affinity data with unit conversions
2. Determining activity types from assay descriptions
3. Analyzing structure-activity relationships
4. Calculating activity statistics
5. Normalizing activity values
6. Handling multiple data sources
"""

import re
import logging
from typing import Dict, List, Optional, Set, Tuple
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from .structure import StructureProcessor


class ActivityProcessor:
    """Handles binding activity data processing."""

    # Combined activity type patterns with enhanced regex
    ACTIVITY_PATTERNS = {
        "superagonist": [
            r"super.?agonist",
            r"high.?efficacy.?agonist",
            r"full.?agonist.+high.?efficacy",
            r"efficacy\s*>\s*100%",
            r"super.?potent.?agonist",
            r"ultra.?potent.?agonist",
        ],
        "full_agonist": [
            r"full.?agonist",
            r"complete.?agonist",
            r"full.?receptor.?activation",
            r"efficacy\s*[~≈≃]\s*100%",
            r"maximal.?response",
            r"full.?efficacy",
        ],
        "partial_agonist": [
            r"partial.?agonist",
            r"submaximal.?activation",
            r"partial.?receptor.?activation",
            r"efficacy\s*[<≈]\s*\d{1,2}%",
            r"partial.?response",
            r"intermediate.?efficacy",
        ],
        "weak_partial_agonist": [
            r"weak.?partial.?agonist",
            r"low.?efficacy.?partial",
            r"weak.?partial.?activation",
            r"efficacy\s*<\s*20%",
            r"minimal.?agonist",
            r"very.?low.?efficacy",
        ],
        "mixed_agonist_antagonist": [
            r"mixed.?agonist.?antagonist",
            r"partial.?agonist.?antagonist",
            r"dual.?activity",
            r"context.?dependent",
            r"tissue.?dependent",
            r"mixed.?profile",
        ],
        "antagonist": [
            r"antagonist",
            r"blocker",
            r"inhibitor",
            r"neutral.?antagonist",
            r"competitive.?antagonist",
            r"receptor.?blocker",
        ],
        "inverse_agonist": [
            r"inverse.?agonist",
            r"negative.?agonist",
            r"inverse.?activity",
            r"negative.?efficacy",
            r"constitutive.?inhibitor",
            r"inverse.?effect",
        ],
        "positive_allosteric_modulator": [
            r"positive.?allosteric",
            r"PAM",
            r"positive.?modulator",
            r"allosteric.?potentiator",
            r"positive.?cooperativity",
            r"allosteric.?enhancer",
        ],
        "negative_allosteric_modulator": [
            r"negative.?allosteric",
            r"NAM",
            r"negative.?modulator",
            r"allosteric.?inhibitor",
            r"negative.?cooperativity",
            r"allosteric.?blocker",
        ],
        "enzyme_inhibitor": [
            r"enzyme.?inhibitor",
            r"inhibits?.?\w+.?enzyme",
            r"inhibits?.?\w+.?activity",
            r"reduces?.?enzyme.?activity",
            r"blocks?.?enzyme.?function",
            r"enzyme.?blocker",
            r"enzymatic.?inhibition",
            r"competitive.?inhibition",
            r"noncompetitive.?inhibition",
            r"uncompetitive.?inhibition",
            r"irreversible.?inhibition",
            r"mechanism.?based.?inhibition",
            r"suicide.?inhibition",
            r"slow.?binding.?inhibition",
            r"tight.?binding.?inhibition",
            r"allosteric.?enzyme.?inhibition",
            r"mixed.?type.?inhibition",
            r"partial.?inhibition",
            r"time.?dependent.?inhibition",
            r"substrate.?analog.?inhibition",
        ],
        "enzyme_inducer": [
            r"enzyme.?inducer",
            r"induces?.?\w+.?enzyme",
            r"increases?.?enzyme.?activity",
            r"enhances?.?enzyme.?function",
            r"upregulates?.?enzyme",
            r"enzyme.?activator",
            r"enzymatic.?induction",
            r"enzyme.?expression",
            r"enzyme.?upregulation",
            r"metabolic.?inducer",
            r"transcriptional.?activation",
            r"post.?translational.?activation",
            r"allosteric.?enzyme.?activation",
            r"enzyme.?stabilization",
            r"cofactor.?mediated.?activation",
            r"substrate.?mediated.?activation",
            r"positive.?enzyme.?regulation",
            r"enzyme.?synthesis.?inducer",
            r"catalytic.?enhancement",
            r"enzyme.?potentiation",
        ],
    }

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

    def __init__(self, structure_processor: Optional[StructureProcessor] = None):
        """
        Initialize activity processor.

        Args:
            structure_processor: Optional StructureProcessor instance for SAR analysis
        """
        self.logger = logging.getLogger(__name__)
        self.structure_processor = structure_processor or StructureProcessor()
        # Pre-compile regex patterns
        self.compiled_patterns = {
            activity: [re.compile(p, re.I) for p in patterns]
            for activity, patterns in self.ACTIVITY_PATTERNS.items()
        }

    def process_affinity_data(
        self, value: str, unit: str = "nM"
    ) -> Tuple[Optional[float], str, Optional[str]]:
        """
        Process binding affinity value and unit.

        Args:
            value: Affinity value as string (e.g. ">100", "<0.1")
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
            self.logger.error(f"Error processing affinity data: {str(e)}")
            return None, unit, None

    def determine_activity_type(self, description: str) -> str:
        """
        Determine activity type from assay description.

        Args:
            description: Assay description text

        Returns:
            Activity type string
        """
        description = description.lower()

        for activity_type, patterns in self.compiled_patterns.items():
            for pattern in patterns:
                if pattern.search(description):
                    return activity_type

        return "unknown"

    def normalize_activities(
        self,
        compounds: List[Dict],
        reference_value: float,
        log_transform: bool = False,
    ) -> List[Dict]:
        """
        Normalize activity values relative to reference.

        Args:
            compounds: List of compound dictionaries
            reference_value: Reference activity value
            log_transform: Whether to log transform values

        Returns:
            List of compounds with normalized activities
        """
        try:
            normalized = []
            for compound in compounds:
                if "activity_value" not in compound:
                    continue

                norm_compound = compound.copy()
                value = float(compound["activity_value"])

                # Handle modifiers
                modifier = compound.get("activity_modifier")
                if modifier == ">":
                    # For '>' values, use 10x for normalization
                    value *= 10
                elif modifier == "<":
                    # For '<' values, use 0.1x for normalization
                    value *= 0.1

                # Calculate normalized value
                if log_transform:
                    norm_value = np.log10(value / reference_value)
                else:
                    norm_value = value / reference_value

                norm_compound["normalized_activity"] = norm_value
                normalized.append(norm_compound)

            return normalized

        except Exception as e:
            self.logger.error(f"Error normalizing activities: {str(e)}")
            return compounds

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

            # Convert compounds to RDKit mols with validation
            mols = []
            activities = []
            valid_compounds = []

            for compound in compounds:
                if not all(k in compound for k in ["smiles", "activity_value"]):
                    continue

                mol = Chem.MolFromSmiles(compound["smiles"])
                if mol is None:
                    continue

                mols.append(mol)
                activities.append(float(compound["activity_value"]))
                valid_compounds.append(compound)

            if not mols:
                return results

            # Calculate all-vs-all similarity matrix
            n_mols = len(mols)
            similarity_matrix = np.zeros((n_mols, n_mols))
            for i in range(n_mols):
                for j in range(i + 1, n_mols):
                    sim = self.structure_processor.calculate_similarity(
                        mols[i], mols[j], method="tanimoto", fp_type="morgan"
                    )
                    similarity_matrix[i, j] = similarity_matrix[j, i] = sim or 0.0

            # Detect activity cliffs
            for i in range(n_mols):
                for j in range(i + 1, n_mols):
                    if similarity_matrix[i, j] >= similarity_threshold:
                        activity_ratio = max(
                            activities[i] / activities[j],
                            activities[j] / activities[i],
                        )
                        if activity_ratio >= 10.0:  # At least 10-fold difference
                            results["activity_cliffs"].append(
                                {
                                    "compound1": valid_compounds[i],
                                    "compound2": valid_compounds[j],
                                    "similarity": similarity_matrix[i, j],
                                    "activity_ratio": activity_ratio,
                                }
                            )

            # Process compounds and analyze structural features
            active_compounds = []
            inactive_compounds = []
            active_features = {}
            inactive_features = {}
            feature_activities = {}
            all_features = set()

            for mol, activity, compound in zip(mols, activities, valid_compounds):
                # Get structural features
                features = self.structure_processor.analyze_structure(mol)
                all_features.update(features.keys())
                compound_data = (compound, features, mol)

                # Categorize compound and update feature statistics
                is_active = activity <= activity_threshold
                target_features = active_features if is_active else inactive_features
                target_compounds = active_compounds if is_active else inactive_compounds

                for feature, count in features.items():
                    target_features[feature] = target_features.get(feature, 0) + count
                    if feature not in feature_activities:
                        feature_activities[feature] = []
                    feature_activities[feature].append(activity)

                target_compounds.append(compound_data)

            # Store feature frequencies
            results["active_features"] = active_features
            results["inactive_features"] = inactive_features

            # Calculate feature contributions
            total_compounds = len(active_compounds) + len(inactive_compounds)
            if total_compounds > 0:
                for feature in all_features:
                    active_count = sum(
                        1 for _, feats, _ in active_compounds if feature in feats
                    )
                    inactive_count = sum(
                        1 for _, feats, _ in inactive_compounds if feature in feats
                    )
                    if active_count + inactive_count > 0:
                        contribution = (
                            active_count / max(len(active_compounds), 1)
                        ) - (inactive_count / max(len(inactive_compounds), 1))
                        results["feature_contributions"][feature] = {
                            "contribution_score": contribution,
                            "active_frequency": active_count
                            / max(len(active_compounds), 1),
                            "inactive_frequency": inactive_count
                            / max(len(inactive_compounds), 1),
                        }

            # Calculate activity ranges and enrichment factors
            for feature, activities in feature_activities.items():
                if activities:
                    results["activity_ranges"][feature] = {
                        "min": float(np.min(activities)),
                        "max": float(np.max(activities)),
                        "mean": float(np.mean(activities)),
                        "median": float(np.median(activities)),
                        "std": float(np.std(activities)),
                        "count": len(activities),
                        "enrichment_factor": (
                            (
                                active_features.get(feature, 0)
                                / max(len(active_features), 1)
                            )
                            / (
                                inactive_features.get(feature, 0)
                                / max(len(inactive_features), 1)
                            )
                            if feature in active_features
                            and feature in inactive_features
                            else 0
                        ),
                    }

            results["active_features"] = active_features
            results["inactive_features"] = inactive_features

            # Perform hierarchical clustering
            distances = 1 - similarity_matrix
            linkage = hierarchy.linkage(squareform(distances), method="complete")
            clusters = hierarchy.fcluster(
                linkage, similarity_threshold, criterion="distance"
            )

            # Organize compounds into clusters
            cluster_dict = {}
            for i, cluster_id in enumerate(clusters):
                if cluster_id not in cluster_dict:
                    cluster_dict[cluster_id] = []
                cluster_dict[cluster_id].append(valid_compounds[i])

            results["structure_clusters"] = [
                {
                    "cluster_id": cluster_id,
                    "compounds": compounds,
                    "mean_activity": float(
                        np.mean([float(c["activity_value"]) for c in compounds])
                    ),
                }
                for cluster_id, compounds in cluster_dict.items()
            ]

            # Generate pharmacophore patterns
            if valid_compounds:
                from rdkit.Chem import ChemicalFeatures

                factory = ChemicalFeatures.BuildFeatureFactory()
                common_features = {}

                for mol in mols:
                    features = factory.GetFeaturesForMol(mol)
                    feature_pattern = tuple(sorted(f.GetFamily() for f in features))
                    common_features[feature_pattern] = (
                        common_features.get(feature_pattern, 0) + 1
                    )

                # Keep patterns present in at least 25% of compounds
                min_support = len(mols) * 0.25
                results["pharmacophores"] = [
                    {
                        "pattern": list(pattern),
                        "frequency": count,
                        "support": count / len(mols),
                    }
                    for pattern, count in common_features.items()
                    if count >= min_support
                ]

            return results

        except Exception as e:
            self.logger.error(f"Error analyzing SAR: {str(e)}")
            return {}

    def calculate_activity_stats(self, compounds: List[Dict]) -> Dict:
        """
        Calculate comprehensive activity statistics.

        Args:
            compounds: List of compound dictionaries

        Returns:
            Dictionary containing activity statistics
        """
        try:
            stats = {
                "activity_types": {},
                "value_ranges": {},
                "correlations": {},
            }

            # Count activity types
            for compound in compounds:
                activity_type = compound.get("activity_type", "unknown")
                stats["activity_types"][activity_type] = (
                    stats["activity_types"].get(activity_type, 0) + 1
                )

            # Get activity values
            values = []
            for compound in compounds:
                if "activity_value" in compound:
                    value = float(compound["activity_value"])
                    modifier = compound.get("activity_modifier")

                    # Handle modifiers
                    if modifier == ">":
                        value *= 10
                    elif modifier == "<":
                        value *= 0.1

                    values.append(value)

            if values:
                # Calculate value ranges
                stats["value_ranges"] = {
                    "min": float(np.min(values)),
                    "max": float(np.max(values)),
                    "mean": float(np.mean(values)),
                    "median": float(np.median(values)),
                    "std": float(np.std(values)),
                }

                # Calculate correlations with descriptors
                descriptor_values = {}
                for compound in compounds:
                    if (
                        "activity_value" not in compound
                        or "descriptors" not in compound
                    ):
                        continue

                    activity = float(compound["activity_value"])
                    for name, value in compound["descriptors"].items():
                        if name not in descriptor_values:
                            descriptor_values[name] = []
                        descriptor_values[name].append((activity, float(value)))

                # Calculate correlation coefficients
                for name, values in descriptor_values.items():
                    if len(values) > 1:
                        activities, descriptors = zip(*values)
                        correlation = float(np.corrcoef(activities, descriptors)[0, 1])
                        stats["correlations"][name] = correlation

            return stats

        except Exception as e:
            self.logger.error(f"Error calculating activity stats: {str(e)}")
            return {}
