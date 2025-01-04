"""Activity type classification and pattern matching.

This module provides comprehensive functionality for:
1. Activity type pattern definitions and matching
2. Activity type classification with confidence scores
3. Activity type hierarchies and relationships
4. Activity type validation and standardization
5. Activity distribution analysis
"""

import logging
import re
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


class ActivityTypes:
    """Handles activity type classification and pattern matching."""

    # Activity type hierarchy with confidence scores
    ACTIVITY_HIERARCHY = {
        "agonist": {
            "superagonist": 1.0,
            "full_agonist": 0.8,
            "partial_agonist": 0.6,
            "weak_partial_agonist": 0.4,
        },
        "antagonist": {
            "competitive_antagonist": 1.0,
            "noncompetitive_antagonist": 0.8,
            "inverse_agonist": 0.6,
            "neutral_antagonist": 0.4,
        },
        "modulator": {
            "positive_allosteric_modulator": 1.0,
            "negative_allosteric_modulator": 1.0,
            "allosteric_agonist": 0.8,
            "allosteric_antagonist": 0.8,
            "complex_modulator": 0.6,
        },
        "enzyme": {
            "enzyme_inhibitor": {
                "competitive_inhibitor": 1.0,
                "noncompetitive_inhibitor": 0.8,
                "uncompetitive_inhibitor": 0.8,
                "mixed_inhibitor": 0.6,
                "irreversible_inhibitor": 1.0,
                "mechanism_based_inhibitor": 1.0,
                "slow_binding_inhibitor": 0.8,
                "tight_binding_inhibitor": 1.0,
            },
            "enzyme_inducer": {
                "transcriptional_inducer": 1.0,
                "post_translational_inducer": 0.8,
                "allosteric_activator": 0.8,
                "cofactor_mediated": 0.6,
            },
        },
    }

    # Activity type categories for grouping related types
    ACTIVITY_CATEGORIES = {
        "agonist": {
            "superagonist",
            "full_agonist",
            "partial_agonist",
            "weak_partial_agonist",
        },
        "antagonist": {
            "competitive_antagonist",
            "noncompetitive_antagonist",
            "inverse_agonist",
            "neutral_antagonist",
        },
        "modulator": {
            "positive_allosteric_modulator",
            "negative_allosteric_modulator",
            "allosteric_agonist",
            "allosteric_antagonist",
            "complex_modulator",
        },
        "enzyme": {
            "enzyme_inhibitor",
            "enzyme_inducer",
            "competitive_inhibitor",
            "noncompetitive_inhibitor",
            "uncompetitive_inhibitor",
            "mixed_inhibitor",
            "irreversible_inhibitor",
            "mechanism_based_inhibitor",
            "slow_binding_inhibitor",
            "tight_binding_inhibitor",
            "transcriptional_inducer",
            "post_translational_inducer",
            "allosteric_activator",
            "cofactor_mediated",
        },
    }

    # Combined activity type patterns with enhanced regex
    ACTIVITY_PATTERNS = {
        # Agonist patterns
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
        # Antagonist patterns
        "competitive_antagonist": [
            r"competitive.?antagonist",
            r"orthosteric.?antagonist",
            r"reversible.?antagonist",
            r"binding.?site.?competition",
        ],
        "noncompetitive_antagonist": [
            r"non.?competitive.?antagonist",
            r"allosteric.?antagonist",
            r"binding.?site.?distinct",
        ],
        "inverse_agonist": [
            r"inverse.?agonist",
            r"negative.?agonist",
            r"inverse.?activity",
            r"negative.?efficacy",
            r"constitutive.?inhibitor",
            r"inverse.?effect",
        ],
        "neutral_antagonist": [
            r"neutral.?antagonist",
            r"silent.?antagonist",
            r"pure.?antagonist",
        ],
        # Modulator patterns
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
        "complex_modulator": [
            r"complex.?modulator",
            r"mixed.?modulator",
            r"complex.?pharmacology",
            r"bitopic",
            r"dual.?mechanism",
            r"complex.?binding",
        ],
        # Enzyme inhibitor patterns
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
        # Enzyme inducer patterns
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

    def __init__(self):
        """Initialize activity types processor."""
        self.logger = logging.getLogger(__name__)
        # Pre-compile regex patterns
        self.compiled_patterns = {
            activity: [re.compile(p, re.I) for p in patterns]
            for activity, patterns in self.ACTIVITY_PATTERNS.items()
        }

    def determine_activity_type(
        self, description: str, confidence_threshold: float = 0.5
    ) -> Tuple[str, float]:
        """
        Determine activity type from assay description.

        Args:
            description: Assay description text
            confidence_threshold: Minimum confidence score

        Returns:
            Tuple of (activity_type, confidence_score)
        """
        try:
            description = description.lower()
            matches = {}

            # Check each activity type pattern
            for activity_type, patterns in self.compiled_patterns.items():
                match_count = 0
                total_patterns = len(patterns)

                for pattern in patterns:
                    if pattern.search(description):
                        match_count += 1

                if match_count > 0:
                    confidence = match_count / total_patterns
                    matches[activity_type] = confidence

            # Get best match above threshold
            if matches:
                best_match = max(matches.items(), key=lambda x: x[1])
                if best_match[1] >= confidence_threshold:
                    return best_match

            return "unknown", 0.0

        except Exception as e:
            self.logger.error(f"Activity type determination error: {str(e)}")
            return "unknown", 0.0

    def analyze_activity_distribution(
        self, descriptions: List[str]
    ) -> Dict[str, Dict[str, float]]:
        """
        Analyze distribution of activity types.

        Args:
            descriptions: List of assay descriptions

        Returns:
            Dictionary mapping activity types to statistics
        """
        try:
            counts = {activity: 0 for activity in self.ACTIVITY_PATTERNS}
            counts["unknown"] = 0
            total = len(descriptions)

            # Count occurrences
            for desc in descriptions:
                activity_type, confidence = self.determine_activity_type(desc)
                counts[activity_type] += 1

            # Calculate statistics
            stats = {}
            for activity_type, count in counts.items():
                category = self.get_activity_category(activity_type)
                stats[activity_type] = {
                    "count": count,
                    "frequency": count / total if total > 0 else 0,
                    "category": category,
                    "related": list(self.get_related_activities(activity_type)),
                }

            return stats

        except Exception as e:
            self.logger.error(f"Error analyzing activity distribution: {str(e)}")
            return {}

    def get_activity_hierarchy(self, activity_type: str) -> List[str]:
        """
        Get hierarchy of activity types (e.g. weak_partial_agonist -> partial_agonist -> agonist).

        Args:
            activity_type: Activity type string

        Returns:
            List of activity types in hierarchy order
        """
        try:
            hierarchy = []

            if activity_type == "unknown":
                return hierarchy

            hierarchy.append(activity_type)

            # Add parent activity types
            if activity_type.startswith("weak_partial_"):
                hierarchy.append(activity_type.replace("weak_", ""))

            if activity_type.startswith("partial_"):
                hierarchy.append(activity_type.replace("partial_", ""))

            if activity_type.endswith("_agonist"):
                hierarchy.append("agonist")

            if activity_type.endswith("_antagonist"):
                hierarchy.append("antagonist")

            if "allosteric" in activity_type:
                hierarchy.append("allosteric_modulator")

            return hierarchy

        except Exception as e:
            logger.error(f"Error getting activity hierarchy: {str(e)}")
            return [activity_type]

    def get_activity_category(self, activity_type: str) -> Optional[str]:
        """
        Get category for an activity type.

        Args:
            activity_type: Activity type string

        Returns:
            Category string or None if not found
        """
        try:
            for category, types in self.ACTIVITY_CATEGORIES.items():
                if activity_type in types:
                    return category
            return None

        except Exception as e:
            self.logger.error(f"Error getting activity category: {str(e)}")
            return None

    def get_parent_activity(self, activity_type: str) -> Optional[str]:
        """
        Get parent activity type from hierarchy.

        Args:
            activity_type: Activity type to look up

        Returns:
            Parent activity type or None
        """
        try:
            # Search through hierarchy
            for parent, children in self.ACTIVITY_HIERARCHY.items():
                if isinstance(children, dict):
                    if activity_type in children:
                        return parent
                    # Check nested dictionaries
                    for child, grandchildren in children.items():
                        if isinstance(grandchildren, dict):
                            if activity_type in grandchildren:
                                return child
            return None

        except Exception as e:
            self.logger.error(f"Parent activity lookup error: {str(e)}")
            return None

    def get_related_activities(
        self, activity_type: str, include_children: bool = True
    ) -> Set[str]:
        """
        Get related activity types.

        Args:
            activity_type: Activity type to find relations for
            include_children: Whether to include child activities

        Returns:
            Set of related activity types
        """
        try:
            related = set()

            # Get activities in same category
            category = self.get_activity_category(activity_type)
            if category:
                related.update(self.ACTIVITY_CATEGORIES[category] - {activity_type})

            # Get parent and siblings
            parent = self.get_parent_activity(activity_type)
            if parent:
                related.add(parent)
                if parent in self.ACTIVITY_HIERARCHY:
                    related.update(self.ACTIVITY_HIERARCHY[parent].keys())

            if include_children:
                # Add child activities
                for parent, children in self.ACTIVITY_HIERARCHY.items():
                    if isinstance(children, dict):
                        if activity_type == parent:
                            related.update(children.keys())
                        # Check nested dictionaries
                        for child, grandchildren in children.items():
                            if isinstance(grandchildren, dict):
                                if activity_type == child:
                                    related.update(grandchildren.keys())

            return related - {activity_type}

        except Exception as e:
            self.logger.error(f"Error getting related activities: {str(e)}")
            return set()

    def validate_activity_type(self, activity_type: str) -> bool:
        """
        Validate activity type against known types.

        Args:
            activity_type: Activity type to validate

        Returns:
            Whether activity type is valid
        """
        try:
            # Check if activity type is in patterns
            if activity_type in self.ACTIVITY_PATTERNS:
                return True

            # Check categories
            for types in self.ACTIVITY_CATEGORIES.values():
                if activity_type in types:
                    return True

            # Check hierarchy
            for parent, children in self.ACTIVITY_HIERARCHY.items():
                if activity_type == parent:
                    return True
                if isinstance(children, dict):
                    if activity_type in children:
                        return True
                    # Check nested dictionaries
                    for child, grandchildren in children.items():
                        if isinstance(grandchildren, dict):
                            if activity_type in grandchildren:
                                return True
            return False

        except Exception as e:
            self.logger.error(f"Activity type validation error: {str(e)}")
            return False

    def standardize_activity_type(self, activity_type: str) -> str:
        """
        Standardize activity type name.

        Args:
            activity_type: Activity type to standardize

        Returns:
            Standardized activity type
        """
        try:
            # Remove special characters and extra whitespace
            activity_type = re.sub(r"[^\w\s-]", "", activity_type)
            activity_type = re.sub(r"\s+", "_", activity_type.strip())
            activity_type = activity_type.lower()

            # Check if standardized type is valid
            if self.validate_activity_type(activity_type):
                return activity_type

            # Try to match against patterns
            for std_type, patterns in self.compiled_patterns.items():
                for pattern in patterns:
                    if pattern.search(activity_type):
                        return std_type

            return "unknown"

        except Exception as e:
            self.logger.error(f"Activity type standardization error: {str(e)}")
            return "unknown"
