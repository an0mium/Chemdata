"""Data analysis functionality for exported compound data.

This module provides functionality to:
1. Analyze compound properties and patterns
2. Generate statistical summaries
3. Identify trends and correlations
4. Detect outliers and anomalies
5. Generate insights and recommendations
"""

import logging
from typing import Dict, List, Optional, Any, Set
import pandas as pd
import numpy as np
from dataclasses import dataclass
from scipy import stats as scipy_stats
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

from ....models.validation import ValidationResult
from .data_enrichment import EnrichedData


@dataclass
class AnalysisResult(ValidationResult):
    """Result of data analysis."""
    
    property_stats: Dict[str, Dict[str, float]]
    correlations: Dict[str, Dict[str, float]]
    clusters: Dict[str, List[str]]
    outliers: Dict[str, List[str]]
    trends: Dict[str, Any]
    insights: List[str]


class DataAnalyzer:
    """Analyzer for exported compound data."""

    # Properties to analyze
    NUMERIC_PROPERTIES = {
        "molecular_weight",
        "logp",
        "psa",
        "hba",
        "hbd",
        "rotatable_bonds",
    }

    # Prediction types to analyze
    PREDICTION_TYPES = {
        "bbb",
        "activity",
        "toxicity",
        "abuse",
    }

    def __init__(
        self,
        log_level: int = logging.INFO,
    ):
        """Initialize data analyzer."""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)

    def analyze_compounds(
        self,
        compounds: List[EnrichedData],
        properties: Optional[Set[str]] = None,
        prediction_types: Optional[Set[str]] = None,
    ) -> AnalysisResult:
        """Analyze compound data."""
        self.logger.debug(f"Analyzing {len(compounds)} compounds")
        
        try:
            # Get properties to analyze
            properties = properties or self.NUMERIC_PROPERTIES
            prediction_types = prediction_types or self.PREDICTION_TYPES
            
            # Calculate statistics
            property_stats = self._calculate_property_stats(
                compounds, properties
            )
            
            # Calculate correlations
            correlations = self._calculate_correlations(
                compounds, properties, prediction_types
            )
            
            # Find clusters
            clusters = self._find_clusters(
                compounds, properties
            )
            
            # Detect outliers
            outliers = self._detect_outliers(
                compounds, properties
            )
            
            # Analyze trends
            trends = self._analyze_trends(
                compounds, properties, prediction_types
            )
            
            # Generate insights
            insights = self._generate_insights(
                property_stats,
                correlations,
                clusters,
                outliers,
                trends,
            )
            
            return AnalysisResult(
                is_valid=True,
                property_stats=property_stats,
                correlations=correlations,
                clusters=clusters,
                outliers=outliers,
                trends=trends,
                insights=insights,
                issues=[],
            )
            
        except Exception as e:
            self.logger.error(
                f"Error analyzing compounds: {str(e)}",
                exc_info=True
            )
            return AnalysisResult(
                is_valid=False,
                property_stats={},
                correlations={},
                clusters={},
                outliers={},
                trends={},
                insights=[],
                issues=[str(e)],
            )

    def _calculate_property_stats(
        self,
        compounds: List[EnrichedData],
        properties: Set[str],
    ) -> Dict[str, Dict[str, float]]:
        """Calculate statistics for numeric properties."""
        stats = {}
        
        for prop in properties:
            values = [
                float(c.properties.get(prop, 0))
                for c in compounds
                if c.properties.get(prop) is not None
            ]
            
            if values:
                stats[prop] = {
                    "mean": np.mean(values),
                    "std": np.std(values),
                    "min": np.min(values),
                    "max": np.max(values),
                    "median": np.median(values),
                    "q1": np.percentile(values, 25),
                    "q3": np.percentile(values, 75),
                    "skew": scipy_stats.skew(values),
                    "kurtosis": scipy_stats.kurtosis(values),
                }
        
        return stats

    def _calculate_correlations(
        self,
        compounds: List[EnrichedData],
        properties: Set[str],
        prediction_types: Set[str],
    ) -> Dict[str, Dict[str, float]]:
        """Calculate correlations between properties and predictions."""
        # Create DataFrame with properties
        data = {
            prop: [
                float(c.properties.get(prop, 0))
                for c in compounds
            ]
            for prop in properties
        }
        
        # Add prediction confidences
        for pred_type in prediction_types:
            data[f"{pred_type}_confidence"] = [
                float(getattr(c, f"{pred_type}_predictions", {})
                      .get("confidence", 0))
                for c in compounds
            ]
        
        df = pd.DataFrame(data)
        
        # Calculate correlations
        correlations = {}
        for col1 in df.columns:
            correlations[col1] = {}
            for col2 in df.columns:
                if col1 != col2:
                    correlations[col1][col2] = df[col1].corr(df[col2])
        
        return correlations

    def _find_clusters(
        self,
        compounds: List[EnrichedData],
        properties: Set[str],
    ) -> Dict[str, List[str]]:
        """Find clusters in property space."""
        # Create feature matrix
        features = []
        valid_compounds = []
        
        for compound in compounds:
            values = [
                compound.properties.get(prop)
                for prop in properties
            ]
            if all(v is not None for v in values):
                features.append(values)
                valid_compounds.append(compound)
        
        if not features:
            return {}
        
        # Scale features
        scaler = StandardScaler()
        X = scaler.fit_transform(features)
        
        # Find clusters
        clustering = DBSCAN(eps=0.5, min_samples=5)
        labels = clustering.fit_predict(X)
        
        # Group compounds by cluster
        clusters = {}
        for label, compound in zip(labels, valid_compounds):
            if label >= 0:  # Ignore noise points (-1)
                if label not in clusters:
                    clusters[str(label)] = []
                clusters[str(label)].append(compound.compound.name)
        
        return clusters

    def _detect_outliers(
        self,
        compounds: List[EnrichedData],
        properties: Set[str],
    ) -> Dict[str, List[str]]:
        """Detect outliers in properties."""
        outliers = {}
        
        for prop in properties:
            values = [
                (c.compound.name, float(c.properties.get(prop, 0)))
                for c in compounds
                if c.properties.get(prop) is not None
            ]
            
            if values:
                names, nums = zip(*values)
                
                # Calculate z-scores
                z_scores = np.abs(scipy_stats.zscore(nums))
                
                # Find outliers (z-score > 3)
                outlier_idx = np.where(z_scores > 3)[0]
                if len(outlier_idx) > 0:
                    outliers[prop] = [
                        names[i] for i in outlier_idx
                    ]
        
        return outliers

    def _analyze_trends(
        self,
        compounds: List[EnrichedData],
        properties: Set[str],
        prediction_types: Set[str],
    ) -> Dict[str, Any]:
        """Analyze trends in data."""
        trends = {
            "property_trends": self._analyze_property_trends(
                compounds, properties
            ),
            "prediction_trends": self._analyze_prediction_trends(
                compounds, prediction_types
            ),
            "temporal_trends": self._analyze_temporal_trends(
                compounds
            ),
        }
        return trends

    def _analyze_property_trends(
        self,
        compounds: List[EnrichedData],
        properties: Set[str],
    ) -> Dict[str, Any]:
        """Analyze trends in properties."""
        trends = {}
        
        for prop in properties:
            values = [
                float(c.properties.get(prop, 0))
                for c in compounds
                if c.properties.get(prop) is not None
            ]
            
            if values:
                # Check for normal distribution
                _, p_value = scipy_stats.normaltest(values)
                is_normal = p_value > 0.05
                
                # Find modality
                kernel = scipy_stats.gaussian_kde(values)
                x = np.linspace(min(values), max(values), 100)
                y = kernel(x)
                peaks = self._find_peaks(y)
                
                trends[prop] = {
                    "distribution": "normal" if is_normal else "non-normal",
                    "modality": len(peaks),
                    "trend": (
                        "increasing" if np.polyfit(
                            range(len(values)), values, 1
                        )[0] > 0
                        else "decreasing"
                    ),
                }
        
        return trends

    def _analyze_prediction_trends(
        self,
        compounds: List[EnrichedData],
        prediction_types: Set[str],
    ) -> Dict[str, Any]:
        """Analyze trends in predictions."""
        trends = {}
        
        for pred_type in prediction_types:
            confidences = [
                float(getattr(c, f"{pred_type}_predictions", {})
                      .get("confidence", 0))
                for c in compounds
            ]
            
            if confidences:
                trends[pred_type] = {
                    "mean_confidence": np.mean(confidences),
                    "confidence_trend": (
                        "improving" if np.polyfit(
                            range(len(confidences)), confidences, 1
                        )[0] > 0
                        else "declining"
                    ),
                }
        
        return trends

    def _analyze_temporal_trends(
        self,
        compounds: List[EnrichedData],
    ) -> Dict[str, Any]:
        """Analyze temporal trends."""
        # Get compounds with timestamps
        dated_compounds = [
            (c, max(c.timestamps.values()))
            for c in compounds
            if c.timestamps
        ]
        
        if not dated_compounds:
            return {}
        
        # Sort by date
        dated_compounds.sort(key=lambda x: x[1])
        compounds_sorted, dates = zip(*dated_compounds)
        
        # Analyze property changes over time
        trends = {}
        for prop in self.NUMERIC_PROPERTIES:
            values = [
                float(c.properties.get(prop, 0))
                for c in compounds_sorted
                if c.properties.get(prop) is not None
            ]
            
            if values:
                slope = np.polyfit(range(len(values)), values, 1)[0]
                trends[prop] = {
                    "trend": "increasing" if slope > 0 else "decreasing",
                    "slope": slope,
                }
        
        return trends

    def _find_peaks(self, y: np.ndarray) -> List[int]:
        """Find peaks in array."""
        peaks = []
        for i in range(1, len(y) - 1):
            if y[i - 1] < y[i] > y[i + 1]:
                peaks.append(i)
        return peaks

    def _generate_insights(
        self,
        property_stats: Dict[str, Dict[str, float]],
        correlations: Dict[str, Dict[str, float]],
        clusters: Dict[str, List[str]],
        outliers: Dict[str, List[str]],
        trends: Dict[str, Any],
    ) -> List[str]:
        """Generate insights from analysis results."""
        insights = []
        
        # Generate different types of insights
        insights.extend(self._generate_property_insights(property_stats))
        insights.extend(self._generate_correlation_insights(correlations))
        insights.extend(self._generate_cluster_insights(clusters))
        insights.extend(self._generate_outlier_insights(outliers))
        insights.extend(self._generate_trend_insights(trends))
        
        return insights

    def _generate_property_insights(
        self,
        property_stats: Dict[str, Dict[str, float]],
    ) -> List[str]:
        """Generate insights about property distributions."""
        insights = []
        
        for prop, stat_values in property_stats.items():
            if stat_values["skew"] > 1:
                insights.append(
                    f"Property {prop} shows significant right skew "
                    f"(skew={stat_values['skew']:.2f})"
                )
            elif stat_values["skew"] < -1:
                insights.append(
                    f"Property {prop} shows significant left skew "
                    f"(skew={stat_values['skew']:.2f})"
                )
        
        return insights

    def _generate_correlation_insights(
        self,
        correlations: Dict[str, Dict[str, float]],
    ) -> List[str]:
        """Generate insights about correlations."""
        insights = []
        strong_correlations = []
        
        # Find strong correlations
        for prop1, corrs in correlations.items():
            for prop2, corr in corrs.items():
                if abs(corr) > 0.7:
                    strong_correlations.append(
                        (prop1, prop2, corr)
                    )
        
        # Format correlation insights
        if strong_correlations:
            insights.append("Strong correlations found:")
            for prop1, prop2, corr in strong_correlations:
                insights.append(
                    f"  {prop1} vs {prop2}: {corr:.2f}"
                )
        
        return insights

    def _generate_cluster_insights(
        self,
        clusters: Dict[str, List[str]],
    ) -> List[str]:
        """Generate insights about clusters."""
        insights = []
        
        if clusters:
            insights.append(
                f"Found {len(clusters)} distinct compound clusters"
            )
            for cluster, compounds in clusters.items():
                insights.append(
                    f"  Cluster {cluster}: {len(compounds)} compounds"
                )
        
        return insights

    def _generate_outlier_insights(
        self,
        outliers: Dict[str, List[str]],
    ) -> List[str]:
        """Generate insights about outliers."""
        insights = []
        
        if outliers:
            insights.append("Outliers detected:")
            for prop, compounds in outliers.items():
                insights.append(
                    f"  {prop}: {len(compounds)} outliers"
                )
        
        return insights

    def _generate_trend_insights(
        self,
        trends: Dict[str, Any],
    ) -> List[str]:
        """Generate insights about trends."""
        insights = []
        
        if "property_trends" in trends:
            for prop, trend_data in trends["property_trends"].items():
                if trend_data["distribution"] == "non-normal":
                    insights.append(
                        f"Property {prop} shows non-normal distribution "
                        f"with {trend_data['modality']} modes"
                    )
        
        return insights
