"""Feature validation and normalization utilities."""

import logging
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_classif

from .base import MLProcessor


class FeatureProcessor(MLProcessor):
    """Process molecular features for ML models."""

    def __init__(self):
        """Initialize feature processor."""
        super().__init__()
        self.scalers = {}
        self.selectors = {}
        self.feature_names = {}

    def validate_features(
        self,
        features: Dict[str, Union[np.ndarray, Dict[str, float]]],
        required_types: Optional[List[str]] = None,
    ) -> Tuple[bool, List[str]]:
        """
        Validate molecular features.

        Args:
            features: Dictionary of features
            required_types: List of required feature types

        Returns:
            Tuple of (is_valid, error_messages)
        """
        try:
            errors = []

            # Check required types
            if required_types:
                missing = [t for t in required_types if t not in features]
                if missing:
                    errors.append(f"Missing required feature types: {missing}")

            # Validate each feature type
            for feat_type, feat_data in features.items():
                # Check for None values
                if feat_data is None:
                    errors.append(f"Feature type {feat_type} has None value")
                    continue

                # Validate numpy arrays
                if isinstance(feat_data, np.ndarray):
                    if not feat_data.size:
                        errors.append(f"Empty array for feature type {feat_type}")
                    if not np.isfinite(feat_data).all():
                        errors.append(f"Non-finite values in feature type {feat_type}")

                # Validate dictionaries (e.g. 3D features)
                elif isinstance(feat_data, dict):
                    if not feat_data:
                        errors.append(f"Empty dictionary for feature type {feat_type}")
                    for key, value in feat_data.items():
                        if not isinstance(value, (int, float)):
                            errors.append(
                                f"Invalid value type for {feat_type}.{key}: {type(value)}"
                            )
                        elif not np.isfinite(value):
                            errors.append(
                                f"Non-finite value for {feat_type}.{key}: {value}"
                            )

                else:
                    errors.append(
                        f"Invalid data type for feature {feat_type}: {type(feat_data)}"
                    )

            return len(errors) == 0, errors

        except Exception as e:
            self.logger.error(f"Error validating features: {str(e)}")
            return False, [str(e)]

    def normalize_features(
        self,
        features: Dict[str, Union[np.ndarray, Dict[str, float]]],
        method: str = "standard",
        fit: bool = True,
    ) -> Dict[str, Union[np.ndarray, Dict[str, float]]]:
        """
        Normalize molecular features.

        Args:
            features: Dictionary of features
            method: Normalization method ('standard', 'minmax', or 'robust')
            fit: Whether to fit scalers on data

        Returns:
            Dictionary of normalized features
        """
        try:
            normalized = {}

            for feat_type, feat_data in features.items():
                if feat_data is None:
                    continue

                # Handle numpy arrays
                if isinstance(feat_data, np.ndarray):
                    if fit or feat_type not in self.scalers:
                        if method == "standard":
                            scaler = StandardScaler()
                        elif method == "minmax":
                            scaler = MinMaxScaler()
                        elif method == "robust":
                            scaler = RobustScaler()
                        else:
                            raise ValueError(f"Unknown normalization method: {method}")

                        # Reshape for 1D arrays
                        reshaped = (
                            feat_data.reshape(-1, 1)
                            if feat_data.ndim == 1
                            else feat_data
                        )
                        normalized[feat_type] = scaler.fit_transform(reshaped)
                        self.scalers[feat_type] = scaler
                    else:
                        reshaped = (
                            feat_data.reshape(-1, 1)
                            if feat_data.ndim == 1
                            else feat_data
                        )
                        normalized[feat_type] = self.scalers[feat_type].transform(
                            reshaped
                        )

                # Handle dictionaries
                elif isinstance(feat_data, dict):
                    if fit or feat_type not in self.scalers:
                        values = np.array(list(feat_data.values())).reshape(-1, 1)
                        if method == "standard":
                            scaler = StandardScaler()
                        elif method == "minmax":
                            scaler = MinMaxScaler()
                        elif method == "robust":
                            scaler = RobustScaler()
                        else:
                            raise ValueError(f"Unknown normalization method: {method}")

                        normalized_values = scaler.fit_transform(values).flatten()
                        self.scalers[feat_type] = scaler
                    else:
                        values = np.array(list(feat_data.values())).reshape(-1, 1)
                        normalized_values = (
                            self.scalers[feat_type].transform(values).flatten()
                        )

                    normalized[feat_type] = dict(
                        zip(feat_data.keys(), normalized_values)
                    )

            return normalized

        except Exception as e:
            self.logger.error(f"Error normalizing features: {str(e)}")
            return features

    def select_features(
        self,
        features: Dict[str, Union[np.ndarray, Dict[str, float]]],
        labels: np.ndarray,
        method: str = "variance",
        threshold: float = 0.0,
        k: Optional[int] = None,
        fit: bool = True,
    ) -> Dict[str, Union[np.ndarray, Dict[str, float]]]:
        """
        Select relevant features.

        Args:
            features: Dictionary of features
            labels: Target labels for supervised selection
            method: Selection method ('variance' or 'kbest')
            threshold: Variance threshold for 'variance' method
            k: Number of features to select for 'kbest' method
            fit: Whether to fit selectors on data

        Returns:
            Dictionary of selected features
        """
        try:
            selected = {}

            for feat_type, feat_data in features.items():
                if feat_data is None:
                    continue

                # Handle numpy arrays
                if isinstance(feat_data, np.ndarray):
                    if fit or feat_type not in self.selectors:
                        if method == "variance":
                            selector = VarianceThreshold(threshold=threshold)
                        elif method == "kbest":
                            if k is None:
                                k = feat_data.shape[1] // 2
                            selector = SelectKBest(score_func=f_classif, k=k)
                        else:
                            raise ValueError(f"Unknown selection method: {method}")

                        # Reshape for 1D arrays
                        reshaped = (
                            feat_data.reshape(-1, 1)
                            if feat_data.ndim == 1
                            else feat_data
                        )
                        selected[feat_type] = selector.fit_transform(reshaped, labels)
                        self.selectors[feat_type] = selector

                        # Store feature names/indices
                        if hasattr(selector, "get_support"):
                            self.feature_names[feat_type] = np.where(
                                selector.get_support()
                            )[0]
                    else:
                        reshaped = (
                            feat_data.reshape(-1, 1)
                            if feat_data.ndim == 1
                            else feat_data
                        )
                        selected[feat_type] = self.selectors[feat_type].transform(
                            reshaped
                        )

                # Handle dictionaries
                elif isinstance(feat_data, dict):
                    if fit or feat_type not in self.selectors:
                        values = np.array(list(feat_data.values())).reshape(-1, 1)
                        if method == "variance":
                            selector = VarianceThreshold(threshold=threshold)
                        elif method == "kbest":
                            if k is None:
                                k = len(feat_data) // 2
                            selector = SelectKBest(score_func=f_classif, k=k)
                        else:
                            raise ValueError(f"Unknown selection method: {method}")

                        selected_values = selector.fit_transform(
                            values, labels
                        ).flatten()
                        self.selectors[feat_type] = selector

                        # Store selected feature names
                        if hasattr(selector, "get_support"):
                            selected_indices = np.where(selector.get_support())[0]
                            self.feature_names[feat_type] = list(feat_data.keys())[
                                selected_indices
                            ]
                            selected[feat_type] = {
                                name: selected_values[i]
                                for i, name in enumerate(self.feature_names[feat_type])
                            }
                    else:
                        values = np.array(list(feat_data.values())).reshape(-1, 1)
                        selected_values = (
                            self.selectors[feat_type].transform(values).flatten()
                        )
                        selected[feat_type] = {
                            name: selected_values[i]
                            for i, name in enumerate(self.feature_names[feat_type])
                        }

            return selected

        except Exception as e:
            self.logger.error(f"Error selecting features: {str(e)}")
            return features

    def get_feature_names(self, feat_type: str) -> Optional[List[str]]:
        """Get names of selected features for a feature type."""
        return self.feature_names.get(feat_type)

    def reset(self):
        """Reset all scalers and selectors."""
        self.scalers = {}
        self.selectors = {}
        self.feature_names = {}
