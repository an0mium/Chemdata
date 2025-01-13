"""B-factor analysis functionality."""

import logging
from typing import Dict, List, Any
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue

logger = logging.getLogger(__name__)


class BFactorAnalyzer:
    """Analyzes protein B-factors."""

    def __init__(self):
        """Initialize B-factor analyzer."""
        self.logger = logging.getLogger(__name__)

    def analyze_bfactors(self, structure: Structure) -> Dict[str, Any]:
        """Analyze B-factors across structure.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of B-factor analysis
        """
        try:
            # Collect B-factors
            bfactors = {
                "backbone": [],
                "sidechain": [],
                "all": [],
            }

            residue_bfactors = {}

            for residue in structure.get_residues():
                res_bfactors = []
                for atom in residue:
                    bfactor = atom.get_bfactor()
                    bfactors["all"].append(bfactor)
                    res_bfactors.append(bfactor)
                    if atom.get_name() in ["N", "CA", "C", "O"]:
                        bfactors["backbone"].append(bfactor)
                    else:
                        bfactors["sidechain"].append(bfactor)

                # Store average B-factor for residue
                if res_bfactors:
                    residue_bfactors[residue.get_id()[1]] = float(np.mean(res_bfactors))

            # Calculate statistics
            stats = {}
            for region, values in bfactors.items():
                if values:
                    stats[region] = {
                        "mean": float(np.mean(values)),
                        "std": float(np.std(values)),
                        "min": float(np.min(values)),
                        "max": float(np.max(values)),
                        "median": float(np.median(values)),
                        "quartiles": [
                            float(np.percentile(values, 25)),
                            float(np.percentile(values, 75)),
                        ],
                    }

            return {
                "bfactors": residue_bfactors,
                "statistics": stats,
                "mobility_profile": self._calculate_mobility_profile(residue_bfactors),
                "regions": self._identify_mobility_regions(residue_bfactors),
                "correlations": self._analyze_bfactor_correlations(structure, residue_bfactors),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing B-factors: {str(e)}")
            return {}

    def get_residue_bfactors(self, residues: List[Residue]) -> Dict[str, Any]:
        """Get B-factors for specific residues.

        Args:
            residues: List of residues to analyze

        Returns:
            Dictionary of B-factor values and statistics
        """
        try:
            values = []
            residue_values = {}
            atom_values = {
                "backbone": [],
                "sidechain": [],
            }

            for residue in residues:
                res_bfactors = []
                for atom in residue:
                    bfactor = atom.get_bfactor()
                    res_bfactors.append(bfactor)
                    if atom.get_name() in ["N", "CA", "C", "O"]:
                        atom_values["backbone"].append(bfactor)
                    else:
                        atom_values["sidechain"].append(bfactor)

                if res_bfactors:
                    avg_bfactor = float(np.mean(res_bfactors))
                    values.append(avg_bfactor)
                    residue_values[residue.get_id()[1]] = avg_bfactor

            if not values:
                return {"mean": 0.0, "std": 0.0, "values": {}}

            # Calculate comprehensive statistics
            stats = {
                "mean": float(np.mean(values)),
                "std": float(np.std(values)),
                "median": float(np.median(values)),
                "min": float(np.min(values)),
                "max": float(np.max(values)),
                "quartiles": [
                    float(np.percentile(values, 25)),
                    float(np.percentile(values, 75)),
                ],
                "backbone_mean": float(np.mean(atom_values["backbone"])) if atom_values["backbone"] else 0.0,
                "sidechain_mean": float(np.mean(atom_values["sidechain"])) if atom_values["sidechain"] else 0.0,
            }

            return {
                **stats,
                "values": residue_values,
                "distribution": self._analyze_bfactor_distribution(values),
            }

        except Exception as e:
            self.logger.error(f"Error getting residue B-factors: {str(e)}")
            return {"mean": 0.0, "std": 0.0, "values": {}}

    def calculate_relative_bfactor(
        self,
        site_bfactors: Dict[int, float],
        structure: Structure,
    ) -> float:
        """Calculate relative B-factor compared to whole structure.

        Args:
            site_bfactors: Dictionary mapping residue IDs to B-factors
            structure: Full structure

        Returns:
            Relative B-factor score
        """
        try:
            # Get all B-factors
            all_bfactors = []
            for residue in structure.get_residues():
                for atom in residue:
                    all_bfactors.append(atom.get_bfactor())

            # Calculate Z-score
            site_mean = np.mean(list(site_bfactors.values()))
            all_mean = np.mean(all_bfactors)
            all_std = np.std(all_bfactors)

            if all_std == 0:
                return 0.0

            return float((site_mean - all_mean) / all_std)

        except Exception as e:
            self.logger.error(f"Error calculating relative B-factor: {str(e)}")
            return 0.0

    def _calculate_mobility_profile(
        self,
        bfactors: Dict[int, float],
    ) -> Dict[str, List[int]]:
        """Calculate mobility profile from B-factors.

        Args:
            bfactors: Dictionary mapping residue IDs to B-factors

        Returns:
            Dictionary categorizing residues by mobility
        """
        try:
            values = list(bfactors.values())
            if not values:
                return {"rigid": [], "moderate": [], "flexible": [], "highly_flexible": []}

            mean = np.mean(values)
            std = np.std(values)

            profile = {
                "rigid": [],  # < mean - std
                "moderate": [],  # between mean ± std
                "flexible": [],  # > mean + std
                "highly_flexible": [],  # > mean + 2*std
            }

            for res_id, bfactor in bfactors.items():
                if bfactor > mean + 2 * std:
                    profile["highly_flexible"].append(res_id)
                elif bfactor > mean + std:
                    profile["flexible"].append(res_id)
                elif bfactor < mean - std:
                    profile["rigid"].append(res_id)
                else:
                    profile["moderate"].append(res_id)

            return profile

        except Exception as e:
            self.logger.error(f"Error calculating mobility profile: {str(e)}")
            return {"rigid": [], "moderate": [], "flexible": [], "highly_flexible": []}

    def _identify_mobility_regions(
        self,
        bfactors: Dict[int, float],
        window_size: int = 5,
    ) -> Dict[str, List[List[int]]]:
        """Identify continuous regions of similar mobility.

        Args:
            bfactors: Dictionary mapping residue IDs to B-factors
            window_size: Window size for smoothing

        Returns:
            Dictionary mapping mobility types to lists of residue ranges
        """
        try:
            # Smooth B-factors using sliding window
            residues = sorted(bfactors.keys())
            smoothed = {}

            for i, res_id in enumerate(residues):
                start = max(0, i - window_size // 2)
                end = min(len(residues), i + window_size // 2 + 1)
                window = [bfactors[residues[j]] for j in range(start, end)]
                smoothed[res_id] = float(np.mean(window))

            # Calculate thresholds
            values = list(smoothed.values())
            mean = np.mean(values)
            std = np.std(values)

            # Identify regions
            regions = {
                "rigid": [],
                "moderate": [],
                "flexible": [],
                "highly_flexible": [],
            }

            current_type = None
            current_region = []

            for res_id in residues:
                bfactor = smoothed[res_id]

                # Determine mobility type
                if bfactor > mean + 2 * std:
                    mobility = "highly_flexible"
                elif bfactor > mean + std:
                    mobility = "flexible"
                elif bfactor < mean - std:
                    mobility = "rigid"
                else:
                    mobility = "moderate"

                # Handle region transitions
                if mobility != current_type:
                    if current_region:
                        regions[current_type].append(current_region)
                    current_type = mobility
                    current_region = [res_id]
                else:
                    current_region.append(res_id)

            # Add final region
            if current_region:
                regions[current_type].append(current_region)

            return regions

        except Exception as e:
            self.logger.error(f"Error identifying mobility regions: {str(e)}")
            return {"rigid": [], "moderate": [], "flexible": [], "highly_flexible": []}

    def _analyze_bfactor_correlations(
        self,
        structure: Structure,
        bfactors: Dict[int, float],
    ) -> Dict[str, float]:
        """Analyze correlations between B-factors and structural properties.

        Args:
            structure: BioPython Structure object
            bfactors: Dictionary mapping residue IDs to B-factors

        Returns:
            Dictionary of correlation coefficients
        """
        try:
            correlations = {}

            # Get structural properties
            properties = {
                "depth": self._calculate_residue_depths(structure),
                "contacts": self._count_residue_contacts(structure),
                "sasa": self._calculate_residue_sasa(structure),
            }

            # Calculate correlations
            for prop_name, prop_values in properties.items():
                common_residues = set(bfactors.keys()) & set(prop_values.keys())
                if len(common_residues) > 2:
                    x = [bfactors[res] for res in common_residues]
                    y = [prop_values[res] for res in common_residues]
                    correlation = float(np.corrcoef(x, y)[0, 1])
                    correlations[f"{prop_name}_correlation"] = correlation

            return correlations

        except Exception as e:
            self.logger.error(f"Error analyzing B-factor correlations: {str(e)}")
            return {}

    def _analyze_bfactor_distribution(self, values: List[float]) -> Dict[str, float]:
        """Analyze statistical distribution of B-factors.

        Args:
            values: List of B-factor values

        Returns:
            Dictionary of distribution properties
        """
        try:
            if len(values) < 3:
                return {}

            return {
                "skewness": float(self._calculate_skewness(values)),
                "kurtosis": float(self._calculate_kurtosis(values)),
                "normality": float(self._test_normality(values)),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing B-factor distribution: {str(e)}")
            return {}

    def _calculate_skewness(self, values: List[float]) -> float:
        """Calculate distribution skewness."""
        try:
            mean = np.mean(values)
            std = np.std(values)
            if std == 0:
                return 0.0
            return float(np.mean(((values - mean) / std) ** 3))
        except Exception:
            return 0.0

    def _calculate_kurtosis(self, values: List[float]) -> float:
        """Calculate distribution kurtosis."""
        try:
            mean = np.mean(values)
            std = np.std(values)
            if std == 0:
                return 0.0
            return float(np.mean(((values - mean) / std) ** 4) - 3)
        except Exception:
            return 0.0

    def _test_normality(self, values: List[float]) -> float:
        """Test for normal distribution using Shapiro-Wilk test."""
        try:
            from scipy import stats

            statistic, _ = stats.shapiro(values)
            return float(statistic)
        except Exception:
            return 0.0

    def _calculate_residue_depths(self, structure: Structure) -> Dict[int, float]:
        """Calculate residue depths from protein surface."""
        try:
            from Bio.PDB.ResidueDepth import get_surface, residue_depth

            surface = get_surface(structure)
            depths = {}

            for residue in structure.get_residues():
                depth = residue_depth(residue, surface)
                depths[residue.get_id()[1]] = float(depth)

            return depths

        except Exception as e:
            self.logger.error(f"Error calculating residue depths: {str(e)}")
            return {}

    def _count_residue_contacts(self, structure: Structure) -> Dict[int, int]:
        """Count contacts for each residue."""
        try:
            contacts = {}
            cutoff = 8.0  # Angstroms

            for res1 in structure.get_residues():
                res1_id = res1.get_id()[1]
                contacts[res1_id] = 0

                for res2 in structure.get_residues():
                    if res1 != res2:
                        # Calculate minimum distance between residues
                        min_dist = float("inf")
                        for atom1 in res1:
                            for atom2 in res2:
                                dist = atom1 - atom2
                                min_dist = min(min_dist, dist)

                        if min_dist < cutoff:
                            contacts[res1_id] += 1

            return contacts

        except Exception as e:
            self.logger.error(f"Error counting residue contacts: {str(e)}")
            return {}

    def _calculate_residue_sasa(self, structure: Structure) -> Dict[int, float]:
        """Calculate solvent accessible surface area for each residue."""
        try:
            from Bio.PDB.SASA import calculate_sasa

            sasa = {}
            for residue in structure.get_residues():
                area = calculate_sasa(residue)
                sasa[residue.get_id()[1]] = float(area)

            return sasa

        except Exception as e:
            self.logger.error(f"Error calculating residue SASA: {str(e)}")
            return {}
