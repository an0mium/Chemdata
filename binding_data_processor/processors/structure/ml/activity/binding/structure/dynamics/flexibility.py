"""Protein flexibility analysis functionality."""

import logging
from typing import Dict, List, Optional, Any
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

logger = logging.getLogger(__name__)


class FlexibilityAnalyzer:
    """Analyzes protein flexibility and dynamics."""

    def __init__(self):
        """Initialize flexibility analyzer."""
        self.logger = logging.getLogger(__name__)
        self.settings = {
            "window_size": 5,  # Residues, for smoothing
            "peak_prominence": 0.5,  # For identifying flexible regions
            "flexibility_threshold": 1.0,  # Standard deviations above mean
            "rigidity_threshold": -1.0,  # Standard deviations below mean
            "min_region_size": 3,  # Minimum residues for a flexible/rigid region
            "secondary_structure_weights": {
                "H": 0.3,  # Alpha helix - more rigid
                "G": 0.4,  # 3-10 helix
                "I": 0.3,  # Pi helix
                "E": 0.3,  # Beta strand - more rigid
                "B": 0.5,  # Beta bridge
                "T": 0.7,  # Turn - more flexible
                "S": 0.8,  # Bend - more flexible
                "-": 1.0,  # Coil - most flexible
            },
        }

    def analyze_flexibility(self, structure: Structure) -> Dict[str, Any]:
        """Analyze protein flexibility using multiple methods.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary of flexibility analysis results
        """
        try:
            # Get B-factors and calculate initial flexibility scores
            flexibility_scores = self._calculate_flexibility_scores(structure)
            if not flexibility_scores:
                return {}

            # Identify flexible and rigid regions
            regions = self._identify_flexibility_regions(flexibility_scores)

            # Calculate additional flexibility metrics
            metrics = self._calculate_flexibility_metrics(structure, flexibility_scores)

            # Analyze secondary structure contribution
            ss_analysis = self._analyze_secondary_structure(structure, flexibility_scores)

            # Analyze local flexibility patterns
            patterns = self._analyze_flexibility_patterns(flexibility_scores)

            # Combine results
            return {
                "residue_flexibility": flexibility_scores,
                "regions": regions,
                "metrics": metrics,
                "secondary_structure": ss_analysis,
                "patterns": patterns,
                "distribution": self._analyze_flexibility_distribution(list(flexibility_scores.values())),
                "correlations": self._analyze_flexibility_correlations(structure, flexibility_scores),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing flexibility: {str(e)}")
            return {}

    def analyze_site_flexibility(
        self,
        site_residues: List[int],
        structure: Structure,
        include_context: bool = True,
    ) -> Dict[str, Any]:
        """Analyze flexibility of specific binding site residues.

        Args:
            site_residues: List of residue numbers in site
            structure: Full structure for context
            include_context: Whether to analyze surrounding context

        Returns:
            Dictionary of site flexibility properties
        """
        try:
            # Get flexibility scores for all residues
            all_scores = self._calculate_flexibility_scores(structure)
            if not all_scores:
                return {}

            # Extract site scores
            site_scores = {res: all_scores[res] for res in site_residues if res in all_scores}
            if not site_scores:
                return {}

            # Calculate site statistics
            site_stats = {
                "mean_flexibility": float(np.mean(list(site_scores.values()))),
                "std_flexibility": float(np.std(list(site_scores.values()))),
                "min_flexibility": float(np.min(list(site_scores.values()))),
                "max_flexibility": float(np.max(list(site_scores.values()))),
            }

            # Calculate relative flexibility compared to whole structure
            all_mean = np.mean(list(all_scores.values()))
            all_std = np.std(list(all_scores.values()))
            site_stats["relative_flexibility"] = float((site_stats["mean_flexibility"] - all_mean) / all_std)

            # Analyze local context if requested
            context = {}
            if include_context:
                context = self._analyze_flexibility_context(site_residues, structure, all_scores)

            return {
                "scores": site_scores,
                "statistics": site_stats,
                "context": context,
                "patterns": self._analyze_flexibility_patterns(site_scores),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing site flexibility: {str(e)}")
            return {}

    def _calculate_flexibility_scores(self, structure: Structure) -> Dict[int, float]:
        """Calculate residue flexibility scores.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary mapping residue numbers to flexibility scores
        """
        try:
            scores = {}

            # Get B-factors
            bfactors = {}
            for residue in structure.get_residues():
                res_id = residue.get_id()[1]
                b_values = [atom.get_bfactor() for atom in residue]
                if b_values:
                    bfactors[res_id] = float(np.mean(b_values))

            if not bfactors:
                return {}

            # Normalize B-factors
            mean_b = np.mean(list(bfactors.values()))
            std_b = np.std(list(bfactors.values()))
            for res_id, value in bfactors.items():
                scores[res_id] = (value - mean_b) / std_b if std_b > 0 else 0.0

            # Apply secondary structure weighting
            try:
                dssp = dssp_dict_from_pdb_file(structure)
                for res_id in scores:
                    if res_id in dssp:
                        ss = dssp[res_id][2]
                        weight = self.settings["secondary_structure_weights"].get(ss, 1.0)
                        scores[res_id] *= weight
            except Exception:
                pass  # Continue without secondary structure weighting

            # Smooth scores
            smoothed = self._smooth_scores(scores)
            return smoothed

        except Exception as e:
            self.logger.error(f"Error calculating flexibility scores: {str(e)}")
            return {}

    def _smooth_scores(self, scores: Dict[int, float]) -> Dict[int, float]:
        """Apply window smoothing to flexibility scores.

        Args:
            scores: Dictionary mapping residue numbers to scores

        Returns:
            Dictionary of smoothed scores
        """
        try:
            smoothed = {}
            residues = sorted(scores.keys())
            window = self.settings["window_size"]

            for i, res_id in enumerate(residues):
                # Get window indices
                start = max(0, i - window // 2)
                end = min(len(residues), i + window // 2 + 1)
                window_ids = residues[start:end]

                # Calculate weighted average
                weights = np.exp(-0.5 * np.abs(np.array(window_ids) - res_id))
                values = [scores[rid] * w for rid, w in zip(window_ids, weights)]
                smoothed[res_id] = float(np.sum(values) / np.sum(weights))

            return smoothed

        except Exception as e:
            self.logger.error(f"Error smoothing scores: {str(e)}")
            return scores

    def _identify_flexibility_regions(self, scores: Dict[int, float]) -> Dict[str, List[List[int]]]:
        """Identify continuous regions of high/low flexibility.

        Args:
            scores: Dictionary mapping residue numbers to flexibility scores

        Returns:
            Dictionary mapping region types to lists of residue ranges
        """
        try:
            regions = {
                "flexible": [],
                "rigid": [],
                "moderate": [],
            }

            # Calculate thresholds
            values = np.array(list(scores.values()))
            mean = np.mean(values)
            std = np.std(values)
            flex_thresh = mean + self.settings["flexibility_threshold"] * std
            rigid_thresh = mean + self.settings["rigidity_threshold"] * std

            # Find continuous regions
            residues = sorted(scores.keys())
            current_region = []
            current_type = None

            for res_id in residues:
                score = scores[res_id]

                # Determine region type
                if score > flex_thresh:
                    region_type = "flexible"
                elif score < rigid_thresh:
                    region_type = "rigid"
                else:
                    region_type = "moderate"

                # Handle region transitions
                if region_type != current_type:
                    if current_region and len(current_region) >= self.settings["min_region_size"]:
                        regions[current_type].append(current_region)
                    current_region = [res_id]
                    current_type = region_type
                else:
                    current_region.append(res_id)

            # Add final region
            if current_region and len(current_region) >= self.settings["min_region_size"]:
                regions[current_type].append(current_region)

            return regions

        except Exception as e:
            self.logger.error(f"Error identifying flexibility regions: {str(e)}")
            return {"flexible": [], "rigid": [], "moderate": []}

    def _calculate_flexibility_metrics(
        self,
        structure: Structure,
        scores: Dict[int, float],
    ) -> Dict[str, float]:
        """Calculate additional flexibility metrics.

        Args:
            structure: BioPython Structure object
            scores: Dictionary mapping residue numbers to flexibility scores

        Returns:
            Dictionary of flexibility metrics
        """
        try:
            values = np.array(list(scores.values()))

            metrics = {
                "mean_flexibility": float(np.mean(values)),
                "std_flexibility": float(np.std(values)),
                "skewness": float(self._calculate_skewness(values)),
                "kurtosis": float(self._calculate_kurtosis(values)),
                "flexibility_index": float(np.sum(np.abs(values)) / len(values)),
            }

            # Calculate domain-level metrics
            try:
                domains = self._identify_domains(structure)
                if domains:
                    domain_flex = []
                    for domain in domains:
                        domain_scores = [scores[res_id] for res_id in domain if res_id in scores]
                        if domain_scores:
                            domain_flex.append(float(np.mean(domain_scores)))

                    metrics["domain_flexibility_variation"] = float(np.std(domain_flex))
                    metrics["max_domain_flexibility"] = float(np.max(domain_flex))
                    metrics["min_domain_flexibility"] = float(np.min(domain_flex))
            except Exception:
                pass

            return metrics

        except Exception as e:
            self.logger.error(f"Error calculating flexibility metrics: {str(e)}")
            return {}

    def _analyze_secondary_structure(
        self,
        structure: Structure,
        scores: Dict[int, float],
    ) -> Dict[str, Any]:
        """Analyze relationship between flexibility and secondary structure.

        Args:
            structure: BioPython Structure object
            scores: Dictionary mapping residue numbers to flexibility scores

        Returns:
            Dictionary of secondary structure analysis
        """
        try:
            ss_scores = {ss: [] for ss in self.settings["secondary_structure_weights"]}

            # Get DSSP assignments
            dssp = dssp_dict_from_pdb_file(structure)

            # Group flexibility scores by secondary structure
            for res_id, score in scores.items():
                if res_id in dssp:
                    ss = dssp[res_id][2]
                    ss_scores[ss].append(score)

            # Calculate statistics for each type
            analysis = {}
            for ss, values in ss_scores.items():
                if values:
                    analysis[ss] = {
                        "mean": float(np.mean(values)),
                        "std": float(np.std(values)),
                        "count": len(values),
                    }

            return analysis

        except Exception as e:
            self.logger.error(f"Error analyzing secondary structure: {str(e)}")
            return {}

    def _analyze_flexibility_patterns(self, scores: Dict[int, float]) -> Dict[str, Any]:
        """Analyze patterns in flexibility distribution.

        Args:
            scores: Dictionary mapping residue numbers to flexibility scores

        Returns:
            Dictionary of pattern analysis results
        """
        try:
            values = np.array(list(scores.values()))
            residues = np.array(list(scores.keys()))

            # Find peaks in flexibility profile
            peaks, properties = find_peaks(
                values,
                prominence=self.settings["peak_prominence"],
                width=self.settings["window_size"],
            )

            # Analyze periodicity
            try:
                from scipy import fft

                frequencies = fft.fftfreq(len(values))
                spectrum = np.abs(fft.fft(values))
                main_freq_idx = np.argmax(spectrum[1:]) + 1
                periodicity = 1.0 / frequencies[main_freq_idx] if frequencies[main_freq_idx] != 0 else 0.0
            except Exception:
                periodicity = 0.0

            return {
                "peaks": {
                    "residues": residues[peaks].tolist(),
                    "scores": values[peaks].tolist(),
                    "prominences": properties["prominences"].tolist(),
                    "widths": properties["widths"].tolist(),
                },
                "periodicity": float(periodicity),
                "autocorrelation": self._calculate_autocorrelation(values),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing flexibility patterns: {str(e)}")
            return {}

    def _analyze_flexibility_distribution(self, values: List[float]) -> Dict[str, Any]:
        """Analyze statistical distribution of flexibility scores.

        Args:
            values: List of flexibility scores

        Returns:
            Dictionary of distribution properties
        """
        try:
            if len(values) < 3:
                return {}

            # Calculate basic statistics
            stats = {
                "skewness": float(self._calculate_skewness(values)),
                "kurtosis": float(self._calculate_kurtosis(values)),
            }

            # Estimate probability density
            try:
                kde = gaussian_kde(values)
                x = np.linspace(min(values), max(values), 100)
                density = kde(x)
                stats["density"] = {
                    "x": x.tolist(),
                    "y": density.tolist(),
                }
            except Exception:
                pass

            return stats

        except Exception as e:
            self.logger.error(f"Error analyzing flexibility distribution: {str(e)}")
            return {}

    def _analyze_flexibility_correlations(
        self,
        structure: Structure,
        scores: Dict[int, float],
    ) -> Dict[str, float]:
        """Analyze correlations between flexibility and structural properties.

        Args:
            structure: BioPython Structure object
            scores: Dictionary mapping residue numbers to flexibility scores

        Returns:
            Dictionary of correlation coefficients
        """
        try:
            correlations = {}

            # Correlation with depth
            try:
                depths = self._calculate_residue_depths(structure)
                common_residues = set(scores.keys()) & set(depths.keys())
                if len(common_residues) > 2:
                    x = [scores[res] for res in common_residues]
                    y = [depths[res] for res in common_residues]
                    correlations["depth"] = float(np.corrcoef(x, y)[0, 1])
            except Exception:
                pass

            # Correlation with contacts
            try:
                contacts = self._count_residue_contacts(structure)
                common_residues = set(scores.keys()) & set(contacts.keys())
                if len(common_residues) > 2:
                    x = [scores[res] for res in common_residues]
                    y = [contacts[res] for res in common_residues]
                    correlations["contacts"] = float(np.corrcoef(x, y)[0, 1])
            except Exception:
                pass

            return correlations

        except Exception as e:
            self.logger.error(f"Error analyzing flexibility correlations: {str(e)}")
            return {}

    def _analyze_flexibility_context(
        self,
        site_residues: List[int],
        structure: Structure,
        scores: Dict[int, float],
    ) -> Dict[str, Any]:
        """Analyze flexibility context around binding site.

        Args:
            site_residues: List of residue numbers in site
            structure: Full structure
            scores: Pre-calculated flexibility scores

        Returns:
            Dictionary of context analysis results
        """
        try:
            # Get surrounding residues
            context_residues = set()
            for res_id in site_residues:
                neighbors = self._get_residue_neighbors(structure, res_id)
                context_residues.update(neighbors)
            context_residues = context_residues - set(site_residues)

            # Calculate context statistics
            if context_residues:
                context_scores = [scores[res] for res in context_residues if res in scores]
                if context_scores:
                    return {
                        "mean_context_flexibility": float(np.mean(context_scores)),
                        "context_residues": sorted(list(context_residues)),
                        "flexibility_gradient": float(np.mean([scores[res] for res in site_residues if res in scores]) - np.mean(context_scores)),
                    }

            return {}

        except Exception as e:
            self.logger.error(f"Error analyzing flexibility context: {str(e)}")
            return {}

    def _calculate_skewness(self, values: np.ndarray) -> float:
        """Calculate distribution skewness."""
        try:
            mean = np.mean(values)
            std = np.std(values)
            if std == 0:
                return 0.0
            return float(np.mean(((values - mean) / std) ** 3))
        except Exception:
            return 0.0

    def _calculate_kurtosis(self, values: np.ndarray) -> float:
        """Calculate distribution kurtosis."""
        try:
            mean = np.mean(values)
            std = np.std(values)
            if std == 0:
                return 0.0
            return float(np.mean(((values - mean) / std) ** 4) - 3)
        except Exception:
            return 0.0

    def _calculate_autocorrelation(self, values: np.ndarray) -> List[float]:
        """Calculate autocorrelation of flexibility profile."""
        try:
            from scipy.signal import correlate

            autocorr = correlate(values, values, mode="full")
            center = len(autocorr) // 2
            return autocorr[center:].tolist()
        except Exception:
            return []

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

    def _get_residue_neighbors(self, structure: Structure, res_id: int) -> List[int]:
        """Get neighboring residues within distance cutoff."""
        try:
            neighbors = set()
            cutoff = 8.0  # Angstroms

            target_res = None
            for residue in structure.get_residues():
                if residue.get_id()[1] == res_id:
                    target_res = residue
                    break

            if target_res is None:
                return []

            for residue in structure.get_residues():
                if residue != target_res:
                    # Calculate minimum distance between residues
                    min_dist = float("inf")
                    for atom1 in target_res:
                        for atom2 in residue:
                            dist = atom1 - atom2
                            min_dist = min(min_dist, dist)

                    if min_dist < cutoff:
                        neighbors.add(residue.get_id()[1])

            return sorted(list(neighbors))

        except Exception as e:
            self.logger.error(f"Error getting residue neighbors: {str(e)}")
            return []

    def _identify_domains(self, structure: Structure) -> List[List[int]]:
        """Identify protein domains using contact topology."""
        try:
            # Build contact network
            import networkx as nx

            G = nx.Graph()
            cutoff = 8.0  # Angstroms

            for res1 in structure.get_residues():
                res1_id = res1.get_id()[1]
                G.add_node(res1_id)

                for res2 in structure.get_residues():
                    res2_id = res2.get_id()[1]
                    if res1_id != res2_id:
                        # Check contact
                        min_dist = float("inf")
                        for atom1 in res1:
                            for atom2 in res2:
                                dist = atom1 - atom2
                                min_dist = min(min_dist, dist)

                        if min_dist < cutoff:
                            G.add_edge(res1_id, res2_id)

            # Use community detection to identify domains
            communities = list(nx.community.louvain_communities(G))
            return [sorted(comm) for comm in communities]

        except Exception as e:
            self.logger.error(f"Error identifying domains: {str(e)}")
            return []
