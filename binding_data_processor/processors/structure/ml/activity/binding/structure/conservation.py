"""Conservation analysis functionality for protein structures."""

import logging
from typing import Dict, List, Optional, Any, Tuple
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Polypeptide import protein_letters_3to1, is_aa
from Bio.Align import substitution_matrices
import requests

logger = logging.getLogger(__name__)


class ConservationAnalyzer:
    """Analyzes sequence conservation in protein structures."""

    def __init__(self):
        """Initialize conservation analyzer."""
        self.logger = logging.getLogger(__name__)
        self.blosum = substitution_matrices.load("BLOSUM62")
        self.conservation_cache = {}

    def analyze_conservation(
        self,
        structure: Structure,
        window_size: int = 5,
        use_pssm: bool = True,
        use_entropy: bool = True,
        use_alphafold: bool = True,
    ) -> Dict[str, Any]:
        """Analyze sequence conservation using multiple methods.

        Args:
            structure: BioPython Structure object
            window_size: Window size for sliding window analysis
            use_pssm: Whether to use position-specific scoring matrices
            use_entropy: Whether to calculate sequence entropy
            use_alphafold: Whether to incorporate AlphaFold confidence scores

        Returns:
            Dictionary of conservation scores and analysis
        """
        try:
            # Get sequence and residue IDs
            sequence = ""
            residue_ids = []
            for residue in structure.get_residues():
                if not is_aa(residue):
                    continue
                sequence += protein_letters_3to1[residue.get_resname()]
                residue_ids.append(residue.get_id()[1])

            # Calculate basic conservation scores
            conservation = {
                "sequence": sequence,
                "residue_ids": residue_ids,
                "window_scores": self._calculate_window_conservation(sequence, window_size),
                "relative_scores": self._calculate_relative_conservation(sequence),
            }

            # Optional PSSM analysis
            if use_pssm:
                conservation["pssm_scores"] = self._calculate_pssm_conservation(sequence)

            # Optional entropy analysis
            if use_entropy:
                conservation["entropy_scores"] = self._calculate_sequence_entropy(sequence, window_size)

            # Optional AlphaFold confidence scores
            if use_alphafold:
                conservation["alphafold_scores"] = self._get_alphafold_confidence(structure)

            # Calculate overall conservation metrics
            conservation["metrics"] = self._calculate_conservation_metrics(conservation)

            return conservation

        except Exception as e:
            self.logger.error(f"Error analyzing conservation: {str(e)}")
            return {}

    def analyze_site_conservation(
        self,
        residues: List[Residue],
        window_size: int = 5,
    ) -> Dict[str, float]:
        """Analyze conservation of binding site residues.

        Args:
            residues: List of residues in binding site
            window_size: Window size for sliding window analysis

        Returns:
            Dictionary of conservation metrics
        """
        try:
            # Get site sequence
            sequence = ""
            for res in residues:
                if not is_aa(res):
                    continue
                sequence += protein_letters_3to1[res.get_resname()]

            # Calculate conservation scores
            window_scores = self._calculate_window_conservation(sequence, window_size)
            entropy_scores = self._calculate_sequence_entropy(sequence, window_size)
            relative_scores = self._calculate_relative_conservation(sequence)

            # Calculate site metrics
            metrics = {
                "average_conservation": float(np.mean(window_scores)),
                "min_conservation": float(np.min(window_scores)),
                "max_conservation": float(np.max(window_scores)),
                "entropy": float(np.mean(entropy_scores)),
                "relative_conservation": float(np.mean(relative_scores)),
                "variability": float(np.std(window_scores)),
                "percentile": self._calculate_percentile(window_scores),
            }

            return metrics

        except Exception as e:
            self.logger.error(f"Error analyzing site conservation: {str(e)}")
            return {}

    def get_residue_conservation(self, residue: Residue) -> float:
        """Get conservation score for single residue.

        Args:
            residue: BioPython Residue object

        Returns:
            Conservation score (0-1)
        """
        try:
            if not is_aa(residue):
                return 0.0

            # Get residue type
            res_type = protein_letters_3to1[residue.get_resname()]

            # Get BLOSUM scores for this residue type
            scores = []
            for aa in self.blosum:
                scores.append(self.blosum[res_type][aa])

            # Normalize score
            score = np.mean(scores)
            max_score = max(self.blosum[aa][aa] for aa in self.blosum)
            return float(score / max_score)

        except Exception as e:
            self.logger.error(f"Error getting residue conservation: {str(e)}")
            return 0.0

    def _calculate_window_conservation(
        self,
        sequence: str,
        window_size: int = 5,
    ) -> List[float]:
        """Calculate conservation scores using sliding window.

        Args:
            sequence: Amino acid sequence
            window_size: Size of sliding window

        Returns:
            List of conservation scores
        """
        try:
            scores = []
            half_window = window_size // 2

            for i in range(len(sequence)):
                # Get window sequence
                start = max(0, i - half_window)
                end = min(len(sequence), i + half_window + 1)
                window = sequence[start:end]

                # Calculate average BLOSUM scores
                window_scores = []
                for j, aa1 in enumerate(window):
                    for aa2 in window[j + 1 :]:
                        score = self.blosum[aa1][aa2]
                        window_scores.append(score)

                if window_scores:
                    avg_score = np.mean(window_scores)
                    max_score = max(self.blosum[aa][aa] for aa in self.blosum)
                    scores.append(float(avg_score / max_score))
                else:
                    scores.append(0.0)

            return scores

        except Exception as e:
            self.logger.error(f"Error calculating window conservation: {str(e)}")
            return [0.0] * len(sequence)

    def _calculate_relative_conservation(self, sequence: str) -> List[float]:
        """Calculate position-specific conservation relative to background.

        Args:
            sequence: Amino acid sequence

        Returns:
            List of relative conservation scores
        """
        try:
            scores = []
            background_freqs = self._get_background_frequencies()

            for aa in sequence:
                # Get BLOSUM scores relative to background
                rel_scores = []
                for bg_aa, freq in background_freqs.items():
                    score = self.blosum[aa][bg_aa] * freq
                    rel_scores.append(score)

                # Normalize score
                score = np.mean(rel_scores)
                max_score = max(self.blosum[aa][aa] for aa in self.blosum)
                scores.append(float(score / max_score))

            return scores

        except Exception as e:
            self.logger.error(f"Error calculating relative conservation: {str(e)}")
            return [0.0] * len(sequence)

    def _calculate_pssm_conservation(self, sequence: str) -> List[float]:
        """Calculate conservation using position-specific scoring matrix.

        Args:
            sequence: Amino acid sequence

        Returns:
            List of PSSM-based conservation scores
        """
        try:
            # Check cache first
            cache_key = hash(sequence)
            if cache_key in self.conservation_cache:
                return self.conservation_cache[cache_key]

            # Calculate position-specific frequencies
            aa_freqs = {}
            for aa in sequence:
                if aa not in aa_freqs:
                    aa_freqs[aa] = 0
                aa_freqs[aa] += 1

            # Normalize frequencies
            total = sum(aa_freqs.values())
            for aa in aa_freqs:
                aa_freqs[aa] /= total

            # Calculate PSSM scores
            scores = []
            background_freqs = self._get_background_frequencies()

            for aa in sequence:
                # Calculate log-odds score
                if aa in aa_freqs and aa in background_freqs:
                    score = np.log2(aa_freqs[aa] / background_freqs[aa])
                    scores.append(float(score))
                else:
                    scores.append(0.0)

            # Cache results
            self.conservation_cache[cache_key] = scores
            return scores

        except Exception as e:
            self.logger.error(f"Error calculating PSSM conservation: {str(e)}")
            return [0.0] * len(sequence)

    def _calculate_sequence_entropy(
        self,
        sequence: str,
        window_size: int = 5,
    ) -> List[float]:
        """Calculate sequence entropy in sliding window.

        Args:
            sequence: Amino acid sequence
            window_size: Size of sliding window

        Returns:
            List of entropy scores
        """
        try:
            scores = []
            half_window = window_size // 2

            for i in range(len(sequence)):
                # Get window sequence
                start = max(0, i - half_window)
                end = min(len(sequence), i + half_window + 1)
                window = sequence[start:end]

                # Calculate amino acid frequencies
                aa_freqs = {}
                for aa in window:
                    if aa not in aa_freqs:
                        aa_freqs[aa] = 0
                    aa_freqs[aa] += 1

                # Calculate entropy
                entropy = 0.0
                for count in aa_freqs.values():
                    p = count / len(window)
                    entropy -= p * np.log2(p)

                scores.append(float(entropy))

            return scores

        except Exception as e:
            self.logger.error(f"Error calculating sequence entropy: {str(e)}")
            return [0.0] * len(sequence)

    def _calculate_conservation_metrics(self, conservation: Dict[str, Any]) -> Dict[str, float]:
        """Calculate overall conservation metrics.

        Args:
            conservation: Dictionary of conservation data

        Returns:
            Dictionary of conservation metrics
        """
        try:
            metrics = {}

            # Window-based metrics
            if "window_scores" in conservation:
                scores = conservation["window_scores"]
                metrics.update(
                    {
                        "average_conservation": float(np.mean(scores)),
                        "conservation_std": float(np.std(scores)),
                        "min_conservation": float(np.min(scores)),
                        "max_conservation": float(np.max(scores)),
                    }
                )

            # Entropy-based metrics
            if "entropy_scores" in conservation:
                entropy = conservation["entropy_scores"]
                metrics.update(
                    {
                        "average_entropy": float(np.mean(entropy)),
                        "entropy_std": float(np.std(entropy)),
                        "min_entropy": float(np.min(entropy)),
                        "max_entropy": float(np.max(entropy)),
                    }
                )

            # PSSM-based metrics
            if "pssm_scores" in conservation:
                pssm = conservation["pssm_scores"]
                metrics.update(
                    {
                        "average_pssm": float(np.mean(pssm)),
                        "pssm_std": float(np.std(pssm)),
                        "min_pssm": float(np.min(pssm)),
                        "max_pssm": float(np.max(pssm)),
                    }
                )

            # AlphaFold confidence metrics
            if "alphafold_scores" in conservation:
                af_scores = conservation["alphafold_scores"]
                if af_scores:
                    metrics.update(
                        {
                            "average_confidence": float(np.mean(list(af_scores.values()))),
                            "min_confidence": float(np.min(list(af_scores.values()))),
                            "max_confidence": float(np.max(list(af_scores.values()))),
                        }
                    )

            return metrics

        except Exception as e:
            self.logger.error(f"Error calculating conservation metrics: {str(e)}")
            return {}

    def _get_alphafold_confidence(self, structure: Structure) -> Dict[int, float]:
        """Get AlphaFold confidence scores.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary mapping residue numbers to confidence scores
        """
        try:
            # Get UniProt ID from structure
            uniprot_id = structure.header.get("idcode", "").split("_")[0]

            if not uniprot_id:
                return {}

            # Query AlphaFold DB
            url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
            response = requests.get(url)

            if not response.ok:
                return {}

            data = response.json()

            # Extract confidence scores
            confidence = {}
            for residue in structure.get_residues():
                res_id = residue.get_id()[1]
                if str(res_id) in data["plddt"]:
                    confidence[res_id] = float(data["plddt"][str(res_id)]) / 100.0

            return confidence

        except Exception as e:
            self.logger.error(f"Error getting AlphaFold confidence: {str(e)}")
            return {}

    def _calculate_percentile(self, scores: List[float]) -> float:
        """Calculate percentile of average score.

        Args:
            scores: List of conservation scores

        Returns:
            Percentile (0-1)
        """
        try:
            if not scores:
                return 0.0

            avg_score = np.mean(scores)
            return float(len([s for s in scores if s < avg_score]) / len(scores))

        except Exception as e:
            self.logger.error(f"Error calculating percentile: {str(e)}")
            return 0.0

    def _get_background_frequencies(self) -> Dict[str, float]:
        """Get background amino acid frequencies.

        Returns:
            Dictionary mapping amino acids to frequencies
        """
        # Standard amino acid frequencies from SwissProt
        return {
            "A": 0.0825,
            "R": 0.0553,
            "N": 0.0406,
            "D": 0.0545,
            "C": 0.0137,
            "Q": 0.0393,
            "E": 0.0675,
            "G": 0.0707,
            "H": 0.0227,
            "I": 0.0595,
            "L": 0.0966,
            "K": 0.0584,
            "M": 0.0242,
            "F": 0.0386,
            "P": 0.0470,
            "S": 0.0657,
            "T": 0.0534,
            "W": 0.0108,
            "Y": 0.0292,
            "V": 0.0687,
        }
