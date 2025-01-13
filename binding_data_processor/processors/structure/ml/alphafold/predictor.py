"""AlphaFold structure prediction."""

import logging
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, Any, Union
import json
import requests
import numpy as np
from Bio.PDB import Structure, PDBParser, PDBIO

from .config import AlphaFoldConfig

logger = logging.getLogger(__name__)


@dataclass
class AlphaFoldResult:
    """Result from AlphaFold prediction."""

    structure: Structure
    confidence: Dict[str, float]  # Contains plddt, pae scores
    plddt: np.ndarray  # Per-residue confidence scores
    pae: Optional[np.ndarray] = None  # Predicted aligned error matrix
    raw_scores: Optional[Dict[str, Any]] = None  # Additional prediction metrics


class AlphaFoldPredictor:
    """Predicts protein structures using AlphaFold."""

    def __init__(
        self,
        config: Optional[Union[Dict, AlphaFoldConfig]] = None,
        cache_dir: Optional[Union[str, Path]] = None,
    ):
        """Initialize predictor.

        Args:
            config: Configuration for AlphaFold
            cache_dir: Directory for caching predictions
        """
        if isinstance(config, dict):
            config = AlphaFoldConfig(**config)
        elif config is None:
            config = AlphaFoldConfig()

        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)

        if cache_dir:
            self.cache_dir = Path(cache_dir)
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        else:
            self.cache_dir = None

        self.pdb_parser = PDBParser(QUIET=True)

    async def predict(self, sequence: str) -> Optional[AlphaFoldResult]:
        """Predict protein structure using AlphaFold.

        Args:
            sequence: Amino acid sequence

        Returns:
            AlphaFoldResult containing structure and confidence metrics,
            or None if prediction fails
        """
        try:
            # Check cache first
            cached = self._load_cached_prediction(sequence)
            if cached:
                return cached

            # Use API by default
            if self.config.use_alphafold_api:
                return await self._predict_with_api(sequence)

            # Fall back to local prediction
            if self.config.use_templates:
                return await self._predict_with_templates(sequence)
            else:
                return await self._predict_no_templates(sequence)

        except Exception as e:
            self.logger.error(f"Error in AlphaFold prediction: {str(e)}")
            return None

    async def _predict_with_api(self, sequence: str) -> Optional[AlphaFoldResult]:
        """Predict structure using AlphaFold API."""
        try:
            # Search for existing prediction
            search_params = {"query": sequence}
            response = requests.get(
                f"{self.config.api_base_url}/prediction/search",
                params=search_params,
                timeout=self.config.api_timeout,
            )
            response.raise_for_status()
            results = response.json()

            if not results:
                self.logger.warning("No AlphaFold structure found for sequence")
                return None

            # Get best model by confidence score
            best_model = max(results, key=lambda x: x["confidence_score"])
            confidence = {
                "plddt": best_model["confidence_score"],
                "pae": best_model.get("pae_score"),
            }

            # Download structure
            download_params = {"id": best_model["id"]}
            response = requests.get(
                f"{self.config.api_base_url}/prediction/download",
                params=download_params,
                timeout=self.config.api_timeout,
            )
            response.raise_for_status()

            # Parse structure
            structure = self.pdb_parser.get_structure("prediction", response.content)

            # Extract pLDDT scores from B-factors
            plddt = np.array([atom.get_bfactor() for atom in structure.get_atoms()])

            # Create result
            result = AlphaFoldResult(
                structure=structure,
                confidence=confidence,
                plddt=plddt,
                pae=best_model.get("pae_matrix"),
                raw_scores=best_model.get("raw_scores"),
            )

            # Cache result
            self._cache_prediction(sequence, result)

            return result

        except Exception as e:
            self.logger.error(f"Error in AlphaFold API prediction: {str(e)}")
            return None

    async def _predict_with_templates(self, sequence: str) -> Optional[AlphaFoldResult]:
        """Predict structure using templates."""
        try:
            if not self.config.model_dir:
                raise ValueError("Model directory not configured for template prediction")

            # TODO: Implement template-based prediction using local AlphaFold
            self.logger.warning("Template-based prediction not yet implemented")
            return None

        except Exception as e:
            self.logger.error(f"Error in template prediction: {str(e)}")
            return None

    async def _predict_no_templates(self, sequence: str) -> Optional[AlphaFoldResult]:
        """Predict structure without templates."""
        try:
            if not self.config.model_dir:
                raise ValueError("Model directory not configured for template-free prediction")

            # TODO: Implement template-free prediction using local AlphaFold
            self.logger.warning("Template-free prediction not yet implemented")
            return None

        except Exception as e:
            self.logger.error(f"Error in template-free prediction: {str(e)}")
            return None

    def _cache_prediction(self, sequence: str, result: AlphaFoldResult) -> None:
        """Cache prediction results."""
        if not self.cache_dir:
            return

        try:
            # Create cache key from sequence
            cache_key = hashlib.sha256(sequence.encode()).hexdigest()

            # Save structure as PDB file
            structure_path = self.cache_dir / f"{cache_key}.pdb"
            io = PDBIO()
            io.set_structure(result.structure)
            io.save(str(structure_path))

            # Save metadata
            meta_path = self.cache_dir / f"{cache_key}.json"
            meta = {
                "sequence": sequence,
                "confidence": result.confidence,
                "plddt": result.plddt.tolist() if result.plddt is not None else None,
                "pae": result.pae.tolist() if result.pae is not None else None,
                "raw_scores": result.raw_scores,
            }
            with open(meta_path, "w") as f:
                json.dump(meta, f)

        except Exception as e:
            self.logger.error(f"Error caching prediction: {str(e)}")

    def _load_cached_prediction(self, sequence: str) -> Optional[AlphaFoldResult]:
        """Load cached prediction if available."""
        if not self.cache_dir:
            return None

        try:
            # Get cache key
            cache_key = hashlib.sha256(sequence.encode()).hexdigest()

            # Check if cache exists
            structure_path = self.cache_dir / f"{cache_key}.pdb"
            meta_path = self.cache_dir / f"{cache_key}.json"

            if not (structure_path.exists() and meta_path.exists()):
                return None

            # Load structure
            structure = self.pdb_parser.get_structure("cached", str(structure_path))

            # Load metadata
            with open(meta_path) as f:
                meta = json.load(f)

            # Verify sequence matches
            if meta["sequence"] != sequence:
                return None

            # Convert arrays back to numpy
            plddt = np.array(meta["plddt"]) if meta.get("plddt") is not None else None
            pae = np.array(meta["pae"]) if meta.get("pae") is not None else None

            return AlphaFoldResult(
                structure=structure,
                confidence=meta["confidence"],
                plddt=plddt,
                pae=pae,
                raw_scores=meta.get("raw_scores"),
            )

        except Exception as e:
            self.logger.error(f"Error loading cached prediction: {str(e)}")
            return None
