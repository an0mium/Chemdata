"""Base classes for activity prediction."""

import abc
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union, Any, Tuple

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import AllChem

from ..base import MLProcessor
from ...base import BaseStructureProcessor
from ..models.graph_utils import mol_to_graph, compute_fingerprints

logger = logging.getLogger(__name__)


@dataclass
class ActivityPredictorConfig:
    """Configuration for activity predictors."""

    # Model configuration
    model_type: str = "graph"  # graph, ensemble, fingerprint
    model_path: Optional[str] = None
    device: str = "cuda" if torch.cuda.is_available() else "cpu"
    cache_dir: Optional[str] = None

    # Feature configuration
    fingerprint_types: List[str] = field(default_factory=lambda: ["morgan", "maccs", "topological"])
    use_3d: bool = False
    add_hs: bool = True
    compute_distances: bool = False

    # Training configuration
    batch_size: int = 32
    num_epochs: int = 100
    learning_rate: float = 0.001
    hidden_size: int = 128
    dropout: float = 0.1
    num_workers: int = 4

    # Optimization
    early_stopping: bool = True
    patience: int = 10
    random_seed: int = 42


class ActivityPredictor(MLProcessor, BaseStructureProcessor):
    """Base class for predicting molecular activity."""

    def __init__(
        self,
        config: Optional[Union[Dict, ActivityPredictorConfig]] = None,
        **kwargs,
    ):
        """Initialize activity predictor.

        Args:
            config: Predictor configuration
            **kwargs: Additional arguments passed to parent classes
        """
        if isinstance(config, dict):
            config = ActivityPredictorConfig(**config)
        elif config is None:
            config = ActivityPredictorConfig()

        super().__init__(
            model_path=config.model_path,
            device=config.device,
            **kwargs,
        )
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)

        if self.config.cache_dir:
            from pathlib import Path

            Path(self.config.cache_dir).mkdir(parents=True, exist_ok=True)

    def preprocess(
        self,
        mol: Chem.Mol,
        generate_3d: bool = False,
    ) -> Dict[str, Any]:
        """Preprocess molecule for prediction.

        Args:
            mol: Input molecule
            generate_3d: Whether to generate 3D conformer if needed

        Returns:
            Dictionary of features
        """
        try:
            if mol is None:
                return None

            # Generate 3D conformer if needed
            if generate_3d and not mol.GetNumConformers():
                AllChem.EmbedMolecule(mol, randomSeed=self.config.random_seed)
                AllChem.MMFFOptimizeMolecule(mol)

            features = {}

            # Graph features
            if self.config.model_type == "graph":
                graph_data = mol_to_graph(
                    mol,
                    add_hs=self.config.add_hs,
                    compute_distances=self.config.compute_distances and mol.GetNumConformers() > 0,
                    return_pyg=True,
                )
                features["graph"] = graph_data

            # Fingerprint features
            if self.config.fingerprint_types:
                fps = compute_fingerprints(
                    mol,
                    fp_types=self.config.fingerprint_types,
                    as_tensor=True,
                )
                features["fingerprints"] = fps

            # 3D features if available and requested
            if self.config.use_3d and mol.GetNumConformers() > 0:
                from ..models.graph_utils import compute_3d_descriptors

                descriptors_3d = compute_3d_descriptors(mol)
                features["3d"] = {k: torch.tensor([v], dtype=torch.float32) for k, v in descriptors_3d.items()}

            return features

        except Exception as e:
            self.logger.error(f"Error preprocessing molecule: {str(e)}")
            return None

    @abc.abstractmethod
    def predict(
        self,
        mol: Union[str, Chem.Mol],
        return_confidence: bool = False,
        **kwargs,
    ) -> Union[Dict[str, float], Tuple[Dict[str, float], Dict[str, float]]]:
        """Make predictions for molecule.

        Args:
            mol: Input molecule or SMILES
            return_confidence: Whether to return prediction confidence
            **kwargs: Additional arguments

        Returns:
            Dictionary of predictions, optionally with confidence scores
        """
        raise NotImplementedError

    def batch_predict(
        self,
        mols: List[Union[str, Chem.Mol]],
        return_confidence: bool = False,
        batch_size: Optional[int] = None,
        **kwargs,
    ) -> Union[List[Dict[str, float]], Tuple[List[Dict[str, float]], List[Dict[str, float]]]]:
        """Make predictions for multiple molecules.

        Args:
            mols: List of molecules or SMILES
            return_confidence: Whether to return prediction confidence
            batch_size: Batch size for predictions
            **kwargs: Additional arguments

        Returns:
            List of prediction dictionaries, optionally with confidence scores
        """
        try:
            if batch_size is None:
                batch_size = self.config.batch_size

            results = []
            confidences = []

            for i in range(0, len(mols), batch_size):
                batch = mols[i : i + batch_size]
                if return_confidence:
                    batch_preds, batch_confs = zip(*[self.predict(mol, return_confidence=True, **kwargs) for mol in batch])
                    results.extend(batch_preds)
                    confidences.extend(batch_confs)
                else:
                    batch_preds = [self.predict(mol, **kwargs) for mol in batch]
                    results.extend(batch_preds)

            if return_confidence:
                return results, confidences
            return results

        except Exception as e:
            self.logger.error(f"Error in batch prediction: {str(e)}")
            if return_confidence:
                return [{} for _ in mols], [{} for _ in mols]
            return [{} for _ in mols]

    def validate(
        self,
        mols: List[Union[str, Chem.Mol]],
        labels: List[Dict[str, float]],
        metrics: Optional[List[str]] = None,
        **kwargs,
    ) -> Dict[str, Dict[str, float]]:
        """Validate predictions against known labels.

        Args:
            mols: List of molecules or SMILES
            labels: List of label dictionaries
            metrics: Metrics to compute
            **kwargs: Additional arguments passed to batch_predict

        Returns:
            Dictionary of validation metrics
        """
        try:
            if metrics is None:
                metrics = ["rmse", "mae", "r2"]

            # Get predictions
            predictions = self.batch_predict(mols, **kwargs)

            # Compute metrics
            results = {}
            for target in labels[0].keys():
                target_results = {}
                y_true = np.array([label[target] for label in labels])
                y_pred = np.array([pred[target] for pred in predictions])

                for metric in metrics:
                    if metric == "rmse":
                        score = np.sqrt(np.mean((y_true - y_pred) ** 2))
                    elif metric == "mae":
                        score = np.mean(np.abs(y_true - y_pred))
                    elif metric == "r2":
                        from sklearn.metrics import r2_score

                        score = r2_score(y_true, y_pred)
                    else:
                        self.logger.warning(f"Unknown metric: {metric}")
                        continue

                    target_results[metric] = float(score)

                results[target] = target_results

            return results

        except Exception as e:
            self.logger.error(f"Error in validation: {str(e)}")
            return {}

    async def analyze(
        self,
        mol: Union[str, Chem.Mol],
        include_confidence: bool = True,
        **kwargs,
    ) -> Dict[str, Any]:
        """Analyze activity predictions.

        Args:
            mol: Input molecule or SMILES
            include_confidence: Whether to include prediction confidence
            **kwargs: Additional arguments

        Returns:
            Dictionary containing analysis results
        """
        try:
            if include_confidence:
                predictions, confidence = self.predict(mol, return_confidence=True, **kwargs)
                return {
                    "predictions": predictions,
                    "confidence": confidence,
                }
            else:
                predictions = self.predict(mol, **kwargs)
                return {"predictions": predictions}

        except Exception as e:
            self.logger.error(f"Error in activity analysis: {str(e)}")
            return {}
