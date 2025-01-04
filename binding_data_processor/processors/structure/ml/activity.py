"""Enhanced activity prediction and analysis using machine learning.

This module provides comprehensive functionality for:
1. Binding affinity prediction using graph neural networks and ensemble methods
2. Activity classification and regression with uncertainty estimation 
3. Target prediction and binding site analysis
4. Structure-activity relationship analysis
5. Psychopharmacological activity prediction
6. Toxicity and abuse potential assessment
7. Web data integration and enrichment
8. Model interpretability and visualization
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import torch
from rdkit import Chem

from .models.gnn import EnhancedGNN
from .models.ensemble import EnsembleModel
from .predictors.affinity import AffinityPredictor
from .predictors.activity import ActivityPredictor
from .predictors.toxicity import ToxicityPredictor
from .predictors.abuse import AbusePotentialPredictor
from .analysis.sar import SARAnalyzer
from .analysis.visualization import PredictionVisualizer
from ..pharmacophore import PharmacophoreGenerator
from ...web_enrichment import WebDataEnricher


class ActivityPredictionPipeline:
    """Comprehensive activity prediction and analysis pipeline."""

    def __init__(
        self,
        model_dir: Optional[str] = None,
        use_ensemble: bool = True,
        uncertainty: bool = True,
        device: Optional[str] = None,
    ):
        """Initialize prediction pipeline.

        Args:
            model_dir: Directory containing pre-trained models
            use_ensemble: Whether to use ensemble models
            uncertainty: Whether to estimate prediction uncertainty
            device: Device to run models on
        """
        self.logger = logging.getLogger(__name__)
        self.model_dir = Path(model_dir) if model_dir else None
        self.use_ensemble = use_ensemble
        self.uncertainty = uncertainty
        self.device = device or "cuda" if torch.cuda.is_available() else "cpu"

        # Initialize predictors
        self.affinity_predictor = AffinityPredictor(
            model_dir=self.model_dir,
            use_ensemble=use_ensemble,
            uncertainty=uncertainty,
            device=self.device,
        )

        self.activity_predictor = ActivityPredictor(
            model_dir=self.model_dir,
            use_ensemble=use_ensemble,
            uncertainty=uncertainty,
            device=self.device,
        )

        self.toxicity_predictor = ToxicityPredictor(
            model_dir=self.model_dir,
            use_ensemble=use_ensemble,
            uncertainty=uncertainty,
            device=self.device,
        )

        self.abuse_predictor = AbusePotentialPredictor(
            model_dir=self.model_dir,
            use_ensemble=use_ensemble,
            uncertainty=uncertainty,
            device=self.device,
        )

        # Initialize analyzers
        self.sar_analyzer = SARAnalyzer()
        self.visualizer = PredictionVisualizer()

        # Initialize enrichment components
        self.pharmacophore_gen = PharmacophoreGenerator()
        self.web_enricher = WebDataEnricher()

    def predict_compound(
        self,
        compound: Union[str, Chem.Mol],
        include_web_data: bool = True,
        confidence_threshold: float = 0.5,
    ) -> Dict:
        """Comprehensive compound prediction and analysis.

        Args:
            compound: SMILES string or RDKit molecule
            include_web_data: Whether to include web data
            confidence_threshold: Minimum confidence threshold

        Returns:
            Dictionary containing:
            - Binding affinity predictions
            - Activity predictions
            - Toxicity predictions
            - Abuse potential predictions
            - Structural analysis
            - Web data (if requested)
            - Visualization data
        """
        try:
            # Convert input to molecule
            mol = (
                compound
                if isinstance(compound, Chem.Mol)
                else Chem.MolFromSmiles(compound)
            )
            if mol is None:
                raise ValueError("Invalid compound input")

            # Get predictions
            results = {}

            # Binding affinity
            affinity_results = self.affinity_predictor.predict(mol)
            results.update(affinity_results)

            # Activity
            activity_results = self.activity_predictor.predict(
                mol, confidence_threshold
            )
            results.update(activity_results)

            # Toxicity
            toxicity_results = self.toxicity_predictor.predict(
                mol, confidence_threshold
            )
            results.update(toxicity_results)

            # Abuse potential
            abuse_results = self.abuse_predictor.predict(mol, confidence_threshold)
            results.update(abuse_results)

            # Add pharmacophore analysis
            results["pharmacophore"] = self.pharmacophore_gen.generate(mol)

            # Add structural analysis
            results["structural_analysis"] = self.sar_analyzer.analyze_compound(mol)

            # Add visualization data
            results["visualization"] = self.visualizer.generate_visualization(
                mol, results
            )

            # Add web data if requested
            if include_web_data:
                web_data = self.web_enricher.enrich_compound(Chem.MolToSmiles(mol))
                results.update(web_data)

            return results

        except Exception as e:
            self.logger.error(f"Error in predict_compound: {str(e)}")
            return {}

    def analyze_compound_series(
        self,
        compounds: List[Union[str, Chem.Mol]],
        activities: Optional[np.ndarray] = None,
        threshold: float = 0.5,
    ) -> Dict:
        """Analyze series of compounds.

        Args:
            compounds: List of compounds
            activities: Optional known activity values
            threshold: Activity threshold

        Returns:
            Dictionary containing:
            - SAR analysis
            - Activity patterns
            - Structural patterns
            - Series visualization
        """
        try:
            # Convert compounds to molecules
            mols = []
            for comp in compounds:
                mol = comp if isinstance(comp, Chem.Mol) else Chem.MolFromSmiles(comp)
                if mol is not None:
                    mols.append(mol)

            if not mols:
                return {}

            # Get predictions if activities not provided
            if activities is None:
                activities = np.array(
                    [
                        float(self.affinity_predictor.predict(mol)["binding_affinity"])
                        for mol in mols
                    ]
                )

            # Analyze SAR
            results = self.sar_analyzer.analyze_series(mols, activities, threshold)

            # Add series visualization
            results["visualization"] = self.visualizer.visualize_series(
                mols, activities
            )

            return results

        except Exception as e:
            self.logger.error(f"Error in analyze_compound_series: {str(e)}")
            return {}

    def batch_predict(
        self,
        compounds: List[Union[str, Chem.Mol]],
        batch_size: int = 32,
        include_web_data: bool = True,
        confidence_threshold: float = 0.5,
    ) -> List[Dict]:
        """Make predictions for multiple compounds.

        Args:
            compounds: List of compounds
            batch_size: Batch size for predictions
            include_web_data: Whether to include web data
            confidence_threshold: Minimum confidence threshold

        Returns:
            List of prediction dictionaries
        """
        try:
            results = []

            # Process in batches
            for i in range(0, len(compounds), batch_size):
                batch = compounds[i : i + batch_size]

                # Get predictions for batch
                batch_results = [
                    self.predict_compound(
                        comp,
                        include_web_data=include_web_data,
                        confidence_threshold=confidence_threshold,
                    )
                    for comp in batch
                ]

                results.extend(batch_results)

            return results

        except Exception as e:
            self.logger.error(f"Error in batch_predict: {str(e)}")
            return []
