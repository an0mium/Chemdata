"""Base psychopharmacological data processor.

This module provides the base PsychopharmProcessor class that orchestrates:
1. Data collection from various sources
2. ML model predictions
3. Data validation and standardization
4. Result aggregation and confidence scoring
5. Batch processing with progress tracking
6. Result analysis and reporting
"""

import logging
import time
from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed

from ...models.psychopharm import (
    PsychoactiveClass,
    NootropicMechanism,
    BBBPermeability,
)
from ...models.core import CompoundData
from ..structure.base import BaseStructureProcessor
from ...utils.cache import CacheManager
from ...utils.progress import ProgressTracker


@dataclass
class PredictionResult:
    """Container for prediction results with confidence scores."""

    value: any
    confidence: float
    supporting_data: Dict = None
    source: str = "predicted"


@dataclass
class BatchResult:
    """Results from batch processing."""

    successful: List[CompoundData]
    failed: List[Tuple[CompoundData, str]]  # (compound, error_message)
    stats: Dict[str, int]
    processing_time: float


class PsychopharmProcessor:
    """Process and predict psychopharmacological properties."""

    def __init__(
        self,
        structure_processor: Optional[BaseStructureProcessor] = None,
        cache_dir: Optional[str] = None,
        n_jobs: int = 4,
    ):
        """Initialize processor with optional dependencies."""
        self.logger = logging.getLogger(__name__)
        self.structure_processor = structure_processor
        self.cache = CacheManager(cache_dir) if cache_dir else None
        self.n_jobs = n_jobs
        self.progress = ProgressTracker()

        # Initialize predictors
        self.predictors = self._initialize_predictors()

        # Track mechanism frequencies
        self.nootropic_mechanism_counts: Dict[NootropicMechanism, int] = {
            mech: 0 for mech in NootropicMechanism
        }
        
        # Track related compounds
        self.similar_compound_groups: Set[Set[str]] = set()

    def _initialize_predictors(self) -> Dict:
        """Initialize ML predictors."""
        from .predictors import (
            BBBPredictor,
            ReceptorProfilePredictor,
            PsychoactiveClassPredictor,
            NootropicPredictor,
            AbusePotentialPredictor,
        )

        return {
            "bbb": BBBPredictor(),
            "receptor_profile": ReceptorProfilePredictor(),
            "psychoactive_class": PsychoactiveClassPredictor(),
            "nootropic": NootropicPredictor(),
            "abuse_potential": AbusePotentialPredictor(),
        }

    def process_batch(
        self,
        compounds: List[CompoundData],
        **kwargs
    ) -> BatchResult:
        """Process a batch of compounds in parallel."""
        start_time = time.time()
        successful = []
        failed = []
        
        # Process compounds in parallel
        with ThreadPoolExecutor(max_workers=self.n_jobs) as executor:
            future_to_compound = {
                executor.submit(self.process_compound, compound, **kwargs): compound
                for compound in compounds
            }
            
            for future in as_completed(future_to_compound):
                compound = future_to_compound[future]
                try:
                    result = future.result()
                    successful.append(result)
                    self._update_statistics(result)
                except Exception as e:
                    failed.append((compound, str(e)))

        # Generate statistics
        stats = self._generate_batch_statistics(successful)
        
        return BatchResult(
            successful=successful,
            failed=failed,
            stats=stats,
            processing_time=time.time() - start_time
        )

    def _update_statistics(self, compound: CompoundData) -> None:
        """Update mechanism frequencies and compound relationships."""
        # Update nootropic mechanism counts
        for mechanism in compound.nootropic_mechanisms:
            self.nootropic_mechanism_counts[mechanism] += 1
            
        # Update similar compound groups
        if self.structure_processor:
            similar_compounds = self.structure_processor.find_similar_compounds(
                compound.smiles
            )
            if similar_compounds:
                self.similar_compound_groups.add(
                    frozenset([compound.name] + similar_compounds)
                )

    def _generate_batch_statistics(self, compounds: List[CompoundData]) -> Dict:
        """Generate statistics from processed compounds."""
        stats = {
            "total_processed": len(compounds),
            "psychoactive_counts": {},
            "nootropic_mechanisms": dict(self.nootropic_mechanism_counts),
            "bbb_permeability": {},
            "receptor_coverage": 0,
            "similar_groups": len(self.similar_compound_groups),
        }
        
        for compound in compounds:
            # Count psychoactive classes
            if compound.psychoactive_class != PsychoactiveClass.UNKNOWN:
                stats["psychoactive_counts"][compound.psychoactive_class] = (
                    stats["psychoactive_counts"].get(compound.psychoactive_class, 0) + 1
                )
            
            # Count BBB permeability classes
            if compound.bbb_permeability != BBBPermeability.UNKNOWN:
                stats["bbb_permeability"][compound.bbb_permeability] = (
                    stats["bbb_permeability"].get(compound.bbb_permeability, 0) + 1
                )
            
            # Calculate receptor coverage
            if compound.receptor_profiles:
                stats["receptor_coverage"] += 1
        
        if compounds:
            stats["receptor_coverage"] /= len(compounds)
            
        return stats

    def process_compound(
        self,
        compound: CompoundData,
        predict_bbb: bool = True,
        predict_receptors: bool = True,
        predict_class: bool = True,
        predict_nootropic: bool = True,
        predict_abuse: bool = True,
        use_cache: bool = True,
    ) -> CompoundData:
        """Process a compound to predict psychopharmacological properties."""
        try:
            if self._try_load_from_cache(compound, use_cache):
                return compound

            self._run_predictions(
                compound,
                predict_bbb=predict_bbb,
                predict_receptors=predict_receptors,
                predict_class=predict_class,
                predict_nootropic=predict_nootropic,
                predict_abuse=predict_abuse,
            )

            if use_cache and self.cache:
                self._save_to_cache(compound)

            return compound

        except Exception as e:
            self.logger.error(f"Error processing compound {compound.name}: {str(e)}")
            return compound

    def _try_load_from_cache(self, compound: CompoundData, use_cache: bool) -> bool:
        """Try to load compound data from cache."""
        if use_cache and self.cache:
            cached = self.cache.get(f"psychopharm_{compound.inchi_key}")
            if cached:
                self._merge_cached_data(compound, cached)
                return True
        return False

    def _run_predictions(
        self,
        compound: CompoundData,
        predict_bbb: bool = True,
        predict_receptors: bool = True,
        predict_class: bool = True,
        predict_nootropic: bool = True,
        predict_abuse: bool = True,
    ) -> None:
        """Run all enabled predictions on compound."""
        with ThreadPoolExecutor(max_workers=self.n_jobs) as executor:
            futures = []

            prediction_tasks = [
                (predict_bbb, self._predict_bbb_properties),
                (predict_receptors, self._predict_receptor_profiles),
                (predict_class, self._predict_psychoactive_class),
                (predict_nootropic, self._predict_nootropic_properties),
                (predict_abuse, self._predict_abuse_potential),
            ]

            for enabled, predictor in prediction_tasks:
                if enabled:
                    futures.append(executor.submit(predictor, compound))

            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    self.logger.error(f"Prediction error: {str(e)}")

    def _save_to_cache(self, compound: CompoundData) -> None:
        """Save compound data to cache."""
        if self.cache:
            self.cache.set(
                f"psychopharm_{compound.inchi_key}",
                self._get_cacheable_data(compound),
            )

    def _predict_bbb_properties(self, compound: CompoundData) -> None:
        """Predict blood-brain barrier properties."""
        try:
            result = self.predictors["bbb"].predict(compound)
            compound.bbb_permeability = result.value
            compound.bbb_score = result.confidence
            compound.p_glycoprotein_substrate = result.supporting_data.get(
                "p_gp_substrate", False
            )
        except Exception as e:
            self.logger.error(f"BBB prediction error: {str(e)}")

    def _predict_receptor_profiles(self, compound: CompoundData) -> None:
        """Predict receptor binding profiles."""
        try:
            result = self.predictors["receptor_profile"].predict(compound)
            compound.receptor_profiles.update(result.value)
        except Exception as e:
            self.logger.error(f"Receptor profile prediction error: {str(e)}")

    def _predict_psychoactive_class(self, compound: CompoundData) -> None:
        """Predict psychoactive classification."""
        try:
            result = self.predictors["psychoactive_class"].predict(compound)
            if result.confidence > 0.5:  # Confidence threshold
                compound.psychoactive_class = result.value
                if result.supporting_data:
                    compound.effect_profile.update(result.supporting_data)
        except Exception as e:
            self.logger.error(f"Psychoactive class prediction error: {str(e)}")

    def _predict_nootropic_properties(self, compound: CompoundData) -> None:
        """Predict nootropic properties."""
        try:
            result = self.predictors["nootropic"].predict(compound)
            if result.confidence > 0.5:  # Confidence threshold
                compound.nootropic_mechanisms.update(result.value)
                if result.supporting_data:
                    compound.cognitive_effects.update(
                        result.supporting_data.get("effects", {})
                    )
                    compound.side_effects.update(
                        result.supporting_data.get("side_effects", {})
                    )
        except Exception as e:
            self.logger.error(f"Nootropic prediction error: {str(e)}")

    def _predict_abuse_potential(self, compound: CompoundData) -> None:
        """Predict abuse potential."""
        try:
            result = self.predictors["abuse_potential"].predict(compound)
            if result.supporting_data:
                compound.tolerance_profile.update(
                    result.supporting_data.get("tolerance", {})
                )
                compound.withdrawal_profile.update(
                    result.supporting_data.get("withdrawal", {})
                )
                compound.cross_tolerance.update(
                    result.supporting_data.get("cross_tolerance", set())
                )
        except Exception as e:
            self.logger.error(f"Abuse potential prediction error: {str(e)}")

    def _merge_cached_data(
        self, compound: CompoundData, cached_data: Dict
    ) -> None:
        """Merge cached psychopharm data into compound."""
        try:
            compound.bbb_permeability = cached_data.get(
                "bbb_permeability", BBBPermeability.UNKNOWN
            )
            compound.bbb_score = cached_data.get("bbb_score", 0.0)
            compound.p_glycoprotein_substrate = cached_data.get(
                "p_glycoprotein_substrate", False
            )
            compound.receptor_profiles.update(
                cached_data.get("receptor_profiles", {})
            )
            compound.psychoactive_class = cached_data.get(
                "psychoactive_class", PsychoactiveClass.UNKNOWN
            )
            compound.effect_profile.update(
                cached_data.get("effect_profile", {})
            )
            compound.nootropic_mechanisms.update(
                cached_data.get("nootropic_mechanisms", set())
            )
            compound.cognitive_effects.update(
                cached_data.get("cognitive_effects", {})
            )
            compound.side_effects.update(
                cached_data.get("side_effects", {})
            )
            compound.tolerance_profile.update(
                cached_data.get("tolerance_profile", {})
            )
            compound.withdrawal_profile.update(
                cached_data.get("withdrawal_profile", {})
            )
            compound.cross_tolerance.update(
                cached_data.get("cross_tolerance", set())
            )
        except Exception as e:
            self.logger.error(f"Error merging cached data: {str(e)}")

    def _get_cacheable_data(self, compound: CompoundData) -> Dict:
        """Extract cacheable psychopharm data from compound."""
        return {
            "bbb_permeability": compound.bbb_permeability,
            "bbb_score": compound.bbb_score,
            "p_glycoprotein_substrate": compound.p_glycoprotein_substrate,
            "receptor_profiles": compound.receptor_profiles,
            "psychoactive_class": compound.psychoactive_class,
            "effect_profile": compound.effect_profile,
            "nootropic_mechanisms": list(compound.nootropic_mechanisms),
            "cognitive_effects": compound.cognitive_effects,
            "side_effects": compound.side_effects,
            "tolerance_profile": compound.tolerance_profile,
            "withdrawal_profile": compound.withdrawal_profile,
            "cross_tolerance": list(compound.cross_tolerance),
        }
