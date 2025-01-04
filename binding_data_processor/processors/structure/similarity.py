"""Enhanced chemical similarity calculations with ML capabilities."""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    AllChem,
    MACCSkeys,
    rdFingerprintGenerator,
    rdMolDescriptors,
    rdFMCS,
    rdShapeHelpers,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Pharm3D import Pharmacophore
from sklearn.ensemble import RandomForestClassifier, IsolationForest
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

from .base import BaseStructureProcessor


class SimilarityProcessor(BaseStructureProcessor):
    """Enhanced molecular similarity calculations with ML capabilities."""

    # Comprehensive fingerprint types with detailed configuration
    FINGERPRINT_TYPES = {
        "morgan": {
            "description": "Morgan (ECFP-like) fingerprints",
            "radius": 2,
            "bits": 2048,
            "generator": lambda mol, radius: AllChem.GetMorganFingerprintAsBitVect(
                mol, radius, nBits=2048
            ),
        },
        "maccs": {
            "description": "MACCS keys (166 bits)",
            "bits": 166,
            "generator": lambda mol, _: MACCSkeys.GenMACCSKeys(mol),
        },
        "topological": {
            "description": "Topological (path-based) fingerprints",
            "bits": 2048,
            "min_path": 1,
            "max_path": 7,
            "generator": lambda mol, _: rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
                mol, nBits=2048
            ),
        },
        "atom_pairs": {
            "description": "Atom pairs fingerprints",
            "bits": 2048,
            "max_length": 30,
            "generator": lambda mol, _: rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
                mol, nBits=2048
            ),
        },
        "torsion": {
            "description": "Topological torsion fingerprints",
            "bits": 2048,
            "generator": lambda mol, _: rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
                mol, nBits=2048
            ),
        },
        "pharmacophore": {
            "description": "3D pharmacophore fingerprints",
            "bits": 2048,
            "generator": lambda mol, _: Pharmacophore.Generate3DFingerprint(mol),
        },
    }

    # Enhanced similarity metrics with descriptions
    SIMILARITY_METRICS = {
        "tanimoto": {
            "description": "Tanimoto/Jaccard similarity",
            "function": DataStructs.TanimotoSimilarity,
        },
        "dice": {
            "description": "Dice/Sørensen similarity",
            "function": DataStructs.DiceSimilarity,
        },
        "cosine": {
            "description": "Cosine similarity",
            "function": DataStructs.CosineSimilarity,
        },
        "sokal": {
            "description": "Sokal-Sneath similarity",
            "function": DataStructs.SokalSimilarity,
        },
        "russel": {
            "description": "Russell-Rao similarity",
            "function": DataStructs.RusselSimilarity,
        },
        "kulczynski": {
            "description": "Kulczynski similarity",
            "function": DataStructs.KulczynskiSimilarity,
        },
        "mcconnaughey": {
            "description": "McConnaughey similarity",
            "function": DataStructs.McConnaugheySimilarity,
        },
        "shape": {
            "description": "3D shape similarity",
            "function": rdShapeHelpers.ShapeTanimotoMol,
        },
    }

    def __init__(self, cache_size: int = 1000):
        """Initialize similarity processor."""
        super().__init__()
        self.logger = logging.getLogger(__name__)
        self._init_fingerprint_generators()
        self._fp_cache = {}
        self._cache_size = cache_size
        self._ml_models = {}
        self._scaler = StandardScaler()

    def _init_fingerprint_generators(self):
        """Initialize optimized fingerprint generators."""
        try:
            self.fp_generators = {
                "morgan": rdFingerprintGenerator.GetMorganGenerator(
                    radius=self.FINGERPRINT_TYPES["morgan"]["radius"],
                    fpSize=self.FINGERPRINT_TYPES["morgan"]["bits"],
                ),
                "topological": rdFingerprintGenerator.GetRDKitFPGenerator(
                    minPath=self.FINGERPRINT_TYPES["topological"]["min_path"],
                    maxPath=self.FINGERPRINT_TYPES["topological"]["max_path"],
                    fpSize=self.FINGERPRINT_TYPES["topological"]["bits"],
                ),
                "atom_pairs": rdFingerprintGenerator.GetAtomPairGenerator(
                    maxLength=self.FINGERPRINT_TYPES["atom_pairs"]["max_length"],
                    fpSize=self.FINGERPRINT_TYPES["atom_pairs"]["bits"],
                ),
                "torsion": rdFingerprintGenerator.GetTopologicalTorsionGenerator(
                    fpSize=self.FINGERPRINT_TYPES["torsion"]["bits"],
                ),
            }
        except Exception as e:
            self.logger.error(f"Error initializing fingerprint generators: {str(e)}")

    def _get_cached_fingerprint(
        self, mol: Chem.Mol, fp_type: str, radius: int
    ) -> Optional[DataStructs.ExplicitBitVect]:
        """Get fingerprint from cache or generate new one."""
        try:
            if mol is None:
                return None

            cache_key = (Chem.MolToSmiles(mol), fp_type, radius)
            if cache_key in self._fp_cache:
                return self._fp_cache[cache_key]

            # Generate fingerprint using appropriate generator
            if fp_type in self.FINGERPRINT_TYPES:
                fp = self.FINGERPRINT_TYPES[fp_type]["generator"](mol, radius)
                if fp is not None and len(self._fp_cache) < self._cache_size:
                    self._fp_cache[cache_key] = fp
                return fp

            return None

        except Exception as e:
            self.logger.error(f"Error getting cached fingerprint: {str(e)}")
            return None

    def calculate_similarity(
        self,
        mol1: Chem.Mol,
        mol2: Chem.Mol,
        method: str = "tanimoto",
        fp_type: str = "morgan",
        radius: int = 2,
        use_3d: bool = False,
        **kwargs,
    ) -> Optional[float]:
        """
        Calculate chemical similarity between two molecules.

        Args:
            mol1: First RDKit molecule
            mol2: Second RDKit molecule
            method: Similarity metric to use
            fp_type: Fingerprint type to use
            radius: Radius for Morgan fingerprints
            use_3d: Use 3D conformers for shape/pharm similarity
            **kwargs: Additional parameters for similarity calculation

        Returns:
            Similarity score or None if error
        """
        try:
            if mol1 is None or mol2 is None:
                return None

            # Handle 3D similarity methods
            if use_3d:
                if method == "shape":
                    return rdShapeHelpers.ShapeTanimotoMol(mol1, mol2)
                elif method == "pharmacophore":
                    return self._calculate_pharmacophore_similarity(mol1, mol2)

            # Get similarity function
            if method not in self.SIMILARITY_METRICS:
                self.logger.warning(f"Unknown similarity method: {method}")
                return None
            sim_func = self.SIMILARITY_METRICS[method]["function"]

            # Generate fingerprints using cache
            fp1 = self._get_cached_fingerprint(mol1, fp_type, radius)
            fp2 = self._get_cached_fingerprint(mol2, fp_type, radius)

            if fp1 is None or fp2 is None:
                return None

            # Calculate similarity
            return float(sim_func(fp1, fp2))

        except Exception as e:
            self.logger.error(f"Error calculating similarity: {str(e)}")
            return None

    def _calculate_pharmacophore_similarity(
        self, mol1: Chem.Mol, mol2: Chem.Mol
    ) -> Optional[float]:
        """Calculate pharmacophore-based similarity."""
        try:
            # Generate 3D conformers if needed
            for mol in [mol1, mol2]:
                if not mol.GetNumConformers():
                    AllChem.EmbedMolecule(mol, randomSeed=42)
                    AllChem.MMFFOptimizeMolecule(mol)

            # Generate pharmacophore fingerprints
            fp1 = Pharmacophore.Generate3DFingerprint(mol1)
            fp2 = Pharmacophore.Generate3DFingerprint(mol2)

            return DataStructs.TanimotoSimilarity(fp1, fp2)

        except Exception as e:
            self.logger.error(f"Error calculating pharmacophore similarity: {str(e)}")
            return None

    def generate_fingerprints_batch(
        self,
        mols: List[Chem.Mol],
        fp_type: str = "morgan",
        radius: int = 2,
        n_jobs: int = 1,
    ) -> List[Optional[DataStructs.ExplicitBitVect]]:
        """
        Generate fingerprints for multiple molecules in parallel.

        Args:
            mols: List of RDKit molecules
            fp_type: Type of fingerprint to generate
            radius: Radius for Morgan fingerprints
            n_jobs: Number of parallel jobs

        Returns:
            List of fingerprints (None for failed molecules)
        """
        try:
            if not mols:
                return []

            # Use cached fingerprints where available
            fingerprints = []
            for mol in mols:
                if mol is not None:
                    fp = self._get_cached_fingerprint(mol, fp_type, radius)
                    fingerprints.append(fp)
                else:
                    fingerprints.append(None)

            return fingerprints

        except Exception as e:
            self.logger.error(f"Error generating batch fingerprints: {str(e)}")
            return [None] * len(mols)

    def calculate_similarity_matrix(
        self,
        mols: List[Chem.Mol],
        method: str = "tanimoto",
        fp_type: str = "morgan",
        radius: int = 2,
        n_jobs: int = 1,
    ) -> Optional[np.ndarray]:
        """
        Calculate similarity matrix for a list of molecules.

        Args:
            mols: List of RDKit molecules
            method: Similarity metric to use
            fp_type: Fingerprint type to use
            radius: Radius for Morgan fingerprints
            n_jobs: Number of parallel jobs

        Returns:
            Similarity matrix as numpy array or None if error
        """
        try:
            if not mols:
                return None

            n_mols = len(mols)
            sim_matrix = np.zeros((n_mols, n_mols))

            # Generate all fingerprints in batch
            fingerprints = []
            for mol in mols:
                fp = self._get_cached_fingerprint(mol, fp_type, radius)
                if fp is None:
                    return None
                fingerprints.append(fp)

            # Get similarity function
            if method not in self.SIMILARITY_METRICS:
                return None
            sim_func = self.SIMILARITY_METRICS[method]["function"]

            # Calculate similarity matrix efficiently
            for i in range(n_mols):
                sim_matrix[i, i] = 1.0  # Self-similarity
                for j in range(i + 1, n_mols):
                    sim = sim_func(fingerprints[i], fingerprints[j])
                    sim_matrix[i, j] = sim_matrix[j, i] = sim

            return sim_matrix

        except Exception as e:
            self.logger.error(f"Error calculating similarity matrix: {str(e)}")
            return None

    def find_similar_compounds(
        self,
        query_mol: Chem.Mol,
        library_mols: List[Chem.Mol],
        threshold: float = 0.7,
        method: str = "tanimoto",
        fp_type: str = "morgan",
        radius: int = 2,
        use_ml: bool = False,
    ) -> List[Tuple[int, float]]:
        """
        Find similar compounds in a library.

        Args:
            query_mol: Query molecule
            library_mols: List of library molecules
            threshold: Similarity threshold
            method: Similarity metric to use
            fp_type: Fingerprint type to use
            radius: Radius for Morgan fingerprints
            use_ml: Use ML model if available

        Returns:
            List of (compound_index, similarity) tuples
        """
        try:
            if query_mol is None or not library_mols:
                return []

            # Use ML model if available and requested
            if use_ml and fp_type in self._ml_models:
                return self.predict_similarity_ml(
                    query_mol, library_mols, fp_type, threshold
                )

            # Traditional fingerprint similarity
            query_fp = self._get_cached_fingerprint(query_mol, fp_type, radius)
            if query_fp is None:
                return []

            # Get similarity function
            if method not in self.SIMILARITY_METRICS:
                return []
            sim_func = self.SIMILARITY_METRICS[method]["function"]

            # Calculate similarities
            similar_compounds = []
            for i, mol in enumerate(library_mols):
                if mol is not None:
                    fp = self._get_cached_fingerprint(mol, fp_type, radius)
                    if fp is not None:
                        sim = sim_func(query_fp, fp)
                        if sim >= threshold:
                            similar_compounds.append((i, float(sim)))

            similar_compounds.sort(key=lambda x: x[1], reverse=True)
            return similar_compounds

        except Exception as e:
            self.logger.error(f"Error finding similar compounds: {str(e)}")
            return []

    def train_similarity_model(
        self,
        active_mols: List[Chem.Mol],
        inactive_mols: List[Chem.Mol],
        fp_type: str = "morgan",
        model_type: str = "rf",
    ) -> bool:
        """
        Train ML model for similarity prediction.

        Args:
            active_mols: List of active molecules
            inactive_mols: List of inactive molecules
            fp_type: Fingerprint type to use
            model_type: ML model type ('rf' or 'knn')

        Returns:
            Success status
        """
        try:
            # Generate fingerprints for all molecules
            X = []
            y = []

            for mol in active_mols:
                fp = self._get_cached_fingerprint(mol, fp_type, radius=2)
                if fp is not None:
                    X.append(list(fp.GetOnBits()))
                    y.append(1)

            for mol in inactive_mols:
                fp = self._get_cached_fingerprint(mol, fp_type, radius=2)
                if fp is not None:
                    X.append(list(fp.GetOnBits()))
                    y.append(0)

            if not X:
                return False

            # Scale features
            X = self._scaler.fit_transform(X)

            # Train model
            if model_type == "rf":
                model = RandomForestClassifier(
                    n_estimators=100,
                    max_depth=10,
                    random_state=42,
                    class_weight="balanced",
                )
            elif model_type == "knn":
                model = NearestNeighbors(
                    n_neighbors=5,
                    algorithm="auto",
                    metric="euclidean",
                )
            else:
                return False

            model.fit(X, y)
            self._ml_models[fp_type] = model
            return True

        except Exception as e:
            self.logger.error(f"Error training similarity model: {str(e)}")
            return False

    def predict_similarity_ml(
        self,
        query_mol: Chem.Mol,
        library_mols: List[Chem.Mol],
        fp_type: str = "morgan",
        threshold: float = 0.5,
    ) -> List[Tuple[int, float]]:
        """
        Predict similar compounds using ML model.

        Args:
            query_mol: Query molecule
            library_mols: List of library molecules
            fp_type: Fingerprint type
            threshold: Similarity threshold

        Returns:
            List of (index, similarity) tuples
        """
        try:
            if query_mol is None or not library_mols or fp_type not in self._ml_models:
                return []

            # Generate fingerprints
            query_fp = self._get_cached_fingerprint(query_mol, fp_type, radius=2)
            if query_fp is None:
                return []

            library_fps = []
            valid_indices = []
            for i, mol in enumerate(library_mols):
                fp = self._get_cached_fingerprint(mol, fp_type, radius=2)
                if fp is not None:
                    library_fps.append(list(fp.GetOnBits()))
                    valid_indices.append(i)

            if not library_fps:
                return []

            # Scale features
            X = self._scaler.transform(library_fps)
            query_fp = self._scaler.transform([list(query_fp.GetOnBits())])

            # Get predictions
            model = self._ml_models[fp_type]
            if isinstance(model, RandomForestClassifier):
                probs = model.predict_proba(X)[:, 1]
            else:  # NearestNeighbors
                distances, _ = model.kneighbors(query_fp)
                probs = 1 / (1 + distances.flatten())

            # Filter and sort results
            similar_compounds = []
            for idx, prob in zip(valid_indices, probs):
                if prob >= threshold:
                    similar_compounds.append((idx, float(prob)))

            similar_compounds.sort(key=lambda x: x[1], reverse=True)
            return similar_compounds

        except Exception as e:
            self.logger.error(f"Error predicting similarities: {str(e)}")
            return []

    def detect_activity_cliffs(
        self,
        mols: List[Chem.Mol],
        activities: List[float],
        similarity_threshold: float = 0.7,
        activity_ratio_threshold: float = 10.0,
        fp_type: str = "morgan",
    ) -> List[Dict]:
        """
        Detect activity cliffs using ML.

        Args:
            mols: List of molecules
            activities: List of activity values
            similarity_threshold: Structural similarity threshold
            activity_ratio_threshold: Activity difference threshold
            fp_type: Fingerprint type

        Returns:
            List of activity cliff information
        """
        try:
            if not mols or len(mols) != len(activities):
                return []

            # Calculate similarity matrix
            sim_matrix = self.calculate_similarity_matrix(mols, fp_type=fp_type)
            if sim_matrix is None:
                return []

            # Use Isolation Forest to detect anomalies
            X = []
            for i in range(len(mols)):
                for j in range(i + 1, len(mols)):
                    if sim_matrix[i, j] >= similarity_threshold:
                        activity_ratio = max(
                            activities[i] / activities[j],
                            activities[j] / activities[i],
                        )
                        X.append([sim_matrix[i, j], activity_ratio])

            if not X:
                return []

            # Fit Isolation Forest
            iso_forest = IsolationForest(contamination=0.1, random_state=42)
            X = np.array(X)
            predictions = iso_forest.fit_predict(X)

            # Collect activity cliffs
            cliffs = []
            idx = 0
            for i in range(len(mols)):
                for j in range(i + 1, len(mols)):
                    if sim_matrix[i, j] >= similarity_threshold:
                        if predictions[idx] == -1:  # Anomaly detected
                            cliffs.append(
                                {
                                    "mol1_idx": i,
                                    "mol2_idx": j,
                                    "similarity": float(sim_matrix[i, j]),
                                    "activity_ratio": float(
                                        max(
                                            activities[i] / activities[j],
                                            activities[j] / activities[i],
                                        )
                                    ),
                                    "activity1": float(activities[i]),
                                    "activity2": float(activities[j]),
                                }
                            )
                        idx += 1

            return cliffs

        except Exception as e:
            self.logger.error(f"Error detecting activity cliffs: {str(e)}")
            return []
