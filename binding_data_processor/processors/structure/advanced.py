"""Advanced chemical structure analysis capabilities."""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union
import numpy as np
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    MACCSkeys,
    rdFingerprintGenerator,
    rdMolDescriptors,
    rdFMCS,
    rdShapeHelpers,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Pharm3D import Pharmacophore, Pharmacophore3D
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

from .base import BaseStructureProcessor


class AdvancedStructureProcessor(BaseStructureProcessor):
    """Advanced molecular structure analysis capabilities."""

    # Pharmacophore feature definitions
    PHARM_FEATURES = {
        "donor": {
            "pattern": "[!#6;!H0]-[!#6]",
            "family": "Donor",
        },
        "acceptor": {
            "pattern": "[$([!#6;+0]);!$([F,Cl,Br,I]);!$([o,s,nX3]);!$([Nv5,Pv5,Sv4,Sv6])]",
            "family": "Acceptor",
        },
        "aromatic": {
            "pattern": "a1aaaaa1",
            "family": "Aromatic",
        },
        "hydrophobe": {
            "pattern": "[#6,#7,#8][CX4H2][#6,#7,#8]",
            "family": "Hydrophobe",
        },
        "poscharge": {
            "pattern": "[+,+2,+3,+4]",
            "family": "PosIonizable",
        },
        "negcharge": {
            "pattern": "[-,-2,-3,-4]",
            "family": "NegIonizable",
        },
    }

    def __init__(self):
        """Initialize advanced structure processor."""
        super().__init__()
        self.logger = logging.getLogger(__name__)
        self._search_index = None
        self._scaler = StandardScaler()
        self._init_feature_factory()

    def _init_feature_factory(self):
        """Initialize pharmacophore feature factory."""
        try:
            from rdkit.Chem import ChemicalFeatures

            factory = ChemicalFeatures.BuildFeatureFactory()
            self._feature_factory = factory
        except Exception as e:
            self.logger.error(f"Error initializing feature factory: {str(e)}")
            self._feature_factory = None

    def analyze_scaffolds(
        self,
        mols: List[Chem.Mol],
        use_murcko: bool = True,
        generic: bool = False,
    ) -> Dict[str, List[int]]:
        """
        Group molecules by scaffold similarity.

        Args:
            mols: List of RDKit molecules
            use_murcko: Use Murcko scaffolds
            generic: Make scaffolds generic

        Returns:
            Dictionary mapping scaffold SMILES to molecule indices
        """
        try:
            scaffold_groups = {}
            for i, mol in enumerate(mols):
                if mol is None:
                    continue

                try:
                    # Generate scaffold
                    if use_murcko:
                        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                    else:
                        scaffold = Chem.Mol(mol)

                    if generic:
                        scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)

                    scaffold_smiles = Chem.MolToSmiles(scaffold)

                    if scaffold_smiles not in scaffold_groups:
                        scaffold_groups[scaffold_smiles] = []
                    scaffold_groups[scaffold_smiles].append(i)

                except Exception as e:
                    self.logger.warning(f"Error processing molecule {i}: {str(e)}")
                    continue

            return scaffold_groups

        except Exception as e:
            self.logger.error(f"Error analyzing scaffolds: {str(e)}")
            return {}

    def optimize_similarity_search(
        self,
        library_mols: List[Chem.Mol],
        fp_type: str = "morgan",
        n_bits: int = 2048,
        radius: int = 2,
    ) -> bool:
        """
        Build optimized search index for library.

        Args:
            library_mols: List of RDKit molecules
            fp_type: Fingerprint type
            n_bits: Number of bits in fingerprint
            radius: Radius for Morgan fingerprints

        Returns:
            Success status
        """
        try:
            # Generate fingerprints
            fps = []
            valid_indices = []

            for i, mol in enumerate(library_mols):
                if mol is None:
                    continue

                try:
                    if fp_type == "morgan":
                        fp = AllChem.GetMorganFingerprintAsBitVect(
                            mol, radius, nBits=n_bits
                        )
                    elif fp_type == "maccs":
                        fp = MACCSkeys.GenMACCSKeys(mol)
                    elif fp_type == "topological":
                        fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
                            mol, nBits=n_bits
                        )
                    else:
                        continue

                    fps.append(list(fp.GetOnBits()))
                    valid_indices.append(i)

                except Exception as e:
                    self.logger.warning(f"Error processing molecule {i}: {str(e)}")
                    continue

            if not fps:
                return False

            # Build search index
            X = np.array(fps)
            X = self._scaler.fit_transform(X)

            model = NearestNeighbors(
                n_neighbors=min(100, len(X)),
                algorithm="ball_tree",
                metric="jaccard",
                n_jobs=-1,
            )
            model.fit(X)

            self._search_index = {
                "model": model,
                "indices": valid_indices,
                "fp_type": fp_type,
                "n_bits": n_bits,
                "radius": radius,
            }

            return True

        except Exception as e:
            self.logger.error(f"Error optimizing search: {str(e)}")
            return False

    def find_similar_fast(
        self,
        query_mol: Chem.Mol,
        n_neighbors: int = 10,
        threshold: float = 0.7,
    ) -> List[Tuple[int, float]]:
        """
        Find similar compounds using optimized index.

        Args:
            query_mol: Query molecule
            n_neighbors: Number of neighbors to return
            threshold: Similarity threshold

        Returns:
            List of (compound_index, similarity) tuples
        """
        try:
            if query_mol is None or self._search_index is None:
                return []

            # Generate query fingerprint
            fp_type = self._search_index["fp_type"]
            n_bits = self._search_index["n_bits"]
            radius = self._search_index["radius"]

            if fp_type == "morgan":
                fp = AllChem.GetMorganFingerprintAsBitVect(
                    query_mol, radius, nBits=n_bits
                )
            elif fp_type == "maccs":
                fp = MACCSkeys.GenMACCSKeys(query_mol)
            elif fp_type == "topological":
                fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
                    query_mol, nBits=n_bits
                )
            else:
                return []

            # Search index
            query_fp = self._scaler.transform([list(fp.GetOnBits())])
            distances, indices = self._search_index["model"].kneighbors(
                query_fp, n_neighbors=n_neighbors
            )

            # Convert distances to similarities
            similarities = 1 / (1 + distances.flatten())

            # Filter by threshold
            results = []
            for idx, sim in zip(indices.flatten(), similarities):
                if sim >= threshold:
                    mol_idx = self._search_index["indices"][idx]
                    results.append((mol_idx, float(sim)))

            results.sort(key=lambda x: x[1], reverse=True)
            return results

        except Exception as e:
            self.logger.error(f"Error searching similar compounds: {str(e)}")
            return []

    def find_pharmacophore_matches(
        self,
        query_mol: Chem.Mol,
        library_mols: List[Chem.Mol],
        min_features: int = 3,
        conformer_rmsd: float = 2.0,
    ) -> List[Tuple[int, float]]:
        """
        Find molecules matching pharmacophore pattern.

        Args:
            query_mol: Query molecule
            library_mols: List of molecules to search
            min_features: Minimum number of matching features
            conformer_rmsd: RMSD threshold for 3D matching

        Returns:
            List of (compound_index, match_score) tuples
        """
        try:
            if query_mol is None or not library_mols or self._feature_factory is None:
                return []

            # Generate query pharmacophore
            query_features = self._get_pharmacophore_features(query_mol)
            if len(query_features) < min_features:
                return []

            # Match against library
            matches = []
            for i, mol in enumerate(library_mols):
                if mol is None:
                    continue

                try:
                    # Generate 3D conformer if needed
                    if not mol.GetNumConformers():
                        AllChem.EmbedMolecule(mol, randomSeed=42)
                        AllChem.MMFFOptimizeMolecule(mol)

                    features = self._get_pharmacophore_features(mol)
                    if len(features) >= min_features:
                        score = self._match_pharmacophores(
                            query_features,
                            features,
                            conformer_rmsd,
                        )
                        if score > 0:
                            matches.append((i, float(score)))

                except Exception as e:
                    self.logger.warning(f"Error processing molecule {i}: {str(e)}")
                    continue

            matches.sort(key=lambda x: x[1], reverse=True)
            return matches

        except Exception as e:
            self.logger.error(f"Error matching pharmacophores: {str(e)}")
            return []

    def _get_pharmacophore_features(self, mol: Chem.Mol) -> List[Dict]:
        """Extract pharmacophore features from molecule."""
        try:
            features = []
            raw_features = self._feature_factory.GetFeaturesForMol(mol)

            for feature in raw_features:
                features.append(
                    {
                        "family": feature.GetFamily(),
                        "type": feature.GetType(),
                        "position": feature.GetPos(),
                        "atoms": feature.GetAtomIds(),
                    }
                )

            return features

        except Exception as e:
            self.logger.error(f"Error extracting features: {str(e)}")
            return []

    def _match_pharmacophores(
        self,
        query_features: List[Dict],
        target_features: List[Dict],
        rmsd_threshold: float = 2.0,
    ) -> float:
        """Calculate pharmacophore matching score."""
        try:
            if not query_features or not target_features:
                return 0.0

            # Count matching feature types
            type_matches = 0
            for qf in query_features:
                for tf in target_features:
                    if qf["family"] == tf["family"]:
                        type_matches += 1
                        break

            # Calculate feature overlap score
            max_features = max(len(query_features), len(target_features))
            type_score = type_matches / max_features

            # Calculate 3D matching score if positions available
            position_score = 0.0
            if all("position" in f for f in query_features + target_features):
                matched_positions = 0
                for qf in query_features:
                    qpos = np.array(qf["position"])
                    for tf in target_features:
                        if qf["family"] != tf["family"]:
                            continue
                        tpos = np.array(tf["position"])
                        if np.linalg.norm(qpos - tpos) <= rmsd_threshold:
                            matched_positions += 1
                            break
                position_score = matched_positions / max_features

            # Combine scores (weight 3D matching more heavily)
            final_score = (type_score + 2 * position_score) / 3
            return final_score

        except Exception as e:
            self.logger.error(f"Error matching features: {str(e)}")
            return 0.0
