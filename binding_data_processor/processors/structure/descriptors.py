"""Molecular descriptor calculation with machine learning enhancements."""

import logging
from typing import Dict, List, Optional, Tuple
import numpy as np
from sklearn.ensemble import RandomForestRegressor, IsolationForest
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Fragments,
    Crippen,
    rdMolDescriptors,
    QED,
    rdRGroupDecomposition,
    MolSurf,
)
from rdkit.Chem.MolSurf import (
    PEOE_VSA1,
    PEOE_VSA2,
    PEOE_VSA3,
    SMR_VSA1,
    SMR_VSA2,
    SMR_VSA3,
    SlogP_VSA1,
    SlogP_VSA2,
    SlogP_VSA3,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

logger = logging.getLogger(__name__)


class MLDescriptorCalculator:
    """ML-enhanced molecular descriptor calculator."""

    def __init__(self):
        """Initialize ML models and descriptor calculator."""
        self.logger = logging.getLogger(__name__)
        self.rf_model = None
        self.scaler = StandardScaler()
        self.pca = PCA(n_components=0.95)  # Keep 95% of variance
        self.isolation_forest = IsolationForest(contamination=0.1)
        self.important_descriptors = None

    def fit(self, mols: List[Chem.Mol], activities: List[float]) -> None:
        """
        Fit ML models to training data.

        Args:
            mols: List of RDKit molecules
            activities: List of activity values
        """
        try:
            # Calculate all descriptors
            X = []
            for mol in mols:
                desc = self.calculate_all_descriptors(mol)
                if desc:
                    X.append(list(desc.values()))

            X = np.array(X)
            y = np.array(activities)

            # Scale features
            X_scaled = self.scaler.fit_transform(X)

            # Reduce dimensionality
            X_pca = self.pca.fit_transform(X_scaled)

            # Train random forest
            self.rf_model = RandomForestRegressor(n_estimators=100)
            self.rf_model.fit(X_pca, y)

            # Get feature importance
            importance = self.rf_model.feature_importances_
            self.important_descriptors = [desc for i, desc in enumerate(self.descriptor_names) if importance[i] > np.mean(importance)]

            # Train outlier detector
            self.isolation_forest.fit(X_scaled)

        except Exception as e:
            self.logger.error(f"Error fitting ML models: {str(e)}")

    def predict_activity(self, mol: Chem.Mol) -> Optional[float]:
        """
        Predict activity using trained model.

        Args:
            mol: RDKit molecule

        Returns:
            Predicted activity value
        """
        try:
            if self.rf_model is None:
                return None

            desc = self.calculate_all_descriptors(mol)
            if not desc:
                return None

            X = np.array(list(desc.values())).reshape(1, -1)
            X_scaled = self.scaler.transform(X)
            X_pca = self.pca.transform(X_scaled)

            return float(self.rf_model.predict(X_pca)[0])

        except Exception as e:
            self.logger.error(f"Error predicting activity: {str(e)}")
            return None

    def is_outlier(self, mol: Chem.Mol) -> bool:
        """
        Check if molecule is an outlier.

        Args:
            mol: RDKit molecule

        Returns:
            True if molecule is an outlier
        """
        try:
            desc = self.calculate_all_descriptors(mol)
            if not desc:
                return True

            X = np.array(list(desc.values())).reshape(1, -1)
            X_scaled = self.scaler.transform(X)

            return self.isolation_forest.predict(X_scaled)[0] == -1

        except Exception as e:
            self.logger.error(f"Error checking outlier: {str(e)}")
            return True

    def get_important_descriptors(self, mol: Chem.Mol) -> Dict[str, Tuple[float, float]]:
        """
        Get most important descriptors and their contributions.

        Args:
            mol: RDKit molecule

        Returns:
            Dictionary of descriptor names to (value, importance) tuples
        """
        try:
            if self.important_descriptors is None:
                return {}

            desc = self.calculate_all_descriptors(mol)
            if not desc:
                return {}

            important_desc = {}
            for name in self.important_descriptors:
                if name in desc:
                    value = desc[name]
                    importance = self.get_descriptor_importance(name)
                    important_desc[name] = (value, importance)

            return important_desc

        except Exception as e:
            self.logger.error(f"Error getting important descriptors: {str(e)}")
            return {}

    def get_descriptor_importance(self, name: str) -> float:
        """Get importance score for a descriptor."""
        try:
            if self.rf_model is None:
                return 0.0

            idx = self.descriptor_names.index(name)
            return float(self.rf_model.feature_importances_[idx])

        except Exception:
            return 0.0


class DescriptorCalculator(MLDescriptorCalculator):
    """Enhanced descriptor calculator with ML capabilities."""

    # SMARTS patterns for structural features
    STRUCTURE_PATTERNS = {
        # Basic functional groups
        "basic_amine": "[NX3;H2,H1;!$(NC=O)]",
        "acidic_oh": "[OH;$(O[#6]);!$(OC=O)]",
        "phenol": "[OH]c1[c,n][c,n][c,n][c,n][c,n]1",
        "carboxylic_acid": "[CX3](=O)[OX2H1]",
        "sulfonamide": "[SX4](=[OX1])(=[OX1])[NX3H2]",
        "amide": "[NX3][CX3](=[OX1])[#6]",
        "ester": "[#6][CX3](=O)[OX2H0][#6]",
        "phosphate": "[PX4](=[OX1])([OX2H1])([OX2H1])[OX2H1]",
        "sulfate": "[SX4](=[OX1])(=[OX1])([OX2H1])[OX2H1]",
        # Ring systems
        "aromatic_ring": "[a;r5,r6]1[a]:[a]:[a]:[a]:[a]1",
        "aliphatic_ring": "[A;R]1[A;R][A;R][A;R][A;R][A;R]1",
        "heterocycle": "[a;!c]",
        # Other features
        "halogen": "[F,Cl,Br,I]",
        "alkene": "[CX3]=[CX3]",
        "alkyne": "[CX2]#[CX2]",
        "ether": "[OX2]([#6])[#6]",
        "ketone": "[#6][CX3](=O)[#6]",
        "alcohol": "[OX2H][CX4]",
        "guanidine": "[N;H2]C(=N)N",
    }

    # Enhanced descriptor categories with descriptions
    DESCRIPTOR_CATEGORIES = {
        "physical": {
            "molecular_weight": Descriptors.ExactMolWt,
            "heavy_atom_count": Descriptors.HeavyAtomCount,
            "rotatable_bonds": Descriptors.NumRotatableBonds,
            "ring_count": Descriptors.RingCount,
            "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings,
            "aliphatic_rings": rdMolDescriptors.CalcNumAliphaticRings,
            "saturated_rings": rdMolDescriptors.CalcNumSaturatedRings,
            "heterocycles": rdMolDescriptors.CalcNumHeterocycles,
        },
        "topological": {
            "description": "2D structure-based descriptors",
            "descriptors": {
                "MolWt": (Descriptors.ExactMolWt, "Molecular weight"),
                "LogP": (Crippen.MolLogP, "Calculated LogP"),
                "TPSA": (Descriptors.TPSA, "Topological polar surface area"),
                "HBA": (rdMolDescriptors.CalcNumHBA, "H-bond acceptors"),
                "HBD": (rdMolDescriptors.CalcNumHBD, "H-bond donors"),
                "RotBonds": (rdMolDescriptors.CalcNumRotatableBonds, "Rotatable bonds"),
                "Rings": (rdMolDescriptors.CalcNumRings, "Number of rings"),
                "AromaticRings": (
                    rdMolDescriptors.CalcNumAromaticRings,
                    "Aromatic rings",
                ),
                "HeteroRings": (
                    rdMolDescriptors.CalcNumAliphaticHeterocycles,
                    "Heterocyclic rings",
                ),
                "BertzCT": (Descriptors.BertzCT, "Bertz complexity index"),
                "BalabanJ": (Descriptors.BalabanJ, "Balaban J index"),
                "WienerIndex": (lambda m: sum(sum(Chem.GetDistanceMatrix(m)[i][j] for j in range(i + 1, m.GetNumAtoms())) for i in range(m.GetNumAtoms())), "Wiener path index"),
            },
        },
        "connectivity": {
            "description": "Molecular connectivity descriptors",
            "descriptors": {
                "Chi0v": (rdMolDescriptors.CalcChi0v, "0th order valence chi"),
                "Chi1v": (rdMolDescriptors.CalcChi1v, "1st order valence chi"),
                "Chi2v": (rdMolDescriptors.CalcChi2v, "2nd order valence chi"),
                "Chi3v": (rdMolDescriptors.CalcChi3v, "3rd order valence chi"),
                "Chi4v": (rdMolDescriptors.CalcChi4v, "4th order valence chi"),
                "Kappa1": (rdMolDescriptors.CalcKappa1, "1st order kappa"),
                "Kappa2": (rdMolDescriptors.CalcKappa2, "2nd order kappa"),
                "Kappa3": (rdMolDescriptors.CalcKappa3, "3rd order kappa"),
                "HallKierAlpha": (Descriptors.HallKierAlpha, "Hall-Kier alpha value"),
            },
        },
        "electronic": {
            "description": "Electronic and charge-based descriptors",
            "descriptors": {
                "PEOE_VSA1": (PEOE_VSA1, "PE/OE charge VSA descriptor 1"),
                "PEOE_VSA2": (PEOE_VSA2, "PE/OE charge VSA descriptor 2"),
                "PEOE_VSA3": (PEOE_VSA3, "PE/OE charge VSA descriptor 3"),
                "SMR_VSA1": (SMR_VSA1, "MR VSA descriptor 1"),
                "SMR_VSA2": (SMR_VSA2, "MR VSA descriptor 2"),
                "SMR_VSA3": (SMR_VSA3, "MR VSA descriptor 3"),
                "SlogP_VSA1": (SlogP_VSA1, "LogP VSA descriptor 1"),
                "SlogP_VSA2": (SlogP_VSA2, "LogP VSA descriptor 2"),
                "SlogP_VSA3": (SlogP_VSA3, "LogP VSA descriptor 3"),
                "MolarRefractivity": (Crippen.MolMR, "Molar refractivity"),
                "MaxPartialCharge": (
                    Descriptors.MaxPartialCharge,
                    "Maximum partial charge",
                ),
                "MinPartialCharge": (
                    Descriptors.MinPartialCharge,
                    "Minimum partial charge",
                ),
            },
        },
        "surface": {
            "description": "Surface and volume descriptors",
            "descriptors": {
                "LabuteASA": (
                    rdMolDescriptors.CalcLabuteASA,
                    "Labute accessible surface area",
                ),
                "TPSA": (Descriptors.TPSA, "Topological polar surface area"),
                "MolVolume": (
                    lambda m: (
                        AllChem.ComputeMolVolume(m)
                        if m.GetNumConformers() > 0
                        else (AllChem.EmbedMolecule(m, randomSeed=42) != -1 and AllChem.MMFFOptimizeMolecule(m) != -1 and AllChem.ComputeMolVolume(m) or 0.0)
                    ),
                    "Molecular volume",
                ),
                "Asphericity": (
                    rdMolDescriptors.CalcAsphericity,
                    "Molecular asphericity",
                ),
                "Eccentricity": (
                    rdMolDescriptors.CalcEccentricity,
                    "Molecular eccentricity",
                ),
                "InertialShapeFactor": (
                    rdMolDescriptors.CalcInertialShapeFactor,
                    "Inertial shape factor",
                ),
                "NPR1": (
                    rdMolDescriptors.CalcNPR1,
                    "Normalized principal moments ratio 1",
                ),
                "NPR2": (
                    rdMolDescriptors.CalcNPR2,
                    "Normalized principal moments ratio 2",
                ),
            },
        },
        "fragment": {
            "description": "Fragment-based descriptors",
            "descriptors": {
                "fr_Al_OH": (Fragments.fr_Al_OH, "Aliphatic hydroxyl groups"),
                "fr_Ar_OH": (Fragments.fr_Ar_OH, "Aromatic hydroxyl groups"),
                "fr_NH2": (Fragments.fr_NH2, "Primary amines"),
                "fr_NH1": (Fragments.fr_NH1, "Secondary amines"),
                "fr_NH0": (Fragments.fr_NH0, "Tertiary amines"),
                "fr_Ar_N": (Fragments.fr_Ar_N, "Aromatic nitrogens"),
                "fr_Al_COO": (Fragments.fr_Al_COO, "Aliphatic carboxylic acids"),
                "fr_Ar_COO": (Fragments.fr_Ar_COO, "Aromatic carboxylic acids"),
                "fr_COO": (Fragments.fr_COO, "Esters"),
                "fr_ketone": (Fragments.fr_ketone, "Ketones"),
                "fr_ether": (Fragments.fr_ether, "Ethers"),
                "fr_phenol": (Fragments.fr_phenol, "Phenols"),
                "fr_aldehyde": (Fragments.fr_aldehyde, "Aldehydes"),
                "fr_amide": (Fragments.fr_amide, "Amides"),
            },
        },
    }

    def __init__(self):
        """Initialize descriptor calculator."""
        super().__init__()
        self.logger = logging.getLogger(__name__)
        # Pre-compile SMARTS patterns
        self.structure_patterns = {name: Chem.MolFromSmarts(smarts) for name, smarts in self.STRUCTURE_PATTERNS.items()}

    def calculate_descriptors(
        self,
        mol: Chem.Mol,
        categories: Optional[List[str]] = None,
        include_3d: bool = True,
        use_ml: bool = True,
    ) -> Dict[str, float]:
        """
        Calculate molecular descriptors with ML enhancements.

        Args:
            mol: RDKit molecule
            categories: Optional list of descriptor categories
            include_3d: Whether to include 3D descriptors
            use_ml: Whether to use ML enhancements

        Returns:
            Dictionary of descriptor names and values
        """
        try:
            # Calculate base descriptors
            results = super().calculate_descriptors(mol, categories, include_3d)

            if use_ml and self.rf_model is not None:
                # Add ML-based predictions
                activity_pred = self.predict_activity(mol)
                if activity_pred is not None:
                    results["predicted_activity"] = activity_pred

                # Add outlier score
                results["outlier_score"] = float(self.isolation_forest.score_samples(np.array(list(results.values())).reshape(1, -1))[0])

                # Add importance-weighted druglikeness
                druglike_score = self.get_druglikeness_score(mol)
                if druglike_score is not None:
                    imp_descriptors = self.get_important_descriptors(mol)
                    weighted_score = sum(druglike_score * imp[1] for imp in imp_descriptors.values()) / len(imp_descriptors)
                    results["ml_druglikeness"] = float(weighted_score)

            return results

        except Exception as e:
            self.logger.error(f"Error calculating descriptors: {str(e)}")
            return {}

    def calculate_descriptors_batch(
        self,
        mols: List[Chem.Mol],
        categories: Optional[List[str]] = None,
        names: Optional[List[str]] = None,
        use_ml: bool = True,
        n_jobs: int = -1,
    ) -> List[Dict[str, float]]:
        """
        Calculate descriptors for multiple molecules in parallel.

        Args:
            mols: List of RDKit molecules
            categories: Optional list of descriptor categories
            names: Optional list of specific descriptors
            use_ml: Whether to use ML enhancements
            n_jobs: Number of parallel jobs (-1 for all cores)

        Returns:
            List of descriptor dictionaries
        """
        try:
            from joblib import Parallel, delayed

            # Filter out None molecules
            valid_mols = [mol for mol in mols if mol is not None]
            if not valid_mols:
                return []

            # Calculate descriptors in parallel
            results = Parallel(n_jobs=n_jobs)(delayed(self.calculate_descriptors)(mol, categories, names, use_ml) for mol in valid_mols)

            # Add QED scores
            for i, mol in enumerate(valid_mols):
                try:
                    qed_score = QED.default(mol)
                    results[i]["qed_score"] = float(qed_score)
                except Exception as e:
                    self.logger.debug(f"Error calculating QED: {str(e)}")

            # Add decomposition analysis
            for i, mol in enumerate(valid_mols):
                try:
                    # Get Murcko scaffold
                    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                    results[i]["scaffold_smiles"] = Chem.MolToSmiles(scaffold)

                    # Get R-group decomposition
                    try:
                        rgroup_decomp = rdRGroupDecomposition.RGroupDecompositionParameters()
                        rgroup_decomp.removeHydrogensPostMatch = True
                        rgroup_decomp.onlyMatchAtRGroups = False
                        decomp = rdRGroupDecomposition.RGroupDecomposition(scaffold, rgroup_decomp)
                        decomp.Add(mol)
                        if decomp.Process():
                            rgroups = decomp.GetRGroupsAsColumns()
                            if rgroups:
                                results[i]["rgroup_count"] = len(rgroups)
                                for rgroup_label, rgroup_mols in rgroups.items():
                                    if rgroup_mols and rgroup_mols[0] is not None:
                                        results[i][f"rgroup_{rgroup_label}_smiles"] = Chem.MolToSmiles(rgroup_mols[0])
                    except Exception as e:
                        self.logger.debug(f"Error in R-group decomposition: {str(e)}")

                    # Get largest ring system
                    ring_info = mol.GetRingInfo()
                    rings = ring_info.AtomRings()
                    if rings:
                        max_ring = max(rings, key=len)
                        results[i]["largest_ring_size"] = len(max_ring)
                except Exception as e:
                    self.logger.debug(f"Error in decomposition: {str(e)}")

            return results

        except Exception as e:
            self.logger.error(f"Error calculating batch descriptors: {str(e)}")
            return []

    def get_druglikeness_score(self, mol: Chem.Mol, use_ml: bool = True) -> Optional[float]:
        """
        Calculate ML-enhanced druglikeness score.

        Args:
            mol: RDKit molecule
            use_ml: Whether to use ML weighting

        Returns:
            Druglikeness score between 0 and 1
        """
        try:
            base_score = super().get_druglikeness_score(mol)
            if base_score is None:
                return None

            if not use_ml or self.rf_model is None:
                return base_score

            # Get descriptor importance weights
            imp_descriptors = self.get_important_descriptors(mol)
            if not imp_descriptors:
                return base_score

            # Calculate weighted score
            weights = np.array([imp[1] for imp in imp_descriptors.values()])
            weighted_score = base_score * np.mean(weights)

            return float(weighted_score)

        except Exception as e:
            self.logger.error(f"Error calculating druglikeness: {str(e)}")
            return None

    def get_complexity_score(self, mol: Chem.Mol, use_ml: bool = True) -> Optional[float]:
        """
        Calculate ML-enhanced complexity score.

        Args:
            mol: RDKit molecule
            use_ml: Whether to use ML weighting

        Returns:
            Complexity score between 0 and 1
        """
        try:
            base_score = super().get_complexity_score(mol)
            if base_score is None:
                return None

            if not use_ml or self.rf_model is None:
                return base_score

            # Get descriptor importance weights
            imp_descriptors = self.get_important_descriptors(mol)
            if not imp_descriptors:
                return base_score

            # Calculate weighted score
            weights = np.array([imp[1] for imp in imp_descriptors.values()])
            weighted_score = base_score * np.mean(weights)

            return float(weighted_score)

        except Exception as e:
            self.logger.error(f"Error calculating complexity: {str(e)}")
            return None

    def _count_lipinski_violations(self, mol: Chem.Mol) -> int:
        """Count Lipinski's Rule of Five violations."""
        try:
            violations = 0
            mw = Descriptors.ExactMolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Descriptors.NumLipinskiHBD(mol)
            hba = Descriptors.NumLipinskiHBA(mol)
            rotatable = Descriptors.NumRotatableBonds(mol)

            if mw > 500:
                violations += 1
            if logp > 5:
                violations += 1
            if hbd > 5:
                violations += 1
            if hba > 10:
                violations += 1
            if rotatable > 10:
                violations += 1

            return violations

        except Exception:
            return 0

    def _count_ghose_violations(self, mol: Chem.Mol) -> int:
        """Count Ghose filter violations."""
        try:
            violations = 0
            mw = Descriptors.ExactMolWt(mol)
            logp = Crippen.MolLogP(mol)
            atoms = mol.GetNumAtoms()

            if not (160 <= mw <= 480):
                violations += 1
            if not (-0.4 <= logp <= 5.6):
                violations += 1
            if not (20 <= atoms <= 70):
                violations += 1

            return violations

        except Exception:
            return 0

    def _count_veber_violations(self, mol: Chem.Mol) -> int:
        """Count Veber filter violations."""
        try:
            violations = 0
            rotatable = Descriptors.NumRotatableBonds(mol)
            tpsa = Descriptors.TPSA(mol)

            if rotatable > 10:
                violations += 1
            if tpsa > 140:
                violations += 1

            return violations

        except Exception:
            return 0

    def calculate_descriptors(
        self,
        mol: Chem.Mol,
        categories: Optional[List[str]] = None,
        names: Optional[List[str]] = None,
    ) -> Dict[str, float]:
        """
        Calculate molecular descriptors.

        Args:
            mol: RDKit molecule
            categories: Optional list of descriptor categories to calculate
            names: Optional list of specific descriptor names to calculate

        Returns:
            Dictionary of descriptor names and values
        """
        try:
            if mol is None:
                return {}

            results = {}

            # Determine which descriptors to calculate
            if names:
                # Calculate specific descriptors
                for name in names:
                    for cat_info in self.DESCRIPTOR_CATEGORIES.values():
                        if name in cat_info["descriptors"]:
                            func = cat_info["descriptors"][name][0]
                            try:
                                value = func(mol)
                                if hasattr(value, "item"):  # Convert numpy types
                                    value = value.item()
                                results[name] = float(value)
                            except Exception as e:
                                self.logger.warning(f"Error calculating descriptor {name}: {str(e)}")
            else:
                # Calculate descriptors by category
                for cat_name, cat_info in self.DESCRIPTOR_CATEGORIES.items():
                    if categories is None or cat_name in categories:
                        for name, (func, _) in cat_info["descriptors"].items():
                            try:
                                value = func(mol)
                                if hasattr(value, "item"):  # Convert numpy types
                                    value = value.item()
                                results[name] = float(value)
                            except Exception as e:
                                self.logger.warning(f"Error calculating descriptor {name}: {str(e)}")

            # Add structural feature counts
            for name, pattern in self.structure_patterns.items():
                if pattern is not None:
                    results[f"count_{name}"] = len(mol.GetSubstructMatches(pattern))

            # Add fragment-based descriptors
            results.update(self.calculate_fragment_descriptors(mol))

            return results

        except Exception as e:
            self.logger.error(f"Error calculating descriptors: {str(e)}")
            return {}

    def calculate_descriptors_batch(
        self,
        mols: List[Chem.Mol],
        categories: Optional[List[str]] = None,
        names: Optional[List[str]] = None,
    ) -> List[Dict[str, float]]:
        """
        Calculate descriptors for multiple molecules.

        Args:
            mols: List of RDKit molecules
            categories: Optional list of descriptor categories
            names: Optional list of specific descriptors

        Returns:
            List of descriptor dictionaries
        """
        try:
            results = []
            for mol in mols:
                desc = self.calculate_descriptors(mol, categories, names)
                results.append(desc)
            return results

        except Exception as e:
            self.logger.error(f"Error calculating batch descriptors: {str(e)}")
            return []

    def calculate_fragment_descriptors(self, mol: Chem.Mol, include_rings: bool = True) -> Dict[str, float]:
        """
        Calculate fragment-based descriptors.

        Args:
            mol: RDKit molecule
            include_rings: Include ring fragments

        Returns:
            Dictionary of fragment counts
        """
        try:
            if mol is None:
                return {}

            fragment_counts = {}

            # Get all available fragment functions
            fragment_functions = [(name, func) for name, func in Fragments.__dict__.items() if name.startswith("fr_") and callable(func)]

            # Calculate all fragment counts
            for name, func in fragment_functions:
                try:
                    count = func(mol)
                    if isinstance(count, (int, float)):
                        fragment_counts[name] = float(count)
                except Exception as e:
                    self.logger.debug(f"Error calculating fragment {name}: {str(e)}")

            # Count ring systems if requested
            if include_rings:
                ring_info = mol.GetRingInfo()
                ring_counts = {}
                for ring in ring_info.AtomRings():
                    size = len(ring)
                    ring_counts[f"ring_{size}"] = float(ring_counts.get(f"ring_{size}", 0) + 1)
                fragment_counts.update(ring_counts)

            return fragment_counts

        except Exception as e:
            self.logger.error(f"Error calculating fragment descriptors: {str(e)}")
            return {}

    def get_descriptor_info(self, category: Optional[str] = None, name: Optional[str] = None) -> Dict:
        """
        Get information about available descriptors.

        Args:
            category: Optional category name to filter by
            name: Optional descriptor name to get specific info for

        Returns:
            Dictionary of descriptor information
        """
        try:
            if name:
                # Find specific descriptor
                for cat_info in self.DESCRIPTOR_CATEGORIES.values():
                    if name in cat_info["descriptors"]:
                        return {
                            "name": name,
                            "description": cat_info["descriptors"][name][1],
                        }
                return {}

            if category:
                # Return info for specific category
                if category in self.DESCRIPTOR_CATEGORIES:
                    return {
                        "category": category,
                        "description": self.DESCRIPTOR_CATEGORIES[category]["description"],
                        "descriptors": {name: desc[1] for name, desc in self.DESCRIPTOR_CATEGORIES[category]["descriptors"].items()},
                    }
                return {}

            # Return all descriptor information
            return {
                cat_name: {
                    "description": cat_info["description"],
                    "descriptors": {name: desc[1] for name, desc in cat_info["descriptors"].items()},
                }
                for cat_name, cat_info in self.DESCRIPTOR_CATEGORIES.items()
            }

        except Exception as e:
            self.logger.error(f"Error getting descriptor info: {str(e)}")
            return {}

    def get_druglikeness_score(self, mol: Chem.Mol, method: str = "combined") -> Optional[float]:
        """
        Calculate druglikeness score using multiple methods.

        Args:
            mol: RDKit molecule
            method: Scoring method ('combined', 'qed', 'rules', or 'ml')

        Returns:
            Druglikeness score between 0 and 1
        """
        try:
            if mol is None:
                return None

            scores = {}

            # Rules-based score
            if method in ["combined", "rules"]:
                mw = Descriptors.ExactMolWt(mol)
                logp = Crippen.MolLogP(mol)
                hba = rdMolDescriptors.CalcNumHBA(mol)
                hbd = rdMolDescriptors.CalcNumHBD(mol)
                tpsa = Descriptors.TPSA(mol)
                rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)

                # Score components (normalized between 0 and 1)
                scores["rules"] = float(
                    sum(
                        [
                            max(0, 1 - abs(mw - 400) / 400) * 1.5,
                            max(0, 1 - abs(logp - 2.5) / 5) * 1.5,
                            max(0, 1 - hba / 10),
                            max(0, 1 - hbd / 5),
                            max(0, 1 - abs(tpsa - 100) / 100),
                            max(0, 1 - rotatable / 10),
                        ]
                    )
                    / 7.0
                )  # Normalize by total weight

            # QED score
            if method in ["combined", "qed"]:
                try:
                    scores["qed"] = float(QED.default(mol))
                except Exception:
                    scores["qed"] = 0.0

            # ML-based score
            if method in ["combined", "ml"] and self.rf_model is not None:
                try:
                    desc = self.calculate_descriptors(mol, use_ml=False)
                    X = np.array(list(desc.values())).reshape(1, -1)
                    X_scaled = self.scaler.transform(X)
                    X_pca = self.pca.transform(X_scaled)
                    scores["ml"] = float(self.rf_model.predict(X_pca)[0])
                except Exception:
                    scores["ml"] = 0.0

            # Return appropriate score
            if method == "combined":
                weights = {"rules": 0.4, "qed": 0.4, "ml": 0.2}
                return float(sum(scores.get(k, 0.0) * w for k, w in weights.items()))
            return scores.get(method, 0.0)

        except Exception as e:
            self.logger.error(f"Error calculating druglikeness score: {str(e)}")
            return None

    def get_complexity_score(self, mol: Chem.Mol) -> Optional[float]:
        """
        Calculate molecular complexity score.

        Args:
            mol: RDKit molecule

        Returns:
            Complexity score between 0 and 1
        """
        try:
            if mol is None:
                return None

            # Calculate complexity-related descriptors
            num_atoms = mol.GetNumAtoms()
            num_bonds = mol.GetNumBonds()
            num_rings = rdMolDescriptors.CalcNumRings(mol)
            num_aromatic = rdMolDescriptors.CalcNumAromaticRings(mol)
            num_hetero = rdMolDescriptors.CalcNumAliphaticHeterocycles(mol)
            num_stereo = len(Chem.FindMolChiralCenters(mol))

            # Normalize and combine scores
            atom_score = min(1.0, num_atoms / 50)
            bond_score = min(1.0, num_bonds / 60)
            ring_score = min(1.0, (num_rings + num_aromatic) / 8)
            hetero_score = min(1.0, num_hetero / 4)
            stereo_score = min(1.0, num_stereo / 6)

            # Weighted average
            weights = [1.0, 1.0, 1.5, 1.2, 1.3]
            total_weight = sum(weights)
            score = sum(
                [
                    atom_score * weights[0],
                    bond_score * weights[1],
                    ring_score * weights[2],
                    hetero_score * weights[3],
                    stereo_score * weights[4],
                ]
            )
            return float(score / total_weight)

        except Exception as e:
            self.logger.error(f"Error calculating complexity score: {str(e)}")
            return None
