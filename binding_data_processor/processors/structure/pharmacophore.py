"""Pharmacophore feature detection functionality."""

import logging
from typing import Dict, List, Optional, Set, Tuple
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem import AllChem

logger = logging.getLogger(__name__)


class PharmacophoreDetector:
    """Handles pharmacophore feature detection."""

    # Enhanced SMARTS patterns for pharmacophore features
    FEATURE_PATTERNS = {
        # Hydrogen bond donors
        "hbd": [
            "[N,O,S;H1,H2]-[!$(*=[O,N,P,S])]",  # Classic HBD
            "[n,o,s;H1;+0]",  # Aromatic HBD
            "[O,S;H1;+0;!$(*-*=[O,N,P,S])]",  # Hydroxyl/thiol
            "[N;H2;+0]",  # Primary amine
            "[N;H1;+0][!$(*=[O,N,P,S])]",  # Secondary amine
            "[N!H0;!$(N[C,S,P]=O)]",  # NH, NH2, NH3
            "[O!H0;!$(O[C,S,P]=O)]",  # OH
            "[S!H0]",  # SH
        ],
        # Hydrogen bond acceptors
        "hba": [
            "[$([O,S;H0;v2]),$([O,S;-])]",  # O, S acceptors
            "[$([N;H0;v3,v4&+1]),$([N;-;v2])]",  # N acceptors
            "[$([O,S;H0;v2])]",  # Carbonyl O, S
            "[N;H0;+0;!$(*-*=[O,N,P,S])]",  # Tertiary amine
            "[O,S;H0;+0]",  # Ether, thioether
            "[N!$(N[C,S,P]=O);!$(N=O);!$(NC=O);!$(NS=O);!$(NP=O)]",  # N
            "[O!$(O[C,S,P]=O);!$(O=N);!$(OC=O);!$(OS=O);!$(OP=O)]",  # O
            "[S!$(S[C,P]=O);!$(S=O);!$(SC=O);!$(SP=O)]",  # S
        ],
        # Positive ionizable
        "pi": [
            "[N;H2,H3;+1]",  # Primary ammonium
            "[N;H2,H1;+1][C]",  # Secondary/tertiary ammonium
            "[N;H0;+1]=[C]",  # Iminium
            "[n;H1;+1]",  # Protonated aromatic N
            "[N;H0;+1]([C])[C]",  # Quaternary N
            "[N+]",  # Quaternary N
            "[NH2;X3;!$(NC=O);!$(NS=O);!$(NP=O)]",  # Primary amine
            "[NH1;X3;!$(NC=O)]",  # Secondary amine
            "[NH0;X3;+]",  # Tertiary amine (protonated)
            "[NH2;X4;+]",  # Primary ammonium
        ],
        # Negative ionizable
        "ni": [
            "[$([O-;!$(*-[N,P,S])]),$([S-;!$(*-[N,P])]),$([O,S;-;!$(*-*=[N,P])])]",  # Carboxylate, sulfonate
            "[C,S](=[O,S])[O-]",  # Carboxylate, sulfonate
            "[P](=[O])[O-]",  # Phosphate
            "[B-](F)(F)(F)[F-]",  # Tetrafluoroborate
            "[C,P,S](=O)[O-,OH]",  # Carboxylate, phosphate, sulfonate
            "[S,P](=O)(=O)[O-,OH]",  # Sulfate, phosphate
            "[O-;X1]",  # Oxide
            "[N-;X2]",  # Azide
        ],
        # Aromatic
        "ar": [
            "a1aaaaa1",  # 6-membered aromatic
            "a1aaaa1",  # 5-membered aromatic
            "[$(c1ccccc1),$(c1cccc1)]",  # Benzene, pyridine
            "a1:a:a:a:a:1",  # 6-membered aromatic
            "a1:a:a:a:1",  # 5-membered aromatic
        ],
        # Hydrophobic
        "hp": [
            "[C;!$(C=[O,N,S]);!$(C[F,Cl,Br,I,N,O,P,S])]",  # Alkyl
            "[C;H3,H2,H1][!N,!O,!S,!P]",  # Terminal alkyl
            "[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I]),$([!#1]);!$([CH3X4,CH2X3,CH1X2][!C]);!$(C(F)(F)F)]",  # General hydrophobic
            "[C;D3,D4;!$(C[N,O,S]);!$(C(N)(N));!$(C(N)(O));!$(C(N)(S))]",  # Alkyl
            "[C;$(C=C);!$(C[N,O,S])]",  # Alkene
            "[CH2;!$(CC[N,O,S])]",  # Methylene
            "[$([CH3][CH2][CH2][CH2]),$([CH3][CH2][CH2][CH3])]",  # Butyl+
        ],
    }

    def __init__(self):
        """Initialize pharmacophore detector."""
        self.logger = logging.getLogger(__name__)

        # Pre-compile SMARTS patterns
        self.feature_patterns = {
            feature: [Chem.MolFromSmarts(p) for p in patterns]
            for feature, patterns in self.FEATURE_PATTERNS.items()
        }

        # Initialize feature factory
        self.factory = ChemicalFeatures.BuildFeatureFactory()
        self._init_feature_factory()

    def _init_feature_factory(self):
        """Initialize chemical feature factory with custom definitions."""
        feature_defs = [
            # Hydrogen bond donors
            """DefineFeature HBD [N,O,S;H1,H2]-[!$(*=[O,N,P,S])]
                Family HBond
                Weights 1.0
            EndFeature""",
            # Hydrogen bond acceptors
            """DefineFeature HBA [$([O,S;H0;v2]),$([O,S;-])]
                Family HBond
                Weights 1.0
            EndFeature""",
            # Positive ionizable
            """DefineFeature PI [N;H2,H3;+1]
                Family PosIon
                Weights 1.0
            EndFeature""",
            # Negative ionizable
            """DefineFeature NI [C,S](=[O,S])[O-]
                Family NegIon
                Weights 1.0
            EndFeature""",
            # Aromatic
            """DefineFeature AR a1aaaaa1
                Family Aromatic
                Weights 1.0,1.0,1.0,1.0,1.0,1.0
            EndFeature""",
            # Hydrophobic
            """DefineFeature HP [C;!$(C=[O,N,S]);!$(C[F,Cl,Br,I,N,O,P,S])]
                Family Hydrophobe
                Weights 1.0
            EndFeature""",
        ]

        for fdef in feature_defs:
            self.factory.ParseFeatureDef(fdef)

    def detect_pharmacophores(
        self, mol: Chem.Mol, include_coords: bool = False
    ) -> Dict[str, List[Tuple[int, ...]]]:
        """
        Detect pharmacophore features in molecule.

        Args:
            mol: RDKit molecule
            include_coords: Include 3D coordinates of features

        Returns:
            Dictionary mapping feature types to lists of atom indices
        """
        try:
            if mol is None:
                return {}

            features = {}

            # Detect features using SMARTS patterns
            for feature, patterns in self.feature_patterns.items():
                matches = []
                for pattern in patterns:
                    if pattern is not None:
                        matches.extend(mol.GetSubstructMatches(pattern))
                if matches:
                    features[feature] = matches

            # Detect features using feature factory
            if include_coords:
                # Generate 3D conformation if needed
                if not mol.GetNumConformers():
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol, randomSeed=42)
                    AllChem.MMFFOptimizeMolecule(mol)

                factory_features = self.factory.GetFeaturesForMol(mol)
                for f in factory_features:
                    feature_type = f.GetFamily().lower()
                    if feature_type not in features:
                        features[feature_type] = []
                    features[feature_type].append((tuple(f.GetAtomIds()), f.GetPos()))

            # Generate 2D pharmacophore fingerprint
            try:
                sigfactory = self._create_sigfactory()
                fp = Generate.Gen2DFingerprint(mol, sigfactory)
                features["fingerprint"] = [
                    tuple([i]) for i, bit in enumerate(fp) if bit
                ]
            except Exception as e:
                self.logger.warning(
                    f"Error generating pharmacophore fingerprint: {str(e)}"
                )

            return features

        except Exception as e:
            self.logger.error(f"Error detecting pharmacophores: {str(e)}")
            return {}

    def _create_sigfactory(self) -> SigFactory:
        """Create pharmacophore signature factory."""
        factory = SigFactory()

        # Add feature definitions
        factory.SetPatternsFromSmarts(
            {
                "HBD": ["[N,O,S;H1,H2]-[!$(*=[O,N,P,S])]"],
                "HBA": ["[$([O,S;H0;v2]),$([O,S;-])]"],
                "PI": ["[N;H2,H3;+1]"],
                "NI": ["[C,S](=[O,S])[O-]"],
                "AR": ["a1aaaaa1"],
                "HP": ["[C;!$(C=[O,N,S]);!$(C[F,Cl,Br,I,N,O,P,S])]"],
            }
        )

        # Set distance bins (in Å)
        factory.SetBins([(0, 2), (2, 4), (4, 6), (6, 8), (8, 10)])

        # Set triangle distance bounds
        factory.SetTrianglePruneBounds(2, 10)

        factory.Init()
        return factory

    def get_3d_pharmacophores(
        self, mol: Chem.Mol, conf_id: int = -1
    ) -> Dict[str, List[Tuple[float, float, float]]]:
        """
        Get 3D coordinates of pharmacophore features.

        Args:
            mol: RDKit molecule with 3D coordinates
            conf_id: Conformer ID to use

        Returns:
            Dictionary mapping feature types to lists of 3D coordinates
        """
        try:
            if mol is None or not mol.GetNumConformers():
                return {}

            features_3d = {}
            conf = mol.GetConformer(conf_id)

            # Get 2D features first
            features_2d = self.detect_pharmacophores(mol)

            # Convert atom indices to 3D coordinates
            for feature_type, matches in features_2d.items():
                if feature_type != "fingerprint":  # Skip fingerprint bits
                    coords = []
                    for match in matches:
                        # Calculate centroid for multi-atom features
                        if len(match) > 1:
                            points = [conf.GetAtomPosition(i) for i in match]
                            centroid = (
                                sum(p.x for p in points) / len(points),
                                sum(p.y for p in points) / len(points),
                                sum(p.z for p in points) / len(points),
                            )
                            coords.append(centroid)
                        else:
                            pos = conf.GetAtomPosition(match[0])
                            coords.append((pos.x, pos.y, pos.z))
                    features_3d[feature_type] = coords

            return features_3d

        except Exception as e:
            self.logger.error(f"Error getting 3D pharmacophores: {str(e)}")
            return {}

    def get_feature_distances(
        self, mol: Chem.Mol, features: Dict[str, List[Tuple[int, ...]]]
    ) -> Dict[str, List[float]]:
        """
        Calculate distances between pharmacophore features.

        Args:
            mol: RDKit molecule
            features: Dictionary of feature atom indices

        Returns:
            Dictionary mapping feature pairs to distances
        """
        try:
            if mol is None or not features:
                return {}

            # Generate 3D conformation if needed
            if not mol.GetNumConformers():
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol)

            distances = {}
            conf = mol.GetConformer()

            # Calculate distances between all feature pairs
            feature_types = list(features.keys())
            for i, type1 in enumerate(feature_types):
                for type2 in feature_types[i:]:
                    if type1 != "fingerprint" and type2 != "fingerprint":
                        key = f"{type1}-{type2}"
                        distances[key] = []

                        for atoms1 in features[type1]:
                            for atoms2 in features[type2]:
                                # Use centroids for multi-atom features
                                pos1 = sum(
                                    conf.GetAtomPosition(a) for a in atoms1
                                ) / len(atoms1)
                                pos2 = sum(
                                    conf.GetAtomPosition(a) for a in atoms2
                                ) / len(atoms2)
                                dist = pos1.Distance(pos2)
                                distances[key].append(float(dist))

            return distances

        except Exception as e:
            self.logger.error(f"Error calculating feature distances: {str(e)}")
            return {}
