"""Pharmacophore feature detection functionality with enhanced RDKit features."""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union, Any
from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures, Draw, AllChem, rdDepictor, rdMolTransforms
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
import numpy as np
from rdkit.Chem.Draw import rdDepictor
from rdkit.Geometry import Point3D
import copy

logger = logging.getLogger(__name__)


@dataclass
class PharmacophoreFeature:
    """Class representing a pharmacophore feature."""

    feature_type: str
    atoms: Tuple[int, ...]
    position: Optional[Tuple[float, float, float]] = None
    radius: float = 1.0
    enabled: bool = True
    weight: float = 1.0
    vector: Optional[Tuple[float, float, float]] = None


class PharmacophoreGenerator:
    """Handles pharmacophore feature detection with enhanced RDKit functionality."""

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
        # Metal binding
        "metal": [
            "[O,N,S;!H0;+0]",  # Metal coordination sites
            "[O,N,S;-]",  # Anionic metal binding
            "[$(O=[C,S,P]);!$([OH,SH,NH])]",  # Carbonyl metal binding
        ],
        # Halogen bonds
        "xbond": [
            "[F,Cl,Br,I]",  # Halogen bond donors
            "[N,O,S;!H0]",  # Halogen bond acceptors
        ],
        # New patterns for additional features
        "stereocenter": [
            "[!H0;$(C([*])([*])([*])[*])]",  # Tetrahedral carbon
            "[!H0;$(N([*])([*])([*]))]",  # Tetrahedral nitrogen
            "[!H0;$(P([*])([*])([*]))]",  # Tetrahedral phosphorus
            "[!H0;$(S([*])([*])([*]))]",  # Tetrahedral sulfur
        ],
        "conjugated": [
            "[$(C=C-C=C),$(C=C-C=N),$(C=C-N=N)]",  # Conjugated double bonds
            "[$(c:c:c:c:c)]",  # Conjugated aromatic systems
            "[$(C=C-C#N),$(C=C-C=O)]",  # Conjugated with electron-withdrawing groups
        ],
        "sulfonamide": [
            "[$([#16X4]([NX3H])([!O])[OX2H0])]",  # Sulfonamide group
        ],
        "phosphate": [
            "[$([PX4](=[OX1])([OX2H0])([OX2H0])[OX2H0])]",  # Phosphate group
        ],
        "guanidine": [
            "[$([NX3H2][CX3](=[NX2H0])[NX3H2])]",  # Guanidine group
        ],
    }

    # Color schemes for visualization
    COLOR_SCHEMES = {
        "default": {
            "hbd": (1, 0, 0),  # Red
            "hba": (0, 0, 1),  # Blue
            "pi": (0, 1, 0),  # Green
            "ni": (1, 0.5, 0),  # Orange
            "ar": (0.5, 0, 0.5),  # Purple
            "hp": (0.5, 0.5, 0.5),  # Gray
            "metal": (1, 1, 0),  # Yellow
            "xbond": (0, 1, 1),  # Cyan
            "stereocenter": (1, 0, 1),  # Magenta
            "conjugated": (0.7, 0.7, 0),  # Yellow-green
            "sulfonamide": (0.5, 0.2, 0.8),  # Purple-blue
            "phosphate": (0.8, 0.4, 0.2),  # Orange-brown
            "guanidine": (0.2, 0.8, 0.4),  # Blue-green
        },
        "colorblind": {
            # Colorblind-friendly palette
            "hbd": (0.9, 0.6, 0),  # Orange
            "hba": (0.35, 0.7, 0.9),  # Light blue
            "pi": (0, 0.6, 0.5),  # Teal
            "ni": (0.95, 0.9, 0.25),  # Yellow
            "ar": (0.8, 0.4, 0),  # Brown
            "hp": (0.8, 0.8, 0.8),  # Gray
            "metal": (0.35, 0.35, 0.35),  # Dark gray
            "xbond": (0.9, 0.6, 0.6),  # Pink
            "stereocenter": (0.5, 0.5, 0.8),  # Purple
            "conjugated": (0.4, 0.7, 0.4),  # Green
            "sulfonamide": (0.7, 0.5, 0.8),  # Light purple
            "phosphate": (0.6, 0.6, 0.2),  # Olive
            "guanidine": (0.2, 0.7, 0.6),  # Turquoise
        },
    }

    def __init__(self):
        """Initialize pharmacophore generator with enhanced features."""
        self.logger = logging.getLogger(__name__)

        # Pre-compile SMARTS patterns
        self.feature_patterns = {feature: [Chem.MolFromSmarts(p) for p in patterns] for feature, patterns in self.FEATURE_PATTERNS.items()}

        # Initialize feature factory with enhanced definitions
        self.factory = ChemicalFeatures.BuildFeatureFactory()
        self._init_feature_factory()

        # Default settings
        self.conformer_settings = {
            "num_confs": 10,
            "num_threads": 0,  # Use all available CPUs
            "random_seed": 42,
            "pruneRmsThresh": 0.5,  # Remove similar conformers
            "enforceChirality": True,
            "useExpTorsionAnglePrefs": True,
            "useBasicKnowledge": True,
            "ETversion": 2,  # Use newer version of experimental torsion preferences
        }

    def _init_feature_factory(self):
        """Initialize chemical feature factory with enhanced definitions."""
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
            # Metal binding
            """DefineFeature MB [O,N,S;!H0;+0]
                Family MetalBinding
                Weights 1.0
            EndFeature""",
            # Halogen bonds
            """DefineFeature XB [F,Cl,Br,I]
                Family HalogenBond
                Weights 1.0
            EndFeature""",
        ]

        # Add new feature definitions
        additional_defs = [
            # Stereocenter definition
            """DefineFeature SC [!H0;$(C([*])([*])([*])[*])]
                Family Stereocenter
                Weights 1.0
            EndFeature""",
            # Conjugated system definition
            """DefineFeature CJ [$(C=C-C=C),$(C=C-C=N),$(C=C-N=N)]
                Family Conjugated
                Weights 1.0,1.0,1.0,1.0
            EndFeature""",
            # Sulfonamide definition
            """DefineFeature SU [$([#16X4]([NX3H])([!O])[OX2H0])]
                Family Sulfonamide
                Weights 1.0
            EndFeature""",
            # Phosphate definition
            """DefineFeature PH [$([PX4](=[OX1])([OX2H0])([OX2H0])[OX2H0])]
                Family Phosphate
                Weights 1.0
            EndFeature""",
            # Guanidine definition
            """DefineFeature GU [$([NX3H2][CX3](=[NX2H0])[NX3H2])]
                Family Guanidine
                Weights 1.0
            EndFeature""",
        ]

        for fdef in additional_defs:
            self.factory.ParseFeatureDef(fdef)

    def prepare_molecule(self, mol: Chem.Mol, generate_conformers: bool = True) -> Optional[Chem.Mol]:
        """
        Prepare molecule for pharmacophore generation with enhanced 3D handling.

        Args:
            mol: Input RDKit molecule
            generate_conformers: Whether to generate multiple conformers

        Returns:
            Prepared molecule or None if preparation fails
        """
        try:
            if mol is None:
                return None

            # Make a copy to avoid modifying input
            mol = Chem.Mol(mol)

            # Clean up structure
            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)

            if generate_conformers:
                # Generate multiple conformers
                params = AllChem.ETKDGv3()
                for key, value in self.conformer_settings.items():
                    if hasattr(params, key):
                        setattr(params, key, value)

                # Generate conformers
                cids = AllChem.EmbedMultipleConfs(
                    mol,
                    numConfs=self.conformer_settings["num_confs"],
                    params=params,
                    numThreads=self.conformer_settings["num_threads"],
                    randomSeed=self.conformer_settings["random_seed"],
                )

                # Optimize all conformers
                for cid in cids:
                    AllChem.MMFFOptimizeMolecule(mol, confId=cid)

                # Prune similar conformers
                if len(cids) > 1:
                    AllChem.AlignMolConformers(mol)
                    rms_list = []
                    for i in range(len(cids)):
                        for j in range(i + 1, len(cids)):
                            rms = AllChem.GetConformerRMS(mol, i, j)
                            rms_list.append((rms, i, j))

                    # Sort by RMS values
                    rms_list.sort()

                    # Remove similar conformers
                    to_remove = set()
                    for rms, i, j in rms_list:
                        if rms < self.conformer_settings["pruneRmsThresh"]:
                            if j not in to_remove:
                                to_remove.add(j)

                    # Create new molecule with unique conformers
                    new_mol = Chem.Mol(mol)
                    new_mol.RemoveAllConformers()
                    for i in range(mol.GetNumConformers()):
                        if i not in to_remove:
                            conf = mol.GetConformer(i)
                            new_mol.AddConformer(conf, assignId=True)
                    mol = new_mol

            return mol

        except Exception as e:
            self.logger.error(f"Error preparing molecule: {str(e)}")
            return None

    def generate_features(self, mol: Chem.Mol, include_factory_features: bool = True, include_vectors: bool = True, confId: int = -1) -> Dict[str, List[PharmacophoreFeature]]:
        """
        Generate pharmacophore features for a molecule with enhanced functionality.

        Args:
            mol: RDKit molecule
            include_factory_features: Whether to include features from feature factory
            include_vectors: Whether to calculate feature vectors (for directed features)
            confId: Conformer ID to use for 3D coordinates

        Returns:
            Dictionary mapping feature types to lists of PharmacophoreFeature objects
        """
        try:
            if mol is None:
                return {}

            features = {}
            conf = mol.GetConformer(confId) if mol.GetNumConformers() > 0 else None

            # Detect features using SMARTS patterns
            for feature_type, patterns in self.feature_patterns.items():
                features[feature_type] = []
                for pattern in patterns:
                    if pattern is not None:
                        matches = mol.GetSubstructMatches(pattern)
                        for match in matches:
                            # Calculate feature position (centroid of matched atoms)
                            position = None
                            if conf is not None:
                                positions = [conf.GetAtomPosition(i) for i in match]
                                position = (
                                    sum(p.x for p in positions) / len(positions),
                                    sum(p.y for p in positions) / len(positions),
                                    sum(p.z for p in positions) / len(positions),
                                )

                            # Calculate feature vector for directed features
                            vector = None
                            if include_vectors and conf is not None:
                                if feature_type in ["hbd", "hba"]:
                                    # Calculate vector based on donor/acceptor direction
                                    if len(match) >= 2:
                                        p1 = conf.GetAtomPosition(match[0])
                                        p2 = conf.GetAtomPosition(match[1])
                                        vector = (p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
                                        # Normalize vector
                                        magnitude = (vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2) ** 0.5
                                        if magnitude > 0:
                                            vector = tuple(v / magnitude for v in vector)

                            feature = PharmacophoreFeature(feature_type=feature_type, atoms=match, position=position, vector=vector)
                            features[feature_type].append(feature)

            # Add ring centers as features
            if mol.GetNumConformers() > 0:
                features["ring_centers"] = []
                ring_info = mol.GetRingInfo()
                for ring in ring_info.AtomRings():
                    if len(ring) > 0:
                        # Calculate ring center and normal vector
                        positions = [conf.GetAtomPosition(i) for i in ring]
                        center = (
                            sum(p.x for p in positions) / len(positions),
                            sum(p.y for p in positions) / len(positions),
                            sum(p.z for p in positions) / len(positions),
                        )

                        # Calculate ring normal vector using cross products
                        if len(positions) >= 3:
                            v1 = Point3D(positions[1].x - positions[0].x, positions[1].y - positions[0].y, positions[1].z - positions[0].z)
                            v2 = Point3D(positions[2].x - positions[0].x, positions[2].y - positions[0].y, positions[2].z - positions[0].z)
                            normal = v1.CrossProduct(v2)
                            magnitude = (normal.x**2 + normal.y**2 + normal.z**2) ** 0.5
                            if magnitude > 0:
                                normal_vector = (normal.x / magnitude, normal.y / magnitude, normal.z / magnitude)
                            else:
                                normal_vector = None
                        else:
                            normal_vector = None

                        feature = PharmacophoreFeature(feature_type="ring_centers", atoms=ring, position=center, vector=normal_vector)
                        features["ring_centers"].append(feature)

            return features

        except Exception as e:
            self.logger.error(f"Error detecting pharmacophores: {str(e)}")
            return {}

    def align_pharmacophores(
        self,
        ref_mol: Chem.Mol,
        probe_mol: Chem.Mol,
        ref_features: Optional[Dict[str, List[PharmacophoreFeature]]] = None,
        probe_features: Optional[Dict[str, List[PharmacophoreFeature]]] = None,
        max_combinations: int = 1000,
    ) -> Tuple[float, Chem.Mol]:
        """
        Align two molecules based on their pharmacophore features.

        Args:
            ref_mol: Reference molecule
            probe_mol: Probe molecule to align
            ref_features: Pre-computed features for reference molecule
            probe_features: Pre-computed features for probe molecule
            max_combinations: Maximum number of feature combinations to try

        Returns:
            Tuple of (RMSD score, aligned probe molecule)
        """
        try:
            # Generate features if not provided
            if ref_features is None:
                ref_features = self.generate_features(ref_mol)
            if probe_features is None:
                probe_features = self.generate_features(probe_mol)

            # Get all feature positions
            ref_positions = []
            probe_positions = []
            ref_weights = []

            for feat_type in ref_features:
                for ref_feat in ref_features[feat_type]:
                    if ref_feat.position is not None:
                        for probe_feat in probe_features.get(feat_type, []):
                            if probe_feat.position is not None:
                                ref_positions.append(ref_feat.position)
                                probe_positions.append(probe_feat.position)
                                # Weight by feature importance
                                ref_weights.append(ref_feat.weight)

            if len(ref_positions) < 3:
                raise ValueError("Need at least 3 matching features for alignment")

            # Convert to numpy arrays
            ref_positions = np.array(ref_positions)
            probe_positions = np.array(probe_positions)
            weights = np.array(ref_weights)

            # Center the coordinates
            ref_centroid = np.average(ref_positions, weights=weights, axis=0)
            probe_centroid = np.average(probe_positions, weights=weights, axis=0)

            ref_centered = ref_positions - ref_centroid
            probe_centered = probe_positions - probe_centroid

            # Calculate optimal rotation matrix
            covariance = np.dot(probe_centered.T * weights[:, None], ref_centered)
            V, S, Wt = np.linalg.svd(covariance)

            # Ensure right-handed coordinate system
            d = np.linalg.det(np.dot(V, Wt))
            if d < 0:
                V[:, -1] *= -1

            # Calculate rotation matrix
            rotation = np.dot(V, Wt)

            # Apply transformation to probe molecule
            aligned_mol = Chem.Mol(probe_mol)
            conf = aligned_mol.GetConformer()

            for i in range(aligned_mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                pos_array = np.array([pos.x, pos.y, pos.z])

                # Apply alignment transformation
                new_pos = np.dot(pos_array - probe_centroid, rotation.T) + ref_centroid

                conf.SetAtomPosition(i, Point3D(*new_pos))

            # Calculate RMSD
            rmsd = np.sqrt(np.mean(np.sum((np.dot(probe_centered, rotation.T) - ref_centered) ** 2, axis=1)))

            return rmsd, aligned_mol

        except Exception as e:
            self.logger.error(f"Error aligning pharmacophores: {str(e)}")
            return float("inf"), probe_mol

    def visualize_features(
        self,
        mol: Chem.Mol,
        features: Dict[str, List[PharmacophoreFeature]],
        size: Tuple[int, int] = (400, 400),
        color_scheme: str = "default",
        show_vectors: bool = True,
        vector_scale: float = 1.0,
        show_labels: bool = True,
        highlight_features: Optional[List[str]] = None,
    ) -> Optional[Draw.Image]:
        """
        Generate an enhanced 2D depiction of the molecule with highlighted pharmacophore features.

        Args:
            mol: RDKit molecule
            features: Dictionary of pharmacophore features
            size: Image size as (width, height)
            color_scheme: Color scheme to use ("default" or "colorblind")
            show_vectors: Whether to show feature vectors
            vector_scale: Scale factor for vector arrows
            show_labels: Whether to show feature labels
            highlight_features: List of feature types to highlight (None for all)

        Returns:
            PIL Image object or None if visualization fails
        """
        try:
            if mol is None or not features:
                return None

            # Make a copy for drawing
            mol = Chem.Mol(mol)

            # Generate 2D coordinates if needed
            if not mol.GetNumConformers():
                rdDepictor.Compute2DCoords(mol)

            # Get color scheme
            colors = self.COLOR_SCHEMES.get(color_scheme, self.COLOR_SCHEMES["default"])

            # Collect atoms to highlight for each feature type
            highlights = {}
            for feat_type, feat_list in features.items():
                if highlight_features is None or feat_type in highlight_features:
                    highlights[feat_type] = set()
                    for feat in feat_list:
                        highlights[feat_type].update(feat.atoms)

            # Create drawing options
            draw_opts = Draw.DrawingOptions()
            draw_opts.includeAtomNumbers = show_labels
            draw_opts.bondLineWidth = 2

            # Generate depiction with all features
            img = Draw.MolToImage(
                mol,
                size=size,
                highlightAtoms=[atom for atoms in highlights.values() for atom in atoms],
                highlightBonds=[],
                highlightColor=colors.get(next(iter(highlights.keys())), (0.7, 0.7, 0.7)),
                drawOptions=draw_opts,
            )

            # Add feature vectors if requested
            if show_vectors:
                drawer = Draw.rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
                drawer.drawOptions().prepareMolsBeforeDrawing = False

                # Draw the molecule first
                drawer.DrawMolecule(mol)

                # Draw feature vectors
                for feat_type, feat_list in features.items():
                    if highlight_features is None or feat_type in highlight_features:
                        for feat in feat_list:
                            if feat.position is not None and feat.vector is not None:
                                # Convert 3D vector to 2D for visualization
                                start = Point3D(*feat.position)
                                end = Point3D(
                                    feat.position[0] + feat.vector[0] * vector_scale,
                                    feat.position[1] + feat.vector[1] * vector_scale,
                                    feat.position[2] + feat.vector[2] * vector_scale,
                                )

                                # Project to 2D
                                start_2d = drawer.GetDrawCoords(start)
                                end_2d = drawer.GetDrawCoords(end)

                                # Draw arrow
                                color = colors.get(feat_type, (0.7, 0.7, 0.7))
                                drawer.SetColour(color)
                                drawer.DrawArrow(start_2d, end_2d)

                drawer.FinishDrawing()
                img = drawer.GetDrawingText()

            return img

        except Exception as e:
            self.logger.error(f"Error visualizing features: {str(e)}")
            return None

    def get_feature_constraints(self, features: Dict[str, List[PharmacophoreFeature]], tolerance: float = 1.0) -> List[Dict[str, Any]]:
        """
        Generate pharmacophore constraints from features.

        Args:
            features: Dictionary of pharmacophore features
            tolerance: Distance tolerance for constraints

        Returns:
            List of constraint dictionaries
        """
        constraints = []

        try:
            # Generate distance constraints between features
            feature_list = []
            for feat_type, feat_items in features.items():
                for feat in feat_items:
                    if feat.position is not None:
                        feature_list.append((feat_type, feat))

            for i, (type1, feat1) in enumerate(feature_list):
                for j, (type2, feat2) in enumerate(feature_list[i + 1 :], i + 1):
                    if feat1.position and feat2.position:
                        # Calculate distance
                        dist = sum((a - b) ** 2 for a, b in zip(feat1.position, feat2.position)) ** 0.5

                        constraint = {
                            "type": "distance",
                            "features": [{"type": type1, "atoms": feat1.atoms}, {"type": type2, "atoms": feat2.atoms}],
                            "min_dist": max(0.0, dist - tolerance),
                            "max_dist": dist + tolerance,
                            "weight": min(feat1.weight, feat2.weight),
                        }
                        constraints.append(constraint)

            # Generate angle constraints for directed features
            for i, (type1, feat1) in enumerate(feature_list):
                if feat1.vector:
                    for j, (type2, feat2) in enumerate(feature_list):
                        if i != j and feat2.vector:
                            # Calculate angle between vectors
                            dot_product = sum(a * b for a, b in zip(feat1.vector, feat2.vector))
                            angle = np.arccos(np.clip(dot_product, -1.0, 1.0))
                            angle_deg = np.degrees(angle)

                            constraint = {
                                "type": "angle",
                                "features": [{"type": type1, "atoms": feat1.atoms}, {"type": type2, "atoms": feat2.atoms}],
                                "min_angle": max(0.0, angle_deg - tolerance),
                                "max_angle": min(180.0, angle_deg + tolerance),
                                "weight": min(feat1.weight, feat2.weight),
                            }
                            constraints.append(constraint)

            return constraints

        except Exception as e:
            self.logger.error(f"Error generating constraints: {str(e)}")
            return []

    def score_alignment(
        self,
        ref_features: Dict[str, List[PharmacophoreFeature]],
        probe_features: Dict[str, List[PharmacophoreFeature]],
        distance_weight: float = 1.0,
        vector_weight: float = 0.5,
        feature_weights: Optional[Dict[str, float]] = None,
    ) -> float:
        """
        Score the alignment between two sets of pharmacophore features.

        Args:
            ref_features: Reference features
            probe_features: Probe features to compare
            distance_weight: Weight for distance contributions
            vector_weight: Weight for vector alignment contributions
            feature_weights: Optional weights for different feature types

        Returns:
            Alignment score (lower is better)
        """
        try:
            if not ref_features or not probe_features:
                return float("inf")

            total_score = 0.0
            total_weight = 0.0

            # Use default weights if none provided
            if feature_weights is None:
                feature_weights = {feat_type: 1.0 for feat_type in ref_features.keys()}

            # Score each feature type
            for feat_type in ref_features:
                if feat_type in probe_features:
                    type_weight = feature_weights.get(feat_type, 1.0)

                    ref_list = ref_features[feat_type]
                    probe_list = probe_features[feat_type]

                    # Find best matches between features
                    for ref_feat in ref_list:
                        best_score = float("inf")

                        for probe_feat in probe_list:
                            score = 0.0

                            # Distance score
                            if ref_feat.position and probe_feat.position:
                                dist = sum((a - b) ** 2 for a, b in zip(ref_feat.position, probe_feat.position)) ** 0.5
                                score += dist * distance_weight

                            # Vector alignment score
                            if ref_feat.vector and probe_feat.vector:
                                dot_product = sum(a * b for a, b in zip(ref_feat.vector, probe_feat.vector))
                                angle = np.arccos(np.clip(dot_product, -1.0, 1.0))
                                score += angle * vector_weight

                            best_score = min(best_score, score)

                        if best_score < float("inf"):
                            total_score += best_score * type_weight
                            total_weight += type_weight

            return total_score / total_weight if total_weight > 0 else float("inf")

        except Exception as e:
            self.logger.error(f"Error scoring alignment: {str(e)}")
            return float("inf")

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
                "MB": ["[O,N,S;!H0;+0]"],
                "XB": ["[F,Cl,Br,I]"],
            }
        )

        # Set distance bins (in Å)
        factory.SetBins([(0, 2), (2, 4), (4, 6), (6, 8), (8, 10)])

        # Set triangle distance bounds
        factory.SetTrianglePruneBounds(2, 10)

        factory.Init()
        return factory

    def get_3d_features(self, mol: Chem.Mol, confId: int = -1) -> Dict[str, List[Tuple[float, float, float]]]:
        """
        Get 3D coordinates of pharmacophore features.

        Args:
            mol: RDKit molecule with 3D conformer
            confId: Conformer ID to use

        Returns:
            Dictionary mapping feature types to lists of 3D coordinates
        """
        try:
            if mol is None or not mol.GetNumConformers():
                return {}

            features = self.generate_features(mol)
            coords = {}
            conf = mol.GetConformer(confId)

            # Get coordinates for each feature
            for name, matches in features.items():
                if name in ["fingerprint", "gobbi_fingerprint"]:
                    continue  # Skip fingerprint bits

                coords[name] = []
                for match in matches:
                    if len(match) > 1:
                        # Calculate centroid for multi-atom features
                        points = [conf.GetAtomPosition(i) for i in match]
                        centroid = (
                            sum(p.x for p in points) / len(points),
                            sum(p.y for p in points) / len(points),
                            sum(p.z for p in points) / len(points),
                        )
                        coords[name].append(centroid)
                    else:
                        # Use atom position for single-atom features
                        pos = conf.GetAtomPosition(match[0])
                        coords[name].append((pos.x, pos.y, pos.z))

            return coords

        except Exception as e:
            self.logger.error(f"Error getting 3D features: {str(e)}")
            return {}

    def get_feature_counts(self, mol: Chem.Mol) -> Dict[str, int]:
        """
        Get counts of each pharmacophore feature type.

        Args:
            mol: RDKit molecule

        Returns:
            Dictionary mapping feature types to counts
        """
        features = self.generate_features(mol)
        return {name: len(matches) for name, matches in features.items() if name not in ["fingerprint", "gobbi_fingerprint"]}

    def get_feature_distances(self, mol: Chem.Mol, features: Dict[str, List[Tuple[int, ...]]]) -> Dict[str, List[float]]:
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
                mol = self.prepare_molecule(mol, optimize_3d=True)
                if mol is None:
                    return {}

            distances = {}
            conf = mol.GetConformer()

            # Calculate distances between all feature pairs
            feature_types = [ft for ft in features.keys() if ft not in ["fingerprint", "gobbi_fingerprint"]]

            for i, type1 in enumerate(feature_types):
                for type2 in feature_types[i:]:
                    key = f"{type1}-{type2}"
                    distances[key] = []

                    for atoms1 in features[type1]:
                        for atoms2 in features[type2]:
                            # Use centroids for multi-atom features
                            pos1 = sum(conf.GetAtomPosition(a) for a in atoms1) / len(atoms1)
                            pos2 = sum(conf.GetAtomPosition(a) for a in atoms2) / len(atoms2)
                            dist = pos1.Distance(pos2)
                            distances[key].append(float(dist))

            return distances

        except Exception as e:
            self.logger.error(f"Error calculating feature distances: {str(e)}")
            return {}

    def visualize_features(self, mol: Chem.Mol, features: Dict[str, List[Tuple[int, ...]]], size: Tuple[int, int] = (400, 400)) -> Optional[Draw.Image]:
        """
        Generate a 2D depiction of the molecule with highlighted pharmacophore features.

        Args:
            mol: RDKit molecule
            features: Dictionary of pharmacophore features
            size: Image size as (width, height)

        Returns:
            PIL Image object or None if visualization fails
        """
        try:
            if mol is None or not features:
                return None

            # Make a copy for drawing
            mol = Chem.Mol(mol)

            # Generate 2D coordinates if needed
            if not mol.GetNumConformers():
                rdDepictor.Compute2DCoords(mol)

            # Define highlight colors for different feature types
            colors = {
                "hbd": (1, 0, 0),  # Red
                "hba": (0, 0, 1),  # Blue
                "pi": (0, 1, 0),  # Green
                "ni": (1, 0.5, 0),  # Orange
                "ar": (0.5, 0, 0.5),  # Purple
                "hp": (0.5, 0.5, 0.5),  # Gray
                "metal": (1, 1, 0),  # Yellow
                "xbond": (0, 1, 1),  # Cyan
            }

            # Collect atoms to highlight for each feature type
            highlights = {}
            for feat_type, matches in features.items():
                if feat_type not in ["fingerprint", "gobbi_fingerprint"]:
                    highlights[feat_type] = set()
                    for match in matches:
                        highlights[feat_type].update(match)

            # Generate depiction
            img = Draw.MolToImage(
                mol, size=size, highlightAtoms=[atom for atoms in highlights.values() for atom in atoms], highlightColor=colors.get(next(iter(highlights.keys())), (0.7, 0.7, 0.7))
            )

            return img

        except Exception as e:
            self.logger.error(f"Error visualizing features: {str(e)}")
            return None
