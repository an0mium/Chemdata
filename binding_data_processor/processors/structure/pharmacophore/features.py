"""Pharmacophore feature definitions and patterns."""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


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
    shape_weight: float = 1.0  # Weight for shape-based alignment
    exclusion_radius: Optional[float] = None  # Radius for exclusion volumes


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
    # Stereocenter
    "stereocenter": [
        "[!H0;$(C([*])([*])([*])[*])]",  # Tetrahedral carbon
        "[!H0;$(N([*])([*])([*]))]",  # Tetrahedral nitrogen
        "[!H0;$(P([*])([*])([*]))]",  # Tetrahedral phosphorus
        "[!H0;$(S([*])([*])([*]))]",  # Tetrahedral sulfur
    ],
    # Conjugated systems
    "conjugated": [
        "[$(C=C-C=C),$(C=C-C=N),$(C=C-N=N)]",  # Conjugated double bonds
        "[$(c:c:c:c:c)]",  # Conjugated aromatic systems
        "[$(C=C-C#N),$(C=C-C=O)]",  # Conjugated with electron-withdrawing groups
    ],
    # Sulfonamide
    "sulfonamide": [
        "[$([#16X4]([NX3H])([!O])[OX2H0])]",  # Sulfonamide group
    ],
    # Phosphate
    "phosphate": [
        "[$([PX4](=[OX1])([OX2H0])([OX2H0])[OX2H0])]",  # Phosphate group
    ],
    # Guanidine
    "guanidine": [
        "[$([NX3H2][CX3](=[NX2H0])[NX3H2])]",  # Guanidine group
    ],
    # Heterocycles by ring size
    "heterocycle_3": [
        "[$(C1[O,N,S]C1)]",  # 3-membered heterocycles
        "[$(O1CC1)]",  # Oxirane
        "[$(N1CC1)]",  # Aziridine
        "[$(S1CC1)]",  # Thiirane
    ],
    "heterocycle_4": [
        "[$(C1[O,N,S]CC1)]",  # 4-membered heterocycles
        "[$(O1CCC1)]",  # Oxetane
        "[$(N1CCC1)]",  # Azetidine
        "[$(S1CCC1)]",  # Thietane
    ],
    "heterocycle_5": [
        "[$(C1[O,N,S]CCC1)]",  # 5-membered heterocycles
        "[$(O1CCCC1)]",  # Tetrahydrofuran
        "[$(N1CCCC1)]",  # Pyrrolidine
        "[$(S1CCCC1)]",  # Tetrahydrothiophene
        "c1[nH]ccc1",  # Pyrrole
        "c1nccc1",  # Pyrazole
        "c1occn1",  # Oxazole
        "c1sccc1",  # Thiophene
    ],
    "heterocycle_6": [
        "[$(C1[O,N,S]CCCC1)]",  # 6-membered heterocycles
        "[$(O1CCCCC1)]",  # Tetrahydropyran
        "[$(N1CCCCC1)]",  # Piperidine
        "[$(S1CCCCC1)]",  # Tetrahydrothiopyran
        "c1ncccc1",  # Pyridine
        "c1ncncc1",  # Pyrimidine
        "c1nnccc1",  # Pyridazine
        "c1cnccn1",  # Pyrazine
    ],
    "heterocycle_7": [
        "[$(C1[O,N,S]CCCCC1)]",  # 7-membered heterocycles
        "[$(O1CCCCCC1)]",  # Oxepane
        "[$(N1CCCCCC1)]",  # Azepane
        "[$(S1CCCCCC1)]",  # Thiepane
    ],
    "heterocycle_8": [
        "[$(C1[O,N,S]CCCCCC1)]",  # 8-membered heterocycles
        "[$(O1CCCCCCC1)]",  # Oxocane
        "[$(N1CCCCCCC1)]",  # Azocane
        "[$(S1CCCCCCC1)]",  # Thiocane
    ],
    # Polycyclic systems
    "bicyclic": [
        "[$(C12CCCCC1CCCC2)]",  # Decalin-like
        "[$(C12CCCCC1CCC2)]",  # Indane-like
        "[$(C12CCCCC1CC2)]",  # Bicyclo[4.2.0]
        "[$(C12CCC1CC2)]",  # Bicyclo[3.1.0]
    ],
    "tricyclic": [
        "[$(C12C3CCCCC1CCCC2CCC3)]",  # Anthracene-like
        "[$(C12C3CCCCC1CCCC2CC3)]",  # Phenalene-like
        "[$(C12C3CCCCC1CCC2CC3)]",  # Fluorene-like
    ],
    "spiro": [
        "[$(C12(CCCC1)CCCC2)]",  # Spiro[4.4]
        "[$(C12(CCCCC1)CCCC2)]",  # Spiro[5.4]
        "[$(C12(CCCCC1)CCCCC2)]",  # Spiro[5.5]
    ],
    # Advanced stereochemistry
    "chiral_carbon": [
        "[C@H]([*])([*])[*]",  # R configuration
        "[C@@H]([*])([*])[*]",  # S configuration
    ],
    "chiral_nitrogen": [
        "[N@H]([*])([*])[*]",  # R configuration
        "[N@@H]([*])([*])[*]",  # S configuration
    ],
    "chiral_sulfur": [
        "[S@]([*])([*])[*]",  # R configuration
        "[S@@]([*])([*])[*]",  # S configuration
    ],
    # Drug-like functional groups
    "amide": [
        "[NX3H2,NX3H1,NX3H0][CX3](=[OX1])[#6]",  # General amide
        "[NX3H2][CX3](=[OX1])[#6]",  # Primary amide
        "[NX3H1][CX3](=[OX1])[#6]",  # Secondary amide
        "[NX3H0][CX3](=[OX1])[#6]",  # Tertiary amide
    ],
    "ester": [
        "[#6][CX3](=[OX1])[OX2H0][#6]",  # General ester
        "[#6]C(=O)OC",  # Simple ester
        "[#6]C(=O)OC(=O)[#6]",  # Anhydride
    ],
    "carbamate": [
        "[NX3][CX3](=[OX1])[OX2H0][#6]",  # General carbamate
        "[NH2]C(=O)O[#6]",  # Primary carbamate
        "[NH1][CX3](=[OX1])[OX2H0][#6]",  # Secondary carbamate
    ],
    "urea": [
        "[NX3][CX3](=[OX1])[NX3]",  # General urea
        "[NH2]C(=O)[NH2]",  # Urea
        "[NH2]C(=O)[NH1][#6]",  # N-substituted urea
    ],
    "sulfonamide": [
        "[SX4](=[OX1])(=[OX1])([NX3])[#6]",  # General sulfonamide
        "[SX4](=[OX1])(=[OX1])([NH2])[#6]",  # Primary sulfonamide
        "[SX4](=[OX1])(=[OX1])([NH1][#6])[#6]",  # Secondary sulfonamide
    ],
    # Extended aromatic systems
    "phenol": [
        "[OH]c1ccccc1",  # Simple phenol
        "[OH]c1cc([#6,#7,#8,#9,#16,#17])ccc1",  # Substituted phenol
    ],
    "aniline": [
        "[NH2]c1ccccc1",  # Simple aniline
        "[NH2]c1cc([#6,#7,#8,#9,#16,#17])ccc1",  # Substituted aniline
    ],
    "indole": [
        "c1ccc2c(c1)[nH]cc2",  # Indole
        "c1ccc2c(c1)n([#6,#1])cc2",  # N-substituted indole
    ],
    # Extended heterocycles
    "benzofused": [
        "c1ccc2[o,n,s]cccc2c1",  # Benzofuran, indole, benzothiophene
        "c1ccc2[o,n,s]ccc2c1",  # Isobenzofuran, isoindole, isobenzothiophene
        "c1ccc2[nH]cnc2c1",  # Benzimidazole
        "c1ccc2ncnc2c1",  # Quinazoline
        "c1ccc2nccn2c1",  # Quinoxaline
    ],
    "bridged_heterocycles": [
        "[$(C12[O,N,S]CCC1CC2)]",  # 2-oxabicyclo[2.2.1]
        "[$(C12[O,N,S]CCCC1CC2)]",  # 2-oxabicyclo[2.2.2]
        "[$(C12[O,N,S]CCC1CCC2)]",  # 2-oxabicyclo[3.2.1]
    ],
    "macrocycle": [
        "[$(C1CCCCCCCCCC1)]",  # 11-membered ring
        "[$(C1CCCCCCCCCCC1)]",  # 12-membered ring
        "[$(C1CCCCCCCCCCCC1)]",  # 13-membered ring
        "[$(C1CCCCCCCCCCCCC1)]",  # 14-membered ring
    ],
    "multiple_heteroatoms": [
        "[$(C1[O,N,S]CC[O,N,S]CC1)]",  # 1,4-heterocycles
        "[$(C1[O,N,S]CCC[O,N,S]CC1)]",  # 1,5-heterocycles
        "[$(C1[O,N,S]C[O,N,S]C[O,N,S]C1)]",  # Triple heteroatom
    ],
    "fused_heterocycles": [
        "c1ccc2c(c1)nccn2",  # Quinoxaline
        "c1ccc2c(c1)ncnc2",  # Quinazoline
        "c1ccc2c(c1)ccnc2",  # Quinoline
        "c1ccc2c(c1)cncc2",  # Isoquinoline
    ],
    "spiro_heterocycles": [
        "[$(C12(CCOC1)CCNC2)]",  # Spiro[oxolane-2,2'-pyrrolidine]
        "[$(C12(CCSC1)CCNC2)]",  # Spiro[thiolane-2,2'-pyrrolidine]
        "[$(C12(CCNC1)CCNC2)]",  # Spiro[pyrrolidine-2,2'-pyrrolidine]
    ],
    "bridged_heterocycles": [
        "[$(C12CCNC1CC2)]",  # 2-azabicyclo[2.2.1]
        "[$(C12CCOC1CC2)]",  # 2-oxabicyclo[2.2.1]
        "[$(C12CCSC1CC2)]",  # 2-thiabicyclo[2.2.1]
    ],
    "acyl_groups": [
        "[CX3](=[OX1])[#6]",  # General acyl
        "[CX3](=[OX1])C",  # Acetyl
        "[CX3](=[OX1])CC",  # Propionyl
        "[CX3](=[OX1])c1ccccc1",  # Benzoyl
    ],
    "alkyl_halides": [
        "C[F,Cl,Br,I]",  # Primary
        "CC[F,Cl,Br,I]",  # Secondary
        "C(C)(C)[F,Cl,Br,I]",  # Tertiary
        "C(F)(F)F",  # Trifluoromethyl
    ],
    "nitriles": [
        "[NX1]#[CX2]",  # General nitrile
        "CC#N",  # Acetonitrile
        "c1ccccc1C#N",  # Benzonitrile
    ],
    "isonitriles": [
        "[CX1-]#[NX2+]",  # General isonitrile
        "C[N+]#[C-]",  # Methyl isocyanide
    ],
    "azides": [
        "[NX1]=[NX2+]=[NX1-]",  # Azide
        "C[NX1]=[NX2+]=[NX1-]",  # Alkyl azide
        "c1ccccc1[NX1]=[NX2+]=[NX1-]",  # Aryl azide
    ],
    "diazo": [
        "[NX1]=[NX1]",  # Diazo group
        "C[NX1]=[NX1]",  # Alkyl diazo
        "c1ccccc1[NX1]=[NX1]",  # Aryl diazo
    ],
    # Stereochemistry patterns
    "ez_alkene": [
        "[CH]=[CH]/[CH3]",  # E-alkene
        "[CH]=[CH]\[CH3]",  # Z-alkene
        "[CH]=[CH]/c1ccccc1",  # E-styrene
        "[CH]=[CH]\c1ccccc1",  # Z-styrene
    ],
    "axial_chirality": [
        "[aX2H0](-[aX2H0])-[aX2H0](-[aX2H0])",  # Biaryl
        "[C@H]1CC[C@H](CC1)c1ccccc1",  # Allene-like
    ],
    "planar_chirality": [
        "[Fe](C1C=CC=C1)(C2C=CC=C2)",  # Ferrocene
        "C1=CC=C(C=C1)[C@]2(C=CC=C2)",  # Paracyclophane
    ],
    "crown_ethers": [
        "[$(C1COCCOCCOCCO1)]",  # 12-crown-4
        "[$(C1COCCOCCOCCOCCO1)]",  # 15-crown-5
        "[$(C1COCCOCCOCCOCCOCCO1)]",  # 18-crown-6
    ],
    "peptide_bonds": [
        "[NX3H1][CX3](=[OX1])[CX4H1][NX3H1]",  # Dipeptide
        "[NX3H1][CX3](=[OX1])[CX4H1][NX3H1][CX3](=[OX1])[CX4H1][NX3H1]",  # Tripeptide
    ],
    "sugar_rings": [
        "[$(C1[C@H]([OH])[C@@H]([OH])[C@H]([OH])[C@H]([OH])O1)]",  # Pyranose
        "[$(C1[C@H]([OH])[C@@H]([OH])[C@H]([OH])CO1)]",  # Furanose
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
        "heterocycle_5": (0.8, 0.2, 0.8),  # Purple
        "heterocycle_6": (0.2, 0.8, 0.8),  # Cyan
        "bicyclic": (0.8, 0.4, 0.4),  # Pink
        "tricyclic": (0.4, 0.8, 0.4),  # Light green
        "spiro": (0.4, 0.4, 0.8),  # Light blue
        "chiral_carbon": (1.0, 0.8, 0.0),  # Gold
        "benzofused": (0.6, 0.3, 0.6),  # Dark purple
        "bridged": (0.3, 0.6, 0.6),  # Dark cyan
        "ez_alkene": (0.8, 0.3, 0.3),  # Red
        "axial_chirality": (0.3, 0.8, 0.3),  # Green
        "planar_chirality": (0.3, 0.3, 0.8),  # Blue
    },
    "colorblind": {
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
        "heterocycle_5": (0.7, 0.5, 0.2),  # Brown
        "heterocycle_6": (0.2, 0.7, 0.5),  # Teal
        "bicyclic": (0.5, 0.2, 0.7),  # Purple
        "tricyclic": (0.2, 0.5, 0.7),  # Blue
        "spiro": (0.7, 0.2, 0.5),  # Magenta
        "chiral_carbon": (0.5, 0.7, 0.2),  # Olive
        "benzofused": (0.4, 0.4, 0.7),  # Blue-gray
        "bridged": (0.7, 0.4, 0.4),  # Red-gray
        "ez_alkene": (0.7, 0.5, 0.3),  # Brown
        "axial_chirality": (0.3, 0.7, 0.5),  # Teal
        "planar_chirality": (0.5, 0.3, 0.7),  # Purple
    },
}

# Update color schemes with new features
COLOR_SCHEMES["default"].update(
    {
        "macrocycle": (0.7, 0.3, 0.7),  # Purple
        "multiple_heteroatoms": (0.3, 0.7, 0.7),  # Cyan
        "fused_heterocycles": (0.7, 0.7, 0.3),  # Yellow
        "spiro_heterocycles": (0.3, 0.3, 0.7),  # Blue
        "bridged_heterocycles": (0.7, 0.3, 0.3),  # Red
        "acyl_groups": (0.5, 0.5, 0.7),  # Blue-gray
        "alkyl_halides": (0.7, 0.5, 0.5),  # Red-gray
        "nitriles": (0.5, 0.7, 0.5),  # Green-gray
        "isonitriles": (0.7, 0.7, 0.5),  # Yellow-gray
        "azides": (0.5, 0.5, 0.3),  # Dark gray
        "diazo": (0.3, 0.5, 0.5),  # Dark cyan
        "crown_ethers": (0.8, 0.6, 0.8),  # Light purple
        "peptide_bonds": (0.6, 0.8, 0.6),  # Light green
        "sugar_rings": (0.8, 0.8, 0.6),  # Light yellow
    }
)

COLOR_SCHEMES["colorblind"].update(
    {
        "macrocycle": (0.6, 0.4, 0.7),  # Purple-blue
        "multiple_heteroatoms": (0.4, 0.7, 0.6),  # Teal-green
        "fused_heterocycles": (0.7, 0.6, 0.4),  # Orange-brown
        "spiro_heterocycles": (0.4, 0.4, 0.7),  # Blue-gray
        "bridged_heterocycles": (0.7, 0.4, 0.4),  # Red-gray
        "acyl_groups": (0.5, 0.6, 0.7),  # Light blue
        "alkyl_halides": (0.7, 0.5, 0.6),  # Pink-purple
        "nitriles": (0.6, 0.7, 0.5),  # Yellow-green
        "isonitriles": (0.7, 0.7, 0.6),  # Light yellow
        "azides": (0.5, 0.5, 0.4),  # Gray-brown
        "diazo": (0.4, 0.5, 0.5),  # Gray-blue
        "crown_ethers": (0.7, 0.5, 0.7),  # Light purple
        "peptide_bonds": (0.5, 0.7, 0.5),  # Light green
        "sugar_rings": (0.7, 0.7, 0.5),  # Light yellow
    }
)

# Feature factory definitions
FEATURE_FACTORY_DEFS = [
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
    # Stereocenter
    """DefineFeature SC [!H0;$(C([*])([*])([*])[*])]
        Family Stereocenter
        Weights 1.0
    EndFeature""",
    # Conjugated system
    """DefineFeature CJ [$(C=C-C=C),$(C=C-C=N),$(C=C-N=N)]
        Family Conjugated
        Weights 1.0,1.0,1.0,1.0
    EndFeature""",
    # Sulfonamide
    """DefineFeature SU [$([#16X4]([NX3H])([!O])[OX2H0])]
        Family Sulfonamide
        Weights 1.0
    EndFeature""",
    # Phosphate
    """DefineFeature PH [$([PX4](=[OX1])([OX2H0])([OX2H0])[OX2H0])]
        Family Phosphate
        Weights 1.0
    EndFeature""",
    # Guanidine
    """DefineFeature GU [$([NX3H2][CX3](=[NX2H0])[NX3H2])]
        Family Guanidine
        Weights 1.0
    EndFeature""",
    # Additional heterocycle definitions
    """DefineFeature HC5 [$(C1[O,N,S]CCC1)]
        Family Heterocycle5
        Weights 1.0,1.0,1.0,1.0,1.0
    EndFeature""",
    """DefineFeature HC6 [$(C1[O,N,S]CCCC1)]
        Family Heterocycle6
        Weights 1.0,1.0,1.0,1.0,1.0,1.0
    EndFeature""",
    # Polycyclic definitions
    """DefineFeature BC [$(C12CCCCC1CCCC2)]
        Family Bicyclic
        Weights 1.0
    EndFeature""",
    """DefineFeature TC [$(C12C3CCCCC1CCCC2CCC3)]
        Family Tricyclic
        Weights 1.0
    EndFeature""",
    # Chiral center definitions
    """DefineFeature CC [C@H,C@@H]([*])([*])[*]
        Family ChiralCarbon
        Weights 1.0
    EndFeature""",
    # Add stereochemistry definitions
    """DefineFeature EZ [CH]=[CH]/[CH3]
        Family EZAlkene
        Weights 1.0
    EndFeature""",
    """DefineFeature AX [aX2H0](-[aX2H0])-[aX2H0](-[aX2H0])
        Family AxialChirality
        Weights 1.0
    EndFeature""",
    """DefineFeature PC [Fe](C1C=CC=C1)(C2C=CC=C2)
        Family PlanarChirality
        Weights 1.0
    EndFeature""",
]

# Add new feature factory definitions
FEATURE_FACTORY_DEFS.extend(
    [
        """DefineFeature MC [$(C1CCCCCCCCCC1)]
        Family Macrocycle
        Weights 1.0
    EndFeature""",
        """DefineFeature MH [$(C1[O,N,S]CC[O,N,S]CC1)]
        Family MultipleHeteroatoms
        Weights 1.0
    EndFeature""",
        """DefineFeature FH c1ccc2c(c1)nccn2
        Family FusedHeterocycles
        Weights 1.0
    EndFeature""",
        """DefineFeature SH [$(C12(CCOC1)CCNC2)]
        Family SpiroHeterocycles
        Weights 1.0
    EndFeature""",
        """DefineFeature BH [$(C12CCNC1CC2)]
        Family BridgedHeterocycles
        Weights 1.0
    EndFeature""",
        """DefineFeature AG [CX3](=[OX1])[#6]
        Family AcylGroups
        Weights 1.0
    EndFeature""",
        """DefineFeature AH C[F,Cl,Br,I]
        Family AlkylHalides
        Weights 1.0
    EndFeature""",
        """DefineFeature NI [NX1]#[CX2]
        Family Nitriles
        Weights 1.0
    EndFeature""",
        """DefineFeature II [CX1-]#[NX2+]
        Family Isonitriles
        Weights 1.0
    EndFeature""",
        """DefineFeature AZ [NX1]=[NX2+]=[NX1-]
        Family Azides
        Weights 1.0
    EndFeature""",
        """DefineFeature DZ [NX1]=[NX1]
        Family Diazo
        Weights 1.0
    EndFeature""",
        """DefineFeature CE [$(C1COCCOCCOCCOCCO1)]
        Family CrownEthers
        Weights 1.0
    EndFeature""",
        """DefineFeature PB [NX3H1][CX3](=[OX1])[CX4H1][NX3H1]
        Family PeptideBonds
        Weights 1.0
    EndFeature""",
        """DefineFeature SR [$(C1[C@H]([OH])[C@@H]([OH])[C@H]([OH])[C@H]([OH])O1)]
        Family SugarRings
        Weights 1.0
    EndFeature""",
    ]
)


# Signature factory patterns for 2D pharmacophore fingerprints
SIGNATURE_FACTORY_PATTERNS = {
    "HBD": ["[N,O,S;H1,H2]-[!$(*=[O,N,P,S])]"],
    "HBA": ["[$([O,S;H0;v2]),$([O,S;-])]"],
    "PI": ["[N;H2,H3;+1]"],
    "NI": ["[C,S](=[O,S])[O-]"],
    "AR": ["a1aaaaa1"],
    "HP": ["[C;!$(C=[O,N,S]);!$(C[F,Cl,Br,I,N,O,P,S])]"],
    "MB": ["[O,N,S;!H0;+0]"],
    "XB": ["[F,Cl,Br,I]"],
    "HC5": ["[$(C1[O,N,S]CCC1)]"],
    "HC6": ["[$(C1[O,N,S]CCCC1)]"],
    "BC": ["[$(C12CCCCC1CCCC2)]"],
    "TC": ["[$(C12C3CCCCC1CCCC2CCC3)]"],
    "CC": ["[C@H,C@@H]([*])([*])[*]"],
    "EZ": ["[CH]=[CH]/[CH3]"],  # E/Z alkenes
    "AX": ["[aX2H0](-[aX2H0])-[aX2H0](-[aX2H0])"],  # Axial chirality
    "PC": ["[Fe](C1C=CC=C1)(C2C=CC=C2)"],  # Planar chirality
}

# Add new signature factory patterns
SIGNATURE_FACTORY_PATTERNS.update(
    {
        "MC": ["[$(C1CCCCCCCCCC1)]"],  # Macrocycle
        "MH": ["[$(C1[O,N,S]CC[O,N,S]CC1)]"],  # Multiple heteroatoms
        "FH": ["c1ccc2c(c1)nccn2"],  # Fused heterocycles
        "SH": ["[$(C12(CCOC1)CCNC2)]"],  # Spiro heterocycles
        "BH": ["[$(C12CCNC1CC2)]"],  # Bridged heterocycles
        "AG": ["[CX3](=[OX1])[#6]"],  # Acyl groups
        "AH": ["C[F,Cl,Br,I]"],  # Alkyl halides
        "NI": ["[NX1]#[CX2]"],  # Nitriles
        "II": ["[CX1-]#[NX2+]"],  # Isonitriles
        "AZ": ["[NX1]=[NX2+]=[NX1-]"],  # Azides
        "DZ": ["[NX1]=[NX1]"],  # Diazo
        "CE": ["[$(C1COCCOCCOCCOCCO1)]"],  # Crown ethers
        "PB": ["[NX3H1][CX3](=[OX1])[CX4H1][NX3H1]"],  # Peptide bonds
        "SR": ["[$(C1[C@H]([OH])[C@@H]([OH])[C@H]([OH])[C@H]([OH])O1)]"],  # Sugar rings
    }
)
