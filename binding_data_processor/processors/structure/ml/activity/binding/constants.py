"""Constants and property scales for binding site prediction and analysis."""

# Residue property scales
HYDROPHOBICITY = {
    "ALA": 1.8,
    "ARG": -4.5,
    "ASN": -3.5,
    "ASP": -3.5,
    "CYS": 2.5,
    "GLN": -3.5,
    "GLU": -3.5,
    "GLY": -0.4,
    "HIS": -3.2,
    "ILE": 4.5,
    "LEU": 3.8,
    "LYS": -3.9,
    "MET": 1.9,
    "PHE": 2.8,
    "PRO": -1.6,
    "SER": -0.8,
    "THR": -0.7,
    "TRP": -0.9,
    "TYR": -1.3,
    "VAL": 4.2,
}

# Residue volume scale (Å³)
VOLUME = {
    "ALA": 88.6,
    "ARG": 173.4,
    "ASN": 114.1,
    "ASP": 111.1,
    "CYS": 108.5,
    "GLN": 143.8,
    "GLU": 138.4,
    "GLY": 60.1,
    "HIS": 153.2,
    "ILE": 166.7,
    "LEU": 166.7,
    "LYS": 168.6,
    "MET": 162.9,
    "PHE": 189.9,
    "PRO": 112.7,
    "SER": 89.0,
    "THR": 116.1,
    "TRP": 227.8,
    "TYR": 193.6,
    "VAL": 140.0,
}

# Residue charge scale (pKa dependent)
CHARGE = {
    "ARG": 1,  # Positive
    "LYS": 1,  # Positive
    "ASP": -1,  # Negative
    "GLU": -1,  # Negative
    "HIS": 0.1,  # Slightly positive (pKa dependent)
}

# Binding site detection parameters
BINDING_SITE_PARAMS = {
    "min_volume": 100.0,  # Minimum pocket volume in Å³
    "max_exposure": 0.5,  # Maximum solvent exposure
    "min_depth": 4.0,  # Minimum pocket depth in Å
    "conservation_cutoff": 0.7,  # Minimum conservation score
    "interaction_distance": 4.5,  # Maximum distance for interactions
    "probe_radius": 1.4,  # Solvent probe radius in Å
    "grid_spacing": 0.5,  # Grid spacing for energy calculations
    "energy_cutoff": -2.0,  # Energy cutoff for favorable regions
    "min_interface_contacts": 3,  # Minimum contacts for interface
    "interface_cutoff": 10.0,  # Cutoff for interface contacts
    "min_community_size": 4,  # Minimum residues for community
    "edge_weight_factor": 1.0,  # Factor for distance-based edge weights
    "pi_stacking_cutoff": 7.0,  # Å, π-stacking interaction cutoff
    "halogen_bond_cutoff": 4.0,  # Å, halogen bond cutoff
    "metal_coord_cutoff": 3.0,  # Å, metal coordination cutoff
    "interface_min_size": 4,  # Minimum interface size
    "interface_core_cutoff": 0.7,  # Core residue burial threshold
    "network_weight_scale": 2.0,  # Scale factor for distance-based weights
    "community_resolution": 1.0,  # Resolution parameter for community detection
    "water_bridge_cutoff": 3.5,  # Å, water-mediated interaction cutoff
    "hotspot_threshold": 0.7,  # Fraction of total interface energy for hotspots
    "allosteric_distance": 8.0,  # Å, distance to consider for allosteric effects
    "path_weight_factor": 2.0,  # Factor for weighting paths by interaction strength
    "community_edge_weight": 1.5,  # Weight factor for community detection
    "betweenness_cutoff": 0.1,  # Cutoff for significant betweenness centrality
}

# Force field parameters
FORCEFIELD_PARAMS = {
    "MMFF94": {"maxIters": 200, "nonBondedThresh": 100.0, "confId": -1, "ignoreInterfragInteractions": False},
    "UFF": {"maxIters": 200, "vdwThresh": 10.0, "confId": -1, "ignoreInterfragInteractions": False},
}

# Conformer generation parameters
CONFORMER_PARAMS = {
    "num_confs": 50,  # Number of conformers to generate
    "num_threads": 0,  # Use all available CPUs
    "random_seed": 42,
    "pruneRmsThresh": 0.5,  # Remove similar conformers
    "enforceChirality": True,
    "useExpTorsionAnglePrefs": True,
    "useBasicKnowledge": True,
    "ETversion": 2,  # Use newer version of experimental torsion preferences
    "useSmallRingTorsions": True,
    "useMacrocycleTorsions": True,
    "forceTransAmides": True,
}

# Feature type weights for scoring
FEATURE_WEIGHTS = {
    "donors": 1.0,
    "acceptors": 1.0,
    "aromatic": 1.2,
    "hydrophobic": 0.8,
    "positive": 1.5,
    "negative": 1.5,
}

# Energy terms for binding energy estimation (kcal/mol)
ENERGY_TERMS = {
    "donors": -2.0,  # H-bond donor
    "acceptors": -2.0,  # H-bond acceptor
    "aromatic": -2.5,  # π-π stacking
    "hydrophobic": -0.8,  # Hydrophobic contact
    "positive": -3.0,  # Salt bridge
    "negative": -3.0,  # Salt bridge
}

# Geometric parameters
GEOMETRY_PARAMS = {
    "hbond_distance": 3.5,  # Maximum H-bond distance in Å
    "hbond_angle": 30.0,  # Maximum H-bond angle deviation in degrees
    "salt_bridge_dist": 4.0,  # Maximum salt bridge distance in Å
    "disulfide_dist": 2.2,  # Maximum disulfide bond distance in Å
    "feature_match_dist": 4.0,  # Maximum pharmacophore feature matching distance in Å
    "pocket_overlap_dist": 5.0,  # Distance for considering pockets as overlapping
}

# Scoring weights for different components
SCORING_WEIGHTS = {
    "volume": 0.25,
    "hydrophobicity": 0.20,
    "conservation": 0.15,
    "shape": 0.15,
    "depth": 0.10,
    "accessibility": 0.10,
    "alphafold": 0.05,
}

# Energy calculation parameters
ENERGY_PARAMS = {
    "vdw_epsilon": 0.1,  # Energy well depth for van der Waals
    "vdw_sigma": 3.5,  # Distance at zero energy
    "dielectric": 80.0,  # Dielectric constant for electrostatics
    "cutoff": 8.0,  # Distance cutoff for interactions
}

# Feature colors (RGB)
FEATURE_COLORS = {
    "hbd": (1, 0, 0),  # Red
    "hba": (0, 0, 1),  # Blue
    "pi": (0, 1, 0),  # Green
    "ni": (1, 0.5, 0),  # Orange
    "ar": (0.5, 0, 0.5),  # Purple
    "hp": (0.5, 0.5, 0.5),  # Gray
    "metal": (1, 1, 0),  # Yellow
    "xbond": (0, 1, 1),  # Cyan
}

# Colorblind-friendly feature colors
COLORBLIND_COLORS = {
    "hbd": (0.9, 0.6, 0),  # Orange
    "hba": (0.35, 0.7, 0.9),  # Light blue
    "pi": (0, 0.6, 0.5),  # Teal
    "ni": (0.95, 0.9, 0.25),  # Yellow
    "ar": (0.8, 0.4, 0),  # Brown
    "hp": (0.8, 0.8, 0.8),  # Gray
    "metal": (0.35, 0.35, 0.35),  # Dark gray
    "xbond": (0.9, 0.6, 0.6),  # Pink
}

# Network analysis parameters
NETWORK_PARAMS = {
    "min_community_size": 4,
    "edge_weight_factor": 1.0,
    "community_resolution": 1.0,
    "path_weight_factor": 2.0,
    "betweenness_cutoff": 0.1,
}

# Interaction parameters
INTERACTION_PARAMS = {
    "hydrophobic_cutoff": 5.0,  # Å
    "hbond_distance": 3.5,  # Å
    "hbond_angle": 30.0,  # degrees
    "ionic_cutoff": 4.0,  # Å
    "aromatic_cutoff": 6.0,  # Å
    "cation_pi_cutoff": 6.0,  # Å
    "disulfide_cutoff": 2.2,  # Å
    "metal_coord_cutoff": 3.0,  # Å
}

# Atom radii (Å)
ATOM_RADII = {
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "S": 1.8,
    "P": 1.8,
    "H": 1.2,
    "F": 1.47,
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
}

# Ring atoms for aromatic residues
RING_ATOMS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}

# Metal coordinating atoms
METAL_COORD_ATOMS = {
    "HIS": ["ND1", "NE2"],
    "CYS": ["SG"],
    "MET": ["SD"],
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}

# H-bond donor-acceptor pairs
HBOND_PAIRS = {
    "donors": {
        "N": ["H"],  # Backbone
        "NE": ["HE"],  # Arg
        "NH1": ["HH11", "HH12"],  # Arg
        "NH2": ["HH21", "HH22"],  # Arg
        "NZ": ["HZ1", "HZ2", "HZ3"],  # Lys
        "ND1": ["HD1"],  # His
        "NE2": ["HE2"],  # His/Gln
        "ND2": ["HD21", "HD22"],  # Asn
        "OG": ["HG"],  # Ser
        "OG1": ["HG1"],  # Thr
        "OH": ["HH"],  # Tyr
        "NE1": ["HE1"],  # Trp
    },
    "acceptors": {"O", "OD1", "OD2", "OE1", "OE2", "ND1", "NE2"},
}
