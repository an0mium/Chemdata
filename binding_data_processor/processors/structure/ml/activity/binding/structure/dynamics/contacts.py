"""Contact network and residue interaction analysis."""

import logging
from typing import Dict, List, Optional, Set, Any, Tuple
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from scipy.spatial.distance import cdist
import networkx as nx

# Add new imports
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio.PDB.vectors import calc_dihedral, calc_angle


class ContactAnalyzer:
    """Analyzes protein contact networks and residue interactions."""

    def __init__(self):
        """Initialize contact analyzer."""
        self.logger = logging.getLogger(__name__)
        self.settings = {
            "contact_cutoff": 8.0,  # Å, distance cutoff for contacts
            "interface_cutoff": 10.0,  # Å, cutoff for interface contacts
            "min_interface_contacts": 3,  # Minimum contacts for interface
            "hydrophobic_cutoff": 5.0,  # Å, cutoff for hydrophobic contacts
            "hbond_distance": 3.5,  # Å, H-bond distance cutoff
            "hbond_angle": 30.0,  # degrees, H-bond angle cutoff
            "ionic_cutoff": 4.0,  # Å, ionic interaction cutoff
            "aromatic_cutoff": 6.0,  # Å, π-π interaction cutoff
            "cation_pi_cutoff": 6.0,  # Å, cation-π interaction cutoff
            "disulfide_cutoff": 2.2,  # Å, disulfide bond cutoff
            "min_community_size": 4,  # Minimum residues for a community
            "edge_weight_factor": 1.0,  # Factor for distance-based edge weights
            "pi_stacking_cutoff": 7.0,  # Å, π-stacking interaction cutoff
            "halogen_bond_cutoff": 4.0,  # Å, halogen bond cutoff
            "metal_coord_cutoff": 3.0,  # Å, metal coordination cutoff
            "interface_min_size": 4,  # Minimum interface size
            "interface_core_cutoff": 0.7,  # Core residue burial threshold
            "network_weight_scale": 2.0,  # Scale factor for distance-based weights
            "community_resolution": 1.0,  # Resolution parameter for community detection
        }

        # Residue properties
        self.hydrophobic = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
        self.charged_pos = {"ARG", "LYS", "HIS"}
        self.charged_neg = {"ASP", "GLU"}
        self.aromatic = {"PHE", "TYR", "TRP", "HIS"}
        self.hbond_donors = {"ARG", "LYS", "ASN", "GLN", "HIS", "SER", "THR", "TYR", "TRP"}
        self.hbond_acceptors = {"ASP", "GLU", "ASN", "GLN", "HIS", "SER", "THR", "TYR"}
        self.pi_stacking = {"PHE", "TYR", "TRP", "HIS"}  # π-stacking capable
        self.halogen_donors = {"PHE", "TYR", "TRP", "LEU", "ILE", "VAL"}  # Can donate halogen bonds
        self.halogen_acceptors = {"ASP", "GLU", "ASN", "GLN", "SER", "THR"}  # Can accept halogen bonds
        self.metal_coordinating = {"HIS", "CYS", "MET", "ASP", "GLU"}  # Can coordinate metals

        # Add new settings
        self.settings.update(
            {
                "water_bridge_cutoff": 3.5,  # Å, water-mediated interaction cutoff
                "hotspot_threshold": 0.7,  # Fraction of total interface energy for hotspots
                "allosteric_distance": 8.0,  # Å, distance to consider for allosteric effects
                "path_weight_factor": 2.0,  # Factor for weighting paths by interaction strength
                "community_edge_weight": 1.5,  # Weight factor for community detection
                "betweenness_cutoff": 0.1,  # Cutoff for significant betweenness centrality
            }
        )

    def analyze_contacts(
        self,
        structure: Structure,
        include_interactions: bool = True,
        include_interfaces: bool = True,
        include_network: bool = True,
    ) -> Dict[str, Any]:
        """Analyze protein contact network and interactions.

        Args:
            structure: BioPython Structure object
            include_interactions: Whether to analyze specific interactions
            include_interfaces: Whether to analyze interfaces
            include_network: Whether to analyze network properties

        Returns:
            Dictionary of contact analysis results
        """
        try:
            # Build contact network
            network = self._build_contact_network(structure)
            if not network:
                return {}

            results = {"contact_network": network}

            # Analyze specific interactions
            if include_interactions:
                interactions = self._analyze_interactions(structure, network)
                results.update(interactions)

            # Analyze interfaces
            if include_interfaces:
                interfaces = self._analyze_interfaces(structure, network)
                results.update({"interfaces": interfaces})

            # Analyze network properties
            if include_network:
                network_props = self._analyze_network_properties(network)
                results.update(network_props)

            return results

        except Exception as e:
            self.logger.error(f"Error analyzing contacts: {str(e)}")
            return {}

    def analyze_site_contacts(
        self,
        site_residues: List[int],
        contact_network: Dict[int, Set[int]],
    ) -> Dict[str, Any]:
        """Analyze contacts for binding site residues.

        Args:
            site_residues: List of site residue numbers
            contact_network: Pre-calculated contact network

        Returns:
            Dictionary of site contact properties
        """
        try:
            if not site_residues or not contact_network:
                return {}

            # Count internal and external contacts
            internal_contacts = 0
            external_contacts = 0
            site_set = set(site_residues)

            for res_id in site_residues:
                if res_id in contact_network:
                    neighbors = contact_network[res_id]
                    internal_contacts += len(neighbors & site_set)
                    external_contacts += len(neighbors - site_set)

            # Adjust for double counting of internal contacts
            internal_contacts //= 2

            # Calculate contact metrics
            n_residues = len(site_residues)
            max_internal = n_residues * (n_residues - 1) / 2

            return {
                "internal_contacts": internal_contacts,
                "external_contacts": external_contacts,
                "contact_density": float(internal_contacts / max_internal) if max_internal > 0 else 0.0,
                "surface_exposure": float(external_contacts / (internal_contacts + external_contacts)) if (internal_contacts + external_contacts) > 0 else 1.0,
                "average_degree": float(sum(len(contact_network.get(res, set())) for res in site_residues) / n_residues) if n_residues > 0 else 0.0,
            }

        except Exception as e:
            self.logger.error(f"Error analyzing site contacts: {str(e)}")
            return {}

    def _build_contact_network(
        self,
        structure: Structure,
    ) -> Dict[int, Set[int]]:
        """Build residue contact network.

        Args:
            structure: BioPython Structure object

        Returns:
            Dictionary mapping residue numbers to sets of contacting residues
        """
        try:
            network = {}
            residues = list(structure.get_residues())

            # Get coordinates for each residue
            coords = {}
            for res in residues:
                res_id = res.get_id()[1]
                res_coords = []
                for atom in res:
                    res_coords.append(atom.get_coord())
                if res_coords:
                    coords[res_id] = np.array(res_coords)

            # Find contacts between residues
            res_ids = list(coords.keys())
            for i, res1_id in enumerate(res_ids):
                network[res1_id] = set()
                for res2_id in res_ids[i + 1 :]:
                    # Calculate minimum distance between residues
                    dist_matrix = cdist(coords[res1_id], coords[res2_id])
                    min_dist = np.min(dist_matrix)

                    if min_dist < self.settings["contact_cutoff"]:
                        network[res1_id].add(res2_id)
                        if res2_id not in network:
                            network[res2_id] = set()
                        network[res2_id].add(res1_id)

            return network

        except Exception as e:
            self.logger.error(f"Error building contact network: {str(e)}")
            return {}

    def _analyze_interactions(
        self,
        structure: Structure,
        network: Dict[int, Set[int]],
    ) -> Dict[str, Any]:
        """Analyze specific residue interactions.

        Args:
            structure: BioPython Structure object
            network: Contact network

        Returns:
            Dictionary of interaction analysis
        """
        try:
            interactions = {
                "hydrophobic": [],
                "hbonds": [],
                "ionic": [],
                "aromatic": [],
                "cation_pi": [],
                "disulfide": [],
            }

            residues = {res.get_id()[1]: res for res in structure.get_residues()}

            # Analyze each contacting pair
            for res1_id, neighbors in network.items():
                if res1_id not in residues:
                    continue
                res1 = residues[res1_id]

                for res2_id in neighbors:
                    if res2_id not in residues or res2_id <= res1_id:
                        continue
                    res2 = residues[res2_id]

                    # Check each interaction type
                    interaction = self._check_interaction(res1, res2)
                    for int_type, pairs in interaction.items():
                        if pairs:
                            interactions[int_type].append(
                                {
                                    "residues": (res1_id, res2_id),
                                    "type": int_type,
                                    "details": pairs,
                                }
                            )

            # Calculate interaction statistics
            stats = {f"{key}_count": len(value) for key, value in interactions.items()}

            return {
                "interactions": interactions,
                "interaction_stats": stats,
            }

        except Exception as e:
            self.logger.error(f"Error analyzing interactions: {str(e)}")
            return {}

    def analyze_water_bridges(
        self,
        structure: Structure,
        network: Dict[int, Set[int]],
    ) -> Dict[str, Any]:
        """Analyze water-mediated interactions.

        Args:
            structure: BioPython Structure object
            network: Contact network

        Returns:
            Dictionary of water bridge analysis
        """
        try:
            water_bridges = []
            cutoff = self.settings["water_bridge_cutoff"]

            # Get water molecules
            waters = []
            for atom in structure.get_atoms():
                if atom.get_name() == "O" and atom.get_parent().get_resname() == "HOH":
                    waters.append(atom)

            # Get protein atoms that can interact with water
            protein_atoms = []
            for residue in structure.get_residues():
                if residue.get_resname() != "HOH":
                    for atom in residue:
                        if atom.get_name() in ["N", "O", "OD1", "OD2", "OE1", "OE2", "NE", "NH1", "NH2", "ND1", "NE2"]:
                            protein_atoms.append(atom)

            # Find water bridges
            for water in waters:
                # Find protein atoms within cutoff of water
                close_atoms = []
                for atom in protein_atoms:
                    if water - atom < cutoff:
                        close_atoms.append(atom)

                # Check for bridges between different residues
                for i, atom1 in enumerate(close_atoms):
                    res1 = atom1.get_parent()
                    for atom2 in close_atoms[i + 1 :]:
                        res2 = atom2.get_parent()
                        if res1 != res2:
                            bridge = {
                                "residues": (res1.get_id()[1], res2.get_id()[1]),
                                "atoms": (atom1.get_name(), atom2.get_name()),
                                "water": water.get_full_id()[3][1],  # Water residue number
                                "distances": (
                                    float(water - atom1),
                                    float(water - atom2),
                                ),
                                "angle": float(
                                    calc_angle(
                                        atom1.get_vector(),
                                        water.get_vector(),
                                        atom2.get_vector(),
                                    )
                                ),
                            }
                            water_bridges.append(bridge)

            return {
                "water_bridges": water_bridges,
                "bridge_count": len(water_bridges),
                "bridged_residues": len(set(res_id for bridge in water_bridges for res_id in bridge["residues"])),
            }

        except Exception as e:
            self.logger.error(f"Error analyzing water bridges: {str(e)}")
            return {}

    def analyze_interface_hotspots(
        self,
        structure: Structure,
        interface_residues: Set[int],
        network: Dict[int, Set[int]],
    ) -> Dict[str, Any]:
        """Analyze interface hotspots using multiple criteria.

        Args:
            structure: BioPython Structure object
            interface_residues: Set of interface residue numbers
            network: Contact network

        Returns:
            Dictionary of hotspot analysis
        """
        try:
            hotspots = []
            residues = {res.get_id()[1]: res for res in structure.get_residues()}

            # Calculate properties for each interface residue
            for res_id in interface_residues:
                if res_id not in residues:
                    continue
                residue = residues[res_id]

                # Get interaction energy
                interaction_energy = 0.0
                if res_id in network:
                    for neighbor_id in network[res_id]:
                        if neighbor_id in residues:
                            neighbor = residues[neighbor_id]
                            energy = self._estimate_pairwise_energy(residue, neighbor)
                            interaction_energy += energy

                # Calculate burial
                burial = self._calculate_burial(residue, structure)

                # Calculate conservation if available
                conservation = self._calculate_conservation([residue])
                conservation_score = float(np.mean(list(conservation.values()))) if conservation else 0.0

                # Combine scores
                total_score = 0.4 * abs(interaction_energy) + 0.3 * burial + 0.3 * conservation_score

                if total_score > self.settings["hotspot_threshold"]:
                    hotspots.append(
                        {
                            "residue": res_id,
                            "score": float(total_score),
                            "energy": float(interaction_energy),
                            "burial": float(burial),
                            "conservation": float(conservation_score),
                        }
                    )

            return {
                "hotspots": sorted(hotspots, key=lambda x: x["score"], reverse=True),
                "total_hotspots": len(hotspots),
                "average_score": float(np.mean([h["score"] for h in hotspots])) if hotspots else 0.0,
            }

        except Exception as e:
            self.logger.error(f"Error analyzing hotspots: {str(e)}")
            return {}

    def analyze_allosteric_sites(
        self,
        structure: Structure,
        binding_site: Set[int],
        network: Dict[int, Set[int]],
    ) -> Dict[str, Any]:
        """Analyze potential allosteric sites using network and dynamics analysis.

        Args:
            structure: BioPython Structure object
            binding_site: Set of binding site residue numbers
            network: Contact network

        Returns:
            Dictionary of allosteric site analysis
        """
        try:
            # Create NetworkX graph with weighted edges
            G = nx.Graph()
            for res1_id, neighbors in network.items():
                for res2_id in neighbors:
                    # Weight edges by interaction strength and distance
                    weight = self._calculate_edge_weight(
                        structure,
                        res1_id,
                        res2_id,
                        binding_site,
                    )
                    G.add_edge(res1_id, res2_id, weight=weight)

            # Find communities
            communities = list(
                nx.community.greedy_modularity_communities(
                    G,
                    weight="weight",
                    resolution=self.settings["community_resolution"],
                )
            )

            # Analyze paths between binding site and other regions
            allosteric_sites = []
            binding_community = None
            for i, comm in enumerate(communities):
                if any(res_id in binding_site for res_id in comm):
                    binding_community = i
                    break

            if binding_community is not None:
                for i, comm in enumerate(communities):
                    if i != binding_community:
                        # Calculate properties of potential allosteric site
                        site_props = self._analyze_potential_site(
                            structure,
                            comm,
                            binding_site,
                            G,
                        )
                        if site_props["score"] > self.settings["hotspot_threshold"]:
                            allosteric_sites.append(
                                {
                                    "residues": sorted(list(comm)),
                                    "score": site_props["score"],
                                    "path_length": site_props["path_length"],
                                    "coupling_strength": site_props["coupling_strength"],
                                    "community_id": i,
                                }
                            )

            return {
                "allosteric_sites": sorted(
                    allosteric_sites,
                    key=lambda x: x["score"],
                    reverse=True,
                ),
                "total_sites": len(allosteric_sites),
                "communities": [sorted(list(c)) for c in communities],
                "binding_community": binding_community,
            }

        except Exception as e:
            self.logger.error(f"Error analyzing allosteric sites: {str(e)}")
            return {}

    def _calculate_edge_weight(
        self,
        structure: Structure,
        res1_id: int,
        res2_id: int,
        binding_site: Set[int],
    ) -> float:
        """Calculate weighted edge for network analysis.

        Args:
            structure: BioPython Structure object
            res1_id: First residue number
            res2_id: Second residue number
            binding_site: Set of binding site residue numbers

        Returns:
            Edge weight combining multiple factors
        """
        try:
            residues = {res.get_id()[1]: res for res in structure.get_residues()}
            if res1_id not in residues or res2_id not in residues:
                return 0.0

            res1 = residues[res1_id]
            res2 = residues[res2_id]

            # Interaction energy
            energy = abs(self._estimate_pairwise_energy(res1, res2))

            # Distance factor
            min_dist = float("inf")
            for atom1 in res1:
                for atom2 in res2:
                    dist = atom1 - atom2
                    min_dist = min(min_dist, dist)
            dist_factor = np.exp(-min_dist / self.settings["path_weight_factor"])

            # Binding site proximity
            site_dist = 0.0
            if res1_id in binding_site or res2_id in binding_site:
                site_dist = 1.0
            else:
                min_site_dist = float("inf")
                for site_id in binding_site:
                    if site_id in residues:
                        site_res = residues[site_id]
                        for atom1 in res1:
                            for atom2 in site_res:
                                dist = atom1 - atom2
                                min_site_dist = min(min_site_dist, dist)
                site_dist = np.exp(-min_site_dist / self.settings["allosteric_distance"])

            # Combine factors
            weight = 0.4 * energy + 0.4 * dist_factor + 0.2 * site_dist
            return float(weight)

        except Exception as e:
            self.logger.error(f"Error calculating edge weight: {str(e)}")
            return 0.0

    def _analyze_potential_site(
        self,
        structure: Structure,
        residues: Set[int],
        binding_site: Set[int],
        network: nx.Graph,
    ) -> Dict[str, float]:
        """Analyze properties of potential allosteric site.

        Args:
            structure: BioPython Structure object
            residues: Set of residue numbers in potential site
            binding_site: Set of binding site residue numbers
            network: NetworkX graph of protein structure

        Returns:
            Dictionary of site properties
        """
        try:
            # Calculate shortest paths to binding site
            path_lengths = []
            coupling_strengths = []

            for res_id in residues:
                min_length = float("inf")
                max_coupling = 0.0

                for site_id in binding_site:
                    if network.has_node(res_id) and network.has_node(site_id):
                        try:
                            path = nx.shortest_path(
                                network,
                                res_id,
                                site_id,
                                weight="weight",
                            )
                            path_lengths.append(len(path))

                            # Calculate coupling strength along path
                            coupling = 1.0
                            for i in range(len(path) - 1):
                                coupling *= network[path[i]][path[i + 1]]["weight"]
                            coupling_strengths.append(coupling)

                        except nx.NetworkXNoPath:
                            continue

            if not path_lengths:
                return {
                    "score": 0.0,
                    "path_length": float("inf"),
                    "coupling_strength": 0.0,
                }

            # Calculate overall score
            avg_length = float(np.mean(path_lengths))
            avg_coupling = float(np.mean(coupling_strengths))

            score = 0.5 * (1.0 / avg_length) + 0.5 * avg_coupling

            return {
                "score": score,
                "path_length": avg_length,
                "coupling_strength": avg_coupling,
            }

        except Exception as e:
            self.logger.error(f"Error analyzing potential site: {str(e)}")
            return {
                "score": 0.0,
                "path_length": float("inf"),
                "coupling_strength": 0.0,
            }

    def _check_interaction(
        self,
        res1: Residue,
        res2: Residue,
    ) -> Dict[str, List[Tuple[str, str]]]:
        """Check all possible interactions between residues.

        Args:
            res1: First residue
            res2: Second residue

        Returns:
            Dictionary mapping interaction types to atom pairs
        """
        try:
            interactions = {
                "hydrophobic": [],
                "hbonds": [],
                "ionic": [],
                "aromatic": [],
                "cation_pi": [],
                "disulfide": [],
                "pi_stacking": [],
                "halogen_bond": [],
                "metal_coord": [],
            }

            res1_name = res1.get_resname()
            res2_name = res2.get_resname()

            # Hydrophobic interactions
            if res1_name in self.hydrophobic and res2_name in self.hydrophobic:
                hydrophobic = self._check_hydrophobic(res1, res2)
                if hydrophobic:
                    interactions["hydrophobic"].extend(hydrophobic)

            # Hydrogen bonds
            if (res1_name in self.hbond_donors and res2_name in self.hbond_acceptors) or (res2_name in self.hbond_donors and res1_name in self.hbond_acceptors):
                hbonds = self._check_hbonds(res1, res2)
                if hbonds:
                    interactions["hbonds"].extend(hbonds)

            # Ionic interactions
            if (res1_name in self.charged_pos and res2_name in self.charged_neg) or (res1_name in self.charged_neg and res2_name in self.charged_pos):
                ionic = self._check_ionic(res1, res2)
                if ionic:
                    interactions["ionic"].extend(ionic)

            # Aromatic interactions
            if res1_name in self.aromatic and res2_name in self.aromatic:
                aromatic = self._check_aromatic(res1, res2)
                if aromatic:
                    interactions["aromatic"].extend(aromatic)

            # Cation-π interactions
            if (res1_name in self.charged_pos and res2_name in self.aromatic) or (res1_name in self.aromatic and res2_name in self.charged_pos):
                cation_pi = self._check_cation_pi(res1, res2)
                if cation_pi:
                    interactions["cation_pi"].extend(cation_pi)

            # Disulfide bonds
            if res1_name == "CYS" and res2_name == "CYS":
                disulfide = self._check_disulfide(res1, res2)
                if disulfide:
                    interactions["disulfide"].extend(disulfide)

            # π-stacking interactions
            if res1_name in self.pi_stacking and res2_name in self.pi_stacking:
                pi_stacking = self._check_pi_stacking(res1, res2)
                if pi_stacking:
                    interactions["pi_stacking"].extend(pi_stacking)

            # Halogen bonds
            if (res1_name in self.halogen_donors and res2_name in self.halogen_acceptors) or (res2_name in self.halogen_donors and res1_name in self.halogen_acceptors):
                halogen = self._check_halogen_bond(res1, res2)
                if halogen:
                    interactions["halogen_bond"].extend(halogen)

            # Metal coordination
            if res1_name in self.metal_coordinating and res2_name in self.metal_coordinating:
                metal = self._check_metal_coordination(res1, res2)
                if metal:
                    interactions["metal_coord"].extend(metal)

            return interactions

        except Exception as e:
            self.logger.error(f"Error checking interactions: {str(e)}")
            return {key: [] for key in interactions.keys()}

    def _check_hydrophobic(
        self,
        res1: Residue,
        res2: Residue,
    ) -> List[Tuple[str, str]]:
        """Check for hydrophobic interactions.

        Args:
            res1: First residue
            res2: Second residue

        Returns:
            List of interacting atom pairs
        """
        try:
            pairs = []
            cutoff = self.settings["hydrophobic_cutoff"]

            # Get carbon atoms
            carbons1 = [atom for atom in res1 if atom.get_name().startswith("C")]
            carbons2 = [atom for atom in res2 if atom.get_name().startswith("C")]

            # Check distances
            for atom1 in carbons1:
                for atom2 in carbons2:
                    if atom1 - atom2 < cutoff:
                        pairs.append((atom1.get_name(), atom2.get_name()))

            return pairs

        except Exception as e:
            self.logger.error(f"Error checking hydrophobic: {str(e)}")
            return []

    def _check_pi_stacking(
        self,
        res1: Residue,
        res2: Residue,
    ) -> List[Tuple[str, str]]:
        """Check for π-stacking interactions.

        Args:
            res1: First residue
            res2: Second residue

        Returns:
            List of stacking atom pairs
        """
        try:
            pairs = []
            cutoff = self.settings["pi_stacking_cutoff"]

            # Get ring atoms
            ring1 = self._get_aromatic_ring(res1)
            ring2 = self._get_aromatic_ring(res2)

            if ring1 and ring2:
                # Calculate ring centers and normals
                center1, normal1 = self._get_ring_geometry(ring1)
                center2, normal2 = self._get_ring_geometry(ring2)

                if center1 is not None and center2 is not None:
                    # Check distance between ring centers
                    dist = np.linalg.norm(center2 - center1)
                    if dist < cutoff:
                        # Calculate angle between ring normals
                        angle = np.degrees(np.arccos(abs(np.dot(normal1, normal2))))

                        # Parallel stacking (angle < 30°) or T-shaped (angle > 60°)
                        if angle < 30 or angle > 60:
                            pairs.extend([(atom1.get_name(), atom2.get_name()) for atom1, atom2 in zip(ring1, ring2)])

            return pairs

        except Exception as e:
            self.logger.error(f"Error checking pi-stacking: {str(e)}")
            return []

    def _check_halogen_bond(
        self,
        res1: Residue,
        res2: Residue,
    ) -> List[Tuple[str, str]]:
        """Check for halogen bonds.

        Args:
            res1: First residue
            res2: Second residue

        Returns:
            List of halogen bond pairs
        """
        try:
            pairs = []
            cutoff = self.settings["halogen_bond_cutoff"]

            # Check both directions
            for donor, acceptor in [(res1, res2), (res2, res1)]:
                if donor.get_resname() in self.halogen_donors and acceptor.get_resname() in self.halogen_acceptors:
                    # Get halogen atoms (C-X where X is halogen)
                    for atom in donor:
                        if atom.get_name().startswith(("CL", "BR", "I")):
                            # Get acceptor atoms (O, N)
                            for acc_atom in acceptor:
                                if acc_atom.get_name().startswith(("O", "N")):
                                    dist = atom - acc_atom
                                    if dist < cutoff:
                                        pairs.append((atom.get_name(), acc_atom.get_name()))

            return pairs

        except Exception as e:
            self.logger.error(f"Error checking halogen bonds: {str(e)}")
            return []

    def _check_metal_coordination(
        self,
        res1: Residue,
        res2: Residue,
    ) -> List[Tuple[str, str]]:
        """Check for metal coordination.

        Args:
            res1: First residue
            res2: Second residue

        Returns:
            List of coordinating atom pairs
        """
        try:
            pairs = []
            cutoff = self.settings["metal_coord_cutoff"]

            # Metal coordinating atoms
            coord_atoms = {
                "HIS": ["ND1", "NE2"],
                "CYS": ["SG"],
                "MET": ["SD"],
                "ASP": ["OD1", "OD2"],
                "GLU": ["OE1", "OE2"],
            }

            res1_name = res1.get_resname()
            res2_name = res2.get_resname()

            if res1_name in coord_atoms and res2_name in coord_atoms:
                for atom1_name in coord_atoms[res1_name]:
                    if atom1_name in res1:
                        for atom2_name in coord_atoms[res2_name]:
                            if atom2_name in res2:
                                dist = res1[atom1_name] - res2[atom2_name]
                                if dist < cutoff:
                                    pairs.append((atom1_name, atom2_name))

            return pairs

        except Exception as e:
            self.logger.error(f"Error checking metal coordination: {str(e)}")
            return []

    def _get_aromatic_ring(self, residue: Residue) -> List["Atom"]:
        """Get aromatic ring atoms from residue.

        Args:
            residue: Residue to analyze

        Returns:
            List of ring atoms or empty list if no ring found
        """
        try:
            # Ring atom patterns
            ring_atoms = {
                "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
                "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
                "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
                "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
            }

            res_name = residue.get_resname()
            if res_name in ring_atoms:
                ring = []
                for atom_name in ring_atoms[res_name]:
                    if atom_name in residue:
                        ring.append(residue[atom_name])
                if len(ring) >= 5:  # Need at least 5 atoms for a ring
                    return ring
            return []

        except Exception as e:
            self.logger.error(f"Error getting aromatic ring: {str(e)}")
            return []

    def _get_ring_geometry(
        self,
        ring_atoms: List["Atom"],
    ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """Calculate ring center and normal vector.

        Args:
            ring_atoms: List of ring atoms

        Returns:
            Tuple of (center coordinates, normal vector) or (None, None)
        """
        try:
            if len(ring_atoms) < 3:
                return None, None

            # Calculate ring center
            coords = np.array([atom.get_coord() for atom in ring_atoms])
            center = np.mean(coords, axis=0)

            # Calculate ring normal using cross product
            v1 = coords[1] - coords[0]
            v2 = coords[2] - coords[0]
            normal = np.cross(v1, v2)
            norm = np.linalg.norm(normal)
            if norm > 0:
                normal = normal / norm
                return center, normal
            return center, None

        except Exception as e:
            self.logger.error(f"Error calculating ring geometry: {str(e)}")
            return None, None

    def _check_hbonds(
        self,
        res1: Residue,
        res2: Residue,
    ) -> List[Tuple[str, str]]:
        """Check for hydrogen bonds.

        Args:
            res1: First residue
            res2: Second residue

        Returns:
            List of donor-acceptor atom pairs
        """
        try:
            pairs = []
            dist_cutoff = self.settings["hbond_distance"]
            angle_cutoff = self.settings["hbond_angle"]

            # Donor-acceptor pairs to check
            donors = {
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
            }

            acceptors = {"O", "OD1", "OD2", "OE1", "OE2", "ND1", "NE2"}

            # Check each donor-acceptor pair
            for donor_name, hydrogens in donors.items():
                if donor_name in res1:
                    donor = res1[donor_name]
                    for acceptor_name in acceptors:
                        if acceptor_name in res2:
                            acceptor = res2[acceptor_name]
                            if donor - acceptor < dist_cutoff:
                                # Check angle if hydrogen present
                                for h_name in hydrogens:
                                    if h_name in res1:
                                        hydrogen = res1[h_name]
                                        angle = self._calc_angle(hydrogen, donor, acceptor)
                                        if abs(180 - angle) < angle_cutoff:
                                            pairs.append((donor_name, acceptor_name))
                                            break

            # Check reverse direction
            for donor_name, hydrogens in donors.items():
                if donor_name in res2:
                    donor = res2[donor_name]
                    for acceptor_name in acceptors:
                        if acceptor_name in res1:
                            acceptor = res1[acceptor_name]
                            if donor - acceptor < dist_cutoff:
                                for h_name in hydrogens:
                                    if h_name in res2:
                                        hydrogen = res2[h_name]
                                        angle = self._calc_angle(hydrogen, donor, acceptor)
                                        if abs(180 - angle) < angle_cutoff:
                                            pairs.append((acceptor_name, donor_name))
                                            break

            return pairs

        except Exception as e:
            self.logger.error(f"Error checking H-bonds: {str(e)}")
            return []

    def _check_ionic(
        self,
        res1: Residue,
        res2: Residue,
    ) -> List[Tuple[str, str]]:
        """Check for ionic interactions.

        Args:
            res1: First residue
            res2: Second residue

        Returns:
            List of charged atom pairs
        """
        try:
            pairs = []
            cutoff = self.settings["ionic_cutoff"]

            # Charged groups
            positive = {
                "ARG": ["NH1", "NH2"],
                "LYS": ["NZ"],
                "HIS": ["ND1", "NE2"],
            }
            negative = {
                "ASP": ["OD1", "OD2"],
                "GLU": ["OE1", "OE2"],
            }

            # Check positive-negative pairs
            res1_name = res1.get_resname()
            res2_name = res2.get_resname()

            if res1_name in positive and res2_name in negative:
                for pos_atom in positive[res1_name]:
                    if pos_atom in res1:
                        for neg_atom in negative[res2_name]:
                            if neg_atom in res2:
                                if res1[pos_atom] - res2[neg_atom] < cutoff:
                                    pairs.append((pos_atom, neg_atom))

            elif res1_name in negative and res2_name in positive:
                for neg_atom in negative[res1_name]:
                    if neg_atom in res1:
                        for pos_atom in positive[res2_name]:
                            if pos_atom in res2:
                                if res1[neg_atom] - res2[pos_atom] < cutoff:
                                    pairs.append((neg_atom, pos_atom))

            return pairs

        except Exception as e:
            self.logger.error(f"Error checking ionic: {str(e)}")
            return []

    def _check_aromatic(
        self,
        res1: Residue,
        res2: Residue,
    ) -> List[Tuple[str, str]]:
        """Check for aromatic interactions.

        Args:
            res1: First residue
            res2: Second residue

        Returns:
            List of ring atom pairs
        """
        try:
            pairs = []
            cutoff = self.settings["aromatic_cutoff"]

            # Ring atoms
            rings = {
                "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
                "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
                "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
                "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
            }

            res1_name = res1.get_resname()
            res2_name = res2.get_resname()

            if res1_name in rings and res2_name in rings:
                ring1 = [atom for name in rings[res1_name] if name in res1]
                ring2 = [atom for name in rings[res2_name] if name in res2]

                if len(ring1) >= 5 and len(ring2) >= 5:
                    # Calculate ring centers
                    center1 = np.mean([atom.get_coord() for atom in ring1], axis=0)
                    center2 = np.mean([atom.get_coord() for atom in ring2], axis=0)

                    # Calculate ring normals
                    normal1 = self._calc_ring_normal(ring1)
                    normal2 = self._calc_ring_normal(ring2)

                    if normal1 is not None and normal2 is not None:
                        # Check distance and angle
                        dist = np.linalg.norm(center2 - center1)
                        angle = np.degrees(np.arccos(np.abs(np.dot(normal1, normal2))))

                        if dist < cutoff and (angle < 30 or angle > 150):
                            pairs.extend([(atom1.get_name(), atom2.get_name()) for atom1, atom2 in zip(ring1, ring2)])

            return pairs

        except Exception as e:
            self.logger.error(f"Error checking aromatic: {str(e)}")
            return []

    def _check_cation_pi(
        self,
        res1: Residue,
        res2: Residue,
    ) -> List[Tuple[str, str]]:
        """Check for cation-π interactions.

        Args:
            res1: First residue
            res2: Second residue

        Returns:
            List of cation-ring atom pairs
        """
        try:
            pairs = []
            cutoff = self.settings["cation_pi_cutoff"]

            # Ring atoms
            rings = {
                "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
                "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
                "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
            }

            # Cation atoms
            cations = {
                "ARG": ["NH1", "NH2"],
                "LYS": ["NZ"],
            }

            res1_name = res1.get_resname()
            res2_name = res2.get_resname()

            # Check ARG/LYS - aromatic
            if res1_name in cations and res2_name in rings:
                cation_atoms = [atom for name in cations[res1_name] if name in res1]
                ring_atoms = [atom for name in rings[res2_name] if name in res2]

                if ring_atoms:
                    ring_center = np.mean([atom.get_coord() for atom in ring_atoms], axis=0)
                    ring_normal = self._calc_ring_normal(ring_atoms)

                    if ring_normal is not None:
                        for cation in cation_atoms:
                            cation_pos = cation.get_coord()
                            dist = np.linalg.norm(cation_pos - ring_center)
                            if dist < cutoff:
                                pairs.extend([(cation.get_name(), atom.get_name()) for atom in ring_atoms])

            # Check aromatic - ARG/LYS
            elif res1_name in rings and res2_name in cations:
                ring_atoms = [atom for name in rings[res1_name] if name in res1]
                cation_atoms = [atom for name in cations[res2_name] if name in res2]

                if ring_atoms:
                    ring_center = np.mean([atom.get_coord() for atom in ring_atoms], axis=0)
                    ring_normal = self._calc_ring_normal(ring_atoms)

                    if ring_normal is not None:
                        for cation in cation_atoms:
                            cation_pos = cation.get_coord()
                            dist = np.linalg.norm(cation_pos - ring_center)
                            if dist < cutoff:
                                pairs.extend([(atom.get_name(), cation.get_name()) for atom in ring_atoms])

            return pairs

        except Exception as e:
            self.logger.error(f"Error checking cation-pi: {str(e)}")
            return []

    def _check_disulfide(
        self,
        res1: Residue,
        res2: Residue,
    ) -> List[Tuple[str, str]]:
        """Check for disulfide bonds.

        Args:
            res1: First residue
            res2: Second residue

        Returns:
            List of sulfur atom pairs
        """
        try:
            pairs = []
            cutoff = self.settings["disulfide_cutoff"]

            if "SG" in res1 and "SG" in res2:
                if res1["SG"] - res2["SG"] < cutoff:
                    pairs.append(("SG", "SG"))

            return pairs

        except Exception as e:
            self.logger.error(f"Error checking disulfide: {str(e)}")
            return []

    def _analyze_interfaces(
        self,
        structure: Structure,
        network: Dict[int, Set[int]],
    ) -> List[Dict[str, Any]]:
        """Analyze interfaces between protein regions.

        Args:
            structure: BioPython Structure object
            network: Contact network

        Returns:
            List of interface properties
        """
        try:
            interfaces = []
            residues = {res.get_id()[1]: res for res in structure.get_residues()}

            # Find connected components
            G = nx.Graph()
            for res1_id, neighbors in network.items():
                for res2_id in neighbors:
                    G.add_edge(res1_id, res2_id)

            components = list(nx.connected_components(G))

            # Analyze interfaces between components
            for i, comp1 in enumerate(components):
                for comp2 in components[i + 1 :]:
                    interface_residues = set()

                    # Find residues at interface
                    for res1_id in comp1:
                        if res1_id in network:
                            neighbors = network[res1_id]
                            if neighbors & comp2:
                                interface_residues.add(res1_id)
                                interface_residues.update(neighbors & comp2)

                    if len(interface_residues) >= self.settings["min_interface_contacts"]:
                        # Analyze interface properties
                        interface = {
                            "residues": sorted(list(interface_residues)),
                            "size": len(interface_residues),
                            "components": (sorted(list(comp1)), sorted(list(comp2))),
                        }

                        # Get interface interactions
                        interactions = []
                        for res1_id in interface_residues & comp1:
                            if res1_id not in residues:
                                continue
                            res1 = residues[res1_id]

                            for res2_id in interface_residues & comp2:
                                if res2_id not in residues:
                                    continue
                                res2 = residues[res2_id]

                                interaction = self._check_interaction(res1, res2)
                                for int_type, pairs in interaction.items():
                                    if pairs:
                                        interactions.append(
                                            {
                                                "residues": (res1_id, res2_id),
                                                "type": int_type,
                                                "details": pairs,
                                            }
                                        )

                        interface["interactions"] = interactions
                        interfaces.append(interface)

            return interfaces

        except Exception as e:
            self.logger.error(f"Error analyzing interfaces: {str(e)}")
            return []

    def _analyze_network_properties(
        self,
        network: Dict[int, Set[int]],
    ) -> Dict[str, Any]:
        """Analyze contact network properties.

        Args:
            network: Contact network

        Returns:
            Dictionary of network properties
        """
        try:
            # Create NetworkX graph
            G = nx.Graph()
            for res1_id, neighbors in network.items():
                for res2_id in neighbors:
                    G.add_edge(res1_id, res2_id)

            properties = {
                "n_nodes": G.number_of_nodes(),
                "n_edges": G.number_of_edges(),
                "density": nx.density(G),
                "average_degree": float(sum(dict(G.degree()).values()) / G.number_of_nodes()),
                "clustering_coefficient": nx.average_clustering(G),
                "average_path_length": float(nx.average_shortest_path_length(G)) if nx.is_connected(G) else float("inf"),
            }

            # Calculate degree distribution
            degrees = dict(G.degree())
            properties["degree_distribution"] = {d: list(degrees.values()).count(d) for d in sorted(set(degrees.values()))}

            # Identify hubs (high degree nodes)
            mean_degree = properties["average_degree"]
            std_degree = float(np.std(list(degrees.values())))
            properties["hubs"] = [node for node, degree in degrees.items() if degree > mean_degree + 2 * std_degree]

            # Calculate centrality measures
            properties["centrality"] = {
                "degree": nx.degree_centrality(G),
                "betweenness": nx.betweenness_centrality(G),
                "closeness": nx.closeness_centrality(G),
                "eigenvector": nx.eigenvector_centrality_numpy(G),
            }

            # Calculate path-based metrics
            properties["paths"] = {
                "diameter": float(nx.diameter(G)),
                "radius": float(nx.radius(G)),
                "center": sorted(nx.center(G)),
                "periphery": sorted(nx.periphery(G)),
                "eccentricity": {node: float(ecc) for node, ecc in nx.eccentricity(G).items()},
            }

            # Identify communities
            communities = list(nx.community.greedy_modularity_communities(G))
            properties["communities"] = [sorted(list(community)) for community in communities]
            properties["modularity"] = float(nx.community.modularity(G, communities))

            # Calculate additional metrics
            properties.update(
                {
                    "transitivity": float(nx.transitivity(G)),
                    "graph_clique_number": float(nx.graph_clique_number(G)),
                    "number_of_cliques": float(nx.graph_number_of_cliques(G)),
                    "degree_assortativity": float(nx.degree_assortativity_coefficient(G)),
                    "global_efficiency": float(nx.global_efficiency(G)),
                    "local_efficiency": float(nx.local_efficiency(G)),
                }
            )

            return properties

        except Exception as e:
            self.logger.error(f"Error analyzing network properties: {str(e)}")
            return {}

    def _calc_angle(
        self,
        atom1: "Atom",
        atom2: "Atom",
        atom3: "Atom",
    ) -> float:
        """Calculate angle between three atoms.

        Args:
            atom1: First atom
            atom2: Central atom
            atom3: Third atom

        Returns:
            Angle in degrees
        """
        try:
            v1 = atom1.get_coord() - atom2.get_coord()
            v2 = atom3.get_coord() - atom2.get_coord()

            cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            return float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))

        except Exception as e:
            self.logger.error(f"Error calculating angle: {str(e)}")
            return 0.0

    def _calc_ring_normal(self, ring_atoms: List["Atom"]) -> Optional[np.ndarray]:
        """Calculate normal vector for aromatic ring.

        Args:
            ring_atoms: List of ring atoms

        Returns:
            Normal vector or None if calculation fails
        """
        try:
            if len(ring_atoms) < 3:
                return None

            # Get coordinates
            coords = np.array([atom.get_coord() for atom in ring_atoms])

            # Calculate normal using cross product of two vectors in ring
            v1 = coords[1] - coords[0]
            v2 = coords[2] - coords[0]
            normal = np.cross(v1, v2)

            # Normalize
            norm = np.linalg.norm(normal)
            if norm > 0:
                normal = normal / norm
                return normal
            return None

        except Exception as e:
            self.logger.error(f"Error calculating ring normal: {str(e)}")
            return None


def _estimate_pairwise_energy(
    self,
    res1: Residue,
    res2: Residue,
) -> float:
    """Estimate interaction energy between residue pair.

    Args:
        res1: First residue
        res2: Second residue

    Returns:
        Estimated interaction energy in arbitrary units
    """
    try:
        energy = 0.0

        # Get interaction types
        interactions = self._check_interaction(res1, res2)

        # Energy contributions
        energy_terms = {
            "hydrophobic": -0.8,  # Hydrophobic contact
            "hbonds": -2.0,  # H-bond
            "ionic": -3.0,  # Salt bridge
            "aromatic": -2.5,  # π-π stacking
            "cation_pi": -2.5,  # Cation-π
            "disulfide": -4.0,  # Disulfide bond
            "pi_stacking": -2.5,  # π-stacking
            "halogen_bond": -1.5,  # Halogen bond
            "metal_coord": -3.5,  # Metal coordination
        }

        # Sum contributions
        for int_type, pairs in interactions.items():
            if pairs:
                base_energy = energy_terms.get(int_type, 0.0)
                energy += base_energy * len(pairs)

        return float(energy)

    except Exception as e:
        self.logger.error(f"Error estimating pairwise energy: {str(e)}")
        return 0.0


def _calculate_burial(
    self,
    residue: Residue,
    structure: Structure,
) -> float:
    """Calculate burial score for residue.

    Args:
        residue: Residue to analyze
        structure: Full structure for context

    Returns:
        Burial score (0-1)
    """
    try:
        # Get residue center
        center = self._get_residue_center(residue)

        # Get surface atoms
        surface_atoms = self._get_surface_atoms(structure)

        if len(surface_atoms) == 0:
            return 0.0

        # Calculate minimum distance to surface
        distances = np.linalg.norm(surface_atoms - center, axis=1)
        min_dist = float(np.min(distances))

        # Calculate burial score (1 = fully buried)
        burial = 1.0 - np.exp(-min_dist / 5.0)  # 5Å characteristic length
        return float(burial)

    except Exception as e:
        self.logger.error(f"Error calculating burial: {str(e)}")
        return 0.0


def _calculate_conservation(
    self,
    residues: List[Residue],
) -> Dict[int, float]:
    """Calculate conservation scores for residues.

    Args:
        residues: List of residues

    Returns:
        Dictionary mapping residue IDs to conservation scores
    """
    try:
        conservation = {}

        # Get sequence window for each residue
        window_size = 5
        for i, residue in enumerate(residues):
            # Get surrounding residues
            window = residues[max(0, i - window_size // 2) : min(len(residues), i + window_size // 2 + 1)]

            # Calculate conservation score
            if window:
                # Get sequence
                sequence = ""
                for res in window:
                    sequence += three_to_one(res.get_resname())

                # Calculate conservation using BLOSUM matrix
                score = self._calculate_window_conservation(sequence)
                conservation[residue.get_id()[1]] = score

        return conservation

    except Exception as e:
        self.logger.error(f"Error calculating conservation: {str(e)}")
        return {}


def _analyze_site_domain_context(
    self,
    residues: List[Residue],
    structure: Structure,
) -> Dict[str, Any]:
    """Analyze domain context of binding site residues.

    Args:
        residues: List of residues in binding site
        structure: Full structure for context

    Returns:
        Dictionary of domain context properties
    """
    try:
        # Identify domains
        domains = self._identify_domains(structure)
        if not domains:
            return {}

        # Find domain membership
        site_domains = set()
        domain_coverage = {}
        for i, domain in enumerate(domains):
            overlap = set(res.get_id()[1] for res in residues) & set(domain)
            if overlap:
                site_domains.add(i)
                domain_coverage[i] = len(overlap) / len(residues)

        # Calculate domain interface properties
        interface_props = {}
        if len(site_domains) > 1:
            for i, j in itertools.combinations(site_domains, 2):
                interface = self._get_domain_interface(
                    domains[i],
                    domains[j],
                    structure,
                )
                if interface:
                    key = f"domain_{i}_{j}"
                    interface_props[key] = self._analyze_interface_properties(
                        interface,
                        structure,
                    )

        return {
            "domains": sorted(list(site_domains)),
            "domain_coverage": domain_coverage,
            "interfaces": interface_props,
            "spans_domains": len(site_domains) > 1,
            "domain_count": len(site_domains),
        }

    except Exception as e:
        self.logger.error(f"Error analyzing domain context: {str(e)}")
        return {}


def _get_domain_interface(
    self,
    domain1: List[int],
    domain2: List[int],
    structure: Structure,
) -> List[int]:
    """Get interface residues between domains.

    Args:
        domain1: First domain residue numbers
        domain2: Second domain residue numbers
        structure: Structure containing domains

    Returns:
        List of interface residue numbers
    """
    try:
        interface = set()
        cutoff = self.settings["interface_cutoff"]

        # Get residues
        residues = {res.get_id()[1]: res for res in structure.get_residues()}

        # Check contacts between domains
        for res1_id in domain1:
            if res1_id not in residues:
                continue
            res1 = residues[res1_id]

            for res2_id in domain2:
                if res2_id not in residues:
                    continue
                res2 = residues[res2_id]

                # Check minimum atom distance
                min_dist = float("inf")
                for atom1 in res1:
                    for atom2 in res2:
                        dist = atom1 - atom2
                        min_dist = min(min_dist, dist)

                if min_dist < cutoff:
                    interface.add(res1_id)
                    interface.add(res2_id)

        return sorted(list(interface))

    except Exception as e:
        self.logger.error(f"Error getting domain interface: {str(e)}")
        return []


def _analyze_interface_properties(
    self,
    interface_residues: List[int],
    structure: Structure,
) -> Dict[str, Any]:
    """Analyze properties of domain interface.

    Args:
        interface_residues: List of interface residue numbers
        structure: Structure containing interface

    Returns:
        Dictionary of interface properties
    """
    try:
        # Get interface residues
        residues = {res.get_id()[1]: res for res in structure.get_residues()}
        interface = [residues[res_id] for res_id in interface_residues if res_id in residues]

        if not interface:
            return {}

        # Calculate properties
        properties = {
            "size": len(interface),
            "hydrophobicity": float(np.mean([self.HYDROPHOBICITY.get(res.get_resname(), 0.0) for res in interface])),
            "charged_residues": len([res for res in interface if res.get_resname() in self.charged_pos | self.charged_neg]),
            "aromatic_residues": len([res for res in interface if res.get_resname() in self.aromatic]),
            "polar_residues": len([res for res in interface if res.get_resname() in ["SER", "THR", "ASN", "GLN"]]),
        }

        # Add interaction analysis
        interactions = []
        for i, res1 in enumerate(interface):
            for res2 in interface[i + 1 :]:
                interaction = self._check_interaction(res1, res2)
                if any(pairs for pairs in interaction.values()):
                    interactions.append(
                        {
                            "residues": (res1.get_id()[1], res2.get_id()[1]),
                            "types": [k for k, v in interaction.items() if v],
                        }
                    )

        properties["interactions"] = interactions
        properties["interaction_density"] = len(interactions) / (len(interface) * (len(interface) - 1) / 2) if len(interface) > 1 else 0.0

        return properties

    except Exception as e:
        self.logger.error(f"Error analyzing interface properties: {str(e)}")
        return {}


def _analyze_allosteric_pathways(
    self,
    structure: Structure,
    site1_residues: Set[int],
    site2_residues: Set[int],
    network: nx.Graph,
) -> Dict[str, Any]:
    """Analyze allosteric communication pathways between two sites.

    Args:
        structure: BioPython Structure object
        site1_residues: First site residue numbers
        site2_residues: Second site residue numbers
        network: NetworkX graph of protein structure

    Returns:
        Dictionary of pathway properties
    """
    try:
        pathways = []
        coupling_scores = []

        # Find all paths between sites
        for res1 in site1_residues:
            for res2 in site2_residues:
                if network.has_node(res1) and network.has_node(res2):
                    try:
                        # Get all simple paths up to certain length
                        paths = list(nx.all_simple_paths(network, res1, res2, cutoff=10))

                        for path in paths:
                            # Calculate pathway coupling
                            coupling = self._calculate_pathway_coupling(structure, path, network)

                            pathways.append(
                                {
                                    "residues": path,
                                    "length": len(path),
                                    "coupling": coupling,
                                }
                            )
                            coupling_scores.append(coupling)

                    except nx.NetworkXNoPath:
                        continue

        if not pathways:
            return {
                "pathways": [],
                "average_coupling": 0.0,
                "critical_residues": [],
            }

        # Identify critical residues
        critical_residues = self._identify_critical_residues(pathways, network)

        return {
            "pathways": sorted(pathways, key=lambda x: x["coupling"], reverse=True),
            "average_coupling": float(np.mean(coupling_scores)),
            "critical_residues": critical_residues,
        }

    except Exception as e:
        self.logger.error(f"Error analyzing allosteric pathways: {str(e)}")
        return {
            "pathways": [],
            "average_coupling": 0.0,
            "critical_residues": [],
        }


def _calculate_pathway_coupling(
    self,
    structure: Structure,
    path: List[int],
    network: nx.Graph,
) -> float:
    """Calculate coupling strength along communication pathway.

    Args:
        structure: BioPython Structure object
        path: List of residue numbers in pathway
        network: NetworkX graph

    Returns:
        Coupling strength score
    """
    try:
        if len(path) < 2:
            return 0.0

        coupling = 1.0
        residues = {res.get_id()[1]: res for res in structure.get_residues()}

        # Calculate pairwise coupling along path
        for i in range(len(path) - 1):
            res1_id = path[i]
            res2_id = path[i + 1]

            if res1_id not in residues or res2_id not in residues:
                continue

            res1 = residues[res1_id]
            res2 = residues[res2_id]

            # Get interaction energy
            energy = abs(self._estimate_pairwise_energy(res1, res2))

            # Get dynamic correlation
            correlation = self._calculate_dynamic_correlation(res1, res2)

            # Combine scores
            pair_coupling = 0.6 * energy + 0.4 * correlation
            coupling *= pair_coupling

        # Normalize by path length
        return float(coupling ** (1.0 / (len(path) - 1)))

    except Exception as e:
        self.logger.error(f"Error calculating pathway coupling: {str(e)}")
        return 0.0


def _identify_critical_residues(
    self,
    pathways: List[Dict[str, Any]],
    network: nx.Graph,
) -> List[Dict[str, Any]]:
    """Identify critical residues in allosteric pathways.

    Args:
        pathways: List of pathway dictionaries
        network: NetworkX graph

    Returns:
        List of critical residue properties
    """
    try:
        # Count residue occurrences in pathways
        residue_counts = {}
        for pathway in pathways:
            for res_id in pathway["residues"]:
                if res_id not in residue_counts:
                    residue_counts[res_id] = {
                        "count": 0,
                        "total_coupling": 0.0,
                    }
                residue_counts[res_id]["count"] += 1
                residue_counts[res_id]["total_coupling"] += pathway["coupling"]

        # Calculate betweenness centrality
        betweenness = nx.betweenness_centrality(network)

        # Identify critical residues
        critical_residues = []
        for res_id, counts in residue_counts.items():
            if counts["count"] > len(pathways) * 0.2:  # In >20% of pathways
                critical_residues.append(
                    {
                        "residue": res_id,
                        "pathway_count": counts["count"],
                        "average_coupling": counts["total_coupling"] / counts["count"],
                        "betweenness": betweenness.get(res_id, 0.0),
                    }
                )

        return sorted(critical_residues, key=lambda x: x["average_coupling"] * x["betweenness"], reverse=True)

    except Exception as e:
        self.logger.error(f"Error identifying critical residues: {str(e)}")
        return []


def _calculate_dynamic_correlation(
    self,
    res1: Residue,
    res2: Residue,
) -> float:
    """Calculate dynamic correlation between residues.

    Args:
        res1: First residue
        res2: Second residue

    Returns:
        Correlation score (0-1)
    """
    try:
        # Get CA atoms
        if "CA" not in res1 or "CA" not in res2:
            return 0.0

        ca1 = res1["CA"]
        ca2 = res2["CA"]

        # Get B-factors
        b1 = ca1.get_bfactor()
        b2 = ca2.get_bfactor()

        # Calculate correlation based on B-factors
        correlation = 1.0 - abs(b1 - b2) / max(b1, b2)

        # Distance weighting
        dist = ca1 - ca2
        dist_weight = np.exp(-dist / 10.0)  # 10Å characteristic length

        return float(correlation * dist_weight)

    except Exception as e:
        self.logger.error(f"Error calculating dynamic correlation: {str(e)}")
        return 0.0


def _analyze_dynamic_coupling(
    self,
    structure: Structure,
    network: nx.Graph,
) -> Dict[str, Any]:
    """Analyze dynamic coupling between residues.

    Args:
        structure: BioPython Structure object
        network: NetworkX graph

    Returns:
        Dictionary of dynamic coupling properties
    """
    try:
        # Calculate correlation matrix
        correlation_matrix = self._calculate_correlation_matrix(structure)

        # Identify dynamic domains
        domains = self._identify_dynamic_domains(correlation_matrix, network)

        # Calculate domain coupling
        domain_coupling = {}
        for i, domain1 in enumerate(domains):
            for j, domain2 in enumerate(domains[i + 1 :], i + 1):
                coupling = self._calculate_domain_coupling(domain1, domain2, correlation_matrix)
                domain_coupling[f"domain_{i+1}_{j+1}"] = coupling

        return {
            "correlation_matrix": correlation_matrix.tolist(),
            "domains": [sorted(list(d)) for d in domains],
            "domain_coupling": domain_coupling,
        }

    except Exception as e:
        self.logger.error(f"Error analyzing dynamic coupling: {str(e)}")
        return {}


def _calculate_correlation_matrix(
    self,
    structure: Structure,
) -> np.ndarray:
    """Calculate residue correlation matrix.

    Args:
        structure: BioPython Structure object

    Returns:
        Correlation matrix
    """
    try:
        residues = list(structure.get_residues())
        n_res = len(residues)
        matrix = np.zeros((n_res, n_res))

        for i, res1 in enumerate(residues):
            for j, res2 in enumerate(residues[i:], i):
                corr = self._calculate_dynamic_correlation(res1, res2)
                matrix[i, j] = corr
                matrix[j, i] = corr

        return matrix

    except Exception as e:
        self.logger.error(f"Error calculating correlation matrix: {str(e)}")
        return np.array([])


def _identify_dynamic_domains(
    self,
    correlation_matrix: np.ndarray,
    network: nx.Graph,
) -> List[Set[int]]:
    """Identify dynamically coupled domains.

    Args:
        correlation_matrix: Residue correlation matrix
        network: NetworkX graph

    Returns:
        List of domain residue sets
    """
    try:
        if len(correlation_matrix) == 0:
            return []

        # Create graph with correlation-weighted edges
        G = nx.Graph()
        n_res = len(correlation_matrix)

        for i in range(n_res):
            for j in range(i + 1, n_res):
                if correlation_matrix[i, j] > 0.5:  # Correlation threshold
                    G.add_edge(i, j, weight=correlation_matrix[i, j])

        # Find communities using correlation weights
        communities = list(nx.community.greedy_modularity_communities(G, weight="weight", resolution=1.2))

        return communities

    except Exception as e:
        self.logger.error(f"Error identifying dynamic domains: {str(e)}")
        return []


def _calculate_domain_coupling(
    self,
    domain1: Set[int],
    domain2: Set[int],
    correlation_matrix: np.ndarray,
) -> float:
    """Calculate coupling between dynamic domains.

    Args:
        domain1: First domain residue indices
        domain2: Second domain residue indices
        correlation_matrix: Residue correlation matrix

    Returns:
        Domain coupling score
    """
    try:
        if len(correlation_matrix) == 0:
            return 0.0

        # Calculate average correlation between domains
        correlations = []
        for i in domain1:
            for j in domain2:
                correlations.append(correlation_matrix[i, j])

        if not correlations:
            return 0.0

        return float(np.mean(correlations))

    except Exception as e:
        self.logger.error(f"Error calculating domain coupling: {str(e)}")
        return 0.0
