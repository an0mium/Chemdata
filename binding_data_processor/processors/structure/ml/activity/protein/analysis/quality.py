"""Protein structure quality analysis module."""

import logging
from typing import Dict, List, Any, Optional, Union, Tuple
import numpy as np
from Bio.PDB import Structure, Residue, Atom
from Bio.PDB.vectors import calc_angle, calc_dihedral
from Bio.PDB.Polypeptide import is_aa

from binding_data_processor.processors.structure.ml.activity.binding.alphafold import AlphaFoldIntegrator

logger = logging.getLogger(__name__)

# Initialize AlphaFold integrator
alphafold = AlphaFoldIntegrator()

# Standard amino acid bond lengths and angles
IDEAL_BOND_LENGTHS = {"N-CA": 1.46, "CA-C": 1.52, "C-N": 1.33, "CA-CB": 1.53}  # Angstroms

IDEAL_BOND_ANGLES = {"N-CA-C": 111.0, "CA-C-N": 116.6, "C-N-CA": 121.9}  # Degrees

RAMACHANDRAN_REGIONS = {
    "alpha": [(-60, -40), (-50, -30)],  # phi, psi ranges for alpha helix
    "beta": [(-140, -110), (130, 150)],  # phi, psi ranges for beta sheet
}


def calculate_quality_metrics(structure: Structure, include_alphafold_metrics: bool = True) -> Dict[str, Any]:
    """Calculate various structure quality metrics.

    Args:
        structure: BioPython Structure object
        include_alphafold_metrics: Whether to include AlphaFold-specific quality metrics

    Returns:
        Dictionary containing quality metrics
    """
    try:
        metrics = {}

        # Calculate geometry metrics
        geometry = _analyze_geometry(structure)
        metrics["geometry"] = geometry

        # Calculate Ramachandran statistics
        rama = _analyze_ramachandran(structure)
        metrics["ramachandran"] = rama

        # Calculate packing quality
        packing = _analyze_packing(structure)
        metrics["packing"] = packing

        # Calculate overall quality score
        metrics["overall_score"] = _calculate_overall_score(geometry, rama, packing)

        # Add AlphaFold-specific metrics if requested
        if include_alphafold_metrics:
            try:
                af_metrics = _calculate_alphafold_metrics(structure)
                metrics["alphafold"] = af_metrics

                # Update overall score to include AlphaFold confidence
                if "mean_plddt" in af_metrics:
                    metrics["overall_score"] = 0.6 * metrics["overall_score"] + 0.4 * (af_metrics["mean_plddt"] / 100)
            except Exception as e:
                logger.warning(f"Could not calculate AlphaFold metrics: {str(e)}")

        return metrics

    except Exception as e:
        logger.error(f"Error calculating quality metrics: {str(e)}")
        return {}


def analyze_clashes(structure: Structure, cutoff: float = 2.0) -> Dict[str, Any]:
    """Analyze steric clashes in structure.

    Args:
        structure: BioPython Structure object
        cutoff: Distance cutoff for clash detection in Angstroms

    Returns:
        Dictionary containing clash analysis results
    """
    try:
        clashes = []
        atoms = list(structure.get_atoms())

        for i, atom1 in enumerate(atoms):
            for atom2 in atoms[i + 1 :]:
                # Skip bonded atoms
                if _are_bonded(atom1, atom2):
                    continue

                distance = atom1 - atom2
                vdw_sum = get_vdw_radius(atom1) + get_vdw_radius(atom2)

                if distance < (vdw_sum - cutoff):
                    clashes.append(
                        {
                            "atom1": f"{atom1.get_parent().get_resname()}{atom1.get_parent().get_id()[1]}.{atom1.get_name()}",
                            "atom2": f"{atom2.get_parent().get_resname()}{atom2.get_parent().get_id()[1]}.{atom2.get_name()}",
                            "distance": float(distance),
                            "overlap": float(vdw_sum - distance),
                        }
                    )

        return {
            "n_clashes": len(clashes),
            "clashes": clashes,
            "max_overlap": max([c["overlap"] for c in clashes]) if clashes else 0.0,
            "mean_overlap": np.mean([c["overlap"] for c in clashes]) if clashes else 0.0,
        }

    except Exception as e:
        logger.error(f"Error analyzing clashes: {str(e)}")
        return {}


def analyze_rotamers(structure: Structure) -> Dict[str, Any]:
    """Analyze side chain rotamer quality.

    Args:
        structure: BioPython Structure object

    Returns:
        Dictionary containing rotamer analysis results
    """
    try:
        rotamer_scores = {}
        outliers = []

        for residue in structure.get_residues():
            if not is_aa(residue):
                continue

            chi_angles = calculate_chi_angles(residue)
            if not chi_angles:
                continue

            score = evaluate_rotamer(residue.get_resname(), chi_angles)
            rotamer_scores[f"{residue.get_resname()}{residue.get_id()[1]}"] = score

            if score < 0.3:  # Arbitrary threshold for outliers
                outliers.append({"residue": f"{residue.get_resname()}{residue.get_id()[1]}", "score": score, "chi_angles": chi_angles})

        return {"mean_score": float(np.mean(list(rotamer_scores.values()))), "n_outliers": len(outliers), "outliers": outliers, "scores": rotamer_scores}

    except Exception as e:
        logger.error(f"Error analyzing rotamers: {str(e)}")
        return {}


def _analyze_geometry(structure: Structure) -> Dict[str, Any]:
    """Analyze geometric quality of structure.

    Args:
        structure: BioPython Structure object

    Returns:
        Dictionary containing geometry analysis results
    """
    try:
        bond_lengths = []
        bond_angles = []
        dihedrals = []

        for residue in structure.get_residues():
            if not is_aa(residue):
                continue

            # Analyze bond lengths
            if all(atom in residue for atom in ["N", "CA", "C"]):
                n_ca = residue["N"] - residue["CA"]
                ca_c = residue["CA"] - residue["C"]
                bond_lengths.extend([abs(n_ca - IDEAL_BOND_LENGTHS["N-CA"]), abs(ca_c - IDEAL_BOND_LENGTHS["CA-C"])])

            # Analyze bond angles
            if all(atom in residue for atom in ["N", "CA", "C"]):
                angle = calc_angle(residue["N"].get_vector(), residue["CA"].get_vector(), residue["C"].get_vector())
                bond_angles.append(abs(np.degrees(angle) - IDEAL_BOND_ANGLES["N-CA-C"]))

            # Analyze backbone dihedrals
            phi, psi = calculate_phi_psi(residue)
            if phi is not None:
                dihedrals.append(phi)
            if psi is not None:
                dihedrals.append(psi)

        return {
            "mean_bond_dev": float(np.mean(bond_lengths)),
            "std_bond_dev": float(np.std(bond_lengths)),
            "max_bond_dev": float(np.max(bond_lengths)),
            "mean_angle_dev": float(np.mean(bond_angles)),
            "std_angle_dev": float(np.std(bond_angles)),
            "max_angle_dev": float(np.max(bond_angles)),
            "mean_dihedral": float(np.mean(dihedrals)) if dihedrals else 0.0,
            "std_dihedral": float(np.std(dihedrals)) if dihedrals else 0.0,
        }

    except Exception as e:
        logger.error(f"Error analyzing geometry: {str(e)}")
        return {}


def _analyze_ramachandran(structure: Structure) -> Dict[str, Any]:
    """Analyze Ramachandran plot statistics.

    Args:
        structure: BioPython Structure object

    Returns:
        Dictionary containing Ramachandran analysis results
    """
    try:
        phi_psi = []
        outliers = []

        for residue in structure.get_residues():
            if not is_aa(residue):
                continue

            phi, psi = calculate_phi_psi(residue)
            if phi is None or psi is None:
                continue

            phi_psi.append((phi, psi))

            # Check if in allowed regions
            in_allowed = False
            for region, (phi_range, psi_range) in RAMACHANDRAN_REGIONS.items():
                if phi_range[0] <= phi <= phi_range[1] and psi_range[0] <= psi <= psi_range[1]:
                    in_allowed = True
                    break

            if not in_allowed:
                outliers.append({"residue": f"{residue.get_resname()}{residue.get_id()[1]}", "phi": phi, "psi": psi})

        return {"n_residues": len(phi_psi), "n_outliers": len(outliers), "outliers": outliers, "percent_favored": 100 * (1 - len(outliers) / len(phi_psi)) if phi_psi else 0.0}

    except Exception as e:
        logger.error(f"Error analyzing Ramachandran plot: {str(e)}")
        return {}


def _analyze_packing(structure: Structure) -> Dict[str, Any]:
    """Analyze packing quality of structure.

    Args:
        structure: BioPython Structure object

    Returns:
        Dictionary containing packing analysis results
    """
    try:
        # Calculate atomic packing densities
        densities = []
        outliers = []
        residue_densities = {}

        atoms = list(structure.get_atoms())
        coords = np.array([atom.get_coord() for atom in atoms])

        for i, atom in enumerate(atoms):
            # Calculate local density in 6A sphere
            dists = np.linalg.norm(coords - coords[i], axis=1)
            local_atoms = np.sum((dists > 0) & (dists < 6.0))
            densities.append(local_atoms)

            # Track density by residue
            res_id = atom.get_parent().get_id()[1]
            if res_id not in residue_densities:
                residue_densities[res_id] = []
            residue_densities[res_id].append(local_atoms)

            if local_atoms < 8:  # Arbitrary threshold for underpacking
                outliers.append({"atom": f"{atom.get_parent().get_resname()}{atom.get_parent().get_id()[1]}.{atom.get_name()}", "density": int(local_atoms)})

        # Calculate per-residue average densities
        avg_residue_densities = {res_id: float(np.mean(densities)) for res_id, densities in residue_densities.items()}

        return {
            "mean_density": float(np.mean(densities)),
            "std_density": float(np.std(densities)),
            "min_density": float(min(densities)),
            "max_density": float(max(densities)),
            "n_outliers": len(outliers),
            "outliers": outliers,
            "residue_densities": avg_residue_densities,
        }

    except Exception as e:
        logger.error(f"Error analyzing packing: {str(e)}")
        return {}


def _calculate_overall_score(geometry: Dict[str, Any], rama: Dict[str, Any], packing: Dict[str, Any]) -> float:
    """Calculate overall quality score from individual metrics.

    Args:
        geometry: Geometry analysis results
        rama: Ramachandran analysis results
        packing: Packing analysis results

    Returns:
        Overall quality score between 0 and 1
    """
    try:
        # Geometry score components
        bond_dev_score = max(0, 1 - (0.5 * geometry.get("mean_bond_dev", 0)))
        angle_dev_score = max(0, 1 - (0.3 * geometry.get("mean_angle_dev", 0) / 10))
        dihedral_score = max(0, 1 - (0.2 * abs(geometry.get("std_dihedral", 0)) / 30))  # Penalize high dihedral variability
        geom_score = (bond_dev_score + angle_dev_score + dihedral_score) / 3

        # Ramachandran score
        rama_score = rama.get("percent_favored", 0) / 100

        # Packing score components
        density_score = max(0, min(1, packing.get("mean_density", 0) / 12))  # Normalize by expected density
        uniformity_score = max(0, 1 - packing.get("std_density", 0) / 4)  # Penalize high density variability
        pack_score = (density_score + uniformity_score) / 2

        # Combine scores with weights
        # Higher weight on geometry and Ramachandran since they're more fundamental
        return 0.4 * geom_score + 0.4 * rama_score + 0.2 * pack_score

    except Exception as e:
        logger.error(f"Error calculating overall score: {str(e)}")
        return 0.0


def _are_bonded(atom1: Atom, atom2: Atom) -> bool:
    """Check if two atoms are covalently bonded.

    Args:
        atom1: First BioPython Atom object
        atom2: Second BioPython Atom object

    Returns:
        True if atoms are likely bonded
    """
    # Simple distance-based check
    cutoff = 2.0  # Maximum covalent bond length
    return atom1 - atom2 < cutoff


def _calculate_alphafold_metrics(structure: Structure) -> Dict[str, Any]:
    """Calculate AlphaFold-specific quality metrics.

    Args:
        structure: BioPython Structure object

    Returns:
        Dictionary containing AlphaFold quality metrics
    """
    try:
        # Get confidence scores from AlphaFold
        confidence = alphafold.get_confidence_scores(structure)

        # Extract pLDDT scores from B-factors
        plddt_scores = [atom.get_bfactor() for atom in structure.get_atoms()]

        metrics = {
            # Overall confidence metrics
            "mean_plddt": float(np.mean(plddt_scores)),
            "min_plddt": float(min(plddt_scores)),
            "max_plddt": float(max(plddt_scores)),
            "ptm_score": confidence.get("ptm", 0.0),
            "iptm_score": confidence.get("iptm", 0.0),
            # Per-residue pLDDT scores from B-factors
            "plddt_per_residue": {int(res.get_id()[1]): float(np.mean([atom.get_bfactor() for atom in res])) for res in structure.get_residues() if is_aa(res)},
            # PAE (Predicted Aligned Error) matrix if available
            "pae_matrix": confidence.get("pae"),
            # Additional raw scores from prediction
            "raw_scores": confidence.get("raw_scores", {}),
        }

        # Classify regions by pLDDT ranges
        regions = {"very_high": [], "high": [], "medium": [], "low": []}
        for res_id, score in metrics["plddt_per_residue"].items():
            if score > 90:
                regions["very_high"].append(res_id)
            elif score > 70:
                regions["high"].append(res_id)
            elif score > 50:
                regions["medium"].append(res_id)
            else:
                regions["low"].append(res_id)
        metrics["confidence_regions"] = regions

        # Calculate PAE-based metrics if available
        if metrics["pae_matrix"] is not None:
            pae = np.array(metrics["pae_matrix"])
            metrics.update(
                {
                    "mean_pae": float(np.mean(pae)),
                    "max_pae": float(np.max(pae)),
                    "pae_long_range": float(np.mean(np.triu(pae, k=12))),  # Average PAE for residues >12 apart
                }
            )

        return metrics

    except Exception as e:
        logger.error(f"Error calculating AlphaFold metrics: {str(e)}")
        return {}


def get_vdw_radius(atom: Atom) -> float:
    """Get van der Waals radius for an atom.

    Args:
        atom: BioPython Atom object

    Returns:
        Van der Waals radius in Angstroms
    """
    # Standard vdW radii
    radii = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "P": 1.8, "H": 1.2}
    return radii.get(atom.element, 1.5)  # Default radius


def is_aa(residue: Residue) -> bool:
    """Check if residue is a standard amino acid.

    Args:
        residue: BioPython Residue object

    Returns:
        True if standard amino acid
    """
    standard_aas = {"ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"}
    return residue.get_resname() in standard_aas


def calculate_phi_psi(residue: Residue) -> Tuple[Optional[float], Optional[float]]:
    """Calculate phi/psi angles for a residue.

    Args:
        residue: BioPython Residue object

    Returns:
        Tuple of (phi, psi) angles in degrees, or (None, None) if can't calculate
    """
    try:
        phi = psi = None

        # Need -C, N, CA, C atoms for phi
        prev = residue.get_parent()[residue.get_id()[1] - 1]
        if prev and "C" in prev and all(a in residue for a in ["N", "CA", "C"]):
            phi = calc_dihedral(prev["C"].get_vector(), residue["N"].get_vector(), residue["CA"].get_vector(), residue["C"].get_vector())
            phi = np.degrees(phi)

        # Need N, CA, C, +N atoms for psi
        next_res = residue.get_parent()[residue.get_id()[1] + 1]
        if next_res and "N" in next_res and all(a in residue for a in ["N", "CA", "C"]):
            psi = calc_dihedral(residue["N"].get_vector(), residue["CA"].get_vector(), residue["C"].get_vector(), next_res["N"].get_vector())
            psi = np.degrees(psi)

        return phi, psi

    except Exception as e:
        logger.error(f"Error calculating phi/psi angles: {str(e)}")
        return None, None


def calculate_chi_angles(residue: Residue) -> List[float]:
    """Calculate side chain chi angles for a residue.

    Args:
        residue: BioPython Residue object

    Returns:
        List of chi angles in degrees
    """
    try:
        chi = []

        # Chi1: N-CA-CB-CG
        if all(a in residue for a in ["N", "CA", "CB", "CG"]):
            angle = calc_dihedral(residue["N"].get_vector(), residue["CA"].get_vector(), residue["CB"].get_vector(), residue["CG"].get_vector())
            chi.append(np.degrees(angle))

        # Add other chi angles as needed...

        return chi

    except Exception as e:
        logger.error(f"Error calculating chi angles: {str(e)}")
        return []


def evaluate_rotamer(resname: str, chi_angles: List[float]) -> float:
    """Evaluate how well chi angles match expected rotamers.

    Args:
        resname: Residue name
        chi_angles: List of chi angles in degrees

    Returns:
        Score between 0 and 1 indicating rotamer quality
    """
    try:
        # Simplified scoring - just check if angles are near multiples of 60°
        scores = []
        for chi in chi_angles:
            # Calculate minimum deviation from 60° rotamer positions
            devs = [abs(chi - x) % 360 for x in range(0, 360, 60)]
            min_dev = min(devs)
            # Convert to score between 0 and 1
            score = max(0, 1 - min_dev / 30)  # Linear falloff up to 30° deviation
            scores.append(score)

        return np.mean(scores) if scores else 0.0

    except Exception as e:
        logger.error(f"Error evaluating rotamer: {str(e)}")
        return 0.0
