"""Examples demonstrating usage of structure depiction classes.

This script shows how to:
1. Create and configure depiction objects
2. Generate 2D and 3D structure visualizations
3. Create grid layouts and highlight substructures
4. Save images in different formats
5. Customize drawing options
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from rdkit import Chem

from ..depiction_2d import Structure2DDepiction
from ..depiction_3d import Structure3DDepiction


def create_example_molecules() -> List[Chem.Mol]:
    """Create example molecules for demonstration."""
    smiles_list = [
        "CC1=C(C(=O)NC(=N1)C2=CC=CC=C2)C3=CC=C(C=C3)S(=O)(=O)N",  # Sildenafil
        "CC(C)(C)NC[C@H](O)C1=CC(O)=CC(O)=C1",  # Salbutamol
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
        "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",  # Testosterone
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
    ]
    return [Chem.MolFromSmiles(s) for s in smiles_list]


def example_2d_depiction():
    """Demonstrate 2D structure depiction."""
    # Create depiction object with custom options
    depict_2d = Structure2DDepiction(
        {
            "size": (400, 400),
            "show_atom_numbers": True,
            "annotate_stereo": True,
        }
    )

    # Get example molecules
    mols = create_example_molecules()
    if not mols:
        print("Error creating example molecules")
        return

    # Create output directory
    output_dir = Path("examples/output")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Basic 2D depiction
    img = depict_2d.depict_2d(mols[0], return_pil=True)
    if img:
        img.save(output_dir / "basic_2d.png")

    # Grid depiction with legends
    legends = ["Sildenafil", "Salbutamol", "Caffeine", "Testosterone", "Aspirin"]
    grid = depict_2d.depict_grid(mols, legends=legends, return_pil=True)
    if grid:
        grid.save(output_dir / "molecule_grid.png")

    # Highlight substructure
    query = Chem.MolFromSmarts("[OH]c1ccccc1")  # Phenol pattern
    if query:
        highlighted = depict_2d.depict_with_highlights(mols[1], query, return_pil=True)
        if highlighted:
            highlighted.save(output_dir / "highlighted_2d.png")


def example_3d_depiction():
    """Demonstrate 3D structure depiction."""
    # Create depiction object with custom options
    depict_3d = Structure3DDepiction(
        {
            "size": (400, 400),
            "add_hydrogens": True,
            "optimize_geometry": True,
        }
    )

    # Get example molecules
    mols = create_example_molecules()
    if not mols:
        print("Error creating example molecules")
        return

    # Create output directory
    output_dir = Path("examples/output")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Basic 3D depiction
    mol_3d = depict_3d.generate_3d_conformer(mols[0])
    if mol_3d:
        img = depict_3d.depict_3d(mol_3d, return_pil=True)
        if img:
            img.save(output_dir / "basic_3d.png")

    # Conformer grid
    mol_3d = depict_3d.generate_conformers(mols[0], n_conf=6)
    if mol_3d:
        grid = depict_3d.depict_conformer_grid(mol_3d, mols_per_row=3, return_pil=True)
        if grid:
            grid.save(output_dir / "conformer_grid.png")

    # Pharmacophore features
    mol_3d = depict_3d.generate_3d_conformer(mols[1])
    if mol_3d:
        # Define pharmacophore features
        features = {
            "donor": [4, 6, 8],  # OH groups
            "acceptor": [3],  # Ether oxygen
            "positive": [1],  # Amine
        }
        img = depict_3d.depict_pharmacophore(mol_3d, features, return_pil=True)
        if img:
            img.save(output_dir / "pharmacophore_3d.png")


def example_advanced_depiction():
    """Demonstrate advanced depiction features."""
    # Create depiction objects
    depict_2d = Structure2DDepiction()
    depict_3d = Structure3DDepiction()

    # Get example molecules
    mols = create_example_molecules()
    if not mols:
        print("Error creating example molecules")
        return

    # Create output directory
    output_dir = Path("examples/output")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Custom atom colors and labels
    mol = mols[2]  # Caffeine
    atom_colors = {
        0: (1, 0, 0),  # Red
        1: (0, 1, 0),  # Green
        2: (0, 0, 1),  # Blue
    }
    atom_labels = {
        0: "Me",
        1: "N",
        2: "C",
    }
    img = depict_2d.depict_2d(
        mol,
        atom_colors=atom_colors,
        atom_labels=atom_labels,
        return_pil=True,
    )
    if img:
        img.save(output_dir / "custom_labels.png")

    # Reaction mechanism
    rxn = Chem.ReactionFromSmarts(
        "[OH:1][C:2]=[O:3].[NH2:4][C:5]>>[O:1]=[C:2][NH:4][C:5].[OH2:3]"
    )
    if rxn:
        steps = [
            {
                "reactants": [mols[4]],  # Aspirin
                "products": [mols[3]],  # Testosterone
                "annotation": "Amide formation",
            }
        ]
        img = depict_3d.depict_reaction_mechanism(rxn, steps, return_pil=True)
        if img:
            img.save(output_dir / "reaction.png")


def main():
    """Run all examples."""
    print("Running 2D depiction examples...")
    example_2d_depiction()

    print("Running 3D depiction examples...")
    example_3d_depiction()

    print("Running advanced depiction examples...")
    example_advanced_depiction()

    print("Done! Check examples/output directory for results.")


if __name__ == "__main__":
    main()
