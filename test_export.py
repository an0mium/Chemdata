#!/usr/bin/env python3

from pathlib import Path
import json
import pandas as pd


def main():
    # Create output directory
    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)

    # Create a test compound
    compound = {
        "name": "Test Compound",
        "cas_number": "123-45-6",
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "compound_type": "OTHER",
        "molecular_weight": 180.16,
        "logp": 1.43,
        "tpsa": 63.6,
        "hbd": 2,
        "hba": 4,
        "rotatable_bonds": 3,
    }

    # Export to different formats
    # TSV
    df = pd.DataFrame([compound])
    df.to_csv(output_dir / "compounds.tsv", sep="\t", index=False)

    # JSON
    with open(output_dir / "compounds.json", "w") as f:
        json.dump([compound], f, indent=2)

    # Excel
    df.to_excel(output_dir / "compounds.xlsx", index=False)


if __name__ == "__main__":
    main()
