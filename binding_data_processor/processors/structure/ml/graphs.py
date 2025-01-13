"""Graph-based feature generation for molecular structures."""

import logging
from typing import Dict, List, Optional, Union

import numpy as np
import torch
from rdkit import Chem
from torch_geometric.data import Data

from .base import MLProcessor
from .utils import mol_to_graph


class GraphFeatureGenerator(MLProcessor):
    """Generate graph-based features from molecular structures."""

    def __init__(
        self,
        use_3d: bool = False,
        use_chirality: bool = True,
        use_features: bool = True,
    ):
        """Initialize graph feature generator.

        Args:
            use_3d: Whether to include 3D coordinates
            use_chirality: Whether to include chirality information
            use_features: Whether to include atomic/bond features
        """
        super().__init__()
        self.use_3d = use_3d
        self.use_chirality = use_chirality
        self.use_features = use_features

    def generate(
        self,
        mols: List[Chem.Mol],
    ) -> List[Data]:
        """Generate graph representations of molecules.

        Args:
            mols: List of RDKit molecules

        Returns:
            List of PyTorch Geometric Data objects
        """
        try:
            graphs = []
            for mol in mols:
                try:
                    graph = mol_to_graph(
                        mol,
                        use_3d=self.use_3d,
                        use_chirality=self.use_chirality,
                        use_features=self.use_features,
                    )
                    graphs.append(graph)
                except Exception as e:
                    self.logger.debug(f"Error converting molecule to graph: {str(e)}")
                    # Create empty graph as placeholder
                    graphs.append(Data())
            return graphs

        except Exception as e:
            self.logger.error(f"Error generating graph features: {str(e)}")
            return []

    def get_info(self) -> Dict[str, Dict]:
        """Get information about graph features."""
        return {
            "graph": {
                "description": "Graph-based molecular representation",
                "use_3d": self.use_3d,
                "use_chirality": self.use_chirality,
                "use_features": self.use_features,
            }
        }
