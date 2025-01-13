"""Types module for structure analysis system.

This module provides:
1. Custom types
2. Type definitions
3. Standard interfaces
4. Type aliases
"""

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import (
    Dict,
    List,
    Set,
    FrozenSet,
    Tuple,
    Optional,
    Union,
    Any,
    TypeVar,
    Generic,
    Protocol,
    runtime_checkable,
    Callable,
    Iterator,
    Mapping,
    Sequence,
    NewType,
)
from pathlib import Path
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom

# Type variables
T = TypeVar("T")
U = TypeVar("U")

# Custom types
AtomID = NewType("AtomID", str)
ResidueID = NewType("ResidueID", str)
ChainID = NewType("ChainID", str)
ModelID = NewType("ModelID", int)
StructureID = NewType("StructureID", str)

# Coordinate types
Coordinate = Tuple[float, float, float]
Vector = np.ndarray  # Shape (3,)
Matrix = np.ndarray  # Shape (N, 3)
Transform = np.ndarray  # Shape (4, 4)

# Structure types
AtomSelection = Set[AtomID]
ResidueSelection = Set[ResidueID]
ChainSelection = Set[ChainID]
StructureSelection = Dict[ChainID, ResidueSelection]

# Analysis types
Distance = float
Angle = float
Energy = float
Score = float
Probability = float

# Result types
ValidationResult = Dict[str, Union[bool, List[str]]]
AnalysisResult = Dict[str, Any]
PredictionResult = Dict[str, Union[float, List[float]]]


@dataclass
class Point3D:
    """3D point with coordinates."""

    x: float
    y: float
    z: float

    def to_array(self) -> np.ndarray:
        """Convert to numpy array."""
        return np.array([self.x, self.y, self.z])

    def distance_to(self, other: "Point3D") -> float:
        """Calculate distance to another point."""
        return np.linalg.norm(self.to_array() - other.to_array())


@dataclass
class BoundingBox:
    """Axis-aligned bounding box."""

    min_point: Point3D
    max_point: Point3D

    @property
    def center(self) -> Point3D:
        """Get center point."""
        return Point3D(
            x=(self.min_point.x + self.max_point.x) / 2,
            y=(self.min_point.y + self.max_point.y) / 2,
            z=(self.min_point.z + self.max_point.z) / 2,
        )

    @property
    def dimensions(self) -> Point3D:
        """Get box dimensions."""
        return Point3D(
            x=self.max_point.x - self.min_point.x,
            y=self.max_point.y - self.min_point.y,
            z=self.max_point.z - self.min_point.z,
        )


@dataclass
class AtomInfo:
    """Atom information."""

    id: AtomID
    element: str
    name: str
    coordinates: Point3D
    occupancy: float = 1.0
    b_factor: float = 0.0
    charge: float = 0.0
    radius: float = 0.0


@dataclass
class ResidueInfo:
    """Residue information."""

    id: ResidueID
    name: str
    number: int
    atoms: Dict[AtomID, AtomInfo] = field(default_factory=dict)
    center: Optional[Point3D] = None


@dataclass
class ChainInfo:
    """Chain information."""

    id: ChainID
    residues: Dict[ResidueID, ResidueInfo] = field(default_factory=dict)
    center: Optional[Point3D] = None


@dataclass
class StructureInfo:
    """Structure information."""

    id: StructureID
    chains: Dict[ChainID, ChainInfo] = field(default_factory=dict)
    center: Optional[Point3D] = None
    box: Optional[BoundingBox] = None


@dataclass
class InteractionInfo:
    """Interaction information."""

    type: str
    atom1: AtomInfo
    atom2: AtomInfo
    distance: float
    angle: Optional[float] = None
    energy: Optional[float] = None


@dataclass
class ContactInfo:
    """Contact information."""

    residue1: ResidueInfo
    residue2: ResidueInfo
    interactions: List[InteractionInfo] = field(default_factory=list)
    score: float = 0.0


@runtime_checkable
class Transformer(Protocol[T, U]):
    """Protocol for structure transformers."""

    def transform(self, data: T) -> U:
        """Transform input data."""
        ...

    def inverse_transform(self, data: U) -> T:
        """Inverse transform data."""
        ...


@runtime_checkable
class Predictor(Protocol[T]):
    """Protocol for structure predictors."""

    def predict(self, structure: Structure) -> T:
        """Make prediction for structure."""
        ...

    def predict_proba(self, structure: Structure) -> Dict[str, float]:
        """Get prediction probabilities."""
        ...


@runtime_checkable
class Analyzer(Protocol):
    """Protocol for structure analyzers."""

    def analyze(self, structure: Structure) -> AnalysisResult:
        """Analyze structure."""
        ...


@runtime_checkable
class Validator(Protocol):
    """Protocol for structure validators."""

    def validate(self, structure: Structure) -> ValidationResult:
        """Validate structure."""
        ...


# Type aliases for common operations
DistanceFunction = Callable[[Point3D, Point3D], float]
AngleFunction = Callable[[Point3D, Point3D, Point3D], float]
EnergyFunction = Callable[[Structure], float]
ScoreFunction = Callable[[Structure], float]
TransformFunction = Callable[[Structure], Structure]

# Type aliases for collections
AtomDict = Dict[AtomID, AtomInfo]
ResidueDict = Dict[ResidueID, ResidueInfo]
ChainDict = Dict[ChainID, ChainInfo]
StructureDict = Dict[StructureID, StructureInfo]

InteractionList = List[InteractionInfo]
ContactList = List[ContactInfo]

# Type aliases for results
ValidationDict = Dict[str, ValidationResult]
AnalysisDict = Dict[str, AnalysisResult]
PredictionDict = Dict[str, PredictionResult]

# Type aliases for configurations
ConfigValue = Union[str, int, float, bool, None]
ConfigDict = Dict[str, ConfigValue]

# Type aliases for data
DataValue = Union[str, int, float, bool, None, List[Any], Dict[str, Any]]
DataDict = Dict[str, DataValue]

# Type aliases for paths
PathLike = Union[str, Path]
