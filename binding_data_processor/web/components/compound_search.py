"""Web component for compound searching.

This module provides components for:
1. Text-based compound search
2. Structure-based search
3. Property-based filtering
4. Advanced search options
"""

import logging
from typing import List, Optional, Dict, Any
from dataclasses import dataclass

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from ...models.compound.enhanced import EnhancedCompound


@dataclass
class CompoundSearchConfig:
    """Configuration for compound search component."""
    
    # Search settings
    enable_text_search: bool = True
    enable_structure_search: bool = True
    enable_property_search: bool = True
    
    # Text search settings
    search_fields: List[str] = None
    min_query_length: int = 2
    max_results: int = 100
    
    # Structure search settings
    similarity_threshold: float = 0.7
    max_structure_results: int = 50
    
    # Property search settings
    property_ranges: Dict[str, Dict[str, float]] = None
    
    def __post_init__(self):
        """Initialize configuration."""
        if self.search_fields is None:
            self.search_fields = ["name", "smiles", "cas_number"]
            
        if self.property_ranges is None:
            self.property_ranges = {
                "molecular_weight": {"min": 0, "max": 1000},
                "logp": {"min": -5, "max": 10},
                "hbd": {"min": 0, "max": 10},
                "hba": {"min": 0, "max": 20},
                "tpsa": {"min": 0, "max": 200},
            }


class CompoundSearch:
    """Component for compound searching."""

    def __init__(
        self,
        config: Optional[CompoundSearchConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize compound search component.
        
        Args:
            config: Optional component configuration
            logger: Optional logger instance
        """
        self.config = config or CompoundSearchConfig()
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        
        # Initialize state
        self.compounds = []
        self.search_results = []
        self.current_query = None
        self.current_filters = {}

    def update_compounds(self, compounds: List[EnhancedCompound]) -> None:
        """Update compound list.
        
        Args:
            compounds: List of compounds to search
        """
        self.compounds = compounds
        self._clear_results()

    def text_search(
        self,
        query: str,
        fields: Optional[List[str]] = None,
    ) -> List[EnhancedCompound]:
        """Search compounds by text.
        
        Args:
            query: Search query
            fields: Optional list of fields to search
            
        Returns:
            List of matching compounds
        """
        if not self.config.enable_text_search:
            return []
            
        if len(query) < self.config.min_query_length:
            return []
            
        try:
            # Normalize query
            query = query.lower().strip()
            
            # Get search fields
            search_fields = fields or self.config.search_fields
            
            # Search compounds
            results = []
            for compound in self.compounds:
                for field in search_fields:
                    if hasattr(compound, field):
                        value = str(getattr(compound, field)).lower()
                        if query in value:
                            results.append(compound)
                            break
                            
                if len(results) >= self.config.max_results:
                    break
            
            self.search_results = results
            self.current_query = query
            return results
            
        except Exception as e:
            self.logger.error(f"Error in text search: {str(e)}")
            return []

    def structure_search(
        self,
        smiles: str,
        threshold: Optional[float] = None,
    ) -> List[EnhancedCompound]:
        """Search compounds by structure similarity.
        
        Args:
            smiles: SMILES string to search
            threshold: Optional similarity threshold
            
        Returns:
            List of similar compounds
        """
        if not self.config.enable_structure_search:
            return []
            
        try:
            # Create query molecule
            query_mol = Chem.MolFromSmiles(smiles)
            if not query_mol:
                return []
            
            # Generate fingerprint
            query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2)
            
            # Get threshold
            sim_threshold = threshold or self.config.similarity_threshold
            
            # Search compounds
            results = []
            for compound in self.compounds:
                try:
                    # Create target molecule
                    target_mol = Chem.MolFromSmiles(compound.smiles)
                    if not target_mol:
                        continue
                    
                    # Generate fingerprint
                    target_fp = AllChem.GetMorganFingerprintAsBitVect(target_mol, 2)
                    
                    # Calculate similarity
                    similarity = DataStructs.TanimotoSimilarity(query_fp, target_fp)
                    
                    # Add if similar enough
                    if similarity >= sim_threshold:
                        compound.similarity = similarity
                        results.append(compound)
                        
                except Exception as e:
                    self.logger.error(
                        f"Error processing compound {compound.name}: {str(e)}"
                    )
                    continue
                    
                if len(results) >= self.config.max_structure_results:
                    break
            
            # Sort by similarity
            results.sort(key=lambda x: x.similarity, reverse=True)
            
            self.search_results = results
            self.current_query = smiles
            return results
            
        except Exception as e:
            self.logger.error(f"Error in structure search: {str(e)}")
            return []

    def property_search(
        self,
        filters: Dict[str, Dict[str, float]],
    ) -> List[EnhancedCompound]:
        """Search compounds by property filters.
        
        Args:
            filters: Dictionary of property filters
            
        Returns:
            List of matching compounds
        """
        if not self.config.enable_property_search:
            return []
            
        try:
            # Validate filters
            valid_filters = {}
            for prop, ranges in filters.items():
                if prop not in self.config.property_ranges:
                    continue
                    
                valid_ranges = {}
                if "min" in ranges:
                    valid_ranges["min"] = max(
                        ranges["min"],
                        self.config.property_ranges[prop]["min"]
                    )
                if "max" in ranges:
                    valid_ranges["max"] = min(
                        ranges["max"],
                        self.config.property_ranges[prop]["max"]
                    )
                    
                if valid_ranges:
                    valid_filters[prop] = valid_ranges
            
            # Search compounds
            results = []
            for compound in self.compounds:
                matches = True
                for prop, ranges in valid_filters.items():
                    if not hasattr(compound, prop):
                        matches = False
                        break
                        
                    value = getattr(compound, prop)
                    if value is None:
                        matches = False
                        break
                        
                    if "min" in ranges and value < ranges["min"]:
                        matches = False
                        break
                        
                    if "max" in ranges and value > ranges["max"]:
                        matches = False
                        break
                        
                if matches:
                    results.append(compound)
            
            self.search_results = results
            self.current_filters = valid_filters
            return results
            
        except Exception as e:
            self.logger.error(f"Error in property search: {str(e)}")
            return []

    def combined_search(
        self,
        text_query: Optional[str] = None,
        structure_query: Optional[str] = None,
        property_filters: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> List[EnhancedCompound]:
        """Perform combined search using multiple criteria.
        
        Args:
            text_query: Optional text search query
            structure_query: Optional structure search query
            property_filters: Optional property filters
            
        Returns:
            List of matching compounds
        """
        try:
            results = set(self.compounds)
            
            # Apply text search
            if text_query:
                text_results = set(self.text_search(text_query))
                results &= text_results
            
            # Apply structure search
            if structure_query:
                structure_results = set(self.structure_search(structure_query))
                results &= structure_results
            
            # Apply property filters
            if property_filters:
                property_results = set(self.property_search(property_filters))
                results &= property_results
            
            # Convert to list and sort
            results = list(results)
            results.sort(key=lambda x: x.name)
            
            self.search_results = results
            return results
            
        except Exception as e:
            self.logger.error(f"Error in combined search: {str(e)}")
            return []

    def get_search_stats(self) -> Dict[str, Any]:
        """Get search statistics.
        
        Returns:
            Dictionary of search statistics
        """
        stats = {
            "total_compounds": len(self.compounds),
            "results": len(self.search_results),
        }
        
        if self.current_query:
            stats["query"] = self.current_query
            
        if self.current_filters:
            stats["filters"] = self.current_filters
            
        return stats

    def _clear_results(self) -> None:
        """Clear search results."""
        self.search_results = []
        self.current_query = None
        self.current_filters = {}
