"""Enhanced compound search component.

This module provides an enhanced web interface for searching compounds with:
- Text search
- Structure search
- Advanced filtering
- Search history
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from datetime import datetime

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from flask import render_template, request, jsonify

from ...models.compound import Compound
from ..base import BaseComponent, ViewResult


class CompoundSearchEnhanced(BaseComponent):
    """Enhanced compound search component."""

    def __init__(
        self,
        template_dir: Optional[Path] = None,
        static_dir: Optional[Path] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize search component.
        
        Args:
            template_dir: Optional template directory
            static_dir: Optional static files directory
            logger: Optional logger instance
        """
        super().__init__(template_dir, static_dir, logger)

        # Initialize tracking
        self.search_stats = {
            "total_searches": 0,
            "text_searches": 0,
            "structure_searches": 0,
            "filter_searches": 0,
            "search_history": [],
        }

    def render_search(
        self,
        compounds: List[Compound],
        query: Optional[str] = None,
        structure: Optional[str] = None,
        filters: Optional[Dict[str, Any]] = None,
    ) -> ViewResult:
        """Render search interface.
        
        Args:
            compounds: List of compounds to search
            query: Optional text query
            structure: Optional structure query (SMILES)
            filters: Optional filters to apply
            
        Returns:
            ViewResult containing rendered HTML and metadata
        """
        try:
            # Perform search
            results = []
            if query:
                results = self._text_search(compounds, query)
                self.search_stats["text_searches"] += 1
            elif structure:
                results = self._structure_search(compounds, structure)
                self.search_stats["structure_searches"] += 1
            else:
                results = compounds

            # Apply filters
            if filters:
                results = self._apply_filters(results, filters)
                self.search_stats["filter_searches"] += 1

            # Update stats
            self.search_stats["total_searches"] += 1
            self.search_stats["search_history"].append({
                "timestamp": datetime.now().isoformat(),
                "query": query,
                "structure": structure,
                "filters": filters,
                "result_count": len(results),
            })

            # Render template
            html = render_template(
                "compound_search.html",
                compounds=results,
                query=query,
                structure=structure,
                filters=filters or {},
                stats=self.search_stats,
            )

            return ViewResult(
                success=True,
                data={
                    "html": html,
                    "compounds": results,
                    "search_stats": self.search_stats,
                },
            )

        except Exception as e:
            self.logger.error(f"Error rendering search: {str(e)}")
            return ViewResult(
                success=False,
                error=str(e),
            )

    def _text_search(
        self,
        compounds: List[Compound],
        query: str,
    ) -> List[Compound]:
        """Perform text search on compounds.
        
        Args:
            compounds: List of compounds to search
            query: Search query
            
        Returns:
            List of matching compounds
        """
        query = query.lower()
        results = []

        for compound in compounds:
            # Search name
            if query in compound.name.lower():
                results.append(compound)
                continue

            # Search CAS number
            if query in compound.cas_number:
                results.append(compound)
                continue

            # Search binding data
            if hasattr(compound, "binding_data"):
                for binding in compound.binding_data:
                    if query in binding["target"].lower():
                        results.append(compound)
                        break

            # Search social data
            if hasattr(compound, "social_data"):
                data = compound.social_data
                if "reddit" in data:
                    for post in data["reddit"].get("posts", []):
                        if query in post["title"].lower():
                            results.append(compound)
                            break
                if "twitter" in data:
                    for tweet in data["twitter"].get("tweets", []):
                        if query in tweet["text"].lower():
                            results.append(compound)
                            break

        return results

    def _structure_search(
        self,
        compounds: List[Compound],
        structure: str,
        threshold: float = 0.7,
    ) -> List[Compound]:
        """Perform structure similarity search.
        
        Args:
            compounds: List of compounds to search
            structure: SMILES query structure
            threshold: Similarity threshold (0-1)
            
        Returns:
            List of similar compounds
        """
        # Parse query structure
        query_mol = Chem.MolFromSmiles(structure)
        if not query_mol:
            raise ValueError(f"Invalid SMILES: {structure}")

        # Generate query fingerprint
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2)

        # Search compounds
        results = []
        for compound in compounds:
            mol = Chem.MolFromSmiles(compound.smiles)
            if not mol:
                continue

            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
            similarity = DataStructs.TanimotoSimilarity(query_fp, fp)

            if similarity >= threshold:
                compound.similarity = similarity
                results.append(compound)

        # Sort by similarity
        results.sort(key=lambda x: x.similarity, reverse=True)
        return results

    def _apply_filters(
        self,
        compounds: List[Compound],
        filters: Dict[str, Any],
    ) -> List[Compound]:
        """Apply filters to compound list.
        
        Args:
            compounds: List of compounds to filter
            filters: Filter parameters
            
        Returns:
            Filtered list of compounds
        """
        filtered = compounds

        # Filter by binding affinity
        if "min_binding" in filters:
            filtered = [
                c for c in filtered
                if hasattr(c, "binding_data")
                and any(
                    float(b["affinity"]) >= filters["min_binding"]
                    for b in c.binding_data
                )
            ]

        # Filter by target
        if "target" in filters:
            filtered = [
                c for c in filtered
                if hasattr(c, "binding_data")
                and any(
                    filters["target"].lower() in b["target"].lower()
                    for b in c.binding_data
                )
            ]

        # Filter by social data
        if "has_social_data" in filters:
            filtered = [
                c for c in filtered
                if hasattr(c, "social_data") == filters["has_social_data"]
            ]

        # Filter by date range
        if "date_from" in filters and "date_to" in filters:
            date_from = datetime.fromisoformat(filters["date_from"])
            date_to = datetime.fromisoformat(filters["date_to"])
            filtered = [
                c for c in filtered
                if hasattr(c, "social_data")
                and any(
                    date_from <= datetime.fromisoformat(post["created_utc"])
                    <= date_to
                    for post in c.social_data.get("reddit", {}).get("posts", [])
                )
            ]

        return filtered

    def get_metrics(self) -> Dict[str, Any]:
        """Get component metrics."""
        return {
            "search_stats": self.search_stats,
        }
