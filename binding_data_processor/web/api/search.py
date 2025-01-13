"""FastAPI endpoints for compound search."""

from typing import Dict, List, Optional, Any
from fastapi import APIRouter, Depends, HTTPException, Query
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

from ...models.compound import CompoundData
from ...models.validation import ValidationError
from ..cache import Cache
from ..rate_limiter import RateLimiter
from .compounds import CompoundResponse, CompoundListResponse


class SearchQuery(BaseModel):
    """Search query model."""

    text: Optional[str] = Field(None, description="Text search query")
    structure: Optional[str] = Field(None, description="Structure search (SMILES/InChI)")
    similarity: Optional[float] = Field(None, ge=0, le=1, description="Similarity threshold")
    property_ranges: Optional[Dict[str, Dict[str, float]]] = Field(
        None,
        description="Property range filters (e.g. {'logP': {'min': 0, 'max': 5}})",
    )
    binding_ranges: Optional[Dict[str, Dict[str, float]]] = Field(
        None,
        description="Binding data range filters (e.g. {'5HT2A': {'min': 7.0}})",
    )
    safety_ranges: Optional[Dict[str, Dict[str, float]]] = Field(
        None,
        description="Safety data range filters",
    )
    quantum_ranges: Optional[Dict[str, Dict[str, float]]] = Field(
        None,
        description="Quantum property range filters",
    )
    tags: Optional[List[str]] = Field(None, description="Tag filters")
    sort_by: Optional[str] = Field(None, description="Sort field")
    sort_order: Optional[str] = Field("asc", description="Sort order (asc/desc)")


class SearchSuggestion(BaseModel):
    """Search suggestion model."""

    text: str = Field(..., description="Suggestion text")
    type: str = Field(..., description="Suggestion type")
    score: float = Field(..., description="Relevance score")
    metadata: Optional[Dict[str, Any]] = Field(None, description="Additional metadata")


class SearchSuggestionsResponse(BaseModel):
    """Response model for search suggestions."""

    suggestions: List[SearchSuggestion]
    total: int = Field(..., description="Total number of suggestions")


class SearchAPI:
    """API endpoints for compound search."""

    def __init__(self):
        """Initialize API."""
        self.router = APIRouter()
        self.cache = Cache()
        self.rate_limiter = RateLimiter()
        self._setup_routes()

    def _setup_routes(self):
        """Set up API routes."""
        self.router.add_api_route(
            "/search",
            self.search,
            methods=["POST"],
            response_model=CompoundListResponse,
            tags=["search"],
        )
        self.router.add_api_route(
            "/search/suggestions",
            self.get_suggestions,
            methods=["GET"],
            response_model=SearchSuggestionsResponse,
            tags=["search"],
        )
        self.router.add_api_route(
            "/search/structure",
            self.search_structure,
            methods=["POST"],
            response_model=CompoundListResponse,
            tags=["search"],
        )
        self.router.add_api_route(
            "/search/similarity",
            self.search_similarity,
            methods=["POST"],
            response_model=CompoundListResponse,
            tags=["search"],
        )

    async def search(
        self,
        query: SearchQuery,
        page: int = Query(1, ge=1, description="Page number"),
        per_page: int = Query(50, ge=1, le=100, description="Items per page"),
        viewport_width: Optional[int] = Query(None, description="Viewport width"),
    ) -> CompoundListResponse:
        """Search compounds with query.

        Args:
            query: Search query parameters
            page: Page number
            per_page: Items per page
            viewport_width: Optional viewport width for responsive data

        Returns:
            Paginated compound list response

        Raises:
            HTTPException: If validation fails
        """
        # Apply rate limiting
        await self.rate_limiter.acquire()

        try:
            # Get cache key
            cache_key = f"search:{hash(str(query))}"

            # Try cache first
            if results := await self.cache.get(cache_key):
                return CompoundListResponse(**results)

            # Execute search
            compounds = await self._execute_search(
                query=query,
                page=page,
                per_page=per_page,
            )

            # Format for viewport if needed
            if viewport_width:
                compounds = [self._format_for_viewport(c, viewport_width) for c in compounds]

            # Get total count
            total = await self._get_search_total(query)

            # Calculate pagination
            has_next = (page * per_page) < total
            has_prev = page > 1

            # Format response
            response = CompoundListResponse(
                compounds=[CompoundResponse(**c.dict()) for c in compounds],
                total=total,
                page=page,
                per_page=per_page,
                has_next=has_next,
                has_prev=has_prev,
            )

            # Cache results
            await self.cache.set(cache_key, response.dict())

            return response

        except ValidationError as e:
            raise HTTPException(status_code=400, detail=str(e))

        finally:
            self.rate_limiter.release()

    async def get_suggestions(
        self,
        query: str = Query(..., description="Search query text"),
        limit: int = Query(10, ge=1, le=50, description="Maximum suggestions"),
    ) -> SearchSuggestionsResponse:
        """Get search suggestions for query.

        Args:
            query: Search query text
            limit: Maximum number of suggestions

        Returns:
            Search suggestions response

        Raises:
            HTTPException: If validation fails
        """
        # Apply rate limiting
        await self.rate_limiter.acquire()

        try:
            # Get cache key
            cache_key = f"suggestions:{query}:{limit}"

            # Try cache first
            if suggestions := await self.cache.get(cache_key):
                return SearchSuggestionsResponse(**suggestions)

            # Get suggestions
            suggestions = await self._get_suggestions(query, limit)

            # Format response
            response = SearchSuggestionsResponse(
                suggestions=[SearchSuggestion(**s) for s in suggestions],
                total=len(suggestions),
            )

            # Cache results
            await self.cache.set(cache_key, response.dict())

            return response

        except ValidationError as e:
            raise HTTPException(status_code=400, detail=str(e))

        finally:
            self.rate_limiter.release()

    async def search_structure(
        self,
        structure: str = Query(..., description="Structure query (SMILES/InChI)"),
        page: int = Query(1, ge=1, description="Page number"),
        per_page: int = Query(50, ge=1, le=100, description="Items per page"),
        viewport_width: Optional[int] = Query(None, description="Viewport width"),
    ) -> CompoundListResponse:
        """Search compounds by structure.

        Args:
            structure: Structure query (SMILES/InChI)
            page: Page number
            per_page: Items per page
            viewport_width: Optional viewport width for responsive data

        Returns:
            Paginated compound list response

        Raises:
            HTTPException: If validation fails
        """
        # Apply rate limiting
        await self.rate_limiter.acquire()

        try:
            # Create structure search query
            query = SearchQuery(structure=structure)

            # Execute search
            return await self.search(
                query=query,
                page=page,
                per_page=per_page,
                viewport_width=viewport_width,
            )

        finally:
            self.rate_limiter.release()

    async def search_similarity(
        self,
        structure: str = Query(..., description="Structure query (SMILES/InChI)"),
        threshold: float = Query(0.7, ge=0, le=1, description="Similarity threshold"),
        page: int = Query(1, ge=1, description="Page number"),
        per_page: int = Query(50, ge=1, le=100, description="Items per page"),
        viewport_width: Optional[int] = Query(None, description="Viewport width"),
    ) -> CompoundListResponse:
        """Search compounds by structural similarity.

        Args:
            structure: Structure query (SMILES/InChI)
            threshold: Similarity threshold
            page: Page number
            per_page: Items per page
            viewport_width: Optional viewport width for responsive data

        Returns:
            Paginated compound list response

        Raises:
            HTTPException: If validation fails
        """
        # Apply rate limiting
        await self.rate_limiter.acquire()

        try:
            # Create similarity search query
            query = SearchQuery(structure=structure, similarity=threshold)

            # Execute search
            return await self.search(
                query=query,
                page=page,
                per_page=per_page,
                viewport_width=viewport_width,
            )

        finally:
            self.rate_limiter.release()

    async def _execute_search(
        self,
        query: SearchQuery,
        page: int,
        per_page: int,
    ) -> List[CompoundData]:
        """Execute search query.

        Args:
            query: Search query
            page: Page number
            per_page: Items per page

        Returns:
            List of compound data
        """
        # TODO: Implement search query
        return []

    async def _get_search_total(self, query: SearchQuery) -> int:
        """Get total count of search results.

        Args:
            query: Search query

        Returns:
            Total count
        """
        # TODO: Implement count query
        return 0

    async def _get_suggestions(
        self,
        query: str,
        limit: int,
    ) -> List[Dict[str, Any]]:
        """Get search suggestions.

        Args:
            query: Search query text
            limit: Maximum suggestions

        Returns:
            List of suggestions
        """
        # TODO: Implement suggestions
        return []

    def _format_for_viewport(
        self,
        compound: CompoundData,
        viewport_width: int,
    ) -> CompoundData:
        """Format compound data for viewport.

        Args:
            compound: Compound data
            viewport_width: Viewport width

        Returns:
            Formatted compound data
        """
        # Mobile viewport
        if viewport_width < 576:
            # Reduce data for mobile
            compound.spectral_data = None
            compound.crystal_data = None
            if compound.quantum_data:
                compound.quantum_data = {k: v for k, v in compound.quantum_data.items() if k in ["energy", "dipole"]}

        # Tablet viewport
        elif viewport_width < 992:
            # Moderate data for tablet
            if compound.quantum_data:
                compound.quantum_data = {k: v for k, v in compound.quantum_data.items() if k not in ["orbital_coefficients"]}

        return compound
