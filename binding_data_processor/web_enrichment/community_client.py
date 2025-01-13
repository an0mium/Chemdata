"""Client for interacting with community data sources."""

from typing import Dict, List, Optional, Set, Tuple

from .base_client import BaseWebClient as BaseClient
from .validation.schema import SchemaValidator, DataSource, ValidationLevel
from .validation.data import clean_community_data


class CommunityClient(BaseClient):
    """Client for fetching and processing community data."""

    def __init__(self, config: Optional[Dict] = None):
        """Initialize the community client.

        Args:
            config: Optional configuration dictionary
        """
        super().__init__(config or {})
        self.base_url = self.config.get("community_api_url", "")
        self.api_key = self.config.get("community_api_key", "")

    async def get_community_data(self, compound_id: str) -> Dict:
        """Get community data for a compound.

        Args:
            compound_id: Compound identifier

        Returns:
            Dictionary containing community data
        """
        endpoint = f"/compounds/{compound_id}/community"
        response = await self._make_request("GET", endpoint)

        if response:
            # Validate and clean the response data
            validator = SchemaValidator(level=ValidationLevel.NORMAL)
            result = validator.validate(response, DataSource.COMMUNITY, "forum")
            if result.is_valid:
                return clean_community_data(response)
        return {}

    async def get_discussions(self, compound_id: str) -> List[Dict]:
        """Get community discussions about a compound.

        Args:
            compound_id: Compound identifier

        Returns:
            List of discussion data dictionaries
        """
        endpoint = f"/compounds/{compound_id}/discussions"
        response = await self._make_request("GET", endpoint)

        if response and isinstance(response, list):
            validator = SchemaValidator(level=ValidationLevel.NORMAL)
            return [clean_community_data(item) for item in response if validator.validate(item, DataSource.COMMUNITY, "forum").is_valid]
        return []

    async def get_experience_reports(self, compound_id: str) -> List[Dict]:
        """Get experience reports for a compound.

        Args:
            compound_id: Compound identifier

        Returns:
            List of experience report dictionaries
        """
        endpoint = f"/compounds/{compound_id}/experiences"
        response = await self._make_request("GET", endpoint)

        if response and isinstance(response, list):
            validator = SchemaValidator(level=ValidationLevel.NORMAL)
            return [clean_community_data(item) for item in response if validator.validate(item, DataSource.COMMUNITY, "forum").is_valid]
        return []

    async def get_safety_notices(self, compound_id: str) -> List[Dict]:
        """Get community safety notices for a compound.

        Args:
            compound_id: Compound identifier

        Returns:
            List of safety notice dictionaries
        """
        endpoint = f"/compounds/{compound_id}/safety"
        response = await self._make_request("GET", endpoint)

        if response and isinstance(response, list):
            validator = SchemaValidator(level=ValidationLevel.NORMAL)
            return [clean_community_data(item) for item in response if validator.validate(item, DataSource.COMMUNITY, "forum").is_valid]
        return []

    async def search_discussions(self, query: str) -> List[Dict]:
        """Search community discussions.

        Args:
            query: Search query string

        Returns:
            List of matching discussion dictionaries
        """
        endpoint = "/discussions/search"
        response = await self._make_request("GET", endpoint, params={"q": query})

        if response and isinstance(response, list):
            validator = SchemaValidator(level=ValidationLevel.NORMAL)
            return [clean_community_data(item) for item in response if validator.validate(item, DataSource.COMMUNITY, "forum").is_valid]
        return []

    async def get_trending_compounds(self) -> List[str]:
        """Get list of trending compounds in community discussions.

        Returns:
            List of compound identifiers
        """
        endpoint = "/compounds/trending"
        response = await self._make_request("GET", endpoint)

        if response and isinstance(response, list):
            return [str(item) for item in response]
        return []

    async def get_community_stats(self, compound_id: str) -> Dict:
        """Get community engagement statistics for a compound.

        Args:
            compound_id: Compound identifier

        Returns:
            Dictionary of statistics
        """
        endpoint = f"/compounds/{compound_id}/stats"
        response = await self._make_request("GET", endpoint)

        if response and isinstance(response, dict):
            return clean_community_data(response)
        return {}
