"""Compound model for web enrichment data."""

from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional, Set


@dataclass
class CompoundMention:
    """Model representing a social media mention of a compound."""

    compound_id: str
    platform: str  # e.g. "reddit", "twitter"
    url: str
    text: str
    timestamp: datetime
    author: str
    sentiment_score: Optional[float] = None
    engagement_score: Optional[float] = None
    metadata: Dict = None

    def __post_init__(self):
        """Initialize optional fields."""
        if self.metadata is None:
            self.metadata = {}


@dataclass
class CompoundDiscussion:
    """Model representing a social media discussion about a compound."""

    compound_id: str
    platform: str
    url: str
    title: str
    content: str
    timestamp: datetime
    author: str
    replies: List["CompoundMention"] = None
    sentiment_score: Optional[float] = None
    engagement_score: Optional[float] = None
    metadata: Dict = None

    def __post_init__(self):
        """Initialize optional fields."""
        if self.replies is None:
            self.replies = []
        if self.metadata is None:
            self.metadata = {}


@dataclass
class CompoundSentiment:
    """Model representing sentiment analysis for a compound."""

    compound_id: str
    sentiment_score: float  # -1.0 to 1.0
    mention_count: int
    positive_mentions: int = 0
    negative_mentions: int = 0
    neutral_mentions: int = 0
    metadata: Dict = None

    def __post_init__(self):
        """Initialize optional fields."""
        if self.metadata is None:
            self.metadata = {}


class Compound:
    """Model representing a chemical compound with web-enriched data."""

    def __init__(
        self,
        name: str,
        cas_number: Optional[str] = None,
    ) -> None:
        """Initialize compound.

        Args:
            name: Chemical compound name
            cas_number: Optional CAS registry number
        """
        self.name = name
        self.cas_number = cas_number

        # Web enrichment data
        self.community_data: Dict = {}
        self.reference_dois: Set[str] = set()
        self.reference_pmids: Set[str] = set()
        self.reference_urls: Dict[str, str] = {}

        # Social media data
        self.social_mentions: List[CompoundMention] = []
        self.social_discussions: List[CompoundDiscussion] = []
        self.social_sentiment: Optional[CompoundSentiment] = None

    def __repr__(self) -> str:
        """String representation."""
        return f"Compound(name='{self.name}', cas_number='{self.cas_number}')"
