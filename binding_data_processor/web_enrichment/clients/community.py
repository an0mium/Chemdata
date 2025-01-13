"""Community data source client.

This module provides functionality to fetch and parse data from community sources:
- PsychonautWiki
- Erowid
- TripSit

Enhanced with:
- Circuit breaker pattern for resilience
- Better error handling and recovery
- Improved metrics collection
- Enhanced text classification
- Data validation and processing
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, TYPE_CHECKING
from datetime import datetime
import json
import re

from bs4 import BeautifulSoup
import pandas as pd
from transformers import pipeline

from .base_client import BaseWebClient, ValidationError
from ..models.compound import Compound
from ...pipeline.infrastructure.circuit_breaker import CircuitBreakerConfig

if TYPE_CHECKING:
    from .http_client_enhanced import HTTPClientEnhanced


class CommunityClient(BaseWebClient):
    """Enhanced client for community data sources."""

    # API endpoints
    PSYCHONAUT_API = "https://api.psychonautwiki.org"
    TRIPSIT_API = "https://tripbot.tripsit.me/api/tripsit/getDrug"
    EROWID_BASE = "https://erowid.org/experiences"

    # Required fields for validation
    REQUIRED_FIELDS = {
        "psychonaut": ["name", "effects", "roas"],
        "tripsit": ["name", "properties", "effects"],
        "erowid": ["title", "substance", "body_text"],
    }

    # Regex patterns for validation
    NAME_PATTERN = re.compile(r"^[a-zA-Z0-9\-\(\)\[\] ]+$")
    URL_PATTERN = re.compile(r"^https?://[^\s/$.?#].[^\s]*$")

    def __init__(
        self,
        http_client: Optional["HTTPClientEnhanced"] = None,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        logger: Optional[logging.Logger] = None,
        circuit_config: Optional[CircuitBreakerConfig] = None,
    ):
        """Initialize community client."""
        super().__init__(http_client, model_dir, cache_dir, logger)

        # Initialize text classifier for experience reports
        if model_dir:
            self.logger.info("Loading text classifier...")
            self.text_classifier = pipeline(
                "text-classification",
                model=str(model_dir / "experience_classifier"),
                device="cuda" if model_dir else "cpu",
            )
        else:
            self.text_classifier = None

        # Initialize tracking
        self.processed_compounds: List[str] = []
        self.failed_compounds: List[str] = []
        self.source_stats = {
            "psychonaut": {"success": 0, "failure": 0},
            "tripsit": {"success": 0, "failure": 0},
            "erowid": {"success": 0, "failure": 0},
        }

        # Initialize data processing
        self.reports_df = pd.DataFrame()

    def process_compounds(
        self,
        compounds: List[Compound],
        skip_predictions: bool = False,
        use_cache: bool = True,
    ) -> None:
        """Process list of compounds."""
        for compound in compounds:
            try:
                self._process_single_compound(compound, skip_predictions, use_cache)
            except ValidationError as e:
                self._handle_compound_error(compound, e, "Validation error")
            except Exception as e:
                self._handle_compound_error(compound, e, "Error processing")

    def _process_single_compound(
        self,
        compound: Compound,
        skip_predictions: bool,
        use_cache: bool,
    ) -> None:
        """Process a single compound."""
        # Validate compound name
        if not self._validate_name(compound.name):
            raise ValidationError(f"Invalid compound name: {compound.name}")

        # Get compound data
        data = self.get_compound_data(
            compound.name,
            compound.cas_number,
            use_cache=use_cache,
        )

        if not data:
            self.failed_compounds.append(compound.name)
            return

        # Update compound with data
        self._update_compound_data(compound, data)
        self.processed_compounds.append(compound.name)

    def _update_compound_data(self, compound: Compound, data: Dict[str, Any]) -> None:
        """Update compound with fetched data."""
        compound.community_data = data

        # Add references
        if "references" in data:
            self._add_compound_references(compound, data["references"])

        # Process Erowid reports if available
        if "erowid" in data and data["erowid"].get("reports"):
            self._process_reports(compound.name, data["erowid"]["reports"])

    def _add_compound_references(
        self,
        compound: Compound,
        references: List[Dict[str, Any]],
    ) -> None:
        """Add references to compound."""
        for ref in references:
            if "doi" in ref:
                compound.reference_dois.add(ref["doi"])
            if "pubmed_id" in ref:
                compound.reference_pmids.add(ref["pubmed_id"])
            if "url" in ref and self._validate_url(ref["url"]):
                compound.reference_urls[ref["title"]] = ref["url"]

    def _handle_compound_error(
        self,
        compound: Compound,
        error: Exception,
        prefix: str,
    ) -> None:
        """Handle errors during compound processing."""
        self.logger.error(f"{prefix} {compound.name}: {str(error)}")
        self.failed_compounds.append(compound.name)

    def _validate_name(self, name: str) -> bool:
        """Validate compound name format."""
        return bool(self.NAME_PATTERN.match(name))

    def _validate_url(self, url: str) -> bool:
        """Validate URL format."""
        return bool(self.URL_PATTERN.match(url))

    def _process_reports(self, compound_name: str, reports: List[Dict[str, Any]]) -> None:
        """Process and analyze experience reports."""
        # Convert reports to DataFrame
        df = pd.DataFrame(reports)
        df["compound"] = compound_name

        # Parse dates
        df["date"] = pd.to_datetime(df["date"], errors="coerce")

        # Add report classification if available
        if "classification" in df.columns:
            df["sentiment"] = df["classification"].apply(lambda x: x.get("label"))
            df["confidence"] = df["classification"].apply(lambda x: x.get("score"))

        # Calculate text statistics
        df["word_count"] = df["body_text"].str.split().str.len()
        df["avg_word_length"] = df["body_text"].str.split().apply(lambda x: sum(len(w) for w in x) / len(x) if x else 0)

        # Append to main DataFrame
        self.reports_df = pd.concat([self.reports_df, df], ignore_index=True)

    def get_compound_data(
        self,
        name: str,
        cas_number: Optional[str] = None,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Get data for a single compound."""
        data = {}

        # Get PsychonautWiki data
        try:
            psychonaut_data = self._get_psychonaut_data(name, use_cache)
            if psychonaut_data:
                data["psychonaut"] = psychonaut_data
                self.source_stats["psychonaut"]["success"] += 1
            else:
                self.source_stats["psychonaut"]["failure"] += 1
        except Exception as e:
            self.logger.error(f"Error getting PsychonautWiki data: {str(e)}")
            self.source_stats["psychonaut"]["failure"] += 1

        # Get TripSit data
        try:
            tripsit_data = self._get_tripsit_data(name, use_cache)
            if tripsit_data:
                data["tripsit"] = tripsit_data
                self.source_stats["tripsit"]["success"] += 1
            else:
                self.source_stats["tripsit"]["failure"] += 1
        except Exception as e:
            self.logger.error(f"Error getting TripSit data: {str(e)}")
            self.source_stats["tripsit"]["failure"] += 1

        # Get Erowid data
        try:
            erowid_data = self._get_erowid_data(name, use_cache)
            if erowid_data:
                data["erowid"] = erowid_data
                self.source_stats["erowid"]["success"] += 1
            else:
                self.source_stats["erowid"]["failure"] += 1
        except Exception as e:
            self.logger.error(f"Error getting Erowid data: {str(e)}")
            self.source_stats["erowid"]["failure"] += 1

        return data if data else None

    def _get_psychonaut_data(
        self,
        name: str,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Get data from PsychonautWiki."""
        query = """
        query Substance($query: String) {
            substances(query: $query) {
                name
                commonNames
                class {
                    chemical
                    psychoactive
                }
                tolerance {
                    full
                    half
                    zero
                }
                roas {
                    name
                    dose {
                        units
                        threshold
                        light {
                            min
                            max
                        }
                        common {
                            min
                            max
                        }
                        strong {
                            min
                            max
                        }
                        heavy
                    }
                    duration {
                        onset {
                            min
                            max
                            units
                        }
                        comeup {
                            min
                            max
                            units
                        }
                        peak {
                            min
                            max
                            units
                        }
                        offset {
                            min
                            max
                            units
                        }
                        total {
                            min
                            max
                            units
                        }
                        afterglow {
                            min
                            max
                            units
                        }
                    }
                    bioavailability {
                        min
                        max
                    }
                }
                effects {
                    name
                    url
                    experience {
                        positive
                        neutral
                        negative
                    }
                }
                interactions {
                    status
                    note
                }
            }
        }
        """

        try:
            response = self.http.post(
                self.PSYCHONAUT_API,
                json_data={"query": query, "variables": {"query": name}},
                use_cache=use_cache,
                fallback=self._get_cached_psychonaut_data,
            )
            data = response.json()

            if "data" in data and "substances" in data["data"]:
                substances = data["data"]["substances"]
                if substances:
                    substance = substances[0]
                    self._validate_response(substance, self.REQUIRED_FIELDS["psychonaut"])
                    return substance

        except Exception as e:
            self.logger.error(f"Error querying PsychonautWiki: {str(e)}")

        return None

    def _get_tripsit_data(
        self,
        name: str,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Get data from TripSit."""
        try:
            response = self.http.get(
                self.TRIPSIT_API,
                params={"name": name},
                use_cache=use_cache,
                fallback=self._get_cached_tripsit_data,
            )
            data = response.json()

            if "data" in data and data["data"]:
                substance = data["data"][0]
                self._validate_response(substance, self.REQUIRED_FIELDS["tripsit"])
                return substance

        except Exception as e:
            self.logger.error(f"Error querying TripSit: {str(e)}")

        return None

    def _get_erowid_data(
        self,
        name: str,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Get data from Erowid."""
        try:
            # Search for experience reports
            search_url = f"{self.EROWID_BASE}/search.php"
            response = self.http.get(
                search_url,
                params={"q": name},
                use_cache=use_cache,
                fallback=self._get_cached_erowid_data,
            )

            # Parse search results
            soup = BeautifulSoup(response.text, "html.parser")
            reports = []

            for result in soup.find_all("div", class_="experience-report"):
                try:
                    report = {
                        "title": result.find("h3").text.strip(),
                        "substance": result.find("div", class_="substance").text.strip(),
                        "author": result.find("div", class_="author").text.strip(),
                        "date": result.find("div", class_="date").text.strip(),
                        "body_text": result.find("div", class_="body").text.strip(),
                    }

                    # Validate report
                    self._validate_response(report, self.REQUIRED_FIELDS["erowid"])

                    # Classify report if model available
                    if self.text_classifier:
                        classification = self.text_classifier(report["body_text"][:512])[0]
                        report["classification"] = {
                            "label": classification["label"],
                            "score": classification["score"],
                        }

                    reports.append(report)

                except Exception as e:
                    self.logger.error(f"Error parsing report: {str(e)}")
                    continue

            if reports:
                return {
                    "reports": reports,
                    "total_reports": len(reports),
                    "last_updated": datetime.now().isoformat(),
                }

        except Exception as e:
            self.logger.error(f"Error querying Erowid: {str(e)}")

        return None

    def _get_cached_psychonaut_data(self, name: str) -> Optional[Dict[str, Any]]:
        """Get cached PsychonautWiki data."""
        if not self.cache_dir:
            return None

        try:
            cache_file = self.cache_dir / f"psychonaut_{name}.json"
            if not cache_file.exists():
                return None

            with cache_file.open() as f:
                return json.load(f)

        except Exception as e:
            self.logger.error(f"Error reading cached PsychonautWiki data: {str(e)}")
            return None

    def _get_cached_tripsit_data(self, name: str) -> Optional[Dict[str, Any]]:
        """Get cached TripSit data."""
        if not self.cache_dir:
            return None

        try:
            cache_file = self.cache_dir / f"tripsit_{name}.json"
            if not cache_file.exists():
                return None

            with cache_file.open() as f:
                return json.load(f)

        except Exception as e:
            self.logger.error(f"Error reading cached TripSit data: {str(e)}")
            return None

    def _get_cached_erowid_data(self, name: str) -> Optional[Dict[str, Any]]:
        """Get cached Erowid data."""
        if not self.cache_dir:
            return None

        try:
            cache_file = self.cache_dir / f"erowid_{name}.json"
            if not cache_file.exists():
                return None

            with cache_file.open() as f:
                return json.load(f)

        except Exception as e:
            self.logger.error(f"Error reading cached Erowid data: {str(e)}")
            return None

    def get_metrics(self) -> Dict[str, Any]:
        """Get client metrics."""
        total = len(self.processed_compounds) + len(self.failed_compounds)
        success_rate = len(self.processed_compounds) / total if total > 0 else 0

        return {
            "processed_compounds": len(self.processed_compounds),
            "failed_compounds": len(self.failed_compounds),
            "success_rate": success_rate,
            "source_stats": self.source_stats,
            "reports_stats": self._get_reports_stats(),
        }

    def _get_reports_stats(self) -> Dict[str, Any]:
        """Get statistics about processed reports."""
        if self.reports_df.empty:
            return {}

        stats = {
            "total_reports": len(self.reports_df),
            "compounds_with_reports": self.reports_df["compound"].nunique(),
            "avg_reports_per_compound": (len(self.reports_df) / self.reports_df["compound"].nunique()),
            "avg_word_count": self.reports_df["word_count"].mean(),
            "avg_word_length": self.reports_df["avg_word_length"].mean(),
        }

        # Add date range if available
        if "date" in self.reports_df.columns:
            stats["date_range"] = {
                "start": self.reports_df["date"].min().isoformat(),
                "end": self.reports_df["date"].max().isoformat(),
            }

        # Add sentiment stats if available
        if "sentiment" in self.reports_df.columns:
            stats["sentiment_stats"] = self.reports_df["sentiment"].value_counts().to_dict()

        return stats
