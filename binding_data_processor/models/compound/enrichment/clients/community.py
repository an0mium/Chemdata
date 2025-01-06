"""Community data source client.

This module provides functionality to fetch and parse data from community sources:
- PsychonautWiki
- Erowid
- TripSit

Features:
- Circuit breaker pattern for resilience
- Caching with fallback
- Metrics collection
- Error tracking
- Enhanced text classification
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, TYPE_CHECKING
from datetime import datetime
import json

from bs4 import BeautifulSoup
from transformers import pipeline

from ....pipeline.infrastructure.circuit_breaker import CircuitConfig
from .base import BaseWebClient
from ...base.core import Compound

if TYPE_CHECKING:
    from .http import HTTPClient


class CommunityClient(BaseWebClient):
    """Client for community data sources."""

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

    def __init__(
        self,
        name: str = "community",
        http_client: Optional["HTTPClient"] = None,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        circuit_config: Optional[CircuitConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize community client.
        
        Args:
            name: Client name for circuit breaker
            http_client: Optional HTTP client to use
            model_dir: Optional directory for ML models
            cache_dir: Optional directory for caching
            circuit_config: Optional circuit breaker configuration
            logger: Optional logger instance
        """
        super().__init__(
            name=name,
            http_client=http_client,
            model_dir=model_dir,
            cache_dir=cache_dir,
            circuit_config=circuit_config,
            logger=logger,
        )

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

    def process_compounds(
        self,
        compounds: List[Compound],
        skip_predictions: bool = False,
        use_cache: bool = True,
    ) -> None:
        """Process list of compounds.
        
        Args:
            compounds: List of compounds to process
            skip_predictions: Whether to skip ML predictions
            use_cache: Whether to use cached results
        """
        for compound in compounds:
            try:
                data = self.get_compound_data(
                    compound.name,
                    compound.cas_number,
                    use_cache=use_cache,
                )
                if data:
                    compound.community_data = data
                    
                    # Add references
                    if "references" in data:
                        for ref in data["references"]:
                            if "doi" in ref:
                                compound.reference_dois.add(ref["doi"])
                            if "pubmed_id" in ref:
                                compound.reference_pmids.add(ref["pubmed_id"])
                            if "url" in ref:
                                compound.reference_urls[ref["title"]] = ref["url"]

                    self.processed_compounds.append(compound.name)
                else:
                    self.failed_compounds.append(compound.name)

            except Exception as e:
                self.logger.error(f"Error processing {compound.name}: {str(e)}")
                self.failed_compounds.append(compound.name)

    def get_compound_data(
        self,
        name: str,
        cas_number: Optional[str] = None,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Get data for a single compound.
        
        Args:
            name: Compound name
            cas_number: Optional CAS number
            use_cache: Whether to use cached results
            
        Returns:
            Dictionary of compound data or None if not found
        """
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
        """Get data from PsychonautWiki.
        
        Args:
            name: Compound name
            use_cache: Whether to use cached results
            
        Returns:
            Dictionary of PsychonautWiki data or None if not found
        """
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
                fallback=lambda: self._get_cached_psychonaut_data(name),
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
        """Get data from TripSit.
        
        Args:
            name: Compound name
            use_cache: Whether to use cached results
            
        Returns:
            Dictionary of TripSit data or None if not found
        """
        try:
            response = self.http.get(
                self.TRIPSIT_API,
                params={"name": name},
                use_cache=use_cache,
                fallback=lambda: self._get_cached_tripsit_data(name),
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
        """Get data from Erowid.
        
        Args:
            name: Compound name
            use_cache: Whether to use cached results
            
        Returns:
            Dictionary of Erowid data or None if not found
        """
        try:
            # Search for experience reports
            search_url = f"{self.EROWID_BASE}/search.php"
            response = self.http.get(
                search_url,
                params={"q": name},
                use_cache=use_cache,
                fallback=lambda: self._get_cached_erowid_data(name),
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
                        classification = self.text_classifier(
                            report["body_text"][:512]
                        )[0]
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
        """Get cached PsychonautWiki data.
        
        Args:
            name: Compound name
            
        Returns:
            Dictionary of cached data or None if not found
        """
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
        """Get cached TripSit data.
        
        Args:
            name: Compound name
            
        Returns:
            Dictionary of cached data or None if not found
        """
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
        """Get cached Erowid data.
        
        Args:
            name: Compound name
            
        Returns:
            Dictionary of cached data or None if not found
        """
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
        metrics = super().get_metrics()
        metrics.update({
            "processed_compounds": len(self.processed_compounds),
            "failed_compounds": len(self.failed_compounds),
            "success_rate": (
                len(self.processed_compounds) /
                (len(self.processed_compounds) + len(self.failed_compounds))
                if self.processed_compounds or self.failed_compounds
                else 0
            ),
            "source_stats": self.source_stats,
        })
        return metrics
