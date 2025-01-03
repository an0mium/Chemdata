"""Swiss* web services integration.

This module handles:
1. SwissTargetPrediction - Target prediction with confidence scores
2. SwissADME - ADME properties and pharmacokinetics
3. SwissSimilarity - Chemical similarity search

Features:
- Progress bars for long-running operations
- Comprehensive error handling
- Detailed property extraction
- Source ID and URL extraction
- Enhanced documentation
"""

import os
import re
import json
import time
from typing import Dict, List, Optional, Any
from urllib.parse import quote
from bs4 import BeautifulSoup
from tqdm import tqdm

from ..http_client import HttpClient
from logger import LogManager

logger = LogManager().get_logger("web_enrichment.data_sources.swiss")


class SwissClient:
    """Client for Swiss* web services with enhanced functionality."""

    # Base URLs for Swiss* services
    BASE_URLS = {
        "target": "http://www.swisstargetprediction.ch",
        "adme": "http://www.swissadme.ch",
        "similarity": "http://www.swisssimilarity.ch",
    }

    def __init__(self, http_client: HttpClient):
        """Initialize Swiss* client."""
        self.http = http_client

    def get_target_predictions(
        self, smiles: str, organism: str = "Homo sapiens"
    ) -> Dict[str, Any]:
        """
        Get target predictions from SwissTargetPrediction.

        Args:
            smiles: SMILES string of compound
            organism: Target organism (default: Homo sapiens)

        Returns:
            Dictionary containing:
            - predictions: List of predicted targets with probabilities
            - url: Results page URL
            Each prediction contains:
            - target: Target protein name
            - common_name: Common name of target
            - uniprot_id: UniProt ID
            - probability: Prediction confidence (0-1)
            - target_class: Protein family/class
            - chembl_id: ChEMBL ID if available
            - gene_name: Gene name if available
        """
        try:
            # Submit prediction request
            submit_url = f"{self.BASE_URLS['target']}/predict/"
            data = {"smiles": smiles, "organism": organism}

            response = self.http.make_request(submit_url, method="POST", data=data)

            if not response:
                return {"predictions": [], "url": None}

            # Get job ID from response
            soup = BeautifulSoup(response.text, "html.parser")
            job_id = soup.find("input", {"name": "job"})["value"]

            # Poll for results with progress bar
            results_url = f"{self.BASE_URLS['target']}/results/{job_id}"
            max_attempts = 30
            progress = tqdm(
                total=max_attempts,
                desc="Waiting for target predictions",
                unit="attempts",
            )

            for attempt in range(max_attempts):
                response = self.http.make_request(results_url)
                if not response:
                    time.sleep(2)
                    progress.update(1)
                    continue

                if "Computation in progress" not in response.text:
                    progress.close()
                    break

                time.sleep(2)
                progress.update(1)
                progress.set_description(f"Attempt {attempt + 1}/{max_attempts}")

            else:
                progress.close()
                logger.warning("Target prediction timed out")
                return {"predictions": [], "url": None}

            # Parse results
            soup = BeautifulSoup(response.text, "html.parser")
            predictions = []

            result_rows = soup.find_all("tr", {"class": "result-row"})
            if result_rows:
                parse_progress = tqdm(
                    result_rows, desc="Parsing predictions", unit="predictions"
                )

                for row in parse_progress:
                    cols = row.find_all("td")
                    if len(cols) >= 4:
                        prediction = {
                            "target": cols[0].text.strip(),
                            "common_name": cols[1].text.strip(),
                            "uniprot_id": cols[2].text.strip(),
                            "probability": float(cols[3].text.strip()),
                            "target_class": (
                                cols[4].text.strip() if len(cols) > 4 else "Unknown"
                            ),
                            "chembl_id": None,
                            "gene_name": None,
                        }

                        # Extract additional IDs from links
                        chembl_link = cols[0].find("a", href=re.compile(r"chembl"))
                        if chembl_link:
                            chembl_match = re.search(r"CHEMBL\d+", chembl_link["href"])
                            if chembl_match:
                                prediction["chembl_id"] = chembl_match.group()

                        uniprot_link = cols[2].find("a", href=re.compile(r"uniprot"))
                        if uniprot_link:
                            prediction["gene_name"] = uniprot_link.get(
                                "title", ""
                            ).split("_")[0]

                        predictions.append(prediction)

                parse_progress.close()

            return {
                "predictions": predictions,
                "url": f"{self.BASE_URLS['target']}/result.php?job={job_id}",
            }

        except Exception as e:
            logger.error(f"Error getting target predictions: {str(e)}")
            return {"predictions": [], "url": None}

    def get_adme_properties(self, smiles: str) -> Dict[str, Any]:
        """
        Get ADME properties from SwissADME.

        Args:
            smiles: SMILES string of compound

        Returns:
            Dictionary containing:
            - physicochemical: Basic physical chemistry properties
            - lipophilicity: LogP values and consensus
            - solubility: Water solubility predictions
            - pharmacokinetics: ADME predictions
            - druglikeness: Drug-likeness scores and rules
            - medicinal_chemistry: Med chem friendliness
            - url: Results page URL
        """
        try:
            # Submit calculation request
            submit_url = f"{self.BASE_URLS['adme']}/calculate/"
            data = {"smiles": smiles}

            response = self.http.make_request(submit_url, method="POST", data=data)

            if not response:
                return {}

            # Get job ID from response
            soup = BeautifulSoup(response.text, "html.parser")
            job_id = soup.find("input", {"name": "job"})["value"]

            # Poll for results with progress bar
            results_url = f"{self.BASE_URLS['adme']}/results/{job_id}"
            max_attempts = 30
            progress = tqdm(
                total=max_attempts, desc="Calculating ADME properties", unit="attempts"
            )

            for attempt in range(max_attempts):
                response = self.http.make_request(results_url)
                if not response:
                    time.sleep(2)
                    progress.update(1)
                    continue

                if "Calculation in progress" not in response.text:
                    progress.close()
                    break

                time.sleep(2)
                progress.update(1)
                progress.set_description(f"Attempt {attempt + 1}/{max_attempts}")

            else:
                progress.close()
                logger.warning("ADME calculation timed out")
                return {}

            # Parse results
            soup = BeautifulSoup(response.text, "html.parser")
            properties = {
                "physicochemical": {
                    "mw": None,
                    "logp": None,
                    "hbd": None,
                    "hba": None,
                    "tpsa": None,
                    "rotatable_bonds": None,
                },
                "lipophilicity": {"logp_values": [], "consensus": None},
                "solubility": {"log_s_values": [], "consensus": None},
                "absorption": {
                    "gi_absorption": None,
                    "bbb_permeant": None,
                    "pgp_substrate": None,
                },
                "metabolism": {
                    "cyp_inhibition": {
                        "cyp1a2": None,
                        "cyp2c19": None,
                        "cyp2c9": None,
                        "cyp2d6": None,
                        "cyp3a4": None,
                    }
                },
                "druglikeness": {},
                "medicinal_chemistry": {},
                "url": f"{self.BASE_URLS['adme']}/result.php?job={job_id}",
            }

            # Extract properties with progress bar
            sections = [
                ("physchem", "physicochemical"),
                ("lipophilicity", "lipophilicity"),
                ("solubility", "solubility"),
                ("pharmacokinetics", "absorption"),
                ("druglikeness", "druglikeness"),
                ("medchem", "medicinal_chemistry"),
            ]

            section_progress = tqdm(
                sections, desc="Parsing property sections", unit="sections"
            )

            for html_id, prop_key in section_progress:
                section = soup.find("div", {"id": html_id})
                if section:
                    for row in section.find_all("tr"):
                        cols = row.find_all("td")
                        if len(cols) >= 2:
                            key = cols[0].text.strip()
                            value = cols[1].text.strip()

                            if prop_key in ["lipophilicity", "solubility"]:
                                if key == "Consensus":
                                    properties[prop_key]["consensus"] = value
                                else:
                                    properties[prop_key][f"{prop_key}_values"].append(
                                        {"method": key, "value": value}
                                    )
                            else:
                                properties[prop_key][key] = value

            section_progress.close()

            return properties

        except Exception as e:
            logger.error(f"Error getting ADME properties: {str(e)}")
            return {}

    def search_similar_compounds(
        self, smiles: str, similarity_threshold: float = 0.7, max_results: int = 100
    ) -> Dict[str, Any]:
        """
        Search for similar compounds using SwissSimilarity.

        Args:
            smiles: SMILES string of query compound
            similarity_threshold: Minimum similarity score (0-1)
            max_results: Maximum number of results to return

        Returns:
            Dictionary containing:
            - similar_compounds: List of similar compounds with scores
            - url: Results page URL
            Each compound contains:
            - name: Compound name
            - smiles: SMILES string
            - similarity: Tanimoto similarity score (0-1)
            - source: Database source
            - source_id: Database ID if available
            - source_url: Database URL if available
        """
        try:
            # Submit search request
            submit_url = f"{self.BASE_URLS['similarity']}/search/"
            data = {
                "smiles": smiles,
                "threshold": similarity_threshold,
                "limit": max_results,
            }

            response = self.http.make_request(submit_url, method="POST", data=data)

            if not response:
                return {"similar_compounds": [], "url": None}

            # Get job ID from response
            soup = BeautifulSoup(response.text, "html.parser")
            job_id = soup.find("input", {"name": "job"})["value"]

            # Poll for results with progress bar
            results_url = f"{self.BASE_URLS['similarity']}/results/{job_id}"
            max_attempts = 30
            progress = tqdm(
                total=max_attempts, desc="Searching similar compounds", unit="attempts"
            )

            for attempt in range(max_attempts):
                response = self.http.make_request(results_url)
                if not response:
                    time.sleep(2)
                    progress.update(1)
                    continue

                if "Search in progress" not in response.text:
                    progress.close()
                    break

                time.sleep(2)
                progress.update(1)
                progress.set_description(f"Attempt {attempt + 1}/{max_attempts}")

            else:
                progress.close()
                logger.warning("Similarity search timed out")
                return {"similar_compounds": [], "url": None}

            # Parse results
            soup = BeautifulSoup(response.text, "html.parser")
            compounds = []

            result_rows = soup.find_all("tr", {"class": "result-row"})
            if result_rows:
                parse_progress = tqdm(
                    result_rows, desc="Parsing similar compounds", unit="compounds"
                )

                for row in parse_progress:
                    cols = row.find_all("td")
                    if len(cols) >= 4:
                        compound = {
                            "name": cols[0].text.strip(),
                            "smiles": cols[1].text.strip(),
                            "similarity": float(cols[2].text.strip()),
                            "source": cols[3].text.strip(),
                            "source_id": None,
                            "source_url": None,
                        }

                        # Extract source IDs and URLs
                        source_link = cols[3].find("a")
                        if source_link:
                            compound["source_url"] = source_link["href"]
                            if "chembl" in compound["source"].lower():
                                chembl_match = re.search(
                                    r"CHEMBL\d+", source_link["href"]
                                )
                                if chembl_match:
                                    compound["source_id"] = chembl_match.group()
                            elif "pubchem" in compound["source"].lower():
                                pubchem_match = re.search(
                                    r"compound/(\d+)", source_link["href"]
                                )
                                if pubchem_match:
                                    compound["source_id"] = pubchem_match.group(1)

                        compounds.append(compound)

                parse_progress.close()

            return {
                "similar_compounds": compounds,
                "url": f"{self.base_urls['similarity']}/result.php?job={job_id}",
            }

        except Exception as e:
            logger.error(f"Error searching similar compounds: {str(e)}")
            return {"similar_compounds": [], "url": None}
