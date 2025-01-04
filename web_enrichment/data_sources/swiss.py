"""Swiss* web services integration.

This module handles:
1. SwissTargetPrediction - Target prediction with confidence scores
2. SwissADME - ADME properties and pharmacokinetics
3. SwissSimilarity - Chemical similarity search
4. SwissDock - Molecular docking predictions
5. SwissBioisostere - Bioisostere suggestions
6. SwissAI - AI-powered predictions
7. SwissParam - Force field parameters
8. SwissSideChain - Side chain modeling

Features:
- Progress bars for long-running operations
- Comprehensive error handling
- Detailed property extraction
- Source ID and URL extraction
- Enhanced documentation
- Parallel processing for batch queries
- Caching of results
- ML integration
- Data enrichment
"""

import os
import re
import json
import time
from typing import Dict, List, Optional, Any, Set, Tuple
from urllib.parse import quote
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from datetime import datetime, timedelta
from dataclasses import asdict

import numpy as np
from bs4 import BeautifulSoup
from tqdm import tqdm
import torch
from transformers import pipeline
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from ..http_client import HttpClient
from logger import LogManager
from binding_data_processor.models.compound import CompoundData
from binding_data_processor.processors.structure.ml.predictors.activity import (
    ActivityPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.toxicity import (
    ToxicityPredictor,
)
from binding_data_processor.processors.structure.ml.predictors.abuse import (
    AbusePotentialPredictor,
)

logger = LogManager().get_logger("web_enrichment.data_sources.swiss")


class SwissClient:
    """Client for Swiss* web services with enhanced functionality."""

    # Base URLs for Swiss* services
    BASE_URLS = {
        "target": "http://www.swisstargetprediction.ch",
        "adme": "http://www.swissadme.ch",
        "similarity": "http://www.swisssimilarity.ch",
        "dock": "http://www.swissdock.ch",
        "bioisostere": "http://www.swissbioisostere.ch",
        "ai": "http://www.swissai.ch",
        "param": "http://www.swissparam.ch",
        "sidechain": "http://www.swisssidechain.ch",
    }

    # Cache settings
    CACHE_DIR = Path(".cache/swiss")
    CACHE_DURATION = timedelta(days=30)  # Cache results for 30 days

    def __init__(
        self,
        http_client: HttpClient,
        model_dir: str = "../models",
        cache_dir: Optional[str] = None,
        n_workers: int = 4,
        max_retries: int = 3,
        device: str = "cuda" if torch.cuda.is_available() else "cpu",
    ):
        """Initialize Swiss* client.

        Args:
            http_client: HTTP client for making requests
            model_dir: Directory containing ML models
            cache_dir: Optional custom cache directory
            n_workers: Number of worker threads
            max_retries: Maximum number of retries for failed operations
            device: Device to run ML models on
        """
        self.http = http_client
        self.n_workers = n_workers
        self.max_retries = max_retries
        self.device = device

        # Initialize cache
        if cache_dir:
            self.CACHE_DIR = Path(cache_dir)
        self.CACHE_DIR.mkdir(parents=True, exist_ok=True)

        # Initialize ML models
        logger.info("Loading ML models...")
        self.activity_predictor = ActivityPredictor(model_dir=f"{model_dir}/activity")
        self.toxicity_predictor = ToxicityPredictor(model_dir=f"{model_dir}/toxicity")
        self.abuse_predictor = AbusePotentialPredictor(model_dir=f"{model_dir}/abuse")

        # Initialize thread pools
        self.prediction_executor = ThreadPoolExecutor(
            max_workers=n_workers, thread_name_prefix="prediction"
        )
        self.web_executor = ThreadPoolExecutor(
            max_workers=n_workers * 2, thread_name_prefix="web"
        )  # More threads for web IO

    def _get_cached_result(self, cache_key: str) -> Optional[Dict[str, Any]]:
        """Get cached result if available and not expired."""
        cache_file = self.CACHE_DIR / f"{cache_key}.json"
        if cache_file.exists():
            try:
                data = json.loads(cache_file.read_text())
                cached_time = datetime.fromisoformat(data["cached_at"])
                if datetime.now() - cached_time < self.CACHE_DURATION:
                    return data["result"]
            except Exception as e:
                logger.error(f"Error reading cache file {cache_file}: {e}")
        return None

    def _cache_result(self, cache_key: str, result: Dict[str, Any]) -> None:
        """Cache result with timestamp."""
        try:
            cache_file = self.CACHE_DIR / f"{cache_key}.json"
            data = {
                "cached_at": datetime.now().isoformat(),
                "result": result,
            }
            cache_file.write_text(json.dumps(data, indent=2))
        except Exception as e:
            logger.error(f"Error writing cache file {cache_file}: {e}")

    def process_compounds(
        self,
        compounds: List[CompoundData],
        skip_predictions: bool = False,
        skip_docking: bool = False,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Process a batch of compounds through Swiss* services.

        Args:
            compounds: List of compounds to process
            skip_predictions: Whether to skip ML predictions
            skip_docking: Whether to skip docking predictions
            use_cache: Whether to use cached results

        Returns:
            Dictionary containing processing statistics
        """
        stats = {
            "total": len(compounds),
            "processed": 0,
            "errors": [],
            "start_time": datetime.now().isoformat(),
        }

        try:
            # Process compounds in parallel
            with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
                futures = []
                for compound in compounds:
                    futures.append(
                        executor.submit(
                            self._process_compound,
                            compound,
                            skip_predictions,
                            skip_docking,
                            use_cache,
                        )
                    )

                # Collect results with progress bar
                with tqdm(
                    total=len(futures),
                    desc="Processing compounds",
                    unit="compounds",
                ) as pbar:
                    for future in as_completed(futures):
                        try:
                            result = future.result()
                            if result:
                                stats["processed"] += 1
                        except Exception as e:
                            stats["errors"].append(str(e))
                        pbar.update(1)

        except Exception as e:
            logger.error(f"Error processing compounds: {str(e)}")
            stats["errors"].append(str(e))

        stats["end_time"] = datetime.now().isoformat()
        return stats

    def _process_compound(
        self,
        compound: CompoundData,
        skip_predictions: bool,
        skip_docking: bool,
        use_cache: bool,
    ) -> bool:
        """Process a single compound through Swiss* services.

        Args:
            compound: Compound to process
            skip_predictions: Whether to skip ML predictions
            skip_docking: Whether to skip docking predictions
            use_cache: Whether to use cached results

        Returns:
            Whether processing was successful
        """
        try:
            if not compound.smiles:
                return False

            # Get target predictions
            target_data = self.get_target_predictions(
                compound.smiles,
                use_cache=use_cache,
            )
            if target_data and target_data["predictions"]:
                compound.target_predictions = target_data["predictions"]
                compound.urls["swiss_target"] = target_data["url"]

            # Get ADME properties
            adme_data = self.get_adme_properties(
                compound.smiles,
                use_cache=use_cache,
            )
            if adme_data:
                compound.adme_properties = adme_data
                compound.urls["swiss_adme"] = adme_data["url"]

            # Get similar compounds
            similarity_data = self.search_similar_compounds(
                compound.smiles,
                use_cache=use_cache,
            )
            if similarity_data and similarity_data["similar_compounds"]:
                compound.similar_compounds = similarity_data["similar_compounds"]
                compound.urls["swiss_similarity"] = similarity_data["url"]

            # Get bioisosteres
            bioisostere_data = self.get_bioisosteres(
                compound.smiles,
                use_cache=use_cache,
            )
            if bioisostere_data and bioisostere_data["replacements"]:
                compound.bioisosteres = bioisostere_data["replacements"]
                compound.urls["swiss_bioisostere"] = bioisostere_data["url"]

            # Get docking predictions if targets available
            if (
                not skip_docking
                and compound.target_predictions
                and len(compound.target_predictions) > 0
            ):
                top_target = compound.target_predictions[0]
                if top_target.get("pdb_ids"):
                    docking_data = self.get_docking_predictions(
                        compound.smiles,
                        top_target["pdb_ids"][0],
                        use_cache=use_cache,
                    )
                    if docking_data and docking_data["poses"]:
                        compound.docking_poses = docking_data["poses"]
                        compound.urls["swiss_dock"] = docking_data["url"]

            # Add ML predictions
            if not skip_predictions:
                try:
                    # Activity predictions
                    activity_pred = self.activity_predictor.predict(compound.smiles)
                    if activity_pred:
                        compound.activity_predictions = activity_pred

                    # Toxicity predictions
                    toxicity_pred = self.toxicity_predictor.predict(compound.smiles)
                    if toxicity_pred:
                        compound.toxicity_predictions = toxicity_pred

                    # Abuse potential predictions
                    abuse_pred = self.abuse_predictor.predict(compound.smiles)
                    if abuse_pred:
                        compound.abuse_predictions = abuse_pred

                except Exception as e:
                    logger.error(f"Error getting predictions: {str(e)}")

            return True

        except Exception as e:
            logger.error(f"Error processing compound {compound.name}: {str(e)}")
            return False

    def get_target_predictions(
        self,
        smiles: str,
        organism: str = "Homo sapiens",
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Get target predictions from SwissTargetPrediction.

        Args:
            smiles: SMILES string of compound
            organism: Target organism (default: Homo sapiens)
            use_cache: Whether to use cached results

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
            - pdb_ids: List of related PDB IDs
            - binding_site: Predicted binding site if available
        """
        # Check cache first
        if use_cache:
            cache_key = f"target_{smiles}_{organism}"
            cached = self._get_cached_result(cache_key)
            if cached:
                return cached

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
                            "pdb_ids": [],
                            "binding_site": None,
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

                            # Get PDB IDs from UniProt
                            try:
                                uniprot_id = prediction["uniprot_id"]
                                pdb_url = (
                                    f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
                                )
                                pdb_response = self.http.make_request(pdb_url)
                                if pdb_response:
                                    pdb_data = pdb_response.json()
                                    for structure in pdb_data.get("structures", []):
                                        prediction["pdb_ids"].append(structure["id"])
                                        if not prediction[
                                            "binding_site"
                                        ] and structure.get("ligands"):
                                            prediction["binding_site"] = structure[
                                                "ligands"
                                            ][0].get("binding_site")
                            except Exception as e:
                                logger.error(f"Error getting PDB data: {str(e)}")

                        predictions.append(prediction)

                parse_progress.close()

            result = {
                "predictions": predictions,
                "url": f"{self.BASE_URLS['target']}/result.php?job={job_id}",
            }

            # Cache result
            if use_cache:
                self._cache_result(cache_key, result)

            return result

        except Exception as e:
            logger.error(f"Error getting target predictions: {str(e)}")
            return {"predictions": [], "url": None}

    def get_adme_properties(
        self, smiles: str, use_cache: bool = True
    ) -> Dict[str, Any]:
        """Get ADME properties from SwissADME.

        Args:
            smiles: SMILES string of compound
            use_cache: Whether to use cached results

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
        # Check cache first
        if use_cache:
            cache_key = f"adme_{smiles}"
            cached = self._get_cached_result(cache_key)
            if cached:
                return cached

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

            # Cache result
            if use_cache:
                self._cache_result(cache_key, properties)

            return properties

        except Exception as e:
            logger.error(f"Error getting ADME properties: {str(e)}")
            return {}

    def search_similar_compounds(
        self,
        smiles: str,
        similarity_threshold: float = 0.7,
        max_results: int = 100,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Search for similar compounds using SwissSimilarity.

        Args:
            smiles: SMILES string of query compound
            similarity_threshold: Minimum similarity score (0-1)
            max_results: Maximum number of results to return
            use_cache: Whether to use cached results

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
            - properties: Basic properties if available
            - bioactivity: Known bioactivity data if available
        """
        # Check cache first
        if use_cache:
            cache_key = f"similarity_{smiles}_{similarity_threshold}_{max_results}"
            cached = self._get_cached_result(cache_key)
            if cached:
                return cached

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
                            "properties": {},
                            "bioactivity": [],
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
                                    # Get ChEMBL bioactivity data
                                    try:
                                        chembl_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{compound['source_id']}/bioactivities"
                                        chembl_response = self.http.make_request(
                                            chembl_url
                                        )
                                        if chembl_response:
                                            bioactivities = chembl_response.json()
                                            compound["bioactivity"] = bioactivities
                                    except Exception as e:
                                        logger.error(
                                            f"Error getting ChEMBL data: {str(e)}"
                                        )
                            elif "pubchem" in compound["source"].lower():
                                pubchem_match = re.search(
                                    r"compound/(\d+)", source_link["href"]
                                )
                                if pubchem_match:
                                    compound["source_id"] = pubchem_match.group(1)
                                    # Get PubChem properties
                                    try:
                                        pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{compound['source_id']}/property/MolecularWeight,XLogP,TPSA/JSON"
                                        pubchem_response = self.http.make_request(
                                            pubchem_url
                                        )
                                        if pubchem_response:
                                            properties = pubchem_response.json()
                                            compound["properties"] = properties
                                    except Exception as e:
                                        logger.error(
                                            f"Error getting PubChem data: {str(e)}"
                                        )

                        compounds.append(compound)

                parse_progress.close()

            result = {
                "similar_compounds": compounds,
                "url": f"{self.BASE_URLS['similarity']}/result.php?job={job_id}",
            }

            # Cache result
            if use_cache:
                self._cache_result(cache_key, result)

            return result

        except Exception as e:
            logger.error(f"Error searching similar compounds: {str(e)}")
            return {"similar_compounds": [], "url": None}

    def get_docking_predictions(
        self,
        smiles: str,
        target_pdb: str,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Get molecular docking predictions from SwissDock.

        Args:
            smiles: SMILES string of compound
            target_pdb: PDB ID of target protein
            use_cache: Whether to use cached results

        Returns:
            Dictionary containing docking predictions and scores
        """
        # Check cache first
        if use_cache:
            cache_key = f"dock_{smiles}_{target_pdb}"
            cached = self._get_cached_result(cache_key)
            if cached:
                return cached

        try:
            # Submit docking request
            submit_url = f"{self.BASE_URLS['dock']}/docking/"
            data = {
                "smiles": smiles,
                "target": target_pdb,
            }

            response = self.http.make_request(submit_url, method="POST", data=data)

            if not response:
                return {}

            # Get job ID and poll for results
            soup = BeautifulSoup(response.text, "html.parser")
            job_id = soup.find("input", {"name": "job"})["value"]

            results_url = f"{self.BASE_URLS['dock']}/results/{job_id}"
            max_attempts = 60  # Docking takes longer
            progress = tqdm(
                total=max_attempts, desc="Running molecular docking", unit="attempts"
            )

            for attempt in range(max_attempts):
                response = self.http.make_request(results_url)
                if not response:
                    time.sleep(5)  # Longer sleep for docking
                    progress.update(1)
                    continue

                if "Docking in progress" not in response.text:
                    progress.close()
                    break

                time.sleep(5)
                progress.update(1)
                progress.set_description(f"Attempt {attempt + 1}/{max_attempts}")

            else:
                progress.close()
                logger.warning("Docking timed out")
                return {}

            # Parse results
            soup = BeautifulSoup(response.text, "html.parser")
            results = {
                "poses": [],
                "url": f"{self.BASE_URLS['dock']}/result.php?job={job_id}",
            }

            pose_divs = soup.find_all("div", {"class": "pose"})
            if pose_divs:
                parse_progress = tqdm(
                    pose_divs, desc="Parsing docking poses", unit="poses"
                )

                for div in parse_progress:
                    pose = {
                        "cluster": int(div.get("data-cluster", -1)),
                        "rank": int(div.get("data-rank", -1)),
                        "score": float(div.get("data-score", 0)),
                        "binding_energy": float(div.get("data-energy", 0)),
                        "interactions": [],
                    }

                    # Extract interaction details
                    interactions_div = div.find("div", {"class": "interactions"})
                    if interactions_div:
                        for interaction in interactions_div.find_all("div"):
                            pose["interactions"].append(
                                {
                                    "type": interaction.get("data-type"),
                                    "residue": interaction.get("data-residue"),
                                    "distance": float(
                                        interaction.get("data-distance", 0)
                                    ),
                                }
                            )

                    results["poses"].append(pose)

                parse_progress.close()

            # Cache result
            if use_cache:
                self._cache_result(cache_key, results)

            return results

        except Exception as e:
            logger.error(f"Error getting docking predictions: {str(e)}")
            return {}

    def get_bioisosteres(
        self,
        smiles: str,
        min_occurrence: int = 5,
        use_cache: bool = True,
    ) -> Dict[str, Any]:
        """Get bioisostere suggestions from SwissBioisostere.

        Args:
            smiles: SMILES string of compound
            min_occurrence: Minimum number of occurrences in database
            use_cache: Whether to use cached results

        Returns:
            Dictionary containing bioisostere suggestions and statistics
        """
        # Check cache first
        if use_cache:
            cache_key = f"bioisostere_{smiles}_{min_occurrence}"
            cached = self._get_cached_result(cache_key)
            if cached:
                return cached

        try:
            # Submit search request
            submit_url = f"{self.BASE_URLS['bioisostere']}/search/"
            data = {
                "smiles": smiles,
                "min_occurrence": min_occurrence,
            }

            response = self.http.make_request(submit_url, method="POST", data=data)

            if not response:
                return {}

            # Get job ID and poll for results
            soup = BeautifulSoup(response.text, "html.parser")
            job_id = soup.find("input", {"name": "job"})["value"]

            results_url = f"{self.BASE_URLS['bioisostere']}/results/{job_id}"
            max_attempts = 30
            progress = tqdm(
                total=max_attempts,
                desc="Searching bioisosteres",
                unit="attempts",
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
                logger.warning("Bioisostere search timed out")
                return {}

            # Parse results
            soup = BeautifulSoup(response.text, "html.parser")
            results = {
                "replacements": [],
                "url": f"{self.BASE_URLS['bioisostere']}/result.php?job={job_id}",
            }

            replacement_divs = soup.find_all("div", {"class": "replacement"})
            if replacement_divs:
                parse_progress = tqdm(
                    replacement_divs, desc="Parsing bioisosteres", unit="replacements"
                )

                for div in parse_progress:
                    replacement = {
                        "smiles": div.get("data-smiles"),
                        "occurrences": int(div.get("data-occurrences", 0)),
                        "similarity": float(div.get("data-similarity", 0)),
                        "activity_delta": float(div.get("data-activity-delta", 0)),
                        "properties": {},
                        "references": [],
                    }

                    # Extract property changes
                    props_div = div.find("div", {"class": "properties"})
                    if props_div:
                        for prop in props_div.find_all("div"):
                            replacement["properties"][prop.get("data-property")] = {
                                "original": float(prop.get("data-original", 0)),
                                "modified": float(prop.get("data-modified", 0)),
                                "delta": float(prop.get("data-delta", 0)),
                            }

                    # Extract references
                    refs_div = div.find("div", {"class": "references"})
                    if refs_div:
                        for ref in refs_div.find_all("div"):
                            replacement["references"].append(
                                {
                                    "doi": ref.get("data-doi"),
                                    "title": ref.get("data-title"),
                                    "year": int(ref.get("data-year", 0)),
                                }
                            )

                    results["replacements"].append(replacement)

                parse_progress.close()

            # Cache result
            if use_cache:
                self._cache_result(cache_key, results)

            return results

        except Exception as e:
            logger.error(f"Error getting bioisosteres: {str(e)}")
            return {}
