"""Community data source functionality.

This module handles:
1. Community data source discovery through search engines
2. Data extraction from discovered sources
3. Validation and filtering of search results
4. Caching of search results and extracted data
5. Error handling and retry logic
6. Progress tracking
7. Model ensembling and confidence scoring
8. Parallel processing and error recovery
"""

import logging
import os
import re
import json
from pathlib import Path
from typing import Dict, List, Optional, Any, Set, Tuple
from datetime import datetime, timedelta
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict

from bs4 import BeautifulSoup
import requests
from tqdm import tqdm
import torch
from transformers import pipeline, AutoTokenizer, AutoModelForTokenClassification
from sentence_transformers import SentenceTransformer, util

from logger import LogManager
from ..http_client import HttpClient
from ..llm_utils import analyze_content_with_llm
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
from binding_data_processor.processors.structure.ml.predictors.psychoactive import (
    PsychoactivePredictor,
)
from binding_data_processor.processors.structure.ml.predictors.nootropic import (
    NootropicPredictor,
)
from binding_data_processor.processors.structure.ml.models.ensemble import (
    EnsembleModel,
)

logger = LogManager().get_logger("web_enrichment.data_sources.community")


class CommunityClient:
    """Client for discovering and extracting data from community sources."""

    # Cache settings
    CACHE_DIR = Path(".cache/community")
    CACHE_DURATIONS = {
        "search": timedelta(days=1),  # Refresh search results daily
        "content": timedelta(days=7),  # Cache content for a week
        "predictions": timedelta(days=30),  # Cache predictions for a month
    }

    # Search settings
    SEARCH_DOMAINS = {
        # Community wikis
        "psychonaut": "psychonautwiki.org",
        "erowid": "erowid.org",
        "tripsit": "drugs.tripsit.me",
        "bluelight": "bluelight.org",
        "drugs-forum": "drugs-forum.com",
        "dmt-nexus": "dmt-nexus.me",
        "lucid": "lucid.me",
        "shroomery": "shroomery.org",
        # Reddit communities
        "reddit_rc": "reddit.com/r/researchchemicals",
        "reddit_drugs": "reddit.com/r/Drugs",
        "reddit_nootropics": "reddit.com/r/Nootropics",
        "reddit_drugnerds": "reddit.com/r/DrugNerds",
        "reddit_psychonaut": "reddit.com/r/Psychonaut",
        "reddit_neuro": "reddit.com/r/neuroscience",
        "reddit_pharm": "reddit.com/r/pharmacology",
        "reddit_medchem": "reddit.com/r/MedicinalChemistry",
        "reddit_theehive": "reddit.com/r/TheeHive",
        "reddit_babybees": "reddit.com/r/BabyBees",
        "reddit_obscuredrugs": "reddit.com/r/ObscureDrugs",
        "reddit_peptides": "reddit.com/r/Peptides",
        "reddit_stackadvice": "reddit.com/r/StackAdvice",
        # Academic sources
        "pubmed": "pubmed.ncbi.nlm.nih.gov",
        "chemrxiv": "chemrxiv.org",
        "biorxiv": "biorxiv.org",
        "sciencedirect": "sciencedirect.com",
        "nature": "nature.com",
        "frontiers": "frontiersin.org",
        "acs": "pubs.acs.org",
        "wiley": "onlinelibrary.wiley.com",
        # Patent databases
        "patents_google": "patents.google.com",
        "patents_uspto": "patft.uspto.gov",
        "patents_espacenet": "worldwide.espacenet.com",
        "patents_wipo": "patentscope.wipo.int",
        # Regulatory databases
        "fda": "fda.gov",
        "ema": "ema.europa.eu",
        "who": "who.int",
        "inchem": "inchem.org",
        "toxnet": "toxnet.nlm.nih.gov",
        "chemidplus": "chem.nlm.nih.gov",
        "pubchem": "pubchem.ncbi.nlm.nih.gov",
        "chemspider": "chemspider.com",
        # Additional sources
        "wikipedia": "wikipedia.org",
        "drugbank": "drugbank.ca",
        "chembl": "ebi.ac.uk/chembl",
        "bindingdb": "bindingdb.org",
        "zinc": "zinc.docking.org",
        "emcdda": "emcdda.europa.eu",
        "unodc": "unodc.org",
        "nida": "nida.nih.gov",
    }

    # Content validation patterns
    VALIDATION_PATTERNS = {
        "chemistry": [
            r"chemical",
            r"compound",
            r"molecule",
            r"substance",
            r"drug",
            r"synthesis",
            r"structure",
            r"formula",
            r"isomer",
            r"analog",
            r"derivative",
            r"salt",
            r"freebase",
            r"crystal",
            r"powder",
            r"solution",
            r"preparation",
            r"purity",
            r"stability",
            r"degradation",
        ],
        "pharmacology": [
            r"pharmacology",
            r"receptor",
            r"binding",
            r"mechanism",
            r"activity",
            r"agonist",
            r"antagonist",
            r"affinity",
            r"potency",
            r"selectivity",
            r"pharmacokinetics?",
            r"metabolism",
            r"half-life",
            r"bioavailability",
            r"distribution",
            r"absorption",
            r"excretion",
            r"enzyme",
            r"transporter",
            r"signaling",
            r"pathway",
        ],
        "effects": [
            r"effects?",
            r"action",
            r"psychoactive",
            r"hallucinogenic?",
            r"psychedelic",
            r"nootropic",
            r"cognitive",
            r"behavioral",
            r"subjective",
            r"experience",
            r"report",
            r"duration",
            r"onset",
            r"dosage",
            r"therapeutic",
            r"clinical",
            r"treatment",
            r"indication",
            r"efficacy",
            r"response",
        ],
        "safety": [
            r"safety",
            r"toxicity",
            r"risk",
            r"danger",
            r"adverse",
            r"contraindication",
            r"interaction",
            r"overdose",
            r"tolerance",
            r"dependence",
            r"addiction",
            r"withdrawal",
            r"lethal",
            r"carcinogenic",
            r"mutagenic",
            r"teratogenic",
            r"neurotoxic",
            r"hepatotoxic",
            r"nephrotoxic",
            r"cardiotoxic",
        ],
        "regulation": [
            r"regulation",
            r"legal",
            r"status",
            r"schedule",
            r"controlled",
            r"approved",
            r"banned",
            r"restricted",
            r"clinical",
            r"trial",
            r"patent",
            r"registration",
            r"authorization",
            r"marketing",
            r"prescription",
            r"otc",
            r"research",
            r"development",
            r"phase",
            r"investigation",
        ],
    }

    def __init__(
        self,
        http_client: HttpClient,
        model_dir: str = "../models",
        cache_dir: Optional[str] = None,
        n_workers: int = 4,
        max_retries: int = 3,
        device: str = "cuda" if torch.cuda.is_available() else "cpu",
    ):
        """Initialize community client.

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
        for cache_type in self.CACHE_DURATIONS.keys():
            (self.CACHE_DIR / cache_type).mkdir(parents=True, exist_ok=True)

        # Initialize ML models
        logger.info("Loading ML models...")

        # Base predictors
        self.activity_predictor = ActivityPredictor(model_dir=f"{model_dir}/activity")
        self.toxicity_predictor = ToxicityPredictor(model_dir=f"{model_dir}/toxicity")
        self.abuse_predictor = AbusePotentialPredictor(model_dir=f"{model_dir}/abuse")
        self.psychoactive_predictor = PsychoactivePredictor(
            model_dir=f"{model_dir}/psychoactive"
        )
        self.nootropic_predictor = NootropicPredictor(
            model_dir=f"{model_dir}/nootropic"
        )

        # Ensemble models
        self.activity_ensemble = EnsembleModel(
            [
                self.activity_predictor,
                self.psychoactive_predictor,
                self.nootropic_predictor,
            ],
            weights=[0.4, 0.3, 0.3],
        )
        self.safety_ensemble = EnsembleModel(
            [self.toxicity_predictor, self.abuse_predictor],
            weights=[0.5, 0.5],
        )

        # NER models
        self.ner_tokenizer = AutoTokenizer.from_pretrained(
            "allenai/scibert_scivocab_uncased"
        )
        self.ner_model = AutoModelForTokenClassification.from_pretrained(
            f"{model_dir}/ner/chemical"
        ).to(device)
        self.ner_bio = AutoModelForTokenClassification.from_pretrained(
            f"{model_dir}/ner/bio"
        ).to(device)
        self.ner_drug = AutoModelForTokenClassification.from_pretrained(
            f"{model_dir}/ner/drug"
        ).to(device)

        # Classification models
        self.text_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/community/compound_classifier",
            device=device,
        )
        self.safety_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/community/safety_classifier",
            device=device,
        )
        self.mechanism_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/community/mechanism_classifier",
            device=device,
        )
        self.effect_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/community/effect_classifier",
            device=device,
        )
        self.relationship_classifier = pipeline(
            "text-classification",
            model=f"{model_dir}/community/relationship_classifier",
            device=device,
        )

        # Similarity models
        self.similarity_model = SentenceTransformer(
            "sentence-transformers/all-MiniLM-L6-v2"
        ).to(device)

    def _get_cached_result(
        self, cache_key: str, cache_type: str
    ) -> Optional[Dict[str, Any]]:
        """Get cached result if available and not expired.

        Args:
            cache_key: Cache key to look up
            cache_type: Type of cached data

        Returns:
            Cached result or None if not found/expired
        """
        cache_file = self.CACHE_DIR / cache_type / f"{cache_key}.json"
        if cache_file.exists():
            try:
                data = json.loads(cache_file.read_text())
                cached_time = datetime.fromisoformat(data["cached_at"])
                if datetime.now() - cached_time < self.CACHE_DURATIONS.get(
                    cache_type, timedelta(days=7)
                ):
                    return data["result"]
            except Exception as e:
                logger.error(f"Error reading cache file {cache_file}: {e}")
        return None

    def _cache_result(
        self, cache_key: str, result: Dict[str, Any], cache_type: str
    ) -> None:
        """Cache result with timestamp.

        Args:
            cache_key: Cache key to store under
            result: Result data to cache
            cache_type: Type of data to cache
        """
        try:
            cache_file = self.CACHE_DIR / cache_type / f"{cache_key}.json"
            data = {
                "cached_at": datetime.now().isoformat(),
                "result": result,
            }
            cache_file.write_text(json.dumps(data, indent=2))
        except Exception as e:
            logger.error(f"Error writing cache file {cache_file}: {e}")

    def get_compound_data(
        self,
        name: str,
        cas_number: Optional[str] = None,
        smiles: Optional[str] = None,
        llm_api_key: Optional[str] = None,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Get compound data from community sources.

        Args:
            name: Compound name
            cas_number: Optional CAS number
            smiles: Optional SMILES structure
            llm_api_key: Optional API key for LLM analysis
            use_cache: Whether to use cached results

        Returns:
            Dictionary containing compound data or None
        """
        # Check cache first
        if use_cache:
            cache_key = f"compound_{name}_{cas_number}_{smiles}"
            cached = self._get_cached_result(cache_key, "content")
            if cached:
                return cached

        data = {
            "chemistry": [],
            "pharmacology": [],
            "effects": [],
            "toxicity": [],
            "regulation": [],
            "sources": [],
            "urls": {},
            "predictions": {},
            "relationships": {},
        }

        # Build search queries
        queries = [
            f"{name} drug",
            f"{name} pharmacology",
            f"{name} effects",
            f"{name} toxicity",
            f"{name} legal status",
            f"{name} mechanism of action",
            f"{name} binding affinity",
            f"{name} structure activity",
            f"{name} clinical trials",
            f"{name} synthesis route",
        ]
        if cas_number:
            queries.append(f"{cas_number} chemical")
        if smiles:
            queries.append(f"{smiles} structure")

        # Search each domain in parallel
        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            futures = {}
            for domain_name, domain in self.SEARCH_DOMAINS.items():
                futures[domain_name] = []
                for query in queries:
                    futures[domain_name].append(
                        executor.submit(
                            self._search_domain,
                            domain,
                            query,
                            name,
                            use_cache,
                        )
                    )

            # Collect results with progress bar
            with tqdm(
                total=len(futures) * len(queries),
                desc="Searching sources",
            ) as pbar:
                for domain_name, domain_futures in futures.items():
                    domain_data = {
                        "chemistry": [],
                        "pharmacology": [],
                        "effects": [],
                        "toxicity": [],
                        "regulation": [],
                        "urls": set(),
                    }

                    for future in domain_futures:
                        try:
                            result = future.result()
                            if result:
                                # Merge data
                                for key in [
                                    "chemistry",
                                    "pharmacology",
                                    "effects",
                                    "toxicity",
                                    "regulation",
                                ]:
                                    if key in result:
                                        domain_data[key].extend(result[key])
                                if "urls" in result:
                                    domain_data["urls"].update(result["urls"])
                                if "relationships" in result:
                                    self._merge_relationships(
                                        data["relationships"],
                                        result["relationships"],
                                    )
                        except Exception as e:
                            logger.error(f"Error processing {domain_name} results: {e}")
                        pbar.update(1)

                    # Add unique data to results
                    if any(v for v in domain_data.values()):
                        for key in [
                            "chemistry",
                            "pharmacology",
                            "effects",
                            "toxicity",
                            "regulation",
                        ]:
                            data[key].extend(
                                [x for x in domain_data[key] if x not in data[key]]
                            )
                        data["urls"][domain_name] = list(domain_data["urls"])
                        if domain_name not in data["sources"]:
                            data["sources"].append(domain_name)

        # Add ML predictions
        try:
            if name:
                # Activity ensemble predictions
                activity_pred = self.activity_ensemble.predict(name)
                if activity_pred:
                    data["predictions"]["activity"] = activity_pred

                # Safety ensemble predictions
                safety_pred = self.safety_ensemble.predict(name)
                if safety_pred:
                    data["predictions"]["safety"] = safety_pred

                # Individual predictions
                for pred_type, predictor in [
                    ("toxicity", self.toxicity_predictor),
                    ("abuse", self.abuse_predictor),
                    ("psychoactive", self.psychoactive_predictor),
                    ("nootropic", self.nootropic_predictor),
                ]:
                    try:
                        pred = predictor.predict(name)
                        if pred:
                            data["predictions"][pred_type] = pred
                    except Exception as e:
                        logger.error(f"Error getting {pred_type} prediction: {e}")

                # Mechanism predictions
                if data["pharmacology"]:
                    mechanism_pred = self.mechanism_classifier(
                        " ".join(data["pharmacology"])
                    )[0]
                    if mechanism_pred["score"] > 0.8:
                        data["predictions"]["mechanism"] = mechanism_pred

                # Effect predictions
                if data["effects"]:
                    effect_pred = self.effect_classifier(" ".join(data["effects"]))[0]
                    if effect_pred["score"] > 0.8:
                        data["predictions"]["effects"] = effect_pred

                # Relationship predictions
                if data["chemistry"]:
                    relationship_pred = self.relationship_classifier(
                        " ".join(data["chemistry"])
                    )[0]
                    if relationship_pred["score"] > 0.8:
                        data["predictions"]["relationships"] = relationship_pred

        except Exception as e:
            logger.error(f"Error getting predictions: {e}")

        # Use LLM to analyze content if available
        if llm_api_key and data["sources"]:
            try:
                llm_analysis = analyze_content_with_llm(
                    {
                        "name": name,
                        "cas": cas_number,
                        "smiles": smiles,
                        "data": data,
                    },
                    llm_api_key,
                )
                if llm_analysis:
                    data["llm_analysis"] = llm_analysis
            except Exception as e:
                logger.error(f"Error analyzing with LLM: {e}")

        # Cache result if we found anything
        if use_cache and data["sources"]:
            self._cache_result(cache_key, data, "content")

        return data if data["sources"] else None

    def _search_domain(
        self,
        domain: str,
        query: str,
        compound_name: str,
        use_cache: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """Search a domain for compound information.

        Args:
            domain: Domain to search
            query: Search query
            compound_name: Name of compound
            use_cache: Whether to use cached results

        Returns:
            Dictionary containing found data or None
        """
        # Check cache first
        if use_cache:
            cache_key = f"search_{domain}_{query}"
            cached = self._get_cached_result(cache_key, "search")
            if cached:
                return cached

        data = {
            "chemistry": [],
            "pharmacology": [],
            "effects": [],
            "toxicity": [],
            "regulation": [],
            "urls": set(),
            "relationships": {},
        }

        try:
            # Use web search to find pages
            search_query = f'site:{domain} "{query}"'
            search_results = self.http.search(search_query, max_results=10)

            # Process each result in parallel
            with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
                futures = []
                for url in search_results:
                    futures.append(
                        executor.submit(
                            self._process_search_result,
                            url,
                            compound_name,
                        )
                    )

                # Collect results
                for future in as_completed(futures):
                    try:
                        result = future.result()
                        if result:
                            # Merge data
                            for key in [
                                "chemistry",
                                "pharmacology",
                                "effects",
                                "toxicity",
                                "regulation",
                            ]:
                                if key in result:
                                    data[key].extend(result[key])
                            if "url" in result:
                                data["urls"].add(result["url"])
                            if "relationships" in result:
                                self._merge_relationships(
                                    data["relationships"],
                                    result["relationships"],
                                )
                    except Exception as e:
                        logger.error(f"Error processing search result: {e}")

            # Cache result if we found anything
            if use_cache and any(v for v in data.values()):
                self._cache_result(cache_key, data, "search")

            return data if any(v for v in data.values()) else None

        except Exception as e:
            logger.error(f"Error searching {domain}: {e}")
            return None

    def _process_search_result(
        self, url: str, compound_name: str
    ) -> Optional[Dict[str, Any]]:
        """Process a search result URL.

        Args:
            url: URL to process
            compound_name: Name of compound to look for

        Returns:
            Dictionary containing extracted data or None
        """
        try:
            # Get page content
            response = self.http.get(url)
            if not response:
                return None

            # Parse HTML
            soup = BeautifulSoup(response, "html.parser")

            # Remove unwanted elements
            for element in soup(["script", "style"]):
                element.decompose()

            # Get text content
            text = soup.get_text()
            text = " ".join(line.strip() for line in text.splitlines() if line.strip())

            # Validate content
            if not self._validate_content(text, compound_name):
                return None

            # Extract data
            data = {
                "chemistry": [],
                "pharmacology": [],
                "effects": [],
                "toxicity": [],
                "regulation": [],
                "url": url,
                "relationships": {},
            }

            # Extract sections using patterns
            for section, patterns in self.VALIDATION_PATTERNS.items():
                for pattern in patterns:
                    matches = re.finditer(f"[^.]*{pattern}[^.]*\.", text, re.I)
                    for match in matches:
                        sentence = match.group(0).strip()
                        if compound_name.lower() in sentence.lower():
                            data[section].append(sentence)

            # Use ML to classify content
            try:
                # Split into chunks to avoid token limits
                chunks = [text[i : i + 512] for i in range(0, len(text), 512)]
                for chunk in chunks:
                    # Classify relevance
                    classification = self.text_classifier(chunk)[0]
                    if (
                        classification["label"] == "relevant"
                        and classification["score"] > 0.8
                    ):
                        # Extract chemical entities using ensemble
                        entities = []
                        for model in [self.ner_model, self.ner_bio, self.ner_drug]:
                            try:
                                model_entities = model(chunk)
                                entities.extend(model_entities)
                            except Exception as e:
                                logger.error(f"Error running NER model: {e}")

                        # Process entities
                        for entity in entities:
                            if entity["entity"].startswith("B-CHEMICAL"):
                                context = self._get_surrounding_text(
                                    chunk,
                                    entity["start"],
                                    entity["end"],
                                )
                                if context:
                                    data["chemistry"].append(context)

                                    # Extract relationships
                                    self._extract_relationships(
                                        data["relationships"],
                                        compound_name,
                                        entity["word"],
                                        context,
                                    )

                        # Classify safety concerns
                        safety = self.safety_classifier(chunk)[0]
                        if safety["label"] == "concern" and safety["score"] > 0.8:
                            data["toxicity"].append(chunk)

                        # Classify mechanisms
                        mechanism = self.mechanism_classifier(chunk)[0]
                        if mechanism["score"] > 0.8:
                            data["pharmacology"].append(chunk)

                        # Classify effects
                        effect = self.effect_classifier(chunk)[0]
                        if effect["score"] > 0.8:
                            data["effects"].append(chunk)

                        # Classify relationships
                        relationship = self.relationship_classifier(chunk)[0]
                        if relationship["score"] > 0.8:
                            self._extract_relationships(
                                data["relationships"],
                                compound_name,
                                None,  # Extract from context
                                chunk,
                            )

            except Exception as e:
                logger.error(f"Error classifying content: {e}")

            return data if any(v for v in data.values()) else None

        except Exception as e:
            logger.error(f"Error processing URL {url}: {e}")
            return None

    def _validate_content(self, text: str, compound_name: str) -> bool:
        """Validate if content is relevant to compound.

        Args:
            text: Content to validate
            compound_name: Name of compound

        Returns:
            Whether content is relevant
        """
        try:
            # Check content length
            if len(text) < 100:
                return False

            # Look for compound name mentions
            name_pattern = re.compile(re.escape(compound_name), re.I)
            if not name_pattern.search(text):
                return False

            # Count context matches
            context_count = 0
            for patterns in self.VALIDATION_PATTERNS.values():
                for pattern in patterns:
                    if re.search(pattern, text, re.I):
                        context_count += 1

            # Use ML to validate content
            try:
                # Check first chunk
                classification = self.text_classifier(text[:512])[0]
                if (
                    classification["label"] == "relevant"
                    and classification["score"] > 0.8
                ):
                    return True

                # Check for chemical entities
                entities = []
                for model in [self.ner_model, self.ner_bio, self.ner_drug]:
                    try:
                        model_entities = model(text[:512])
                        entities.extend(model_entities)
                    except Exception as e:
                        logger.error(f"Error running NER model: {e}")

                if any(
                    e["entity"].startswith("B-CHEMICAL") and e["score"] > 0.8
                    for e in entities
                ):
                    return True

            except Exception as e:
                logger.error(f"Error classifying content: {e}")

            # Require at least 3 context matches
            return context_count >= 3

        except Exception as e:
            logger.error(f"Error validating content: {e}")
            return False

    def _get_surrounding_text(
        self, text: str, start: int, end: int, context: int = 200
    ) -> Optional[str]:
        """Get text surrounding a match.

        Args:
            text: Full text
            start: Start index of match
            end: End index of match
            context: Number of characters of context

        Returns:
            Text with surrounding context
        """
        try:
            context_start = max(0, start - context)
            context_end = min(len(text), end + context)
            return text[context_start:context_end].strip()
        except Exception as e:
            logger.error(f"Error getting surrounding text: {e}")
            return None

    def _extract_relationships(
        self,
        relationships: Dict[str, Dict[str, Any]],
        compound1: str,
        compound2: Optional[str],
        context: str,
    ) -> None:
        """Extract relationships between compounds.

        Args:
            relationships: Dictionary to store relationships
            compound1: First compound name
            compound2: Optional second compound name (extracted from context if None)
            context: Context containing relationship
        """
        try:
            # Extract compound names if needed
            if compound2 is None:
                # Use NER ensemble to find compounds
                entities = []
                for model in [self.ner_model, self.ner_bio, self.ner_drug]:
                    try:
                        model_entities = model(context)
                        entities.extend(model_entities)
                    except Exception as e:
                        logger.error(f"Error running NER model: {e}")

                # Get high confidence chemical entities
                compounds = [
                    e["word"]
                    for e in entities
                    if e["entity"].startswith("B-CHEMICAL") and e["score"] > 0.8
                ]

                # Remove compound1 and duplicates
                compounds = list(
                    set(c for c in compounds if c.lower() != compound1.lower())
                )

                # Process each found compound
                for compound2 in compounds:
                    self._add_relationship(
                        relationships,
                        compound1,
                        compound2,
                        context,
                    )
            else:
                # Process single compound pair
                self._add_relationship(
                    relationships,
                    compound1,
                    compound2,
                    context,
                )

        except Exception as e:
            logger.error(f"Error extracting relationships: {e}")

    def _add_relationship(
        self,
        relationships: Dict[str, Dict[str, Any]],
        compound1: str,
        compound2: str,
        context: str,
    ) -> None:
        """Add a relationship between two compounds.

        Args:
            relationships: Dictionary to store relationships
            compound1: First compound name
            compound2: Second compound name
            context: Context containing relationship
        """
        try:
            # Calculate text similarity
            similarity = self.similarity_model.encode(
                [compound1, compound2],
                convert_to_tensor=True,
            )
            similarity_score = float(
                util.pytorch_cos_sim(similarity[0], similarity[1])[0][0]
            )

            # Extract relationship type
            relationship_type = "related"
            if re.search(r"analog|derivative|similar", context, re.I):
                relationship_type = "analog"
            elif re.search(r"metabolite|breakdown", context, re.I):
                relationship_type = "metabolite"
            elif re.search(r"precursor|synthetic", context, re.I):
                relationship_type = "precursor"

            # Store relationship
            key = f"{compound1.lower()}-{compound2.lower()}"
            if key not in relationships:
                relationships[key] = {
                    "compound1": compound1,
                    "compound2": compound2,
                    "type": relationship_type,
                    "similarity": similarity_score,
                    "contexts": [],
                }
            relationships[key]["contexts"].append(context)

        except Exception as e:
            logger.error(f"Error adding relationship: {e}")

    def _merge_relationships(
        self,
        relationships1: Dict[str, Dict[str, Any]],
        relationships2: Dict[str, Dict[str, Any]],
    ) -> None:
        """Merge two relationship dictionaries.

        Args:
            relationships1: First relationship dictionary
            relationships2: Second relationship dictionary
        """
        try:
            for key, rel2 in relationships2.items():
                if key in relationships1:
                    # Update existing relationship
                    rel1 = relationships1[key]
                    rel1["contexts"].extend(rel2["contexts"])
                    # Use highest similarity score
                    rel1["similarity"] = max(rel1["similarity"], rel2["similarity"])
                    # Use most specific relationship type
                    if rel2["type"] != "related":
                        rel1["type"] = rel2["type"]
                else:
                    # Add new relationship
                    relationships1[key] = rel2

        except Exception as e:
            logger.error(f"Error merging relationships: {e}")
