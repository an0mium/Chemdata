"""BindingDB data source integration.

This module provides comprehensive functionality for:
1. BindingDB data downloading and processing
2. Target and activity classification
3. Structure validation and standardization
4. Property calculation and enrichment
5. Integration with other data sources
6. ML predictions and analysis

The module combines efficient data loading and caching with
comprehensive target pattern matching and activity classification.
"""

import logging
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Any
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict

import pandas as pd
import numpy as np
import requests
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from ..models.compound import Compound, CompoundType, TargetData
from ..web_enrichment.data_sources.swiss import SwissClient
from ..web_enrichment.data_sources.community import CommunityClient
from ..web_enrichment.data_sources.social import SocialDataHarvester
from ..processors.structure.ml.predictors import (
    ActivityPredictor,
    ToxicityPredictor,
    AbusePotentialPredictor,
)


class BindingDBSource:
    """BindingDB data source handler with comprehensive processing."""

    # BindingDB download URLs
    BINDINGDB_URLS = {
        "all": "https://bindingdb.org/bind/downloads/BindingDB_All.tsv",
        "psychoactive": "https://bindingdb.org/bind/downloads/BindingDB_Psychoactive.tsv",
    }

    # Target patterns for filtering compounds
    TARGET_PATTERNS = {
        # Serotonin receptors and transporters
        "serotonin": (
            r"5-HT\d*[A-Z]?|serotonin|SLC6A4|HTR[12][A-Z]|"
            r"5-hydroxytryptamine|tryptamine"
        ),
        # Dopamine receptors and transporters
        "dopamine": (
            r"D\d+|dopamine|DAT|SLC6A3|DRD[1-5]|" r"dopaminergic|catecholamine"
        ),
        # Norepinephrine system
        "norepinephrine": (
            r"norepinephrine|NET|SLC6A2|ADRA\d[A-Z]|" r"adrenergic|noradrenergic"
        ),
        # GABA system
        "gaba": (
            r"GABA[A-Z]?\d*|SLC6A1|GABR[A-Z]\d|" r"gamma-aminobutyric|benzodiazepine"
        ),
        # Glutamate system
        "glutamate": (
            r"glutamate|NMDA|AMPA|mGluR|GluR|" r"GRM\d|GRIN[12][A-Z]|kainate"
        ),
        # Opioid system
        "opioid": (r"[μμκδ]?-?opioid|MOR|KOR|DOR|OPRM1|" r"OPRK1|OPRD1|endorphin"),
        # Cannabinoid system
        "cannabinoid": r"cannabinoid|CB[12]|CNR[12]|endocannabinoid",
        # Psychedelic-related
        "psychedelic": (
            r"psychedelic|hallucinogen|entheogen|"
            r"5-HT2A|HTR2A|DMT|LSD|psilocybin|"
            r"mescaline|ayahuasca|ibogaine"
        ),
        # Nootropic-related
        "nootropic": (
            r"nootropic|cognitive enhancer|smart drug|"
            r"acetylcholine|nicotinic|muscarinic|CHRN[A-Z]\d|"
            r"racetam|ampakine|eugeroic"
        ),
        # Dissociative-related
        "dissociative": (
            r"dissociative|NMDA|ketamine|PCP|"
            r"GRIN[12][A-Z]|glutamate|arylcyclohexylamine"
        ),
        # Stimulant-related
        "stimulant": (
            r"stimulant|amphetamine|cocaine|methylphenidate|"
            r"DAT|NET|SERT|monoamine|cathinone"
        ),
    }

    # Activity type patterns
    ACTIVITY_PATTERNS = {
        "Ki": r"(?:K_?i|inhibition constant|binding constant)",
        "IC50": r"IC_?50|inhibitory concentration|half maximal inhibition",
        "EC50": r"EC_?50|effective concentration|half maximal effect",
        "Kd": r"K_?d|dissociation constant|binding affinity",
        "potency": r"potency|activity|binding|affinity|efficacy",
        "ED50": r"ED_?50|effective dose|half maximal dose",
        "Emax": r"E_?max|maximal effect|maximal response",
        "Hill": r"Hill|slope|cooperativity|n_H",
    }

    def __init__(
        self,
        data_dir: Optional[Path] = None,
        model_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        n_workers: int = 4,
        batch_size: int = 1000,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize BindingDB source.
        
        Args:
            data_dir: Optional directory for data files
            model_dir: Optional directory for ML models
            cache_dir: Optional directory for caching
            n_workers: Number of worker threads
            batch_size: Batch size for processing
            logger: Optional logger instance
        """
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        self.data_dir = data_dir
        self.model_dir = model_dir
        self.cache_dir = cache_dir
        self.n_workers = n_workers
        self.batch_size = batch_size

        # Create directories
        for d in [data_dir, model_dir, cache_dir]:
            if d:
                d.mkdir(parents=True, exist_ok=True)

        # Initialize standardizer
        self.standardizer = rdMolStandardize.Standardizer()

        # Initialize ML models if model_dir provided
        if model_dir:
            self.logger.info("Loading ML models...")
            self.activity_predictor = ActivityPredictor(
                model_dir=model_dir / "activity"
            )
            self.toxicity_predictor = ToxicityPredictor(
                model_dir=model_dir / "toxicity"
            )
            self.abuse_predictor = AbusePotentialPredictor(
                model_dir=model_dir / "abuse"
            )

        # Initialize thread pools
        self.prediction_executor = ThreadPoolExecutor(
            max_workers=n_workers,
            thread_name_prefix="prediction",
        )
        self.web_executor = ThreadPoolExecutor(
            max_workers=n_workers * 2,
            thread_name_prefix="web",
        )

    def download_data(
        self,
        dataset: str = "all",
        force: bool = False,
    ) -> Path:
        """Download BindingDB dataset.
        
        Args:
            dataset: Which dataset to download ("all" or "psychoactive")
            force: Whether to force download even if cached
            
        Returns:
            Path to downloaded file
        """
        if dataset not in self.BINDINGDB_URLS:
            raise ValueError(f"Unknown dataset: {dataset}")

        url = self.BINDINGDB_URLS[dataset]
        filename = f"bindingdb_{dataset}.tsv"
        
        if self.cache_dir:
            cache_file = self.cache_dir / filename
            if not force and cache_file.exists():
                self.logger.info(f"Using cached file: {cache_file}")
                return cache_file
            
            self.logger.info(f"Downloading {url} to {cache_file}")
            self._download_file(url, cache_file)
            return cache_file
        
        # No cache dir, download to temp file
        import tempfile
        tmp_file = Path(tempfile.mktemp(suffix=".tsv"))
        self.logger.info(f"Downloading {url} to {tmp_file}")
        self._download_file(url, tmp_file)
        return tmp_file

    def _download_file(self, url: str, output_file: Path) -> None:
        """Download file with progress bar."""
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        total_size = int(response.headers.get("content-length", 0))
        block_size = 8192
        
        with output_file.open("wb") as f, tqdm(
            total=total_size,
            unit="B",
            unit_scale=True,
            desc=output_file.name,
        ) as pbar:
            for chunk in response.iter_content(block_size):
                f.write(chunk)
                pbar.update(len(chunk))

    def load_compounds(
        self,
        input_file: Path,
        target_patterns: Optional[Dict[str, str]] = None,
        skip_predictions: bool = False,
        skip_web_data: bool = False,
        use_cache: bool = True,
    ) -> List[Compound]:
        """Load and process compounds from BindingDB file.
        
        Args:
            input_file: Path to BindingDB TSV file
            target_patterns: Optional patterns to filter targets
            skip_predictions: Whether to skip ML predictions
            skip_web_data: Whether to skip web data enrichment
            use_cache: Whether to use cached results
            
        Returns:
            List of processed CompoundData instances
        """
        self.logger.info(f"Loading compounds from {input_file}")
        
        # Use provided patterns or defaults
        patterns = target_patterns or self.TARGET_PATTERNS
        
        # Read TSV file in chunks
        compounds = {}
        chunks = pd.read_csv(
            input_file,
            sep="\t",
            chunksize=self.batch_size,
            usecols=[
                "Ligand Name",
                "Ligand SMILES",
                "Ligand InChI",
                "Ligand InChI Key",
                "Target Name",
                "Target Source Organism",
                "Ki (nM)",
                "IC50 (nM)",
                "Kd (nM)",
                "EC50 (nM)",
                "DOI",
                "PubMed ID",
                "Patent Number",
                "Authors",
                "Institution",
                "Article Title",
                "Journal",
                "Publication Year",
                "Measurement Type",
                "Assay Description",
                "Assay Type",
                "pH",
                "Temperature",
                "Curation Level",
            ],
        )
        
        for chunk in tqdm(chunks, desc="Processing chunks"):
            # Process chunk in parallel
            with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
                futures = []
                for _, row in chunk.iterrows():
                    futures.append(
                        executor.submit(
                            self._process_compound,
                            row,
                            patterns,
                            skip_predictions,
                            skip_web_data,
                            use_cache,
                        )
                    )
                
                # Collect results
                for future in as_completed(futures):
                    try:
                        result = future.result()
                        if result:
                            name, compound = result
                            if name not in compounds:
                                compounds[name] = compound
                            else:
                                compounds[name].merge(compound)
                    except Exception as e:
                        self.logger.error(f"Error processing compound: {str(e)}")
        
        self.logger.info(f"Loaded {len(compounds)} compounds")
        return list(compounds.values())

    def _process_compound(
        self,
        row: pd.Series,
        target_patterns: Dict[str, str],
        skip_predictions: bool,
        skip_web_data: bool,
        use_cache: bool,
    ) -> Optional[tuple[str, Compound]]:
        """Process a single compound row.
        
        Args:
            row: DataFrame row containing compound data
            target_patterns: Patterns to filter targets
            skip_predictions: Whether to skip ML predictions
            skip_web_data: Whether to skip web data enrichment
            use_cache: Whether to use cached results
            
        Returns:
            Tuple of (compound name, CompoundData) or None if invalid
        """
        try:
            name = row["Ligand Name"]
            smiles = row["Ligand SMILES"]
            
            if not (name and smiles):
                return None
            
            # Validate and standardize structure
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None
                
            mol = self.standardizer.standardize(mol)
            standard_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            
            # Calculate properties
            mw = Descriptors.ExactMolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            rotatable = Descriptors.NumRotatableBonds(mol)
            
            # Create compound
            compound = Compound(
                name=name,
                smiles=standard_smiles,
                inchi=row["Ligand InChI"],
                inchi_key=row["Ligand InChI Key"],
                compound_type=CompoundType.OTHER,
                molecular_weight=mw,
                logp=logp,
                tpsa=tpsa,
                hbd=hbd,
                hba=hba,
                rotatable_bonds=rotatable,
            )
            
            # Add target data
            target = TargetData(
                common_name=row["Target Name"],
                organism=row["Target Source Organism"],
            )
            
            # Add binding data
            if pd.notna(row["Ki (nM)"]):
                target.affinity_value = row["Ki (nM)"]
                target.affinity_type = "Ki"
            elif pd.notna(row["IC50 (nM)"]):
                target.affinity_value = row["IC50 (nM)"]
                target.affinity_type = "IC50"
            elif pd.notna(row["Kd (nM)"]):
                target.affinity_value = row["Kd (nM)"]
                target.affinity_type = "Kd"
            elif pd.notna(row["EC50 (nM)"]):
                target.affinity_value = row["EC50 (nM)"]
                target.affinity_type = "EC50"
            
            # Add confidence score
            target.confidence = self._calculate_confidence(row)
            
            # Add references
            if pd.notna(row["DOI"]):
                compound.reference_dois.add(row["DOI"])
            if pd.notna(row["PubMed ID"]):
                compound.reference_pmids.add(str(row["PubMed ID"]))
            if pd.notna(row["Patent Number"]):
                if "patents" not in compound.patent_data:
                    compound.patent_data["patents"] = []
                compound.patent_data["patents"].append({
                    "number": row["Patent Number"],
                    "authors": row["Authors"] if pd.notna(row["Authors"]) else None,
                    "institution": row["Institution"] if pd.notna(row["Institution"]) else None,
                    "title": row["Article Title"] if pd.notna(row["Article Title"]) else None,
                    "journal": row["Journal"] if pd.notna(row["Journal"]) else None,
                    "year": int(row["Publication Year"]) if pd.notna(row["Publication Year"]) else None,
                })
            
            # Filter by target patterns
            if target_patterns:
                for pattern_name, pattern in target_patterns.items():
                    if pd.notna(row["Target Name"]) and re.search(pattern, row["Target Name"], re.I):
                        target.is_primary = True
                        compound.compound_type = CompoundType.PSYCHOACTIVE
                        break
            
            compound.targets.append(target)
            
            # Add ML predictions
            if not skip_predictions and self.model_dir:
                try:
                    # Activity predictions
                    activity_pred = self.activity_predictor.predict(standard_smiles)
                    if activity_pred:
                        compound.activity_predictions = activity_pred

                    # Toxicity predictions
                    toxicity_pred = self.toxicity_predictor.predict(standard_smiles)
                    if toxicity_pred:
                        compound.toxicity_predictions = toxicity_pred

                    # Abuse potential predictions
                    abuse_pred = self.abuse_predictor.predict(standard_smiles)
                    if abuse_pred:
                        compound.abuse_predictions = abuse_pred

                except Exception as e:
                    self.logger.error(f"Error getting predictions: {str(e)}")

            # Add web data
            if not skip_web_data and self.cache_dir:
                try:
                    # Get Swiss data
                    swiss_client = SwissClient(
                        http_client=None,
                        model_dir=str(self.model_dir),
                        cache_dir=str(self.cache_dir),
                    )
                    swiss_client.process_compounds(
                        [compound],
                        skip_predictions=skip_predictions,
                        use_cache=use_cache,
                    )

                    # Get community data
                    community_client = CommunityClient(
                        http_client=None,
                        model_dir=str(self.model_dir),
                        cache_dir=str(self.cache_dir),
                    )
                    community_data = community_client.get_compound_data(
                        name,
                        cas_number=compound.cas_number,
                        use_cache=use_cache,
                    )
                    if community_data:
                        compound.community_data = community_data

                    # Get social data
                    social_harvester = SocialDataHarvester(
                        model_dir=str(self.model_dir),
                        cache_dir=str(self.cache_dir),
                    )
                    social_data = social_harvester.enrich_compound_data(
                        compound,
                        use_cache=use_cache,
                    )
                    if social_data:
                        compound.social_data = social_data

                except Exception as e:
                    self.logger.error(f"Error getting web data: {str(e)}")
            
            return name, compound
            
        except Exception as e:
            self.logger.error(f"Error processing compound row: {str(e)}")
            return None

    def _calculate_confidence(self, row: pd.Series) -> float:
        """Calculate confidence score for binding data.

        Args:
            row: DataFrame row containing compound data

        Returns:
            Confidence score between 0 and 1
        """
        score = 0.0
        total = 0.0

        # Check for key fields
        if row.get("Measurement Value") is not None:
            score += 1
            total += 1
        if row.get("Measurement Unit"):
            score += 1
            total += 1
        if row.get("Measurement Type"):
            score += 1
            total += 1
        if row.get("Assay Description"):
            score += 1
            total += 1
        if row.get("Target Name"):
            score += 1
            total += 1
        if row.get("Target Source Organism"):
            score += 1
            total += 1
        if row.get("Target Sequence"):
            score += 1
            total += 1
        if row.get("Article DOI"):
            score += 2  # Weight references more heavily
            total += 2

        # Check data quality
        if row.get("pH") is not None:
            score += 0.5
            total += 0.5
        if row.get("Temperature") is not None:
            score += 0.5
            total += 0.5
        if row.get("Curation Level", "").lower() == "expert":
            score += 1
            total += 1

        return score / total if total > 0 else 0.0
