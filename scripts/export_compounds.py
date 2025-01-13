#!/usr/bin/env python3
"""Export compounds of interest with comprehensive data.

This script exports:
1. 5-HT2 agonists
2. NMDA antagonists
3. Known psychoactive compounds
4. Specified proteins/peptides
5. Basic biomolecules
6. Emerging abuse threats
7. Longevity compounds
8. Physical enhancement compounds

The output includes:
- TSV format with essential data (CAS numbers, structures, properties)
- JSON format with complete data (including all metadata)
- Excel format for easier viewing
- SDF format for chemical structure software
- MOL format for individual structures
"""

import logging
import json
import csv
import time
import os
import sys
import asyncio
import argparse
import multiprocessing
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Set, Optional, Tuple, Any, Union
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import torch
from tqdm.asyncio import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdDeprotect, SDWriter
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Crippen
from transformers import AutoTokenizer, AutoModelForTokenClassification, pipeline
from sentence_transformers import SentenceTransformer

from binding_data_processor.data_sources.bindingdb_enhanced import BindingDBSourceEnhanced
from binding_data_processor.data_sources.pubchem_enhanced import PubChemClientEnhanced
from binding_data_processor.data_sources.chembl_enhanced import ChEMBLClientEnhanced
from binding_data_processor.web_enrichment.swiss_client_enhanced import SwissClientEnhanced
from binding_data_processor.web_enrichment.community_client_enhanced import CommunityClientEnhanced
from binding_data_processor.web_enrichment.social_client_enhanced import SocialClientEnhanced
from binding_data_processor.web_enrichment.clients.patents import EnhancedPatentClient
from binding_data_processor.data_sources.pubmed_enhanced import PubMedClientEnhanced
from binding_data_processor.web_enrichment.clients.sciencedirect_enhanced import ScienceDirectClientEnhanced
from binding_data_processor.web_enrichment.clients.reddit_enhanced import RedditClientEnhanced
from binding_data_processor.web_enrichment.clients.bluelight_enhanced import BluelightClientEnhanced
from binding_data_processor.models.compound import Compound, CompoundType, TargetData
from binding_data_processor.models.compound.enrichment.validation import validate_compound_data
from binding_data_processor.models.compound.base.types import CompoundData
from binding_data_processor.pipeline.infrastructure.circuit_breaker import CircuitConfig

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", handlers=[logging.StreamHandler(), logging.FileHandler("compound_export.log")]
)
logger = logging.getLogger(__name__)

# Project directories
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
CACHE_DIR = PROJECT_ROOT / "cache"
OUTPUT_DIR = PROJECT_ROOT / "output"
MODEL_DIR = PROJECT_ROOT / "models"

# Ensure directories exist
for directory in [DATA_DIR, CACHE_DIR, OUTPUT_DIR, MODEL_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# Circuit breaker configuration
CIRCUIT_CONFIG = CircuitConfig(failure_threshold=5, recovery_timeout=60, reset_timeout=300, concurrency_limit=10)

# Global variables for models
tokenizer = None
model = None
device = None
ner_model = None
analyzer = None
relevance_model = None
nootropic_model = None

# Global locks for model initialization
_model_init_lock = multiprocessing.Lock()
_models_initialized = multiprocessing.Value("i", 0)


def initialize_models() -> bool:
    """Initialize ML models with proper error handling."""
    global tokenizer, model, device, ner_model, analyzer, relevance_model, nootropic_model

    # Check if models are already initialized
    if _models_initialized.value:
        return True

    # Use lock to ensure thread-safe initialization
    with _model_init_lock:
        try:
            # Double-check initialization after acquiring lock
            if _models_initialized.value:
                return True

            # Set device with proper error handling
            try:
                if torch.backends.mps.is_available():
                    device = torch.device("mps")
                elif torch.cuda.is_available():
                    device = torch.device("cuda")
                else:
                    device = torch.device("cpu")
                logger.info(f"Device set to use {device}")
            except Exception as e:
                logger.warning(f"Error setting device, falling back to CPU: {str(e)}")
                device = torch.device("cpu")

            # Initialize tokenizer first
            tokenizer = AutoTokenizer.from_pretrained(
                "allenai/scibert_scivocab_uncased",
                use_fast=True,
                local_files_only=False,
                model_max_length=512,
                trust_remote_code=False,
                cache_dir=MODEL_DIR / "tokenizer",
            )

            # Initialize model with PyTorch weights
            model = AutoModelForTokenClassification.from_pretrained(
                "allenai/scibert_scivocab_uncased",
                num_labels=2,
                local_files_only=False,
                use_auth_token=None,
                torch_dtype=torch.float32,
                low_cpu_mem_usage=True,
                from_pt=True,  # Use PyTorch weights
                trust_remote_code=False,
                model_type="bert",
                _fast_init=False,
                framework="pt",
                device_map="auto" if torch.cuda.is_available() else None,
                cache_dir=MODEL_DIR / "model",
            )

            # Move model to device if not using CUDA device_map="auto"
            if device.type != "cuda":
                model = model.to(device)
            model.eval()  # Set to evaluation mode

            # Create pipeline with initialized components
            ner_model = pipeline(
                "token-classification",
                model=model,
                tokenizer=tokenizer,
                aggregation_strategy="simple",
                device=device,
            )

            # Initialize sentence transformers with proper settings
            analyzer = SentenceTransformer(
                "allenai/scibert_scivocab_uncased",
                device=device,
                cache_folder=str(MODEL_DIR / "sentence_transformers"),
                from_pt=True,  # Use PyTorch weights
            )

            relevance_model = SentenceTransformer(
                "pritamdeka/S-PubMedBert-MS-MARCO",
                device=device,
                cache_folder=str(MODEL_DIR / "sentence_transformers"),
                from_pt=True,  # Use PyTorch weights
            )

            nootropic_model = SentenceTransformer(
                "microsoft/BiomedNLP-PubMedBert-base-uncased-abstract",
                device=device,
                cache_folder=str(MODEL_DIR / "sentence_transformers"),
                from_pt=True,  # Use PyTorch weights
            )

            # Verify models loaded successfully
            if not all([ner_model, analyzer, relevance_model, nootropic_model]):
                raise RuntimeError("One or more models failed to initialize")

            # Set models to evaluation mode
            analyzer.eval()
            relevance_model.eval()
            nootropic_model.eval()

            # Mark initialization as complete
            _models_initialized.value = 1
            return True

        except Exception as e:
            logger.error(f"Error initializing models: {str(e)}")
            # Reset initialization flag on failure
            _models_initialized.value = 0
            return False


def initialize_models_multi():
    """Initialize ML models with multiprocessing safety and proper error handling."""
    global tokenizer, model, device

    # Check if models are already initialized in this process
    if _models_initialized.value:
        logger.debug("Models already initialized in this process")
        return True

    # Use lock to ensure thread-safe initialization
    with _model_init_lock:
        try:
            # Double-check initialization after acquiring lock
            if _models_initialized.value:
                logger.debug("Models initialized by another thread")
                return True

            logger.info("Initializing scibert model...")

            # Set device with proper error handling
            try:
                if torch.backends.mps.is_available():
                    device = torch.device("mps")
                elif torch.cuda.is_available():
                    device = torch.device("cuda")
                else:
                    device = torch.device("cpu")
                logger.info(f"Device set to use {device}")
            except Exception as e:
                logger.warning(f"Error setting device, falling back to CPU: {str(e)}")
                device = torch.device("cpu")

            # Initialize tokenizer with retry logic
            max_retries = 3
            for attempt in range(max_retries):
                try:
                    tokenizer = AutoTokenizer.from_pretrained(
                        "allenai/scibert_scivocab_uncased",
                        use_fast=True,
                        local_files_only=False,
                        model_max_length=512,
                        trust_remote_code=False,
                        cache_dir=MODEL_DIR / "tokenizer",
                    )
                    break
                except Exception as e:
                    if attempt == max_retries - 1:
                        raise RuntimeError(f"Failed to initialize tokenizer after {max_retries} attempts: {str(e)}")
                    logger.warning(f"Tokenizer initialization attempt {attempt + 1} failed: {str(e)}")
                    time.sleep(1)

            # Initialize model with retry logic and safe settings
            for attempt in range(max_retries):
                try:
                    model = AutoModelForTokenClassification.from_pretrained(
                        "allenai/scibert_scivocab_uncased",
                        num_labels=2,
                        local_files_only=False,
                        use_auth_token=None,
                        torch_dtype=torch.float32,
                        # Multiprocessing safety settings
                        use_cache=False,
                        low_cpu_mem_usage=True,
                        from_pt=True,
                        trust_remote_code=False,
                        model_type="bert",
                        _fast_init=False,
                        framework="pt",
                        device_map="auto" if torch.cuda.is_available() else None,
                        is_parallel=False,
                        offload_folder=None,
                        gradient_checkpointing=False,
                        cache_dir=MODEL_DIR / "model",
                    )
                    break
                except Exception as e:
                    if attempt == max_retries - 1:
                        raise RuntimeError(f"Failed to initialize model after {max_retries} attempts: {str(e)}")
                    logger.warning(f"Model initialization attempt {attempt + 1} failed: {str(e)}")
                    time.sleep(1)

            # Move model to device with error handling
            try:
                if device.type != "cuda":  # Only move if not using CUDA device_map="auto"
                    model = model.to(device)
                model.eval()  # Set to evaluation mode
            except Exception as e:
                logger.error(f"Error moving model to device {device}: {str(e)}")
                raise

            # Verify initialization
            if model is None or tokenizer is None:
                raise RuntimeError("Model initialization failed - components are None")

            # Mark initialization as complete
            _models_initialized.value = 1
            logger.info("Model initialization successful")
            return True

        except Exception as e:
            logger.error(f"Error initializing models: {str(e)}")
            # Reset initialization flag on failure
            _models_initialized.value = 0
            return False


# Known compounds and biomolecules
PROTEINS = {
    "Follistatin-288": "P19883",  # UniProt ID
    "Follistatin-315": "P19883-2",  # UniProt ID
    "alpha-Klotho": "Q9UEF7",
    "Myoglobin": "P02144",
    "Hemoglobin": ["P69905", "P68871"],  # Alpha and beta chains
    "Profilin": "P07737",
    "Human apolipoprotein E": "P02649",
    "Apolipoprotein A-I Milano": "P02647",
    "Ferritin": ["P02794", "P02792"],  # Heavy and light chains
    "Tubulin alpha": ["P68363", "P68366", "Q71U36", "Q9BQE3", "Q13748"],  # 5 types
    "Tubulin beta": ["P07437", "Q13885", "Q13509", "P68371", "Q9BVA1"],
    "Actin": "P60709",
    "Troponin": ["P19429", "P45379", "P45378"],  # I, T, C
    "Myosin": ["P13533", "P12883", "P13535"],  # Heavy chains
    "BDNF": "P23560",  # Brain-derived neurotrophic factor
    "NGF": "P01138",  # Nerve growth factor
    "GDNF": "P39905",  # Glial cell-derived neurotrophic factor
    "IGF1": "P05019",  # Insulin-like growth factor 1
    "CNTF": "P26441",  # Ciliary neurotrophic factor
}

BIOMOLECULES = {
    "Creatinine": "64-13-7",  # CAS number
    "Creatine": "57-00-1",
    "ATP": "56-65-5",
    "ADP": "58-64-0",
    "AMP": "61-19-8",
    "NAD+": "53-84-9",
    "NADH": "606-68-8",
    "NADP+": "53-59-8",
    "NADPH": "2646-71-1",
    "FAD": "146-14-5",
    "CoA": "85-61-0",
    "Acetyl-CoA": "72-89-9",
    "Glutamate": "56-86-0",
    "GABA": "56-12-2",
    "Dopamine": "51-61-6",
    "Serotonin": "50-67-9",
    "Norepinephrine": "51-41-2",
    "Acetylcholine": "51-84-3",
    "Melatonin": "73-31-4",
    "Histamine": "51-45-6",
}

# Custom compounds dictionary remains unchanged
CUSTOM_COMPOUNDS = {
    # Previous entries remain unchanged
    "EMERGING_THREATS": [
        # Novel synthetic opioids
        {"name": "Protonitazene", "cas": "2276831-46-9", "smiles": "CC(C)N(C)CC1=CC=C(C=C1)N2C(=O)C3=C(N=CN3C)C2=O"},
        {"name": "Isotonitazene", "cas": "14188-81-9", "smiles": "CC(C)N(C)CC1=CC=C(C=C1)N2C(=O)C3=C(N=CN3C)C2=O"},
        {"name": "Metonitazene", "cas": "1239943-76-0", "smiles": "CC(C)N(C)CC1=CC=C(C=C1)N2C(=O)C3=C(N=CN3C)C2=O"},
        # Novel synthetic cannabinoids
        {"name": "ADB-BUTINACA", "cas": "2365460-08-8", "smiles": "CCCC(=O)N1CCC(CC1)N(C)C(=O)NC(C)(C)C1=NC=C2C(=NN(C2=C1)C1=CC=CC=C1)C(=O)OC"},
        {"name": "MDMB-4en-PINACA", "cas": "2316744-97-7", "smiles": "CC(C)(C)NC(=O)C(CC=C)NC(=O)C1=NN(C2=C1C=CN=C2)C1=CC=CC=C1"},
        # Novel designer benzodiazepines
        {"name": "Bromazolam", "cas": "86386-73-4", "smiles": "CC1=NN=C2N1C(=NC(=O)CN2C)C1=CC=C(Br)C=C1"},
        {"name": "Flubromazolam", "cas": "612526-40-6", "smiles": "CC1=NN=C2N1C(=NC(=O)CN2C)C1=CC=C(F)C=C1Br"},
        # Novel psychostimulants
        {"name": "3-MMC", "cas": "1246816-62-5", "smiles": "CC(NC)C(=O)C1=CC=CC(OC)=C1"},
        {"name": "4-CMC", "cas": "1225622-14-9", "smiles": "CNC(C)C(=O)C1=CC=C(Cl)C=C1"},
        # Novel dissociatives
        {"name": "3-HO-PCP", "cas": "943294-78-0", "smiles": "OC1=CC=CC(=C1)C1(N2CCCCC2)CCCCC1"},
        {"name": "3-MeO-PCE", "cas": "72242-03-6", "smiles": "COC1=CC=CC(=C1)C1(CCN)CCCCC1"},
    ],
    "5-HT2_AGONISTS": [
        {"name": "DOI", "cas": "83619-26-1", "smiles": "CC(NC)C(C1=CC(=C(C=C1)OC)I)C"},
        {"name": "2C-B", "cas": "66142-81-2", "smiles": "CC(NC)CC1=CC(=C(C=C1)BR)OC"},
        {"name": "DOM", "cas": "15588-95-1", "smiles": "CC(C)NC(C)C1=CC(=C(C=C1)OC)CC"},
        {"name": "Bromo-DragonFLY", "cas": "502759-67-3", "smiles": "CC(NC)C1=C2C=C(BR)C=CC2=C3C=CC=CC3=C1"},
        {"name": "25I-NBOMe", "cas": "919797-19-6", "smiles": "COc1cc(I)c(OC)cc1CCNCc1ccccc1OC"},
        {"name": "2C-I", "cas": "64584-32-3", "smiles": "CCNCCc1cc(I)c(OC)cc1OC"},
        {"name": "2C-T-7", "cas": "207740-26-9", "smiles": "CCNCCc1cc(SC(C)CC)c(OC)cc1OC"},
        {"name": "TCB-2", "cas": "1404-55-3", "smiles": "CN(C)CCc1c[nH]c2ccc(Br)cc12"},
        {"name": "AL-LAD", "cas": "1352-79-4", "smiles": "CCN(CC)C(=O)C1CN(C)CCc2c1[nH]c3ccc(CC)cc23"},
    ],
    "NMDA_ANTAGONISTS": [
        {"name": "Ketamine", "cas": "6740-88-1", "smiles": "CN1C(=O)CCC(C1=O)(C2=CC=CC=C2)N(C)C"},
        {"name": "Memantine", "cas": "19982-08-2", "smiles": "CC12CC3CC(C1)(CC(C3)(C2)N)C"},
        {"name": "Dextromethorphan", "cas": "125-71-3", "smiles": "CN1CCC23CCCCC2C1CC4=C3C=C(OC)C=C4"},
        {"name": "MK-801", "cas": "77086-22-7", "smiles": "CN1C2CCCCC2C2=C1C=CC=C2"},
        {"name": "PCP", "cas": "77-10-1", "smiles": "c1ccc(C2(N3CCCCC3)CCCCC2)cc1"},
        {"name": "3-MeO-PCP", "cas": "91164-58-8", "smiles": "COc1cccc(C2(N3CCCCC3)CCCCC2)c1"},
        {"name": "DCK", "cas": "111982-50-4", "smiles": "O=C1CC2(c3ccccc3Cl)CC1N2"},
        {"name": "2-FDCK", "cas": "111982-50-4", "smiles": "O=C1CC2(c3ccccc3F)CC1N2"},
        {"name": "Lanicemine", "cas": "102627-32-5", "smiles": "O=C1CC2(c3ccccc3)CC1N2"},
    ],
    "NOOTROPICS": [
        {"name": "Piracetam", "cas": "7491-74-9", "smiles": "O=C1NC(=O)CC1N1CCCC1"},
        {"name": "Aniracetam", "cas": "72432-10-1", "smiles": "O=C1NC(=O)CC1N1CCN(CC1)C(=O)C1=CC=CC=C1"},
        {"name": "Oxiracetam", "cas": "62613-82-5", "smiles": "CC(O)C1=NC(=O)CC(=O)N1"},
        {"name": "Pramiracetam", "cas": "68497-62-1", "smiles": "CC(C)CN1C(=O)CC(N2CCCC2)C(=O)N1"},
        {"name": "Phenylpiracetam", "cas": "77472-70-9", "smiles": "O=C1NC(=O)CC1N1CCCC1C1=CC=CC=C1"},
        {"name": "Noopept", "cas": "157115-85-0", "smiles": "CC(C)CC(=O)N1CCCC1C(=O)NCC(=O)N1CCN(Cc2ccccc2)CC1"},
        {"name": "Semax", "cas": "80714-61-0", "smiles": "CC(C)CC(NC(=O)CNC(=O)C(N)Cc1ccc(O)cc1)C(=O)N1CCCC1C(=O)O"},
        {"name": "Selank", "cas": "129954-34-3", "smiles": "CC(C)CC(NC(=O)CNC(=O)C(N)CCCNC(=O)C(N)CC(=O)N)C(=O)N1CCCC1C(=O)O"},
        {"name": "Coluracetam", "cas": "135463-81-9", "smiles": "CC1(C)CC(=O)N(CCCCN2C(=O)CC3(CC3)C2=O)C1=O"},
        {"name": "Sunifiram", "cas": "314728-85-3", "smiles": "O=C1NC(=O)CC1N1CCN(C(=O)c2ccncc2)CC1"},
    ],
    "ANTI_ADDICTIVE": [
        {"name": "Ibogaine", "cas": "83-74-9", "smiles": "CC1C2CC3C(CN2CCc2c1[nH]c4ccccc24)CC1NCCc2c3[nH]c3ccccc23"},
        {"name": "18-MC", "cas": "2746-81-8", "smiles": "CC1C2CC3C(CN2CCc2c1nc4ccccc24)CC1NCCc2c3nc3ccccc23"},
        {"name": "Naltrexone", "cas": "16590-41-3", "smiles": "Oc1ccc2CC3C4CCC(O)C5Oc1c2C45CCN3"},
        {"name": "Acamprosate", "cas": "77337-76-9", "smiles": "CC(=O)NCCCS(=O)(=O)O"},
        {"name": "Baclofen", "cas": "1134-47-0", "smiles": "O=C(O)CCCC1CCCCC1"},
        {"name": "Varenicline", "cas": "249296-44-4", "smiles": "c1cc2c(cc1O)C(=O)C1CNCC2C1"},
        {"name": "Bupropion", "cas": "34841-39-9", "smiles": "CC(NC(C)(C)C)C(=O)c1cccc(Cl)c1"},
        {"name": "Topiramate", "cas": "97240-79-4", "smiles": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1"},
        {"name": "Gabapentin", "cas": "60142-96-3", "smiles": "NCC1(CC(=O)O)CCCCC1"},
    ],
    "LONGEVITY": [
        {"name": "Rapamycin", "cas": "53123-88-9", "smiles": "CC1CCC2CC(=O)C(=C(O)C2(O1)C(=O)C1=C(O)CC(OC)C(=O)C=C1OC)C(=O)OC"},
        {"name": "Metformin", "cas": "657-24-9", "smiles": "CN(C)C(=N)NC(=N)N"},
        {"name": "NMN", "cas": "1094-61-7", "smiles": "NC(=O)C1=CN(C=CC1=O)[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O"},
        {"name": "NR", "cas": "23111-00-4", "smiles": "NC(=O)C1=CN([C@H]2O[C@H](CO)[C@@H](O)[C@H]2O)C=CC1=O"},
        {"name": "Spermidine", "cas": "124-20-9", "smiles": "NCCCNCCCCNCCCN"},
        {"name": "Resveratrol", "cas": "501-36-0", "smiles": "OC1=CC(O)=CC(C=CC2=CC(O)=CC(O)=C2)=C1"},
        {"name": "Fisetin", "cas": "528-48-3", "smiles": "O=C1C=C(O)C2=C(O)C=C(O)C=C2O1"},
        {"name": "Quercetin", "cas": "117-39-5", "smiles": "O=C1C(O)=C(C(=O)C2=C1C=C(O)C(O)=C2)C1=CC=C(O)C(O)=C1"},
    ],
    "PHYSICAL_ENHANCEMENT": [
        {"name": "Cardarine", "cas": "317318-70-0", "smiles": "CC1=C(C2=CC=C(OCC3=CC=C(NC(=O)C4=CC=C(Cl)C=C4)C=C3)C=C2)C(=O)N(CC(=O)O)C1=O"},
        {"name": "Ostarine", "cas": "841205-47-8", "smiles": "CC1=C(C2=CC=C(OCC3=CC=C(NC(=O)C4=CC=C(F)C=C4)C=C3)C=C2)C(=O)N(CC(=O)O)C1=O"},
        {"name": "SR9009", "cas": "1379686-30-2", "smiles": "CC1=C(C2=CC=C(OCC3=CC=C(NC(=O)C4=CC=C(F)C=C4)C=C3)C=C2)C(=O)N(CC(=O)O)C1=O"},
        {
            "name": "BPC-157",
            "cas": "137525-51-0",
            "smiles": "CC(C)CC(NC(=O)CNC(=O)C(N)CC(=O)N)C(=O)NC(CC(=O)O)C(=O)NC(CC(C)C)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC1=CN=CN1)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CC(=O)O)C(=O)NC(CC1=CC=C(O)C=O",
        },
    ],
}


class CompoundExporter:
    """Handles export of compound data with enrichment from multiple sources."""

    def __init__(self, cache_dir: Optional[Path] = None, n_workers: int = 4):
        self.cache_dir = cache_dir or Path("cache")
        self.n_workers = n_workers
        self.setup_clients()
        self.setup_standardizer()

    def setup_standardizer(self):
        """Initialize structure standardization parameters."""
        self.standardizer = rdMolStandardize.Standardizer()
        self.uncharger = rdMolStandardize.Uncharger()
        self.normalizer = rdMolStandardize.Normalizer()
        self.reionizer = rdMolStandardize.Reionizer()

    def setup_clients(self):
        """Initialize all data source clients with enhanced versions."""
        self.bindingdb = BindingDBSourceEnhanced(cache_dir=self.cache_dir, model_dir=MODEL_DIR, data_dir=DATA_DIR, n_workers=self.n_workers, circuit_config=CIRCUIT_CONFIG)
        self.pubchem = PubChemClientEnhanced(cache_dir=self.cache_dir, circuit_config=CIRCUIT_CONFIG)
        self.chembl = ChEMBLClientEnhanced(cache_dir=self.cache_dir, circuit_config=CIRCUIT_CONFIG)
        self.swiss = SwissClientEnhanced(cache_dir=self.cache_dir, circuit_config=CIRCUIT_CONFIG)
        self.community = CommunityClientEnhanced(cache_dir=self.cache_dir, circuit_config=CIRCUIT_CONFIG)
        self.social = SocialClientEnhanced(cache_dir=self.cache_dir, circuit_config=CIRCUIT_CONFIG)
        self.patents = PatentClientEnhanced(cache_dir=self.cache_dir, circuit_config=CIRCUIT_CONFIG)
        self.pubmed = PubMedClientEnhanced(cache_dir=self.cache_dir, circuit_config=CIRCUIT_CONFIG)
        self.sciencedirect = ScienceDirectClientEnhanced(cache_dir=self.cache_dir, circuit_config=CIRCUIT_CONFIG)
        self.reddit = RedditClientEnhanced(cache_dir=self.cache_dir, circuit_config=CIRCUIT_CONFIG)
        self.bluelight = BluelightClientEnhanced(cache_dir=self.cache_dir, circuit_config=CIRCUIT_CONFIG)

    def setup_directories(self) -> Dict[str, Path]:
        """Create necessary directories."""
        base_dir = Path.cwd()
        data_dir = base_dir / "data"
        cache_dir = base_dir / "cache"
        output_dir = base_dir / "output"
        log_dir = base_dir / "logs"

        for d in [data_dir, cache_dir, output_dir, log_dir]:
            d.mkdir(parents=True, exist_ok=True)

        return {"data": data_dir, "cache": cache_dir, "output": output_dir, "logs": log_dir}

    async def use_patent_search(self, compound: CompoundData) -> List[CompoundData]:
        """Use patent-search MCP server and enhanced patent client to find additional compounds."""
        try:
            compounds = []
            seen_inchi_keys = set()

            # 1. Search patents using compound structure and name with enhanced parameters
            results = await self.use_mcp_tool(
                server_name="patent-search",
                tool_name="search_patents",
                arguments={
                    "query": f"{compound.name} OR {compound.cas_number} OR {compound.inchi_key}",
                    "chemical_structure": compound.smiles,
                    "from_date": "2000-01-01",  # Include historical patents
                },
            )

            if results and "patents" in results:
                for patent in results["patents"]:
                    # 2. Extract compounds from each patent with enhanced extraction
                    extracted = await self.use_mcp_tool(
                        server_name="patent-search",
                        tool_name="extract_compounds",
                        arguments={"patent_number": patent["number"], "include_similar": True, "similarity_threshold": 0.7},  # Get structurally similar compounds
                    )

                    if extracted and "compounds" in extracted:
                        for compound_data in extracted["compounds"]:
                            try:
                                # 3. Use enhanced patent client to get additional data
                                enriched_data = await self.patents.get_compound_details(smiles=compound_data["smiles"], include_similar=True, include_references=True)

                                mol = Chem.MolFromSmiles(compound_data["smiles"])
                                if mol:
                                    mol = self.standardize_structure(mol)
                                    if mol:
                                        # Generate InChI key for deduplication
                                        inchi_key = Chem.MolToInchiKey(mol)

                                        # Skip if we've already seen this compound
                                        if inchi_key in seen_inchi_keys:
                                            continue

                                        seen_inchi_keys.add(inchi_key)

                                        # Create compound with enhanced data
                                        new_compound = CompoundData(
                                            name=compound_data.get("name", enriched_data.get("name", "Unknown")),
                                            cas_number=compound_data.get("cas", enriched_data.get("cas_number", "")),
                                            smiles=Chem.MolToSmiles(mol, isomericSmiles=True),
                                            inchi=Chem.MolToInchi(mol),
                                            inchi_key=inchi_key,
                                            compound_type=compound.compound_type,
                                            molecular_weight=Descriptors.ExactMolWt(mol),
                                            logp=Crippen.MolLogP(mol),
                                            tpsa=Descriptors.TPSA(mol),
                                            hbd=Descriptors.NumHDonors(mol),
                                            hba=Descriptors.NumHAcceptors(mol),
                                            rotatable_bonds=Descriptors.NumRotatableBonds(mol),
                                            patent_refs=[patent["number"]] + (enriched_data.get("additional_patents", []) or []),
                                            source="patent",
                                            activity_predictions=enriched_data.get("activity_predictions"),
                                            toxicity_predictions=enriched_data.get("toxicity_predictions"),
                                            patent_confidence=enriched_data.get("confidence"),
                                            patent_properties=enriched_data.get("properties", {}),
                                            literature_refs=enriched_data.get("literature_refs", []),
                                            target_data=enriched_data.get("target_data", []),
                                        )

                                        # 4. Cross-reference with other data sources
                                        try:
                                            # Check PubChem for additional data
                                            pubchem_data = await self.pubchem.get_compound_data(smiles=new_compound.smiles, cas_number=new_compound.cas_number)
                                            if pubchem_data:
                                                new_compound.merge(pubchem_data)

                                            # Check ChEMBL for bioactivity data
                                            chembl_data = await self.chembl.search_similar_compounds(
                                                new_compound.smiles, similarity=0.9  # High similarity threshold for direct matches
                                            )
                                            if chembl_data:
                                                for match in chembl_data:
                                                    if match.inchi_key == new_compound.inchi_key:
                                                        new_compound.merge(match)
                                                        break

                                        except Exception as e:
                                            logger.warning(f"Error enriching patent compound with additional data: {str(e)}")

                                        compounds.append(new_compound)

                            except Exception as e:
                                logger.error(f"Error processing patent compound: {str(e)}")
                                continue

            # 5. Find similar compounds in related patents
            try:
                similar_results = await self.patents.find_similar_compounds(compound.smiles, similarity_threshold=0.7, include_references=True)

                for similar in similar_results:
                    if similar.inchi_key not in seen_inchi_keys:
                        seen_inchi_keys.add(similar.inchi_key)
                        compounds.append(similar)

            except Exception as e:
                logger.warning(f"Error finding similar compounds: {str(e)}")

            return compounds

        except Exception as e:
            logger.error(f"Patent search failed: {str(e)}")
            return []

    async def find_similar_compounds(self, compound: CompoundData, similarity_threshold: float = 0.7) -> List[CompoundData]:
        """Find structurally similar compounds using multiple sources."""
        similar_compounds = []

        # Search ChEMBL
        try:
            chembl_similar = await self.chembl.search_similar_compounds(compound.smiles, similarity=similarity_threshold)
            similar_compounds.extend(chembl_similar)
        except Exception as e:
            logger.error(f"ChEMBL similarity search failed: {str(e)}")

        # Search PubChem
        try:
            pubchem_similar = await self.pubchem.search_similar_compounds(compound.smiles, similarity=similarity_threshold)
            similar_compounds.extend(pubchem_similar)
        except Exception as e:
            logger.error(f"PubChem similarity search failed: {str(e)}")

        # Remove duplicates based on InChI Key
        unique_compounds = {}
        for compound in similar_compounds:
            if compound.inchi_key and compound.inchi_key not in unique_compounds:
                unique_compounds[compound.inchi_key] = compound

        return list(unique_compounds.values())

    async def collect_compounds(self, custom_compounds_file: Optional[Path] = None) -> List[CompoundData]:
        """Collect compounds from all sources including optional custom file."""
        tasks = [self.collect_5ht2_agonists(), self.collect_nmda_antagonists(), self.collect_proteins(), self.collect_biomolecules(), self.collect_custom_compounds()]

        if custom_compounds_file:
            tasks.append(self.load_custom_compounds_file(custom_compounds_file))

        async with tqdm.create() as progress:
            results = await asyncio.gather(*tasks)

        all_compounds = []
        for compounds in results:
            all_compounds.extend(compounds)

        # Use enhanced BindingDB to find additional compounds
        logger.info("Searching for additional compounds using enhanced BindingDB...")
        try:
            additional = await self.bindingdb.process_all_sources()
            if additional:
                all_compounds.extend(additional)
                logger.info(f"Found {len(additional)} additional compounds")
        except Exception as e:
            logger.error(f"Error getting additional compounds: {str(e)}")

        return self.merge_duplicates(all_compounds)

    async def collect_5ht2_agonists(self) -> List[CompoundData]:
        """Collect 5-HT2 agonist compounds using enhanced sources."""
        logger.info("Collecting 5-HT2 agonists...")

        compounds = []

        # 1. Get initial compounds from enhanced BindingDB
        initial_compounds = await self.bindingdb.search_compounds(target="5-HT2", activity_type="agonist", min_activity=7.0, include_patents=True, include_similar=True)

        # 2. Search patents for each compound to find related compounds
        for compound in initial_compounds:
            try:
                # Get enhanced data
                enriched = await self.bindingdb.enrich_compound_data(compound, include_patents=True, include_community=True, include_social=True)

                if enriched:
                    compounds.append(enriched)

                    # Use patent search to find additional compounds
                    patent_compounds = await self.use_patent_search(enriched)
                    if patent_compounds:
                        compounds.extend(patent_compounds)

                        # For each new compound found, also search for similar compounds
                        for patent_compound in patent_compounds:
                            similar = await self.find_similar_compounds(patent_compound, similarity_threshold=0.7)
                            compounds.extend(similar)

            except Exception as e:
                logger.error(f"Error processing compound {compound.name}: {str(e)}")

        # Remove duplicates
        return self.merge_duplicates(compounds)

    async def collect_nmda_antagonists(self) -> List[CompoundData]:
        """Collect NMDA antagonist compounds using enhanced sources."""
        logger.info("Collecting NMDA antagonists...")

        compounds = []

        # 1. Get initial compounds from enhanced BindingDB
        initial_compounds = await self.bindingdb.search_compounds(target="NMDA", activity_type="antagonist", min_activity=6.5, include_patents=True, include_similar=True)

        # 2. Search patents and find similar compounds for each initial compound
        for compound in initial_compounds:
            try:
                # Get enhanced data
                enriched = await self.bindingdb.enrich_compound_data(compound, include_patents=True, include_community=True, include_social=True)

                if enriched:
                    compounds.append(enriched)

                    # Use patent search to find additional compounds
                    patent_compounds = await self.use_patent_search(enriched)
                    if patent_compounds:
                        compounds.extend(patent_compounds)

                        # For each new compound found, also search for similar compounds
                        for patent_compound in patent_compounds:
                            similar = await self.find_similar_compounds(patent_compound, similarity_threshold=0.7)
                            compounds.extend(similar)

                            # Additional search for NMDA-specific analogs
                            try:
                                nmda_analogs = await self.chembl.search_compounds(
                                    target="NMDA",
                                    activity_type="antagonist",
                                    similar_to=patent_compound.smiles,
                                    similarity_threshold=0.6,  # Lower threshold to catch more potential analogs
                                )
                                if nmda_analogs:
                                    compounds.extend(nmda_analogs)
                            except Exception as e:
                                logger.warning(f"Error searching NMDA analogs: {str(e)}")

            except Exception as e:
                logger.error(f"Error processing compound {compound.name}: {str(e)}")

        # Remove duplicates and return
        return self.merge_duplicates(compounds)

    async def collect_proteins(self) -> List[CompoundData]:
        """Collect protein/peptide data using enhanced Swiss client."""
        logger.info("Collecting proteins...")
        compounds = []

        for name, uniprot_ids in PROTEINS.items():
            if isinstance(uniprot_ids, str):
                uniprot_ids = [uniprot_ids]

            for uniprot_id in uniprot_ids:
                try:
                    # Get initial protein data with enhanced features
                    data = await self.swiss.get_protein_data(
                        uniprot_id, include_patents=True, include_literature=True, include_interactions=True  # Get protein-protein interactions
                    )

                    if data:
                        # Search for related compounds in patents
                        try:
                            patent_compounds = await self.use_patent_search(data)
                            if patent_compounds:
                                compounds.extend(patent_compounds)

                                # For each patent compound, look for similar compounds
                                for patent_compound in patent_compounds:
                                    similar = await self.find_similar_compounds(patent_compound, similarity_threshold=0.7)
                                    compounds.extend(similar)

                        except Exception as e:
                            logger.warning(f"Error searching patents for protein {name}: {str(e)}")

                        # Add the original protein data
                        compounds.append(data)

                        # Get protein family members
                        try:
                            family_members = await self.swiss.get_protein_family(uniprot_id)
                            if family_members:
                                for member in family_members:
                                    member_data = await self.swiss.get_protein_data(member, include_patents=True, include_literature=True, include_interactions=True)
                                    if member_data:
                                        compounds.append(member_data)
                        except Exception as e:
                            logger.warning(f"Error getting protein family data for {name}: {str(e)}")

                except Exception as e:
                    logger.error(f"Error getting protein data for {uniprot_id}: {str(e)}")

        # Remove duplicates and return
        return self.merge_duplicates(compounds)

    async def collect_biomolecules(self) -> List[CompoundData]:
        """Collect basic biomolecule data with enhanced enrichment."""
        logger.info("Collecting biomolecules...")
        compounds = []

        for name, cas in BIOMOLECULES.items():
            try:
                # Get initial data from PubChem
                data = await self.pubchem.get_compound_data(name=name, cas_number=cas)
                if data:
                    # Enrich with enhanced sources
                    enriched = await self.bindingdb.enrich_compound_data(data, include_patents=True, include_community=True, include_social=True)

                    if enriched:
                        compounds.append(enriched)

                        # Search patents for related compounds
                        try:
                            patent_compounds = await self.use_patent_search(enriched)
                            if patent_compounds:
                                compounds.extend(patent_compounds)

                                # For each patent compound, look for similar compounds
                                for patent_compound in patent_compounds:
                                    similar = await self.find_similar_compounds(patent_compound, similarity_threshold=0.7)
                                    compounds.extend(similar)

                                    # Look for biochemical analogs
                                    try:
                                        biochem_analogs = await self.chembl.search_compounds(
                                            target=name,  # Use biomolecule name as target
                                            similar_to=patent_compound.smiles,
                                            similarity_threshold=0.6,  # Lower threshold for biochemical similarity
                                            include_metabolites=True,  # Include metabolic products
                                            include_natural_variants=True,  # Include natural variants
                                        )
                                        if biochem_analogs:
                                            compounds.extend(biochem_analogs)
                                    except Exception as e:
                                        logger.warning(f"Error searching biochemical analogs for {name}: {str(e)}")

                        except Exception as e:
                            logger.warning(f"Error searching patents for biomolecule {name}: {str(e)}")

                        # Search for metabolic derivatives
                        try:
                            metabolites = await self.chembl.get_metabolites(enriched.smiles)
                            if metabolites:
                                compounds.extend(metabolites)
                        except Exception as e:
                            logger.warning(f"Error getting metabolites for {name}: {str(e)}")

                        # Search for natural variants
                        try:
                            variants = await self.chembl.get_natural_variants(enriched.smiles)
                            if variants:
                                compounds.extend(variants)
                        except Exception as e:
                            logger.warning(f"Error getting natural variants for {name}: {str(e)}")

            except Exception as e:
                logger.error(f"Error getting biomolecule data for {name}: {str(e)}")

        # Remove duplicates and return
        return self.merge_duplicates(compounds)

    async def collect_custom_compounds(self) -> List[CompoundData]:
        """Collect custom compound entries with enhanced enrichment."""
        logger.info("Collecting custom compounds...")
        compounds = []

        with ThreadPoolExecutor() as executor:
            futures = []
            for category, entries in CUSTOM_COMPOUNDS.items():
                compound_type = CompoundType.PSYCHOACTIVE
                if "NOOTROPIC" in category:
                    compound_type = CompoundType.NOOTROPIC

                for entry in entries:
                    futures.append(executor.submit(self.create_custom_compound, entry["name"], entry["cas"], entry["smiles"], compound_type))

            for future in as_completed(futures):
                try:
                    compound = future.result()
                    if compound:
                        # Enrich with enhanced sources
                        enriched = await self.bindingdb.enrich_compound_data(compound, include_patents=True, include_community=True, include_social=True)
                        if enriched:
                            compounds.append(enriched)

                            # Search patents for related compounds
                            try:
                                patent_compounds = await self.use_patent_search(enriched)
                                if patent_compounds:
                                    compounds.extend(patent_compounds)

                                    # For each patent compound, look for similar compounds
                                    for patent_compound in patent_compounds:
                                        similar = await self.find_similar_compounds(patent_compound, similarity_threshold=0.7)
                                        compounds.extend(similar)

                                        # Look for category-specific analogs
                                        try:
                                            if "NOOTROPIC" in category:
                                                analogs = await self.chembl.search_compounds(
                                                    target="nootropic",
                                                    similar_to=patent_compound.smiles,
                                                    similarity_threshold=0.6,
                                                    include_cognitive_enhancers=True,
                                                )
                                                if analogs:
                                                    compounds.extend(analogs)
                                            elif "ANTI_ADDICTIVE" in category:
                                                analogs = await self.chembl.search_compounds(
                                                    target="addiction",
                                                    similar_to=patent_compound.smiles,
                                                    similarity_threshold=0.6,
                                                    include_therapeutic_agents=True,
                                                )
                                                if analogs:
                                                    compounds.extend(analogs)
                                            elif "LONGEVITY" in category:
                                                analogs = await self.chembl.search_compounds(
                                                    target="aging",
                                                    similar_to=patent_compound.smiles,
                                                    similarity_threshold=0.6,
                                                    include_metabolic_modulators=True,
                                                )
                                                if analogs:
                                                    compounds.extend(analogs)
                                            elif "PHYSICAL_ENHANCEMENT" in category:
                                                analogs = await self.chembl.search_compounds(
                                                    target="performance",
                                                    similar_to=patent_compound.smiles,
                                                    similarity_threshold=0.6,
                                                    include_ergogenic_aids=True,
                                                )
                                                if analogs:
                                                    compounds.extend(analogs)
                                            elif "EMERGING_THREATS" in category:
                                                analogs = await self.chembl.search_compounds(
                                                    target="abuse",
                                                    similar_to=patent_compound.smiles,
                                                    similarity_threshold=0.6,
                                                    include_novel_substances=True,
                                                )
                                                if analogs:
                                                    compounds.extend(analogs)
                                        except Exception as e:
                                            logger.warning(f"Error searching category-specific analogs for {category}: {str(e)}")

                            except Exception as e:
                                logger.warning(f"Error searching patents for custom compound {compound.name}: {str(e)}")

                except Exception as e:
                    logger.error(f"Error creating/enriching custom compound: {str(e)}")

        # Remove duplicates and return
        return self.merge_duplicates(compounds)

    async def load_custom_compounds_file(self, file_path: Path) -> List[CompoundData]:
        """Load compounds from custom TSV/JSON file with enhanced enrichment."""
        if not file_path.exists():
            logger.error(f"Custom compounds file not found: {file_path}")
            return []

        compounds = []
        entries = []

        # Load entries from file
        if file_path.suffix.lower() == ".json":
            with open(file_path) as f:
                entries = json.load(f)
        elif file_path.suffix.lower() in [".tsv", ".csv"]:
            df = pd.read_csv(file_path, sep="\t" if file_path.suffix.lower() == ".tsv" else ",")
            entries = df.to_dict("records")

        # Process each entry
        for entry in entries:
            try:
                compound = await self.create_compound_from_dict(entry)
                if compound:
                    # Enrich with enhanced sources
                    enriched = await self.bindingdb.enrich_compound_data(compound, include_patents=True, include_community=True, include_social=True)
                    if enriched:
                        compounds.append(enriched)

                        # Search patents for related compounds
                        try:
                            patent_compounds = await self.use_patent_search(enriched)
                            if patent_compounds:
                                compounds.extend(patent_compounds)

                                # For each patent compound, look for similar compounds
                                for patent_compound in patent_compounds:
                                    similar = await self.find_similar_compounds(patent_compound, similarity_threshold=0.7)
                                    compounds.extend(similar)

                                    # Look for category-specific analogs based on compound type
                                    try:
                                        if compound.compound_type == CompoundType.NOOTROPIC:
                                            analogs = await self.chembl.search_compounds(
                                                target="nootropic",
                                                similar_to=patent_compound.smiles,
                                                similarity_threshold=0.6,
                                                include_cognitive_enhancers=True,
                                            )
                                            if analogs:
                                                compounds.extend(analogs)
                                        elif compound.compound_type == CompoundType.PSYCHOACTIVE:
                                            analogs = await self.chembl.search_compounds(
                                                target="psychoactive",
                                                similar_to=patent_compound.smiles,
                                                similarity_threshold=0.6,
                                                include_novel_substances=True,
                                            )
                                            if analogs:
                                                compounds.extend(analogs)
                                    except Exception as e:
                                        logger.warning(f"Error searching category-specific analogs for {compound.name}: {str(e)}")

                        except Exception as e:
                            logger.warning(f"Error searching patents for custom compound {compound.name}: {str(e)}")

                        # Search for metabolic derivatives
                        try:
                            metabolites = await self.chembl.get_metabolites(enriched.smiles)
                            if metabolites:
                                compounds.extend(metabolites)
                        except Exception as e:
                            logger.warning(f"Error getting metabolites for {compound.name}: {str(e)}")

                        # Search for natural variants
                        try:
                            variants = await self.chembl.get_natural_variants(enriched.smiles)
                            if variants:
                                compounds.extend(variants)
                        except Exception as e:
                            logger.warning(f"Error getting natural variants for {compound.name}: {str(e)}")

            except Exception as e:
                logger.error(f"Error processing custom compound: {str(e)}")

        # Remove duplicates and return
        return self.merge_duplicates(compounds)

    def create_custom_compound(self, name: str, cas: str, smiles: str, compound_type: CompoundType = CompoundType.OTHER, retry_count: int = 3) -> Optional[CompoundData]:
        """Create a custom compound with retry logic."""
        for attempt in range(retry_count):
            try:
                mol = Chem.MolFromSmiles(smiles)
                if not mol:
                    raise ValueError(f"Invalid SMILES for {name}")

                mol = self.standardize_structure(mol)
                if not mol:
                    raise ValueError(f"Failed to standardize {name}")

                return CompoundData(
                    name=name,
                    cas_number=cas,
                    smiles=Chem.MolToSmiles(mol, isomericSmiles=True),
                    inchi=Chem.MolToInchi(mol),
                    inchi_key=Chem.MolToInchiKey(mol),
                    compound_type=compound_type,
                    molecular_weight=Descriptors.ExactMolWt(mol),
                    logp=Crippen.MolLogP(mol),
                    tpsa=Descriptors.TPSA(mol),
                    hbd=Descriptors.NumHDonors(mol),
                    hba=Descriptors.NumHAcceptors(mol),
                    rotatable_bonds=Descriptors.NumRotatableBonds(mol),
                )

            except Exception as e:
                if attempt == retry_count - 1:
                    logger.error(f"Failed to create compound {name}: {str(e)}")
                    return None
                time.sleep(1)

        return None

    def standardize_structure(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Enhanced structure standardization with retry logic."""
        try:
            if not mol:
                return None

            # Remove salt
            mol = self.uncharger.uncharge(mol)

            # Normalize structure
            mol = self.normalizer.normalize(mol)

            # Reionize
            mol = self.reionizer.reionize(mol)

            # Standardize
            mol = self.standardizer.standardize(mol)

            # Generate 3D conformation if needed
            if not mol.GetConformer().Is3D():
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol)

            return mol

        except Exception as e:
            logger.error(f"Error standardizing structure: {str(e)}")
            return None

    async def enrich_compounds(self, compounds: List[CompoundData], retry_count: int = 3) -> List[CompoundData]:
        """Enrich compound data with focus on abuse potential monitoring."""
        enriched = []
        for compound in compounds:
            try:
                # Patent data
                patent_data = await self.patents.search_compound(smiles=compound.smiles, name=compound.name)
                if patent_data:
                    compound.patent_refs = patent_data.references
                    compound.patent_properties = patent_data.properties

                # Literature data
                lit_data = await self.pubmed.search_compound(name=compound.name, cas_number=compound.cas_number)
                if lit_data:
                    compound.literature_refs = lit_data.references

                # Community data
                community_data = await self.community.get_compound_mentions(name=compound.name, smiles=compound.smiles)
                if community_data:
                    compound.community_refs = community_data.references
                    compound.safety_reports = community_data.safety_reports

                # Enhanced social media monitoring for abuse patterns
                social_data = await self.social.get_compound_mentions(name=compound.name, monitor_abuse_patterns=True, track_regional_trends=True, alert_on_clusters=True)
                if social_data:
                    compound.social_refs = social_data.references
                    compound.abuse_patterns = social_data.abuse_patterns
                    compound.regional_trends = social_data.regional_trends
                    compound.alert_clusters = social_data.alert_clusters

                # Calculate abuse threat score (0-100)
                compound.abuse_threat_score = self._calculate_abuse_threat_score(
                    mentions_velocity=social_data.mentions_velocity,
                    abuse_patterns=social_data.abuse_patterns,
                    regional_spread=len(social_data.regional_trends),
                    alert_clusters=len(social_data.alert_clusters),
                )

                # Enhanced validation with abuse focus
                if validate_compound_data(compound):
                    enriched.append(compound)

            except Exception as e:
                logger.error(f"Error enriching {compound.name}: {str(e)}")

        return enriched

    def _calculate_abuse_threat_score(self, mentions_velocity: float, abuse_patterns: List[str], regional_spread: int, alert_clusters: int) -> float:
        """Calculate abuse threat score from monitoring data."""
        # Velocity score (0-25)
        velocity_score = min(25, mentions_velocity * 5)

        # Patterns score (0-25)
        patterns_score = min(25, len(abuse_patterns) * 5)

        # Spread score (0-25)
        spread_score = min(25, regional_spread * 2)

        # Alert score (0-25)
        alert_score = min(25, alert_clusters * 5)

        return velocity_score + patterns_score + spread_score + alert_score

    def merge_duplicates(self, compounds: List[CompoundData]) -> List[CompoundData]:
        """Merge duplicate compounds based on InChI Key."""
        logger.info("Merging duplicate compounds...")
        merged = {}
        for compound in compounds:
            if not compound.inchi_key:
                continue
            if compound.inchi_key in merged:
                merged[compound.inchi_key].merge(compound)
            else:
                merged[compound.inchi_key] = compound
        return list(merged.values())

    def export_data(self, compounds: List[CompoundData], output_dir: Path, formats: Optional[List[str]] = None) -> None:
        """Export compounds in multiple formats with progress tracking."""
        formats = formats or ["tsv", "json", "excel", "sdf", "mol"]
        output_dir.mkdir(parents=True, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        with tqdm(total=len(formats), desc="Exporting formats") as pbar:
            for fmt in formats:
                try:
                    if fmt == "tsv":
                        self.export_tsv(compounds, output_dir / f"compounds_{timestamp}.tsv")
                    elif fmt == "json":
                        self.export_json(compounds, output_dir / f"compounds_{timestamp}.json")
                    elif fmt == "excel":
                        self.export_excel(compounds, output_dir / f"compounds_{timestamp}.xlsx")
                    elif fmt == "sdf":
                        self.export_sdf(compounds, output_dir / f"compounds_{timestamp}.sdf")
                    elif fmt == "mol":
                        self.export_mol_files(compounds, output_dir / "mol")
                    pbar.update(1)
                except Exception as e:
                    logger.error(f"Error exporting {fmt} format: {str(e)}")

    def export_tsv(self, compounds: List[CompoundData], output_path: Path) -> None:
        """Export compounds to TSV format."""
        fieldnames = [
            "name",
            "cas_number",
            "smiles",
            "inchi",
            "inchi_key",
            "compound_type",
            "molecular_weight",
            "logp",
            "tpsa",
            "hbd",
            "hba",
            "rotatable_bonds",
            "target",
            "activity_type",
            "activity_value",
            "patent_count",
            "community_mentions",
            "safety_score",
            "data_quality",
            "doi_refs",
            "pmid_refs",
            "patent_refs",
            "source",
            "prediction_confidence",
            "activity_predictions",
            "toxicity_predictions",
            "abuse_predictions",
            "community_safety_score",
            "literature_confidence",
            "patent_confidence",
            "social_confidence",
            "total_references",
            "validation_status",
            "last_updated",
        ]

        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for compound in compounds:
                row = compound.to_dict()

                # Add reference counts
                row["doi_refs"] = ";".join(compound.literature_refs or [])
                row["pmid_refs"] = ";".join(compound.pubmed_refs or [])
                row["patent_refs"] = ";".join(compound.patent_refs or [])
                row["total_references"] = len((compound.literature_refs or []) + (compound.pubmed_refs or []) + (compound.patent_refs or []))

                # Add prediction data
                if hasattr(compound, "activity_predictions"):
                    row["activity_predictions"] = json.dumps(compound.activity_predictions)
                if hasattr(compound, "toxicity_predictions"):
                    row["toxicity_predictions"] = json.dumps(compound.toxicity_predictions)
                if hasattr(compound, "abuse_predictions"):
                    row["abuse_predictions"] = json.dumps(compound.abuse_predictions)

                # Add confidence scores
                row["prediction_confidence"] = getattr(compound, "prediction_confidence", None)
                row["literature_confidence"] = getattr(compound, "literature_confidence", None)
                row["patent_confidence"] = getattr(compound, "patent_confidence", None)
                row["social_confidence"] = getattr(compound, "social_confidence", None)

                # Add metadata
                row["validation_status"] = getattr(compound, "validation_status", "unvalidated")
                row["last_updated"] = datetime.now().isoformat()

                writer.writerow(row)

    def export_json(self, compounds: List[CompoundData], output_path: Path) -> None:
        """Export compounds to JSON format with full data."""
        with open(output_path, "w") as f:
            json.dump([c.to_dict(include_all=True) for c in compounds], f, indent=2)

    def export_excel(self, compounds: List[CompoundData], output_path: Path) -> None:
        """Export compounds to Excel format for easier viewing."""
        df = pd.DataFrame([c.to_dict() for c in compounds])
        df.to_excel(output_path, index=False, engine="openpyxl")

    def export_sdf(self, compounds: List[CompoundData], output_path: Path) -> None:
        """Export compounds to SDF format."""
        writer = SDWriter(str(output_path))
        for compound in compounds:
            if compound.smiles:
                mol = Chem.MolFromSmiles(compound.smiles)
                if mol:
                    mol = self.standardize_structure(mol)
                    if mol:
                        for key, value in compound.to_dict().items():
                            if value is not None:
                                mol.SetProp(key, str(value))
                        writer.write(mol)
        writer.close()

    def export_mol_files(self, compounds: List[CompoundData], output_dir: Path) -> None:
        """Export individual MOL files for each compound."""
        output_dir.mkdir(parents=True, exist_ok=True)
        for compound in compounds:
            if compound.smiles:
                mol = Chem.MolFromSmiles(compound.smiles)
                if mol:
                    mol = self.standardize_structure(mol)
                    if mol:
                        # Create a safe filename from compound name
                        safe_name = "".join(c for c in compound.name if c.isalnum() or c in (" ", "-", "_")).rstrip()
                        mol_path = output_dir / f"{safe_name}.mol"
                        try:
                            Chem.MolToMolFile(mol, str(mol_path))
                        except Exception as e:
                            logger.error(f"Error saving MOL file for {compound.name}: {str(e)}")


async def main():
    """Main async function to run the export process."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--custom-compounds", type=Path, help="Path to custom compounds file (TSV/JSON)")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Output directory")
    parser.add_argument("--formats", nargs="+", default=["tsv", "json", "excel", "sdf", "mol"], help="Output formats")
    args = parser.parse_args()

    # Initialize exporter
    exporter = CompoundExporter(cache_dir=CACHE_DIR)

    # Collect compounds
    compounds = await exporter.collect_compounds(args.custom_compounds)
    logger.info(f"Collected {len(compounds)} compounds")

    # Enrich compounds with additional data
    enriched = await exporter.enrich_compounds(compounds)
    logger.info(f"Enriched {len(enriched)} compounds")

    # Export in requested formats
    exporter.export_data(enriched, args.output_dir, args.formats)
    logger.info(f"Export complete to {args.output_dir}")


def run_main():
    """Run the main function with proper error handling."""
    try:
        # Initialize models in the main process
        if not initialize_models():
            logger.error("Failed to initialize models")
            sys.exit(1)

        # Run the main async function
        try:
            asyncio.run(main())
            sys.exit(0)
        except asyncio.CancelledError:
            logger.warning("Operation cancelled")
            sys.exit(130)
        except Exception as e:
            logger.error(f"Error in main execution: {str(e)}")
            sys.exit(1)
    except KeyboardInterrupt:
        logger.warning("Export interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Export failed: {str(e)}")
        sys.exit(1)


def main_wrapper():
    """Wrapper function to handle multiprocessing initialization."""
    try:
        # Initialize multiprocessing support
        multiprocessing.freeze_support()

        # Set multiprocessing start method to 'spawn' for all platforms
        if sys.platform != "win32":  # Windows already uses 'spawn'
            try:
                multiprocessing.set_start_method("spawn")
            except RuntimeError:
                # Method may already be set
                pass

        # Initialize models in the main process
        if not initialize_models():
            logger.error("Failed to initialize models")
            sys.exit(1)

        # Run the main async function
        try:
            asyncio.run(main())
            sys.exit(0)
        except asyncio.CancelledError:
            logger.warning("Operation cancelled")
            sys.exit(130)
        except Exception as e:
            logger.error(f"Error in main execution: {str(e)}")
            sys.exit(1)
    except KeyboardInterrupt:
        logger.warning("Export interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Export failed: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    try:
        # Initialize multiprocessing support
        multiprocessing.freeze_support()

        # Set multiprocessing start method to 'spawn' for all platforms
        if sys.platform != "win32":  # Windows already uses 'spawn'
            try:
                multiprocessing.set_start_method("spawn")
            except RuntimeError:
                # Method may already be set
                pass

        # Set event loop policy for Windows
        if sys.platform == "win32":
            asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

        # Initialize models in the main process before any multiprocessing
        if not initialize_models():
            logger.error("Failed to initialize models")
            sys.exit(1)

        # Run the main async function
        try:
            asyncio.run(main())
            sys.exit(0)
        except asyncio.CancelledError:
            logger.warning("Operation cancelled")
            sys.exit(130)
        except Exception as e:
            logger.error(f"Error in main execution: {str(e)}")
            sys.exit(1)
    except KeyboardInterrupt:
        logger.warning("Export interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Export failed: {str(e)}")
        sys.exit(1)
