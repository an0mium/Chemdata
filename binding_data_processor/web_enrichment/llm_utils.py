"""LLM utilities for text analysis and extraction.

This module provides enhanced text analysis capabilities using LLMs:
1. Extract chemical information and structures
2. Analyze scientific and patent text
3. Extract effects and mechanisms
4. Extract safety information
5. Extract community insights
6. Validate and cross-reference information
7. Process PDF documents
8. Analyze transporter interactions
"""

import logging
import re
import json
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field
from datetime import datetime
from io import BytesIO

from transformers import pipeline
from sentence_transformers import SentenceTransformer
import torch
from crawl4ai import AsyncWebCrawler, BrowserConfig
from crawl4ai.extraction_strategy import LLMExtractionStrategy
import requests
import pypdf
from rdkit import Chem

from binding_data_processor.logger import LogManager

logger = LogManager().get_logger("web_enrichment.llm_utils")

__all__ = [
    "extract_chemical_info",
    "analyze_patent_text",
    "extract_biological_targets",
    "analyze_activity_data",
    "NootropicLLMExtractor",
    "LLMProcessor",  # Ensure this is listed
    "Effect",
    "Mechanism",
    "SafetyData",
    "CommunityInsight",
    "ExtractionResult",
    "analyze_content_with_llm",
    "extract_patent_compound",
    "extract_transporter_data",
    "process_pdf_content",
]


# Initialize models flag
MODELS_LOADED = False


def initialize_models():
    """Initialize required models."""
    global MODELS_LOADED, ner_model, analyzer, relevance_model, nootropic_model

    try:
        if MODELS_LOADED:
            return True

        # Use BERT-based model for chemical entity recognition
        ner_model = pipeline(
            "token-classification",
            model="allenai/scibert_scivocab_uncased",
            aggregation_strategy="simple",
            from_pt=True,  # Use PyTorch weights
        )

        # Use SciBERT for text analysis
        analyzer = SentenceTransformer(
            "allenai/scibert_scivocab_uncased",
            from_pt=True,  # Use PyTorch weights
        )

        # Use domain-specific model for relevance scoring
        relevance_model = SentenceTransformer(
            "pritamdeka/S-PubMedBert-MS-MARCO",
            from_pt=True,  # Use PyTorch weights
        )

        # Use PubMedBERT for nootropic analysis
        nootropic_model = SentenceTransformer(
            "microsoft/BiomedNLP-PubMedBert-base-uncased-abstract",
            from_pt=True,  # Use PyTorch weights
        )

        MODELS_LOADED = True
        return True
    except Exception as e:
        logging.error(f"Error loading models: {str(e)}")
        MODELS_LOADED = False
        return False


class LLMProcessor:
    """Base class for LLM-based text processing."""

    def __init__(
        self,
        model: str = "gpt-4",
        api_key: Optional[str] = None,
        temperature: float = 0.0,
        max_tokens: int = 1000,
        cache_dir: Optional[str] = None,
        device: Optional[str] = None,
    ):
        """Initialize LLM processor.

        Args:
            model: LLM model to use
            api_key: Optional API key
            temperature: Sampling temperature
            max_tokens: Maximum tokens to generate
            cache_dir: Optional cache directory
            device: Optional device to use (cpu, cuda, mps)
        """
        # Initialize models if not already loaded
        if not MODELS_LOADED:
            if not initialize_models():
                raise RuntimeError("Failed to initialize required models")

        self.model = model
        self.api_key = api_key
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.cache_dir = cache_dir

        # Set device
        if device:
            self.device = device
        else:
            if torch.backends.mps.is_available():
                self.device = "mps"
            elif torch.cuda.is_available():
                self.device = "cuda"
            else:
                self.device = "cpu"

        # Configure extraction strategy
        self.strategy = LLMExtractionStrategy(
            input_format="text",
            provider=model,
            api_token=api_key,
            temperature=temperature,
            max_tokens=max_tokens,
            chunking={
                "type": "sliding_window",
                "window_size": 1000,
                "step_size": 500,
                "overlap": 0.5,
                "min_chunk_size": 100,
            },
        )

    async def process_text(self, text: str, instruction: str) -> Dict[str, Any]:
        """Process text with LLM using given instruction.

        Args:
            text: Text to process
            instruction: Instruction for LLM processing

        Returns:
            Dictionary containing processed results
        """
        try:
            # Update strategy instruction
            self.strategy.instruction = instruction

            # Process with LLM
            async with AsyncWebCrawler() as crawler:
                result = await crawler.arun(
                    text=text,
                    extraction_strategy=self.strategy,
                )

                if result.success and result.extracted_content:
                    return {
                        "content": result.extracted_content,
                        "confidence": self._calculate_confidence(result.extracted_content),
                        "metadata": {
                            "model": self.model,
                            "timestamp": datetime.now().isoformat(),
                        },
                    }

        except Exception as e:
            logger.error(f"Error processing text with LLM: {str(e)}")

        return {}

    def _calculate_confidence(self, content: Any) -> float:
        """Calculate confidence score for extracted content."""
        if isinstance(content, (list, tuple)):
            if not content:
                return 0.0
            return sum(self._calculate_confidence(item) for item in content) / len(content)

        if isinstance(content, dict):
            if "confidence" in content:
                return float(content["confidence"])
            if not content:
                return 0.0

            return sum(self._calculate_confidence(v) for v in content.values()) / len(content)

        return 1.0 if content else 0.0


# Helper functions for preprocessing
def _preprocess_chemical_content(content: str) -> str:
    """Preprocess content to highlight chemical information."""
    try:
        # Remove HTML tags
        content = re.sub(r"<[^>]+>", " ", content)

        # Normalize whitespace
        content = re.sub(r"\s+", " ", content)

        # Highlight chemical formulas
        content = re.sub(r"([A-Z][a-z]?\d*)+", r" \1 ", content)

        # Highlight CAS numbers
        content = re.sub(r"\b\d{1,7}-\d{2}-\d\b", r" \0 ", content)

        return content.strip()
    except Exception as e:
        logger.error(f"Error preprocessing chemical content: {str(e)}")
        return content


def _preprocess_patent_text(text: str) -> str:
    """Preprocess patent text for analysis."""
    try:
        # Remove patent-specific formatting
        text = re.sub(r"\[\d+\]", " ", text)
        text = re.sub(r"\((\d+)\)", r" \1 ", text)

        # Clean up whitespace
        text = re.sub(r"\s+", " ", text)

        return text.strip()
    except Exception as e:
        logger.error(f"Error preprocessing patent text: {str(e)}")
        return text


def _preprocess_transporter_content(content: str) -> str:
    """Preprocess content for transporter analysis."""
    try:
        # Highlight transporter names
        for transporter in TRANSPORTER_PATTERNS:
            for name in TRANSPORTER_PATTERNS[transporter]["names"]:
                content = re.sub(rf"\b{re.escape(name)}\b", f" {name} ", content, flags=re.IGNORECASE)

        # Clean up whitespace
        content = re.sub(r"\s+", " ", content)

        return content.strip()
    except Exception as e:
        logger.error(f"Error preprocessing transporter content: {str(e)}")
        return content


def _split_into_sections(text: str) -> List[str]:
    """Split text into logical sections."""
    try:
        # Split on common section headers
        sections = re.split(r"\n\s*(?:EXAMPLE|DESCRIPTION|CLAIMS|ABSTRACT|BACKGROUND|SUMMARY|DETAILED DESCRIPTION|EXPERIMENTAL|REFERENCES)\s*\n", text)

        # Remove empty sections and clean whitespace
        sections = [s.strip() for s in sections if s.strip()]

        return sections
    except Exception as e:
        logger.error(f"Error splitting into sections: {str(e)}")
        return [text]


def _validate_chemical_data(data: Dict[str, Any]) -> Dict[str, Any]:
    """Validate and standardize extracted chemical data."""
    try:
        validated = {}

        # Required fields
        required = ["name", "cas_number", "smiles", "inchi"]
        for field in required:
            validated[field] = data.get(field, "")

        # Optional fields with defaults
        validated["properties"] = data.get("properties", {})
        validated["references"] = data.get("references", [])
        validated["confidence"] = float(data.get("confidence", 0.0))

        # Validate SMILES if present
        if validated["smiles"]:
            mol = Chem.MolFromSmiles(validated["smiles"])
            if mol:
                validated["smiles"] = Chem.MolToSmiles(mol, isomericSmiles=True)

        return validated
    except Exception as e:
        logger.error(f"Error validating chemical data: {str(e)}")
        return data


def _validate_patent_data(data: Dict[str, Any]) -> Dict[str, Any]:
    """Validate and standardize extracted patent data."""
    try:
        validated = {}

        # Required fields
        required = ["chemical_identifiers", "chemical_names", "structure_details"]
        for field in required:
            validated[field] = data.get(field, {})

        # Optional fields
        validated["physical_properties"] = data.get("physical_properties", {})
        validated["preparation_method"] = data.get("preparation_method", {})
        validated["analytical_data"] = data.get("analytical_data", {})
        validated["biological_activity"] = data.get("biological_activity", {})

        return validated
    except Exception as e:
        logger.error(f"Error validating patent data: {str(e)}")
        return data


def _validate_transporter_data(data: Dict[str, Any]) -> Dict[str, Any]:
    """Validate and standardize extracted transporter data."""
    try:
        validated = {}

        # Required fields
        required = ["lat1_interactions", "other_transporters", "transport_mechanisms"]
        for field in required:
            validated[field] = data.get(field, {})

        # Optional fields
        validated["structure_activity"] = data.get("structure_activity", {})
        validated["supporting_evidence"] = data.get("supporting_evidence", [])
        validated["confidence"] = float(data.get("confidence", 0.0))

        return validated
    except Exception as e:
        logger.error(f"Error validating transporter data: {str(e)}")
        return data


def _deduplicate_compounds(compounds: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Remove duplicate compounds based on InChI Key."""
    try:
        unique = {}
        for compound in compounds:
            # Generate InChI Key if SMILES is present
            if compound.get("smiles"):
                mol = Chem.MolFromSmiles(compound["smiles"])
                if mol:
                    inchi_key = Chem.MolToInchiKey(mol)
                    if inchi_key not in unique:
                        unique[inchi_key] = compound
                    else:
                        # Merge data from duplicate
                        for key, value in compound.items():
                            if value and not unique[inchi_key].get(key):
                                unique[inchi_key][key] = value

        return list(unique.values())
    except Exception as e:
        logger.error(f"Error deduplicating compounds: {str(e)}")
        return compounds


# Initialize models flag
MODELS_LOADED = False


def initialize_models():
    """Initialize required models."""
    global MODELS_LOADED, ner_model, analyzer, relevance_model, nootropic_model

    try:
        if MODELS_LOADED:
            return True

        # Use BERT-based model for chemical entity recognition
        ner_model = pipeline(
            "token-classification",
            model="allenai/scibert_scivocab_uncased",
            aggregation_strategy="simple",
            from_pt=True,  # Use PyTorch weights
        )

        # Use SciBERT for text analysis
        analyzer = SentenceTransformer(
            "allenai/scibert_scivocab_uncased",
            from_pt=True,  # Use PyTorch weights
        )

        # Use domain-specific model for relevance scoring
        relevance_model = SentenceTransformer(
            "pritamdeka/S-PubMedBert-MS-MARCO",
            from_pt=True,  # Use PyTorch weights
        )

        # Use PubMedBERT for nootropic analysis
        nootropic_model = SentenceTransformer(
            "microsoft/BiomedNLP-PubMedBert-base-uncased-abstract",
            from_pt=True,  # Use PyTorch weights
        )

        MODELS_LOADED = True
        return True
    except Exception as e:
        logging.error(f"Error loading models: {str(e)}")
        MODELS_LOADED = False
        return False


# Do not initialize models at module import
# initialize_models()


# Amino acid transporter patterns
TRANSPORTER_PATTERNS = {
    "lat1": {
        "names": {"lat1", "lat-1", "l-type amino acid transporter 1", "slc7a5", "solute carrier family 7 member 5"},
        "substrates": {"leucine", "isoleucine", "valine", "phenylalanine", "tyrosine", "tryptophan", "methionine"},
    },
    "lat2": {
        "names": {"lat2", "lat-2", "l-type amino acid transporter 2", "slc7a8", "solute carrier family 7 member 8"},
        "substrates": {"alanine", "serine", "cysteine", "threonine", "asparagine", "glutamine", "histidine"},
    },
    "asct1": {
        "names": {"asct1", "asct-1", "alanine serine cysteine transporter 1", "slc1a4", "solute carrier family 1 member 4"},
        "substrates": {"alanine", "serine", "cysteine", "threonine"},
    },
    "asct2": {
        "names": {"asct2", "asct-2", "alanine serine cysteine transporter 2", "slc1a5", "solute carrier family 1 member 5"},
        "substrates": {"alanine", "serine", "cysteine", "threonine", "glutamine", "asparagine"},
    },
}


@dataclass
class Effect:
    """Schema for extracted effect data."""

    name: str
    description: str
    category: str
    confidence: float
    evidence: List[str]
    metadata: Dict[str, Any]


@dataclass
class Mechanism:
    """Schema for extracted mechanism data."""

    name: str
    description: str
    category: str
    confidence: float
    evidence: List[str]
    metadata: Dict[str, Any]


@dataclass
class SafetyData:
    """Schema for extracted safety data."""

    warning: str
    severity: str
    category: str
    confidence: float
    evidence: List[str]
    metadata: Dict[str, Any]


@dataclass
class CommunityInsight:
    """Schema for extracted community insight data."""

    insight: str
    sentiment: str
    confidence: float
    sources: List[str]
    metadata: Dict[str, Any]


@dataclass
class ExtractionResult:
    """Result of LLM extraction."""

    text: str
    confidence: float
    metadata: Dict[str, Any]
    timestamp: datetime = field(default_factory=datetime.now)


def extract_chemical_info(text: str) -> List[Dict[str, str]]:
    """Extract chemical compound information from text.

    Args:
        text: Text to extract chemical information from

    Returns:
        List of dictionaries containing extracted chemical information with keys:
        - name: Chemical name
        - smiles: SMILES string if found
        - inchi: InChI string if found
        - cas: CAS number if found
        - description: Any descriptive text about the compound
    """
    try:
        if not MODELS_LOADED:
            logger.warning("Models not loaded, returning empty results")
            return []

        # Extract chemical entities
        entities = ner_model(text)
        chemicals = [e["word"] for e in entities if e["entity_group"] in ["CHEMICAL", "COMPOUND"]]

        # Extract SMILES/InChI patterns
        smiles_pattern = r"(?:SMILES|smiles)[:=]\s*([^\s;]+)"
        inchi_pattern = r"(?:InChI|inchi)[:=]\s*([^\s;]+)"
        cas_pattern = r"\b\d{1,7}-\d{2}-\d\b"

        smiles = re.findall(smiles_pattern, text)
        inchi = re.findall(inchi_pattern, text)
        cas_numbers = re.findall(cas_pattern, text)

        results = []
        for chem in chemicals:
            result = {"name": chem, "smiles": "", "inchi": "", "cas": "", "description": ""}

            # Find description context
            sentences = re.split(r"[.!?]+", text)
            for sentence in sentences:
                if chem in sentence:
                    result["description"] = sentence.strip()
                    break

            # Add identifiers if found
            if smiles:
                result["smiles"] = smiles[0]
            if inchi:
                result["inchi"] = inchi[0]
            if cas_numbers:
                result["cas"] = cas_numbers[0]

            results.append(result)

        return results

    except Exception as e:
        logger.error(f"Error extracting chemical info: {str(e)}")
        return []


def analyze_patent_text(text: str) -> Dict[str, List[str]]:
    """Analyze patent text to extract relevant chemical and biological information.

    Args:
        text: Patent text to analyze

    Returns:
        Dictionary containing extracted information with keys:
        - compounds: List of chemical compound names
        - targets: List of biological targets
        - assays: List of assay descriptions
        - activities: List of activity descriptions
    """
    try:
        if not MODELS_LOADED:
            logger.warning("Models not loaded, returning empty results")
            return {"compounds": [], "targets": [], "assays": [], "activities": []}

        # Extract chemical compounds
        compounds = [info["name"] for info in extract_chemical_info(text)]

        # Extract biological targets
        targets = []
        target_patterns = [
            r"(?:receptor|enzyme|protein|transporter)\s+([A-Za-z0-9\-]+)",
            r"([A-Z][A-Za-z0-9\-]+(?:\s+receptor|\s+enzyme|\s+protein|\s+transporter))",
        ]
        for pattern in target_patterns:
            targets.extend(re.findall(pattern, text))

        # Extract assay information
        assay_patterns = [
            r"(?:assay|test|analysis|measurement)\s+(?:of|for|to)\s+([^.]+)",
            r"(?:in\s+vitro|in\s+vivo)\s+(?:assay|test|study)\s+([^.]+)",
        ]
        assays = []
        for pattern in assay_patterns:
            assays.extend(re.findall(pattern, text))

        # Extract activity data
        activity_patterns = [
            r"(?:IC50|EC50|Ki|Kd|affinity).{0,50}(?:\d+\.?\d*)\s*(?:nM|uM|mM)",
            r"(?:inhibition|activation|binding).{0,50}(?:\d+\.?\d*)\s*(?:%|percent)",
        ]
        activities = []
        for pattern in activity_patterns:
            activities.extend(re.findall(pattern, text))

        return {
            "compounds": list(set(compounds)),
            "targets": list(set(targets)),
            "assays": list(set(assays)),
            "activities": list(set(activities)),
        }

    except Exception as e:
        logger.error(f"Error analyzing patent text: {str(e)}")
        return {"compounds": [], "targets": [], "assays": [], "activities": []}


def extract_biological_targets(text: str) -> List[Dict[str, str]]:
    """Extract biological target information from text.

    Args:
        text: Text to extract target information from

    Returns:
        List of dictionaries containing extracted target information with keys:
        - name: Target name
        - type: Target type (e.g. receptor, enzyme)
        - organism: Source organism
        - description: Any descriptive text about the target
    """
    try:
        if not MODELS_LOADED:
            logger.warning("Models not loaded, returning empty results")
            return []

        # Extract targets using NER model
        entities = ner_model(text)
        targets = []

        for entity in entities:
            if entity["entity_group"] in ["PROTEIN", "GENE"]:
                target = {"name": entity["word"], "type": "unknown", "organism": "unknown", "description": ""}

                # Determine target type
                if "receptor" in text.lower():
                    target["type"] = "receptor"
                elif "enzyme" in text.lower():
                    target["type"] = "enzyme"
                elif "transporter" in text.lower():
                    target["type"] = "transporter"
                elif "channel" in text.lower():
                    target["type"] = "channel"

                # Find organism context
                organism_pattern = r"(?:human|rat|mouse|mammalian|bacterial)\s+[A-Za-z0-9\-]+"
                organisms = re.findall(organism_pattern, text)
                if organisms:
                    target["organism"] = organisms[0]

                # Find description context
                sentences = re.split(r"[.!?]+", text)
                for sentence in sentences:
                    if entity["word"] in sentence:
                        target["description"] = sentence.strip()
                        break

                targets.append(target)

        return targets

    except Exception as e:
        logger.error(f"Error extracting biological targets: {str(e)}")
        return []


def analyze_activity_data(text: str) -> List[Dict[str, str]]:
    """Extract and analyze activity data from text.

    Args:
        text: Text containing activity data

    Returns:
        List of dictionaries containing activity information with keys:
        - compound: Compound name/ID
        - target: Target name/ID
        - activity_type: Type of activity (e.g. IC50, Ki)
        - value: Activity value
        - unit: Unit of measurement
    """
    try:
        if not MODELS_LOADED:
            logger.warning("Models not loaded, returning empty results")
            return []

        # Extract compounds and targets
        compounds = extract_chemical_info(text)
        targets = extract_biological_targets(text)

        # Extract activity values
        activity_pattern = r"(IC50|EC50|Ki|Kd|affinity)\s*(?:=|:|\s+of\s+)?\s*(\d+\.?\d*)\s*(nM|uM|mM|M)"
        activities = []

        for match in re.finditer(activity_pattern, text):
            activity = {
                "compound": "",
                "target": "",
                "activity_type": match.group(1),
                "value": match.group(2),
                "unit": match.group(3),
            }

            # Find closest compound and target
            sentence = re.split(r"[.!?]+", text[max(0, match.start() - 200) : match.end() + 200])[0]

            for compound in compounds:
                if compound["name"] in sentence:
                    activity["compound"] = compound["name"]
                    break

            for target in targets:
                if target["name"] in sentence:
                    activity["target"] = target["name"]
                    break

            activities.append(activity)

        return activities

    except Exception as e:
        logger.error(f"Error analyzing activity data: {str(e)}")
        return []


def analyze_content_with_llm(content: str, compound_data: Dict[str, Any], llm_api_key: str) -> Dict[str, Any]:
    """
    Use LLM to analyze webpage content and extract structured data.

    Args:
        content: Webpage text content
        compound_data: Known compound identifiers
        llm_api_key: API key for LLM service

    Returns:
        Dictionary of extracted data
    """
    try:
        # Pre-process content to highlight chemical information
        processed_content = _preprocess_chemical_content(content)

        # Construct enhanced prompt for LLM
        prompt = f"""Analyze this chemical compound webpage content and extract relevant information:

Known identifiers:
Name: {compound_data.get('name', 'Unknown')}
CAS: {compound_data.get('cas', 'Unknown')}
SMILES: {compound_data.get('smiles', 'Unknown')}
InChI: {compound_data.get('inchi', 'Unknown')}

Webpage content:
{processed_content[:4000]}

Please extract the following information if present:
1. Alternative names/synonyms (including systematic names, trade names, and common names)
2. Chemical identifiers:
   - CAS numbers
   - SMILES strings
   - InChI strings
   - Registry numbers
3. Chemical properties:
   - Molecular weight
   - Melting point
   - Boiling point
   - LogP
   - pKa
4. Pharmacological properties:
   - Mechanism of action
   - Target receptors
   - Binding affinities
   - Activity type (agonist/antagonist/inhibitor/allosteric modulator)
   - Efficacy data
5. Transporter interactions:
   - LAT1/LAT2 substrate status
   - Other amino acid transporters
   - Transport kinetics (Km, Vmax)
   - Competition data
6. BBB permeability:
   - Experimental evidence
   - Transport mechanisms
   - Rate/extent of BBB crossing
7. Safety data
   - Toxicity (LD50, etc.)
   - Side effects
   - Handling precautions
   - Safety classifications
8. Regulatory status
   - Scheduling status
   - Legal status by country
   - Control measures
9. References
   - Scientific papers
   - Patents
   - Database entries
   - DOIs

Format the response as a JSON object with these fields."""

        # Make LLM API request
        response = requests.post(
            "https://api.openai.com/v1/chat/completions",
            headers={"Authorization": f"Bearer {llm_api_key}", "Content-Type": "application/json"},
            json={"model": "gpt-4", "messages": [{"role": "user", "content": prompt}], "temperature": 0.7, "max_tokens": 2000},
        )

        if response.status_code == 200:
            result = response.json()
            try:
                extracted_data = json.loads(result["choices"][0]["message"]["content"])
                # Post-process and validate extracted data
                validated_data = _validate_chemical_data(extracted_data)
                return validated_data
            except json.JSONDecodeError:
                logger.error("Failed to parse LLM response as JSON")

    except Exception as e:
        logger.error(f"Error analyzing content with LLM: {str(e)}")

    return {}


def extract_patent_compound(example_text: str, llm_api_key: str, include_synthesis: bool = True) -> Dict[str, Any]:
    """
    Extract compound information from patent using LLM.

    Args:
        example_text: Patent example text
        llm_api_key: API key for LLM service
        include_synthesis: Whether to extract synthesis details

    Returns:
        Dictionary containing extracted identifiers and data
    """
    try:
        # Pre-process patent text
        processed_text = _preprocess_patent_text(example_text)

        # Enhanced prompt for chemical information
        prompt = f"""Extract detailed chemical information from this patent example:

Example text: {processed_text}

Please provide:
1. Chemical identifiers:
   - CAS number
   - SMILES string
   - InChI string
   - Registry numbers
2. Chemical names:
   - Systematic name (IUPAC)
   - Common names
   - Trade names
3. Structure details:
   - Molecular formula
   - Molecular weight
   - Stereochemistry
   - Structural features
4. Physical properties:
   - Melting point
   - Boiling point
   - Solubility
   - Physical form
5. Preparation method:
   - Synthetic route
   - Reaction conditions
   - Yields
   - Purification
6. Analytical data:
   - NMR data
   - Mass spec data
   - Elemental analysis
7. Biological activity:
   - Target receptors
   - Activity values
   - Assay conditions

Format the response as a JSON object with these fields."""

        # Make LLM API request with increased context
        response = requests.post(
            "https://api.openai.com/v1/chat/completions",
            headers={"Authorization": f"Bearer {llm_api_key}", "Content-Type": "application/json"},
            json={"model": "gpt-4", "messages": [{"role": "user", "content": prompt}], "temperature": 0.7, "max_tokens": 2000},
        )

        if response.status_code == 200:
            result = response.json()
            try:
                extracted_data = json.loads(result["choices"][0]["message"]["content"])
                # Validate and standardize extracted data
                validated_data = _validate_patent_data(extracted_data)
                return validated_data
            except json.JSONDecodeError:
                logger.error("Failed to parse LLM response as JSON")

    except Exception as e:
        logger.error(f"Error extracting patent compound: {str(e)}")

    return {}


def extract_transporter_data(content: str, compound_data: Dict[str, Any], llm_api_key: str) -> Dict[str, Any]:
    """Extract amino acid transporter information from text."""
    try:
        # Pre-process content
        processed_content = _preprocess_transporter_content(content)

        # Construct transporter-specific prompt
        prompt = f"""Extract amino acid transporter information from this text:

Compound:
Name: {compound_data.get('name', 'Unknown')}
SMILES: {compound_data.get('smiles', 'Unknown')}

Text content:
{processed_content[:4000]}

Please extract:
1. LAT1 interactions:
   - Substrate evidence
   - Binding affinity
   - Transport kinetics
   - Competition data
2. Other transporters:
   - LAT2 interactions
   - ASCT1/2 interactions
   - Other amino acid transporters
3. Transport mechanisms:
   - Active vs passive
   - Energy dependence
   - pH dependence
4. Structure-activity relationships:
   - Key structural features
   - SAR patterns
   - Optimization strategies
5. Supporting evidence:
   - Experimental methods
   - Literature citations
   - Confidence assessment

Format as JSON with these fields."""

        # Process with LLM
        response = requests.post(
            "https://api.openai.com/v1/chat/completions",
            headers={"Authorization": f"Bearer {llm_api_key}", "Content-Type": "application/json"},
            json={"model": "gpt-4", "messages": [{"role": "user", "content": prompt}], "temperature": 0.7, "max_tokens": 2000},
        )

        if response.status_code == 200:
            result = response.json()
            try:
                extracted_data = json.loads(result["choices"][0]["message"]["content"])
                # Validate transporter data
                validated_data = _validate_transporter_data(extracted_data)
                return validated_data
            except json.JSONDecodeError:
                logger.error("Failed to parse transporter LLM response")

    except Exception as e:
        logger.error(f"Error extracting transporter data: {str(e)}")

    return {}


def process_pdf_content(pdf_content: bytes, llm_api_key: str) -> List[Dict[str, Any]]:
    """
    Process PDF content to extract chemical information.

    Args:
        pdf_content: PDF file content as bytes
        llm_api_key: API key for LLM service

    Returns:
        List of dictionaries containing extracted compound data
    """
    try:
        # Read PDF content
        pdf_file = BytesIO(pdf_content)
        pdf_reader = pypdf.PdfReader(pdf_file)

        # Extract text from each page
        full_text = ""
        for page in pdf_reader.pages:
            full_text += page.extract_text() + "\n"

        # Pre-process text
        processed_text = _preprocess_chemical_content(full_text)
        processed_text = _preprocess_transporter_content(processed_text)

        # Split into sections
        sections = _split_into_sections(processed_text)

        compounds = []
        for section in sections:
            # Extract compounds from each section
            prompt = f"""Extract chemical compound information from this text section:

Text:
{section}

Please identify and extract information about any chemical compounds mentioned, including:
1. Chemical names and identifiers
2. Structural information
3. Properties and characteristics
4. Biological activity
5. Transporter interactions
6. BBB permeability data
7. Synthesis details if present

Format each compound as a JSON object."""

            # Process with LLM
            response = requests.post(
                "https://api.openai.com/v1/chat/completions",
                headers={"Authorization": f"Bearer {llm_api_key}", "Content-Type": "application/json"},
                json={"model": "gpt-4", "messages": [{"role": "user", "content": prompt}], "temperature": 0.7, "max_tokens": 2000},
            )

            if response.status_code == 200:
                result = response.json()
                try:
                    section_compounds = json.loads(result["choices"][0]["message"]["content"])
                    if isinstance(section_compounds, list):
                        compounds.extend(section_compounds)
                    elif isinstance(section_compounds, dict):
                        compounds.append(section_compounds)
                except json.JSONDecodeError:
                    continue

        # Deduplicate and validate compounds
        unique_compounds = _deduplicate_compounds(compounds)
        validated_compounds = [_validate_chemical_data(comp) for comp in unique_compounds]

        return validated_compounds

    except Exception as e:
        logger.error(f"Error processing PDF content: {str(e)}")
        return []


class NootropicLLMExtractor:
    """Enhanced text extraction with LLMs."""

    def __init__(
        self,
        model: str = "gpt-4",
        api_key: Optional[str] = None,
        temperature: float = 0.0,
        max_tokens: int = 1000,
        cache_dir: Optional[str] = None,
    ):
        """Initialize extractor.

        Args:
            model: LLM model to use
            api_key: Optional API key
            temperature: Sampling temperature
            max_tokens: Maximum tokens to generate
            cache_dir: Optional cache directory
        """
        if not MODELS_LOADED:
            raise RuntimeError("Required models not loaded")

        self.model = model
        self.api_key = api_key
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.cache_dir = cache_dir

        # Effect categories
        self.effect_categories = {
            "cognitive": [
                "memory enhancement",
                "focus improvement",
                "learning enhancement",
                "mental clarity",
                "cognitive processing",
            ],
            "mood": [
                "anxiety reduction",
                "mood enhancement",
                "stress reduction",
                "motivation increase",
                "emotional stability",
            ],
            "physical": [
                "neuroprotection",
                "neuroplasticity",
                "brain health",
                "neural repair",
                "neurogenesis",
            ],
        }

        # Mechanism categories
        self.mechanism_categories = {
            "neurotransmitter": [
                "acetylcholine",
                "dopamine",
                "serotonin",
                "glutamate",
                "GABA",
            ],
            "receptor": [
                "nicotinic",
                "muscarinic",
                "AMPA",
                "NMDA",
                "5-HT",
            ],
            "cellular": [
                "neuroplasticity",
                "neurogenesis",
                "neuroprotection",
                "anti-inflammation",
                "antioxidant",
            ],
        }

        # Safety categories
        self.safety_categories = {
            "side_effects": [
                "headache",
                "insomnia",
                "anxiety",
                "nausea",
                "fatigue",
            ],
            "interactions": [
                "drug interactions",
                "supplement interactions",
                "contraindications",
                "precautions",
                "warnings",
            ],
            "tolerance": [
                "tolerance development",
                "dependence risk",
                "withdrawal effects",
                "addiction potential",
                "habituation",
            ],
        }

        # Configure extraction strategy
        self.strategy = LLMExtractionStrategy(
            input_format="html",
            provider=model,
            api_token=api_key,
            temperature=temperature,
            max_tokens=max_tokens,
            chunking={
                "type": "sliding_window",
                "window_size": 1000,
                "step_size": 500,
                "overlap": 0.5,
                "min_chunk_size": 100,
                "preprocessing": {
                    "remove_html": True,
                    "normalize_whitespace": True,
                    "preserve_sentences": True,
                },
                "postprocessing": {
                    "merge_strategy": "union",
                    "dedup_threshold": 0.9,
                    "min_confidence": 0.7,
                },
            },
        )

        # Configure crawler
        self.crawler_config = BrowserConfig(
            browser_type="firefox",
            headless=True,
            viewport_width=1920,
            viewport_height=1080,
            user_agent_mode="random",
            ignore_https_errors=True,
            java_script_enabled=True,
            accept_downloads=False,
            sleep_on_close=False,
            verbose=True,
            wait_for_selectors=[
                ".docsum-content",
                ".authors-list",
                ".abstract-content",
            ],
            stealth_mode=True,
        )

    def extract_effects(self, text: str) -> List[Effect]:
        """Extract nootropic effects from text.

        Args:
            text: Text to analyze

        Returns:
            List of Effect objects
        """
        try:
            if not MODELS_LOADED:
                return []

            effects = []
            sentences = re.split(r"[.!?]+", text)

            for sentence in sentences:
                # Look for effect-related keywords
                if any(
                    keyword in sentence.lower()
                    for keyword in [
                        "effect",
                        "impact",
                        "improve",
                        "enhance",
                        "increase",
                        "decrease",
                        "reduce",
                    ]
                ):
                    # Analyze sentence
                    embedding = nootropic_model.encode(sentence)
                    confidence = float(torch.sigmoid(torch.tensor(embedding.mean())))

                    effects.append(
                        Effect(
                            name=sentence.strip(),
                            description=sentence.strip(),
                            category=self._categorize_effect(sentence),
                            confidence=confidence,
                            evidence=[sentence.strip()],
                            metadata={},
                        )
                    )

            return effects

        except Exception as e:
            logger.error(f"Error extracting effects: {str(e)}")
            return []

    async def extract_effects_async(self, text: str) -> List[Effect]:
        """Extract nootropic effects asynchronously.

        Args:
            text: Text to analyze

        Returns:
            List of Effect objects
        """
        try:
            # Extract with crawl4ai
            async with AsyncWebCrawler() as crawler:
                self.strategy.instruction = """
                Extract nootropic effects from the text:
                - Cognitive effects (memory, focus, etc.)
                - Mood effects (anxiety, depression, etc.)
                - Physical effects (energy, fatigue, etc.)
                - Onset, duration, and intensity
                - Evidence quality and confidence
                """

                result = await crawler.arun(
                    text=text,
                    config=self.crawler_config,
                    extraction_strategy=self.strategy,
                )

                effects = []
                if result.success and result.extracted_content:
                    for effect_data in result.extracted_content:
                        effect = Effect(
                            name=effect_data["name"],
                            description=effect_data["description"],
                            category=effect_data.get("category", "unknown"),
                            confidence=self._calculate_confidence(effect_data),
                            evidence=effect_data.get("evidence", []),
                            metadata={
                                "source": "crawl4ai",
                                "timestamp": datetime.now().isoformat(),
                                "raw_response": effect_data,
                            },
                        )
                        effects.append(effect)

                return sorted(effects, key=lambda x: x.confidence, reverse=True)

        except Exception as e:
            logger.error(f"Error extracting effects asynchronously: {str(e)}")
            return []

    def extract_mechanisms(self, text: str) -> List[Mechanism]:
        """Extract mechanisms of action from text.

        Args:
            text: Text to analyze

        Returns:
            List of Mechanism objects
        """
        try:
            if not MODELS_LOADED:
                return []

            mechanisms = []
            sentences = re.split(r"[.!?]+", text)

            for sentence in sentences:
                # Look for mechanism-related keywords
                if any(
                    keyword in sentence.lower()
                    for keyword in [
                        "mechanism",
                        "pathway",
                        "receptor",
                        "enzyme",
                        "neurotransmitter",
                        "signaling",
                        "modulate",
                    ]
                ):
                    # Analyze sentence
                    embedding = nootropic_model.encode(sentence)
                    confidence = float(torch.sigmoid(torch.tensor(embedding.mean())))

                    mechanisms.append(
                        Mechanism(
                            name=sentence.strip(),
                            description=sentence.strip(),
                            category=self._categorize_mechanism(sentence),
                            confidence=confidence,
                            evidence=[sentence.strip()],
                            metadata={},
                        )
                    )

            return mechanisms

        except Exception as e:
            logger.error(f"Error extracting mechanisms: {str(e)}")
            return []

    async def extract_mechanisms_async(self, text: str) -> List[Mechanism]:
        """Extract mechanisms of action asynchronously.

        Args:
            text: Text to analyze

        Returns:
            List of Mechanism objects
        """
        try:
            # Extract with crawl4ai
            async with AsyncWebCrawler() as crawler:
                self.strategy.instruction = """
                Extract mechanisms of action:
                - Neurotransmitter systems affected
                - Receptor interactions
                - Cellular/molecular mechanisms
                - Pharmacokinetics
                - Evidence quality and confidence
                """

                result = await crawler.arun(
                    text=text,
                    config=self.crawler_config,
                    extraction_strategy=self.strategy,
                )

                mechanisms = []
                if result.success and result.extracted_content:
                    for mech_data in result.extracted_content:
                        mechanism = Mechanism(
                            name=mech_data["name"],
                            description=mech_data["description"],
                            category=mech_data.get("category", "unknown"),
                            confidence=self._calculate_confidence(mech_data),
                            evidence=mech_data.get("evidence", []),
                            metadata={
                                "source": "crawl4ai",
                                "timestamp": datetime.now().isoformat(),
                                "raw_response": mech_data,
                            },
                        )
                        mechanisms.append(mechanism)

                return sorted(mechanisms, key=lambda x: x.confidence, reverse=True)

        except Exception as e:
            logger.error(f"Error extracting mechanisms asynchronously: {str(e)}")
            return []

    def extract_safety_data(self, text: str) -> List[SafetyData]:
        """Extract safety information from text.

        Args:
            text: Text to analyze

        Returns:
            List of SafetyData objects
        """
        try:
            if not MODELS_LOADED:
                return []

            safety_data = []
            sentences = re.split(r"[.!?]+", text)

            for sentence in sentences:
                # Look for safety-related keywords
                if any(
                    keyword in sentence.lower()
                    for keyword in [
                        "safety",
                        "risk",
                        "warning",
                        "caution",
                        "danger",
                        "adverse",
                        "side effect",
                        "toxicity",
                    ]
                ):
                    # Analyze sentence
                    embedding = nootropic_model.encode(sentence)
                    confidence = float(torch.sigmoid(torch.tensor(embedding.mean())))

                    safety_data.append(
                        SafetyData(
                            warning=sentence.strip(),
                            severity=self._assess_severity(sentence),
                            category=self._categorize_safety_issue(sentence),
                            confidence=confidence,
                            evidence=[sentence.strip()],
                            metadata={},
                        )
                    )

            return safety_data

        except Exception as e:
            logger.error(f"Error extracting safety data: {str(e)}")
            return []

    async def extract_safety_data_async(self, text: str) -> List[SafetyData]:
        """Extract safety information asynchronously.

        Args:
            text: Text to analyze

        Returns:
            List of SafetyData objects
        """
        try:
            # Extract with crawl4ai
            async with AsyncWebCrawler() as crawler:
                self.strategy.instruction = """
                Extract safety information:
                - Side effects and adverse reactions
                - Drug interactions
                - Contraindications
                - Risk factors
                - Long-term effects
                - Tolerance and dependence
                - Evidence quality and confidence
                """

                result = await crawler.arun(
                    text=text,
                    config=self.crawler_config,
                    extraction_strategy=self.strategy,
                )

                safety_items = []
                if result.success and result.extracted_content:
                    for safety_data in result.extracted_content:
                        safety = SafetyData(
                            warning=safety_data["warning"],
                            severity=safety_data["severity"],
                            category=safety_data.get("category", "unknown"),
                            confidence=self._calculate_confidence(safety_data),
                            evidence=safety_data.get("evidence", []),
                            metadata={
                                "source": "crawl4ai",
                                "timestamp": datetime.now().isoformat(),
                                "raw_response": safety_data,
                            },
                        )
                        safety_items.append(safety)

                return sorted(safety_items, key=lambda x: x.confidence, reverse=True)

        except Exception as e:
            logger.error(f"Error extracting safety data asynchronously: {str(e)}")
            return []

    def extract_community_insights(self, text: str) -> List[CommunityInsight]:
        """Extract community insights from text.

        Args:
            text: Text to analyze

        Returns:
            List of CommunityInsight objects
        """
        try:
            if not MODELS_LOADED:
                return []

            insights = []
            sentences = re.split(r"[.!?]+", text)

            for sentence in sentences:
                # Look for experience-related keywords
                if any(
                    keyword in sentence.lower()
                    for keyword in [
                        "experience",
                        "report",
                        "anecdotal",
                        "user",
                        "community",
                        "forum",
                        "discussion",
                    ]
                ):
                    # Analyze sentence
                    embedding = nootropic_model.encode(sentence)
                    confidence = float(torch.sigmoid(torch.tensor(embedding.mean())))

                    insights.append(
                        CommunityInsight(
                            insight=sentence.strip(),
                            sentiment=self._analyze_sentiment(sentence),
                            confidence=confidence,
                            sources=[],
                            metadata={},
                        )
                    )

            return insights

        except Exception as e:
            logger.error(f"Error extracting community insights: {str(e)}")
            return []

    async def extract_community_insights_async(self, text: str) -> List[CommunityInsight]:
        """Extract community insights asynchronously.

        Args:
            text: Text to analyze

        Returns:
            List of CommunityInsight objects
        """
        try:
            # Extract with crawl4ai
            async with AsyncWebCrawler() as crawler:
                self.strategy.instruction = """
                Extract community insights:
                - User experiences and reports
                - Dosage patterns
                - Common combinations
                - Risk mitigation strategies
                - Perceived benefits and drawbacks
                - Source credibility assessment
                """

                result = await crawler.arun(
                    text=text,
                    config=self.crawler_config,
                    extraction_strategy=self.strategy,
                )

                insights = []
                if result.success and result.extracted_content:
                    for insight_data in result.extracted_content:
                        insight = CommunityInsight(
                            insight=insight_data["insight"],
                            sentiment=insight_data["sentiment"],
                            confidence=self._calculate_confidence(insight_data),
                            sources=insight_data.get("sources", []),
                            metadata={
                                "source": "crawl4ai",
                                "timestamp": datetime.now().isoformat(),
                                "raw_response": insight_data,
                            },
                        )
                        insights.append(insight)

                return sorted(insights, key=lambda x: x.confidence, reverse=True)

        except Exception as e:
            logger.error(f"Error extracting community insights asynchronously: {str(e)}")
            return []

    def _extract_evidence(self, text: str, keyword: str) -> Optional[str]:
        """Extract supporting evidence from text.

        Args:
            text: Source text
            keyword: Keyword to find evidence for

        Returns:
            Extracted evidence text or None
        """
        # Find sentences containing the keyword
        sentences = re.split(r"[.!?]+", text)
        relevant = [s for s in sentences if keyword.lower() in s.lower()]

        if not relevant:
            return None

        # Return most relevant sentence
        return max(relevant, key=len).strip()

    def _calculate_confidence(self, data: Dict[str, Any]) -> float:
        """Calculate confidence score for extracted data.

        Args:
            data: Extracted data dictionary

        Returns:
            Confidence score between 0 and 1
        """
        score = 0.0
        total = 0

        # Check required fields presence
        required_fields = ["name", "description"] if "name" in data else ["warning", "severity"]
        for field in required_fields:
            total += 1
            if data.get(field):
                score += 1.0

        # Check evidence quality
        if data.get("evidence"):
            score += len(data["evidence"]) * 0.2
            total += 1

        # Check source credibility
        if data.get("sources"):
            score += len(data["sources"]) * 0.2
            total += 1

        # Calculate final score
        return min(score / total if total > 0 else 0.0, 1.0)

    def _categorize_effect(self, text: str) -> str:
        """Categorize the type of effect described in text."""
        categories = {
            "cognitive": ["memory", "focus", "attention", "learning", "cognition"],
            "mood": ["mood", "anxiety", "depression", "stress", "emotion"],
            "energy": ["energy", "fatigue", "alertness", "wakefulness"],
            "physical": ["strength", "endurance", "performance", "motor"],
        }

        text_lower = text.lower()
        for category, keywords in categories.items():
            if any(keyword in text_lower for keyword in keywords):
                return category

        return "other"

    def _categorize_mechanism(self, text: str) -> str:
        """Categorize the type of mechanism described in text."""
        categories = {
            "receptor": ["receptor", "binding", "agonist", "antagonist"],
            "enzyme": ["enzyme", "inhibitor", "substrate", "metabolism"],
            "neurotransmitter": ["neurotransmitter", "serotonin", "dopamine", "acetylcholine"],
            "signaling": ["pathway", "cascade", "signaling", "transduction"],
        }

        text_lower = text.lower()
        for category, keywords in categories.items():
            if any(keyword in text_lower for keyword in keywords):
                return category

        return "other"

    def _assess_severity(self, text: str) -> str:
        """Assess the severity level of a safety warning."""
        severity_levels = {
            "high": ["severe", "serious", "dangerous", "fatal", "death"],
            "medium": ["moderate", "significant", "concerning", "caution"],
            "low": ["mild", "minor", "slight", "minimal"],
        }

        text_lower = text.lower()
        for level, keywords in severity_levels.items():
            if any(keyword in text_lower for keyword in keywords):
                return level

        return "unknown"

    def _categorize_safety_issue(self, text: str) -> str:
        """Categorize the type of safety issue described in text."""
        categories = {
            "side_effect": ["side effect", "adverse", "reaction"],
            "interaction": ["interaction", "contraindication", "combine"],
            "toxicity": ["toxic", "poisoning", "overdose"],
            "dependency": ["dependence", "addiction", "withdrawal"],
        }

        text_lower = text.lower()
        for category, keywords in categories.items():
            if any(keyword in text_lower for keyword in keywords):
                return category

        return "other"

    def _analyze_sentiment(self, text: str) -> str:
        """Analyze the sentiment of text."""
        positive_words = ["positive", "good", "effective", "helpful", "beneficial", "success"]
        negative_words = ["negative", "bad", "ineffective", "harmful", "adverse", "failure"]

        text_lower = text.lower()
        positive_count = sum(1 for word in positive_words if word in text_lower)
        negative_count = sum(1 for word in negative_words if word in text_lower)

        if positive_count > negative_count:
            return "positive"
        elif negative_count > positive_count:
            return "negative"
        return "neutral"


class LLMProcessor:
    """Base class for LLM-based text processing."""

    def __init__(
        self,
        model: str = "gpt-4",
        api_key: Optional[str] = None,
        temperature: float = 0.0,
        max_tokens: int = 1000,
        cache_dir: Optional[str] = None,
        device: Optional[str] = None,
    ):
        """Initialize LLM processor.

        Args:
            model: LLM model to use
            api_key: Optional API key
            temperature: Sampling temperature
            max_tokens: Maximum tokens to generate
            cache_dir: Optional cache directory
            device: Optional device to use (cpu, cuda, mps)
        """
        # Initialize models if not already loaded
        if not MODELS_LOADED:
            if not initialize_models():
                raise RuntimeError("Failed to initialize required models")

        self.model = model
        self.api_key = api_key
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.cache_dir = cache_dir

        # Set device
        if device:
            self.device = device
        else:
            if torch.backends.mps.is_available():
                self.device = "mps"
            elif torch.cuda.is_available():
                self.device = "cuda"
            else:
                self.device = "cpu"

        # Configure extraction strategy
        self.strategy = LLMExtractionStrategy(
            input_format="text",
            provider=model,
            api_token=api_key,
            temperature=temperature,
            max_tokens=max_tokens,
            chunking={
                "type": "sliding_window",
                "window_size": 1000,
                "step_size": 500,
                "overlap": 0.5,
                "min_chunk_size": 100,
            },
        )


if __name__ == "__main__":
    # Initialize models when module is run directly
    initialize_models()
