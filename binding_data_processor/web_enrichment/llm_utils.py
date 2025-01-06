"""LLM utilities for text analysis and extraction.

This module provides enhanced text analysis capabilities using LLMs:
1. Extract chemical information and structures
2. Analyze scientific and patent text
3. Extract effects and mechanisms
4. Extract safety information
5. Extract community insights
6. Validate and cross-reference information
"""

import logging
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from datetime import datetime
import re

from transformers import pipeline
from sentence_transformers import SentenceTransformer
import torch
from crawl4ai import Config


# Initialize models
try:
    # Use BERT-based model for chemical entity recognition
    ner_model = pipeline(
        "token-classification",
        model="allenai/scibert_scivocab_uncased",
        aggregation_strategy="simple",
    )

    # Use SciBERT for text analysis
    analyzer = SentenceTransformer("allenai/scibert_scivocab_uncased")

    # Use domain-specific model for relevance scoring
    relevance_model = SentenceTransformer("pritamdeka/S-PubMedBert-MS-MARCO")

    # Use PubMedBERT for nootropic analysis
    nootropic_model = SentenceTransformer("microsoft/BiomedNLP-PubMedBert-base-uncased-abstract")

    MODELS_LOADED = True
except Exception as e:
    logging.error(f"Error loading models: {str(e)}")
    MODELS_LOADED = False


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


class NootropicLLMExtractor:
    """Enhanced text extraction with LLMs."""

    def __init__(
        self,
        llm_provider: str = "ollama/llama2",
        api_token: Optional[str] = None,
        **kwargs,
    ):
        """Initialize extractor.

        Args:
            llm_provider: LLM provider for text extraction
            api_token: Optional API token for LLM provider
            **kwargs: Additional arguments
        """
        if not MODELS_LOADED:
            raise RuntimeError("Required models not loaded")

        self.llm_provider = llm_provider
        self.api_token = api_token

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

        # Configure LLM with specialized prompts
        self.config = Config.LLM(
            provider=llm_provider,
            api_token=api_token,
            prompts={
                "effects": """
                Extract nootropic effects from the text:
                - Cognitive effects (memory, focus, etc.)
                - Mood effects (anxiety, depression, etc.)
                - Physical effects (energy, fatigue, etc.)
                - Onset, duration, and intensity
                - Evidence quality and confidence
                """,
                "mechanisms": """
                Extract mechanisms of action:
                - Neurotransmitter systems affected
                - Receptor interactions
                - Cellular/molecular mechanisms
                - Pharmacokinetics
                - Evidence quality and confidence
                """,
                "safety": """
                Extract safety information:
                - Side effects and adverse reactions
                - Drug interactions
                - Contraindications
                - Risk factors
                - Long-term effects
                - Tolerance and dependence
                - Evidence quality and confidence
                """,
                "community": """
                Extract community insights:
                - User experiences and reports
                - Dosage patterns
                - Common combinations
                - Risk mitigation strategies
                - Perceived benefits and drawbacks
                - Source credibility assessment
                """,
            },
        )

    def extract_chemical_info(
        self,
        description: str,
        claims: str,
    ) -> Dict[str, Any]:
        """Extract chemical information from text.

        Args:
            description: Description text
            claims: Claims text

        Returns:
            Dictionary containing extracted chemical information
        """
        try:
            # Combine text
            text = f"{description}\n{claims}"

            # Extract chemical entities
            entities = ner_model(text)
            chemicals = [e["word"] for e in entities if e["entity_group"] in ["CHEMICAL", "COMPOUND"]]

            # Extract SMILES/InChI patterns
            smiles_pattern = r"(?:SMILES|smiles)[:=]\s*([^\s;]+)"
            inchi_pattern = r"(?:InChI|inchi)[:=]\s*([^\s;]+)"

            smiles = re.findall(smiles_pattern, text)
            inchi = re.findall(inchi_pattern, text)

            # Extract molecular weights
            mw_pattern = r"(?:molecular weight|MW|mw)[:=]\s*([\d.]+)"
            weights = re.findall(mw_pattern, text)

            return {
                "chemicals": list(set(chemicals)),
                "smiles": list(set(smiles)),
                "inchi": list(set(inchi)),
                "molecular_weights": [float(w) for w in weights],
            }

        except Exception as e:
            logging.error(f"Error extracting chemical info: {str(e)}")
            return {}

    def analyze_text(
        self,
        title: str,
        abstract: str,
        content: str,
    ) -> Dict[str, Any]:
        """Analyze text for relevance and insights.

        Args:
            title: Title text
            abstract: Abstract text
            content: Main content text

        Returns:
            Dictionary containing text analysis results
        """
        try:
            # Combine text sections
            text = f"{title}\n{abstract}\n{content}"

            # Generate embeddings
            text_embedding = analyzer.encode(text)

            # Calculate relevance scores for key topics
            topics = [
                "psychoactive compounds",
                "receptor binding",
                "brain penetration",
                "pharmacological activity",
                "toxicity",
                "abuse potential",
                "therapeutic use",
                "synthesis methods",
                "drug formulation",
            ]

            topic_embeddings = relevance_model.encode(topics)
            relevance_scores = torch.nn.functional.cosine_similarity(
                torch.tensor(text_embedding).unsqueeze(0),
                torch.tensor(topic_embeddings),
                dim=1,
            )

            # Extract key findings
            findings = []

            # Look for activity data
            if re.search(r"(?:IC50|EC50|Ki|Kd|affinity).{0,50}(?:\d+\.?\d*)\s*(?:nM|uM|mM)", text):
                findings.append("Contains binding/activity data")

            # Look for BBB data
            if re.search(r"blood.?brain.?barrier|BBB|brain.?penetration", text):
                findings.append("Contains BBB permeability data")

            # Look for safety data
            if re.search(r"toxicity|safety|adverse|side.?effects", text):
                findings.append("Contains safety/toxicity data")

            # Look for abuse potential
            if re.search(r"abuse|dependence|addiction|withdrawal", text):
                findings.append("Contains abuse potential data")

            # Build scores dictionary
            scores_dict = {}
            for topic, score in zip(topics, relevance_scores):
                scores_dict[topic] = float(score)

            return {
                "relevance_scores": scores_dict,
                "key_findings": findings,
            }

        except Exception as e:
            logging.error(f"Error analyzing text: {str(e)}")
            return {}

    async def extract_effects(self, text: str) -> List[Effect]:
        """Extract nootropic effects.

        Args:
            text: Text to analyze

        Returns:
            List of extracted effects
        """
        try:
            # Extract effects using both LLM and embeddings
            llm_response = await self._extract_with_llm(text, "effects")
            text_embedding = nootropic_model.encode(text)

            effects = []
            for category, keywords in self.effect_categories.items():
                keyword_embeddings = nootropic_model.encode(keywords)
                scores = torch.nn.functional.cosine_similarity(
                    torch.tensor(text_embedding).unsqueeze(0),
                    torch.tensor(keyword_embeddings),
                    dim=1,
                )

                for keyword, score in zip(keywords, scores):
                    if score > 0.5:  # Confidence threshold
                        effect = Effect(
                            name=keyword,
                            description=self._extract_evidence(text, keyword) or "",
                            category=category,
                            confidence=float(score),
                            evidence=(
                                [self._extract_evidence(text, keyword)] if self._extract_evidence(text, keyword) else []
                            ),
                            metadata={
                                "source": "embeddings",
                                "timestamp": datetime.now().isoformat(),
                                "raw_score": float(score),
                            },
                        )
                        effects.append(effect)

            # Add LLM-extracted effects
            for effect_data in llm_response:
                effect = Effect(
                    name=effect_data["name"],
                    description=effect_data["description"],
                    category=effect_data.get("category", "unknown"),
                    confidence=self._calculate_confidence(effect_data),
                    evidence=effect_data.get("evidence", []),
                    metadata={
                        "source": "llm",
                        "timestamp": datetime.now().isoformat(),
                        "raw_response": effect_data,
                    },
                )
                effects.append(effect)

            return sorted(effects, key=lambda x: x.confidence, reverse=True)

        except Exception as e:
            logging.error(f"Error extracting effects: {str(e)}")
            return []

    async def extract_mechanisms(self, text: str) -> List[Mechanism]:
        """Extract mechanisms of action.

        Args:
            text: Text to analyze

        Returns:
            List of extracted mechanisms
        """
        try:
            # Extract mechanisms using both LLM and embeddings
            llm_response = await self._extract_with_llm(text, "mechanisms")
            text_embedding = nootropic_model.encode(text)

            mechanisms = []
            for category, keywords in self.mechanism_categories.items():
                keyword_embeddings = nootropic_model.encode(keywords)
                scores = torch.nn.functional.cosine_similarity(
                    torch.tensor(text_embedding).unsqueeze(0),
                    torch.tensor(keyword_embeddings),
                    dim=1,
                )

                for keyword, score in zip(keywords, scores):
                    if score > 0.5:  # Confidence threshold
                        mechanism = Mechanism(
                            name=keyword,
                            description=self._extract_evidence(text, keyword) or "",
                            category=category,
                            confidence=float(score),
                            evidence=(
                                [self._extract_evidence(text, keyword)] if self._extract_evidence(text, keyword) else []
                            ),
                            metadata={
                                "source": "embeddings",
                                "timestamp": datetime.now().isoformat(),
                                "raw_score": float(score),
                            },
                        )
                        mechanisms.append(mechanism)

            # Add LLM-extracted mechanisms
            for mech_data in llm_response:
                mechanism = Mechanism(
                    name=mech_data["name"],
                    description=mech_data["description"],
                    category=mech_data.get("category", "unknown"),
                    confidence=self._calculate_confidence(mech_data),
                    evidence=mech_data.get("evidence", []),
                    metadata={
                        "source": "llm",
                        "timestamp": datetime.now().isoformat(),
                        "raw_response": mech_data,
                    },
                )
                mechanisms.append(mechanism)

            return sorted(mechanisms, key=lambda x: x.confidence, reverse=True)

        except Exception as e:
            logging.error(f"Error extracting mechanisms: {str(e)}")
            return []

    async def extract_safety(self, text: str) -> List[SafetyData]:
        """Extract safety information.

        Args:
            text: Text to analyze

        Returns:
            List of extracted safety data
        """
        try:
            # Extract safety info using both LLM and embeddings
            llm_response = await self._extract_with_llm(text, "safety")
            text_embedding = nootropic_model.encode(text)

            safety_items = []
            for category, keywords in self.safety_categories.items():
                keyword_embeddings = nootropic_model.encode(keywords)
                scores = torch.nn.functional.cosine_similarity(
                    torch.tensor(text_embedding).unsqueeze(0),
                    torch.tensor(keyword_embeddings),
                    dim=1,
                )

                for keyword, score in zip(keywords, scores):
                    if score > 0.5:  # Confidence threshold
                        safety = SafetyData(
                            warning=keyword,
                            severity="unknown",
                            category=category,
                            confidence=float(score),
                            evidence=(
                                [self._extract_evidence(text, keyword)] if self._extract_evidence(text, keyword) else []
                            ),
                            metadata={
                                "source": "embeddings",
                                "timestamp": datetime.now().isoformat(),
                                "raw_score": float(score),
                            },
                        )
                        safety_items.append(safety)

            # Add LLM-extracted safety data
            for safety_data in llm_response:
                safety = SafetyData(
                    warning=safety_data["warning"],
                    severity=safety_data["severity"],
                    category=safety_data.get("category", "unknown"),
                    confidence=self._calculate_confidence(safety_data),
                    evidence=safety_data.get("evidence", []),
                    metadata={
                        "source": "llm",
                        "timestamp": datetime.now().isoformat(),
                        "raw_response": safety_data,
                    },
                )
                safety_items.append(safety)

            return sorted(safety_items, key=lambda x: x.confidence, reverse=True)

        except Exception as e:
            logging.error(f"Error extracting safety data: {str(e)}")
            return []

    async def extract_community_insights(self, text: str) -> List[CommunityInsight]:
        """Extract community insights.

        Args:
            text: Text to analyze

        Returns:
            List of extracted community insights
        """
        try:
            # Extract community insights using LLM
            response = await self._extract_with_llm(text, "community")

            # Parse and validate insights
            insights = []
            for insight_data in response:
                insight = CommunityInsight(
                    insight=insight_data["insight"],
                    sentiment=insight_data["sentiment"],
                    confidence=self._calculate_confidence(insight_data),
                    sources=insight_data.get("sources", []),
                    metadata={
                        "source": "llm",
                        "timestamp": datetime.now().isoformat(),
                        "raw_response": insight_data,
                    },
                )
                insights.append(insight)

            return sorted(insights, key=lambda x: x.confidence, reverse=True)

        except Exception as e:
            logging.error(f"Error extracting community insights: {str(e)}")
            return []

    async def _extract_with_llm(self, text: str, prompt_key: str) -> List[Dict[str, Any]]:
        """Extract information using LLM.

        Args:
            text: Text to analyze
            prompt_key: Key for prompt template

        Returns:
            List of extracted data dictionaries
        """
        try:
            # Get prompt template
            prompt = self.config.prompts[prompt_key]

            # Call LLM with prompt and text
            response = await self._call_llm(prompt, text)

            # Parse and validate response
            return self._parse_llm_response(response)

        except Exception as e:
            logging.error(f"Error in LLM extraction: {str(e)}")
            return []

    async def _call_llm(self, prompt: str, text: str) -> str:
        """Call LLM with prompt and text.

        Args:
            prompt: Prompt template
            text: Text to analyze

        Returns:
            LLM response text
        """
        # TODO: Implement actual LLM call
        # This is a placeholder that should be replaced with actual LLM integration
        return "[]"

    def _parse_llm_response(self, response: str) -> List[Dict[str, Any]]:
        """Parse LLM response into structured data.

        Args:
            response: Raw LLM response text

        Returns:
            List of parsed data dictionaries
        """
        # TODO: Implement actual response parsing
        # This is a placeholder that should be replaced with actual parsing logic
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
