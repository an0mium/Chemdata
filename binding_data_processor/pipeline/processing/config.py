"""Processing pipeline configuration.

This module provides configuration classes for:
1. Target pattern matching
2. Web source integration
3. Social media monitoring
4. Processing settings
"""

from typing import Dict, List
import re


class TargetPatterns:
    """Target patterns for filtering compounds."""
    
    # Serotonin system
    SEROTONIN_PATTERNS = [
        r"5-HT\d*[A-Z]?",
        r"serotonin",
        r"SLC6A4",
        r"HTR[12][A-Z]",
        r"5-hydroxytryptamine",
        r"tryptamine",
    ]
    
    # Dopamine system
    DOPAMINE_PATTERNS = [
        r"D\d+",
        r"dopamine",
        r"DAT",
        r"SLC6A3",
        r"DRD[1-5]",
        r"dopaminergic",
        r"catecholamine",
    ]
    
    # Norepinephrine system
    NOREPINEPHRINE_PATTERNS = [
        r"norepinephrine",
        r"NET",
        r"SLC6A2",
        r"ADRA\d[A-Z]",
        r"adrenergic",
        r"noradrenergic",
    ]
    
    # GABA system
    GABA_PATTERNS = [
        r"GABA[A-Z]?\d*",
        r"SLC6A1",
        r"GABR[A-Z]\d",
        r"gamma-aminobutyric",
        r"benzodiazepine",
    ]
    
    # Glutamate system
    GLUTAMATE_PATTERNS = [
        r"glutamate",
        r"NMDA",
        r"AMPA",
        r"mGluR",
        r"GluR",
        r"GRM\d",
        r"GRIN[12][A-Z]",
        r"kainate",
    ]
    
    # Opioid system
    OPIOID_PATTERNS = [
        r"[μμκδ]?-?opioid",
        r"MOR",
        r"KOR",
        r"DOR",
        r"OPRM1",
        r"OPRK1",
        r"OPRD1",
        r"endorphin",
    ]
    
    # Cannabinoid system
    CANNABINOID_PATTERNS = [
        r"cannabinoid",
        r"CB[12]",
        r"CNR[12]",
        r"endocannabinoid",
    ]
    
    # Psychedelic-related
    PSYCHEDELIC_PATTERNS = [
        r"psychedelic",
        r"hallucinogen",
        r"entheogen",
        r"5-HT2A",
        r"HTR2A",
        r"DMT",
        r"LSD",
        r"psilocybin",
        r"mescaline",
        r"ayahuasca",
        r"ibogaine",
    ]
    
    # Nootropic-related
    NOOTROPIC_PATTERNS = [
        r"nootropic",
        r"cognitive enhancer",
        r"smart drug",
        r"acetylcholine",
        r"nicotinic",
        r"muscarinic",
        r"CHRN[A-Z]\d",
        r"racetam",
        r"ampakine",
        r"eugeroic",
    ]
    
    # Dissociative-related
    DISSOCIATIVE_PATTERNS = [
        r"dissociative",
        r"NMDA",
        r"ketamine",
        r"PCP",
        r"GRIN[12][A-Z]",
        r"glutamate",
        r"arylcyclohexylamine",
    ]
    
    # Stimulant-related
    STIMULANT_PATTERNS = [
        r"stimulant",
        r"amphetamine",
        r"cocaine",
        r"methylphenidate",
        r"DAT",
        r"NET",
        r"SERT",
        r"monoamine",
        r"cathinone",
    ]
    
    @classmethod
    def get_all_patterns(cls) -> Dict[str, str]:
        """Get all target patterns.
        
        Returns:
            Dictionary mapping target types to regex patterns
        """
        return {
            "serotonin": "|".join(cls.SEROTONIN_PATTERNS),
            "dopamine": "|".join(cls.DOPAMINE_PATTERNS),
            "norepinephrine": "|".join(cls.NOREPINEPHRINE_PATTERNS),
            "gaba": "|".join(cls.GABA_PATTERNS),
            "glutamate": "|".join(cls.GLUTAMATE_PATTERNS),
            "opioid": "|".join(cls.OPIOID_PATTERNS),
            "cannabinoid": "|".join(cls.CANNABINOID_PATTERNS),
            "psychedelic": "|".join(cls.PSYCHEDELIC_PATTERNS),
            "nootropic": "|".join(cls.NOOTROPIC_PATTERNS),
            "dissociative": "|".join(cls.DISSOCIATIVE_PATTERNS),
            "stimulant": "|".join(cls.STIMULANT_PATTERNS),
        }


class WebSources:
    """Web sources for compound data."""
    
    # Primary databases
    DATABASES = [
        "pubchem",
        "chembl",
        "drugbank",
        "bindingdb",
    ]
    
    # Community sources
    COMMUNITY = [
        "wikipedia",
        "erowid",
        "psychonautwiki",
        "tripsit",
        "bluelight",
    ]
    
    # Scientific sources
    SCIENTIFIC = [
        "pubmed",
        "patents",
        "thesis",
        "google_scholar",
    ]
    
    # Regulatory sources
    REGULATORY = [
        "fda",
        "ema",
        "who",
        "unscheduler",
    ]
    
    @classmethod
    def get_all_sources(cls) -> List[str]:
        """Get all web sources.
        
        Returns:
            List of all web sources
        """
        return (
            cls.DATABASES +
            cls.COMMUNITY +
            cls.SCIENTIFIC +
            cls.REGULATORY
        )


class SocialMedia:
    """Social media sources for compound monitoring."""
    
    # Reddit communities
    SUBREDDITS = [
        "researchchemicals",
        "DrugNerds",
        "Nootropics",
        "psychopharmacology",
        "Drugs",
        "pharmacology",
        "chemistry",
        "MedicalChem",
        "neuroscience",
    ]
    
    # Twitter search queries
    TWITTER_QUERIES = [
        "research chemical",
        "novel compound",
        "new synthesis",
        "receptor binding",
        "pharmacology",
        "drug discovery",
        "medicinal chemistry",
        "psychoactive",
        "nootropic",
    ]
    
    # Bluesky search queries
    BLUESKY_QUERIES = TWITTER_QUERIES
    
    # Discord servers
    DISCORD_SERVERS = [
        "Chemistry",
        "Pharmacology",
        "Drug Discovery",
        "Medicinal Chemistry",
        "Neuroscience",
    ]
