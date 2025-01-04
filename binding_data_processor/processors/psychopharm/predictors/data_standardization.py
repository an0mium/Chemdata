"""Data standardization for BBB permeability prediction.

This module provides functionality to:
1. Standardize data formats across sources
2. Convert units and measurements
3. Normalize text and categorical data
4. Handle missing and inconsistent values
5. Merge data from multiple sources
"""

import logging
from typing import Dict, List, Optional, Any, Set
import pandas as pd
import numpy as np
from dataclasses import dataclass
from datetime import datetime


@dataclass
class StandardizedData:
    """Container for standardized compound data."""
    
    # Basic info
    compound_id: str
    name: str
    smiles: str
    cas: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    
    # Properties
    molecular_weight: Optional[float] = None
    logp: Optional[float] = None
    psa: Optional[float] = None
    hba: Optional[int] = None
    hbd: Optional[int] = None
    rotatable_bonds: Optional[int] = None
    
    # Activity data
    activities: Dict[str, List[Dict[str, Any]]] = None
    targets: Dict[str, List[Dict[str, Any]]] = None
    mechanisms: Dict[str, List[Dict[str, Any]]] = None
    
    # Safety data
    toxicity: Dict[str, Any] = None
    side_effects: Dict[str, Any] = None
    warnings: List[str] = None
    
    # References
    references: Dict[str, List[str]] = None
    timestamps: Dict[str, datetime] = None
    sources: Set[str] = None


class DataStandardizer:
    """Standardizer for compound data from multiple sources."""

    def __init__(
        self,
        cache_dir: Optional[str] = None,
        log_level: int = logging.INFO,
    ):
        """Initialize data standardizer."""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(log_level)
        self.cache_dir = cache_dir

        # Unit conversions
        self.unit_conversions = {
            "activity": {
                "nm": 1.0,
                "um": 1000.0,
                "mm": 1000000.0,
                "m": 1000000000.0,
            },
            "weight": {
                "g": 1.0,
                "mg": 0.001,
                "ug": 0.000001,
                "ng": 0.000000001,
            },
            "time": {
                "s": 1.0,
                "min": 60.0,
                "h": 3600.0,
                "d": 86400.0,
            },
        }

        # Text normalization patterns
        self.text_patterns = {
            "target": {
                r"5-HT(\d[A-Za-z]?)": r"serotonin_\1",
                r"D(\d)": r"dopamine_\1",
                r"NMDA": "nmda",
                r"AMPA": "ampa",
                r"GABA[_-]?([AB])": r"gaba_\1",
            },
            "activity": {
                r"IC50": "ic50",
                r"EC50": "ec50",
                r"Ki": "ki",
                r"Kd": "kd",
                r"pIC50": "pic50",
                r"pEC50": "pec50",
                r"pKi": "pki",
            },
            "effect": {
                r"antidepress\w+": "antidepressant",
                r"anxiolytic\w*": "anxiolytic",
                r"psychedelic\w*": "psychedelic",
                r"stimulant\w*": "stimulant",
                r"sedativ\w+": "sedative",
            },
        }

    def standardize_compound_data(
        self,
        raw_data: Dict[str, Dict[str, Any]],
        compound_id: str,
        name: str,
        smiles: str,
    ) -> StandardizedData:
        """Standardize compound data from multiple sources."""
        self.logger.debug(f"Standardizing data for {name}")
        
        # Initialize standardized data
        std_data = StandardizedData(
            compound_id=compound_id,
            name=name,
            smiles=smiles,
            activities={},
            targets={},
            mechanisms={},
            toxicity={},
            side_effects={},
            warnings=[],
            references={},
            timestamps={},
            sources=set(),
        )
        
        # Process each source
        for source, data in raw_data.items():
            try:
                self._process_source_data(source, data, std_data)
                std_data.sources.add(source)
            except Exception as e:
                self.logger.error(
                    f"Error processing {source} data for {name}: {str(e)}",
                    exc_info=True
                )
        
        # Post-process standardized data
        self._post_process_data(std_data)
        
        return std_data

    def _process_source_data(
        self,
        source: str,
        data: Dict[str, Any],
        std_data: StandardizedData,
    ) -> None:
        """Process data from a specific source."""
        # Update basic info
        self._update_basic_info(data, std_data)
        
        # Update properties
        self._update_properties(data, std_data)
        
        # Update activity data
        self._update_activities(source, data, std_data)
        
        # Update target data
        self._update_targets(source, data, std_data)
        
        # Update mechanism data
        self._update_mechanisms(source, data, std_data)
        
        # Update safety data
        self._update_safety_data(source, data, std_data)
        
        # Update references
        self._update_references(source, data, std_data)
        
        # Update timestamp
        if "timestamp" in data:
            try:
                std_data.timestamps[source] = pd.to_datetime(
                    data["timestamp"]
                ).to_pydatetime()
            except (ValueError, TypeError):
                self.logger.warning(
                    f"Invalid timestamp in {source} data"
                )

    def _update_basic_info(
        self, data: Dict[str, Any], std_data: StandardizedData
    ) -> None:
        """Update basic compound information."""
        # Update CAS number
        if "cas" in data and not std_data.cas:
            std_data.cas = str(data["cas"])
        elif "cas_number" in data and not std_data.cas:
            std_data.cas = str(data["cas_number"])
        
        # Update InChI
        if "inchi" in data and not std_data.inchi:
            std_data.inchi = str(data["inchi"])
        
        # Update InChIKey
        if "inchikey" in data and not std_data.inchikey:
            std_data.inchikey = str(data["inchikey"])

    def _update_properties(
        self, data: Dict[str, Any], std_data: StandardizedData
    ) -> None:
        """Update compound properties."""
        props = data.get("properties", {})
        
        # Update numeric properties
        self._update_numeric_properties(props, std_data)
        
        # Update integer properties
        self._update_integer_properties(props, std_data)

    def _update_numeric_properties(
        self, props: Dict[str, Any], std_data: StandardizedData
    ) -> None:
        """Update numeric properties."""
        # Update molecular weight
        self._update_float_property(
            props, std_data, "molecular_weight", "molecular_weight"
        )
        
        # Update LogP
        self._update_float_property(
            props, std_data, "logp", "logp"
        )
        
        # Update PSA
        self._update_float_property(
            props, std_data, "psa", "psa"
        )

    def _update_integer_properties(
        self, props: Dict[str, Any], std_data: StandardizedData
    ) -> None:
        """Update integer properties."""
        # Update HBA
        self._update_int_property(
            props, std_data, "hba", "hba"
        )
        
        # Update HBD
        self._update_int_property(
            props, std_data, "hbd", "hbd"
        )
        
        # Update rotatable bonds
        self._update_int_property(
            props, std_data, "rotatable_bonds", "rotatable_bonds"
        )

    def _update_float_property(
        self,
        props: Dict[str, Any],
        std_data: StandardizedData,
        source_key: str,
        target_attr: str,
    ) -> None:
        """Update a float property if not already set."""
        if source_key in props and not getattr(std_data, target_attr):
            try:
                setattr(std_data, target_attr, float(props[source_key]))
            except (ValueError, TypeError):
                pass

    def _update_int_property(
        self,
        props: Dict[str, Any],
        std_data: StandardizedData,
        source_key: str,
        target_attr: str,
    ) -> None:
        """Update an integer property if not already set."""
        if source_key in props and not getattr(std_data, target_attr):
            try:
                setattr(std_data, target_attr, int(props[source_key]))
            except (ValueError, TypeError):
                pass

    def _update_activities(
        self,
        source: str,
        data: Dict[str, Any],
        std_data: StandardizedData,
    ) -> None:
        """Update activity data."""
        activities = data.get("activity_data", [])
        if not isinstance(activities, list):
            activities = [activities]
        
        for activity in activities:
            if not isinstance(activity, dict):
                continue
            
            # Standardize activity type
            activity_type = self._normalize_text(
                activity.get("type", ""),
                "activity",
            )
            if not activity_type:
                continue
            
            # Convert value and unit
            try:
                value = float(activity.get("value", 0))
                unit = activity.get("unit", "")
                std_value = self._convert_unit(value, unit, "activity")
            except (ValueError, TypeError):
                continue
            
            # Create standardized activity entry
            std_activity = {
                "value": std_value,
                "original_value": value,
                "original_unit": unit,
                "source": source,
                "confidence": float(activity.get("confidence", 1.0)),
                "reference": activity.get("reference", ""),
            }
            
            # Add to standardized data
            if activity_type not in std_data.activities:
                std_data.activities[activity_type] = []
            std_data.activities[activity_type].append(std_activity)

    def _update_targets(
        self,
        source: str,
        data: Dict[str, Any],
        std_data: StandardizedData,
    ) -> None:
        """Update target data."""
        targets = data.get("target_data", [])
        if not isinstance(targets, list):
            targets = [targets]
        
        for target in targets:
            if not isinstance(target, dict):
                continue
            
            # Standardize target name
            target_name = self._normalize_text(
                target.get("name", ""),
                "target",
            )
            if not target_name:
                continue
            
            # Create standardized target entry
            std_target = {
                "type": target.get("type", "unknown"),
                "organism": target.get("organism", "human"),
                "confidence": float(target.get("confidence", 1.0)),
                "source": source,
                "reference": target.get("reference", ""),
            }
            
            # Add to standardized data
            if target_name not in std_data.targets:
                std_data.targets[target_name] = []
            std_data.targets[target_name].append(std_target)

    def _update_mechanisms(
        self,
        source: str,
        data: Dict[str, Any],
        std_data: StandardizedData,
    ) -> None:
        """Update mechanism data."""
        mechanisms = data.get("mechanism_data", [])
        if not isinstance(mechanisms, list):
            mechanisms = [mechanisms]
        
        for mechanism in mechanisms:
            if not isinstance(mechanism, dict):
                continue
            
            # Standardize mechanism type
            mech_type = mechanism.get("type", "").lower()
            if not mech_type:
                continue
            
            # Create standardized mechanism entry
            std_mechanism = {
                "description": mechanism.get("description", ""),
                "confidence": float(mechanism.get("confidence", 1.0)),
                "source": source,
                "reference": mechanism.get("reference", ""),
            }
            
            # Add to standardized data
            if mech_type not in std_data.mechanisms:
                std_data.mechanisms[mech_type] = []
            std_data.mechanisms[mech_type].append(std_mechanism)

    def _update_safety_data(
        self,
        source: str,
        data: Dict[str, Any],
        std_data: StandardizedData,
    ) -> None:
        """Update safety data."""
        safety = data.get("safety", {})
        if not isinstance(safety, dict):
            return
        
        # Update each type of safety data
        self._update_toxicity_data(source, safety, std_data)
        self._update_side_effects(source, safety, std_data)
        self._update_warnings(safety, std_data)

    def _update_toxicity_data(
        self,
        source: str,
        safety: Dict[str, Any],
        std_data: StandardizedData,
    ) -> None:
        """Update toxicity data."""
        if "toxicity" not in safety:
            return
            
        tox_data = safety["toxicity"]
        if not isinstance(tox_data, dict):
            return
            
        confidence = safety.get("confidence", 1.0)
        for tox_type, tox_value in tox_data.items():
            if tox_type not in std_data.toxicity:
                std_data.toxicity[tox_type] = {
                    "value": tox_value,
                    "source": source,
                    "confidence": confidence,
                }

    def _update_side_effects(
        self,
        source: str,
        safety: Dict[str, Any],
        std_data: StandardizedData,
    ) -> None:
        """Update side effects data."""
        if "side_effects" not in safety:
            return
            
        effects = safety["side_effects"]
        if not isinstance(effects, dict):
            return
            
        for effect, effect_data in effects.items():
            if effect not in std_data.side_effects:
                std_data.side_effects[effect] = {
                    "frequency": effect_data.get("frequency", "unknown"),
                    "severity": effect_data.get("severity", "unknown"),
                    "source": source,
                }

    def _update_warnings(
        self,
        safety: Dict[str, Any],
        std_data: StandardizedData,
    ) -> None:
        """Update warnings data."""
        if "warnings" in safety:
            warnings = safety["warnings"]
            if isinstance(warnings, list):
                std_data.warnings.extend(warnings)

    def _update_references(
        self,
        source: str,
        data: Dict[str, Any],
        std_data: StandardizedData,
    ) -> None:
        """Update reference data."""
        references = data.get("references", [])
        if isinstance(references, list):
            std_data.references[source] = references

    def _normalize_text(
        self, text: str, pattern_type: str
    ) -> str:
        """Normalize text using defined patterns."""
        if not isinstance(text, str):
            return ""
        
        text = text.strip().lower()
        
        # Apply patterns
        if pattern_type in self.text_patterns:
            for pattern, replacement in self.text_patterns[pattern_type].items():
                text = pd.Series([text]).str.replace(
                    pattern,
                    replacement,
                    regex=True,
                )[0]
        
        return text

    def _convert_unit(
        self,
        value: float,
        unit: str,
        unit_type: str,
    ) -> float:
        """Convert value to standard unit."""
        if unit_type not in self.unit_conversions:
            return value
        
        unit = unit.lower()
        conversions = self.unit_conversions[unit_type]
        
        if unit in conversions:
            return value * conversions[unit]
        
        return value

    def _post_process_data(self, std_data: StandardizedData) -> None:
        """Post-process standardized data."""
        # Sort activities by value
        for activity_type in std_data.activities:
            std_data.activities[activity_type].sort(
                key=lambda x: x["value"]
            )
        
        # Remove duplicate warnings
        if std_data.warnings:
            std_data.warnings = list(set(std_data.warnings))
        
        # Calculate average properties if missing
        self._calculate_missing_properties(std_data)
        
        # Validate data consistency
        self._validate_data_consistency(std_data)

    def _calculate_missing_properties(
        self, std_data: StandardizedData
    ) -> None:
        """Calculate missing properties from available data."""
        # Calculate molecular weight if missing
        if not std_data.molecular_weight and std_data.inchi:
            try:
                from rdkit import Chem
                from rdkit.Chem import Descriptors
                mol = Chem.MolFromInchi(std_data.inchi)
                if mol:
                    std_data.molecular_weight = Descriptors.ExactMolWt(mol)
            except ImportError:
                pass
        
        # Calculate LogP if missing
        if not std_data.logp and std_data.inchi:
            try:
                from rdkit import Chem
                from rdkit.Chem import Crippen
                mol = Chem.MolFromInchi(std_data.inchi)
                if mol:
                    std_data.logp = Crippen.MolLogP(mol)
            except ImportError:
                pass

    def _validate_data_consistency(
        self, std_data: StandardizedData
    ) -> None:
        """Validate consistency of standardized data."""
        # Check activity values
        for activity_type, activities in std_data.activities.items():
            values = [a["value"] for a in activities]
            if len(values) > 1:
                mean = np.mean(values)
                std = np.std(values)
                # Flag outliers
                for activity in activities:
                    if abs(activity["value"] - mean) > 2 * std:
                        activity["is_outlier"] = True
        
        # Check property consistency
        if std_data.molecular_weight and std_data.molecular_weight < 0:
            std_data.molecular_weight = None
        
        if std_data.psa and std_data.psa < 0:
            std_data.psa = None
        
        if std_data.hba and std_data.hba < 0:
            std_data.hba = None
        
        if std_data.hbd and std_data.hbd < 0:
            std_data.hbd = None
        
        if std_data.rotatable_bonds and std_data.rotatable_bonds < 0:
            std_data.rotatable_bonds = None
