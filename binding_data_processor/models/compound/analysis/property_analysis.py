"""Property analysis functionality for compound data.

This module provides the PropertyAnalysisMixin class that adds property analysis capabilities:
- Physicochemical property analysis
- Drug-likeness analysis
- ADMET predictions
- Property alerts
- Bioavailability analysis
"""

from dataclasses import dataclass, field
from typing import Dict, List, Tuple


@dataclass
class PropertyAnalysisMixin:
    """Mixin class adding property analysis capabilities."""

    _property_analysis: Dict = field(default_factory=dict)

    def analyze_properties(self) -> Dict:
        """Analyze physicochemical properties and drug-likeness."""
        analysis = {
            "physicochemical": self._analyze_physicochemical(),
            "drug_likeness": self._analyze_drug_likeness(),
            "admet_predictions": self._analyze_admet(),
            "property_alerts": self._analyze_property_alerts(),
            "bioavailability": self._analyze_bioavailability(),
        }
        self._property_analysis = analysis
        return analysis

    def _analyze_physicochemical(self) -> Dict:
        """Analyze physicochemical properties."""
        return {
            "molecular_weight": self.molecular_weight,
            "logp": self.logp,
            "hbd": self.hbd,
            "hba": self.hba,
            "tpsa": self.tpsa,
            "rotatable_bonds": self.rotatable_bonds,
            "charge": self.charge,
            "stereocenter_count": self.stereocenter_count,
            "ring_count": self.ring_count,
        }

    def _analyze_drug_likeness(self) -> Dict:
        """Analyze drug-likeness."""
        # Lipinski's Rule of 5
        lipinski = {
            "mw_ok": self.molecular_weight <= 500,
            "logp_ok": self.logp <= 5,
            "hbd_ok": self.hbd <= 5,
            "hba_ok": self.hba <= 10,
            "violations": 0,
        }
        
        lipinski["violations"] = sum(
            1 for ok in [
                lipinski["mw_ok"],
                lipinski["logp_ok"],
                lipinski["hbd_ok"],
                lipinski["hba_ok"],
            ]
            if not ok
        )

        # Veber's rules
        veber = {
            "rotatable_ok": self.rotatable_bonds <= 10,
            "tpsa_ok": self.tpsa <= 140,
        }

        # Ghose filter
        ghose = {
            "mw_ok": 160 <= self.molecular_weight <= 480,
            "logp_ok": -0.4 <= self.logp <= 5.6,
            "atom_count_ok": 20 <= self.atom_count <= 70,
        }

        # Muegge filter
        muegge = {
            "mw_ok": 200 <= self.molecular_weight <= 600,
            "logp_ok": -2 <= self.logp <= 5,
            "tpsa_ok": self.tpsa <= 150,
            "ring_count_ok": self.ring_count <= 7,
            "hba_ok": self.hba <= 10,
            "hbd_ok": self.hbd <= 5,
        }

        return {
            "lipinski": lipinski,
            "veber": veber,
            "ghose": ghose,
            "muegge": muegge,
            "overall": self._calculate_drug_likeness_score(),
        }

    def _calculate_drug_likeness_score(self) -> Dict:
        """Calculate overall drug-likeness score."""
        score = 1.0
        rules_passed = 0
        total_rules = 4  # Lipinski, Veber, Ghose, Muegge

        # Lipinski contribution
        if self._check_lipinski():
            score *= 1.0
            rules_passed += 1
        else:
            score *= 0.7

        # Veber contribution
        if self._check_veber():
            score *= 1.0
            rules_passed += 1
        else:
            score *= 0.8

        # Ghose contribution
        if self._check_ghose():
            score *= 1.0
            rules_passed += 1
        else:
            score *= 0.8

        # Muegge contribution
        if self._check_muegge():
            score *= 1.0
            rules_passed += 1
        else:
            score *= 0.8

        # Classification
        if rules_passed == total_rules:
            classification = "highly drug-like"
        elif rules_passed >= total_rules * 0.75:
            classification = "moderately drug-like"
        elif rules_passed >= total_rules * 0.5:
            classification = "weakly drug-like"
        else:
            classification = "non-drug-like"

        return {
            "score": score,
            "classification": classification,
            "rules_passed": rules_passed,
            "total_rules": total_rules,
        }

    def _check_lipinski(self) -> bool:
        """Check Lipinski's Rule of 5."""
        return all([
            self.molecular_weight <= 500,
            self.logp <= 5,
            self.hbd <= 5,
            self.hba <= 10,
        ])

    def _check_veber(self) -> bool:
        """Check Veber's rules."""
        return all([
            self.rotatable_bonds <= 10,
            self.tpsa <= 140,
        ])

    def _check_ghose(self) -> bool:
        """Check Ghose filter."""
        return all([
            160 <= self.molecular_weight <= 480,
            -0.4 <= self.logp <= 5.6,
            20 <= self.atom_count <= 70,
        ])

    def _check_muegge(self) -> bool:
        """Check Muegge filter."""
        return all([
            200 <= self.molecular_weight <= 600,
            -2 <= self.logp <= 5,
            self.tpsa <= 150,
            self.ring_count <= 7,
            self.hba <= 10,
            self.hbd <= 5,
        ])

    def _analyze_admet(self) -> Dict:
        """Analyze ADMET predictions."""
        if not hasattr(self, "adme_properties"):
            return {}

        return {
            "absorption": self.adme_properties.get("absorption", {}),
            "distribution": self.adme_properties.get("distribution", {}),
            "metabolism": self.adme_properties.get("metabolism", {}),
            "excretion": self.adme_properties.get("excretion", {}),
            "toxicity": self.adme_properties.get("toxicity", {}),
        }

    def _analyze_property_alerts(self) -> List[Dict]:
        """Analyze property-based alerts."""
        alerts = []

        # Molecular weight alerts
        if self.molecular_weight > 500:
            alerts.append({
                "type": "molecular_weight",
                "description": "High molecular weight may reduce bioavailability",
                "value": self.molecular_weight,
                "threshold": 500,
                "severity": "moderate",
            })

        # LogP alerts
        if self.logp > 5:
            alerts.append({
                "type": "logp",
                "description": "High LogP may cause poor solubility",
                "value": self.logp,
                "threshold": 5,
                "severity": "moderate",
            })
        elif self.logp < -0.4:
            alerts.append({
                "type": "logp",
                "description": "Low LogP may cause poor membrane permeability",
                "value": self.logp,
                "threshold": -0.4,
                "severity": "moderate",
            })

        # TPSA alerts
        if self.tpsa > 140:
            alerts.append({
                "type": "tpsa",
                "description": "High TPSA may reduce membrane permeability",
                "value": self.tpsa,
                "threshold": 140,
                "severity": "moderate",
            })

        # Rotatable bonds alerts
        if self.rotatable_bonds > 10:
            alerts.append({
                "type": "rotatable_bonds",
                "description": "High flexibility may reduce oral bioavailability",
                "value": self.rotatable_bonds,
                "threshold": 10,
                "severity": "moderate",
            })

        return alerts

    def _calculate_bioavailability_score(self) -> float:
        """Calculate basic bioavailability score."""
        score = 1.0

        # Apply penalties based on property thresholds
        if self.molecular_weight > 500:
            score *= 0.8
        if self.logp > 5:
            score *= 0.8
        if self.tpsa > 140:
            score *= 0.8
        if self.rotatable_bonds > 10:
            score *= 0.9
        if self.hbd > 5:
            score *= 0.9
        if self.hba > 10:
            score *= 0.9

        return score

    def _calculate_bbb_score(self) -> Tuple[float, str]:
        """Calculate BBB permeability score and classification."""
        score = 1.0

        # Apply penalties based on BBB-specific thresholds
        if self.molecular_weight > 400:
            score *= 0.8
        if self.logp < 0 or self.logp > 6:
            score *= 0.7
        if self.tpsa > 90:
            score *= 0.6
        if self.rotatable_bonds > 8:
            score *= 0.9
        if self.hbd + self.hba > 8:
            score *= 0.8

        # Classify BBB permeability
        if score >= 0.80:
            classification = "high"
        elif score >= 0.60:
            classification = "moderate"
        else:
            classification = "low"

        return score, classification

    def _classify_bioavailability(self, score: float) -> str:
        """Classify bioavailability based on score."""
        if score >= 0.85:
            return "high"
        elif score >= 0.70:
            return "moderate"
        elif score >= 0.50:
            return "low"
        return "very low"

    def _analyze_bioavailability(self) -> Dict:
        """Analyze predicted bioavailability."""
        # Calculate basic bioavailability score
        score = self._calculate_bioavailability_score()
        classification = self._classify_bioavailability(score)

        # Calculate BBB permeability
        bbb_score, bbb_class = self._calculate_bbb_score()

        # Calculate confidence
        confidence = min(
            1.0,
            0.9 * (1.0 if hasattr(self, "experimental_data") else 0.7)
        )

        return {
            "score": score,
            "classification": classification,
            "limiting_factors": self._get_bioavailability_limiting_factors(),
            "bbb_permeability": {
                "score": bbb_score,
                "classification": bbb_class,
            },
            "predictions": {
                "oral": score >= 0.70,
                "intestinal": score >= 0.60,
                "bbb": bbb_score >= 0.60,
            },
            "confidence": confidence,
        }

    def _get_bioavailability_limiting_factors(self) -> List[str]:
        """Get factors limiting bioavailability."""
        factors = []
        if self.molecular_weight > 500:
            factors.append("high molecular weight")
        if self.logp > 5:
            factors.append("high lipophilicity")
        if self.tpsa > 140:
            factors.append("high polar surface area")
        if self.rotatable_bonds > 10:
            factors.append("high flexibility")
        if self.hbd > 5:
            factors.append("many H-bond donors")
        if self.hba > 10:
            factors.append("many H-bond acceptors")
        return factors
