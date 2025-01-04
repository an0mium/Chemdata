"""Community data integration functionality."""

from dataclasses import field
from typing import Dict, List, Optional, Set
from datetime import datetime

from .base import DoseRange, TimeRange, EffectScore, RiskLevel


class CommunityDataMixin:
    """Mixin class providing community data integration."""

    # Experience reports
    experience_reports: Dict[str, Dict] = field(default_factory=dict)
    report_sources: Set[str] = field(default_factory=set)
    
    # Dosage information
    dosage_data: Dict[str, DoseRange] = field(default_factory=dict)  # ROA -> DoseRange
    duration_data: Dict[str, TimeRange] = field(default_factory=dict)  # Phase -> TimeRange
    
    # Effect reports
    reported_effects: Dict[str, List[EffectScore]] = field(default_factory=dict)
    effect_frequencies: Dict[str, float] = field(default_factory=dict)
    
    # Safety information
    reported_risks: Dict[str, RiskLevel] = field(default_factory=dict)
    contraindications: Set[str] = field(default_factory=set)
    interactions: Dict[str, Dict] = field(default_factory=dict)
    
    # Usage statistics
    total_reports: int = field(default=0)
    last_report_date: Optional[datetime] = field(default=None)
    popularity_score: float = field(default=0.0)

    def get_community_dict(self) -> Dict:
        """Get dictionary of community data."""
        return {
            "reports": {
                "total": self.total_reports,
                "sources": list(self.report_sources),
                "last_date": self.last_report_date.isoformat() if self.last_report_date else None,
                "popularity": self.popularity_score,
            },
            "dosage": {
                roa: {
                    "min": min_dose,
                    "max": max_dose,
                    "recommended": rec_dose,
                }
                for roa, (min_dose, max_dose, rec_dose) in self.dosage_data.items()
            },
            "duration": {
                phase: {
                    "onset": onset,
                    "duration": duration,
                }
                for phase, (onset, duration) in self.duration_data.items()
            },
            "effects": {
                effect: {
                    "scores": [
                        {"magnitude": mag, "confidence": conf}
                        for mag, conf in scores
                    ],
                    "frequency": self.effect_frequencies.get(effect, 0.0),
                }
                for effect, scores in self.reported_effects.items()
            },
            "safety": {
                "risks": {
                    risk: level.value
                    for risk, level in self.reported_risks.items()
                },
                "contraindications": list(self.contraindications),
                "interactions": self.interactions,
            },
        }

    def _update_report_tracking(
        self,
        source: str,
        report_data: Dict,
        report_date: Optional[datetime] = None
    ) -> str:
        """Update report tracking information."""
        report_id = f"{source}_{len(self.experience_reports)}"
        self.experience_reports[report_id] = report_data
        self.report_sources.add(source)
        self.total_reports += 1
        
        if report_date:
            if not self.last_report_date or report_date > self.last_report_date:
                self.last_report_date = report_date
                
        return report_id

    def _process_dosage_timeline(self, report_data: Dict) -> None:
        """Process dosage and timeline information."""
        if "dosage" in report_data:
            roa = report_data["dosage"].get("route", "oral")
            dose = report_data["dosage"].get("amount")
            if dose:
                self._update_dosage_data(roa, dose)
                
        if "timeline" in report_data:
            for phase, times in report_data["timeline"].items():
                if "onset" in times and "duration" in times:
                    self._update_duration_data(
                        phase,
                        times["onset"],
                        times["duration"]
                    )

    def _process_effects_safety(self, report_data: Dict) -> None:
        """Process effects and safety information."""
        if "effects" in report_data:
            for effect, data in report_data["effects"].items():
                magnitude = data.get("magnitude", 0.0)
                confidence = data.get("confidence", 0.5)
                self._update_effect_data(effect, magnitude, confidence)
                
        if "risks" in report_data:
            for risk, level in report_data["risks"].items():
                self._update_risk_data(risk, level)
                
        if "interactions" in report_data:
            for substance, interaction in report_data["interactions"].items():
                self._update_interaction_data(substance, interaction)

    def add_experience_report(
        self,
        source: str,
        report_data: Dict,
        report_date: Optional[datetime] = None,
    ) -> None:
        """Add a new experience report."""
        # Update report tracking
        self._update_report_tracking(source, report_data, report_date)
        
        # Process dosage and timeline
        self._process_dosage_timeline(report_data)
        
        # Process effects and safety
        self._process_effects_safety(report_data)
        
        # Update popularity score
        self._update_popularity_score()

    def _update_dosage_data(self, roa: str, dose: float) -> None:
        """Update dosage data with new report."""
        if roa not in self.dosage_data:
            # Initialize with conservative range
            self.dosage_data[roa] = (dose, dose, dose)
        else:
            min_dose, max_dose, rec_dose = self.dosage_data[roa]
            # Update range
            min_dose = min(min_dose, dose)
            max_dose = max(max_dose, dose)
            # Update recommended dose (moving average)
            rec_dose = (rec_dose + dose) / 2
            self.dosage_data[roa] = (min_dose, max_dose, rec_dose)

    def _update_duration_data(
        self,
        phase: str,
        onset: float,
        duration: float
    ) -> None:
        """Update duration data with new report."""
        if phase not in self.duration_data:
            self.duration_data[phase] = (onset, duration)
        else:
            old_onset, old_duration = self.duration_data[phase]
            # Update with moving average
            new_onset = (old_onset + onset) / 2
            new_duration = (old_duration + duration) / 2
            self.duration_data[phase] = (new_onset, new_duration)

    def _update_effect_data(
        self,
        effect: str,
        magnitude: float,
        confidence: float
    ) -> None:
        """Update effect data with new report."""
        if effect not in self.reported_effects:
            self.reported_effects[effect] = []
            self.effect_frequencies[effect] = 0
            
        # Add new score
        self.reported_effects[effect].append((magnitude, confidence))
        
        # Update frequency
        self.effect_frequencies[effect] = (
            len(self.reported_effects[effect]) / self.total_reports
        )

    def _update_risk_data(self, risk: str, level: str) -> None:
        """Update risk data with new report."""
        try:
            risk_level = RiskLevel[level.upper()]
            if risk not in self.reported_risks:
                self.reported_risks[risk] = risk_level
            else:
                # Keep highest risk level
                self.reported_risks[risk] = max(
                    self.reported_risks[risk],
                    risk_level,
                    key=lambda x: x.value
                )
        except KeyError:
            pass

    def _update_interaction_data(
        self,
        substance: str,
        interaction: Dict
    ) -> None:
        """Update interaction data with new report."""
        if substance not in self.interactions:
            self.interactions[substance] = interaction
        else:
            # Update risk level if higher
            old_level = self.interactions[substance].get("risk_level")
            new_level = interaction.get("risk_level")
            if old_level and new_level:
                try:
                    old_risk = RiskLevel[old_level.upper()]
                    new_risk = RiskLevel[new_level.upper()]
                    if new_risk.value > old_risk.value:
                        self.interactions[substance] = interaction
                except KeyError:
                    pass

    def _update_popularity_score(self) -> None:
        """Update popularity score based on report metrics."""
        # Factors that influence popularity:
        # 1. Number of reports
        # 2. Recency of reports
        # 3. Diversity of sources
        
        # Base score from number of reports (logarithmic scale)
        import math
        report_score = math.log(self.total_reports + 1) / 10  # 0.0 - 1.0
        
        # Recency score
        recency_score = 0.0
        if self.last_report_date:
            days_since = (datetime.now() - self.last_report_date).days
            recency_score = 1.0 / (1.0 + days_since/365)  # Decay over a year
            
        # Source diversity score
        source_score = len(self.report_sources) / 10  # Assume max 10 sources
        
        # Combine scores with weights
        self.popularity_score = (
            0.4 * report_score +
            0.4 * recency_score +
            0.2 * source_score
        )

    def _merge_report_data(self, other: "CommunityDataMixin") -> None:
        """Merge basic report data."""
        self.experience_reports.update(other.experience_reports)
        self.report_sources.update(other.report_sources)
        self.total_reports += other.total_reports
        
        if (
            other.last_report_date and
            (not self.last_report_date or other.last_report_date > self.last_report_date)
        ):
            self.last_report_date = other.last_report_date

    def _merge_dosage_duration(self, other: "CommunityDataMixin") -> None:
        """Merge dosage and duration data."""
        # Merge dosage data
        for roa, dose_range in other.dosage_data.items():
            if roa not in self.dosage_data:
                self.dosage_data[roa] = dose_range
            else:
                min_dose, max_dose, rec_dose = self.dosage_data[roa]
                other_min, other_max, other_rec = dose_range
                self.dosage_data[roa] = (
                    min(min_dose, other_min),
                    max(max_dose, other_max),
                    (rec_dose + other_rec) / 2,
                )
                
        # Merge duration data
        for phase, time_range in other.duration_data.items():
            if phase not in self.duration_data:
                self.duration_data[phase] = time_range
            else:
                onset, duration = self.duration_data[phase]
                other_onset, other_duration = time_range
                self.duration_data[phase] = (
                    (onset + other_onset) / 2,
                    (duration + other_duration) / 2,
                )

    def _merge_effects_risks(self, other: "CommunityDataMixin") -> None:
        """Merge effects and risks data."""
        # Merge effect data
        for effect, scores in other.reported_effects.items():
            if effect not in self.reported_effects:
                self.reported_effects[effect] = scores
            else:
                self.reported_effects[effect].extend(scores)
                
        # Update effect frequencies
        for effect in self.reported_effects:
            self.effect_frequencies[effect] = (
                len(self.reported_effects[effect]) / self.total_reports
            )
            
        # Merge risk data
        for risk, level in other.reported_risks.items():
            if risk not in self.reported_risks:
                self.reported_risks[risk] = level
            else:
                self.reported_risks[risk] = max(
                    self.reported_risks[risk],
                    level,
                    key=lambda x: x.value
                )
                
        # Merge interaction data
        for substance, interaction in other.interactions.items():
            self._update_interaction_data(substance, interaction)

    def merge_community_data(self, other: "CommunityDataMixin") -> None:
        """Merge community data from another instance."""
        # Merge report data
        self._merge_report_data(other)
        
        # Merge dosage and duration data
        self._merge_dosage_duration(other)
        
        # Merge effects and risks data
        self._merge_effects_risks(other)
        
        # Update popularity score
        self._update_popularity_score()
