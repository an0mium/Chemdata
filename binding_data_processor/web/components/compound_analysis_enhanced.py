"""Enhanced compound analysis component.

This module provides an enhanced web interface for analyzing compound data with:
- Binding analysis
- Activity analysis
- Safety analysis
- SAR analysis
- Property analysis
- Community data analysis
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from datetime import datetime

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Fragments
from flask import render_template, request, jsonify

from ...models.compound import Compound
from ..base import BaseComponent, ViewResult


class CompoundAnalysisEnhanced(BaseComponent):
    """Enhanced compound analysis component."""

    def __init__(
        self,
        template_dir: Optional[Path] = None,
        static_dir: Optional[Path] = None,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize analysis component.
        
        Args:
            template_dir: Optional template directory
            static_dir: Optional static files directory
            logger: Optional logger instance
        """
        super().__init__(template_dir, static_dir, logger)

        # Initialize tracking
        self.analysis_stats = {
            "total_analyses": 0,
            "binding_analyses": 0,
            "activity_analyses": 0,
            "safety_analyses": 0,
            "sar_analyses": 0,
            "property_analyses": 0,
            "community_analyses": 0,
            "analysis_history": [],
        }

    def render_analysis(
        self,
        compound: Compound,
        analysis_types: Optional[List[str]] = None,
    ) -> ViewResult:
        """Render analysis interface.
        
        Args:
            compound: Compound to analyze
            analysis_types: Optional list of analysis types to perform
            
        Returns:
            ViewResult containing rendered HTML and metadata
        """
        try:
            # Generate analyses
            analyses = {}

            # Analysis types
            if not analysis_types:
                analysis_types = [
                    "binding",
                    "activity",
                    "safety",
                    "sar",
                    "properties",
                    "community",
                ]

            # Perform requested analyses
            for analysis_type in analysis_types:
                if analysis_type == "binding" and hasattr(compound, "binding_data"):
                    analysis = self._analyze_binding(compound)
                    analyses["binding_analysis"] = analysis
                    self.analysis_stats["binding_analyses"] += 1

                elif analysis_type == "activity" and hasattr(compound, "predictions"):
                    analysis = self._analyze_activity(compound)
                    analyses["activity_analysis"] = analysis
                    self.analysis_stats["activity_analyses"] += 1

                elif analysis_type == "safety" and hasattr(compound, "predictions"):
                    analysis = self._analyze_safety(compound)
                    analyses["safety_analysis"] = analysis
                    self.analysis_stats["safety_analyses"] += 1

                elif analysis_type == "sar":
                    analysis = self._analyze_sar(compound)
                    analyses["sar_analysis"] = analysis
                    self.analysis_stats["sar_analyses"] += 1

                elif analysis_type == "properties":
                    analysis = self._analyze_properties(compound)
                    analyses["property_analysis"] = analysis
                    self.analysis_stats["property_analyses"] += 1

                elif analysis_type == "community" and hasattr(compound, "social_data"):
                    analysis = self._analyze_community(compound)
                    analyses["community_analysis"] = analysis
                    self.analysis_stats["community_analyses"] += 1

            # Update stats
            self.analysis_stats["total_analyses"] += 1
            self.analysis_stats["analysis_history"].append({
                "timestamp": datetime.now().isoformat(),
                "compound": compound.name,
                "analysis_types": analysis_types,
            })

            # Render template
            html = render_template(
                "compound_analysis.html",
                compound=compound,
                analyses=analyses,
                analysis_types=analysis_types,
                stats=self.analysis_stats,
            )

            return ViewResult(
                success=True,
                data={
                    "html": html,
                    "analyses": analyses,
                    "analysis_types": analysis_types,
                    "analysis_stats": self.analysis_stats,
                },
            )

        except Exception as e:
            self.logger.error(f"Error rendering analysis: {str(e)}")
            return ViewResult(
                success=False,
                error=str(e),
            )

    def _analyze_binding(
        self,
        compound: Compound,
    ) -> Dict[str, Any]:
        """Analyze binding data.
        
        Args:
            compound: Compound to analyze
            
        Returns:
            Analysis results
        """
        # Extract data
        targets = []
        affinities = []
        confidences = []

        for binding in compound.binding_data:
            targets.append(binding["target"])
            affinities.append(float(binding["affinity"]))
            confidences.append(binding.get("confidence", 1.0))

        # Calculate statistics
        stats = {
            "strongest_target": targets[np.argmin(affinities)],
            "strongest_affinity": min(affinities),
            "mean_affinity": np.mean(affinities),
            "mean_confidence": np.mean(confidences),
            "target_count": len(targets),
        }

        # Analyze selectivity
        selectivity = {}
        for i, target in enumerate(targets):
            for j, other_target in enumerate(targets):
                if i != j:
                    ratio = affinities[j] / affinities[i]
                    selectivity[f"{target} vs {other_target}"] = ratio

        return {
            "statistics": stats,
            "selectivity": selectivity,
            "targets": targets,
            "affinities": affinities,
            "confidences": confidences,
        }

    def _analyze_activity(
        self,
        compound: Compound,
    ) -> Dict[str, Any]:
        """Analyze activity predictions.
        
        Args:
            compound: Compound to analyze
            
        Returns:
            Analysis results
        """
        # Extract data
        activities = compound.predictions.get("activity", {})

        # Calculate statistics
        stats = {
            "primary_activity": max(activities.items(), key=lambda x: x[1])[0],
            "primary_score": max(activities.values()),
            "mean_score": np.mean(list(activities.values())),
            "activity_count": len(activities),
        }

        # Analyze patterns
        patterns = {}
        for activity, score in activities.items():
            if score > 0.7:
                patterns["high"] = patterns.get("high", []) + [activity]
            elif score > 0.3:
                patterns["moderate"] = patterns.get("moderate", []) + [activity]
            else:
                patterns["low"] = patterns.get("low", []) + [activity]

        return {
            "statistics": stats,
            "patterns": patterns,
            "activities": activities,
        }

    def _analyze_safety(
        self,
        compound: Compound,
    ) -> Dict[str, Any]:
        """Analyze safety predictions.
        
        Args:
            compound: Compound to analyze
            
        Returns:
            Analysis results
        """
        # Extract data
        safety = compound.predictions.get("safety", {})

        # Calculate statistics
        stats = {
            "highest_risk": max(safety.items(), key=lambda x: x[1])[0],
            "highest_score": max(safety.values()),
            "mean_score": np.mean(list(safety.values())),
            "risk_count": len(safety),
        }

        # Analyze patterns
        patterns = {}
        for risk, score in safety.items():
            if score > 0.7:
                patterns["high"] = patterns.get("high", []) + [risk]
            elif score > 0.3:
                patterns["moderate"] = patterns.get("moderate", []) + [risk]
            else:
                patterns["low"] = patterns.get("low", []) + [risk]

        return {
            "statistics": stats,
            "patterns": patterns,
            "risks": safety,
        }

    def _analyze_sar(
        self,
        compound: Compound,
    ) -> Dict[str, Any]:
        """Analyze structure-activity relationships.
        
        Args:
            compound: Compound to analyze
            
        Returns:
            Analysis results
        """
        # Parse structure
        mol = Chem.MolFromSmiles(compound.smiles)
        if not mol:
            raise ValueError(f"Invalid SMILES: {compound.smiles}")

        # Calculate fingerprints
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)

        # Calculate descriptors
        descriptors = {
            "MW": Descriptors.ExactMolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "TPSA": Descriptors.TPSA(mol),
            "HBA": Descriptors.NumHAcceptors(mol),
            "HBD": Descriptors.NumHDonors(mol),
            "RotBonds": Descriptors.NumRotatableBonds(mol),
            "Rings": Descriptors.RingCount(mol),
            "AromaticRings": Descriptors.NumAromaticRings(mol),
        }

        # Detect fragments
        fragments = {
            "Aromatic": Fragments.fr_Ar_N(mol),
            "Alkyl": Fragments.fr_alkyl_halide(mol),
            "Carbonyl": Fragments.fr_C_O(mol),
            "Carboxyl": Fragments.fr_COO(mol),
            "Amine": Fragments.fr_NH2(mol),
            "Amide": Fragments.fr_amide(mol),
        }

        return {
            "fingerprint": list(fp.GetOnBits()),
            "descriptors": descriptors,
            "fragments": fragments,
        }

    def _analyze_properties(
        self,
        compound: Compound,
    ) -> Dict[str, Any]:
        """Analyze physicochemical properties.
        
        Args:
            compound: Compound to analyze
            
        Returns:
            Analysis results
        """
        # Parse structure
        mol = Chem.MolFromSmiles(compound.smiles)
        if not mol:
            raise ValueError(f"Invalid SMILES: {compound.smiles}")

        # Calculate properties
        properties = {
            "MW": Descriptors.ExactMolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "TPSA": Descriptors.TPSA(mol),
            "HBA": Descriptors.NumHAcceptors(mol),
            "HBD": Descriptors.NumHDonors(mol),
            "RotBonds": Descriptors.NumRotatableBonds(mol),
            "Rings": Descriptors.RingCount(mol),
            "AromaticRings": Descriptors.NumAromaticRings(mol),
            "HeavyAtoms": mol.GetNumHeavyAtoms(),
            "Charge": Chem.GetFormalCharge(mol),
        }

        # Calculate drug-likeness
        lipinski = {
            "MW_ok": properties["MW"] <= 500,
            "LogP_ok": properties["LogP"] <= 5,
            "HBA_ok": properties["HBA"] <= 10,
            "HBD_ok": properties["HBD"] <= 5,
            "RotBonds_ok": properties["RotBonds"] <= 10,
        }

        return {
            "properties": properties,
            "lipinski": lipinski,
            "violations": 5 - sum(lipinski.values()),
        }

    def _analyze_community(
        self,
        compound: Compound,
    ) -> Dict[str, Any]:
        """Analyze community data.
        
        Args:
            compound: Compound to analyze
            
        Returns:
            Analysis results
        """
        # Extract data
        social_data = compound.social_data
        reddit_posts = social_data.get("reddit", {}).get("posts", [])
        twitter_tweets = social_data.get("twitter", {}).get("tweets", [])

        # Calculate statistics
        stats = {
            "reddit_posts": len(reddit_posts),
            "twitter_mentions": len(twitter_tweets),
            "total_mentions": len(reddit_posts) + len(twitter_tweets),
        }

        # Analyze dates
        dates = []
        for post in reddit_posts:
            dates.append(datetime.fromisoformat(post["created_utc"]))
        for tweet in twitter_tweets:
            dates.append(datetime.fromisoformat(tweet["created_at"]))

        if dates:
            date_stats = {
                "first_mention": min(dates).isoformat(),
                "last_mention": max(dates).isoformat(),
                "mention_days": (max(dates) - min(dates)).days,
            }
        else:
            date_stats = {}

        return {
            "statistics": stats,
            "date_stats": date_stats,
            "reddit_posts": reddit_posts,
            "twitter_tweets": twitter_tweets,
        }

    def get_metrics(self) -> Dict[str, Any]:
        """Get component metrics."""
        return {
            "analysis_stats": self.analysis_stats,
        }
