"""Storage system for Bluelight monitoring data.

This module provides persistent storage for:
1. Monitored Bluelight posts
2. Extracted compounds
3. Safety alerts
4. Analysis results
5. Screenshots and JavaScript logs
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict

from ...processors.psychopharm.predictors.bbb.predictors import BBBPredictor
from ..validation.schema import BaseSchema


@dataclass
class BluelightAlert(BaseSchema):
    """Alert data for monitored posts."""

    post_id: str
    subforum: str
    title: str
    url: str
    compounds: List[str]
    safety_notes: List[str]
    severity: str  # high, medium, low
    created_at: str
    bbb_predictions: Dict[str, float]
    metadata: Dict[str, Any]
    screenshot: Optional[str] = None
    javascript_logs: Optional[List[str]] = None


class BluelightStorage:
    """Storage manager for Bluelight monitoring data."""

    def __init__(self, storage_dir: str = "data/bluelight"):
        """Initialize storage.

        Args:
            storage_dir: Directory for storing data
        """
        self.storage_dir = Path(storage_dir)
        self.storage_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(__name__)

        # Initialize BBB predictor
        self.bbb_predictor = BBBPredictor()

        # Load existing data
        self.posts_file = self.storage_dir / "posts.json"
        self.alerts_file = self.storage_dir / "alerts.json"
        self.compounds_file = self.storage_dir / "compounds.json"
        self.errors_file = self.storage_dir / "errors.json"

        self.posts = self._load_json(self.posts_file, default={})
        self.alerts = self._load_json(self.alerts_file, default=[])
        self.compounds = self._load_json(self.compounds_file, default={})
        self.errors = self._load_json(self.errors_file, default=[])

    def _load_json(self, path: Path, default: Any = None) -> Any:
        """Load JSON data from file.

        Args:
            path: File path
            default: Default value if file doesn't exist

        Returns:
            Loaded data or default
        """
        try:
            if path.exists():
                with open(path, "r") as f:
                    return json.load(f)
        except Exception as e:
            self.logger.error(f"Error loading {path}: {str(e)}")
        return default

    def _save_json(self, path: Path, data: Any):
        """Save data to JSON file.

        Args:
            path: File path
            data: Data to save
        """
        try:
            with open(path, "w") as f:
                json.dump(data, f, indent=2)
        except Exception as e:
            self.logger.error(f"Error saving {path}: {str(e)}")

    def store_post(self, post_id: str, data: Dict[str, Any]):
        """Store Bluelight post data.

        Args:
            post_id: Unique post identifier
            data: Post data to store
        """
        self.posts[post_id] = {
            **data,
            "stored_at": datetime.now().isoformat(),
        }
        self._save_json(self.posts_file, self.posts)

        # Extract and store compounds
        if "compounds" in data:
            for compound in data["compounds"]:
                if compound not in self.compounds:
                    self.compounds[compound] = {
                        "first_seen": datetime.now().isoformat(),
                        "posts": [],
                        "bbb_prediction": None,
                    }
                self.compounds[compound]["posts"].append(post_id)
            self._save_json(self.compounds_file, self.compounds)

    def store_error(self, error: str, metadata: Optional[Dict[str, Any]] = None):
        """Store error information.

        Args:
            error: Error message
            metadata: Optional error context
        """
        self.errors.append(
            {
                "error": error,
                "timestamp": datetime.now().isoformat(),
                "metadata": metadata or {},
            }
        )
        self._save_json(self.errors_file, self.errors)

    def create_alert(
        self,
        post_id: str,
        compounds: List[str],
        safety_notes: List[str],
        severity: str = "medium",
        metadata: Optional[Dict[str, Any]] = None,
    ) -> BluelightAlert:
        """Create safety alert for post.

        Args:
            post_id: Post identifier
            compounds: Detected compounds
            safety_notes: Safety concerns
            severity: Alert severity
            metadata: Optional additional data

        Returns:
            Created alert
        """
        if post_id not in self.posts:
            raise ValueError(f"Post {post_id} not found")

        post = self.posts[post_id]

        # Get BBB predictions for compounds
        bbb_predictions = {}
        for compound in compounds:
            try:
                if compound in self.compounds:
                    # Use cached prediction
                    pred = self.compounds[compound].get("bbb_prediction")
                    if pred is not None:
                        bbb_predictions[compound] = pred
                        continue

                # Get new prediction
                pred = self.bbb_predictor.predict(compound)
                bbb_predictions[compound] = pred

                # Cache prediction
                if compound in self.compounds:
                    self.compounds[compound]["bbb_prediction"] = pred
                    self._save_json(self.compounds_file, self.compounds)

            except Exception as e:
                self.logger.error(f"BBB prediction failed for {compound}: {str(e)}")

        alert = BluelightAlert(
            post_id=post_id,
            subforum=post["subforum"],
            title=post["title"],
            url=post.get("url", ""),
            compounds=compounds,
            safety_notes=safety_notes,
            severity=severity,
            created_at=datetime.now().isoformat(),
            bbb_predictions=bbb_predictions,
            metadata=metadata or {},
            screenshot=post.get("screenshot"),
            javascript_logs=post.get("javascript_logs"),
        )

        self.alerts.append(asdict(alert))
        self._save_json(self.alerts_file, self.alerts)
        return alert

    def get_alerts(
        self,
        severity: Optional[str] = None,
        subforum: Optional[str] = None,
        compound: Optional[str] = None,
        days: Optional[int] = None,
    ) -> List[BluelightAlert]:
        """Get filtered alerts.

        Args:
            severity: Filter by severity
            subforum: Filter by subforum
            compound: Filter by compound
            days: Filter by days (from now)

        Returns:
            List of matching alerts
        """
        filtered = []
        now = datetime.now()

        for alert_dict in self.alerts:
            alert = BluelightAlert(**alert_dict)

            if severity and alert.severity != severity:
                continue
            if subforum and alert.subforum != subforum:
                continue
            if compound and compound not in alert.compounds:
                continue
            if days:
                created = datetime.fromisoformat(alert.created_at)
                if (now - created).days > days:
                    continue

            filtered.append(alert)

        return filtered

    def get_compound_stats(self, compound: str) -> Dict[str, Any]:
        """Get statistics for compound.

        Args:
            compound: Compound name

        Returns:
            Compound statistics
        """
        if compound not in self.compounds:
            return {}

        data = self.compounds[compound]
        posts = [self.posts[pid] for pid in data["posts"] if pid in self.posts]
        alerts = [a for a in self.alerts if compound in a["compounds"]]

        return {
            "name": compound,
            "first_seen": data["first_seen"],
            "post_count": len(posts),
            "alert_count": len(alerts),
            "subforums": list({p["subforum"] for p in posts}),
            "bbb_prediction": data.get("bbb_prediction"),
            "latest_alerts": sorted(
                alerts,
                key=lambda x: x["created_at"],
                reverse=True,
            )[:5],
        }

    def get_trending_compounds(self, days: int = 7) -> List[Dict[str, Any]]:
        """Get trending compounds by post frequency.

        Args:
            days: Number of days to analyze

        Returns:
            List of trending compounds with stats
        """
        now = datetime.now()
        trending = []

        for compound, data in self.compounds.items():
            recent_posts = []
            for pid in data["posts"]:
                if pid in self.posts:
                    post = self.posts[pid]
                    created = datetime.fromisoformat(post["stored_at"])
                    if (now - created).days <= days:
                        recent_posts.append(post)

            if recent_posts:
                trending.append(
                    {
                        "name": compound,
                        "recent_posts": len(recent_posts),
                        "total_posts": len(data["posts"]),
                        "subforums": list({p["subforum"] for p in recent_posts}),
                        "bbb_prediction": data.get("bbb_prediction"),
                    }
                )

        return sorted(trending, key=lambda x: x["recent_posts"], reverse=True)

    def get_safety_summary(self, days: Optional[int] = None) -> Dict[str, Any]:
        """Get safety monitoring summary.

        Args:
            days: Optional number of days to analyze

        Returns:
            Safety summary statistics
        """
        alerts = self.get_alerts(days=days)
        compounds = set()
        subforums = set()
        severity_counts = {"high": 0, "medium": 0, "low": 0}
        bbb_concerns = []

        for alert in alerts:
            compounds.update(alert.compounds)
            subforums.add(alert.subforum)
            severity_counts[alert.severity] += 1

            # Check BBB predictions
            for compound, pred in alert.bbb_predictions.items():
                if pred > 0.7:  # High BBB permeability
                    bbb_concerns.append(
                        {
                            "compound": compound,
                            "prediction": pred,
                            "alert": alert,
                        }
                    )

        return {
            "total_alerts": len(alerts),
            "unique_compounds": len(compounds),
            "affected_subforums": list(subforums),
            "severity_distribution": severity_counts,
            "high_bbb_compounds": sorted(
                bbb_concerns,
                key=lambda x: x["prediction"],
                reverse=True,
            ),
        }

    def get_error_summary(self, days: Optional[int] = None) -> Dict[str, Any]:
        """Get error monitoring summary.

        Args:
            days: Optional number of days to analyze

        Returns:
            Error summary statistics
        """
        now = datetime.now()
        filtered_errors = []

        for error in self.errors:
            if days:
                created = datetime.fromisoformat(error["timestamp"])
                if (now - created).days > days:
                    continue
            filtered_errors.append(error)

        return {
            "total_errors": len(filtered_errors),
            "recent_errors": sorted(
                filtered_errors,
                key=lambda x: x["timestamp"],
                reverse=True,
            )[:10],
            "error_types": list({e["error"].split(":")[0] for e in filtered_errors}),
        }
