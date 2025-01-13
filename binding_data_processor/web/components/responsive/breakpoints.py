"""Breakpoint manager implementation."""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from .base import Breakpoint


@dataclass
class BreakpointConfig:
    """Breakpoint configuration settings."""

    name: str
    min_width: int
    max_width: Optional[int]
    columns: int
    gap: int
    margin: int
    container_width: Optional[int]


class BreakpointManager:
    """Manages responsive breakpoints and related functionality."""

    def __init__(self):
        self.breakpoints: Dict[str, BreakpointConfig] = {
            "xs": BreakpointConfig(
                name="xs",
                min_width=0,
                max_width=575,
                columns=4,
                gap=8,
                margin=8,
                container_width=None,
            ),
            "sm": BreakpointConfig(
                name="sm",
                min_width=576,
                max_width=767,
                columns=8,
                gap=12,
                margin=12,
                container_width=540,
            ),
            "md": BreakpointConfig(
                name="md",
                min_width=768,
                max_width=991,
                columns=12,
                gap=16,
                margin=16,
                container_width=720,
            ),
            "lg": BreakpointConfig(
                name="lg",
                min_width=992,
                max_width=1199,
                columns=12,
                gap=24,
                margin=24,
                container_width=960,
            ),
            "xl": BreakpointConfig(
                name="xl",
                min_width=1200,
                max_width=None,
                columns=12,
                gap=32,
                margin=32,
                container_width=1140,
            ),
        }
        self.ordered_breakpoints = ["xs", "sm", "md", "lg", "xl"]

    def get_breakpoint(self, width: int) -> str:
        """Get breakpoint name for viewport width."""
        for name in self.ordered_breakpoints:
            config = self.breakpoints[name]
            if config.min_width <= width and (config.max_width is None or width <= config.max_width):
                return name
        return "xs"  # Default to smallest breakpoint

    def get_config(self, breakpoint: str) -> BreakpointConfig:
        """Get configuration for breakpoint."""
        return self.breakpoints[breakpoint]

    def get_container_width(self, breakpoint: str) -> int:
        """Get container width for breakpoint."""
        config = self.breakpoints[breakpoint]
        return config.container_width or 100  # 100% for xs

    def get_grid_columns(self, breakpoint: str) -> int:
        """Get number of grid columns for breakpoint."""
        return self.breakpoints[breakpoint].columns

    def get_spacing(self, breakpoint: str) -> int:
        """Get spacing unit for breakpoint."""
        return self.breakpoints[breakpoint].gap

    def get_margin(self, breakpoint: str) -> int:
        """Get margin for breakpoint."""
        return self.breakpoints[breakpoint].margin

    def is_mobile(self, breakpoint: str) -> bool:
        """Check if breakpoint is considered mobile."""
        return breakpoint in ["xs", "sm"]

    def is_tablet(self, breakpoint: str) -> bool:
        """Check if breakpoint is considered tablet."""
        return breakpoint == "md"

    def is_desktop(self, breakpoint: str) -> bool:
        """Check if breakpoint is considered desktop."""
        return breakpoint in ["lg", "xl"]

    def get_media_query(self, breakpoint: str, type_: str = "min") -> str:
        """Get media query for breakpoint."""
        config = self.breakpoints[breakpoint]
        if type_ == "min":
            return f"@media (min-width: {config.min_width}px)"
        else:  # max
            if config.max_width is None:
                return f"@media (min-width: {config.min_width}px)"
            return f"@media (max-width: {config.max_width}px)"

    def get_breakpoint_range(self, start: str, end: str) -> str:
        """Get media query for breakpoint range."""
        start_config = self.breakpoints[start]
        end_config = self.breakpoints[end]
        return f"@media (min-width: {start_config.min_width}px) and " f"(max-width: {end_config.max_width or 9999}px)"

    def get_next_breakpoint(self, breakpoint: str) -> Optional[str]:
        """Get next larger breakpoint."""
        try:
            idx = self.ordered_breakpoints.index(breakpoint)
            if idx < len(self.ordered_breakpoints) - 1:
                return self.ordered_breakpoints[idx + 1]
        except ValueError:
            pass
        return None

    def get_prev_breakpoint(self, breakpoint: str) -> Optional[str]:
        """Get next smaller breakpoint."""
        try:
            idx = self.ordered_breakpoints.index(breakpoint)
            if idx > 0:
                return self.ordered_breakpoints[idx - 1]
        except ValueError:
            pass
        return None

    def get_breakpoints_up(self, breakpoint: str) -> List[str]:
        """Get list of breakpoints from given one up."""
        try:
            idx = self.ordered_breakpoints.index(breakpoint)
            return self.ordered_breakpoints[idx:]
        except ValueError:
            return []

    def get_breakpoints_down(self, breakpoint: str) -> List[str]:
        """Get list of breakpoints from given one down."""
        try:
            idx = self.ordered_breakpoints.index(breakpoint)
            return self.ordered_breakpoints[: idx + 1]
        except ValueError:
            return []

    def get_breakpoint_bounds(self, breakpoint: str) -> Tuple[int, Optional[int]]:
        """Get min and max width for breakpoint."""
        config = self.breakpoints[breakpoint]
        return config.min_width, config.max_width
