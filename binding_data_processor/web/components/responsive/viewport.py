"""Enhanced viewport management system optimized for scientific data visualization."""

from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, Dict, List, Optional, Set, Union

from .base import Breakpoint


class ColorSpace(Enum):
    """Color space capabilities for scientific visualization."""

    SRGB = "srgb"
    P3 = "p3"  # Wide color gamut
    REC2020 = "rec2020"  # Ultra-wide color gamut


class PointerType(Enum):
    """Input device pointer characteristics."""

    COARSE = "coarse"  # Touch input
    FINE = "fine"  # Mouse/stylus
    NONE = "none"  # No pointing device


@dataclass
class ViewportCapabilities:
    """Device capability configuration."""

    color_space: ColorSpace = ColorSpace.SRGB
    pixel_density: float = 1.0
    pointer_type: PointerType = PointerType.FINE
    supports_hover: bool = True
    supports_touch: bool = False
    supports_3d: bool = False
    max_texture_size: int = 2048
    memory_class: str = "default"  # low, default, high
    performance_class: str = "default"  # low, default, high


@dataclass
class ViewportMetrics:
    """Enhanced viewport metrics with scientific optimization."""

    width: int
    height: int
    device_pixel_ratio: float
    orientation: str
    color_scheme: str
    reduced_motion: bool
    safe_area_insets: Dict[str, int] = field(default_factory=dict)
    capabilities: ViewportCapabilities = field(default_factory=ViewportCapabilities)
    performance_metrics: Dict[str, float] = field(default_factory=dict)


@dataclass
class ViewportPreferences:
    """User viewport preferences with scientific display options."""

    zoom_level: float = 1.0
    font_size: int = 16
    color_scheme: str = "light"
    reduced_motion: bool = False
    high_contrast: bool = False
    safe_area_inset: bool = True
    data_precision: int = 4  # Decimal places for numerical display
    scientific_notation: bool = False
    use_si_units: bool = True
    enable_animations: bool = True
    enable_transitions: bool = True


class EnhancedViewportManager:
    """Sophisticated viewport management system with scientific optimization."""

    def __init__(self):
        """Initialize viewport manager with advanced features."""
        self.metrics = ViewportMetrics(
            width=1024,
            height=768,
            device_pixel_ratio=1.0,
            orientation="landscape",
            color_scheme="light",
            reduced_motion=False,
            safe_area_insets={"top": 0, "right": 0, "bottom": 0, "left": 0},
            performance_metrics={"fps": 60.0, "frame_time": 16.67, "cpu_load": 0.0, "gpu_load": 0.0, "memory_usage": 0.0},
        )
        self.preferences = ViewportPreferences()
        self.breakpoints: Dict[str, Breakpoint] = {
            "xs": Breakpoint("xs", 0, 575),
            "sm": Breakpoint("sm", 576, 767),
            "md": Breakpoint("md", 768, 991),
            "lg": Breakpoint("lg", 992, 1199),
            "xl": Breakpoint("xl", 1200),
        }
        self._current_breakpoint: str = "md"
        self._media_features: Dict[str, str] = {}
        self._viewport_listeners: List[Callable[[ViewportMetrics], None]] = []
        self._performance_monitors: Set[Callable[[], Dict[str, float]]] = set()
        self._active_optimizations: Set[str] = set()

    async def setup(self) -> None:
        """Set up viewport with comprehensive initialization."""
        await self._detect_device_capabilities()
        await self._configure_viewport()
        await self._setup_media_features()
        await self._setup_scientific_display()
        await self._initialize_performance_monitoring()
        await self._apply_optimizations()

    async def _detect_device_capabilities(self) -> None:
        """Detect and configure device capabilities."""
        # This would integrate with actual browser detection
        capabilities = ViewportCapabilities(
            color_space=self._detect_color_space(),
            pixel_density=self.metrics.device_pixel_ratio,
            pointer_type=self._detect_pointer_type(),
            supports_hover=not self._is_touch_device(),
            supports_touch=self._is_touch_device(),
            supports_3d=self._detect_3d_support(),
            max_texture_size=self._get_max_texture_size(),
            memory_class=self._detect_memory_class(),
            performance_class=self._detect_performance_class(),
        )
        self.metrics.capabilities = capabilities

    def _detect_color_space(self) -> ColorSpace:
        """Detect supported color space capabilities."""
        # This would integrate with actual browser detection
        return ColorSpace.P3 if self._supports_p3() else ColorSpace.SRGB

    def _detect_pointer_type(self) -> PointerType:
        """Detect primary pointing device type."""
        if self._is_touch_device():
            return PointerType.COARSE
        return PointerType.FINE

    def _is_touch_device(self) -> bool:
        """Detect touch capability."""
        return self.metrics.width <= self.breakpoints["md"].max_width

    async def _configure_viewport(self) -> None:
        """Configure viewport with comprehensive settings."""
        self._media_features.update(
            {
                "width": "device-width",
                "initial-scale": str(self.preferences.zoom_level),
                "viewport-fit": "cover" if self.preferences.safe_area_inset else "contain",
                "user-scalable": "yes",
                "minimum-scale": "1.0",
                "maximum-scale": "5.0",
            }
        )

        if self.metrics.device_pixel_ratio > 1:
            self._media_features["resolution"] = f"{self.metrics.device_pixel_ratio}dppx"

        # Scientific display optimizations
        if self.metrics.capabilities.color_space != ColorSpace.SRGB:
            self._media_features["color-gamut"] = self.metrics.capabilities.color_space.value

        # Accessibility features
        if self.preferences.high_contrast:
            self._media_features["prefers-contrast"] = "high"
        if self.preferences.reduced_motion:
            self._media_features["prefers-reduced-motion"] = "reduce"

    async def _setup_scientific_display(self) -> None:
        """Configure advanced scientific display features."""
        if self.metrics.capabilities.color_space == ColorSpace.P3:
            self._active_optimizations.add("wide_gamut_color")

        if self.metrics.device_pixel_ratio > 1:
            self._active_optimizations.add("high_dpi_rendering")

        if self.metrics.capabilities.supports_3d:
            self._active_optimizations.add("hardware_acceleration")

    async def _initialize_performance_monitoring(self) -> None:
        """Initialize performance monitoring systems."""
        self._performance_monitors.add(self._monitor_frame_rate)
        self._performance_monitors.add(self._monitor_memory_usage)
        self._performance_monitors.add(self._monitor_gpu_utilization)

    async def update_metrics(self, width: int, height: int) -> None:
        """Update viewport metrics with performance optimization."""
        new_breakpoint = self._get_breakpoint(width)
        metrics_changed = width != self.metrics.width or height != self.metrics.height or new_breakpoint != self._current_breakpoint

        if metrics_changed:
            self.metrics.width = width
            self.metrics.height = height
            self._current_breakpoint = new_breakpoint
            await self._reconfigure_for_breakpoint()
            await self._notify_listeners()
            await self._update_performance_metrics()

    def _get_breakpoint(self, width: int) -> str:
        """Get appropriate breakpoint for viewport width."""
        for name, bp in self.breakpoints.items():
            if bp.min_width <= width and (bp.max_width is None or width <= bp.max_width):
                return name
        return "xs"

    async def _reconfigure_for_breakpoint(self) -> None:
        """Reconfigure viewport for new breakpoint."""
        await self._configure_viewport()
        await self._optimize_for_breakpoint()
        await self._update_scientific_display()

    async def _optimize_for_breakpoint(self) -> None:
        """Apply breakpoint-specific optimizations."""
        if self._current_breakpoint in ["xs", "sm"]:
            self._active_optimizations.add("mobile_optimization")
            self._active_optimizations.add("touch_optimization")
        else:
            self._active_optimizations.discard("mobile_optimization")
            if not self._is_touch_device():
                self._active_optimizations.discard("touch_optimization")

    def get_scientific_display_config(self) -> Dict[str, Union[str, bool, int]]:
        """Get scientific display configuration."""
        return {
            "precision": self.preferences.data_precision,
            "notation": "scientific" if self.preferences.scientific_notation else "standard",
            "units": "si" if self.preferences.use_si_units else "imperial",
            "color_space": self.metrics.capabilities.color_space.value,
            "high_contrast": self.preferences.high_contrast,
            "animations_enabled": self.preferences.enable_animations,
            "transitions_enabled": self.preferences.enable_transitions,
        }

    def get_optimization_status(self) -> Dict[str, bool]:
        """Get status of active optimizations."""
        return {opt: opt in self._active_optimizations for opt in ["wide_gamut_color", "high_dpi_rendering", "hardware_acceleration", "mobile_optimization", "touch_optimization"]}

    async def _update_performance_metrics(self) -> None:
        """Update performance metrics from all monitors."""
        for monitor in self._performance_monitors:
            self.metrics.performance_metrics.update(await monitor())

    async def _monitor_frame_rate(self) -> Dict[str, float]:
        """Monitor frame rate and timing metrics."""
        # This would integrate with actual performance monitoring
        return {"fps": 60.0, "frame_time": 16.67}

    async def _monitor_memory_usage(self) -> Dict[str, float]:
        """Monitor memory usage metrics."""
        return {"memory_usage": 0.0}

    async def _monitor_gpu_utilization(self) -> Dict[str, float]:
        """Monitor GPU utilization metrics."""
        return {"gpu_load": 0.0}

    def get_meta_tag(self) -> str:
        """Get complete viewport meta tag."""
        features = [f"{key}={value}" for key, value in self._media_features.items()]
        return f'<meta name="viewport" content="{", ".join(features)}">'
