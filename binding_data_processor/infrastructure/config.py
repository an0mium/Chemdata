"""Configuration module for infrastructure and processing components.

This module provides a unified configuration system for:
1. Infrastructure components (cache, rate limiting, etc.)
2. Structure processing (pharmacophore, binding, etc.)
3. ML components (AlphaFold, GNN, etc.)
4. Data enrichment (web, social, literature)
5. Monitoring and observability
6. Psychopharm-specific components
7. Web interface components

The configuration classes are designed to:
- Enable seamless integration between components
- Provide consistent validation
- Support monitoring and observability
- Handle domain-specific requirements
"""

from dataclasses import dataclass, field
from datetime import timedelta
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Set, Tuple, Union


@dataclass
class CacheConfig:
    """Cache configuration."""

    # Storage settings
    backend: Literal["memory", "file", "redis"] = "memory"
    ttl: int = 300  # Default TTL in seconds
    max_size: int = 1000  # Max items for memory cache
    cache_dir: Optional[Path] = None  # Directory for file cache

    # Redis settings
    redis_host: str = "localhost"
    redis_port: int = 6379
    redis_db: int = 0

    # Performance settings
    compression: bool = False
    serializer: Literal["json", "pickle", "msgpack"] = "json"

    # Monitoring
    track_stats: bool = True
    stats_window: int = 3600  # 1 hour

    def __post_init__(self):
        """Validate configuration."""
        if self.ttl < 0:
            raise ValueError("TTL must be >= 0")
        if self.max_size < 1:
            raise ValueError("max_size must be >= 1")
        if self.backend == "file" and not self.cache_dir:
            raise ValueError("cache_dir required for file backend")
        if isinstance(self.cache_dir, str):
            self.cache_dir = Path(self.cache_dir)


@dataclass
class RateLimitConfig:
    """Rate limit configuration."""

    # Limit settings
    algorithm: Literal["token", "window"] = "token"
    limit: int = 100  # Requests per window
    window: int = 60  # Window size in seconds
    burst: Optional[int] = None  # Max burst size

    # Storage settings
    backend: Literal["memory", "redis"] = "memory"
    redis_host: str = "localhost"
    redis_port: int = 6379
    redis_db: int = 0

    # Monitoring
    track_stats: bool = True
    stats_window: int = 3600  # 1 hour

    def __post_init__(self):
        """Validate configuration."""
        if self.limit < 1:
            raise ValueError("limit must be >= 1")
        if self.window < 1:
            raise ValueError("window must be >= 1")
        if self.burst is not None and self.burst < self.limit:
            raise ValueError("burst must be >= limit")


@dataclass
class CircuitBreakerConfig:
    """Circuit breaker configuration."""

    # Failure settings
    failure_threshold: int = 5  # Failures before opening
    failure_timeout: int = 60  # Failure counting window
    min_throughput: int = 10  # Min requests before opening

    # Recovery settings
    reset_timeout: int = 60  # Time before half-open
    half_open_timeout: int = 30  # Time between test requests
    test_requests: int = 1  # Requests to try in half-open

    # Monitoring settings
    track_metrics: bool = True
    metric_window: int = 60  # Window size in seconds
    metric_buckets: int = 10  # Number of time buckets

    def __post_init__(self):
        """Validate configuration."""
        if self.failure_threshold < 1:
            raise ValueError("failure_threshold must be >= 1")
        if self.failure_timeout < 1:
            raise ValueError("failure_timeout must be >= 1")
        if self.min_throughput < 1:
            raise ValueError("min_throughput must be >= 1")


@dataclass
class StructureConfig:
    """Structure analysis configuration."""

    # Conformer generation
    max_conformers: int = 100
    energy_window: float = 20.0  # kcal/mol
    rmsd_threshold: float = 0.5  # Angstroms
    optimize_conformers: bool = True

    # Geometry optimization
    force_field: Literal["MMFF94", "UFF"] = "MMFF94"
    max_iterations: int = 1000
    energy_tolerance: float = 0.1
    gradient_tolerance: float = 0.1

    # Surface analysis
    surface_type: Literal["SAS", "VDW", "MS"] = "SAS"
    probe_radius: float = 1.4  # Angstroms
    grid_spacing: float = 0.5  # Angstroms

    # Pocket detection
    min_pocket_volume: float = 100.0  # Å³
    max_pockets: int = 10
    pocket_score_cutoff: float = 0.5

    def __post_init__(self):
        """Validate configuration."""
        if self.max_conformers < 1:
            raise ValueError("max_conformers must be >= 1")
        if self.energy_window <= 0:
            raise ValueError("energy_window must be > 0")
        if self.rmsd_threshold <= 0:
            raise ValueError("rmsd_threshold must be > 0")
        if self.max_iterations < 1:
            raise ValueError("max_iterations must be >= 1")
        if self.energy_tolerance <= 0:
            raise ValueError("energy_tolerance must be > 0")
        if self.gradient_tolerance <= 0:
            raise ValueError("gradient_tolerance must be > 0")
        if self.probe_radius <= 0:
            raise ValueError("probe_radius must be > 0")
        if self.grid_spacing <= 0:
            raise ValueError("grid_spacing must be > 0")
        if self.min_pocket_volume <= 0:
            raise ValueError("min_pocket_volume must be > 0")
        if self.max_pockets < 1:
            raise ValueError("max_pockets must be >= 1")
        if not 0 <= self.pocket_score_cutoff <= 1:
            raise ValueError("pocket_score_cutoff must be in [0, 1]")


@dataclass
class PharmacophoreConfig:
    """Pharmacophore detection and analysis configuration."""

    # Feature detection
    feature_types: Set[str] = field(default_factory=lambda: {"donor", "acceptor", "aromatic", "hydrophobic"})
    max_features: int = 10
    feature_radius: float = 1.5  # Angstroms

    # Alignment settings
    align_method: Literal["rigid", "flexible"] = "flexible"
    max_alignments: int = 10
    rmsd_cutoff: float = 2.0  # Angstroms

    # Scoring settings
    scoring_function: Literal["overlap", "pharmacophore", "hybrid"] = "hybrid"
    score_weights: Dict[str, float] = field(default_factory=lambda: {"overlap": 0.4, "pharmacophore": 0.6})

    # Performance
    cache_results: bool = True
    n_jobs: int = -1  # Number of parallel jobs (-1 for all cores)

    def __post_init__(self):
        """Validate configuration."""
        if not self.feature_types:
            raise ValueError("At least one feature type required")
        if self.max_features < 1:
            raise ValueError("max_features must be >= 1")
        if self.feature_radius <= 0:
            raise ValueError("feature_radius must be > 0")
        if self.max_alignments < 1:
            raise ValueError("max_alignments must be >= 1")
        if self.rmsd_cutoff <= 0:
            raise ValueError("rmsd_cutoff must be > 0")


@dataclass
class BindingConfig:
    """Binding site prediction and analysis configuration."""

    # Site detection
    detection_method: Literal["geometric", "energetic", "ml"] = "ml"
    min_site_size: int = 100  # Minimum site volume in Å³
    max_sites: int = 5  # Maximum number of sites to consider

    # Scoring settings
    scoring_functions: List[str] = field(default_factory=lambda: ["vina", "docking", "conservation"])
    score_weights: Dict[str, float] = field(default_factory=lambda: {"vina": 0.4, "docking": 0.4, "conservation": 0.2})

    # Interaction analysis
    interaction_types: Set[str] = field(default_factory=lambda: {"hbond", "hydrophobic", "pi_stacking", "salt_bridge"})
    distance_cutoffs: Dict[str, float] = field(default_factory=lambda: {"hbond": 3.5, "hydrophobic": 4.0, "pi_stacking": 5.0, "salt_bridge": 4.0})
    angle_cutoffs: Dict[str, float] = field(default_factory=lambda: {"hbond": 30.0, "pi_stacking": 30.0})

    # ML model settings
    model_path: Optional[Path] = None
    batch_size: int = 32
    use_gpu: bool = True

    # Performance
    cache_results: bool = True
    n_jobs: int = -1

    def __post_init__(self):
        """Validate configuration."""
        if self.min_site_size < 1:
            raise ValueError("min_site_size must be >= 1")
        if self.max_sites < 1:
            raise ValueError("max_sites must be >= 1")
        if not self.scoring_functions:
            raise ValueError("At least one scoring function required")
        if sum(self.score_weights.values()) != 1.0:
            raise ValueError("Score weights must sum to 1.0")
        if not self.interaction_types:
            raise ValueError("At least one interaction type required")
        if not all(v > 0 for v in self.distance_cutoffs.values()):
            raise ValueError("All distance cutoffs must be > 0")
        if not all(v > 0 for v in self.angle_cutoffs.values()):
            raise ValueError("All angle cutoffs must be > 0")


@dataclass
class AlphaFoldConfig:
    """AlphaFold integration configuration."""

    # Model settings
    model_preset: Literal["monomer", "monomer_casp14", "multimer"] = "monomer"
    model_path: Optional[Path] = None
    max_seq_len: int = 2500

    # Prediction settings
    num_ensemble: int = 1
    num_recycle: int = 3
    use_templates: bool = True

    # Performance
    use_gpu: bool = True
    cache_results: bool = True
    batch_size: int = 1

    # Output settings
    save_intermediates: bool = False
    output_dir: Optional[Path] = None

    def __post_init__(self):
        """Validate configuration."""
        if self.max_seq_len < 1:
            raise ValueError("max_seq_len must be >= 1")
        if self.num_ensemble < 1:
            raise ValueError("num_ensemble must be >= 1")
        if self.num_recycle < 0:
            raise ValueError("num_recycle must be >= 0")
        if self.batch_size < 1:
            raise ValueError("batch_size must be >= 1")


@dataclass
class MLConfig:
    """Machine learning configuration."""

    # Model architecture
    model_type: Literal["gnn", "cnn", "transformer"] = "gnn"
    hidden_dim: int = 256
    num_layers: int = 6
    dropout: float = 0.1
    activation: Literal["relu", "gelu", "swish"] = "gelu"

    # Training settings
    batch_size: int = 32
    learning_rate: float = 0.001
    weight_decay: float = 0.0001
    num_epochs: int = 100
    early_stopping: int = 10
    gradient_clip: float = 1.0

    # Feature generation
    atom_features: List[str] = field(default_factory=lambda: ["atomic_num", "degree", "formal_charge", "chiral_tag", "hybridization", "aromatic", "ring_size"])
    bond_features: List[str] = field(default_factory=lambda: ["bond_type", "conjugated", "rotatable", "ring"])
    graph_features: List[str] = field(default_factory=lambda: ["num_atoms", "num_bonds", "num_rings", "molecular_weight"])

    # Performance
    num_workers: int = 4
    prefetch_factor: int = 2
    pin_memory: bool = True
    use_amp: bool = True  # Automatic mixed precision
    compile_model: bool = True  # Using torch.compile

    def __post_init__(self):
        """Validate configuration."""
        if self.hidden_dim < 1:
            raise ValueError("hidden_dim must be >= 1")
        if self.num_layers < 1:
            raise ValueError("num_layers must be >= 1")
        if not 0 <= self.dropout < 1:
            raise ValueError("dropout must be in [0, 1)")
        if self.batch_size < 1:
            raise ValueError("batch_size must be >= 1")
        if self.learning_rate <= 0:
            raise ValueError("learning_rate must be > 0")
        if self.weight_decay < 0:
            raise ValueError("weight_decay must be >= 0")
        if self.num_epochs < 1:
            raise ValueError("num_epochs must be >= 1")
        if self.early_stopping < 1:
            raise ValueError("early_stopping must be >= 1")
        if self.gradient_clip <= 0:
            raise ValueError("gradient_clip must be > 0")
        if not self.atom_features:
            raise ValueError("At least one atom feature required")
        if not self.bond_features:
            raise ValueError("At least one bond feature required")
        if self.num_workers < 0:
            raise ValueError("num_workers must be >= 0")
        if self.prefetch_factor < 1:
            raise ValueError("prefetch_factor must be >= 1")


@dataclass
class EnrichmentConfig:
    """Data enrichment configuration."""

    # Web enrichment
    web_sources: Set[str] = field(default_factory=lambda: {"pubchem", "chembl", "bindingdb", "drugbank"})
    max_retries: int = 3
    timeout: int = 30
    concurrent_requests: int = 5

    # Social enrichment
    social_platforms: Set[str] = field(default_factory=lambda: {"reddit", "twitter", "bluelight"})
    max_posts: int = 1000
    max_comments: int = 5000
    min_relevance: float = 0.5

    # Literature enrichment
    literature_sources: Set[str] = field(default_factory=lambda: {"pubmed", "patents", "sciencedirect"})
    max_papers: int = 1000
    min_year: int = 1990
    require_full_text: bool = False

    def __post_init__(self):
        """Validate configuration."""
        if not self.web_sources:
            raise ValueError("At least one web source required")
        if self.max_retries < 1:
            raise ValueError("max_retries must be >= 1")
        if self.timeout < 1:
            raise ValueError("timeout must be >= 1")
        if self.concurrent_requests < 1:
            raise ValueError("concurrent_requests must be >= 1")
        if not self.social_platforms:
            raise ValueError("At least one social platform required")
        if self.max_posts < 1:
            raise ValueError("max_posts must be >= 1")
        if self.max_comments < 1:
            raise ValueError("max_comments must be >= 1")
        if not 0 <= self.min_relevance <= 1:
            raise ValueError("min_relevance must be in [0, 1]")
        if not self.literature_sources:
            raise ValueError("At least one literature source required")
        if self.max_papers < 1:
            raise ValueError("max_papers must be >= 1")


@dataclass
class MonitorConfig:
    """Monitoring and observability configuration."""

    # Core settings
    metrics_enabled: bool = True
    alerts_enabled: bool = True
    log_dir: Optional[Path] = None
    report_format: Literal["json", "prometheus"] = "json"

    # Metric settings
    metric_window: int = 60  # Window size in seconds
    metric_buckets: int = 10  # Number of time buckets
    default_tags: Dict[str, str] = field(default_factory=dict)

    # Alert settings
    alert_check_interval: int = 60  # Seconds between checks
    alert_history_size: int = 1000  # Max alerts to keep

    # Storage settings
    backend: Literal["memory", "file", "redis"] = "memory"
    storage_dir: Optional[Path] = None  # For file backend
    redis_host: str = "localhost"  # For redis backend
    redis_port: int = 6379
    redis_db: int = 0

    def __post_init__(self):
        """Validate configuration."""
        if self.metric_window < 1:
            raise ValueError("metric_window must be >= 1")
        if self.metric_buckets < 1:
            raise ValueError("metric_buckets must be >= 1")
        if self.alert_check_interval < 1:
            raise ValueError("alert_check_interval must be >= 1")
        if self.alert_history_size < 1:
            raise ValueError("alert_history_size must be >= 1")
        if self.backend == "file" and not self.storage_dir:
            raise ValueError("storage_dir required for file backend")


@dataclass
class DocumentConfig:
    """Document processing configuration."""

    # Input settings
    allowed_formats: Set[str] = field(default_factory=lambda: {"pdf", "docx", "txt", "html", "xml"})
    max_file_size: int = 10 * 1024 * 1024  # 10MB
    encoding: str = "utf-8"

    # Processing settings
    extract_tables: bool = True
    extract_images: bool = True
    ocr_enabled: bool = True
    language: str = "en"

    # Output settings
    output_format: Literal["json", "xml", "txt"] = "json"
    preserve_formatting: bool = True
    include_metadata: bool = True

    def __post_init__(self):
        """Validate configuration."""
        if not self.allowed_formats:
            raise ValueError("At least one format must be allowed")
        if self.max_file_size < 1:
            raise ValueError("max_file_size must be >= 1")


@dataclass
class StructureAnalysisConfig:
    """Structure analysis configuration integrating all structure-related components."""

    # Core structure settings
    structure: StructureConfig = field(default_factory=StructureConfig)
    pharmacophore: PharmacophoreConfig = field(default_factory=PharmacophoreConfig)
    binding: BindingConfig = field(default_factory=BindingConfig)

    # Integration settings
    shared_cache: bool = True  # Share cache between components
    coordinate_analysis: bool = True  # Coordinate analysis between components
    validate_results: bool = True  # Cross-validate results between components

    # Performance
    parallel_processing: bool = True
    max_parallel_jobs: int = -1  # -1 for all cores
    memory_limit: Optional[int] = None  # MB, None for no limit

    def __post_init__(self):
        """Validate configuration."""
        if self.max_parallel_jobs < -1:
            raise ValueError("max_parallel_jobs must be >= -1")
        if self.memory_limit is not None and self.memory_limit < 1:
            raise ValueError("memory_limit must be >= 1")

    def validate_integration(self):
        """Validate integration between structure components."""
        # Ensure consistent cache settings
        if self.shared_cache:
            self.pharmacophore.cache_results = True
            self.binding.cache_results = True

        # Align parallel processing settings
        if self.parallel_processing:
            self.pharmacophore.n_jobs = self.max_parallel_jobs
            self.binding.n_jobs = self.max_parallel_jobs


@dataclass
class MLSystemConfig:
    """Machine learning system configuration integrating all ML components."""

    # Core ML settings
    ml: MLConfig = field(default_factory=MLConfig)
    alphafold: AlphaFoldConfig = field(default_factory=AlphaFoldConfig)

    # Integration settings
    shared_gpu: bool = True  # Share GPU between components
    coordinate_training: bool = True  # Coordinate training between models
    validate_predictions: bool = True  # Cross-validate predictions

    # Resource management
    gpu_memory_fraction: float = 0.9  # Fraction of GPU memory to use
    cpu_threads: int = -1  # -1 for all threads
    batch_scheduling: bool = True  # Enable batch scheduling

    def __post_init__(self):
        """Validate configuration."""
        if not 0 < self.gpu_memory_fraction <= 1:
            raise ValueError("gpu_memory_fraction must be in (0, 1]")
        if self.cpu_threads < -1:
            raise ValueError("cpu_threads must be >= -1")

    def validate_integration(self):
        """Validate integration between ML components."""
        # Validate resource settings
        if not 0 < self.gpu_memory_fraction <= 1:
            raise ValueError("gpu_memory_fraction must be in (0, 1]")
        if self.cpu_threads < -1:
            raise ValueError("cpu_threads must be >= -1")

        # Ensure consistent GPU settings
        if self.shared_gpu:
            self.ml.use_gpu = True
            self.alphafold.use_gpu = True

        # Align batch processing settings
        if self.batch_scheduling:
            self.ml.batch_size = min(self.ml.batch_size, 32)
            self.alphafold.batch_size = min(self.alphafold.batch_size, 1)

        # Validate CPU thread allocation
        if self.cpu_threads > 0:
            self.ml.num_workers = min(self.ml.num_workers, self.cpu_threads)


@dataclass
class WebInterfaceConfig:
    """Web interface configuration."""

    # Server settings
    host: str = "localhost"
    port: int = 8000
    workers: int = 4
    debug: bool = False

    # API settings
    api_prefix: str = "/api/v1"
    rate_limit: int = 100  # requests per minute
    require_auth: bool = True
    token_expiry: int = 3600  # seconds

    # UI settings
    page_size: int = 20
    max_items: int = 1000
    enable_export: bool = True
    export_formats: Set[str] = field(default_factory=lambda: {"csv", "json", "sdf", "mol2"})

    # Visualization
    enable_3d: bool = True
    max_structure_size: int = 1000  # atoms
    image_width: int = 400
    image_height: int = 400

    def __post_init__(self):
        """Validate configuration."""
        if self.port < 1 or self.port > 65535:
            raise ValueError("port must be in range [1, 65535]")
        if self.workers < 1:
            raise ValueError("workers must be >= 1")
        if self.rate_limit < 1:
            raise ValueError("rate_limit must be >= 1")
        if self.token_expiry < 1:
            raise ValueError("token_expiry must be >= 1")
        if self.page_size < 1:
            raise ValueError("page_size must be >= 1")
        if self.max_items < 1:
            raise ValueError("max_items must be >= 1")
        if not self.export_formats:
            raise ValueError("At least one export format required")
        if self.image_width < 1 or self.image_height < 1:
            raise ValueError("Image dimensions must be >= 1")


@dataclass
class DataProcessingConfig:
    """Data processing configuration."""

    # Input validation
    validate_input: bool = True
    allowed_formats: Set[str] = field(default_factory=lambda: {"smiles", "inchi", "mol", "sdf"})
    max_input_size: int = 1000000  # bytes

    # Structure processing
    standardize_structures: bool = True
    generate_3d: bool = True
    optimize_3d: bool = True
    max_conformers: int = 100

    # Property calculation
    calc_descriptors: bool = True
    calc_fingerprints: bool = True
    descriptor_set: str = "rdkit"
    fingerprint_type: str = "morgan"
    fingerprint_radius: int = 2
    fingerprint_bits: int = 2048

    # Performance
    parallel_processing: bool = True
    max_parallel_jobs: int = -1  # -1 for all cores
    memory_limit: Optional[int] = None  # MB, None for no limit

    # Output validation
    validate_output: bool = True
    require_properties: Set[str] = field(default_factory=lambda: {"molecular_weight", "logp", "hbd", "hba", "tpsa"})

    def __post_init__(self):
        """Validate configuration."""
        if not self.allowed_formats:
            raise ValueError("At least one input format required")
        if self.max_input_size < 1:
            raise ValueError("max_input_size must be >= 1")
        if self.max_conformers < 1:
            raise ValueError("max_conformers must be >= 1")
        if self.fingerprint_radius < 1:
            raise ValueError("fingerprint_radius must be >= 1")
        if self.fingerprint_bits < 1:
            raise ValueError("fingerprint_bits must be >= 1")
        if not self.require_properties:
            raise ValueError("At least one required property needed")
        if self.max_parallel_jobs < -1:
            raise ValueError("max_parallel_jobs must be >= -1")
        if self.memory_limit is not None and self.memory_limit < 1:
            raise ValueError("memory_limit must be >= 1")


@dataclass
class DataSourceConfig:
    """Data source configuration."""

    # Source settings
    enabled_sources: Set[str] = field(default_factory=lambda: {"bindingdb", "pubchem", "chembl", "drugbank"})
    cache_results: bool = True
    cache_ttl: int = 86400  # 24 hours

    # Query settings
    max_results: int = 1000
    timeout: int = 30
    retry_count: int = 3
    batch_size: int = 100

    # Validation
    validate_input: bool = True
    validate_output: bool = True
    require_structures: bool = True

    def __post_init__(self):
        """Validate configuration."""
        if not self.enabled_sources:
            raise ValueError("At least one source must be enabled")
        if self.cache_ttl < 0:
            raise ValueError("cache_ttl must be >= 0")
        if self.max_results < 1:
            raise ValueError("max_results must be >= 1")
        if self.timeout < 1:
            raise ValueError("timeout must be >= 1")
        if self.retry_count < 0:
            raise ValueError("retry_count must be >= 0")
        if self.batch_size < 1:
            raise ValueError("batch_size must be >= 1")


@dataclass
class WebEnrichmentConfig:
    """Web enrichment configuration."""

    # Client settings
    user_agent: str = "ChemData/1.0"
    timeout: int = 30
    max_retries: int = 3
    concurrent_requests: int = 5

    # Rate limiting
    requests_per_second: float = 1.0
    burst_size: int = 10
    cooldown_period: int = 60

    # Caching
    cache_enabled: bool = True
    cache_ttl: int = 3600  # 1 hour
    cache_size: int = 10000  # Number of items

    # Content processing
    extract_tables: bool = True
    extract_images: bool = True
    follow_links: bool = False
    max_depth: int = 1

    def __post_init__(self):
        """Validate configuration."""
        if self.timeout < 1:
            raise ValueError("timeout must be >= 1")
        if self.max_retries < 0:
            raise ValueError("max_retries must be >= 0")
        if self.concurrent_requests < 1:
            raise ValueError("concurrent_requests must be >= 1")
        if self.requests_per_second <= 0:
            raise ValueError("requests_per_second must be > 0")
        if self.burst_size < 1:
            raise ValueError("burst_size must be >= 1")
        if self.cooldown_period < 0:
            raise ValueError("cooldown_period must be >= 0")
        if self.cache_ttl < 0:
            raise ValueError("cache_ttl must be >= 0")
        if self.cache_size < 1:
            raise ValueError("cache_size must be >= 1")
        if self.max_depth < 0:
            raise ValueError("max_depth must be >= 0")


@dataclass
class PsychopharmConfig:
    """Psychopharm-specific configuration."""

    # Prediction types
    prediction_types: Set[str] = field(default_factory=lambda: {"bbb", "abuse", "toxicity", "nootropic", "psychoactive"})
    confidence_threshold: float = 0.7
    require_explanation: bool = True

    # Model settings
    model_dir: Optional[Path] = None
    ensemble_size: int = 5
    use_gpu: bool = True

    # Validation
    validation_fraction: float = 0.2
    min_samples: int = 100
    cross_validation_folds: int = 5

    # Integration
    integrate_literature: bool = True
    integrate_social: bool = True
    integrate_clinical: bool = True

    def __post_init__(self):
        """Validate configuration."""
        if not self.prediction_types:
            raise ValueError("At least one prediction type required")
        if not 0 <= self.confidence_threshold <= 1:
            raise ValueError("confidence_threshold must be in [0, 1]")
        if self.ensemble_size < 1:
            raise ValueError("ensemble_size must be >= 1")
        if not 0 < self.validation_fraction < 1:
            raise ValueError("validation_fraction must be in (0, 1)")
        if self.min_samples < 1:
            raise ValueError("min_samples must be >= 1")
        if self.cross_validation_folds < 2:
            raise ValueError("cross_validation_folds must be >= 2")


@dataclass
class AppConfig:
    """Application configuration combining all components."""

    # Infrastructure
    cache: CacheConfig = field(default_factory=CacheConfig)
    rate_limit: RateLimitConfig = field(default_factory=RateLimitConfig)
    circuit_breaker: CircuitBreakerConfig = field(default_factory=CircuitBreakerConfig)
    monitor: MonitorConfig = field(default_factory=MonitorConfig)

    # Core processing
    structure_analysis: StructureAnalysisConfig = field(default_factory=StructureAnalysisConfig)
    ml_system: MLSystemConfig = field(default_factory=MLSystemConfig)
    data_processing: DataProcessingConfig = field(default_factory=DataProcessingConfig)
    document: DocumentConfig = field(default_factory=DocumentConfig)

    # Data sources and enrichment
    data_sources: DataSourceConfig = field(default_factory=DataSourceConfig)
    web_enrichment: WebEnrichmentConfig = field(default_factory=WebEnrichmentConfig)
    enrichment: EnrichmentConfig = field(default_factory=EnrichmentConfig)

    # Domain-specific
    psychopharm: PsychopharmConfig = field(default_factory=PsychopharmConfig)

    # Interface
    web: WebInterfaceConfig = field(default_factory=WebInterfaceConfig)

    # Global settings
    debug: bool = False
    log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "INFO"
    data_dir: Optional[Path] = None
    output_dir: Optional[Path] = None

    def __post_init__(self):
        """Validate configuration."""
        if isinstance(self.data_dir, str):
            self.data_dir = Path(self.data_dir)
        if isinstance(self.output_dir, str):
            self.output_dir = Path(self.output_dir)

    def validate_integration(self):
        """Validate integration between all components."""
        # Validate component integrations
        self.structure_analysis.validate_integration()
        self.ml_system.validate_integration()

        # Ensure cache settings are consistent
        if self.cache.backend == "file":
            if not self.data_dir:
                raise ValueError("data_dir required for file cache")
            self.cache.cache_dir = self.data_dir / "cache"

        # Validate rate limits
        if self.web.rate_limit > self.rate_limit.limit:
            raise ValueError("Web rate limit cannot exceed global rate limit")
        if self.enrichment.concurrent_requests > self.rate_limit.limit:
            raise ValueError("Enrichment concurrent requests cannot exceed rate limit")
        if self.web_enrichment.concurrent_requests > self.rate_limit.limit:
            raise ValueError("Web enrichment concurrent requests cannot exceed rate limit")

        # Ensure monitoring settings are consistent
        if self.monitor.metrics_enabled:
            if not self.monitor.storage_dir and self.monitor.backend == "file":
                self.monitor.storage_dir = self.data_dir / "metrics"

        # Validate ML and structure analysis compatibility
        if self.ml_system.ml.model_type == "gnn":
            if not self.structure_analysis.structure.optimize_conformers:
                raise ValueError("GNN models require optimized conformers")

        # Ensure psychopharm settings are compatible
        if self.psychopharm.use_gpu and not self.ml_system.shared_gpu:
            raise ValueError("Psychopharm GPU usage requires shared GPU setting")

        # Validate data processing settings
        if self.data_processing.parallel_processing:
            if self.data_processing.max_parallel_jobs > self.structure_analysis.max_parallel_jobs:
                self.data_processing.max_parallel_jobs = self.structure_analysis.max_parallel_jobs

        # Validate document processing settings
        if self.document.extract_images and not self.document.output_format == "json":
            raise ValueError("Image extraction requires JSON output format")

        # Ensure data source settings are compatible
        if self.data_sources.cache_results:
            if not self.cache.backend == "file":
                raise ValueError("Data source caching requires file cache backend")

        # Validate web enrichment settings
        if self.web_enrichment.cache_enabled:
            if not self.cache.backend == "file":
                raise ValueError("Web enrichment caching requires file cache backend")
