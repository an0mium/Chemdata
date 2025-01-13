"""AlphaFold configuration."""

from dataclasses import dataclass
from typing import Optional, Dict, Any
from pathlib import Path


@dataclass
class AlphaFoldConfig:
    """Configuration for AlphaFold prediction."""

    # Model paths
    model_dir: Optional[str] = None
    data_dir: Optional[str] = None
    params_dir: Optional[str] = None

    # Runtime settings
    use_gpu: bool = True
    num_ensemble: int = 1
    num_recycle: int = 3
    max_seq_len: int = 512

    # MSA options
    use_templates: bool = True
    max_templates: int = 4

    # Performance
    batch_size: int = 1
    chunk_size: int = 128

    # Output options
    save_raw: bool = False
    save_debug: bool = False

    # API settings
    api_url: str = "https://api.alphafold.ebi.ac.uk/v1"
    api_key: Optional[str] = None
    timeout: int = 3600

    # Cache settings
    cache_dir: Optional[Path] = None
    use_cache: bool = True

    def __post_init__(self):
        """Convert string paths to Path objects."""
        if self.cache_dir and isinstance(self.cache_dir, str):
            self.cache_dir = Path(self.cache_dir)
