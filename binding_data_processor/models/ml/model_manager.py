"""Module for managing ML model initialization and access."""

import logging
import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM
from typing import Optional, Tuple

logger = logging.getLogger(__name__)


class MLModels:
    """Manages ML model initialization and access."""

    def __init__(self):
        self.tokenizer = None
        self.model = None
        self.device = None

    def initialize(self) -> bool:
        """Initialize ML models.

        Returns:
            bool: True if initialization successful, False otherwise
        """
        try:
            logger.info("Initializing scibert model...")

            # Set device - prefer MPS on Mac, then CUDA, then CPU
            if torch.backends.mps.is_available():
                device = torch.device("mps")
            elif torch.cuda.is_available():
                device = torch.device("cuda")
            else:
                device = torch.device("cpu")
            logger.info(f"Device set to use {device}")

            # Initialize tokenizer
            self.tokenizer = AutoTokenizer.from_pretrained("allenai/scibert_scivocab_uncased")

            # Initialize model for masked language modeling
            self.model = AutoModelForMaskedLM.from_pretrained("allenai/scibert_scivocab_uncased", trust_remote_code=True, from_pt=True)  # Load from PyTorch weights
            self.model = self.model.to(device)
            self.model.eval()  # Set to evaluation mode

            # Store device
            self.device = device

            # Verify model loaded successfully
            if self.model is None or self.tokenizer is None:
                raise RuntimeError("Failed to initialize models")

            logger.info("Model initialization successful")
            return True

        except Exception as e:
            logger.error(f"Error initializing models: {str(e)}")
            return False

    def is_initialized(self) -> bool:
        """Check if models are initialized.

        Returns:
            bool: True if models are initialized, False otherwise
        """
        return self.tokenizer is not None and self.model is not None and self.device is not None

    def get_device(self) -> Optional[torch.device]:
        """Get the current device.

        Returns:
            Optional[torch.device]: The current device (mps, cuda, or cpu), or None if not initialized
        """
        return self.device

    def get_models(self) -> Tuple[Optional[AutoTokenizer], Optional[AutoModelForMaskedLM]]:
        """Get the tokenizer and model.

        Returns:
            Tuple[Optional[AutoTokenizer], Optional[AutoModelForMaskedLM]]:
                The tokenizer and model, or (None, None) if not initialized
        """
        return self.tokenizer, self.model

    def cleanup(self) -> None:
        """Clean up resources."""
        if self.model is not None:
            # Move model to CPU before deletion to free GPU memory
            self.model.to("cpu")
            del self.model
            self.model = None

        if self.tokenizer is not None:
            del self.tokenizer
            self.tokenizer = None

        # Clear CUDA cache if available
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        logger.info("ML models cleaned up")
