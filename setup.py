#!/usr/bin/env python
"""Setup script for binding_data_processor package."""

import os
from setuptools import setup, find_packages

# Read version from __init__.py
with open(os.path.join("binding_data_processor", "__init__.py"), "r") as f:
    for line in f:
        if line.startswith("__version__"):
            version = line.split("=")[1].strip().strip('"').strip("'")
            break

# Read long description from README
with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

# Core dependencies
CORE_DEPS = [
    "numpy>=1.21.0",
    "pandas>=1.3.0",
    "scipy>=1.7.0",
    "scikit-learn>=1.0.0",
    "pyarrow>=6.0.0",
]

# Chemistry tools
CHEM_DEPS = [
    "rdkit>=2022.3.1",
    "mordred>=1.2.0",
    "openbabel>=3.1.0",
]

# Machine learning
ML_DEPS = [
    "torch>=1.10.0",
    "torch-geometric>=2.0.0",
    "transformers>=4.15.0",
    "pytorch-lightning>=1.5.0",
    "xgboost>=1.5.0",
    "lightgbm>=3.3.0",
    "optuna>=2.10.0",
]

# Web scraping and APIs
WEB_DEPS = [
    "requests>=2.26.0",
    "beautifulsoup4>=4.10.0",
    "selenium>=4.1.0",
    "aiohttp>=3.8.0",
    "praw>=7.5.0",
    "tweepy>=4.4.0",
    "discord.py>=2.0.0",
    "atproto>=0.0.20",
]

# Web interface
WEB_INTERFACE_DEPS = [
    "dash>=2.0.0",
    "dash-bio>=1.0.0",
    "dash-bootstrap-components>=1.0.0",
    "plotly>=5.3.0",
    "flask>=2.0.0",
    "flask-caching>=1.10.0",
    "fastapi>=0.70.0",
    "uvicorn>=0.15.0",
    "pydantic>=1.8.0",
    "starlette>=0.16.0",
]

# Data processing
DATA_DEPS = [
    "tqdm>=4.62.0",
    "joblib>=1.1.0",
    "dask>=2022.1.0",
    "distributed>=2022.1.0",
    "networkx>=2.6.0",
    "tenacity>=8.0.1",
    "cachetools>=4.2.4",
]

# Visualization
VIZ_DEPS = [
    "matplotlib>=3.5.0",
    "seaborn>=0.11.0",
]

# Patent analysis
PATENT_DEPS = [
    "patent-client>=2.2.0",
    "scispacy>=0.5.1",
    "chemdataextractor>=2.1.0",
]

# Development dependencies
DEV_DEPS = [
    "pytest>=6.2.5",
    "pytest-cov>=2.12.0",
    "pytest-asyncio>=0.16.0",
    "pytest-mock>=3.6.1",
    "black>=21.12b0",
    "flake8>=4.0.1",
    "mypy>=0.910",
    "isort>=5.10.1",
    "pre-commit>=2.16.0",
    "bandit>=1.7.0",
]

# Documentation dependencies
DOC_DEPS = [
    "sphinx>=4.3.0",
    "sphinx-rtd-theme>=1.0.0",
    "sphinx-autodoc-typehints>=1.12.0",
    "nbsphinx>=0.8.0",
    "jupyter>=1.0.0",
]

# Type checking dependencies
TYPE_DEPS = [
    "types-requests>=2.26.1",
    "types-beautifulsoup4>=4.9.3",
    "types-tqdm>=4.62.2",
    "types-PyYAML>=6.0.1",
]

# Utility dependencies
UTIL_DEPS = [
    "python-dotenv>=0.19.0",
    "pyyaml>=6.0.0",
    "click>=8.0.0",
    "loguru>=0.5.0",
    "rich>=10.16.0",
]

setup(
    name="binding_data_processor",
    version=version,
    description="A comprehensive pipeline for processing and analyzing chemical compound data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/yourusername/chemdata",
    packages=find_packages(exclude=["tests*", "docs*"]),
    python_requires=">=3.8",
    install_requires=CORE_DEPS,
    extras_require={
        "chemistry": CHEM_DEPS,
        "ml": ML_DEPS,
        "web": WEB_DEPS,
        "web-interface": WEB_INTERFACE_DEPS,
        "data": DATA_DEPS,
        "viz": VIZ_DEPS,
        "patent": PATENT_DEPS,
        "dev": DEV_DEPS,
        "docs": DOC_DEPS,
        "types": TYPE_DEPS,
        "utils": UTIL_DEPS,
        "all": (
            CHEM_DEPS + ML_DEPS + WEB_DEPS + WEB_INTERFACE_DEPS +
            DATA_DEPS + VIZ_DEPS + PATENT_DEPS + UTIL_DEPS
        ),
        "complete": (
            CHEM_DEPS + ML_DEPS + WEB_DEPS + WEB_INTERFACE_DEPS +
            DATA_DEPS + VIZ_DEPS + PATENT_DEPS + DEV_DEPS +
            DOC_DEPS + TYPE_DEPS + UTIL_DEPS
        ),
    },
    entry_points={
        "console_scripts": [
            "chemdata=binding_data_processor.cli:main",
            "chemdata-web=binding_data_processor.web.app:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    keywords=[
        "chemistry",
        "machine learning",
        "data processing",
        "binding data",
        "psychopharmacology",
        "drug discovery",
    ],
    project_urls={
        "Documentation": "https://chemdata.readthedocs.io/",
        "Source": "https://github.com/yourusername/chemdata",
        "Issues": "https://github.com/yourusername/chemdata/issues",
    },
    include_package_data=True,
    zip_safe=False,
)
