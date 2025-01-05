"""Setup script for binding_data_processor examples."""

from pathlib import Path
from setuptools import setup, find_namespace_packages


# Read requirements
def read_requirements(filename: str) -> list[str]:
    """Read requirements from file."""
    with open(filename) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


# Read README
readme = Path("README.md").read_text(encoding="utf-8")


setup(
    name="binding-data-processor-examples",
    version="0.1.0",
    description="Examples for binding_data_processor package",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/yourusername/binding-data-processor",
    packages=find_namespace_packages(include=["examples.*"]),
    python_requires=">=3.9",
    install_requires=read_requirements("requirements-dev.txt"),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    entry_points={
        "console_scripts": [
            "enrich-compounds=examples.scripts.enrich_compounds:main",
        ],
    },
)
