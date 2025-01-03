# ChemData

A comprehensive chemical data collection and analysis tool that gathers information from multiple sources including PubChem, ChEMBL, PubMed, patent databases, and more.

## Features

- Collects chemical compound data from multiple authoritative sources
- Provides a web interface for browsing and searching compounds
- Includes structure drawing and similarity search capabilities
- Enriches data with pharmacological, toxicological, and legal information
- Supports export to TSV format with comprehensive column fields

## Installation

ChemData uses the [uv](https://github.com/astral/uv) package manager for fast, reliable dependency management. uv is significantly faster than pip and provides better dependency resolution.

### Quick Start

```bash
# Clone the repository
git clone https://github.com/yourusername/chemdata.git
cd chemdata

# Run the setup script (installs uv and sets up the development environment)
./scripts/setup_dev.sh
```

### Manual Setup

If you prefer to set up manually:

1. Install uv:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

2. Create and activate a virtual environment:
```bash
uv venv -p python3.12 .venv
source .venv/bin/activate
```

3. Install dependencies:
```bash
uv pip install -e ".[dev]"  # Includes development dependencies
# or
uv pip install -e .         # Production dependencies only
```

## Development Scripts

ChemData includes several utility scripts to help with development:

### Development Environment
- `./scripts/dev.sh`: Start the development environment with auto-reloading
- `./scripts/setup_dev.sh`: Set up the development environment
- `./scripts/update_deps.sh`: Update dependencies to their latest compatible versions

### Testing
- `./scripts/test.sh`: Run tests
  - `--coverage`: Generate coverage report
  - `--watch`: Watch for changes and run tests automatically

### Maintenance
- `./scripts/clean.sh`: Clean up development artifacts
  - `--all`: Clean everything
  - `--venv`: Clean virtual environment only
  - `--cache`: Clean uv cache only
  - `--build`: Clean build artifacts only

### Examples
- `./scripts/example.sh`: Run example workflows
- `./scripts/benchmark.sh`: Compare uv performance against pip

## Usage

### Web Interface

Start the web interface:
```bash
python -m chemdata.web
```

This will launch a web server at http://localhost:5000 where you can:
- Browse the compound database
- Search by name, structure, or properties
- View detailed compound information
- Export search results

### Command Line Interface

Process chemical data:
```bash
python -m chemdata.cli --target "5-HT2A" --output results.tsv
```

Available options:
- `--target`: Target receptor name (e.g., "5-HT2A")
- `--compounds`: File containing list of compounds
- `--chemical-class`: Chemical class to search for
- `--output`: Output TSV file path
- `--max-compounds`: Maximum number of compounds to process
- `--skip-sources`: Skip specific data sources
- `--verbose`: Enable verbose logging

## Development

### Environment

The development environment includes:
- pytest for testing
- black for code formatting
- isort for import sorting
- mypy for type checking
- flake8 for linting
- pre-commit hooks for code quality

### Running Tests

```bash
./scripts/test.sh
```

### Code Quality

The project uses pre-commit hooks to maintain code quality. They are automatically installed by the setup script, but you can manually install them:

```bash
uv pip install pre-commit
pre-commit install
```

## Configuration

Environment variables:
- `OPENAI_API_KEY`: OpenAI/OpenRouter API key for enhanced data extraction
- `SERP_API_KEY`: SerpAPI key for web searching

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Run tests (`./scripts/test.sh`)
4. Commit your changes (`git commit -m 'Add amazing feature'`)
5. Push to the branch (`git push origin feature/amazing-feature`)
6. Open a Pull Request

## Performance

ChemData uses uv for dependency management, which provides several benefits:
- Faster dependency resolution
- Better caching
- More reliable installations
- Improved reproducibility

You can benchmark uv's performance against pip:
```bash
./scripts/benchmark.sh
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
