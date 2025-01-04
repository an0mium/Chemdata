# Contributing to ChemData

Thank you for your interest in contributing to ChemData! This document provides guidelines and instructions for contributing to the project.

## Code of Conduct

This project and everyone participating in it is governed by our Code of Conduct. By participating, you are expected to uphold this code.

## Getting Started

1. Fork the repository
2. Clone your fork:
```bash
git clone https://github.com/yourusername/chemdata.git
cd chemdata
```

3. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

4. Install development dependencies:
```bash
# Install base requirements
pip install -r requirements-dev.txt

# Install special dependencies (ML libraries, chemical toolkits)
./scripts/install_special_deps.sh

# Install optional dependencies for development
./scripts/setup_dev.sh --all
```

5. Set up pre-commit hooks:
```bash
pre-commit install
```

6. Configure data sources:
```bash
# Set up API keys and credentials
cp .env.example .env
# Edit .env with your API keys

# Download initial datasets
./scripts/get_bindingdb.sh
```

## Development Process

1. Create a new branch for your feature:
```bash
git checkout -b feature/your-feature-name
```

2. Make your changes, following our coding standards

3. Write or update tests as needed

4. Run the test suite:
```bash
# Run all tests
pytest

# Run specific test categories
pytest tests/test_ml_predictions.py  # ML tests
pytest tests/test_web_enrichment.py  # Web enrichment tests
pytest tests/test_pipeline.py        # Pipeline tests
```

5. Run code quality checks:
```bash
# Style checks
flake8 binding_data_processor tests

# Type checking
mypy binding_data_processor

# Security checks
bandit -r binding_data_processor

# ML model validation
./scripts/validate_models.sh
```

6. Commit your changes:
```bash
git add .
git commit -m "feat: add your feature description"
```

7. Push to your fork:
```bash
git push origin feature/your-feature-name
```

8. Submit a pull request

## Coding Standards

### Python Style Guide

- Follow PEP 8 style guide
- Use type hints for function arguments and return values
- Maximum line length is 100 characters
- Use double quotes for strings
- Sort imports using isort
- Format code using black
- Keep functions focused and under 50 lines
- Use descriptive variable names
- Add comments for complex algorithms
- Break files longer than 700 lines into modules

### Documentation

- Use Google-style docstrings
- Document all public functions, classes, and methods
- Include examples in docstrings where appropriate
- Update README.md if adding new features
- Add doctest examples for API functions
- Document ML model parameters and assumptions
- Include performance metrics for ML models
- Add data validation rules to docstrings
- Document API rate limits and requirements

### Testing

- Write unit tests for all new functionality
- Maintain test coverage above 80%
- Use pytest fixtures for test setup
- Mock external API calls in tests
- Include integration tests for complex features
- Add ML model validation tests
- Test edge cases and error conditions
- Include performance benchmarks
- Test data validation rules
- Verify API error handling

### ML Model Development

- Document model architecture and hyperparameters
- Include training and validation metrics
- Save model checkpoints
- Version training data
- Document feature engineering steps
- Include uncertainty estimates
- Test model robustness
- Validate predictions
- Monitor for drift
- Include interpretability analysis

### Commit Messages

Follow the Conventional Commits specification:

- feat: New feature
- fix: Bug fix
- docs: Documentation changes
- style: Code style changes (formatting, etc)
- refactor: Code refactoring
- test: Test updates
- chore: Maintenance tasks
- ml: Machine learning updates
- data: Data processing changes
- web: Web interface updates

Example:
```
feat(ml): add receptor binding prediction model

Add deep learning model for predicting receptor binding profiles:
- Implement graph neural network architecture
- Add feature extraction pipeline
- Include uncertainty estimation
- Add model validation suite
- Document performance metrics

Closes #456
```

## Project Structure

```
binding_data_processor/
├── data_sources/              # Data source integrations
│   ├── bindingdb.py
│   ├── chembl.py
│   └── pubchem.py
├── models/                    # Data models and ML
│   ├── compound/             # Compound data models
│   └── psychopharm/          # Psychopharmacology models
├── pipeline/                  # Processing pipeline
│   ├── analysis/             # Analysis components
│   ├── enrichment/           # Data enrichment
│   └── validation/           # Data validation
├── processors/               # Data processors
│   ├── structure/            # Structure processing
│   └── psychopharm/          # Psychopharm analysis
├── web_enrichment/          # Web data enrichment
│   ├── community/            # Community data
│   └── social/               # Social media
└── web/                     # Web interface
    ├── api/                  # REST API
    ├── components/           # UI components
    └── pages/                # Web pages
```

### Module Guidelines

1. Data Sources
- One module per external data source
- Include rate limiting and error handling
- Cache responses when appropriate
- Document API requirements
- Handle authentication securely
- Validate responses
- Retry on failures
- Log API usage

2. Models
- Keep models focused and single-purpose
- Use dataclasses or Pydantic models
- Include validation
- Document all fields
- Add type hints
- Include serialization
- Handle versioning
- Add data migrations

3. ML Models
- Document architecture
- Include training pipeline
- Add validation suite
- Monitor performance
- Handle versioning
- Include interpretability
- Document limitations
- Add uncertainty estimates

4. Processors
- Make processors configurable
- Support batch processing
- Include progress reporting
- Handle errors gracefully
- Add validation
- Support cancellation
- Include logging
- Monitor resources

5. Web Interface
- Follow React-like component structure
- Use TypeScript for frontend code
- Make components reusable
- Include responsive design
- Add error boundaries
- Include loading states
- Support offline mode
- Add analytics

## Pull Request Process

1. Update documentation
2. Add or update tests
3. Run full test suite
4. Update CHANGELOG.md
5. Request review from maintainers
6. Address review feedback
7. Ensure CI passes
8. Squash commits if requested

## Release Process

1. Update version in setup.py
2. Update CHANGELOG.md
3. Create release branch
4. Run full test suite
5. Create GitHub release
6. Upload to PyPI
7. Update documentation
8. Notify users

## Getting Help

- Open an issue for bugs
- Use discussions for questions
- Join our Discord server
- Check the wiki for guides
- Read the API documentation
- Review example notebooks
- Check FAQ section

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
