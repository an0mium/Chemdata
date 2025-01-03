#!/bin/bash
set -e

# Install uv if not already installed
if ! command -v uv &> /dev/null; then
    echo "Installing uv package manager..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
fi

# Add uv to PATH if not already there
if [[ ":$PATH:" != *":$HOME/.local/bin:"* ]]; then
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc
    source ~/.zshrc
fi

# Create and activate virtual environment
echo "Creating virtual environment..."
rm -rf .venv
uv venv -p python3.12 .venv

# Activate virtual environment
source .venv/bin/activate

# Install standard dependencies with uv
echo "Installing standard dependencies..."
uv pip install -r requirements.txt

# Make special deps script executable
chmod +x scripts/install_special_deps.sh

# Install special dependencies
echo "Installing special dependencies..."
./scripts/install_special_deps.sh

# Install package in development mode
echo "Installing package in development mode..."
uv pip install -e .

# Install pre-commit hooks
echo "Setting up pre-commit hooks..."
cat > .pre-commit-config.yaml << EOL
repos:
-   repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v4.5.0
    hooks:
    -   id: trailing-whitespace
    -   id: end-of-file-fixer
    -   id: check-yaml
    -   id: check-added-large-files

-   repo: https://github.com/psf/black
    rev: 23.12.1
    hooks:
    -   id: black

-   repo: https://github.com/pycqa/isort
    rev: 5.13.2
    hooks:
    -   id: isort

-   repo: https://github.com/pycqa/flake8
    rev: 7.0.0
    hooks:
    -   id: flake8
EOL

uv pip install pre-commit
pre-commit install

echo "Development environment setup complete!"
echo "To activate the environment, run: source .venv/bin/activate"
