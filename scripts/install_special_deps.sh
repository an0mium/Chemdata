#!/bin/bash
set -e

echo "Installing special dependencies..."

# Create a temporary directory for building packages
TEMP_DIR=$(mktemp -d)
cd $TEMP_DIR

# Install MCP SDK
echo "Installing MCP SDK..."
npm install @modelcontextprotocol/sdk@1.0.4
cd node_modules/@modelcontextprotocol/sdk
pip install -e .
cd ../../..

# Install Ketcher and related packages
echo "Installing structure drawing packages..."

# 1. Ketcher Python
git clone https://github.com/epam/ketcher-python.git
cd ketcher-python
python setup.py install
cd ..

# 2. RDKit-Ketcher
git clone https://github.com/rdkit/rdkit-ketcher.git
cd rdkit-ketcher
python setup.py install
cd ..

# 3. Indigo-Ketcher
git clone https://github.com/epam/indigo-ketcher.git
cd indigo-ketcher
python setup.py install
cd ..

# 4. JSME Wrapper
git clone https://github.com/peter-ertl/jsme-wrapper.git
cd jsme-wrapper
python setup.py install
cd ..

# Clean up
cd ..
rm -rf $TEMP_DIR

echo "Special dependencies installation complete!"
