#!/bin/bash
set -e

# This script downloads and prepares the BindingDB data

# Ensure we're in the project root
cd "$(dirname "$0")/.."

# Create data directory if it doesn't exist
mkdir -p data

# Function to show progress
show_progress() {
    local pid=$1
    local delay=0.5
    local spinstr='|/-\'
    while [ "$(ps a | awk '{print $1}' | grep $pid)" ]; do
        local temp=${spinstr#?}
        printf " [%c]  " "$spinstr"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "\b\b\b\b\b\b"
    done
    printf "    \b\b\b\b"
}

# Function to download with progress and resume capability
download_with_progress() {
    local url=$1
    local output=$2
    echo "Downloading $(basename $output)..."
    if [ -f "$output.tmp" ]; then
        echo "Resuming previous download..."
        wget -c "$url" -O "$output.tmp" 2>&1 & show_progress $!
    else
        wget "$url" -O "$output.tmp" 2>&1 & show_progress $!
    fi
    mv "$output.tmp" "$output"
}

# Function to check file size
check_file_size() {
    local file=$1
    local min_size=$2  # in bytes
    if [ -f "$file" ]; then
        local size=$(stat -f%z "$file")
        if [ $size -ge $min_size ]; then
            return 0
        fi
    fi
    return 1
}

# Download BindingDB data
BINDINGDB_URL="https://bindingdb.org/bind/downloads/BindingDB_All_202501_tsv.zip"
BINDINGDB_ZIP="data/BindingDB_All.zip"
BINDINGDB_TSV="data/BindingDB_All.tsv"

if [ ! -f "$BINDINGDB_TSV" ] || ! check_file_size "$BINDINGDB_TSV" 1000000000; then
    echo "Downloading BindingDB data..."
    download_with_progress "$BINDINGDB_URL" "$BINDINGDB_ZIP"
    
    echo "Extracting BindingDB data..."
    unzip -o "$BINDINGDB_ZIP" -d data/
    rm "$BINDINGDB_ZIP"
else
    echo "BindingDB data already downloaded and extracted"
fi

# Create index for faster searching
echo "Creating search index..."
if command -v sqlite3 &> /dev/null; then
    echo "Creating SQLite database for faster searching..."
    sqlite3 data/bindingdb.db <<EOF
    DROP TABLE IF EXISTS compounds;
    CREATE TABLE compounds (
        id INTEGER PRIMARY KEY,
        smiles TEXT,
        inchi TEXT,
        inchi_key TEXT,
        monomer_id INTEGER,
        name TEXT,
        target_name TEXT,
        organism TEXT,
        ki REAL,
        ic50 REAL,
        kd REAL,
        ec50 REAL,
        doi TEXT,
        pubmed_id TEXT,
        patent_number TEXT
    );
    CREATE INDEX idx_smiles ON compounds(smiles);
    CREATE INDEX idx_inchi_key ON compounds(inchi_key);
    CREATE INDEX idx_name ON compounds(name);
    CREATE INDEX idx_target ON compounds(target_name);
    .mode tabs
    .import $BINDINGDB_TSV compounds
EOF
    echo "SQLite database created"
else
    echo "SQLite not found, skipping database creation"
    echo "Install SQLite for faster searching"
fi

echo "BindingDB data preparation complete!"
echo "Data file: $BINDINGDB_TSV"
if [ -f "data/bindingdb.db" ]; then
    echo "SQLite database: data/bindingdb.db"
fi
