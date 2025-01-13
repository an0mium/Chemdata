#!/bin/bash

# Setup script for Chemdata database
# This script will:
# 1. Create a new PostgreSQL database
# 2. Apply the consolidated schema
# 3. Add reference data
# 4. Verify installation
# 5. Create backup

set -e # Exit on error

# Configuration
DB_NAME="chemdata"
DB_USER="$USER" # Current user
SCHEMA_DIR="database/schema/consolidated"
REFERENCE_DATA_DIR="database/schema/reference_data"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging functions
log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')] $1${NC}"
}

error() {
    echo -e "${RED}[$(date +'%Y-%m-%d %H:%M:%S')] ERROR: $1${NC}"
    exit 1
}

warn() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')] WARNING: $1${NC}"
}

# Function to execute SQL with proper error handling
execute_sql() {
    if ! psql -v ON_ERROR_STOP=1 "$@"; then
        error "Failed to execute SQL command"
    fi
}

# Check if PostgreSQL is installed
if ! command -v psql &> /dev/null; then
    error "PostgreSQL is not installed. Please install PostgreSQL first."
fi

# Check if database exists
if psql -lqt | cut -d \| -f 1 | grep -qw "$DB_NAME"; then
    warn "Database $DB_NAME already exists."
    read -p "Do you want to drop and recreate it? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        log "Dropping existing database..."
        dropdb "$DB_NAME" || error "Failed to drop database"
    else
        error "Aborting to prevent data loss."
    fi
fi

# Create database
log "Creating database $DB_NAME..."
createdb "$DB_NAME" || error "Failed to create database"

# Apply consolidated schema
log "Applying consolidated schema..."
execute_sql "$DB_NAME" -f "$SCHEMA_DIR/00_main.sql"

# Add reference data
log "Adding reference data..."
for sql_file in "$REFERENCE_DATA_DIR"/*.sql; do
    if [ -f "$sql_file" ]; then
        log "Loading reference data from $(basename "$sql_file")..."
        execute_sql "$DB_NAME" -f "$sql_file"
    fi
done

# Verify installation
log "Verifying schema installation..."

# Core tables to verify
TABLES=(
    "compounds"
    "molecular_descriptors"
    "binding_data"
    "quantum_properties"
    "receptor_families"
    "toxicity_endpoints"
    "therapeutic_classes"
    "research_findings"
    "ml_models"
    "web_templates"
)

for table in "${TABLES[@]}"; do
    if ! psql -tAc "SELECT to_regclass('public.$table');" "$DB_NAME" | grep -q "$table"; then
        error "Table $table was not created properly"
    fi
done

# Create backup
BACKUP_FILE="${DB_NAME}_initial_backup_$(date +%Y%m%d_%H%M%S).sql"
log "Creating initial backup..."
if ! pg_dump "$DB_NAME" > "$BACKUP_FILE"; then
    warn "Failed to create backup file $BACKUP_FILE"
else
    log "Backup created: $BACKUP_FILE"
fi

# Print summary
echo
log "Setup Summary:"
echo "-------------"
echo "Database name: $DB_NAME"
echo "User: $DB_USER"
echo "Backup file: $BACKUP_FILE"
echo
echo "Tables created:"
psql -d "$DB_NAME" -c "\dt" | tail -n +3
echo
echo "To start using the database:"
echo "1. Connect: psql $DB_NAME"
echo "2. List tables: \dt"
echo "3. Describe table: \d table_name"
echo

log "Database setup completed successfully!"
