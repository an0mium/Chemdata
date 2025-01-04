# Use multi-stage build for smaller final image
FROM python:3.8-slim as builder

# Install system dependencies required for building packages
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libboost-all-dev \
    libopenbabel-dev \
    librdkit-dev \
    python3-dev \
    graphviz \
    libgraphviz-dev \
    pkg-config \
    libcairo2-dev \
    libjpeg-dev \
    libgif-dev \
    libpng-dev \
    libtiff-dev \
    libxml2-dev \
    libxslt-dev \
    libffi-dev \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libreadline-dev \
    libsqlite3-dev \
    libfreetype6-dev \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install uv
RUN curl -LsSf https://astral.sh/uv/install.sh | sh

# Create virtual environment
RUN /root/.cargo/bin/uv venv /opt/venv

# Install Python dependencies
COPY requirements.txt pyproject.toml setup.cfg ./
RUN /root/.cargo/bin/uv pip install --no-cache -r requirements.txt

# Install additional dependencies
RUN /root/.cargo/bin/uv pip install --no-cache \
    rdkit \
    openbabel-wheel \
    epam.indigo \
    torch torchvision torchaudio \
    "dgl>=1.0.0" \
    "dgllife>=0.3.0" \
    "deepchem>=2.7.0" \
    "mordred>=1.2.0" \
    "descriptastorus>=2.6.0" \
    "scikit-learn>=1.0.0" \
    "xgboost>=1.7.0" \
    "lightgbm>=4.0.0" \
    "optuna>=3.0.0" \
    "ray[tune]>=2.0.0" \
    "plotly>=5.0.0" \
    "dash>=2.0.0" \
    "dash-bio>=1.0.0" \
    "dash-bootstrap-components>=1.0.0" \
    "networkx>=3.0" \
    "graphviz>=0.20.0" \
    "py3Dmol>=2.0.0" \
    "nglview>=3.0.0" \
    "selenium>=4.0.0" \
    "beautifulsoup4>=4.10.0" \
    "requests>=2.28.0" \
    "aiohttp>=3.8.0" \
    "playwright>=1.30.0"

# Final stage
FROM python:3.8-slim

# Install system dependencies required at runtime
RUN apt-get update && apt-get install -y \
    libopenbabel-dev \
    librdkit-dev \
    graphviz \
    libgraphviz-dev \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Copy virtual environment from builder stage
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Create non-root user
RUN useradd -m -u 1000 chemdata
USER chemdata

# Set working directory
WORKDIR /app

# Copy application code
COPY --chown=chemdata:chemdata . .

# Create necessary directories with correct permissions
RUN mkdir -p \
    /app/data/{raw,processed,interim,external} \
    /app/models/{toxicity,abuse,activity,affinity} \
    /app/logs \
    /app/cache \
    /app/output \
    /app/reports/{coverage,test-results,profiling} \
    /app/web/{static,templates} \
    && chown -R chemdata:chemdata \
        /app/data \
        /app/models \
        /app/logs \
        /app/cache \
        /app/output \
        /app/reports \
        /app/web

# Environment variables
ENV PYTHONPATH=/app
ENV PYTHONUNBUFFERED=1
ENV CHEMDATA_CACHE_DIR=/app/cache
ENV CHEMDATA_DATA_DIR=/app/data
ENV CHEMDATA_LOG_DIR=/app/logs
ENV CHEMDATA_OUTPUT_DIR=/app/output
ENV CHEMDATA_MODEL_DIR=/app/models
ENV FLASK_APP=binding_data_processor.web.app
ENV FLASK_ENV=production

# Expose ports
EXPOSE 8000  # Web interface
EXPOSE 8080  # API

# Health check
HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

# Default command
CMD ["python", "-m", "binding_data_processor.web.app"]

# Usage instructions in comments:
# Build:
#   docker build -t chemdata .
#
# Run web interface:
#   docker run -p 8000:8000 -p 8080:8080 \
#     -v $(pwd)/data:/app/data \
#     -v $(pwd)/logs:/app/logs \
#     -v $(pwd)/cache:/app/cache \
#     -v $(pwd)/output:/app/output \
#     -v $(pwd)/models:/app/models \
#     --env-file .env \
#     chemdata
#
# Run pipeline:
#   docker run --rm \
#     -v $(pwd)/data:/app/data \
#     -v $(pwd)/output:/app/output \
#     -v $(pwd)/models:/app/models \
#     chemdata python -m binding_data_processor.cli process-all \
#     --input data/raw/BindingDB_All.tsv \
#     --output data/processed/compounds.tsv
#
# Development shell:
#   docker run -it --rm \
#     -v $(pwd):/app \
#     chemdata bash
