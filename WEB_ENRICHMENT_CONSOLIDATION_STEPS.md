# Web Enrichment Consolidation Steps

## Overview

The web enrichment code currently has:
1. Base implementations in web_enrichment/
2. Enhanced versions with circuit breaker pattern
3. Transition to new clients/ directory structure

## Current Structure

### Base Implementations
```
web_enrichment/
├── http_client.py
├── base_client.py
├── swiss_client.py
├── community_client.py
├── social_client.py
└── manager.py
```

### Enhanced Versions
```
web_enrichment/
├── http_client_enhanced.py
├── swiss_client_enhanced.py
├── community_client_enhanced.py
├── social_client_enhanced.py
└── manager_enhanced.py
```

## Target Structure

```
web_enrichment/
├── clients/
│   ├── __init__.py
│   ├── base.py           # Base HTTP client
│   ├── swiss.py          # Swiss tools client
│   ├── community.py      # Community data client
│   └── social.py         # Social media client
├── validation/
│   ├── __init__.py
│   ├── schema.py         # Schema validation
│   └── data.py          # Data validation
├── processing/
│   ├── __init__.py
│   ├── nlp.py           # Text processing
│   └── trends.py        # Trend analysis
└── manager.py           # Enhanced manager
```

## Step-by-Step Plan

### 1. Create New Structure
```bash
# Create directories
mkdir -p web_enrichment/{clients,validation,processing}

# Create __init__.py files
touch web_enrichment/{clients,validation,processing}/__init__.py
```

### 2. Consolidate HTTP Client
```python
# In web_enrichment/clients/base.py
class HTTPClient:
    """HTTP client with circuit breaker and caching."""
    def __init__(self):
        self.circuit = CircuitBreaker()
        self.cache = TTLCache()
```

### 3. Consolidate Data Source Clients
```python
# In web_enrichment/clients/swiss.py
class SwissClient(BaseWebClient):
    """Swiss tools client with enhanced features."""
    def __init__(self):
        self.circuit = CircuitBreaker()
        self.validator = SchemaValidator()
```

### 4. Add Validation
```python
# In web_enrichment/validation/schema.py
class SchemaValidator:
    """Schema validation for web data."""
    def validate(self):
        pass

# In web_enrichment/validation/data.py
class DataValidator:
    """Data validation for web responses."""
    def validate(self):
        pass
```

### 5. Add Processing
```python
# In web_enrichment/processing/nlp.py
class TextProcessor:
    """NLP processing for web data."""
    def process(self):
        pass

# In web_enrichment/processing/trends.py
class TrendAnalyzer:
    """Trend analysis for web data."""
    def analyze(self):
        pass
```

### 6. Update Manager
```python
# In web_enrichment/manager.py
class WebEnrichmentManager:
    """Enhanced manager with all features."""
    def __init__(self):
        self.http = HTTPClient()
        self.validator = SchemaValidator()
        self.processor = TextProcessor()
```

### 7. Update Tests
```python
# In tests/web_enrichment/test_clients.py
def test_http_client():
    client = HTTPClient()
    assert client.circuit is not None

# In tests/web_enrichment/test_validation.py
def test_schema_validator():
    validator = SchemaValidator()
    assert validator.validate() is not None
```

## Migration Steps

### Day 1: Setup
1. Create new directory structure
2. Move HTTP client
3. Add validation framework

### Day 2: Clients
1. Move Swiss client
2. Move community client
3. Move social client

### Day 3: Processing
1. Add NLP processing
2. Add trend analysis
3. Add metrics collection

### Day 4: Integration
1. Update manager
2. Update imports
3. Add tests

## Success Criteria

### Code Quality
- [ ] All files under 500 lines
- [ ] Test coverage > 80%
- [ ] No duplicate code
- [ ] Consistent style

### Functionality
- [ ] All features preserved
- [ ] Enhanced validation
- [ ] Better error handling
- [ ] Improved metrics

### Documentation
- [ ] Updated docstrings
- [ ] Added examples
- [ ] Added guides
- [ ] Added API docs

## Next Steps

1. Create directory structure
2. Move and merge clients
3. Add validation framework
4. Add processing tools
5. Update manager
6. Add tests
7. Update documentation
