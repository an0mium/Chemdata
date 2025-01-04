# Model Consolidation Plan

## Current Structure

### Psychopharm Models (Target)
- psychopharm/base.py - Core base mixin
- psychopharm/binding.py - Receptor binding profiles
- psychopharm/activity.py - Activity analysis
- psychopharm/safety.py - Safety assessment
- psychopharm/enrichment.py - Web data enrichment
- psychopharm/compound.py - Main compound class

### Legacy Models (To Consolidate)
- compound.py - Duplicate compound model
- compound_base.py - Duplicate base functionality
- compound_ml.py - ML prediction functionality
- compound_enrichment.py - Duplicate enrichment
- compound_analysis.py - Duplicate analysis

## Consolidation Steps

### 1. Base Functionality
- Move core dataclass fields from compound_base.py to psychopharm/base.py
- Merge validation methods
- Update type hints and docstrings
- Remove compound_base.py

### 2. ML Integration
- Move ML prediction methods from compound_ml.py to psychopharm/binding.py
- Enhance receptor profile with ML capabilities
- Add uncertainty estimation
- Remove compound_ml.py

### 3. Web Enrichment
- Move enrichment methods from compound_enrichment.py to psychopharm/enrichment.py
- Enhance web data processing
- Add new data sources
- Remove compound_enrichment.py

### 4. Analysis Capabilities
- Move analysis methods from compound_analysis.py to:
  - psychopharm/activity.py for activity analysis
  - psychopharm/safety.py for safety analysis
- Enhance prediction capabilities
- Remove compound_analysis.py

### 5. Main Compound Class
- Move any unique functionality from compound.py to psychopharm/compound.py
- Update class inheritance
- Enhance documentation
- Remove compound.py

## Implementation Order

1. Base Functionality
- [ ] Audit base.py and compound_base.py
- [ ] Identify unique features
- [ ] Merge functionality
- [ ] Update tests
- [ ] Remove old file

2. ML Integration
- [ ] Audit binding.py and compound_ml.py
- [ ] Merge prediction capabilities
- [ ] Enhance receptor profiles
- [ ] Update tests
- [ ] Remove old file

3. Web Enrichment
- [ ] Audit enrichment.py and compound_enrichment.py
- [ ] Merge enrichment features
- [ ] Add new data sources
- [ ] Update tests
- [ ] Remove old file

4. Analysis Capabilities
- [ ] Audit activity.py, safety.py, and compound_analysis.py
- [ ] Split analysis features
- [ ] Enhance predictions
- [ ] Update tests
- [ ] Remove old file

5. Main Class
- [ ] Audit compound.py files
- [ ] Merge unique features
- [ ] Update inheritance
- [ ] Update tests
- [ ] Remove old file

## Testing Strategy

1. Create Backup
- [ ] Branch: backup/models-original
- [ ] Commit all current files

2. Feature Branch
- [ ] Branch: feature/model-consolidation
- [ ] Implement changes incrementally

3. Testing
- [ ] Run existing tests after each merge
- [ ] Add new tests for enhanced features
- [ ] Verify no functionality loss

4. Documentation
- [ ] Update API documentation
- [ ] Add migration guide
- [ ] Update examples

## Validation

For each consolidated file:
1. Verify all functionality preserved
2. Run full test suite
3. Check import statements
4. Update documentation
5. Remove old file

## Next Steps

1. Create backup branch
2. Create feature branch
3. Start with base functionality
4. Follow implementation order
5. Run tests frequently
6. Update documentation
7. Clean up old files
