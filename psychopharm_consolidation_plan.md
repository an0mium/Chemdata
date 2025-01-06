# Psychopharm Model Consolidation Plan

## Current Status

### Completed Components ✓
1. Base Infrastructure ✓
   - HTTP Client with caching ✓
   - Rate limiting ✓
   - Error handling ✓
   - Response validation ✓

2. Scientific Sources ✓
   - BindingDB integration complete ✓
   - ChEMBL integration complete ✓
   - PubChem integration complete ✓
   - PubMed integration complete ✓
   - Swiss* services complete ✓

3. Patent Integration ✓
   - Structure search complete ✓
   - Family lookup complete ✓
   - Legal status tracking complete ✓
   - Analytics complete ✓

## Priority Components

### 1. Base Model Enhancement
```python
@dataclass
class PsychopharmBase:
    """Enhanced base class for psychopharmacological compounds."""
    
    identifiers: Dict[str, str] = field(default_factory=dict)
    properties: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
    community_data: Dict[str, Any] = field(default_factory=dict)
    safety_data: Dict[str, Any] = field(default_factory=dict)
    
    def validate(self) -> bool:
        """Enhanced validation with safety checks."""
        return all([
            self._validate_identifiers(),
            self._validate_properties(),
            self._validate_metadata(),
            self._validate_community_data(),
            self._validate_safety_data()
        ])
    
    def to_dict(self) -> Dict[str, Any]:
        """Enhanced serialization with all data."""
        return {
            "identifiers": self.identifiers,
            "properties": self.properties,
            "metadata": self.metadata,
            "predictions": self.get_all_predictions(),
            "web_data": self.get_all_web_data(),
            "community_data": self.community_data,
            "safety_data": self.safety_data
        }
```

### 2. Community Integration
```python
class CommunityEnrichment:
    """Enhanced community data integration."""
    
    async def enrich_from_community(self,
                                  sources: Optional[List[str]] = None
                                 ) -> Dict[str, Any]:
        """Enrich with community data."""
        sources = sources or [
            "reddit",
            "bluelight",
            "psychonautwiki",
            "erowid"
        ]
        
        async with aiohttp.ClientSession() as session:
            tasks = []
            for source in sources:
                task = self._fetch_community_data(session, source)
                tasks.append(task)
            
            results = await asyncio.gather(*tasks)
            analyzed = self._analyze_community_data(results)
            safety = self._assess_safety(analyzed)
            
            return {
                'data': analyzed,
                'safety': safety,
                'trends': self._detect_trends(analyzed),
                'alerts': self._generate_alerts(safety)
            }
```

### 3. Safety Analysis
```python
class SafetyAnalyzer:
    """Enhanced safety analysis with LLM support."""
    
    def analyze_safety(self,
                      analysis_types: Optional[List[str]] = None
                     ) -> Dict[str, Any]:
        """Comprehensive safety analysis."""
        analysis_types = analysis_types or [
            "toxicity",
            "interactions",
            "abuse_potential",
            "contraindications"
        ]
        
        results = {}
        for analysis_type in analysis_types:
            method = getattr(self, f"_analyze_{analysis_type}")
            results[analysis_type] = method()
            
        alerts = self._generate_alerts(results)
        recommendations = self._make_recommendations(results)
        
        return {
            'analysis': results,
            'alerts': alerts,
            'recommendations': recommendations,
            'risk_level': self._assess_risk_level(results)
        }
```

### 4. Content Analysis
```python
class ContentAnalyzer:
    """Enhanced content analysis with LLM support."""
    
    def analyze_content(self,
                       content_types: Optional[List[str]] = None
                      ) -> Dict[str, Any]:
        """Analyze content with LLM."""
        content_types = content_types or [
            "experience_reports",
            "safety_reports",
            "research_papers",
            "discussions"
        ]
        
        results = {}
        for content_type in content_types:
            method = getattr(self, f"_analyze_{content_type}")
            results[content_type] = method()
            
        entities = self._extract_entities(results)
        patterns = self._detect_patterns(results)
        insights = self._generate_insights(patterns)
        
        return {
            'analysis': results,
            'entities': entities,
            'patterns': patterns,
            'insights': insights
        }
```

## Implementation Steps

### Phase 1: Community Integration (Priority)
1. Reddit Integration
   - [ ] Complete OAuth flow
   - [ ] Add content monitoring
   - [ ] Add safety analysis
   - [ ] Add trend detection

2. Safety Analysis
   - [ ] Add content analysis
   - [ ] Add risk assessment
   - [ ] Add alert system
   - [ ] Add reporting

3. Content Analysis
   - [ ] Add LLM integration
   - [ ] Add trend detection
   - [ ] Add visualization
   - [ ] Add reporting

### Phase 2: Model Enhancement
1. Base Model
   - [ ] Add community fields
   - [ ] Add safety fields
   - [ ] Update validation
   - [ ] Update serialization

2. Analysis Tools
   - [ ] Add safety analysis
   - [ ] Add trend detection
   - [ ] Add visualization
   - [ ] Add reporting

### Phase 3: Integration
1. Data Integration
   - [ ] Add cross-validation
   - [ ] Add data merging
   - [ ] Add conflict resolution
   - [ ] Add reporting

2. System Integration
   - [ ] Add monitoring
   - [ ] Add logging
   - [ ] Add metrics
   - [ ] Add dashboards

## Testing Structure

### Unit Tests
```python
# tests/models/psychopharm/test_base.py
def test_enhanced_validation():
    """Test enhanced validation with safety."""
    compound = PsychopharmBase()
    compound.safety_data = {"risk_level": "low"}
    assert compound.validate()

# tests/models/psychopharm/test_community.py
def test_community_enrichment():
    """Test community data enrichment."""
    enricher = CommunityEnrichment()
    result = await enricher.enrich_from_community()
    assert "safety" in result
    assert "trends" in result

# tests/models/psychopharm/test_safety.py
def test_safety_analysis():
    """Test safety analysis."""
    analyzer = SafetyAnalyzer()
    result = analyzer.analyze_safety()
    assert "risk_level" in result
    assert "recommendations" in result
```

## Success Criteria

### Code Quality
- [x] No duplicate code in scientific integration ✓
- [x] Clear inheritance in scientific clients ✓
- [x] Complete type hints in scientific code ✓
- [x] Full docstrings in scientific code ✓
- [ ] No duplicate code in community integration
- [ ] Clear inheritance in community clients
- [ ] Complete type hints in community code
- [ ] Full docstrings in community code

### Functionality
- [x] Scientific integration complete ✓
- [x] Patent integration complete ✓
- [ ] Community integration working
- [ ] Safety analysis working
- [ ] Content analysis working
- [ ] Trend detection working

### Testing
- [x] Scientific tests passing ✓
- [x] Patent tests passing ✓
- [ ] Community tests passing
- [ ] Safety tests passing
- [ ] Content tests passing
- [ ] Integration tests passing

### Documentation
- [x] Scientific docs complete ✓
- [x] Patent docs complete ✓
- [ ] Community docs complete
- [ ] Safety docs complete
- [ ] Content docs complete
- [ ] Integration docs complete

## Next Steps

1. Complete Reddit Integration
   - Implement OAuth flow
   - Add content monitoring
   - Add safety analysis
   - Add trend detection

2. Add Safety Analysis
   - Implement content analysis
   - Add risk assessment
   - Add alert system
   - Add reporting

3. Enhance Content Analysis
   - Add LLM integration
   - Add trend detection
   - Add visualization
   - Add reporting

4. Update Documentation
   - Add community docs
   - Add safety docs
   - Add content docs
   - Add integration docs
