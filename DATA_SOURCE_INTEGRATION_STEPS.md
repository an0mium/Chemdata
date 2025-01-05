# Data Source Integration Steps

## Overview

The project needs to integrate multiple data sources:
1. BindingDB (completed)
2. ChEMBL API
3. PubChem API
4. Community Sources (PsychonautWiki, Erowid, TripSit)
5. Social Media (Reddit, Twitter)

## Current Structure

```
data_sources/
└── bindingdb.py    # BindingDB processing
```

## Target Structure

```
data_sources/
├── core/
│   ├── __init__.py
│   ├── base.py        # Base client
│   ├── cache.py       # Caching
│   └── rate.py        # Rate limiting
├── scientific/
│   ├── __init__.py
│   ├── bindingdb.py   # BindingDB
│   ├── chembl.py      # ChEMBL
│   └── pubchem.py     # PubChem
├── community/
│   ├── __init__.py
│   ├── psychonaut.py  # PsychonautWiki
│   ├── erowid.py      # Erowid
│   └── tripsit.py     # TripSit
└── social/
    ├── __init__.py
    ├── reddit.py      # Reddit API
    └── twitter.py     # Twitter API
```

## Step-by-Step Plan

### 1. Core Infrastructure

```python
# In data_sources/core/base.py
class BaseClient:
    """Base API client with shared functionality."""
    def __init__(self):
        self.cache = Cache()
        self.rate_limiter = RateLimiter()
        self.session = Session()

    async def get(self, url: str, **params) -> Dict:
        """Make rate-limited GET request with caching."""
        cache_key = self._make_cache_key(url, params)
        
        # Check cache
        if cached := await self.cache.get(cache_key):
            return cached
            
        # Rate limit
        await self.rate_limiter.acquire()
        
        try:
            # Make request
            async with self.session.get(url, params=params) as response:
                data = await response.json()
                
            # Cache response
            await self.cache.set(cache_key, data)
            
            return data
            
        finally:
            self.rate_limiter.release()
```

### 2. Scientific Sources

```python
# In data_sources/scientific/chembl.py
class ChEMBLClient(BaseClient):
    """ChEMBL API client."""
    def __init__(self):
        super().__init__()
        self.base_url = "https://www.ebi.ac.uk/chembl/api/data"
        
    async def get_compound(self, chembl_id: str) -> CompoundData:
        """Get compound data from ChEMBL."""
        data = await self.get(f"{self.base_url}/molecule/{chembl_id}")
        return self._parse_compound(data)
        
    async def search_compounds(self, query: str) -> List[CompoundData]:
        """Search compounds by text."""
        data = await self.get(
            f"{self.base_url}/molecule/search",
            q=query
        )
        return [self._parse_compound(item) for item in data["molecules"]]
```

### 3. Community Sources

```python
# In data_sources/community/psychonaut.py
class PsychonautClient(BaseClient):
    """PsychonautWiki API client."""
    def __init__(self):
        super().__init__()
        self.base_url = "https://api.psychonautwiki.org"
        
    async def get_substance(self, name: str) -> Dict:
        """Get substance data from PsychonautWiki."""
        data = await self.get(
            f"{self.base_url}/substances/search",
            name=name
        )
        return self._parse_substance(data)
```

### 4. Social Sources

```python
# In data_sources/social/reddit.py
class RedditClient(BaseClient):
    """Reddit API client."""
    def __init__(self):
        super().__init__()
        self.base_url = "https://oauth.reddit.com"
        self._authenticate()
        
    async def search_subreddit(
        self,
        subreddit: str,
        query: str,
        **params
    ) -> List[Dict]:
        """Search posts in subreddit."""
        data = await self.get(
            f"{self.base_url}/r/{subreddit}/search",
            q=query,
            **params
        )
        return [self._parse_post(post) for post in data["data"]["children"]]
```

### 5. Integration Layer

```python
# In data_sources/integration.py
class DataSourceManager:
    """Manages multiple data sources."""
    def __init__(self):
        # Scientific
        self.chembl = ChEMBLClient()
        self.pubchem = PubChemClient()
        
        # Community
        self.psychonaut = PsychonautClient()
        self.erowid = ErowidClient()
        self.tripsit = TripSitClient()
        
        # Social
        self.reddit = RedditClient()
        self.twitter = TwitterClient()
        
    async def search_all(self, query: str) -> Dict[str, List]:
        """Search across all data sources."""
        results = {}
        
        # Scientific
        results["chembl"] = await self.chembl.search_compounds(query)
        results["pubchem"] = await self.pubchem.search_compounds(query)
        
        # Community
        results["psychonaut"] = await self.psychonaut.search_substances(query)
        results["erowid"] = await self.erowid.search_substances(query)
        
        # Social
        results["reddit"] = await self.reddit.search_all_sources(query)
        results["twitter"] = await self.twitter.search_tweets(query)
        
        return results
```

## Implementation Steps

### Day 1: Core Infrastructure
1. Create directory structure
2. Implement base client
3. Add caching
4. Add rate limiting

### Day 2: Scientific Sources
1. Implement ChEMBL client
2. Implement PubChem client
3. Add data parsing
4. Add validation

### Day 3: Community Sources
1. Implement PsychonautWiki client
2. Implement Erowid client
3. Implement TripSit client
4. Add data parsing

### Day 4: Social Sources
1. Implement Reddit client
2. Implement Twitter client
3. Add authentication
4. Add rate limiting

### Day 5: Integration
1. Implement manager
2. Add error handling
3. Add logging
4. Add monitoring

## Validation Steps

### 1. API Integration
- [ ] Test API connections
- [ ] Test rate limiting
- [ ] Test caching
- [ ] Test error handling

### 2. Data Quality
- [ ] Validate responses
- [ ] Check data types
- [ ] Handle missing data
- [ ] Handle errors

### 3. Performance
- [ ] Check response times
- [ ] Monitor rate limits
- [ ] Test concurrency
- [ ] Test recovery

## Success Criteria

### 1. Functionality
- All APIs accessible
- Data properly parsed
- Errors handled
- Rate limits respected

### 2. Performance
- Fast response times
- Efficient caching
- Resource management
- Error recovery

### 3. Integration
- Clean interfaces
- Type safety
- Good documentation
- Easy to use

## Next Steps

1. Set up infrastructure
2. Add scientific sources
3. Add community sources
4. Add social sources
5. Test integration
6. Document usage
