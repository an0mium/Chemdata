"""Tests for enhanced PubMed client combining API and web scraping functionality."""

import pytest
from unittest.mock import AsyncMock, MagicMock, Mock, patch
from datetime import datetime
from pathlib import Path
from bs4 import BeautifulSoup

from ..clients.pubmed import EnhancedPubMedClient, PubMedArticle
from ...models.core import CompoundData


@pytest.fixture
def mock_api_client():
    """Mock PubMed API client."""
    with patch("...data_sources.pubmed.PubMedClient") as mock:
        client = MagicMock()
        client._make_request = AsyncMock()
        client.ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        client.EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        mock.return_value = client
        yield client


@pytest.fixture
def mock_crawler():
    """Mock AsyncWebCrawler with comprehensive test data."""
    with patch("crawl4ai.AsyncWebCrawler") as mock:
        crawler = AsyncMock()
        crawler.arun.return_value = MagicMock(
            success=True,
            extracted_content=[
                {
                    "title": ["Test Article"],
                    "abstract": ["Test abstract mentioning caffeine"],
                    "authors": [["Author 1", "Author 2"]],
                    "journal": ["Test Journal"],
                    "year": ["2024"],
                    "doi": ["10.1234/test"],
                    "pmid": ["12345678"],
                    "citations": ["42"],
                    "compounds": ["caffeine", "adenosine"],
                    "effects": ["improved memory", "increased focus"],
                    "mechanisms": ["adenosine receptor antagonism"],
                    "safety": ["well tolerated", "mild side effects"],
                    "methods": ["double-blind study"],
                    "results": ["significant improvement"],
                    "conclusions": ["effective nootropic"],
                }
            ],
            screenshot=Path("test.png"),
            javascript_logs=["Page loaded successfully"],
            chunking_stats={"chunks": 5, "overlap": 100},
        )
        mock.return_value.__aenter__.return_value = crawler
        mock.return_value.__aexit__.return_value = None
        yield crawler


@pytest.fixture
def mock_session():
    """Mock aiohttp ClientSession with realistic API responses."""
    with patch("aiohttp.ClientSession") as mock:
        session = AsyncMock()
        response = AsyncMock()
        response.ok = True
        response.json.return_value = {
            "esearchresult": {
                "idlist": ["12345678", "87654321"],
                "count": "2",
                "retmax": "2",
                "retstart": "0",
                "querytranslation": "caffeine[All Fields]",
            }
        }
        response.text.return_value = """<?xml version="1.0"?>
        <PubmedArticleSet>
            <PubmedArticle>
                <MedlineCitation Status="MEDLINE">
                    <PMID Version="1">12345678</PMID>
                    <Article>
                        <ArticleTitle>Test API Article</ArticleTitle>
                        <Abstract>
                            <AbstractText>Test API abstract studying caffeine effects.</AbstractText>
                        </Abstract>
                        <Journal>
                            <Title>Test API Journal</Title>
                            <JournalIssue>
                                <PubDate>
                                    <Year>2024</Year>
                                </PubDate>
                            </JournalIssue>
                        </Journal>
                        <AuthorList CompleteYN="Y">
                            <Author ValidYN="Y">
                                <LastName>Smith</LastName>
                                <ForeName>John</ForeName>
                                <Affiliation>Test University</Affiliation>
                            </Author>
                        </AuthorList>
                    </Article>
                    <ChemicalList>
                        <Chemical>
                            <NameOfSubstance UI="D002110">Caffeine</NameOfSubstance>
                            <RegistryNumber>58-08-2</RegistryNumber>
                        </Chemical>
                        <Chemical>
                            <NameOfSubstance UI="D000241">Adenosine</NameOfSubstance>
                            <RegistryNumber>58-61-7</RegistryNumber>
                        </Chemical>
                    </ChemicalList>
                    <MeshHeadingList>
                        <MeshHeading>
                            <DescriptorName UI="D002110">Caffeine</DescriptorName>
                            <QualifierName UI="PK">pharmacokinetics</QualifierName>
                        </MeshHeading>
                    </MeshHeadingList>
                </MedlineCitation>
                <PubmedData>
                    <ArticleIdList>
                        <ArticleId IdType="pubmed">12345678</ArticleId>
                        <ArticleId IdType="doi">10.1234/test</ArticleId>
                    </ArticleIdList>
                    <ReferenceList>
                        <Reference>
                            <Citation>Related Article 1</Citation>
                            <ArticleIdList>
                                <ArticleId IdType="pubmed">11111111</ArticleId>
                            </ArticleIdList>
                        </Reference>
                    </ReferenceList>
                </PubmedData>
            </PubmedArticle>
        </PubmedArticleSet>
        """
        session.get.return_value.__aenter__.return_value = response
        mock.return_value.__aenter__.return_value = session
        mock.return_value.__aexit__.return_value = None
        yield session


@pytest.fixture
def mock_llm_extractor():
    """Mock NootropicLLMExtractor with realistic extraction results."""
    extractor = MagicMock()
    extractor.extract_chemical_info = AsyncMock(
        return_value=[
            "caffeine",
            "adenosine",
            "theobromine",
        ]
    )
    extractor.extract_mechanisms = AsyncMock(
        return_value=[
            "adenosine receptor antagonism",
            "phosphodiesterase inhibition",
        ]
    )
    extractor.extract_safety = AsyncMock(
        return_value=[
            "well tolerated in healthy adults",
            "caution in sensitive individuals",
            "monitor cardiovascular effects",
        ]
    )
    return extractor


@pytest.fixture
def mock_validator():
    """Mock validator with configurable validation results."""
    with patch("..validation.enhanced.EnhancedValidator") as mock:
        validator = AsyncMock()
        validator.validate.return_value = MagicMock(
            is_valid=True,
            confidence=0.95,
            validation_details={
                "schema_valid": True,
                "cross_source_valid": True,
                "temporal_valid": True,
                "community_valid": True,
            },
        )
        mock.return_value = validator
        yield validator


@pytest.fixture
def client(mock_crawler, mock_session, mock_validator, mock_llm_extractor):
    """Create test client with all mocked dependencies."""
    return EnhancedPubMedClient(
        api_key="test_key",
        validator=mock_validator,
        llm_provider="test_llm",
        api_token="test_token",
        cache_dir=None,
        proxy_enabled=True,
        proxy_rotation=True,
        proxy_retry_count=3,
        rate_limit=10,
        rate_limit_delay=60,
        logger=MagicMock(),
    )


@pytest.fixture
def test_compound():
    """Create test compound with comprehensive metadata."""
    return CompoundData(
        name="Caffeine",
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        cas_number="58-08-2",
        iupac_name="1,3,7-trimethylpurine-2,6-dione",
        common_name_1="caffeine",
        common_name_2="trimethylxanthine",
        common_name_3="N/A",
        molecular_weight=194.19,
        molecular_formula="C8H10N4O2",
        known_effects=["stimulant", "nootropic"],
        known_mechanisms=["adenosine antagonist"],
        safety_class="GRAS",
    )


@pytest.mark.asyncio
async def test_search_articles_success(client, mock_crawler, mock_session):
    """Test successful article search with both API and web data."""
    # Search articles
    articles = await client.search_articles(
        query="caffeine",
        max_results=10,
        from_year=2020,
        to_year=2024,
        sort="relevance",
    )

    # Verify results
    assert len(articles) > 0
    article = articles[0]
    assert isinstance(article, PubMedArticle)
    assert article.title == "Test Article"
    assert "caffeine" in article.abstract.lower()
    assert len(article.authors) == 2
    assert article.journal == "Test Journal"
    assert article.year == 2024
    assert article.doi == "10.1234/test"
    assert article.pmid == "12345678"
    assert article.citations == 42
    assert len(article.compounds) == 2
    assert len(article.effects) == 2
    assert len(article.mechanisms) == 1
    assert len(article.safety_notes) == 2
    assert article.confidence > 0.5
    assert article.metadata["source"] == "pubmed"
    assert "extraction_date" in article.metadata
    assert "screenshot" in article.metadata
    assert "javascript_logs" in article.metadata
    assert "chunking_stats" in article.metadata


@pytest.mark.asyncio
async def test_search_articles_web_only(client, mock_crawler, mock_session):
    """Test article search with only web data when API fails."""
    # Mock API error
    mock_session.get.return_value.__aenter__.return_value.ok = False

    # Search articles
    articles = await client.search_articles("caffeine")

    # Verify results
    assert len(articles) > 0
    article = articles[0]
    assert article.title == "Test Article"
    assert article.confidence < 0.8  # Lower confidence with only web data
    assert "web_only" in article.metadata


@pytest.mark.asyncio
async def test_search_articles_api_only(client, mock_crawler, mock_session):
    """Test article search with only API data when web scraping fails."""
    # Mock web scraping failure
    mock_crawler.arun.return_value.success = False

    # Search articles
    articles = await client.search_articles("caffeine")

    # Verify results
    assert len(articles) > 0
    article = articles[0]
    assert article.title == "Test API Article"
    assert article.confidence < 0.8  # Lower confidence with only API data
    assert "api_only" in article.metadata


@pytest.mark.asyncio
async def test_search_articles_partial_data(client, mock_crawler, mock_session):
    """Test article search with partial data from both sources."""
    # Mock crawler response with missing fields
    mock_crawler.arun.return_value = MagicMock(
        success=True,
        extracted_content=[
            {
                "title": ["Test Title"],
                "abstract": [],
                "authors": [],
                "journal": [],
                "year": [],
                "citations": [],
                "pmid": [],
                "doi": [],
                "compounds": [],
                "effects": [],
                "mechanisms": [],
                "safety": [],
            }
        ],
        screenshot=Path("test.png"),
        javascript_logs=[],
    )

    # Mock API response with missing fields
    mock_session.get.return_value.__aenter__.return_value.text.return_value = """<?xml version="1.0"?>
    <PubmedArticleSet>
        <PubmedArticle>
            <MedlineCitation>
                <Article>
                    <ArticleTitle>Test API Title</ArticleTitle>
                </Article>
            </MedlineCitation>
        </PubmedArticle>
    </PubmedArticleSet>
    """

    # Search articles
    articles = await client.search_articles("test")

    # Verify results
    assert len(articles) == 1
    article = articles[0]
    assert article.title == "Test Title"  # Web data takes precedence
    assert article.abstract is None
    assert article.authors == []
    assert article.journal is None
    assert article.year is None
    assert article.citations is None
    assert article.pmid is None
    assert article.doi is None
    assert article.compounds == []
    assert article.effects == []
    assert article.mechanisms == []
    assert article.safety_notes == []
    assert article.confidence < 0.5  # Lower confidence with missing fields
    assert "partial_data" in article.metadata


@pytest.mark.asyncio
async def test_search_articles_caching(client, mock_crawler, mock_session):
    """Test article search caching behavior."""
    # First search
    articles1 = await client.search_articles("caffeine")
    assert len(articles1) > 0

    # Mock crawler and API failure
    mock_crawler.arun.side_effect = Exception("Crawl error")
    mock_session.get.return_value.__aenter__.return_value.ok = False

    # Second search should use cache
    articles2 = await client.search_articles("caffeine")
    assert len(articles2) > 0
    assert articles2[0].title == articles1[0].title
    assert "cache_hit" in articles2[0].metadata


@pytest.mark.asyncio
async def test_get_article_success(client, mock_crawler, mock_session):
    """Test successful single article retrieval."""
    # Get article
    article = await client.get_article("12345678")

    # Verify results
    assert article is not None
    assert isinstance(article, PubMedArticle)
    assert article.pmid == "12345678"
    assert article.title == "Test Article"
    assert "caffeine" in article.abstract.lower()
    assert len(article.authors) == 2
    assert article.journal == "Test Journal"
    assert article.year == 2024
    assert article.doi == "10.1234/test"
    assert article.citations == 42
    assert len(article.compounds) == 2
    assert len(article.effects) == 2
    assert len(article.mechanisms) == 1
    assert len(article.safety_notes) == 2
    assert article.confidence > 0.5
    assert article.metadata["source"] == "pubmed"
    assert "extraction_date" in article.metadata


@pytest.mark.asyncio
async def test_get_article_not_found(client, mock_crawler, mock_session):
    """Test article not found handling."""
    # Mock crawler and API failure
    mock_crawler.arun.return_value.success = False
    mock_session.get.return_value.__aenter__.return_value.ok = False

    # Get article
    article = await client.get_article("invalid")

    # Verify results
    assert article is None


@pytest.mark.asyncio
async def test_get_compound_articles_success(client, mock_crawler, mock_session, test_compound):
    """Test successful compound article search."""
    # Get articles
    articles = await client.get_compound_articles(
        compound=test_compound,
        max_results=10,
    )

    # Verify results
    assert len(articles) > 0
    article = articles[0]
    assert isinstance(article, PubMedArticle)
    assert article.title == "Test Article"
    assert "caffeine" in article.compounds
    assert len(article.effects) > 0
    assert len(article.mechanisms) > 0
    assert article.confidence > 0.5
    assert "compound_search" in article.metadata


@pytest.mark.asyncio
async def test_get_compound_articles_caching(client, mock_crawler, mock_session, test_compound):
    """Test compound article search caching."""
    # First search
    articles1 = await client.get_compound_articles(test_compound)
    assert len(articles1) > 0

    # Mock crawler and API failure
    mock_crawler.arun.side_effect = Exception("Crawl error")
    mock_session.get.return_value.__aenter__.return_value.ok = False

    # Second search should use cache
    articles2 = await client.get_compound_articles(test_compound)
    assert len(articles2) > 0
    assert articles2[0].title == articles1[0].title
    assert "cache_hit" in articles2[0].metadata


@pytest.mark.asyncio
async def test_error_handling(client, mock_crawler, mock_session):
    """Test comprehensive error handling."""
    # Test network errors
    mock_crawler.arun.side_effect = Exception("Network error")
    mock_session.get.side_effect = Exception("API error")

    articles = await client.search_articles("caffeine")
    assert articles == []

    article = await client.get_article("12345678")
    assert article is None

    # Test parsing errors
    mock_crawler.arun.side_effect = None
    mock_session.get.side_effect = None
    mock_session.get.return_value.__aenter__.return_value.text.return_value = "Invalid XML"

    articles = await client.search_articles("caffeine")
    assert len(articles) > 0  # Should still get web results
    assert "parse_error" in articles[0].metadata

    # Test validation errors
    mock_crawler.arun.return_value.success = True
    mock_session.get.return_value.__aenter__.return_value.ok = True
    client.validator.validate.return_value = MagicMock(is_valid=False)

    articles = await client.search_articles("caffeine")
    assert articles == []


@pytest.mark.asyncio
async def test_validation(client, mock_crawler, mock_session, mock_validator):
    """Test data validation scenarios."""
    # Test schema validation
    mock_validator.validate.return_value = MagicMock(is_valid=False, validation_details={"schema_valid": False})
    articles = await client.search_articles("caffeine")
    assert len(articles) == 0

    # Test cross-source validation
    mock_validator.validate.return_value = MagicMock(is_valid=False, validation_details={"cross_source_valid": False})
    article = await client.get_article("12345678")
    assert article is None

    # Test temporal validation
    mock_validator.validate.return_value = MagicMock(is_valid=False, validation_details={"temporal_valid": False})
    articles = await client.get_compound_articles(test_compound())
    assert len(articles) == 0


def test_build_compound_query(client, test_compound):
    """Test compound query building with various inputs."""
    query = client._build_compound_query(test_compound)
    assert '"Caffeine"' in query
    assert '"58-08-2"' in query
    assert '"1,3,7-trimethylpurine-2,6-dione"' in query
    assert '"caffeine"' in query
    assert '"trimethylxanthine"' in query
    assert "N/A" not in query
    assert "C8H10N4O2" in query
    assert "stimulant" in query.lower()
    assert "nootropic" in query.lower()


@pytest.mark.asyncio
async def test_relevance_check(client, test_compound):
    """Test article relevance checking."""
    # Create test article
    article = PubMedArticle(
        pmid="12345678",
        title="Test Article",
        abstract="Study of caffeine effects",
        authors=["Author 1"],
        journal="Test Journal",
        year=2024,
        compounds=["caffeine"],
        effects=["improved memory"],
        mechanisms=["adenosine antagonism"],
        binding_data=[{"compound": "caffeine", "registry_number": "58-08-2"}],
    )

    # Check relevance
    assert await client._is_relevant(article, test_compound)

    # Test irrelevant article
    irrelevant = PubMedArticle(
        pmid="87654321",
        title="Unrelated Article",
        abstract="Study of glucose metabolism",
        authors=["Author 2"],
        journal="Test Journal",
        year=2024,
        compounds=[],
        effects=[],
        mechanisms=[],
        binding_data=[],
    )
    assert not await client._is_relevant(irrelevant, test_compound)


def test_calculate_confidence(client):
    """Test confidence scoring with various data combinations."""
    # Full data from both web and API
    web_data = {
        "title": "Title",
        "abstract": "Abstract",
        "authors": ["Author"],
        "journal": "Journal",
        "year": 2024,
        "doi": "10.1234/test",
        "compounds": ["compound"],
        "effects": ["effect"],
        "mechanisms": ["mechanism"],
        "safety": ["safe"],
        "methods": ["method"],
        "results": ["result"],
        "conclusions": ["conclusion"],
    }
    api_data = {
        "title": "Title",
        "abstract": "Abstract",
        "authors": ["Author"],
        "journal": "Journal",
        "year": 2024,
        "doi": "10.1234/test",
        "binding_data": [{"compound": "test", "registry_number": "123"}],
        "mesh_terms": ["Term 1", "Term 2"],
    }
    assert client._calculate_confidence(web_data, api_data) == 1.0

    # Partial data from web only
    web_data = {
        "title": "Title",
        "abstract": "Abstract",
        "authors": ["Author"],
    }
    api_data = {}
    assert 0.2 < client._calculate_confidence(web_data, api_data) < 0.5

    # Partial data from API only
    web_data = {}
    api_data = {
        "title": "Title",
        "abstract": "Abstract",
        "authors": ["Author"],
    }
    assert 0.2 < client._calculate_confidence(web_data, api_data) < 0.5

    # Conflicting data
    web_data = {
        "title": "Title 1",
        "abstract": "Abstract 1",
    }
    api_data = {
        "title": "Title 2",
        "abstract": "Abstract 2",
    }
    assert client._calculate_confidence(web_data, api_data) < 0.3

    # No data
    assert client._calculate_confidence({}, {}) == 0.0


def test_extract_article_metadata():
    """Test XML metadata extraction with various input formats."""
    xml = """<?xml version="1.0"?>
    <PubmedArticle>
        <MedlineCitation Status="MEDLINE">
            <PMID Version="1">12345678</PMID>
            <Article>
                <ArticleTitle>Test Title</ArticleTitle>
                <Abstract>
                    <AbstractText>Test Abstract</AbstractText>
                </Abstract>
                <Journal>
                    <Title>Test Journal</Title>
                    <JournalIssue>
                        <PubDate>
                            <Year>2024</Year>
                        </PubDate>
                    </JournalIssue>
                </Journal>
                <AuthorList CompleteYN="Y">
                    <Author ValidYN="Y">
                        <LastName>Smith</LastName>
                        <ForeName>John</ForeName>
                        <Affiliation>Test University</Affiliation>
                    </Author>
                </AuthorList>
            </Article>
            <ChemicalList>
                <Chemical>
                    <NameOfSubstance UI="D002110">Caffeine</NameOfSubstance>
                    <RegistryNumber>58-08-2</RegistryNumber>
                </Chemical>
                <Chemical>
                    <NameOfSubstance UI="D000241">Adenosine</NameOfSubstance>
                    <RegistryNumber>58-61-7</RegistryNumber>
                </Chemical>
            </ChemicalList>
            <MeshHeadingList>
                <MeshHeading>
                    <DescriptorName UI="D002110">Caffeine</DescriptorName>
                    <QualifierName UI="PK">pharmacokinetics</QualifierName>
                </MeshHeading>
            </MeshHeadingList>
        </MedlineCitation>
        <PubmedData>
            <ArticleIdList>
                <ArticleId IdType="pubmed">12345678</ArticleId>
                <ArticleId IdType="doi">10.1234/test</ArticleId>
            </ArticleIdList>
            <ReferenceList>
                <Reference>
                    <Citation>Related Article 1</Citation>
                    <ArticleIdList>
                        <ArticleId IdType="pubmed">11111111</ArticleId>
                    </ArticleIdList>
                </Reference>
            </ReferenceList>
        </PubmedData>
    </PubmedArticle>
    """

    client = EnhancedPubMedClient()
    soup = BeautifulSoup(xml, "xml")
    metadata = client._extract_article_metadata(soup.find("PubmedArticle"))

    assert metadata["title"] == "Test Title"
    assert metadata["abstract"] == "Test Abstract"
    assert metadata["journal"] == "Test Journal"
    assert metadata["authors"] == ["John Smith"]
    assert metadata["year"] == 2024
    assert metadata["doi"] == "10.1234/test"
    assert metadata["chemicals"] == ["Caffeine", "Adenosine"]
    assert metadata["binding_data"] == [
        {"compound": "Caffeine", "registry_number": "58-08-2"},
        {"compound": "Adenosine", "registry_number": "58-61-7"},
    ]
    assert metadata["mesh_terms"] == ["Caffeine pharmacokinetics"]
    assert metadata["references"] == ["11111111"]


def test_article_schema():
    """Test PubMedArticle schema with comprehensive metadata."""
    article = PubMedArticle(
        pmid="12345678",
        title="Test Title",
        abstract="Test Abstract",
        authors=["Test Author"],
        journal="Test Journal",
        year=2024,
        doi="10.1234/test",
        citations=42,
        compounds=["Test Compound"],
        effects=["Test Effect"],
        mechanisms=["Test Mechanism"],
        safety_notes=["Test Safety"],
        methods=["Test Method"],
        results=["Test Result"],
        conclusions=["Test Conclusion"],
        binding_data=[{"compound": "test", "registry_number": "123"}],
        confidence=0.9,
        metadata={
            "source": "test",
            "extraction_date": datetime.now().isoformat(),
            "screenshot": "test.png",
            "javascript_logs": ["test log"],
            "chunking_stats": {"chunks": 5},
            "mesh_terms": ["Term 1"],
            "references": ["11111111"],
        },
    )

    assert article.pmid == "12345678"
    assert article.title == "Test Title"
    assert article.abstract == "Test Abstract"
    assert article.authors == ["Test Author"]
    assert article.journal == "Test Journal"
    assert article.year == 2024
    assert article.doi == "10.1234/test"
    assert article.citations == 42
    assert article.compounds == ["Test Compound"]
    assert article.effects == ["Test Effect"]
    assert article.mechanisms == ["Test Mechanism"]
    assert article.safety_notes == ["Test Safety"]
    assert article.methods == ["Test Method"]
    assert article.results == ["Test Result"]
    assert article.conclusions == ["Test Conclusion"]
    assert article.binding_data == [{"compound": "test", "registry_number": "123"}]
    assert article.confidence == 0.9
    assert article.metadata["source"] == "test"
    assert "extraction_date" in article.metadata
    assert "screenshot" in article.metadata
    assert "javascript_logs" in article.metadata
    assert "chunking_stats" in article.metadata
    assert "mesh_terms" in article.metadata
    assert "references" in article.metadata
