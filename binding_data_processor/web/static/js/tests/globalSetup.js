// Import required utilities
const path = require('path');
const fs = require('fs-extra');
const { MongoMemoryServer } = require('mongodb-memory-server');
const { startTestServer } = require('../server/testUtils');

// Define test environment paths
const TEST_DATA_DIR = path.join(__dirname, '../../../test_data');
const TEST_CACHE_DIR = path.join(__dirname, '../../../test_cache');
const TEST_EXPORT_DIR = path.join(__dirname, '../../../test_exports');
const TEST_UPLOAD_DIR = path.join(__dirname, '../../../test_uploads');

// Test database configuration
const TEST_DB_CONFIG = {
    compounds: 'test_compounds',
    predictions: 'test_predictions',
    analysis: 'test_analysis',
    web_data: 'test_web_data'
};

// Test server configuration
const TEST_SERVER_CONFIG = {
    port: 3001,
    host: 'localhost',
    apiPrefix: '/api/v1'
};

async function setup() {
    try {
        // Create test directories
        await Promise.all([
            fs.ensureDir(TEST_DATA_DIR),
            fs.ensureDir(TEST_CACHE_DIR),
            fs.ensureDir(TEST_EXPORT_DIR),
            fs.ensureDir(TEST_UPLOAD_DIR)
        ]);

        // Copy test fixtures
        await fs.copy(
            path.join(__dirname, 'fixtures'),
            TEST_DATA_DIR,
            { overwrite: true }
        );

        // Start in-memory MongoDB
        const mongod = await MongoMemoryServer.create();
        const mongoUri = mongod.getUri();
        process.env.MONGO_URI = mongoUri;

        // Initialize test collections
        const { MongoClient } = require('mongodb');
        const client = new MongoClient(mongoUri);
        await client.connect();
        const db = client.db('test_db');

        // Create collections with indexes
        await Promise.all([
            db.createCollection(TEST_DB_CONFIG.compounds),
            db.createCollection(TEST_DB_CONFIG.predictions),
            db.createCollection(TEST_DB_CONFIG.analysis),
            db.createCollection(TEST_DB_CONFIG.web_data)
        ]);

        await Promise.all([
            db.collection(TEST_DB_CONFIG.compounds).createIndex({ cas_number: 1 }, { unique: true }),
            db.collection(TEST_DB_CONFIG.compounds).createIndex({ smiles: 1 }),
            db.collection(TEST_DB_CONFIG.predictions).createIndex({ compound_id: 1 }),
            db.collection(TEST_DB_CONFIG.analysis).createIndex({ compound_id: 1 }),
            db.collection(TEST_DB_CONFIG.web_data).createIndex({ compound_id: 1 })
        ]);

        await client.close();

        // Start test server
        const server = await startTestServer(TEST_SERVER_CONFIG);

        // Store configuration in global
        global.__TEST_CONFIG__ = {
            dirs: {
                data: TEST_DATA_DIR,
                cache: TEST_CACHE_DIR,
                exports: TEST_EXPORT_DIR,
                uploads: TEST_UPLOAD_DIR
            },
            db: {
                uri: mongoUri,
                config: TEST_DB_CONFIG
            },
            server: {
                ...TEST_SERVER_CONFIG,
                instance: server
            },
            mongod
        };

        // Set environment variables
        process.env.NODE_ENV = 'test';
        process.env.TEST_DATA_DIR = TEST_DATA_DIR;
        process.env.TEST_CACHE_DIR = TEST_CACHE_DIR;
        process.env.TEST_EXPORT_DIR = TEST_EXPORT_DIR;
        process.env.TEST_UPLOAD_DIR = TEST_UPLOAD_DIR;

        console.log('Test environment setup complete');
    } catch (error) {
        console.error('Test environment setup failed:', error);
        throw error;
    }
}

module.exports = setup;
