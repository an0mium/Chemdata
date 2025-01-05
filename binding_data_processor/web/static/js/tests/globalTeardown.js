// Import required utilities
const fs = require('fs-extra');

async function teardown() {
    try {
        // Get test configuration
        const config = global.__TEST_CONFIG__;
        if (!config) {
            console.warn('No test configuration found during teardown');
            return;
        }

        // Stop test server
        if (config.server?.instance) {
            await new Promise((resolve) => {
                config.server.instance.close(() => {
                    console.log('Test server stopped');
                    resolve();
                });
            });
        }

        // Stop MongoDB instance
        if (config.mongod) {
            await config.mongod.stop();
            console.log('MongoDB memory server stopped');
        }

        // Clean up test directories
        await Promise.all([
            fs.remove(config.dirs.data),
            fs.remove(config.dirs.cache),
            fs.remove(config.dirs.exports),
            fs.remove(config.dirs.uploads)
        ]);
        console.log('Test directories cleaned');

        // Clean up environment variables
        delete process.env.NODE_ENV;
        delete process.env.MONGO_URI;
        delete process.env.TEST_DATA_DIR;
        delete process.env.TEST_CACHE_DIR;
        delete process.env.TEST_EXPORT_DIR;
        delete process.env.TEST_UPLOAD_DIR;

        // Clean up global configuration
        delete global.__TEST_CONFIG__;

        console.log('Test environment teardown complete');
    } catch (error) {
        console.error('Test environment teardown failed:', error);
        throw error;
    }
}

module.exports = teardown;
