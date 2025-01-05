module.exports = {
    // Test Environment
    testEnvironment: 'jsdom',
    setupFilesAfterEnv: ['<rootDir>/tests/setup.js'],
    globalSetup: '<rootDir>/tests/globalSetup.js',
    globalTeardown: '<rootDir>/tests/globalTeardown.js',

    // Test File Patterns
    testMatch: [
        '**/tests/**/*.test.js',
        '**/__tests__/**/*.js',
        '**/tests/**/*.spec.js',
        '**/binding_data_processor/**/*.test.js'
    ],
    testPathIgnorePatterns: [
        '/node_modules/',
        '/dist/',
        '/coverage/',
        '/vendor/'
    ],

    // Module Configuration
    moduleFileExtensions: ['js', 'json'],
    moduleNameMapper: {
        // Static assets
        '\\.(css|less|sass|scss)$': '<rootDir>/tests/mocks/styleMock.js',
        '\\.(gif|ttf|eot|svg|png|jpg|jpeg)$': '<rootDir>/tests/mocks/fileMock.js',
        // Project aliases
        '^@/(.*)$': '<rootDir>/src/$1',
        '^@binding_data_processor/(.*)$': '<rootDir>/binding_data_processor/$1',
        '^@models/(.*)$': '<rootDir>/binding_data_processor/models/$1',
        '^@processors/(.*)$': '<rootDir>/binding_data_processor/processors/$1',
        '^@pipeline/(.*)$': '<rootDir>/binding_data_processor/pipeline/$1',
        '^@web/(.*)$': '<rootDir>/binding_data_processor/web/$1'
    },

    // Coverage Configuration
    collectCoverage: true,
    collectCoverageFrom: [
        'binding_data_processor/**/*.js',
        '!**/node_modules/**',
        '!**/vendor/**',
        '!**/tests/**',
        '!**/coverage/**',
        '!**/__tests__/**',
        '!**/dist/**',
        '!jest.config.js'
    ],
    coverageDirectory: 'coverage',
    coverageReporters: ['text', 'lcov', 'html'],
    coverageThreshold: {
        global: {
            branches: 80,
            functions: 80,
            lines: 80,
            statements: 80
        }
    },

    // Performance Configuration
    maxWorkers: '50%',
    timers: 'modern',
    testTimeout: 10000,

    // Mocking Configuration
    clearMocks: true,
    resetMocks: false,
    restoreMocks: true,

    // Reporting Configuration
    verbose: true,
    notify: true,
    notifyMode: 'failure-change',
    reporters: [
        'default',
        [
            'jest-junit',
            {
                outputDirectory: 'reports/junit',
                outputName: 'js-test-results.xml',
                classNameTemplate: '{classname}',
                titleTemplate: '{title}',
                ancestorSeparator: ' › ',
                usePathForSuiteName: true
            }
        ]
    ],

    // Watch Configuration
    watchPlugins: [
        'jest-watch-typeahead/filename',
        'jest-watch-typeahead/testname'
    ],

    // Error Handling
    bail: 1,
    errorOnDeprecated: true,

    // Transform Configuration
    transform: {
        '^.+\\.jsx?$': 'babel-jest'
    },
    transformIgnorePatterns: [
        '/node_modules/',
        '\\.pnp\\.[^\\/]+$'
    ],

    // Global Settings
    globals: {
        __DEV__: true,
        BINDING_DATA_PROCESSOR_VERSION: require('./package.json').version
    }
};
