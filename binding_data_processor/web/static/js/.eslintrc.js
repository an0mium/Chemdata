module.exports = {
    root: true,
    env: {
        browser: true,
        es2021: true,
        node: true,
        jest: true
    },
    extends: [
        'eslint:recommended',
        'plugin:jest/recommended',
        'plugin:jest/style',
        'plugin:import/errors',
        'plugin:import/warnings',
        'prettier'
    ],
    parserOptions: {
        ecmaVersion: 12,
        sourceType: 'module'
    },
    plugins: [
        'jest',
        'import',
        'prettier'
    ],
    settings: {
        'import/resolver': {
            node: {
                extensions: ['.js']
            }
        }
    },
    rules: {
        // Error Prevention
        'no-console': ['warn', { allow: ['warn', 'error'] }],
        'no-debugger': 'warn',
        'no-unused-vars': ['error', {
            argsIgnorePattern: '^_',
            varsIgnorePattern: '^_'
        }],
        'no-use-before-define': ['error', {
            functions: false,
            classes: true,
            variables: true
        }],

        // Best Practices
        'array-callback-return': 'error',
        'consistent-return': 'error',
        'default-case': 'error',
        'eqeqeq': ['error', 'always'],
        'no-eval': 'error',
        'no-implied-eval': 'error',
        'no-param-reassign': 'error',
        'no-return-await': 'error',
        'no-throw-literal': 'error',
        'prefer-promise-reject-errors': 'error',
        'require-await': 'error',

        // ES6+ Features
        'arrow-body-style': ['error', 'as-needed'],
        'no-var': 'error',
        'prefer-const': 'error',
        'prefer-destructuring': ['error', {
            array: true,
            object: true
        }],
        'prefer-rest-params': 'error',
        'prefer-spread': 'error',
        'prefer-template': 'error',

        // Import Rules
        'import/first': 'error',
        'import/no-duplicates': 'error',
        'import/no-unresolved': 'error',
        'import/no-cycle': 'error',
        'import/order': ['error', {
            groups: [
                'builtin',
                'external',
                'internal',
                'parent',
                'sibling',
                'index'
            ],
            'newlines-between': 'always',
            alphabetize: {
                order: 'asc',
                caseInsensitive: true
            }
        }],

        // Jest Rules
        'jest/consistent-test-it': ['error', {
            fn: 'test',
            withinDescribe: 'test'
        }],
        'jest/expect-expect': 'error',
        'jest/no-disabled-tests': 'warn',
        'jest/no-focused-tests': 'error',
        'jest/no-identical-title': 'error',
        'jest/prefer-to-be': 'error',
        'jest/prefer-to-contain': 'error',
        'jest/prefer-to-have-length': 'error',
        'jest/valid-expect': 'error',

        // Style Rules (not covered by Prettier)
        'no-nested-ternary': 'error',
        'no-unneeded-ternary': 'error',
        'spaced-comment': ['error', 'always'],
        'capitalized-comments': ['error', 'always', {
            ignorePattern: 'pragma|ignore',
            ignoreInlineComments: true
        }],

        // Documentation
        'valid-jsdoc': ['error', {
            requireReturn: false,
            requireReturnType: false,
            requireParamType: false,
            prefer: {
                returns: 'return'
            }
        }]
    },
    overrides: [
        {
            files: ['**/*.test.js', '**/__tests__/**/*.js'],
            env: {
                jest: true
            },
            rules: {
                'max-lines': 'off',
                'max-nested-callbacks': 'off',
                'max-statements': 'off'
            }
        }
    ]
};
