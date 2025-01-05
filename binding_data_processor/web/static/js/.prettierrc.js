module.exports = {
    // Basic Formatting
    printWidth: 100,
    tabWidth: 4,
    useTabs: false,
    semi: true,
    singleQuote: true,
    quoteProps: 'as-needed',

    // JSX Formatting
    jsxSingleQuote: false,
    jsxBracketSameLine: false,

    // Punctuation
    trailingComma: 'es5',
    bracketSpacing: true,

    // Special Formatting
    arrowParens: 'always',
    endOfLine: 'lf',

    // Wrapping
    proseWrap: 'preserve',
    htmlWhitespaceSensitivity: 'css',

    // File Patterns
    overrides: [
        {
            files: '*.{json,yml,yaml,md}',
            options: {
                tabWidth: 2
            }
        },
        {
            files: '*.md',
            options: {
                proseWrap: 'always'
            }
        }
    ],

    // Plugin Options
    plugins: [
        'prettier-plugin-organize-imports'
    ],

    // Import Organization
    importOrder: [
        '^@binding_data_processor/(.*)$',
        '^@testing-library/(.*)$',
        '^[a-zA-Z](.*)$',
        '^[./]'
    ],
    importOrderSeparation: true,
    importOrderSortSpecifiers: true,

    // Editor Integration
    requirePragma: false,
    insertPragma: false,

    // Error Handling
    rangeStart: 0,
    rangeEnd: Infinity,

    // Parser Selection
    parser: 'babel'
};
