module.exports = {
    // JavaScript/TypeScript files
    '**/*.{js,jsx,ts,tsx}': [
        'prettier --write',
        'eslint --fix',
        'jest --bail --findRelatedTests',
    ],

    // Style files
    '**/*.{css,scss,sass}': [
        'prettier --write',
        'stylelint --fix',
    ],

    // JSON, YAML, Markdown files
    '**/*.{json,yml,yaml,md}': [
        'prettier --write',
    ],

    // HTML files
    '**/*.html': [
        'prettier --write',
        'htmlhint',
    ],

    // Configuration files
    '.{prettierrc,eslintrc,stylelintrc,babelrc}': [
        'prettier --write',
    ],

    // Package files
    'package.json': [
        'prettier --write',
        'sort-package-json',
    ],

    // All files
    '*': [
        // Prevent large files from being committed
        (files) => files.map((file) => {
            const maxSize = 500 * 1024; // 500KB
            const stats = require('fs').statSync(file);
            if (stats.size > maxSize) {
                throw new Error(
                    `File ${file} is larger than ${maxSize / 1024}KB. ` +
                    'Please check if this file should be committed.'
                );
            }
            return null;
        }),
        // Prevent binary files from being committed
        (files) => files.map((file) => {
            const isBinary = require('is-binary-file').isBinaryFileSync;
            if (isBinary(file)) {
                throw new Error(
                    `File ${file} appears to be binary. ` +
                    'Please check if this file should be committed.'
                );
            }
            return null;
        }),
    ],
};
