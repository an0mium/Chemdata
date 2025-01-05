#!/bin/sh

# Install Husky
npx husky install

# Add pre-commit hook
npx husky add .husky/pre-commit "npx lint-staged"

# Add pre-push hook
npx husky add .husky/pre-push "npm run test && npm run type-check && npm run build"

# Make hooks executable
chmod +x .husky/pre-commit
chmod +x .husky/pre-push

# Install dependencies
npm install --save-dev \
  husky \
  lint-staged \
  prettier \
  eslint \
  stylelint \
  htmlhint \
  sort-package-json \
  is-binary-file

echo "Husky hooks installed successfully!"
