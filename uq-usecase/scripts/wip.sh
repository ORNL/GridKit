#!/bin/bash
# Stage and push work-in-progress changes from uq-usecase/ to fork.
# Usage: bash wip.sh "your commit message"
#        bash wip.sh          (uses default "wip" message)

set -e

GRIDKIT_DIR=~/gridkit
MSG="${1:-wip}"

cd "$GRIDKIT_DIR"

git add uq-usecase/
# Stage .gitignore only if it has changes
git diff --cached --quiet .gitignore 2>/dev/null || true
git diff --name-only .gitignore | grep -q . && git add .gitignore || true

git commit -m "$MSG"
git push

echo ""
echo "Pushed to origin/isatkaus/uq-usecase"
