#!/bin/bash
# Stage and push work-in-progress changes from uq-usecase/ to isatkaus/uq-usecase branch.
# Mac version - pushes directly to ORNL/GridKit
# Usage: bash wip_mac.sh "your commit message"
#        bash wip_mac.sh          (uses default "wip" message)

set -e

GRIDKIT_DIR=/Users/isatkaus/projects/scidac/repos/gridkit
BRANCH=isatkaus/uq-usecase
MSG="${1:-wip}"

cd "$GRIDKIT_DIR"

# ── Stage and commit ──────────────────────────────────────────────────────────
git add uq-usecase/
# Stage .gitignore only if it has changes
git diff --name-only .gitignore | grep -q . && git add .gitignore || true

git commit -m "$MSG"
git push origin "$BRANCH"

echo ""
echo "✓ Pushed to origin/$BRANCH"
