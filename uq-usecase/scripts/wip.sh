#!/bin/bash
# Stage and push work-in-progress changes from uq-usecase/ to fork.
# Usage: bash wip.sh "your commit message"
#        bash wip.sh          (uses default "wip" message)
#
# Before pushing, checks if origin has new commits. If they only touch .md files,
# rebases automatically. If they touch anything else, prints a warning and exits
# so you can review and pull manually.

set -e

GRIDKIT_DIR=~/gridkit
BRANCH=isatkaus/uq-usecase
MSG="${1:-wip}"

cd "$GRIDKIT_DIR"

# ── Check for remote changes ──────────────────────────────────────────────────
git fetch origin "$BRANCH" --quiet

LOCAL=$(git rev-parse HEAD)
REMOTE=$(git rev-parse "origin/$BRANCH")

if [[ "$LOCAL" != "$REMOTE" ]]; then
  # Find commits on remote not yet in local
  NON_MD=$(git diff --name-only "$LOCAL" "$REMOTE" | grep -v '\.md$' || true)
  if [[ -n "$NON_MD" ]]; then
    echo "ERROR: origin/$BRANCH has new commits that touch non-.md files:"
    echo "$NON_MD" | sed 's/^/  /'
    echo ""
    echo "Pull manually before running wip.sh:"
    echo "  git pull --rebase origin $BRANCH"
    exit 1
  fi
  echo "Remote has .md-only changes; rebasing before push..."
  git pull --rebase --autostash origin "$BRANCH" --quiet
fi

# ── Stage and commit ──────────────────────────────────────────────────────────
git add uq-usecase/
# Stage .gitignore only if it has changes
git diff --cached --quiet .gitignore 2>/dev/null || true
git diff --name-only .gitignore | grep -q . && git add .gitignore || true

git commit -m "$MSG"
git push

echo ""
echo "Pushed to origin/$BRANCH"
