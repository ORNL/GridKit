#!/bin/bash
# Check and display GridKit version status.
# Shows current local branch, pinned version, and distance from both develop and v0.2.0.
# Usage: bash check_gridkit_version.sh

set -e

GRIDKIT_DIR=~/gridkit
cd "$GRIDKIT_DIR"

echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║         GridKit Version & Branch Status                           ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""

# Current branch
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
CURRENT_COMMIT=$(git rev-parse --short HEAD)
echo "📍 Local Branch:  $CURRENT_BRANCH @ $CURRENT_COMMIT"

# Fetch latest tags and branches
git fetch origin tag v0.2.0 2>/dev/null || true
git fetch origin develop 2>/dev/null || true

echo ""
echo "────────────────────────────────────────────────────────────────────"
echo ""

# v0.2.0 status
V020_COMMIT=$(git rev-parse --short v0.2.0)
V020_DATE=$(git log -1 --format=%ai v0.2.0 | awk '{print $1}')
COMMITS_AHEAD_V020=$(git log --oneline v0.2.0..HEAD 2>/dev/null | wc -l)

echo "🏷️  v0.2.0 Release"
echo "   Commit:       $V020_COMMIT"
echo "   Date:         $V020_DATE"
echo "   Your branch:  $COMMITS_AHEAD_V020 commits ahead"

if [[ $COMMITS_AHEAD_V020 -eq 0 ]]; then
  echo "   Status:       ✅ AT PINNED VERSION"
else
  echo "   Status:       ⚠️  Diverged (use rebuild_if_updated_v0.2.0.sh to resync)"
fi

echo ""
echo "────────────────────────────────────────────────────────────────────"
echo ""

# develop status
DEVELOP_COMMIT=$(git rev-parse --short origin/develop)
DEVELOP_DATE=$(git log -1 --format=%ai origin/develop | awk '{print $1}')
COMMITS_AHEAD_DEVELOP=$(git log --oneline origin/develop..HEAD 2>/dev/null | wc -l)
COMMITS_BEHIND_DEVELOP=$(git log --oneline HEAD..origin/develop 2>/dev/null | wc -l)

echo "🌳 Origin/Develop (Unreleased)"
echo "   Commit:       $DEVELOP_COMMIT"
echo "   Date:         $DEVELOP_DATE"
if [[ $COMMITS_AHEAD_DEVELOP -gt 0 ]]; then
  echo "   Your branch:  $COMMITS_AHEAD_DEVELOP commits ahead"
elif [[ $COMMITS_BEHIND_DEVELOP -gt 0 ]]; then
  echo "   Your branch:  $COMMITS_BEHIND_DEVELOP commits behind"
else
  echo "   Your branch:  in sync"
fi
echo "   Status:       🚫 NOT PULLED (staying at v0.2.0)"

echo ""
echo "────────────────────────────────────────────────────────────────────"
echo ""
echo "ℹ️  Workflow:"
echo "   • Save work:               bash wip.sh 'your message'"
echo "   • Rebuild (if v0.2.0 moves): bash rebuild_if_updated_v0.2.0.sh"
echo "   • Force rebuild:           bash rebuild_if_updated_v0.2.0.sh --force"
echo ""
