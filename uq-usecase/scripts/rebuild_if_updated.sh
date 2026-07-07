#!/bin/bash
# Rebuild GridKit and MkDocs docs only if develop branch has new commits.
# Usage: bash rebuild_if_updated.sh [--force]
#rebuild_if_updated.sh # skip if no new commits
#rebuild_if_updated.sh --force  # rebuild regardless

set -e

GRIDKIT_DIR=~/gridkit
SCRIPTS_DIR=$GRIDKIT_DIR/uq-usecase/scripts
STAMP_FILE=$GRIDKIT_DIR/.last_built_commit

# ── Check for --force flag ────────────────────────────────────────────────────
FORCE=0
if [[ "${1}" == "--force" ]]; then
  FORCE=1
fi

# ── Fetch latest refs without changing working tree ───────────────────────────
echo "=== Checking for new commits on upstream/develop ==="
cd "$GRIDKIT_DIR"
git fetch upstream develop 2>&1

REMOTE_COMMIT=$(git rev-parse upstream/develop)
LAST_BUILT=$(cat "$STAMP_FILE" 2>/dev/null || echo "none")

echo "  Remote:     $REMOTE_COMMIT"
echo "  Last built: $LAST_BUILT"

if [[ $FORCE -eq 0 && "$REMOTE_COMMIT" == "$LAST_BUILT" ]]; then
  echo ""
  echo "No new commits since last build. Nothing to do."
  echo "(Use --force to rebuild anyway.)"
  exit 0
fi

echo ""
echo "New commits detected (or --force). Merging develop and rebuilding..."
echo ""

# ── Merge develop into personal branch (safe: uq-usecase/ not on develop) ─────
git merge upstream/develop

# ── Rebuild GridKit ───────────────────────────────────────────────────────────
echo ""
echo "=== [1/3] Building GridKit ==="
bash "$SCRIPTS_DIR/6_build_gridkit.sh"

# ── Rebuild MkDocs docs ───────────────────────────────────────────────────────
echo ""
echo "=== [2/3] Building MkDocs docs ==="
bash "$SCRIPTS_DIR/build_mkdocs.sh"

# ── Run tests ────────────────────────────────────────────────────────────────
echo ""
echo "=== [3/3] Running tests ==="
ctest --test-dir "$GRIDKIT_DIR/build" --output-on-failure

# ── Record the commit we just built ──────────────────────────────────────────
# Record the upstream/develop commit we just built against
echo "$REMOTE_COMMIT" > "$STAMP_FILE"

echo ""
echo "=== Done ==="
echo "  GridKit built:  upstream/develop @ $REMOTE_COMMIT"
echo "  MkDocs built:   $GRIDKIT_DIR/uq-usecase/mkdocs-site/"
echo "  Tests ran:      all passed (stamp updated)"
echo ""
echo "To rsync docs to local machine (e.g. Mac), run from your local machine:"
echo "  cd /Users/isatkaus/projects/scidac/gridkit/"
echo "  ./rsync_and_open_mkdocs.sh"
