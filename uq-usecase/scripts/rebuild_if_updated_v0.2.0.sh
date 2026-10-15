#!/bin/bash
# Rebuild GridKit and MkDocs docs only if v0.2.0 tag has new commits.
# STABLE VERSION: This script pins to GridKit release v0.2.0 and does NOT advance
# beyond it, preventing incompatibilities from unreleased development work.
#
# Usage: bash rebuild_if_updated_v0.2.0.sh [--force]
#        bash rebuild_if_updated_v0.2.0.sh          # skip if v0.2.0 unchanged
#        bash rebuild_if_updated_v0.2.0.sh --force  # rebuild regardless

set -e

GRIDKIT_DIR=~/gridkit
SCRIPTS_DIR=$GRIDKIT_DIR/uq-usecase/scripts
STAMP_FILE=$GRIDKIT_DIR/.last_built_commit_v0.2.0
PINNED_VERSION="v0.2.0"

# ── Check for --force flag ────────────────────────────────────────────────────
FORCE=0
if [[ "${1}" == "--force" ]]; then
  FORCE=1
fi

# ── Fetch latest v0.2.0 tag (does not change working tree) ──────────────────
echo "=== Checking for updates to tag $PINNED_VERSION ==="
cd "$GRIDKIT_DIR"
git fetch origin tag "$PINNED_VERSION" 2>&1

REMOTE_COMMIT=$(git rev-parse "$PINNED_VERSION")
LAST_BUILT=$(cat "$STAMP_FILE" 2>/dev/null || echo "none")

echo "  Remote $PINNED_VERSION:  $REMOTE_COMMIT"
echo "  Last built:          $LAST_BUILT"

if [[ $FORCE -eq 0 && "$REMOTE_COMMIT" == "$LAST_BUILT" ]]; then
  echo ""
  echo "No updates to $PINNED_VERSION since last build. Nothing to do."
  echo "(Use --force to rebuild anyway.)"
  exit 0
fi

echo ""
echo "New version or --force detected. Rebasing onto $PINNED_VERSION and rebuilding..."
echo ""

# ── Rebase onto v0.2.0 tag (safe: uq-usecase/ is not on v0.2.0) ────────────
git checkout isatkaus/uq-usecase
git rebase "$PINNED_VERSION"

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

# ── Record the commit we just built against ──────────────────────────────────
echo "$REMOTE_COMMIT" > "$STAMP_FILE"

# ── Push rebased state to personal branch ────────────────────────────────────
echo ""
echo "=== Pushing rebased state to origin/isatkaus/uq-usecase ==="
git push -f origin isatkaus/uq-usecase

echo ""
echo "=== Done (version-pinned to $PINNED_VERSION) ==="
echo "  Pinned to:     $PINNED_VERSION @ $REMOTE_COMMIT"
echo "  Commit hash:   $REMOTE_COMMIT"
echo "  GridKit built: $GRIDKIT_DIR/build"
echo "  MkDocs built:  $GRIDKIT_DIR/uq-usecase/mkdocs-site/"
echo "  Tests:         all passed (stamp updated)"
echo ""
echo "To rsync docs to local machine (e.g. Mac), run from your local machine:"
echo "  rsync -avz --delete kestrel:~/gridkit/uq-usecase/mkdocs-site/ ~/projects/scidac/gridkit/gridkit-mkdocs/ && open ~/projects/scidac/gridkit/gridkit-mkdocs/index.html"
