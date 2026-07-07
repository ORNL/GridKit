#!/bin/bash
# Build GridKit MkDocs documentation.
#
# Requirements:
#   conda env f-py312-basic (Python 3.12), or any Python 3.8+ env.
#   3 packages auto-installed on first run (idempotent):
#     mkdocs  mkdocs-material  mkdocs-awesome-pages-plugin
#   hooks.py uses only stdlib (re).
#   Module: conda  (loaded by this script)
#
# $GRIDKIT_DIR/ layout:
#   mkdocs/            ← MKDOCS_DIR
#     mkdocs.yml       ← generated each run
#     hooks.py         ← generated each run
#     docs/            ← docs_dir
#       index.md       -> GridKit/README.md
#       GridKit/
#         GridKit/     -> GridKit/GridKit/   (model READMEs)
#         examples/    -> GridKit/examples/
#         docs/        -> GridKit/docs/      (figures)
#         kestrel-usecase/ -> sh/            (kestrel_install.md)
#         *.md         -> GridKit/*.md
#       javascripts/   ← generated each run
#   mkdocs-site/       ← HTML output (rsync to local, e.g. Mac)
#
# GridKit/build/ is excluded — never symlinked (large, no .md content).
# To include more repo dirs, add a symlink in the setup block below.

echo "=== Building GridKit MkDocs ==="

GRIDKIT_DIR=~/gridkit
GRIDKIT_SRC=$GRIDKIT_DIR
MKDOCS_DIR=$GRIDKIT_DIR/uq-usecase/mkdocs

module load conda
conda activate f-py312-basic

# Install mkdocs packages if not present
pip install -q mkdocs mkdocs-material mkdocs-awesome-pages-plugin

# Nuke docs/ completely each run to avoid stale symlinks
rm -rf $MKDOCS_DIR/docs
mkdir -p $MKDOCS_DIR/docs/GridKit
ln -sfn $GRIDKIT_SRC/README.md $MKDOCS_DIR/docs/index.md
ln -sfn $GRIDKIT_SRC/GridKit   $MKDOCS_DIR/docs/GridKit/GridKit
ln -sfn $GRIDKIT_SRC/examples  $MKDOCS_DIR/docs/GridKit/examples
ln -sfn $GRIDKIT_SRC/docs      $MKDOCS_DIR/docs/GridKit/docs
mkdir -p $MKDOCS_DIR/docs/GridKit/uq-usecase
# Auto-generate .pages for uq-usecase and all its subdirs to preserve exact names.
# Runs after every build so new dirs are picked up automatically.
echo "title: uq-usecase" > $MKDOCS_DIR/docs/GridKit/uq-usecase/.pages
ln -sfn $GRIDKIT_DIR/uq-usecase/notebooks  $MKDOCS_DIR/docs/GridKit/uq-usecase/notebooks
ln -sfn $GRIDKIT_DIR/uq-usecase/cases      $MKDOCS_DIR/docs/GridKit/uq-usecase/cases
ln -sfn $GRIDKIT_DIR/uq-usecase/work-notes $MKDOCS_DIR/docs/GridKit/uq-usecase/work-notes
for f in kestrel_install.md build_mkdocs.sh 1_build_llvm.sh 2_build_enzyme.sh \
          3_build_suitesparse.sh 4_build_sundials.sh 5_build_ipopt.sh \
          6_build_gridkit.sh rebuild_if_updated.sh; do
    [[ -f $GRIDKIT_DIR/uq-usecase/scripts/$f ]] && ln -sfn $GRIDKIT_DIR/uq-usecase/scripts/$f $MKDOCS_DIR/docs/GridKit/uq-usecase/$f
done
[[ -f $GRIDKIT_DIR/uq-usecase/kestrel_install.md ]] && ln -sfn $GRIDKIT_DIR/uq-usecase/kestrel_install.md $MKDOCS_DIR/docs/GridKit/uq-usecase/kestrel_install.md
# Auto-generate .pages for every subdir so exact dir names appear in nav.
# Written into the real source dirs (symlink targets) since mkdocs follows symlinks.
for d in $GRIDKIT_DIR/uq-usecase/cases $GRIDKIT_DIR/uq-usecase/notebooks $GRIDKIT_DIR/uq-usecase/work-notes; do
    echo "title: $(basename "$d")" > "$d/.pages"
done
for f in README.md INSTALL.md CHANGELOG.md CONTRIBUTING.md; do
    ln -sfn $GRIDKIT_SRC/$f $MKDOCS_DIR/docs/GridKit/$f
done

# Write mkdocs.yml
cat > $MKDOCS_DIR/mkdocs.yml << 'EOF'
site_name: GridKit
docs_dir: docs
use_directory_urls: false
hooks:
  - hooks.py
theme:
  name: material
  features:
    - navigation.sections
    - navigation.indexes
plugins:
  - awesome-pages
  - search
markdown_extensions:
  - attr_list
  - pymdownx.arithmatex:
      generic: true
  - tables
  - toc:
      permalink: true
extra_javascript:
  - javascripts/mathjax.js
  - https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js
EOF

# Write MathJax config script (must load before the MathJax CDN script)
mkdir -p $MKDOCS_DIR/docs/javascripts
cat > $MKDOCS_DIR/docs/javascripts/mathjax.js << 'JSEOF'
window.MathJax = {
  tex: {
    macros: {
      degree: ['^{\\circ}', 0]
    }
  }
};
JSEOF

# Hook: preprocess Markdown before Python-Markdown sees it
cat > $MKDOCS_DIR/hooks.py << 'EOF'
import re

def on_page_markdown(markdown, page, **kwargs):
    # If every non-empty line starts with whitespace (uniformly indented file — source
    # formatting bug), strip the common leading indent. Affects Illinois README.
    lines = markdown.splitlines()
    non_empty = [l for l in lines if l.strip()]
    if non_empty and all(l[0] in (' ', '\t') for l in non_empty):
        import textwrap
        markdown = textwrap.dedent(markdown)

    # Prepend the repo-relative source path so readers know which file to edit.
    # src_path is already relative to docs_dir (mkdocs/docs/), so:
    #   index.md                           -> GridKit/README.md  (home page)
    #   GridKit/GridKit/Model/.../README.md -> as-is (repo-relative)
    #   GridKit/examples/.../README.md      -> as-is
    src = page.file.src_path
    display_path = 'GridKit/README.md' if src == 'index.md' else src
    markdown = f'*`{display_path}`*\n\n---\n\n' + markdown

    # Fix links in index.md (GridKit/README.md) that reference sibling files like
    # INSTALL.md — in our docs tree those live at GridKit/INSTALL.md, not at root.
    if src == 'index.md':
        markdown = re.sub(r'\]\(INSTALL\.md', '](GridKit/INSTALL.md', markdown)
        markdown = re.sub(r'\]\(CHANGELOG\.md', '](GridKit/CHANGELOG.md', markdown)
        markdown = re.sub(r'\]\(CONTRIBUTING\.md', '](GridKit/CONTRIBUTING.md', markdown)

    # Convert GitHub-style backtick-dollar inline math ($`...`$) to standard $...$
    # which arithmatex understands.
    markdown = re.sub(r'\$`(.*?)`\$', r'$\1$', markdown)

    # Python-Markdown requires a blank line before a list or table that follows
    # non-list/non-table text. Apply BEFORE math fence conversion so that | inside
    # ```math blocks is not mistaken for a table row.
    # Insert blank line before list items (require space after -/*/+ to avoid
    # matching table separators like ---|---):
    markdown = re.sub(
        r'(?m)^([ \t]*(?![-*+]\s|\d+[.)]\s|#+\s|\|)(?:\S[^\n]*))\n([ \t]*(?:[-*+]\s|\d+[.)]\s))',
        r'\1\n\n\2',
        markdown
    )
    # Insert blank line before table rows (lines starting with |):
    markdown = re.sub(
        r'(?m)^([ \t]*(?![-*+]\s|\d+[.)]\s|#+\s|\|)(?:\S[^\n]*))\n([ \t]*\|)',
        r'\1\n\n\2',
        markdown
    )

    # Fix display math blocks (handles both ```math and ``` math, including unclosed at EOF)
    def convert_math(m):
        content = m.group(1)
        # Collapse internal blank lines - they cause Python-Markdown to split the block
        content = re.sub(r'\n\n+', '\n', content)
        # If content uses \\ line breaks but has no \begin{}, wrap in aligned env so
        # Python-Markdown doesn't consume \\ as hard line breaks before arithmatex sees it.
        # Also strip trailing \\ from last content line (redundant in LaTeX, but would
        # otherwise split the closing $$ into a separate paragraph).
        if '\\\\' in content and '\\begin{' not in content:
            content = re.sub(r'\\\\\s*$', '', content)
            content = '\\begin{aligned}\n' + content + '\n\\end{aligned}'
        # Blank lines before/after ensure each block is a separate block element
        return '\n$$\n' + content + '\n$$\n'

    def convert_math_unclosed(m):
        content = m.group(1).rstrip('\n')
        content = re.sub(r'\n\n+', '\n', content)
        if '\\\\' in content and '\\begin{' not in content:
            content = re.sub(r'\\\\\s*$', '', content)
            content = '\\begin{aligned}\n' + content + '\n\\end{aligned}'
        return '\n$$\n' + content + '\n$$\n'

    markdown = re.sub(r'```\s*math\n(.*?)\n```', convert_math, markdown, flags=re.DOTALL)
    # Handle unclosed math fence at end of file (source bug in some READMEs)
    markdown = re.sub(r'```\s*math\n([^`]+)\Z', convert_math_unclosed, markdown, flags=re.DOTALL)

    # Fix figure paths.
    # In the repo, READMEs reference figures as e.g. ../../../docs/Figures/fig.png
    # In the new docs tree (docs_dir=mkdocs/docs/, GridKit/=repo root symlink),
    # figures are served at GridKit/docs/Figures/ from the site root.
    # Rewrite any (\.\./)+docs/Figures/ to the correct depth-relative path.
    depth = page.file.src_path.count('/')
    rel = '../' * depth
    markdown = re.sub(r'(\.\./)+docs/Figures/', rel + 'GridKit/docs/Figures/', markdown)

    # Fix stale src/Model/ links — the model code lives at GridKit/Model/ not src/Model/.
    # e.g. ../../../src/Model/PowerFlow/Branch/README.md -> ../../../GridKit/Model/...
    markdown = re.sub(r'((?:\.\./)+)src/Model/', r'\1GridKit/Model/', markdown)

    # Fix GENSAL image filenames — README uses GENSAL.png / GENSAL_ERROR.png but
    # the actual files are Gensal_validation.png / Gensal_validation_error.png.
    markdown = re.sub(r'\bGENSAL_ERROR\.png\b', 'Gensal_validation_error.png', markdown)
    markdown = re.sub(r'\bGENSAL\.png\b', 'Gensal_validation.png', markdown)

    return markdown

def on_nav(nav, config, files, **kwargs):
    """For pages under uq-usecase/, use the raw filename stem as the nav title."""
    def _fix(items):
        for item in items:
            if hasattr(item, 'file') and item.file is not None:
                if 'GridKit/uq-usecase/' in item.file.src_path:
                    stem = item.file.src_path.rsplit('/', 1)[-1].rsplit('.', 1)[0]
                    item.title = stem
            if hasattr(item, 'children') and item.children:
                _fix(item.children)
    _fix(nav.items)
    return nav
EOF

# Build
cd $MKDOCS_DIR
mkdocs build --site-dir $GRIDKIT_DIR/uq-usecase/mkdocs-site

echo ""
echo "Docs built at: $GRIDKIT_DIR/uq-usecase/mkdocs-site/index.html"
echo ""
echo "To rsync to local (e.g. Mac), run from your local machine:"
# echo "  rsync -avz --delete kestrel:~/gridkit/mkdocs-site/ ~/gridkit-mkdocs/ && open ~/gridkit-mkdocs/index.html"

echo "cd /Users/isatkaus/projects/scidac/gridkit/"
echo "./rsync_and_open_mkdocs.sh"

# Generate Copilot context prompt file listing all GridKit .md files.
# Re-generated each build so it stays current as docs are added/removed.
PROMPT_FILE=/home/isatkaus/.github/prompts/gridkit_docs.prompt.md
mkdir -p "$(dirname "$PROMPT_FILE")"
{
  echo "---"
  echo "name: gridkit_docs"
  echo "description: Answer questions using GridKit documentation and source READMEs"
  echo "mode: ask"
  echo "---"
  echo "You are answering questions about GridKit. The following files are the complete"
  echo "GridKit documentation and source READMEs. Read them to answer questions accurately."
  echo ""
  find "$GRIDKIT_SRC" -name '*.md' ! -path '*/build/*' ! -path '*/third-party/*' | sort | while read -r f; do
    # Use path relative to GRIDKIT_SRC as the label
    label="${f#$GRIDKIT_SRC/}"
    echo "- [$label]($f)"
  done
  echo ""
  echo "User Question: \${input}"
} > "$PROMPT_FILE"
echo "Copilot context prompt written to: $PROMPT_FILE"
echo "  $(grep -c '^- \[' "$PROMPT_FILE") .md files listed."

