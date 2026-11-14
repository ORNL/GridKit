# Markdown-to-LaTeX Report Conversion: convert_pilot_report.py

**Purpose:** Automatically convert `datacenter_pilot_report.md` (Markdown) → `pilot_uq_usecase.tex` (LaTeX) for inclusion in the SciDAC technical report.

**Status:** Phase 1 (MVP) implemented. Converts §1 "What this is" (text only).

---

## Quick Start

### Run the Script

```bash
cd /Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase
python convert_pilot_report.py
```

**Output:** 
```
[Phase 1] Converting datacenter_pilot_report.md → pilot_uq_usecase.tex
Source: /Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/datacenter_pilot_report.md
Target: /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections/pilot_uq_usecase.tex
✓ Read markdown file
✓ Extracted §1 'What this is'
✓ Generated LaTeX
✓ Target directory ready
✓ Wrote LaTeX file
[SUCCESS] Phase 1 conversion complete!
```

### Quick Test (Without Full Report)

A minimal test document is provided for rapid iteration **without rebuilding the full SciDAC report**:

```bash
cd /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections
pdflatex test_main.tex
open test_main.pdf  # or view in your PDF reader
```

This generates a standalone PDF with your section included. Use this for testing during development.

**Files:**
- `test_main.tex` — Minimal test harness (includes `pilot_uq_usecase.tex`)
- `test_main.pdf` — Compiled output (regenerated each run)

### Integrate into Full Report

After running the script and verifying via `test_main.tex`, add this line to `/Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq.tex` at the appropriate location (likely after the existing aleatoric section):

```latex
\input{chapters/uq-subsections/pilot_uq_usecase}
```

Then compile the full report:

```bash
cd /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report
pdflatex main.tex
```

---

## Configuration

All paths are defined as variables at the top of `convert_pilot_report.py`. Adjust these if the file locations change:

```python
SOURCE_MD = "/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/datacenter_pilot_report.md"
FIGURES_DIR = "/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/figures/"
TARGET_PATH = "/Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections/pilot_uq_usecase.tex"
```

| Variable | Purpose | Notes |
|---|---|---|
| `SOURCE_MD` | Path to input markdown file | Phase 1 uses §1 section only; future phases expand to §2–10 |
| `FIGURES_DIR` | Directory containing PNG figures | Used by Phase 3+ to locate and copy figures |
| `TARGET_PATH` | Path where LaTeX output is written | Creates directories if needed |

---

## How It Works

### Design Philosophy

**Idempotent and stateless:** Each run reads the full markdown file and regenerates the entire LaTeX output from scratch. No intermediate state files, checksums, or dependencies on prior runs.

**Why this approach:**
- Simple: No state tracking complexity
- Safe: Can re-run anytime; markdown changes handled automatically
- Transparent: Script logic is straightforward, easy to debug
- Future-proof: Phases 2–5 simply add new conversion logic without changing the pattern

**Key benefit:** When the markdown report is updated (new paragraph, caption change, new figure, etc.), just run the script again and all changes flow through to LaTeX automatically.

### Conversion Pipeline

```
datacenter_pilot_report.md
        ↓
   [Extract §1 by regex]
        ↓
   [Markdown→LaTeX formatting]
        ↓
   [Wrap in \subsection]
        ↓
pilot_uq_usecase.tex
```

### Phase 1: Text-Only Conversion (§1 "What this is")

**Input:** Section §1 of markdown  
**Output:** LaTeX `\subsection{Pilot Aleatoric UQ Usecase}` with formatted text  
**Conversions:**
- `## 1. What this is` → Removed (wrapped in `\subsection{}` LaTeX command)
- Lines starting with `>` (blockquotes) → `\begin{quote}...\end{quote}`
- `**bold**` → `\textbf{bold}`
- `*italic*` → `\textit{italic}`
- Blank lines preserved (standard LaTeX paragraph separation)
- No tables, figures, code blocks, or links (handled in later phases)

**Verification:**
1. Script runs without errors
2. Output file created at target path
3. Manual inspection: Open generated `.tex` file, verify §1 text appears correctly
4. Compile test: Add `\input{}` to `uq.tex`, run `pdflatex`, check for errors
5. Visual inspection: PDF shows formatted text, bold/italic correct

---

## Implementation Phases

### Phase 1: Text-Only (DONE ✓)
- **Sections:** §1 "What this is"
- **Features:** Basic text, bold, italic, blockquotes
- **Time:** 1.5–2.5 hours
- **Status:** ✅ Implemented, tested, blockquotes working

### Phase 2: Tables (TODO)
- **Sections:** §2 Glossary (21-row table)
- **Features:** Markdown table → LaTeX `tabular` environment
- **Conversion logic:** Parse markdown table syntax, generate `\begin{tabular}...\end{tabular}`
- **Estimated time:** 1–1.5 hours
- **Dependencies:** None; can be done independently
- **Verification:** Check table formatting, column alignment, row heights

### Phase 3: Figures & Mermaid Diagrams (TODO)
- **Sections:** §3–7 (with embedded figures and diagrams)
- **Features:** 
  - PNG figures: `![caption](figures/file.png)` → `\includegraphics[width=...]{figs/ignas/...}`
  - Mermaid diagrams: Extract code blocks, render to PNG via `mmdc`, include as figures
- **Figure files:** All in `/gridkit/uq-usecase/work-notes/figures/`
- **Dependencies:** `mmdc` (mermaid-cli) for diagram rendering
- **Estimated time:** 2–3 hours
- **Verification:** All figures appear in output, captions correct, cross-references (`\ref{fig:...}`) work

### Phase 4: Math, Code, Footnotes (TODO)
- **Sections:** §5–10 (math formulas, inline code, cross-document references)
- **Features:**
  - Math: `$...$` and `$$...$$` → LaTeX math mode (mostly compatible)
  - Code: `` `code` `` → `\texttt{code}`
  - Links to details: `[Details](datacenter_load_uq_plan.md#anchor)` → Footnote with gridkit repo reference
- **Estimated time:** 2–2.5 hours
- **Verification:** Math renders correctly in PDF, code formatting visible, footnotes appear

### Phase 5: Full Integration & Testing (TODO)
- **Steps:**
  1. Add `\input{}` to `uq.tex`
  2. Run `pdflatex main.tex` from report root
  3. Visual inspection: figures, tables, math, page breaks
  4. Author review (Ignas)
- **Estimated time:** 0.5–1 hour
- **Verification:** Full report compiles, all content appears correctly

---

## Process Notes

### Why Idempotent?

The markdown source will change over time (new content, clarifications, figure updates, etc.). The script is designed to handle all changes without manual intervention:

- **New paragraph?** Re-run script → new text added automatically
- **Changed caption?** Re-run script → updated caption in LaTeX
- **Added figure?** Re-run script → new figure included
- **Link to details changed?** Re-run script → footnote updated

No manual merging, no state tracking, no lost updates.

### Markdown Source Structure

The markdown file `datacenter_pilot_report.md` is organized as:

```
- Front matter (title, status, table of contents)
- ## 1. What this is          ← Phase 1
- ## 2. Glossary              ← Phase 2
- ## 3. The system            ← Phase 3
- ## 4. Building one condition ← Phase 3
- ## 5. CCT by bisection       ← Phase 4
- ## 6. Parallel execution     ← Phase 3
- ## 7. Results                ← Phase 3–4
- ## 8. GridKit feedback       ← Phase 4
- ## 9. Limitations            ← Phase 4
- ## 10. ACTIVSg2000           ← Phase 4
```

Each phase extracts and converts one or more sections. The script uses regex to identify section boundaries.

### Cross-References to Plan Document

The markdown contains links to `datacenter_load_uq_plan.md` marked as **Details**. Example:

```markdown
[Details](datacenter_load_uq_plan.md#where-the-wind-timeseries-came-from-and-what-it-looks-like)
```

**Phase 4** will convert these to LaTeX footnotes:

```latex
\footnote{See detailed analysis in gridkit uq-usecase branch, section: 
``where the wind timeseries came from and what it looks like''}
```

This keeps the main text clean while allowing readers to follow links to the repository if interested.

---

## Troubleshooting

### Script won't find the markdown file

**Error:** `ERROR: Source markdown file not found: ...`

**Solution:** Verify the path in `SOURCE_MD` variable. Check:
- File exists: `ls /Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/datacenter_pilot_report.md`
- Correct branch: Ensure you're on the `uq-usecase` branch

### Section §1 not found

**Error:** `ERROR: Could not extract §1: Could not find '## 1. What this is' section`

**Solution:** Check the markdown file for the exact heading text. Markdown is whitespace-sensitive:
- Should be exactly `## 1. What this is` (note the space after `##`)
- If heading was changed, update the regex in `extract_section_1()` function

### Generated LaTeX has formatting issues

**Symptoms:** Bold/italic not appearing, extra backslashes, malformed LaTeX

**Solution:** 
1. Check raw markdown has proper formatting: `**bold**`, `*italic*`
2. Review the regex patterns in `markdown_to_latex_text()` function
3. Run a small test with a simplified markdown snippet
4. Compare generated `.tex` file with expected output

### Can't write to target directory

**Error:** `ERROR: Could not create target directory: Permission denied`

**Solution:** Check permissions:
```bash
ls -ld /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/
chmod u+w /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/
```

Or ensure the directory exists:
```bash
mkdir -p /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections
```

### LaTeX compilation fails after integration

**Error:** `\input{chapters/uq-subsections/pilot_uq_usecase}` causes pdflatex errors

**Solution:**
1. Check path is relative to `uq.tex` location (should resolve from `chapters/` directory)
2. Ensure the generated `.tex` file is syntactically valid:
   ```bash
   pdflatex pilot_uq_usecase.tex  # Test file independently
   ```
3. Check for special characters in markdown that need escaping (e.g., `_`, `&`, `%`, `#`)

---

## Future Improvements

### Code Quality
- Add unit tests for markdown→LaTeX conversions
- Test with various markdown edge cases (nested emphasis, escaped characters, etc.)
- Add command-line argument parsing (`--phase 1`, `--source`, `--target`, etc.)
- Add logging levels (verbose, debug, quiet)

### Robustness
- Add checksum/logging to detect if markdown has changed since last run
- Optional: cache intermediate results (still deterministic, just faster)
- Better error messages with line numbers pointing to problematic markdown

### Usability
- Add `--dry-run` flag to preview LaTeX without writing
- Add `--validate` flag to check generated LaTeX without integration
- Add `--diff` flag to show changes from previous version
- Auto-detect and warn about unhandled markdown features (code blocks, blockquotes, etc.)

### Features for Later Phases
- Phase 2: Auto-detect table column widths from content
- Phase 3: Auto-render Mermaid diagrams in parallel (faster for many diagrams)
- Phase 3: Auto-copy figure files to `figs/ignas/` with name collision detection
- Phase 4: Auto-generate footnote anchors and validate cross-references
- All phases: Support for markdown footnotes and citations (pandoc AST integration)

---

## Developer Notes

### Regex Patterns Used

| Pattern | Purpose |
|---|---|
| `^## 1\.\s+What this is\s*$` | Locate start of §1 |
| `^## 2\.` | Locate end of §1 (start of §2) |
| `\*\*([^\*]+?)\*\*` | Match `**bold**` text |
| `(?<!\*)\*([^\*]+?)\*(?!\*)` | Match `*italic*` text (not `**`) |

### Adding a New Phase

1. Create a new function `extract_section_N()` (similar to `extract_section_1()`)
2. Create a new function `generate_latex_section_N()` (similar to `generate_latex_section_1()`)
3. Add conversion logic to `markdown_to_latex_text()` or create a specialized function
4. Update `convert_pilot_report()` to call the new phase
5. Update this documentation with phase details
6. Test with a subset of the markdown, then the full file

Example for Phase 2 (tables):

```python
def extract_section_2(md_content, section_1_end):
    """Extract section §2 Glossary (which is a table)."""
    # Use similar regex pattern to extract §2...
    pass

def convert_markdown_table_to_latex(table_markdown):
    """Parse markdown table, generate LaTeX tabular."""
    # Parse rows and cells, generate \begin{tabular}...\end{tabular}
    pass
```

### Testing a Phase

Small test before running on full file:

```bash
# Extract just §1 from markdown and save to test file
head -120 datacenter_pilot_report.md > test_section_1.md

# Modify script to read from test_section_1.md, output to /tmp/test.tex
# Run and check output
python convert_pilot_report.py
cat /tmp/test.tex
```

---

## Contact & History

- **Created:** Phase 1 MVP, November 2026
- **Author:** Converted by Copilot from markdown source
- **Owner:** Ignas Satkauskas (original report author)
- **Repository:** `/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/` (uq-usecase branch)
- **Report root:** `/Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/`

For questions or updates to this process, refer to the plan document at `/memories/session/plan.md`.
