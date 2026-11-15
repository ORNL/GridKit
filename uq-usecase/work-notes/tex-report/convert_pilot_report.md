# LaTeX Conversion Workflow: Markdown → LaTeX

## Overview

The `convert_pilot_report.py` script incrementally converts sections of the markdown pilot UQ report into publication-quality LaTeX for inclusion in the SciDAC technical report.

**Key design principle**: Section-by-section conversion with idempotent, deterministic output. Each run parses the full markdown and regenerates from scratch—no state files or checksums.

## Dual-Path Output Architecture

Generated files are saved to **two locations**:

### Primary Output (SciDAC Report Integration)
```
/Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/
├── chapters/uq-subsections/
│   ├── pilot_uq_usecase.tex      ← Main LaTeX file
│   ├── figs/                      ← Figures (.png)
│   └── uq-usecase_refs.bib       ← Bibliography
```

**Purpose**: Live editing and integration into main SciDAC report

### Secondary Output (Repository Version Control)
```
/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/tex-report/uq-subsections/
├── pilot_uq_usecase.tex          ← LaTeX version in repo
├── figs/                          ← Figures (mirrored)
└── uq-usecase_refs.bib           ← Bibliography in repo
```

**Purpose**: Track LaTeX generation history; reproducible builds from repo

## Phase Structure

Phases are incremental; running a higher phase generates all prior sections:

| Phase | Sections | Description |
|-------|----------|-------------|
| 1     | §1       | "What this is" (text only) |
| 2     | §1-2     | + Glossary table + footnotes |
| 3     | §1-3     | + The system (with diagram) |
| 4     | §1-4     | + Building condition (diagram + table) |
| 5     | §1-4.3   | + Wind becomes negative load |
| 6     | §1-4.4   | + The datacenter (hunt summary) |
| 7-13  | §5-12    | Future phases (not yet implemented) |

## Usage

### Generate Phase 3 (default: sections 1-3)
```bash
cd /Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/tex-report
python convert_pilot_report.py
```

### Generate Specific Phase
```bash
python convert_pilot_report.py phase=6
```

### Output
```
✓ Read markdown file (67795 bytes)
✓ Extracted §1 'What this is' (824 bytes)
✓ Generated §1 LaTeX (987 bytes)
...
✓ Wrote LaTeX file: /Users/isatkaus/projects/scidac/docs/report/.../pilot_uq_usecase.tex
✓ Wrote LaTeX file: /Users/isatkaus/projects/scidac/repos/gridkit/.../pilot_uq_usecase.tex
✓ Copied figure: /Users/isatkaus/projects/scidac/docs/report/.../figs/fig:3_...png
✓ Copied figure: /Users/isatkaus/projects/scidac/repos/gridkit/.../figs/fig:3_...png
✓ Copied bibliography file to repo: /Users/isatkaus/projects/scidac/repos/gridkit/.../uq-usecase_refs.bib

[SUCCESS] Phase 6 conversion complete!
```

## Configuration

Paths are hardcoded in `convert_pilot_report.py` (lines 25-37):

```python
# Primary output (SciDAC report)
TARGET_PATH = "...chapters/uq-subsections/pilot_uq_usecase.tex"
FIGURES_OUTPUT_DIR = "...chapters/uq-subsections/figs/"

# Secondary output (repo)
REPO_OUTPUT_PATH = "...tex-report/uq-subsections/pilot_uq_usecase.tex"
REPO_FIGURES_OUTPUT_DIR = "...tex-report/uq-subsections/figs/"
REPO_BIBTEX_PATH = "...tex-report/uq-subsections/uq-usecase_refs.bib"
```

## Citation Mapping

The script maps numeric markdown citations `[1]`, `[2]`, etc. to BibTeX keys:

**Markdown**:
```markdown
...controllable current source [1] on a problem...
```

**Converted LaTeX**:
```latex
...controllable current source \cite{Elliott2021} on a problem...
```

**Compiled PDF** (with BibTeX):
```
...controllable current source [1] on a problem...
```

Map is defined in `CITATION_MAP` (lines 52-75) with 23 BibTeX entries.

## PDF Compilation

### Standalone Test Document
Test the LaTeX output with the included `test_main.tex`:

```bash
cd /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections

# Full BibTeX cycle
rm -f test_main.pdf test_main.aux test_main.bbl test_main.blg
pdflatex -interaction=nonstopmode test_main.tex
bibtex test_main
pdflatex -interaction=nonstopmode test_main.tex
pdflatex -interaction=nonstopmode test_main.tex

# View
open test_main.pdf
```

### Integration into Main Report
Once `pilot_uq_usecase.tex` is complete, include in main SciDAC report via:

```latex
% In /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq.tex
\input{chapters/uq-subsections/pilot_uq_usecase}
```

## Design Patterns

### 1. Markdown Section Extraction
Each section has dedicated extraction function:
- `extract_section_1()` → `extract_section_4_4()`
- Handles subsections, tables, blockquotes, lists, figures

### 2. LaTeX Generation
Each extraction pairs with a generation function:
- `generate_latex_section_1()` → `generate_latex_section_4_4()`
- Wraps extracted content in LaTeX commands
- Applies formatting transformations (markdown → LaTeX)

### 3. Markdown-to-LaTeX Conversion
Helper function `markdown_to_latex_text()` handles inline transformations:
- `**bold**` → `\textbf{bold}`
- `*italic*` → `\textit{italic}`
- `` `code` `` → `\texttt{code}` (or `$code$` in math mode)
- `[link](url)` → `\href{url}{link}`
- `[1]` → `\cite{Elliott2021}` (via `CITATION_MAP`)

### 4. Figure Handling
Pre-generated Mermaid diagrams (.png) are copied to both output directories:
- `extract_and_generate_mermaid.py` generates PNG files in `figures/`
- `convert_pilot_report.py` copies them to both `figs/` directories

### 5. Bibliography
- Single source: `uq-usecase_refs.bib` (23 entries)
- Copied to both output locations for compilation
- Style: `plain.bst` (numeric [1], [2], etc.)

## Workflow Example: Adding New Section

To add section 5 (Phase 7):

1. **Extract function** (lines ~1900+):
   ```python
   def extract_section_5(md_content):
       # Find section, extract text/tables/blockquotes
       # Return (text, tables, figures)
   ```

2. **Generation function** (lines ~2000+):
   ```python
   def generate_latex_section_5(text, tables, figures):
       # Wrap in \section{}, \subsection{}, etc.
       # Apply markdown_to_latex_text() to inline content
       # Include figures with \includegraphics
       return latex_output
   ```

3. **Phase 7 handler** (lines ~2100+):
   ```python
   elif phase_int == 7:
       # Extract section 5
       text, tables, figures = extract_section_5(md_content)
       latex_section_5 = generate_latex_section_5(...)
       # Add to latex_sections list
   ```

4. **Update phase guidance** (lines ~2250+):
   ```python
   elif phase_int < 8:
       print("2. Run again: python convert_pilot_report.py phase=8 (add Section 6...)")
   ```

## Troubleshooting

### PDF citations show as "?"
- Verify BibTeX file path in LaTeX `\bibliography{...}` command
- Run full 4-pass BibTeX cycle (see PDF Compilation section)
- Check bibtex log: `cat test_main.blg | grep -i error`

### Figure not found
- Ensure `extract_and_generate_mermaid.py` has been run
- Check `figures/` directory for `.png` file
- Run: `python extract_and_generate_mermaid.py`

### LaTeX compilation errors
- Check log file: `tail -50 test_main.log | grep -i error`
- Common issues: unmatched braces `{}`, invalid Unicode, missing packages
- Validate markdown formatting before conversion

### Margins too large/small
- `test_main.tex` uses: `\usepackage[margin=1in]{geometry}`
- Matches main SciDAC report style
- To adjust: edit line in `test_main.tex` (search `margin=`)

## Files Modified/Created

### Core Scripts
- `convert_pilot_report.py` — Main conversion engine (~2300 lines)
- `extract_and_generate_mermaid.py` — Mermaid diagram generator

### Generated Outputs
- `pilot_uq_usecase.tex` — Primary LaTeX (SciDAC report)
- `uq-subsections/pilot_uq_usecase.tex` — Repo version (git)
- `pilot_uq_usecase.pdf` — Compiled output (for review)
- `figs/*.png` — Figures (mirrored to both locations)
- `uq-usecase_refs.bib` — Bibliography (mirrored to both locations)

### Test/Support
- `test_main.tex` — Standalone test wrapper
- `test_main.pdf` — Test output (for validation)
- `convert_pilot_report.md` — This file (documentation)

## References

- Source markdown: `/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/datacenter_pilot_report.md`
- Branch: `isatkaus/uq-usecase`
- Repository: `https://github.com/ORNL/GridKit`
