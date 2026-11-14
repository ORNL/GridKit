# LaTeX Report Conversion: Pilot UQ Usecase

This directory contains scripts and documentation for converting the markdown pilot UQ report (`datacenter_pilot_report.md`) into LaTeX format for inclusion in the SciDAC technical report.

## Two-Script Architecture

The conversion process uses **two separate scripts** with a clear separation of concerns:

### 1. **extract_and_generate_mermaid.py** — Diagram Generation
- **Purpose**: Extract all mermaid diagram blocks from markdown, convert to PNG
- **Output**: Pre-generated PNG files in `/work-notes/figures/`
- **When to run**: Once before conversion, or when diagrams change
- **Idempotent**: Re-running regenerates all diagrams from scratch

### 2. **convert_pilot_report.py** — LaTeX Conversion
- **Purpose**: Convert markdown sections to LaTeX, reference pre-generated figures
- **Input**: Markdown + pre-generated PNGs
- **Output**: LaTeX file + copied figures to report directory
- **Workflow**: Read markdown → extract text/captions → reference pre-generated PNG → generate LaTeX

## Directory Structure

```
tex-report/
├── convert_pilot_report.py              # Main LaTeX conversion script
├── extract_and_generate_mermaid.py      # Mermaid diagram extraction (NEW)
├── README.md                             # This file
└── PHASES.md                             # Phase documentation

work-notes/
├── datacenter_pilot_report.md           # Source markdown
├── tex-report/                          # Scripts (this directory)
└── figures/                             # Generated PNG files (from mermaid extraction)

SciDAC report output:
├── chapters/uq-subsections/
│   ├── pilot_uq_usecase.tex             # Generated LaTeX (input via \input{})
│   ├── figs/                            # Copied figures for report
│   │   └── fig:3_one_time_seed_building_blue_and_the_pe.png
│   ├── test_main.tex                    # Test harness
│   └── test_main.pdf                    # Generated PDF (for testing)
```

## Usage

### Step 1: Extract and Generate Mermaid Diagrams

```bash
cd /Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/tex-report
python extract_and_generate_mermaid.py
```

**Output**:
```
[Mermaid Extraction] Scanning: .../datacenter_pilot_report.md
Output directory: .../figures/

Found 4 mermaid diagram(s):

  [3] fig:3_one_time_seed_building_blue_and_the_pe
      Caption: One-time seed building (blue) and the per-condition build (grey).
      → .../figures/fig:3_one_time_seed_building_blue_and_the_pe.png
      ✓ Generated from mermaid-cli
  [4] fig:4_per_condition_pipeline_from_hour_and_co
      ...
```

### Step 2: Convert Markdown to LaTeX

```bash
python convert_pilot_report.py              # Run all phases
python convert_pilot_report.py phase=1      # Run specific phase
```

**Output**:
```
[Phase 3] Converting datacenter_pilot_report.md → pilot_uq_usecase.tex
...
  ✓ Copied figure: fig:3_one_time_seed_building_blue_and_the_pe.png
✓ Generated §3 LaTeX (2935 bytes)
...
[SUCCESS] Phase 3 conversion complete!
```

### Step 3: Test PDF Compilation (Optional)

```bash
cd /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections
pdflatex test_main.tex
open test_main.pdf
```

## Configuration

Edit configuration in each script:

### extract_and_generate_mermaid.py (lines ~20-30)
```python
SOURCE_MD = "...datacenter_pilot_report.md"
FIGURES_OUTPUT_DIR = ".../work-notes/figures/"
OUTPUT_FORMAT = "png"  # or "tikz" (future)
```

### convert_pilot_report.py (lines ~20-40)
```python
SOURCE_MD = "...datacenter_pilot_report.md"
FIGURES_SOURCE_DIR = ".../work-notes/figures/"  # Read pre-generated PNGs
TARGET_PATH = ".../pilot_uq_usecase.tex"        # LaTeX output
FIGURES_OUTPUT_DIR = ".../uq-subsections/figs/" # Copy figures here
INCLUDE_FOOTNOTES = True
REPO_BRANCH = "isatkaus/uq-usecase"
REPO_BASE = "https://github.com/ORNL/GridKit"
```

## Phases

- **Phase 1**: §1 "What this is" - Text with blockquotes, bold, italic
- **Phase 2**: §2 "Glossary" - Table → definition-list with footnotes
- **Phase 3**: §3 "The system" - Text with mermaid diagram (PNG) and details boxes
- **Phase 4+**: Future sections §4-10

## Features

### Markdown to LaTeX Conversion
- **Backticks** (`` `code` ``) → math mode (`$code$`) with underscore escaping
- **Bold** (`**text**`) → `\textbf{text}`
- **Italic** (`*text*`) → `\textit{text}`
- **Blockquotes** (`> text`) → `\begin{quote}...\end{quote}`
- **Subsections** (`### text`) → `\subsubsection{text}`
- **Markdown links** (`[text](url)`) → `\url{url}` (text)
- **Unicode** (§, →, •, –, —) → LaTeX equivalents

### Figure Handling

**Figure label generation** (consistent across scripts):
```
Caption: "One-time seed building (blue) and the per-condition build (grey)."
Label: "fig:3_one_time_seed_building_blue_and_the_pe"
```

- **Extraction**: `extract_and_generate_mermaid.py` finds all mermaid blocks and generates PNGs
- **Conversion**: `convert_pilot_report.py` matches captions to PNGs using same label logic
- **References**: LaTeX uses relative paths `\includegraphics{figs/fig:3_...png}`

### Special Handling
- **Tables**: Markdown tables → LaTeX definition-lists (better for long content)
- **Pipe handling**: Correctly handles `|` in backticks and math mode
- **Multiple links**: Details boxes support multiple section links
- **Flexible formatting**: Robust pattern handling for edge cases

## Design Principles

1. **Idempotent**: Each run regenerates all output from scratch
2. **Deterministic**: Same input always produces identical output
3. **Configuration-driven**: All paths/settings as module variables
4. **Separated concerns**: Diagram generation ≠ LaTeX conversion
5. **Preprocessing-first**: Fix markdown issues before scripts process
6. **Simple**: Rely on standard LaTeX, no fancy packages

## Troubleshooting

### Mermaid diagram issues

**Problem**: Placeholder PNG instead of real diagram
- **Solution**: Install mermaid-cli with `npm install -g @mermaid-js/mermaid-cli`
- Re-run: `python extract_and_generate_mermaid.py`

**Problem**: PIL not available for placeholder
- **Solution**: Install with `pip install Pillow`

### LaTeX conversion issues

**Problem**: Figure not found during conversion
- **Solution**: Run mermaid extraction first: `python extract_and_generate_mermaid.py`
- Then re-run: `python convert_pilot_report.py`

**Problem**: PDF compilation errors
- Check: Do both `figures/` and `figs/` directories exist?
- Verify: Are PNG files valid? (`file *.png`)
- Ensure: Markdown source is valid UTF-8

## Dependencies

- **Python 3** (standard library only: `re`, `sys`, `pathlib`, `subprocess`)
- **pdflatex** (TeX Live 2026)
- **mermaid-cli** (optional, for PNG generation: `npm install -g @mermaid-js/mermaid-cli`)
- **Pillow** (optional, for better placeholders: `pip install Pillow`)

## Integration with Main Report

To include the converted content in the main SciDAC report:

```latex
% In uq.tex or chapters file:
\input{chapters/uq-subsections/pilot_uq_usecase}
```

## Development Notes

- Script uses only standard library for portability
- PNG generation optional (falls back gracefully)
- All paths use absolute paths to avoid confusion
- Test harness: `test_main.tex` provides standalone compilation
- Future: Can extend to handle §4-10 with same architecture

## Next Steps (Phase 4+)

When ready to implement §4 and beyond:
1. Add new mermaid blocks to markdown
2. Run: `python extract_and_generate_mermaid.py` (auto-detects all sections)
3. Run: `python convert_pilot_report.py` with desired phase
4. Script will automatically handle new figures and sections

