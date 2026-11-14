# Two-Script Architecture for Mermaid→LaTeX Conversion

## Overview

The conversion from markdown to LaTeX now uses **two separate, independent scripts** with clear separation of concerns:

1. **Mermaid Extraction**: `extract_and_generate_mermaid.py` — Diagram generation pipeline
2. **LaTeX Conversion**: `convert_pilot_report.py` — Markdown to LaTeX conversion

## Architecture Diagram

```
datacenter_pilot_report.md
    ↓
    ├─→ extract_and_generate_mermaid.py ──→ PNG files (work-notes/figures/)
    │                                           ↓
    └─→ convert_pilot_report.py ──(figures)──→ pilot_uq_usecase.tex
                                 ──(figures)──→ figs/
                                    ↑
                                    └─ copies PNG to report output dir
```

## Script Responsibilities

### extract_and_generate_mermaid.py

**What it does:**
1. Scans markdown document for ALL mermaid diagram blocks (```mermaid...```)
2. Extracts caption from nearby italic text (e.g., `*Figure 3.1. ...* `)
3. Generates consistent label from section number + caption slug
4. Converts mermaid code to PNG using mermaid-cli
5. Saves PNG to FIGURES_SOURCE_DIR with consistent naming

**Output:**
- PNG files in: `/Users/.../work-notes/figures/`
- Named: `fig:3_one_time_seed_building_blue_and_the_pe.png`
- One file per mermaid diagram (4 found in current document)

**Key Features:**
- ✅ Finds ALL mermaid blocks (handles §1-10, not just current section)
- ✅ Generates consistent file naming based on section + caption
- ✅ Idempotent: Re-running regenerates all PNGs from scratch
- ✅ Fallback placeholders if mermaid-cli unavailable

### convert_pilot_report.py

**What it does:**
1. Reads markdown document sections
2. For each section:
   - Extracts text, captions, blockquotes
   - Converts markdown formatting to LaTeX
   - **Looks up pre-generated PNG file** using same label logic
   - **Copies PNG to output figures directory**
   - Generates LaTeX includegraphics reference with relative path
3. Writes complete LaTeX file

**Changes from Previous Version:**
- ❌ NO LONGER generates mermaid→PNG (that's `extract_and_generate_mermaid.py`'s job)
- ❌ NO LONGER runs mermaid-cli
- ✅ ONLY references pre-generated PNG files
- ✅ ONLY copies PNGs to output directory
- Cleaner, simpler, more focused

**Figure Reference Flow:**
```
Caption in markdown: "*Figure 3.1. One-time seed building (blue)...*"
    ↓
Label generation: "fig:3_one_time_seed_building_blue_and_the_pe"
    ↓
File lookup: FIGURES_SOURCE_DIR / "fig:3_one_time_seed_building_blue_and_the_pe.png"
    ↓
Copy to: FIGURES_OUTPUT_DIR / same filename
    ↓
LaTeX ref: \includegraphics{figs/fig:3_one_time_seed_building_blue_and_the_pe.png}
```

## Label Generation (Consistent)

Both scripts use identical logic to generate labels:

**Input:** Caption text (with "Figure N.N. " prefix and special characters)
```
"Figure 3.1. One-time seed building (blue) and the per-condition build (grey)."
```

**Processing:**
1. Remove "Figure N.N. " prefix → "One-time seed building (blue)..."
2. Limit to first ~40 chars → "One-time seed building (blue) and th"
3. Convert to lowercase → "one-time seed building (blue) and th"
4. Remove special chars, spaces→underscores → "one_time_seed_building_blue_and_th"
5. Combine with section number → "fig:3_one_time_seed_building_blue_and_th"

**Output:** `fig:3_one_time_seed_building_blue_and_the_pe.png`

This ensures BOTH scripts generate identical labels → PNGs found/referenced correctly.

## Workflow (Day-to-Day Usage)

### Initial Setup (One Time)
```bash
cd /Users/.../work-notes/tex-report/
python extract_and_generate_mermaid.py    # Generate all PNGs
python convert_pilot_report.py             # Convert to LaTeX
```

### When Markdown Changes
```bash
# Option 1: If only text changed (no diagram changes)
python convert_pilot_report.py             # Re-convert to LaTeX

# Option 2: If diagrams or captions changed
python extract_and_generate_mermaid.py    # Regenerate PNGs
python convert_pilot_report.py             # Re-convert LaTeX
```

### Testing
```bash
cd /Users/.../uq-subsections/
pdflatex test_main.tex
open test_main.pdf
```

## Benefits of This Architecture

### 1. **Separation of Concerns**
- Diagram generation isolated from markdown→LaTeX conversion
- Each script has a single, clear responsibility
- Easier to test, debug, maintain

### 2. **Scalability**
- Adding more mermaid blocks? Just add to markdown
- `extract_and_generate_mermaid.py` auto-detects ALL sections
- No changes needed to conversion script
- Scales to §1-10 without modification

### 3. **Reusability**
- PNGs generated once, used forever
- Can regenerate diagrams without reconverting markdown
- Can update LaTeX style without regenerating PNGs
- Independent workflows for different tasks

### 4. **Debugging**
- Can inspect generated PNG files directly
- Can verify mermaid→PNG conversion independently
- Can test LaTeX conversion without mermaid issues
- Clear error messages if figures missing

### 5. **Idempotency**
- Both scripts fully idempotent (same output every run)
- No state files or checksums needed
- Can safely re-run multiple times
- Deterministic: identical markdown → identical output

### 6. **Future-Proof**
- Ready to handle §4-10 with no code changes
- Easy to add alternative formats (TikZ, SVG, etc.)
- Can add preprocessing/post-processing passes as needed
- Modular design allows plugin architectures

## Data Flow Example (Section 3)

```
Markdown Source:
  ```mermaid
  graph TD
    A[...]
  ```
  *Figure 3.1. One-time seed building (blue)...*

  ↓ [extract_and_generate_mermaid.py]

Generated PNG:
  /work-notes/figures/fig:3_one_time_seed_building_blue_and_the_pe.png (23 KB)

  ↓ [convert_pilot_report.py]

LaTeX Output:
  \begin{figure}[htb]
    \centering
    \includegraphics[width=0.9\textwidth]{figs/fig:3_one_time_seed_building_blue_and_the_pe.png}
    \caption{One-time seed building (blue)...}
    \label{fig:3_one_time_seed_building_blue_and_the_pe}
  \end{figure}

  ↓ [pdflatex]

Final PDF:
  6 pages, 241 KB, figure rendered correctly
```

## Error Handling

### Mermaid Extraction Fails
→ Script reports which diagram failed
→ Tries fallback: PIL placeholder PNG
→ If PIL unavailable: clear error message
→ User action: Install mermaid-cli or Pillow

### Figure Not Found During Conversion
→ Script reports: "Pre-generated figure not found: fig:..."
→ Clear instruction: "Run: python extract_and_generate_mermaid.py"
→ No silent failures; always clear feedback

### LaTeX Compilation Fails
→ Check figure files exist: `ls figs/`
→ Verify PNG format: `file figs/*.png`
→ Check paths in LaTeX file: `grep includegraphics pilot_uq_usecase.tex`
→ Standard pdflatex debugging applies

## Files Involved

```
Scripts:
  /Users/.../tex-report/extract_and_generate_mermaid.py    (NEW)
  /Users/.../tex-report/convert_pilot_report.py            (UPDATED)

Source Data:
  /Users/.../work-notes/datacenter_pilot_report.md

Intermediate:
  /Users/.../work-notes/figures/fig:*.png                  (Generated)

Output:
  /Users/.../uq-subsections/pilot_uq_usecase.tex
  /Users/.../uq-subsections/figs/fig:*.png                 (Copied)
  /Users/.../uq-subsections/test_main.pdf                  (For testing)
```

## Test Results

✅ **Mermaid extraction**: 4/4 diagrams found and converted to PNG
✅ **Label generation**: Consistent across both scripts
✅ **Figure copying**: PNG successfully copied to output directory
✅ **LaTeX compilation**: PDF generated (6 pages, 241 KB)
✅ **Figure rendering**: Diagram visible in PDF

## Future Enhancements

### Possible Extensions
1. **Caching**: Skip re-generation if source hasn't changed
2. **Formats**: Add TikZ, SVG output options
3. **Validation**: Pre-check mermaid syntax before mermaid-cli
4. **Logging**: Detailed JSON output for CI/CD integration
5. **Parallelization**: Generate multiple PNGs concurrently

### For Phase 4+
- Add §4 section with new mermaid block to markdown
- Run: `python extract_and_generate_mermaid.py` (auto-detects)
- Run: `python convert_pilot_report.py phase=4`
- New diagram automatically integrated!

## Conclusion

This two-script architecture provides:
- ✅ Clear separation of diagram generation and LaTeX conversion
- ✅ Scalable to any number of sections and figures
- ✅ Idempotent and deterministic workflows
- ✅ Easy debugging and maintenance
- ✅ Future-proof for phase 4-10 implementation
