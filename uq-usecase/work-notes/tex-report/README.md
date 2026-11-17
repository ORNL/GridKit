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

Two different roots are involved — keep them straight:

**1. Repo root (source of truth, version-controlled).** The scripts live here, and
this is where the clean **upload copy** is written.
Root: `/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/`

```
work-notes/
├── datacenter_pilot_report.md            # SOURCE markdown — this is what you edit
├── figures/                              # generated mermaid PNGs (from script 1)
└── tex-report/                           # scripts + docs (this directory)
    ├── convert_pilot_report.py           # markdown -> LaTeX converter (script 2)
    ├── extract_and_generate_mermaid.py   # mermaid -> PNG (script 1)
    ├── merge_bib.py                       # merge our refs into the main bib
    ├── README.md                          # this file
    └── uq-subsections/                    # *** UPLOAD COPY: only final files ***
        ├── pilot_uq_usecase.tex           #   -> Overleaf  chapters/uq-subsections/
        ├── uq-usecase_refs.bib            #   refs to merge into the main bib
        └── figs/uq-pilot/                 #   -> Overleaf  <root>/figs/uq-pilot/
            └── *.png                      #   9 figures
```

**2. Report root (local copy of the Overleaf project, for compiling/previewing).**
Extracted from the SciDAC Overleaf zip. We own **only** the
`chapters/uq-subsections/` folder inside it.
Root: `/Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/`

```
SciDAC_OE_Technical_Report/                # report root: main.tex, common/, figs/ (all from the zip)
└── chapters/
    └── uq-subsections/                    # <- the only folder we own here
        ├── pilot_uq_usecase.tex           # generated content (copy of the upload file)
        ├── uq-usecase_refs.bib            # refs (reference copy)
        ├── figs/uq-pilot/*.png            # generated figures
        ├── test_uq_pilot.tex              # faithful local test harness
        ├── build_test.sh                  # one-command build -> _build/test_uq_pilot.pdf
        └── _build/                        # compile artifacts (disposable)
```

Each `convert_pilot_report.py` run writes **both** copies: the repo upload copy
(version control / Overleaf) and the report copy (local compiling). Nothing of ours
is written to the report root — `build_test.sh` compiles *from* the root (so paths
resolve like `main.tex`) but sends all artifacts into `_build/`.

## Quick edit loop (minor text changes)

For a small wording or content change, edit the markdown and preview:

```bash
# 1. Edit the SOURCE markdown
#    /Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/datacenter_pilot_report.md

# 2. Regenerate the LaTeX (writes BOTH the upload copy and the report copy)
cd /Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/tex-report
python convert_pilot_report.py

# 3. Rebuild the preview PDF and open it
bash /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections/build_test.sh
open /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections/_build/test_uq_pilot.pdf
```

If you added or changed a **mermaid diagram**, run
`python extract_and_generate_mermaid.py` once before step 2 to regenerate its PNG.

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

  [3] fig_3_one_time_seed_building_blue_and_the_pe
      Caption: One-time seed building (blue) and the per-condition build (grey).
      → .../figures/fig_3_one_time_seed_building_blue_and_the_pe.png
      ✓ Generated from mermaid-cli
  [4] fig_4_per_condition_pipeline_from_hour_and_co
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
  ✓ Copied figure: fig_3_one_time_seed_building_blue_and_the_pe.png
✓ Generated §3 LaTeX (2935 bytes)
...
[SUCCESS] Phase 3 conversion complete!
```

### Step 3: Test PDF Compilation (Optional)

Use the **faithful harness** `test_uq_pilot.tex`. It loads the report's real
preamble (`common/header.tex`, i.e. `titlesec`, `natbib`, `glossaries`, …) and the
real chapter/section context, so it reproduces Overleaf's **section numbering and
heading formatting exactly** (this is what catches `titlesec` stretching and the
`9.1.x.y` numbering — a generic standalone preamble would not).

Build it with the helper script — one command, and every artifact stays in
`_build/` (nothing touches the report root):

```bash
bash /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections/build_test.sh
open /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections/_build/test_uq_pilot.pdf
```

Full paths (all under our folder `chapters/uq-subsections/`):

- Harness source: `.../chapters/uq-subsections/test_uq_pilot.tex`
- Build script: `.../chapters/uq-subsections/build_test.sh`
- **Output PDF (this is the one to inspect):** `.../chapters/uq-subsections/_build/test_uq_pilot.pdf`

The script `cd`s to the report root (so `common/`, `figs/`, `references.bib` resolve
like `main.tex`) and runs the full `pdflatex → bibtex → pdflatex → pdflatex` cycle
with `-output-directory=.../_build`. The harness wraps the generated content in
`\chapter`(9) → `\section`(9.1) so the pilot renders as `9.1.1`, `9.1.2`, …,
`9.1.4.2`, exactly as on Overleaf, and uses `\graphicspath` to find the figures
under `chapters/uq-subsections/figs/uq-pilot/`.

> **Numbering note:** every inner heading is emitted as a *numbered* command
> (`\subsection{…}`, not `\subsection*{…}`), and the file brackets its content with
> `\setcounter{secnumdepth}{3}` … `{2}` so subsubsections get `9.1.x.y` numbers
> *only within this block*. Manual numbers from the markdown (e.g. "§7.6") become
> `\label{uqp:7.6}`, and prose cross-references resolve through `\ref{…}`, so the
> printed number always matches regardless of where the block sits in the chapter.
> (Starred headings are avoided because `titlesec` both un-numbers them and stretches
> their titles to full width.)

## Configuration

Edit configuration in each script:

### extract_and_generate_mermaid.py (lines ~20-30)
```python
SOURCE_MD = "...datacenter_pilot_report.md"
FIGURES_OUTPUT_DIR = ".../work-notes/figures/"
OUTPUT_FORMAT = "png"  # or "tikz" (future)
```

### convert_pilot_report.py — `CONFIGURATION` block (top of file, ~lines 25-66)
```python
SOURCE_MD = "...datacenter_pilot_report.md"      # markdown you edit
FIGURES_SOURCE_DIR = ".../work-notes/figures/"   # read pre-generated PNGs
FIG_TEX_DIR = "figs/uq-pilot"                    # figure sub-path in \includegraphics
TARGET_PATH = ".../chapters/uq-subsections/pilot_uq_usecase.tex"        # report copy
REPO_OUTPUT_PATH = ".../tex-report/uq-subsections/pilot_uq_usecase.tex" # upload copy
INCLUDE_FOOTNOTES = True

# Sections to OMIT from the report (the markdown keeps them). Inbound "§N"
# references to a skipped section are cleaned up automatically.
SKIP_SECTIONS = {10, 11, 12}   # feedback, limitations, implementation

# "> **Details:** [...](plan.md#...)" link boxes: kept in the markdown, dropped here.
INCLUDE_DETAILS = False

# Prose terms rendered as \ac{...} (first use spelled out, hyperlinked).
ACRONYMS = [("DC OPF", r"\ac{dcopf}"), ("AC PF", r"\ac{acpf}"), ("CCT", r"\ac{cct}"), ...]
```

**Markdown-only text.** Wrap anything that should stay in the markdown but not
appear in the report in skip markers:

```markdown
<!-- tex:skip -->
Run identifiers in the plan: aleatoric campaign 1 = `aleatoric1_ap08_qlim`, ...
<!-- /tex:skip -->
```

**Acronyms.** Every key in `ACRONYMS` must exist in the report's
`common/acronyms.tex`. Ours that the report lacked (`dcopf`, `acpf`, `cct`, `ups`)
are in `chapters/uq-subsections/uq_pilot_acronyms.tex`; the local harness loads it,
and on Overleaf its four `\newacronym` lines go into `common/acronyms.tex`.
Headings and glossary labels use the plain short form (DCOPF, ACPF) instead.

**Skip / un-skip a section.** Edit `SKIP_SECTIONS` near the top of
`convert_pilot_report.py`, then re-run the converter. Examples:
`{10, 11, 12}` drops all three; `{10, 12}` brings Limitations back; `set()` keeps everything.

**Adjust figure sizes.** Mermaid diagram widths are picked automatically from each
PNG's aspect ratio in `_mermaid_width_for()` (search for that function):

```python
if ratio >= 2.0:  return 1.0               # wide/horizontal  -> full text width
if ratio >= 1.3:  return max(fallback, 0.95)
return fallback                            # portrait/square   -> caller's default
```

To change a size, either tune those thresholds, or change the per-figure default
width passed by the generators:
- mermaid diagrams (all sections): `_emit_mermaid_from_caption(...)` passes `0.72`;
  wide ones are then widened to full width and portrait ones capped at `0.6`.
- matplotlib plots (`main_*` KM figures): `generate_latex_figure_image(..., width=0.85)`.
  (These are not aspect-scaled — set their width directly.)

## Phases

Each run converts the requested section **and all prior sections**. The default
(no argument) is `phase=12`, which converts the whole report (§1 through §12).

- **Phase 1**: §1 "What this is" - Text with blockquotes, bold, italic
- **Phase 2**: §2 "Glossary" - Table → definition-list with footnotes
- **Phase 3**: §3 "The system" - Text with mermaid diagram (PNG) and details boxes
- **Phase 4**: §4 "Building one operating condition" (incl. §4.1–§4.4)
- **Phase 5–9**: §5–§9 via the generic robust block converter
  (tables, `$$` display math, images, fenced code/mermaid, headings, details boxes)
- **Phase 10–12**: §10 feedback, §11 limitations, §12 implementation

Citations: `[N]` markers in the markdown map through `CITATION_MAP` to
`\cite{key}`; keys are defined in `uq-usecase_refs.bib`. The generated
`pilot_uq_usecase.tex` is `\input` into a parent document that owns the
bibliography (the SciDAC report, or the faithful harness `test_uq_pilot.tex`).

## Features

### Markdown to LaTeX Conversion
- **Backticks** (`` `code` ``) → `\texttt{code}` (monospace); LaTeX specials escaped and long identifiers made breakable at underscores
- **Bold / italic** (`**b**`, `*i*`) → `\textbf{}` / `\textit{}`
- **Headings** (`##`/`###`/`####`) → starred `\subsection*`/`\subsubsection*`/`\paragraph*` with the **manual section number kept** (e.g. `7.0 …`), so the report shows our numbering, not LaTeX's
- **Blockquotes** → `\begin{quote}`; **Details** and **Implementation note** boxes → boxed `\fbox{\parbox{…}}`
- **Ordered / bullet lists** → `enumerate` / `itemize`
- **Links** (`[text](url)`) → plain text (URL dropped), specials in the text escaped
- **Inline / display math** (`$…$`, `$$…$$`) → protected from text formatting; Unicode mapped per math vs text context (°, ×, ·, −, ≤, ≥, Δ, …)
- **Citations** (`[N]`) → `\cite{key}` via `CITATION_MAP`

### Figure Handling

**Figure label generation** (consistent across scripts):
```
Caption: "One-time seed building (blue) and the per-condition build (grey)."
Label: "fig_3_one_time_seed_building_blue_and_the_pe"
```

- **Extraction**: `extract_and_generate_mermaid.py` finds all mermaid blocks and generates PNGs
- **Conversion**: `convert_pilot_report.py` matches captions to PNGs using same label logic
- **References**: LaTeX uses `\includegraphics{figs/uq-pilot/fig_3_...png}` — a path under the project-root `figs/` tree (underscore, not colon — colons break on Overleaf)

### Special Handling
- **Tables**: the glossary → a `description` list; every other table → `booktabs` `tabularx` (no vertical rules, content-sized first column, wrapping `X` columns)
- **Pipe handling**: Correctly handles `|` in backticks and math mode
- **Multiple links**: Details boxes support multiple section links
- **Flexible formatting**: Robust pattern handling for edge cases

## Design Principles

1. **Idempotent**: Each run regenerates all output from scratch
2. **Deterministic**: Same input always produces identical output
3. **Configuration-driven**: All paths/settings as module variables
4. **Separated concerns**: Diagram generation ≠ LaTeX conversion
5. **Preprocessing-first**: Fix markdown issues before scripts process
6. **Standard packages only**: `tabularx`, `booktabs`, `array`, `amsmath`/`amssymb`, `graphicx` — all already in the SciDAC report preamble

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

## Porting to Overleaf

The main SciDAC report lives on Overleaf (`report` document class, plain BibTeX,
`\bibliography{references.bib}`). The pilot is a self-contained chunk that is
`\input` into `chapters/uq.tex`, under §9.1 (Modeling Aleatoric Uncertainties).

**Numbering.** Every heading is a *numbered* command, so the pilot renders as
`9.1.1`, `9.1.2`, …, `9.1.7.3` with the report's own numbers. The markdown's manual
numbers become labels (`§7.6` → `\label{uqp:7.6}`), and in-text references are
printed as `Section~\ref{uqp:7.6}`, so they always match the printed heading.

### What to copy

Run the pipeline first so the clean copy is current:

```bash
cd /Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/tex-report
python extract_and_generate_mermaid.py   # only if diagrams changed
python convert_pilot_report.py           # default = whole report, §1–§12
```

Each run (re)writes a clean, self-contained copy to this folder:

```
/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/tex-report/uq-subsections/
├── pilot_uq_usecase.tex     # → Overleaf  chapters/uq-subsections/
├── figs/uq-pilot/           # → Overleaf  figs/uq-pilot/   (PROJECT ROOT, not chapters/)
└── uq-usecase_refs.bib      # the 26 references (merged into the main bib, step 3)
```

> **Why figures go to the project-root `figs/`, not `chapters/uq-subsections/figs/`:**
> `\includegraphics{figs/uq-pilot/...}` is resolved relative to the **main document**
> (the Overleaf project root where `main.tex` lives), *not* relative to the `\input`-ed
> file. The report already keeps its images in a top-level `figs/` tree (`figs/ignas/...`),
> so ours must live at `<project-root>/figs/uq-pilot/`. Putting them under
> `chapters/uq-subsections/figs/` is why they render as empty boxes. (Change `FIG_TEX_DIR`
> at the top of `convert_pilot_report.py` if you want a different sub-folder.)

### Steps

1. **Upload the content via the Overleaf web UI**, to **two** destinations (both
   live under
   `/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/tex-report/uq-subsections/`):

   - `pilot_uq_usecase.tex`  →  Overleaf `chapters/uq-subsections/`
   - the 9 PNGs in `figs/uq-pilot/`  →  Overleaf `figs/uq-pilot/` **at the project root**
     (make a `uq-pilot` folder inside the existing top-level `figs/`; do **not** put
     them under `chapters/`).

   (The `uq-usecase_refs.bib` next to them is only a reference copy — you merge it
   into the main bibliography in step 3, so don't upload it separately.)

2. **Include it in `uq.tex`.** One line, under `\section{Modeling Aleatoric Uncertainties}`:
   ```latex
   \input{chapters/uq-subsections/pilot_uq_usecase}
   ```
   The file opens with its own numbered `\subsection{Pilot Study: Datacenter Load
   and Transient Security}`, so it appears as `9.1.x` in the main ToC with no extra
   wrapper; everything inside keeps our manual numbers.

3. **Merge the bibliography.** In Overleaf, download the project's `references.bib`
   (right-click the file → Download) and save it over the copy in your local unzipped
   report folder, i.e. at
   `/Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/references.bib`.
   Then augment it with `merge_bib.py` — it appends only the entries whose key is
   **not** already present (nothing existing is touched), writes a `.bak` backup,
   and aborts if it would introduce a duplicate key:
   ```bash
   cd /Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/tex-report
   python merge_bib.py /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/references.bib --dry-run   # preview
   python merge_bib.py /Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/references.bib             # do it
   ```
   Typical output: `appended 25 entries; skipped (present): ['wachter2006']; duplicate keys (MUST be empty): []`.
   Then upload that updated `references.bib` back to Overleaf (overwrite the existing one).

   Why `wachter2006` is skipped: our Ipopt entry deliberately uses the report's
   existing key, so `\cite{wachter2006}` resolves to the entry already there — no
   duplicate. **SUNDIALS** (`hindmarsh2005sundials`) and **GridKit**
   (`peles2025gridkit`) exist under *different* keys than ours, so they are appended
   as a (harmless but redundant) second entry for the same work; repoint those via
   `CITATION_MAP` and reconvert if you prefer a single entry. **Never** put a
   `\bibliography` inside `pilot_uq_usecase.tex` — `main.tex` owns the single
   bibliography.

4. **Packages — nothing to add.** `common/header.tex` already loads `amsmath`,
   `amssymb`, `graphicx`, `booktabs`, `tabularx` (which pulls in `array`),
   `xcolor`, `placeins`, `lmodern`, `microtype`.

   **Acronyms — one-time.** Paste the four `\newacronym` lines from
   `chapters/uq-subsections/uq_pilot_acronyms.tex` (without the `\ifglsentryexists`
   guards) into Overleaf's `common/acronyms.tex`. Without them, `\ac{cct}` etc. fail.

5. **Compile.** Overleaf's latexmk runs the BibTeX passes automatically. Locally
   it is the usual four-pass cycle: `pdflatex → bibtex → pdflatex → pdflatex`.

### Re-porting after edits

Edit the markdown, re-run `convert_pilot_report.py`, then repeat step 1 (and step 3
only if citations changed). The output is deterministic, so re-uploading the
folder is safe.

## Dependencies

- **Python 3** (standard library only: `re`, `sys`, `pathlib`, `subprocess`)
- **pdflatex + bibtex** (TeX Live 2026)
- **mermaid-cli** (for diagram PNGs: `npm install -g @mermaid-js/mermaid-cli`)
- **Pillow** (optional, for placeholder PNGs if mermaid-cli is missing)

