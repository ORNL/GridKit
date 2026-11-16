#!/usr/bin/env python3
"""
Convert datacenter_pilot_report.md (Markdown) → pilot_uq_usecase.tex (LaTeX)

This script incrementally converts sections of the markdown pilot UQ report
into LaTeX format for inclusion in the SciDAC technical report.

Phase 1 (MVP): Converts §1 "What this is" (text only, no figures/tables).
Phases 2-6: Add tables (§2), figures/Mermaid (§3-4.4), math/footnotes.
Phases 7-13: Add remaining sections (§5-12) with text/math/tables.

Design: Idempotent, deterministic conversion. Each run parses the full markdown
and regenerates the LaTeX output from scratch. No state files or checksums.
Handles all incremental markdown changes automatically.
"""

import re
import sys
from pathlib import Path

# ============================================================================
# CONFIGURATION — Adjust these paths as needed
# ============================================================================

SOURCE_MD = "/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/datacenter_pilot_report.md"
FIGURES_SOURCE_DIR = (
    "/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/figures/"
)

# Figures are referenced as \includegraphics{figs/uq-pilot/<name>}. LaTeX resolves
# that path from the MAIN document's directory (the Overleaf PROJECT ROOT, where
# main.tex lives), NOT from the \input'd file. So on Overleaf the PNGs must be
# uploaded to <project-root>/figs/uq-pilot/ — the same top-level figs/ tree the
# report already uses (figs/ignas/...), not chapters/uq-subsections/figs/.
# Change FIG_TEX_DIR if you want a different sub-folder (e.g. "figs/ignas/uq-pilot").
FIG_TEX_DIR = "figs/uq-pilot"

# PRIMARY OUTPUT: SciDAC technical report. The faithful harness (test_uq_pilot.tex,
# built via chapters/uq-subsections/build_test.sh) finds these figures through
# \graphicspath{{chapters/uq-subsections/}}, so they must sit at
# chapters/uq-subsections/figs/uq-pilot/.
TARGET_PATH = "/Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections/pilot_uq_usecase.tex"
FIGURES_OUTPUT_DIR = (
    "/Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections/"
    + FIG_TEX_DIR
    + "/"
)

# SECONDARY OUTPUT: Keep generated LaTeX in repo for version control
REPO_OUTPUT_PATH = "/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/tex-report/uq-subsections/pilot_uq_usecase.tex"
REPO_FIGURES_OUTPUT_DIR = (
    "/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/tex-report/uq-subsections/"
    + FIG_TEX_DIR
    + "/"
)
REPO_BIBTEX_PATH = "/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/tex-report/uq-subsections/uq-usecase_refs.bib"

# Phase 2+ options
INCLUDE_FOOTNOTES = True  # Set to False to omit cross-document references as footnotes
REPO_BRANCH = "isatkaus/uq-usecase"  # Git branch where details are located
REPO_BASE = "https://github.com/ORNL/GridKit"  # GitHub repository base URL

# Sections to omit from the report (the markdown keeps them; they are just not
# emitted here). Inbound "§N" cross-references to a skipped section are cleaned up
# automatically by apply_numbered_headings_and_refs so nothing dangles.
SKIP_SECTIONS = {10, 11}  # 10 = GridKit end-user feedback, 11 = Limitations

# Phase 3+ options (Mermaid diagrams)
MERMAID_FORMAT = "png"  # Options: "tikz" (inline, styled), "png" (requires mermaid-cli, creates placeholder if unavailable)

# Inline code formatting options
# INLINE_CODE_FORMAT = "math"   # `code` → $code$ (inline math mode, current behavior)
# INLINE_CODE_FORMAT = "texttt"  # `code` → \texttt{code} (monospace font, simpler, preserves spacing)
INLINE_CODE_FORMAT = "texttt"  # Choose preferred rendering: "math" or "texttt"

# Sections to exclude from conversion
# Useful for optional subsections that may not be needed in final report
# Example: EXCLUDE_SECTIONS = ["4.2"]  # Skip a section (4.2a and 4.2b are now paragraphs within 4.2)
EXCLUDE_SECTIONS = []  # Currently empty; set to ["4.2"] or other sections to exclude

# ============================================================================
# CITATION MAPPING — Maps numeric citations [1], [2], etc. to BibTeX keys
# ============================================================================
# Usage in markdown: (Author et al., YEAR [1], details)
# Converts to LaTeX: (Author et al., YEAR \cite{BibtexKey}, details)

CITATION_MAP = {
    "1": "Elliott2021",  # Elliott et al., 2021 - OSTI 1819998 (Advanced Control and Monitoring of IBRs)
    "2": "PolitisRomano1994",  # Politis and Romano, 1994 - Subsampling under minimal assumptions
    "3": "GridKit2025",  # GridKit v0.2.0 - the simulation framework used
    "4": "Birchfield2017",  # ACTIVSg synthetic grids (validation criteria)
    "5": "TAMUrepo",  # Texas A&M Electric Grid Test Case Repository
    "6": "Coffrin2018",  # PowerModels.jl
    "7": "wachter2006",  # Ipopt interior-point NLP solver (matches main report key)
    "8": "Hindmarsh2005",  # SUNDIALS / IDA DAE integrator
    "9": "Yoo2003",  # SLURM resource manager
    "10": "Bezanson2017",  # Julia language
    "11": "DerKiureghian2009",  # aleatory vs epistemic uncertainty
    "12": "Milano2018",  # low-inertia systems (renewables reduce inertia)
    "13": "IEEELoad1993",  # ZIP load model (IEEE Task Force)
    "14": "Milanovic2013",  # international industry load-modeling practice
    "15": "Pourbeik2017",  # generic dynamic models for IBRs
    "16": "Shehabi2024",  # datacenter energy / grid impact
    "17": "Kundur2004",  # definition and classification of power system stability
    "18": "Kundur1994",  # Power System Stability and Control (textbook)
    "19": "SauerPai1998",  # Power System Dynamics and Stability (COI)
    "20": "Stott1974",  # review of load-flow methods (PV/PQ switching)
    "21": "KaplanMeier1958",  # Kaplan-Meier estimator
    "22": "Greenwood1926",  # Greenwood's variance formula
    "23": "Giles2015",  # multilevel Monte Carlo
    "24": "Battaglia2018",  # graph networks / relational inductive bias
    "25": "Donon2020",  # GNN surrogate for power flow
    "26": "Settles2009",  # active learning survey
}

# ============================================================================
# UNICODE HANDLING — text-mode and math-mode replacements
# ============================================================================
# The report uses Unicode symbols that pdflatex cannot typeset directly.
# Inside math ($...$ and $$...$$) they map to math commands; in prose they map
# to text-mode equivalents (wrapped in $...$ where only a math glyph exists).

# Applied INSIDE math spans (content is already LaTeX math).
MATH_UNICODE = {
    "°": r"^{\circ}",
    "×": r"\times ",
    "·": r"\cdot ",
    "−": "-",
    "∼": r"\sim ",
    "≤": r"\le ",
    "≥": r"\ge ",
    "≈": r"\approx ",
    "≠": r"\ne ",
    "±": r"\pm ",
    "²": r"^{2}",
    "³": r"^{3}",
    "…": r"\dots ",
    "→": r"\to ",
    "∞": r"\infty ",
    "∈": r"\in ",
    "≡": r"\equiv ",
    "∝": r"\propto ",
    "–": "-",
    "—": "---",
    "Δ": r"\Delta ",
    "δ": r"\delta ",
    "α": r"\alpha ",
    "β": r"\beta ",
    "σ": r"\sigma ",
    "μ": r"\mu ",
    "θ": r"\theta ",
    "λ": r"\lambda ",
    "π": r"\pi ",
    "ρ": r"\rho ",
    "ω": r"\omega ",
    "γ": r"\gamma ",
    "Σ": r"\Sigma ",
}

# Applied in text (prose) mode. Order does not matter (all keys are distinct).
TEXT_UNICODE = [
    ("§", r"{\textsection}"),
    ("•", r"\textbullet{}"),
    ("–", "--"),
    ("—", "---"),
    ("→", r"$\to$"),
    ("←", r"$\leftarrow$"),
    ("↔", r"$\leftrightarrow$"),
    ("×", r"$\times$"),
    ("·", r"$\cdot$"),
    ("−", r"$-$"),
    ("∼", r"$\sim$"),
    ("≈", r"$\approx$"),
    ("≤", r"$\le$"),
    ("≥", r"$\ge$"),
    ("≠", r"$\ne$"),
    ("±", r"$\pm$"),
    ("°", r"$^{\circ}$"),
    ("²", r"$^{2}$"),
    ("³", r"$^{3}$"),
    ("…", r"\dots "),
    ("∞", r"$\infty$"),
    ("∈", r"$\in$"),
    ("≡", r"$\equiv$"),
    ("∝", r"$\propto$"),
    ("✓", r"$\checkmark$"),
    ("“", "``"),
    ("”", "''"),
    ("‘", "`"),
    ("’", "'"),
    ("Δ", r"$\Delta$"),
    ("α", r"$\alpha$"),
    ("β", r"$\beta$"),
    ("δ", r"$\delta$"),
    ("σ", r"$\sigma$"),
    ("μ", r"$\mu$"),
    ("θ", r"$\theta$"),
    ("λ", r"$\lambda$"),
    ("π", r"$\pi$"),
    ("ρ", r"$\rho$"),
    ("ω", r"$\omega$"),
    ("γ", r"$\gamma$"),
    ("Σ", r"$\Sigma$"),
]


def latexify_math(s):
    """Fix Unicode symbols inside a LaTeX math span (content between $...$)."""
    for uni, rep in MATH_UNICODE.items():
        s = s.replace(uni, rep)
    return s


def latex_text_safe(s):
    """Sanitize a short plain-text string for LaTeX text mode (e.g. link labels).

    Drops backticks, maps Unicode, and escapes LaTeX specials _ & # %.
    """
    s = s.replace("`", "")
    for uni, rep in TEXT_UNICODE:
        s = s.replace(uni, rep)
    s = re.sub(r"(?<!\\)([_&#%])", r"\\\1", s)
    return s


# ============================================================================
# PHASE 1: Extract §1 and convert to LaTeX
# ============================================================================


def extract_section_1(md_content):
    """
    Extract section §1 "What this is" from markdown.
    Assumes section 1 starts with "## 1. What this is" and ends before "## 2. ".
    Returns the raw markdown text of §1 (without the ## heading).
    """
    # Find the start of §1
    pattern_start = r"^## 1\.\s+What this is\s*$"
    match_start = re.search(pattern_start, md_content, re.MULTILINE)
    if not match_start:
        raise ValueError("Could not find '## 1. What this is' section")

    start_pos = match_start.end()

    # Find the start of §2 (section ends at next ## heading)
    pattern_end = r"^## 2\."
    match_end = re.search(pattern_end, md_content[start_pos:], re.MULTILINE)

    if match_end:
        end_pos = start_pos + match_end.start()
    else:
        # If §2 not found, take everything remaining (shouldn't happen)
        end_pos = len(md_content)

    section_1_text = md_content[start_pos:end_pos].strip()
    return section_1_text


def markdown_to_latex_text(text):
    r"""
    Convert basic markdown text formatting to LaTeX.
    Phase 1 handles:
    - Blockquotes (lines starting with >) → \begin{quote}...\end{quote}
    - Bullet lists (lines starting with * or -) → \begin{itemize}...\end{itemize}
      * Continuation lines indented with spaces are joined to the item
    - `code` → formatted per INLINE_CODE_FORMAT setting (math mode or texttt font)
    - **bold** → \textbf{bold}
    - *italic* → \textit{italic}
    - Markdown links [text](url) → \url{url} with link text
    - Subsections (### text) → \subsubsection{text}
    - Preserves paragraphs (double newlines)
    - Unicode characters escaped to LaTeX
    """
    # Handle blockquotes and bullet lists
    lines = text.split("\n")
    result_lines = []
    in_quote = False
    in_list = False
    quote_buffer = []
    list_buffer = []
    list_kind = None  # "itemize" or "enumerate"
    quote_started = False  # Track if we've output \begin{quote}

    def close_quote():
        nonlocal in_quote, quote_buffer, quote_started
        if in_quote:
            processed_quote = "\n".join(quote_buffer)

            # Boxed note: Details and Implementation-note blockquotes use the fbox
            if (
                "Details:" in processed_quote
                or "Implementation note:" in processed_quote
            ) and INCLUDE_FOOTNOTES:
                # Format as Details blockquote with fbox (no quote tags)
                # Close quote tag if it was started (for mixed blockquotes)
                if quote_started:
                    result_lines.append("\\end{quote}")

                blockquote_md = "> " + processed_quote  # Reconstruct markdown format
                details_latex = format_details_blockquote(blockquote_md)
                if details_latex:
                    result_lines.append(r"\par\vspace{0.4em}")
                    result_lines.append(r"\noindent\fbox{\parbox[t]{0.98\textwidth}{%")
                    result_lines.append(r"\setlength{\parskip}{0pt}%")
                    result_lines.append(details_latex)
                    result_lines.append(r"}}")
                    result_lines.append(r"\vspace{0.2em}")
            else:
                # Regular blockquote
                if not quote_started:
                    result_lines.append("\\begin{quote}")
                    quote_started = True

                processed_quote = markdown_inline_formatting(processed_quote)

                # Convert markdown links in blockquote
                def link_to_url(m):
                    link_text = m.group(1)
                    link_url = m.group(2)
                    return "\\url{" + link_url + "} (" + link_text + ")"

                processed_quote = re.sub(
                    r"\[([^\]]+)\]\(([^)]+)\)", link_to_url, processed_quote
                )
                result_lines.extend(processed_quote.split("\n"))
                result_lines.append("\\end{quote}")
                quote_started = False

            quote_buffer = []
            in_quote = False

    def close_list():
        nonlocal in_list, list_buffer, list_kind
        if in_list:
            env = "enumerate" if list_kind == "enumerate" else "itemize"
            result_lines.append("\\begin{" + env + "}")
            for list_item in list_buffer:
                # Join multi-line items with spaces and apply formatting
                formatted_item = markdown_inline_formatting(list_item)
                result_lines.append(f"  \\item {formatted_item}")
            result_lines.append("\\end{" + env + "}")
            list_buffer = []
            in_list = False
            list_kind = None

    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Check for blockquote line
        if stripped.startswith("> "):
            close_list()
            quote_line = stripped[2:].strip()
            quote_buffer.append(quote_line)
            if not in_quote:
                in_quote = True
                # Don't add \begin{quote} yet - wait until close_quote() to decide
            i += 1

        # Check for bullet list line (starts with * or -)
        elif (
            stripped
            and len(stripped) > 1
            and stripped[0] in ("*", "-")
            and stripped[1] == " "
        ):
            close_quote()
            if in_list and list_kind != "itemize":
                close_list()
            # Extract bullet item text
            item_text = stripped[2:].strip()

            # Look ahead for continuation lines (lines starting with whitespace)
            i += 1
            while i < len(lines):
                next_line = lines[i]
                next_stripped = next_line.strip()

                # Stop if we hit a blank line or new item
                if not next_stripped:
                    break
                if (
                    next_stripped
                    and len(next_stripped) > 1
                    and next_stripped[0] in ("*", "-")
                    and next_stripped[1] == " "
                ):
                    break
                if next_line and next_line[0] not in (" ", "\t"):  # Not indented
                    break

                # This is a continuation line - join it to the item
                item_text += " " + next_stripped
                i += 1

            list_buffer.append(item_text)
            if not in_list:
                in_list = True
                list_kind = "itemize"
            # i is already incremented by while loop
            continue

        # Check for ordered list line (starts with "N. ")
        elif re.match(r"^\d+\.\s", stripped):
            close_quote()
            if in_list and list_kind != "enumerate":
                close_list()
            item_text = re.sub(r"^\d+\.\s+", "", stripped)

            # Look ahead for continuation lines (indented, non-item lines)
            i += 1
            while i < len(lines):
                next_line = lines[i]
                next_stripped = next_line.strip()
                if not next_stripped:
                    break
                if re.match(r"^\d+\.\s", next_stripped):
                    break
                if (
                    len(next_stripped) > 1
                    and next_stripped[0] in ("*", "-")
                    and next_stripped[1] == " "
                ):
                    break
                if next_line and next_line[0] not in (" ", "\t"):  # Not indented
                    break
                item_text += " " + next_stripped
                i += 1

            list_buffer.append(item_text)
            if not in_list:
                in_list = True
                list_kind = "enumerate"
            continue

        else:
            # Regular text or empty line
            close_quote()
            close_list()

            # Handle headings: ## / ### / #### → section levels (strip manual numbers)
            heading_match = re.match(r"^(#{2,4})\s+(.*)$", stripped)
            if heading_match:
                hashes = heading_match.group(1)
                # Keep the manual section number (e.g. "7.0") so "§7.0" references
                # in the prose resolve to the heading. Use starred (unnumbered)
                # commands so LaTeX does not add its own conflicting number.
                htext = markdown_inline_formatting(heading_match.group(2).strip())
                if hashes == "##":
                    result_lines.append(f"\\subsection*{{{htext}}}")
                    result_lines.append(
                        f"\\addcontentsline{{toc}}{{subsection}}{{{htext}}}"
                    )
                elif hashes == "###":
                    result_lines.append(f"\\subsubsection*{{{htext}}}")
                    result_lines.append(
                        f"\\addcontentsline{{toc}}{{subsubsection}}{{{htext}}}"
                    )
                else:
                    result_lines.append(f"\\paragraph*{{{htext}}}")
            elif stripped == "---":
                pass  # Markdown horizontal rule → ignore
            elif stripped:  # Non-empty line
                # Collect consecutive non-empty lines as a paragraph
                # This allows bold/italic formatting to work across line breaks
                paragraph_lines = [line]
                i += 1
                while i < len(lines):
                    next_line = lines[i]
                    next_stripped = next_line.strip()

                    # Stop if we hit an empty line, blockquote, or list item
                    if not next_stripped:
                        break
                    if next_stripped.startswith(">"):
                        break
                    if (
                        next_stripped
                        and len(next_stripped) > 1
                        and next_stripped[0] in ("*", "-")
                        and next_stripped[1] == " "
                    ):
                        break
                    if re.match(r"^\d+\.\s", next_stripped):
                        break
                    if re.match(r"^#{2,4}\s", next_stripped):
                        break
                    if next_stripped == "---":
                        break

                    # This is part of the same paragraph
                    paragraph_lines.append(next_line)
                    i += 1

                # Join paragraph lines with spaces and apply formatting
                paragraph_text = " ".join(
                    line.strip() for line in paragraph_lines if line.strip()
                )
                formatted_para = markdown_inline_formatting(paragraph_text)
                result_lines.append(formatted_para)
                continue  # Skip the normal i += 1 at the end
            elif (
                result_lines and result_lines[-1]
            ):  # Preserve blank lines between paragraphs
                result_lines.append(line)

            i += 1

    # Close any open structures at end of text
    close_quote()
    close_list()

    text = "\n".join(result_lines)
    return text


def generate_latex_section_1(md_section_1):
    """
    Convert extracted §1 markdown text to LaTeX.
    Wraps in \subsection and applies markdown→LaTeX conversions.
    """
    # Convert markdown formatting
    latex_text = markdown_to_latex_text(md_section_1)

    # Numbered subsection so the pilot gets a report-numbered ToC entry (9.1.x)
    # when \input into uq.tex. All inner headings are likewise numbered (see the
    # apply_numbered_headings_and_refs post-processor); uqp:1 lets prose "§1"
    # references resolve to this heading's auto-number.
    latex_output = r"""\subsection{Pilot Study: Datacenter Load and Transient Security}
\label{sec:pilot_uq_usecase}
\label{uqp:1}

"""
    latex_output += latex_text
    latex_output += "\n"

    return latex_output


# ============================================================================
# PHASE 2: Extract §2 (Glossary table) and convert to LaTeX
# ============================================================================


def extract_section_2(md_content):
    """
    Extract section §2 "Glossary" from markdown.
    Assumes section 2 starts with "## 2. Glossary" and ends before "---" or "## 3. ".
    Returns tuple: (table_markdown, footnote_markdown)
    """
    pattern_start = r"^## 2\.\s+Glossary\s*$"
    match_start = re.search(pattern_start, md_content, re.MULTILINE)
    if not match_start:
        raise ValueError("Could not find '## 2. Glossary' section")

    start_pos = match_start.end()

    # Find the end of §2 (before "---" or "## 3.")
    pattern_end = r"^(?:---|## 3\.)"
    match_end = re.search(pattern_end, md_content[start_pos:], re.MULTILINE)

    if match_end:
        end_pos = start_pos + match_end.start()
    else:
        end_pos = len(md_content)

    section_2_text = md_content[start_pos:end_pos].strip()

    # Split into table and footnotes (blockquotes)
    # Table comes first, footnotes after
    table_lines = []
    footnote_lines = []
    in_footnote = False

    for line in section_2_text.split("\n"):
        if line.strip().startswith("> "):
            in_footnote = True
        if in_footnote:
            footnote_lines.append(line)
        else:
            table_lines.append(line)

    return "\n".join(table_lines).strip(), "\n".join(footnote_lines).strip()


def split_on_pipe_outside_math(text):
    """
    Split text on | character, but ignore | inside math mode ($...$) or backticks (`...`).
    Preserves pipes inside $ or ` delimiters as part of the content.

    Returns list of cell contents (including outer | characters).
    """
    cells = []
    current_cell = ""
    in_math = False
    in_code = False

    for char in text:
        if char == "$":
            in_math = not in_math
            current_cell += char
        elif char == "`":
            in_code = not in_code
            current_cell += char
        elif char == "|" and not in_math and not in_code:
            cells.append(current_cell)
            current_cell = ""
        else:
            current_cell += char

    if current_cell or text.endswith("|"):
        cells.append(current_cell)

    return cells


def parse_markdown_table(table_md):
    """
    Parse a markdown table and return list of (header, rows).
    Handles wrapped rows (continuation lines that don't start with |).
    Expects format:
    | col1 | col2 |
    |---|---|
    | val1 | val2 |
    """
    # First, reconstruct wrapped rows by joining continuation lines
    all_lines = table_md.split("\n")
    reconstructed_lines = []
    current_line = ""

    for line in all_lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith("|"):
            if current_line:
                reconstructed_lines.append(current_line)
            current_line = line
        else:
            # Continuation of previous line
            if current_line:
                current_line += " " + line

    if current_line:
        reconstructed_lines.append(current_line)

    if len(reconstructed_lines) < 3:
        raise ValueError(
            "Invalid markdown table (need header, separator, at least 1 row)"
        )

    # Parse header
    header = [
        cell.strip()
        for cell in split_on_pipe_outside_math(reconstructed_lines[0])[1:-1]
    ]

    # Skip separator line (line 1)
    # Parse rows (from line 2 onwards)
    rows = []
    for line in reconstructed_lines[2:]:
        row = [cell.strip() for cell in split_on_pipe_outside_math(line)[1:-1]]
        if len(row) == len(header):
            rows.append(row)

    return header, rows


def markdown_glossary_to_latex_deflist(table_md):
    """
    Convert markdown glossary table to LaTeX definition list format.
    Uses bold terms followed by definitions as flowing text (no table).
    Much better for long glossaries that span multiple pages.
    """
    header, rows = parse_markdown_table(table_md)

    latex = r"\begin{description}" + "\n"

    for row in rows:
        if len(row) >= 2:
            term = row[0].strip()
            definition = row[1].strip()

            # Convert markdown formatting in both term and definition
            term = markdown_inline_formatting(term)
            definition = markdown_inline_formatting(definition)

            # Use \item[\textbf{term}] for definition list format
            latex += r"  \item[" + term + r"] " + definition + "\n"

    latex += r"\end{description}" + "\n"

    return latex


def markdown_table_to_latex_tabular(table_md, float=True, compact=False):
    """
    Convert markdown table to LaTeX tabular environment.

    Args:
        table_md: Markdown table string
        float: If True, wrap in table environment (can float). If False, use bare tabularx (inline).
               Use float=False when table should stay exactly where it appears in text.
        compact: If True, use simpler format without borders (better for small 2-3 row tables).
    """
    header, rows = parse_markdown_table(table_md)

    # For small inline tables (like 2-row commitment table), use compact format
    if compact:
        # Simple compact format without borders or fancy styling.
        # The centering/size changes are kept inside a group ({...\par}) so
        # \parfillskip set by \centering cannot leak out and make titlesec
        # stretch every following heading to full width.
        latex = r"\vspace{0.15cm}" + "\n"
        latex += r"{\centering\small" + "\n"
        latex += r"\setlength{\tabcolsep}{8pt}" + "\n"
        latex += r"\begin{tabular}{lll}" + "\n"  # 3 columns, left-aligned, no borders

        # Add header
        header_cells = []
        for cell in header:
            cell = markdown_inline_formatting(cell)
            header_cells.append(cell)
        latex += " & ".join(header_cells) + r" \\" + "\n"
        latex += r"\midrule" + "\n"  # Thin line under header

        # Add rows
        for row in rows:
            row_cells = []
            for cell in row:
                cell = markdown_inline_formatting(cell)
                row_cells.append(cell)
            latex += " & ".join(row_cells) + r" \\" + "\n"

        latex += r"\end{tabular}" + "\n"
        latex += r"\par}" + "\n"
        latex += r"\vspace{0.15cm}" + "\n"
        return latex

    # Clean, borderless layout (matching the compact "commitment" tables):
    # a content-sized first column so long identifiers never overlap the next
    # column, and ragged, wrapping X columns for the rest so wide cells stay
    # within \textwidth instead of overflowing. If the first column holds long
    # text (not a short label), make it a wrapping X column too.
    n_cols = len(header)
    x_col = r">{\raggedright\arraybackslash}X"
    col0_cells = [header[0]] + [r[0] for r in rows if r]
    col0_max = max(
        (len(c.replace("*", "").replace("`", "")) for c in col0_cells), default=0
    )
    if n_cols <= 1:
        col_spec = x_col
    elif col0_max <= 24:
        col_spec = "l " + " ".join([x_col] * (n_cols - 1))
    else:
        col_spec = " ".join([x_col] * n_cols)

    # Start with table wrapper (floating) or direct tabularx (inline)
    if float:
        latex = r"\begin{table}[htb]" + "\n"
        latex += r"  \centering\small" + "\n"
        latex += r"  \renewcommand{\arraystretch}{1.25}" + "\n"
        latex += r"  \begin{tabularx}{\textwidth}{" + col_spec + "}" + "\n"
    else:
        # Inline table: no float wrapper, just tabularx with spacing. Group the
        # centering/size so \parfillskip from \centering does not leak (which
        # would make titlesec stretch subsequent headings to full width).
        latex = r"\vspace{0.25cm}" + "\n"
        latex += r"{\centering\small" + "\n"
        latex += r"\renewcommand{\arraystretch}{1.25}" + "\n"
        latex += r"\begin{tabularx}{\textwidth}{" + col_spec + "}" + "\n"

    latex += r"\toprule" + "\n"

    # Add header
    header_cells = [markdown_inline_formatting(cell) for cell in header]
    latex += "    " + " & ".join(header_cells) + r" \\" + "\n"
    latex += r"\midrule" + "\n"

    # Add rows
    for row in rows:
        row_cells = [markdown_inline_formatting(cell) for cell in row]
        latex += "    " + " & ".join(row_cells) + r" \\" + "\n"

    latex += r"\bottomrule" + "\n"
    latex += r"\end{tabularx}" + "\n"

    if float:
        latex += r"\end{table}" + "\n"
    else:
        latex += r"\par}" + "\n"
        latex += r"\vspace{0.25cm}" + "\n"

    return latex


def markdown_inline_formatting(text):
    """
    Convert inline markdown formatting: **bold**, *italic*, `code`, $math$.
    Escapes Unicode characters FIRST, then backticks, then bold/italic.
    This prevents nested formatting issues (e.g., bold with code inside).
    Used within table cells.

    INLINE_CODE_FORMAT controls backtick rendering:
    - "math": `code` → $code$ (inline math mode, preserves operators but can lose spacing)
    - "texttt": `code` → \texttt{code} (monospace font, preserves spacing, simpler rendering)
    """
    # Protect inline math $...$ from text formatting; fix Unicode inside it
    math_spans = []

    def _stash_math(m):
        math_spans.append(latexify_math(m.group(1)))
        return f"\x00M{len(math_spans) - 1}\x00"

    text = re.sub(r"\$([^$]+?)\$", _stash_math, text)

    # `code` → format based on INLINE_CODE_FORMAT setting
    def convert_code(m):
        code = m.group(1)

        # Normalize Unicode symbols that appear inside code to ASCII equivalents
        for u, r in (
            ("−", "-"),
            ("–", "-"),
            ("—", "--"),
            ("∼", "~"),
            ("×", "x"),
            ("·", "."),
            ("≤", "<="),
            ("≥", ">="),
            ("°", "deg"),
        ):
            code = code.replace(u, r)

        if INLINE_CODE_FORMAT == "math":
            # Math mode: $code$
            code = code.replace("\\|", "|")  # Remove escapes from pipes
            code = code.replace("_", "\\_")  # Escape underscores in math mode
            return "$" + code + "$"
        else:  # texttt mode
            # Monospace font: \texttt{code}
            # Remove markdown escape sequences (\| → |)
            code = code.replace("\\|", "|")  # Remove pipe escapes
            code = code.replace("\\{", "{")  # Remove brace escapes
            code = code.replace("\\}", "}")  # Remove brace escapes
            # Escape LaTeX special characters in texttt mode
            code = code.replace(
                "\\", r"\textbackslash{}"
            )  # Backslash → \textbackslash{}
            code = code.replace("&", r"\&")  # Ampersand
            code = code.replace("#", r"\#")  # Hash
            code = code.replace("$", r"\$")  # Dollar
            code = code.replace("%", r"\%")  # Percent
            code = code.replace(
                "_", r"\_"
            )  # Underscore (still needs escaping in texttt)
            code = code.replace("~", r"\textasciitilde{}")  # Tilde
            code = code.replace("^", r"\textasciicircum{}")  # Caret
            # Allow long identifiers to wrap after underscores (no overflow/overlap)
            code = code.replace(r"\_", r"\_\allowbreak{}")
            return r"\texttt{" + code + "}"

    # Stash code spans so text-mode substitutions below do not touch them
    code_spans = []

    def _stash_code(m):
        code_spans.append(convert_code(m))
        return f"\x00C{len(code_spans) - 1}\x00"

    text = re.sub(r"`([^`]+?)`", _stash_code, text)

    # Escape Unicode and angle brackets for LaTeX (text mode)
    for uni, rep in TEXT_UNICODE:
        text = text.replace(uni, rep)
    text = text.replace("<", r"\textless{}").replace(">", r"\textgreater{}")
    # Escape remaining LaTeX specials in prose (code and math are stashed above)
    text = text.replace("%", r"\%").replace("&", r"\&").replace("#", r"\#")

    # [text](url) → text (markdown links to LaTeX)
    # Strip URL and keep only link text. Escape LaTeX specials the link text may
    # carry (e.g. underscores in a file name like dc_campaign_pipeline.md).
    #
    # Section references ("§N") are resolved later by apply_numbered_headings_and_refs,
    # which turns "{\textsection}N" into "\ref{uqp:N}". That is only correct for
    # *internal* references (bare "§N", or links whose target is a same-document
    # "#anchor"). A link to another file (e.g. "ida_ic_issue.md#...") that happens to
    # contain "§3a" points at a DIFFERENT document's section, so we render its "§"
    # as "\S{}" (same glyph) to hide it from the ref pass and keep it literal.
    def _convert_link(m):
        link_text = re.sub(r"(?<!\\)([_&#%])", r"\\\1", m.group(1))
        target = m.group(2).strip()
        if not target.startswith("#"):
            link_text = link_text.replace(r"{\textsection}", r"\S{}")
        return link_text

    text = re.sub(r"\[([^\]]+?)\]\(([^)]+?)\)", _convert_link, text)

    # [1], [2], etc. → \cite{BibtexKey} (numeric citations to LaTeX)
    # Convert numeric citations like "[1]" to "\cite{Elliott2021}"
    text = re.sub(
        r"\[(\d+)\]",
        lambda m: r"\cite{" + CITATION_MAP.get(m.group(1), f"cite{m.group(1)}") + "}",
        text,
    )

    # **bold** → \textbf{bold}
    # Use DOTALL flag to match across newlines and normalize whitespace
    text = re.sub(
        r"\*\*(.+?)\*\*",
        lambda m: r"\textbf{" + m.group(1).replace("\n", " ") + "}",
        text,
        flags=re.DOTALL,
    )

    # *italic* → \textit{italic}
    # Use DOTALL flag to match across newlines and normalize whitespace
    # Also escape underscores inside italic text so they render correctly
    def convert_italic(m):
        text_content = m.group(1).replace("\n", " ")
        # Escape underscores in italic text
        text_content = text_content.replace("_", r"\_")
        return r"\textit{" + text_content + "}"

    text = re.sub(
        r"(?<!\*)\*(.+?)\*(?!\*)",
        convert_italic,
        text,
        flags=re.DOTALL,
    )

    # Restore protected code and inline math
    text = re.sub(r"\x00C(\d+)\x00", lambda m: code_spans[int(m.group(1))], text)
    text = re.sub(
        r"\x00M(\d+)\x00", lambda m: "$" + math_spans[int(m.group(1))] + "$", text
    )

    return text


def extract_footnotes_from_blockquote(footnote_md):
    """
    Extract markdown links from blockquotes like:
    > **Details:** [physical intuition](datacenter_load_uq_plan.md#physical-intuition) covering ...

    Converts to plain text GitHub URL that users can easily navigate to.
    Returns footnote text or None if no link found.
    """
    if not footnote_md or not INCLUDE_FOOTNOTES:
        return None

    # Remove blockquote markers and flatten into single line
    lines = [
        l.strip()[2:].strip() if l.strip().startswith("> ") else l.strip()
        for l in footnote_md.split("\n")
        if l.strip()
    ]
    text = " ".join(lines)

    # Find markdown link: [link text](url#anchor)
    match = re.search(r"\[([^\]]+)\]\(([^)]+)\)", text)
    if not match:
        return None

    link_text = match.group(1)
    link_url = match.group(2)

    # Extract filename and anchor
    filename = link_url.split("#")[0] if "#" in link_url else link_url
    section_anchor = link_url.split("#")[-1] if "#" in link_url else ""

    # Convert anchor to readable format: physical-intuition → Physical Intuition
    if section_anchor:
        section_name = section_anchor.replace("-", " ").title()
    else:
        section_name = link_text

    # Get the rest of the text after the link
    remaining = text[match.end() :].lstrip()
    # Remove leading words like "covering"
    remaining = re.sub(
        r"^(covering|discussing|see|see also)\s+",
        "",
        remaining,
        flags=re.IGNORECASE,
    )

    # Build full GitHub URL as plain text (easy to copy/paste)
    # Path: gridkit/uq-usecase/work-notes/datacenter_load_uq_plan.md#physical-intuition
    github_url = f"{REPO_BASE}/blob/{REPO_BRANCH}/uq-usecase/work-notes/{filename}"
    if section_anchor:
        github_url += f"#{section_anchor}"

    # Build plain-text reference with full GitHub URL and section name
    # Use \url{} command for better line breaking of long URLs
    footnote_text = f"\\url{{{github_url}}} ({section_name} section)"
    if remaining:
        footnote_text += f": {remaining}"
    footnote_text += "."

    return footnote_text


def markdown_blockquote_to_latex(blockquote_md):
    """
    Convert markdown blockquote to LaTeX quote environment.
    Handles markdown formatting (bold, italic, links) inside blockquotes.

    Input: Lines starting with '> ' (markdown blockquote format)
    Output: LaTeX quote environment with converted formatting
    """
    if not blockquote_md:
        return None

    # Remove blockquote markers (> ) and join into single text
    lines = []
    for line in blockquote_md.split("\n"):
        if line.strip().startswith("> "):
            # Remove '> ' prefix and strip leading/trailing whitespace
            content = line.strip()[2:].strip()
            lines.append(content)
        elif line.strip():
            # Continuation line without marker
            lines.append(line.strip())

    if not lines:
        return None

    # Join all lines and convert markdown formatting
    text = " ".join(lines)
    latex_text = markdown_to_latex_text(text)

    # Wrap in LaTeX quote environment
    quote_latex = r"\begin{quote}" + "\n"
    quote_latex += latex_text + "\n"
    quote_latex += r"\end{quote}"

    return quote_latex


def format_details_blockquote(blockquote_md):
    """
    Format Details blockquotes cleanly: one full GitHub URL, then quoted section names.
    Also handles Implementation note blockquotes.

    Input: Blockquote markdown like:
    > **Details:** [link text 1](file.md#anchor1) ·
    > **Implementation note:** text about implementation

    Output: LaTeX formatted as:
    \\textbf{Details:} \\url{...} or \\textbf{Implementation note:} text
    """
    if not blockquote_md:
        return None

    # Check if this is a Details or Implementation note blockquote
    has_details = "Details:" in blockquote_md
    has_impl_note = (
        "Implementation note:" in blockquote_md or "note:" in blockquote_md.lower()
    )

    if not (has_details or has_impl_note):
        return None

    # Remove blockquote markers and flatten into single text
    lines = [
        l.strip()[2:].strip() if l.strip().startswith("> ") else l.strip()
        for l in blockquote_md.split("\n")
        if l.strip()
    ]
    text = " ".join(lines)

    # Handle Implementation note format (just plain text, no links)
    if has_impl_note and not has_details:
        # Remove **Implementation note:** prefix
        text = re.sub(r"\*\*Implementation note:\*\*\s*", "", text, flags=re.IGNORECASE)
        # Remove generic "**note:**" if present
        text = re.sub(r"\*\*[Nn]ote:\*\*\s*", "", text)
        return f"\\textbf{{Implementation note:}} {markdown_inline_formatting(text)}"

    # Handle Details format (with links)
    if has_details:
        # Remove **Details:** prefix to extract just the links
        text = re.sub(r"\*\*Details:\*\*\s*", "", text)

        # Find all markdown links: [text](url)
        link_pattern = r"\[([^\]]+)\]\(([^)]+)\)"
        matches = list(re.finditer(link_pattern, text))

        if not matches:
            # No links found, just return the text with Details label
            return f"\\textbf{{Details:}} {markdown_inline_formatting(text)}"

        # Extract filename from FIRST link
        first_link_url = matches[0].group(2)
        filename = (
            first_link_url.split("#")[0] if "#" in first_link_url else first_link_url
        )

        # Build full GitHub URL (one URL for all sections from same file)
        github_url = f"{REPO_BASE}/blob/{REPO_BRANCH}/uq-usecase/work-notes/{filename}"

        # Extract all link texts (section names as they appear in markdown)
        section_names = [latex_text_safe(match.group(1)) for match in matches]

        # Format: Details: [full URL]. Sections: "section1", "section2", "section3"
        # Section names are quoted so they stand out as section references
        sections_str = ", ".join(f'"{name}"' for name in section_names)

        details_latex = (
            f"\\textbf{{Details:}} \\url{{{github_url}}}. Sections: {sections_str}"
        )

        return details_latex

    return None


def generate_figure_label(section_num, caption_text=None):
    """
    Generate figure label matching the format used by extract_and_generate_mermaid.py.
    Examples:
        section=3, caption="One-time seed building (blue)..." → "fig_3_one_time_seed_building_blue_..."

    Returns: Label for looking up pre-generated PNG file
    """
    if caption_text:
        # Limit to first ~40 chars, same logic as mermaid extraction script
        short_caption = caption_text.strip()[:40].strip()

        # Convert to slug: remove special chars, spaces→underscores
        slug = re.sub(r"[^a-z0-9\s_-]", "", short_caption.lower())
        slug = re.sub(r"[\s-]+", "_", slug)  # Spaces/hyphens → underscores
        slug = re.sub(r"_+", "_", slug)  # Remove consecutive underscores
        slug = slug.strip("_")  # Remove leading/trailing underscores

        # Underscore, not colon: colons in filenames break on Overleaf/Windows.
        return f"fig_{section_num}_{slug}"
    else:
        return f"fig_{section_num}"


def find_figure_file(fig_label):
    """
    Find a pre-generated figure file by label.
    Looks in FIGURES_SOURCE_DIR for files matching fig_label.png.

    Returns: Path object if found, None otherwise.
    """
    fig_path = Path(FIGURES_SOURCE_DIR) / f"{fig_label}.png"
    if fig_path.exists():
        return fig_path
    return None


def extract_mermaid_from_markdown(text):
    """
    Extract mermaid diagram code block and caption from markdown.
    Note: mermaid_code is extracted for validation but figures are pre-generated.
    Returns (before_text, mermaid_code, caption, after_text) or (text, None, None, None) if no diagram.
    """
    # Match ```mermaid ... ``` block
    pattern = r"```mermaid\n(.*?)\n```"
    match = re.search(pattern, text, re.DOTALL)

    if not match:
        return text, None, None, None

    mermaid_code = match.group(1)
    before = text[: match.start()]
    after = text[match.end() :]

    # Extract caption (typically in italic after diagram: *Figure 3.1. ...*)
    caption_pattern = r"\*\*?Figure\s+[\d.]+\..*?\*\*?"
    caption_match = re.search(caption_pattern, after)
    caption = None
    if caption_match:
        caption = caption_match.group(0)
        # Remove caption from after_text to avoid duplication
        after = after[: caption_match.start()] + after[caption_match.end() :]

    return before.strip(), mermaid_code.strip(), caption, after.strip()


def mermaid_to_tikz(mermaid_code, diagram_label="diagram"):
    """
    Convert mermaid flowchart to TikZ code for inline LaTeX inclusion.
    Supports simple flowcharts with rectangular boxes, edges, and basic styling.

    This is a basic converter for simple flowcharts. Complex diagrams may need manual adjustment.
    """
    # Basic TikZ template for flowchart
    tikz_code = r"""\begin{tikzpicture}[
  node distance=2cm and 1cm,
  every node/.style={draw, rounded corners, minimum width=3cm, minimum height=0.8cm, align=center},
  box/.style={draw, rounded corners, fill=blue!20, minimum width=3cm, minimum height=0.8cm},
  once/.style={draw, rounded corners, fill=blue!10, minimum width=3cm, minimum height=0.8cm},
]
"""

    # Parse simple flowchart nodes and edges
    # Look for patterns like: NODEID["text"] and connections like NODEID --> NODEID

    # For now, return a placeholder with the mermaid code embedded as comment
    # User can then manually refine or use the PNG option
    tikz_code += f"  % Converted from Mermaid:\n"
    for line in mermaid_code.split("\n"):
        if line.strip():
            tikz_code += f"  % {line}\n"

    tikz_code += r"""\end{tikzpicture}"""

    return tikz_code


def mermaid_to_png(mermaid_code, output_path):
    """
    Convert mermaid diagram to PNG using mermaid-cli.
    Falls back to creating a placeholder PNG if mermaid-cli is not available.

    Requires: npm install -g @mermaid-js/mermaid-cli

    Args:
        mermaid_code: Mermaid diagram code (string)
        output_path: Path to save PNG file

    Returns: True if successful, False otherwise
    """
    import subprocess
    import tempfile

    try:
        # Write mermaid code to temp file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".mmd", delete=False) as f:
            f.write(mermaid_code)
            temp_file = f.name

        # Run mermaid-cli to generate PNG
        subprocess.run(
            ["mmdc", "-i", temp_file, "-o", output_path],
            check=True,
            capture_output=True,
        )

        # Clean up temp file
        Path(temp_file).unlink()
        return True
    except subprocess.CalledProcessError as e:
        print(f"Warning: mermaid-cli failed, creating placeholder PNG", file=sys.stderr)
        return create_placeholder_png(output_path)
    except FileNotFoundError:
        print(
            "Warning: mermaid-cli not found, creating placeholder PNG", file=sys.stderr
        )
        return create_placeholder_png(output_path)
    except Exception as e:
        print(
            f"Warning: Could not generate PNG, creating placeholder: {e}",
            file=sys.stderr,
        )
        return create_placeholder_png(output_path)


def create_placeholder_png(output_path):
    """
    Create a simple placeholder PNG for diagrams not yet rendered.
    Uses PIL if available, otherwise creates a minimal PNG.
    """
    try:
        from PIL import Image, ImageDraw

        # Create simple placeholder image
        width, height = 600, 400
        image = Image.new("RGB", (width, height), color="white")
        draw = ImageDraw.Draw(image)

        # Draw border
        draw.rectangle([10, 10, width - 10, height - 10], outline="gray", width=2)

        # Add placeholder text
        text = "[Diagram placeholder]\n\nFull diagram to be generated\nwith mermaid-cli or manual TikZ"
        draw.text((30, 180), text, fill="gray")

        # Save PNG
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        image.save(output_path, "PNG")
        return True
    except ImportError:
        # Fallback: create minimal PNG without PIL
        print(
            f"Note: PIL not available, using minimal PNG placeholder for {output_path}",
            file=sys.stderr,
        )
        # For now, just indicate failure - LaTeX will show file missing
        return False
    except Exception as e:
        print(f"Could not create placeholder PNG: {e}", file=sys.stderr)
        return False


# ============================================================================
# PHASE 3: Extract §3 (The system) and convert to LaTeX
# ============================================================================


def extract_section_3(md_content):
    """
    Extract section §3 "The system" from markdown.
    Assumes section 3 starts with "## 3. The system" and ends before "## 4. ".
    Returns tuple: (text_before_diagram, mermaid_code, caption, text_after_diagram, details_blockquotes)
    where details_blockquotes is a list of blockquote strings (may have 0, 1, or 2)
    """
    pattern_start = r"^## 3\.\s+The system\s*$"
    match_start = re.search(pattern_start, md_content, re.MULTILINE)
    if not match_start:
        raise ValueError("Could not find '## 3. The system' section")

    start_pos = match_start.end()

    # Find the end of §3 (before "## 4." or "---")
    pattern_end = r"^(?:---|## 4\.)"
    match_end = re.search(pattern_end, md_content[start_pos:], re.MULTILINE)

    if match_end:
        end_pos = start_pos + match_end.start()
    else:
        end_pos = len(md_content)

    section_3_text = md_content[start_pos:end_pos].strip()

    # Extract ALL Details blockquotes from entire section
    # Find all lines starting with > that contain "Details:"
    details_blockquotes = []
    lines = section_3_text.split("\n")
    i = 0

    while i < len(lines):
        if lines[i].strip().startswith(">") and "Details:" in lines[i]:
            # Found a Details blockquote - collect all consecutive > lines
            start_idx = i
            while i < len(lines) and lines[i].strip().startswith(">"):
                i += 1
            # Store the blockquote
            blockquote = "\n".join(lines[start_idx:i])
            details_blockquotes.append(blockquote)
        else:
            i += 1

    # Remove ALL blockquote lines from section and reconstruct
    lines_to_keep = []
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith(">") and "Details:" in lines[i]:
            # Skip this blockquote
            while i < len(lines) and lines[i].strip().startswith(">"):
                i += 1
        else:
            lines_to_keep.append(lines[i])
            i += 1

    section_3_text = "\n".join(lines_to_keep).strip()

    # Now extract mermaid diagram and caption from remaining text
    text_before, mermaid_code, caption, text_after = extract_mermaid_from_markdown(
        section_3_text
    )

    return text_before, mermaid_code, caption, text_after, details_blockquotes


def _png_size(path):
    """Return (width, height) in pixels from a PNG header, or None. Stdlib only."""
    import struct

    try:
        with open(path, "rb") as f:
            head = f.read(24)
        if len(head) < 24 or head[:8] != b"\x89PNG\r\n\x1a\n":
            return None
        return struct.unpack(">II", head[16:24])
    except Exception:
        return None


def _mermaid_width_for(png_source, fallback):
    """Pick an \\includegraphics width fraction from a mermaid PNG's aspect ratio.

    Mermaid flowcharts render at a fixed small font, so wide (horizontal) diagrams
    are illegible when shrunk. Give landscape diagrams the full text width; keep
    portrait/square ones at the caller's smaller default so they don't overflow.
    """
    size = _png_size(str(png_source))
    if not size or not size[1]:
        return fallback
    ratio = size[0] / size[1]
    if ratio >= 2.0:  # clearly horizontal → use all available width
        return 1.0
    if ratio >= 1.3:  # mildly landscape
        return max(fallback, 0.95)
    return fallback  # portrait/square → keep the narrower default


def generate_latex_figure_with_mermaid(
    mermaid_code, caption, fig_label="fig_system_overview", width=0.6
):
    """
    Generate LaTeX figure environment with pre-generated mermaid diagram.

    Looks for pre-generated PNG file in FIGURES_SOURCE_DIR.
    Copies it to FIGURES_OUTPUT_DIR and creates LaTeX includegraphics reference.

    Args:
        mermaid_code: (unused, for compatibility)
        caption: Figure caption text
        fig_label: Figure label for file lookup
        width: Width as fraction of textwidth (default 0.6, range 0.0-1.0)

    Note: Actual mermaid→PNG conversion is handled by extract_and_generate_mermaid.py
    """
    import shutil

    # Find pre-generated figure
    png_source = find_figure_file(fig_label)
    if not png_source:
        print(
            f"WARNING: Pre-generated figure not found: {fig_label}.png", file=sys.stderr
        )
        print(f"  → Run: python extract_and_generate_mermaid.py", file=sys.stderr)
        print(f"  → Then re-run this script", file=sys.stderr)
        return ""

    # Copy to both output directories
    png_filename = f"{fig_label}.png"

    # Size the figure from its aspect ratio: wide mermaid diagrams get full width.
    width = _mermaid_width_for(png_source, width)

    for figures_dir in [FIGURES_OUTPUT_DIR, REPO_FIGURES_OUTPUT_DIR]:
        Path(figures_dir).mkdir(parents=True, exist_ok=True)
        png_output_path = Path(figures_dir) / png_filename

        try:
            shutil.copy(str(png_source), str(png_output_path))
            print(f"  ✓ Copied figure: {png_output_path}")
        except Exception as e:
            print(
                f"WARNING: Could not copy figure to {figures_dir}: {e}", file=sys.stderr
            )
            return ""

    # Create LaTeX figure with a manual (non-floating-numbered) caption so the
    # figure numbers match the "Figure N.N" references in the prose.
    return _emit_figure_latex(png_filename, caption, width)


def _emit_figure_latex(png_filename, caption_md, width=0.75):
    """Emit a figure environment with a manual caption line (keeps 'Figure N.N').

    Using plain text rather than \\caption avoids LaTeX auto-numbering, which does
    not match the manual 'Figure 7.2' style used throughout the report.
    """
    latex = r"\begin{figure}[htb]" + "\n"
    latex += r"  \centering" + "\n"
    latex += (
        r"  \includegraphics[width="
        + f"{width:.2f}"
        + r"\textwidth]{"
        + FIG_TEX_DIR
        + "/"
        + png_filename
        + "}\n"
    )
    if caption_md:
        cap = caption_md.replace("*", "").strip()
        cap = markdown_inline_formatting(cap)
        latex += r"  \par\vspace{0.4em}" + "\n"
        latex += r"  {\small " + cap + "}\n"
    latex += r"\end{figure}" + "\n"
    return latex


def generate_latex_figure_image(image_rel_path, caption_md, width=0.85):
    """Convert a Markdown image (![alt](figures/foo.png)) to a LaTeX figure.

    Copies the referenced PNG from FIGURES_SOURCE_DIR into the output figs/
    directories and emits a figure with a manual caption. Returns "" if the
    source image cannot be found.
    """
    import shutil

    basename = Path(image_rel_path).name
    png_source = Path(FIGURES_SOURCE_DIR) / basename
    if not png_source.exists():
        print(f"WARNING: Image not found: {png_source}", file=sys.stderr)
        return ""

    for figures_dir in [FIGURES_OUTPUT_DIR, REPO_FIGURES_OUTPUT_DIR]:
        Path(figures_dir).mkdir(parents=True, exist_ok=True)
        try:
            shutil.copy(str(png_source), str(Path(figures_dir) / basename))
            print(f"  ✓ Copied figure: {basename}")
        except Exception as e:
            print(
                f"WARNING: Could not copy image to {figures_dir}: {e}", file=sys.stderr
            )
            return ""

    return _emit_figure_latex(basename, caption_md, width)


def generate_latex_section_3(
    text_before, mermaid_code, caption, text_after, details_md_list, section_num=3
):
    """
    Convert extracted §3 markdown to LaTeX.
    Includes: intro text, mermaid figure with caption, followup text, details boxes.

    Args:
        section_num: Section number (for figure label generation)
        details_md_list: List of blockquote strings (may be empty)
    """
    latex_output = r"""\subsection*{3. The system}
\addcontentsline{toc}{subsection}{3. The system}
\label{sec:the_system}

"""

    # Add text before diagram
    if text_before:
        latex_text = markdown_to_latex_text(text_before)
        latex_output += latex_text + "\n\n"

    # Add mermaid figure
    if mermaid_code and caption:
        # Extract caption text and remove asterisks and "Figure N.N. " prefix
        caption_text = caption.replace("*", "").strip()
        # Remove "Figure N.N. " prefix (e.g., "Figure 3.1. " → "")
        caption_text = re.sub(r"^Figure\s+[\d.]+\.\s*", "", caption_text)

        fig_label = generate_figure_label(section_num, caption_text)

        fig_latex = generate_latex_figure_with_mermaid(mermaid_code, caption, fig_label)
        if fig_latex:  # Only add if figure was successfully found/copied
            latex_output += fig_latex + "\n"

    # Add text after diagram
    if text_after:
        latex_text = markdown_to_latex_text(text_after)
        latex_output += latex_text + "\n"

    # Add all details boxes (if present)
    if details_md_list and INCLUDE_FOOTNOTES:
        for details_md in details_md_list:
            details_latex = format_details_blockquote(details_md)
            if details_latex:
                latex_output += "\n"
                latex_output += r"\par\vspace{0.8em}" + "\n"
                latex_output += r"\noindent\fbox{\parbox[t]{0.95\textwidth}{%" + "\n"
                latex_output += r"\setlength{\parskip}{0.3em}%" + "\n"
                latex_output += r"\setlength{\baselineskip}{1.15\baselineskip}%" + "\n"
                latex_output += details_latex + "\n"
                latex_output += r"}}" + "\n"
                latex_output += r"\vspace{0.3em}" + "\n"

    return latex_output


def generate_latex_section_2(table_md, footnote_md):
    """
    Convert extracted §2 markdown to LaTeX.
    Uses definition-list format (flowing text with bold terms) instead of table.
    Much better for long glossaries.
    """
    latex_output = r"""\subsection*{2. Glossary}
\addcontentsline{toc}{subsection}{2. Glossary}
\label{sec:glossary}

Main terms used in our pilot study that follows.

"""

    # Add glossary as definition list (flowing text format, not table)
    latex_output += markdown_glossary_to_latex_deflist(table_md)
    latex_output += "\n"

    # Add note box after glossary (if requested)
    if INCLUDE_FOOTNOTES:
        details_latex = format_details_blockquote(footnote_md)
        if details_latex:
            # Create a compact styled details box with minimal spacing
            latex_output += r"\vspace{-0.3em}" + "\n"
            latex_output += r"\noindent\fbox{\parbox[t]{0.98\textwidth}{%" + "\n"
            latex_output += r"\setlength{\parskip}{0pt}%" + "\n"
            latex_output += details_latex + "\n"
            latex_output += r"}}" + "\n"
            latex_output += r"\vspace{-0.2em}" + "\n"

    return latex_output


# PHASE 4: Extract §4.1 (Building one operating condition - Commitment subsection)


def extract_section_4_1(md_content):
    """
    Extract subsection §4.1 "Commitment: a feasible band, not an economic optimum" from markdown.
    Assumes 4.1 starts with "### 4.1 " and ends before "### 4.2 " or end of section 4.
    Returns tuple: (before_table_text, table_md, after_table_text, details_blockquotes)
    This preserves the original order: text before table, table, text after table, then blockquotes.
    """
    pattern_start = r"^### 4\.1\s+Commitment.*$"
    match_start = re.search(pattern_start, md_content, re.MULTILINE)
    if not match_start:
        raise ValueError("Could not find '### 4.1 Commitment' section")

    # Skip past the heading line to the next line
    start_pos = match_start.end()
    # Find the next newline to skip the rest of the heading line
    next_newline = md_content.find("\n", start_pos)
    if next_newline != -1:
        start_pos = next_newline + 1

    # Find the end of §4.1 (before "### 4.2" or "## 5.")
    pattern_end = r"^(?:### 4\.2|## 5\.)"
    match_end = re.search(pattern_end, md_content[start_pos:], re.MULTILINE)

    if match_end:
        end_pos = start_pos + match_end.start()
    else:
        end_pos = len(md_content)

    section_4_1_text = md_content[start_pos:end_pos].strip()

    # Extract table and blockquotes, preserving order
    table_md = ""
    details_blockquotes = []
    before_table_lines = []
    after_table_lines = []
    found_table = False
    lines = section_4_1_text.split("\n")
    i = 0

    while i < len(lines):
        line = lines[i]

        # Check for blockquote (Details:)
        if line.strip().startswith(">") and "Details:" in line:
            # Collect all consecutive > lines
            blockquote_lines = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                blockquote_lines.append(lines[i])
                i += 1
            details_blockquotes.append("\n".join(blockquote_lines))
        # Check for table (lines with pipes)
        elif "|" in line and i + 1 < len(lines) and "---" in lines[i + 1]:
            # Found a markdown table - collect header, separator, and rows
            table_lines = [line]
            i += 1
            table_lines.append(lines[i])  # separator line
            i += 1
            while i < len(lines) and "|" in lines[i]:
                table_lines.append(lines[i])
                i += 1
            table_md = "\n".join(table_lines)
            found_table = True
        else:
            # Add to before_table or after_table based on whether we've seen the table yet
            if not found_table:
                before_table_lines.append(line)
            else:
                after_table_lines.append(line)
            i += 1

    before_table = "\n".join(before_table_lines).strip() if before_table_lines else ""
    after_table = "\n".join(after_table_lines).strip() if after_table_lines else ""

    return before_table, table_md, after_table, details_blockquotes


def generate_latex_section_4_1(
    before_table, table_md, after_table, details_blockquotes
):
    """
    Convert extracted §4.1 markdown to LaTeX.
    Handles text before table, inline table, text after table, and details blockquotes.
    Preserves the original markdown order.
    """
    latex_output = r"""\subsubsection*{4.1 Commitment: a feasible band, not an economic optimum}
\addcontentsline{toc}{subsubsection}{4.1 Commitment: a feasible band, not an economic optimum}
\label{sec:commitment}

"""

    # Convert text BEFORE table (handles bullet lists, blockquotes, bold, italic, etc.)
    if before_table.strip():
        latex_text = markdown_to_latex_text(before_table)
        latex_output += latex_text + "\n"

    # Convert table to LaTeX tabular if present (use float=False to keep inline, not floating)
    # For small 2-row tables, use compact format
    if table_md.strip():
        latex_output += "\n"
        latex_table = markdown_table_to_latex_tabular(
            table_md, float=False, compact=True
        )
        latex_output += latex_table
        latex_output += "\n"

    # Convert text AFTER table (handles bullet lists, blockquotes, bold, italic, etc.)
    if after_table.strip():
        latex_text = markdown_to_latex_text(after_table)
        latex_output += latex_text + "\n"

    # Add details blockquote if present
    if details_blockquotes and INCLUDE_FOOTNOTES:
        latex_output += r"\par\vspace{0.8em}" + "\n"
        latex_output += r"\noindent\fbox{\parbox[t]{0.98\textwidth}{%" + "\n"
        latex_output += r"\setlength{\parskip}{0pt}%" + "\n"

        # Process all blockquotes - use specialized Details handler for cleaner formatting
        for blockquote in details_blockquotes:
            # Try to format as Details blockquote first
            details_latex = format_details_blockquote(blockquote)
            if details_latex:
                latex_output += details_latex + "\n"
            else:
                # Fallback to generic blockquote formatting
                blockquote_latex = markdown_blockquote_to_latex(blockquote)
                if blockquote_latex:
                    latex_output += blockquote_latex + "\n"

        latex_output += r"}}" + "\n"
        latex_output += r"\vspace{0.2em}" + "\n"

    return latex_output


# PHASE 4: Extract §4 (Building one operating condition) with mermaid diagram
# and §4.1 (Commitment subsection)


def extract_section_4(md_content):
    """
    Extract section §4 "Building one operating condition" from markdown.
    Includes intro text (if any) and mermaid diagram with caption.
    Extracts from "## 4. Building one operating condition" to "### 4.1 ".
    Returns tuple: (text_before_diagram, mermaid_code, caption, text_after_diagram)
    """
    pattern_start = r"^## 4\.\s+Building one operating condition\s*$"
    match_start = re.search(pattern_start, md_content, re.MULTILINE)
    if not match_start:
        raise ValueError(
            "Could not find '## 4. Building one operating condition' section"
        )

    start_pos = match_start.end()

    # Find the end of §4 (before "### 4.1")
    pattern_end = r"^### 4\.1\s+"
    match_end = re.search(pattern_end, md_content[start_pos:], re.MULTILINE)

    if match_end:
        end_pos = start_pos + match_end.start()
    else:
        end_pos = len(md_content)

    section_4_text = md_content[start_pos:end_pos].strip()

    # Use the existing mermaid extraction function
    text_before, mermaid_code, caption, text_after = extract_mermaid_from_markdown(
        section_4_text
    )

    return text_before, mermaid_code, caption, text_after


def generate_latex_section_4(text_before, mermaid_code, caption, text_after):
    """
    Convert extracted §4 markdown to LaTeX.
    Includes: intro text, mermaid figure with caption, followup text.
    """
    latex_output = r"""\subsection*{4. Building one operating condition}
\addcontentsline{toc}{subsection}{4. Building one operating condition}
\label{sec:building_condition}

"""

    # Add text before diagram (usually empty for section 4)
    if text_before and text_before.strip():
        latex_text = markdown_to_latex_text(text_before)
        latex_output += latex_text + "\n\n"

    # Add mermaid figure (larger size for section 4)
    if mermaid_code and caption:
        # Extract caption text and remove asterisks and "Figure N.N. " prefix
        caption_text = caption.replace("*", "").strip()
        # Remove "Figure N.N. " prefix (e.g., "Figure 4.1. " → "")
        caption_text = re.sub(r"^Figure\s+[\d.]+\.\s*", "", caption_text)

        fig_label = generate_figure_label(4, caption_text)

        fig_latex = generate_latex_figure_with_mermaid(
            mermaid_code, caption, fig_label, width=0.95
        )
        if fig_latex:  # Only add if figure was successfully found/copied
            latex_output += fig_latex + "\n"

    # Add text after diagram (usually empty for section 4)
    if text_after and text_after.strip():
        latex_text = markdown_to_latex_text(text_after)
        latex_output += latex_text + "\n"

    return latex_output


# PHASE 4.2: Extract §4.2, §4.2a, §4.2b (Dispatch and reactive limits)


def extract_section_4_2(md_content):
    """
    Extract section §4.2 "Dispatch: DC OPF then AC PF" from markdown.
    Returns: (text_before, mermaid_code, caption, text_after, details_blockquotes)
    Note: details_blockquotes includes both Details and Implementation note blockquotes
    """
    pattern_start = r"^### 4\.2\s+Dispatch"
    match_start = re.search(pattern_start, md_content, re.MULTILINE)
    if not match_start:
        return ("", "", "", "", [])

    start_pos = match_start.end()
    next_newline = md_content.find("\n", start_pos)
    if next_newline != -1:
        start_pos = next_newline + 1

    # Section ends at 4.2a, 4.3, or ## 5
    pattern_end = r"^(?:### 4\.2a|### 4\.2b|### 4\.3|## 5\.)"
    match_end = re.search(pattern_end, md_content[start_pos:], re.MULTILINE)
    end_pos = start_pos + match_end.start() if match_end else len(md_content)

    section_text = md_content[start_pos:end_pos].strip()

    # Extract text before diagram, diagram, caption, text after, and blockquotes (Details or Implementation note)
    text_before_lines = []
    mermaid_code = ""
    caption = ""
    text_after_lines = []
    details_blockquotes = (
        []
    )  # Collects both "Details:" and "Implementation note:" blockquotes
    found_mermaid = False
    found_caption = False

    lines = section_text.split("\n")
    i = 0

    while i < len(lines):
        line = lines[i]

        if line.strip().startswith("```mermaid"):
            found_mermaid = True
            mermaid_lines = []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("```"):
                mermaid_lines.append(lines[i])
                i += 1
            mermaid_code = "\n".join(mermaid_lines).strip()
            i += 1
        elif found_mermaid and line.strip().startswith("*Figure"):
            # Caption line
            caption = line.strip().strip("*").strip()
            found_caption = True
            i += 1
        elif line.strip().startswith(">") and (
            "Details:" in line or "note:" in line.lower()
        ):
            # Blockquote with Details or Implementation note
            blockquote_lines = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                blockquote_lines.append(lines[i])
                i += 1
            details_blockquotes.append("\n".join(blockquote_lines))
        else:
            if not found_mermaid:
                text_before_lines.append(line)
            else:
                text_after_lines.append(line)
            i += 1

    text_before = "\n".join(text_before_lines).strip()
    text_after = "\n".join(text_after_lines).strip()

    return text_before, mermaid_code, caption, text_after, details_blockquotes


def extract_section_4_2a(md_content):
    """
    Extract section §4.2a "Reactive limits: PV to PQ switching" from markdown.
    Returns: (text, details_blockquotes)
    """
    pattern_start = r"^### 4\.2a\s+Reactive limits"
    match_start = re.search(pattern_start, md_content, re.MULTILINE)
    if not match_start:
        return ("", [])

    start_pos = match_start.end()
    next_newline = md_content.find("\n", start_pos)
    if next_newline != -1:
        start_pos = next_newline + 1

    # Section ends at 4.2b, 4.3, or ## 5
    pattern_end = r"^(?:### 4\.2b|### 4\.3|## 5\.)"
    match_end = re.search(pattern_end, md_content[start_pos:], re.MULTILINE)
    end_pos = start_pos + match_end.start() if match_end else len(md_content)

    section_text = md_content[start_pos:end_pos].strip()

    # Extract main text and details blockquotes
    text_lines = []
    details_blockquotes = []

    lines = section_text.split("\n")
    i = 0

    while i < len(lines):
        line = lines[i]

        if line.strip().startswith(">") and "Details:" in line:
            blockquote_lines = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                blockquote_lines.append(lines[i])
                i += 1
            details_blockquotes.append("\n".join(blockquote_lines))
        else:
            text_lines.append(line)
            i += 1

    text = "\n".join(text_lines).strip()
    return text, details_blockquotes


def extract_section_4_2b(md_content):
    """
    Extract section §4.2b "AC OPF as a second Q-limit route" from markdown.
    Returns: (text, details_blockquotes)
    """
    pattern_start = r"^### 4\.2b\s+AC OPF"
    match_start = re.search(pattern_start, md_content, re.MULTILINE)
    if not match_start:
        return ("", [])

    start_pos = match_start.end()
    next_newline = md_content.find("\n", start_pos)
    if next_newline != -1:
        start_pos = next_newline + 1

    # Section ends at 4.3 or ## 5
    pattern_end = r"^(?:### 4\.3|## 5\.)"
    match_end = re.search(pattern_end, md_content[start_pos:], re.MULTILINE)
    end_pos = start_pos + match_end.start() if match_end else len(md_content)

    section_text = md_content[start_pos:end_pos].strip()

    # Extract main text and details blockquotes
    text_lines = []
    details_blockquotes = []

    lines = section_text.split("\n")
    i = 0

    while i < len(lines):
        line = lines[i]

        if line.strip().startswith(">") and "Details:" in line:
            blockquote_lines = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                blockquote_lines.append(lines[i])
                i += 1
            details_blockquotes.append("\n".join(blockquote_lines))
        else:
            text_lines.append(line)
            i += 1

    text = "\n".join(text_lines).strip()
    return text, details_blockquotes


def generate_latex_section_4_2(
    text_before, mermaid_code, caption, text_after, details_blockquotes
):
    """Convert extracted §4.2 markdown to LaTeX."""
    latex_output = r"\subsubsection*{4.2 Dispatch: DC OPF then AC PF}" + "\n"
    latex_output += (
        r"\addcontentsline{toc}{subsubsection}{4.2 Dispatch: DC OPF then AC PF}" + "\n"
    )
    latex_output += r"\label{sec:dispatch_dc_opf_ac_pf}" + "\n\n"

    # Separate blockquotes: Implementation note goes at beginning, Details goes at end
    impl_notes = []
    details_notes = []
    for blockquote in details_blockquotes:
        if "note:" in blockquote.lower():
            impl_notes.append(blockquote)
        elif "Details:" in blockquote:
            details_notes.append(blockquote)

    # Add Implementation note blockquotes at the beginning (as boxed content)
    if impl_notes and INCLUDE_FOOTNOTES:
        latex_output += r"\noindent\fbox{\parbox[t]{0.98\textwidth}{%" + "\n"
        latex_output += r"\setlength{\parskip}{0pt}%" + "\n"
        for blockquote in impl_notes:
            note_latex = format_details_blockquote(blockquote)
            if note_latex:
                latex_output += note_latex + "\n"
        latex_output += r"}}" + "\n"
        latex_output += r"\vspace{0.4em}" + "\n\n"

    if text_before:
        latex_text = markdown_to_latex_text(text_before)
        latex_output += latex_text + "\n\n"

    if mermaid_code and caption:
        caption_text = caption.replace("*", "").strip()
        caption_text = re.sub(r"^Figure\s+[\d.]+\.\s*", "", caption_text)
        fig_label = generate_figure_label(4, caption_text)
        fig_latex = generate_latex_figure_with_mermaid(
            mermaid_code, caption, fig_label, width=0.95
        )
        if fig_latex:
            latex_output += fig_latex + "\n"

    if text_after:
        latex_text = markdown_to_latex_text(text_after)
        latex_output += latex_text + "\n"

    # Add Details blockquotes at the end (as boxed content)
    if details_notes and INCLUDE_FOOTNOTES:
        latex_output += r"\par\vspace{0.8em}" + "\n"
        latex_output += r"\noindent\fbox{\parbox[t]{0.98\textwidth}{%" + "\n"
        latex_output += r"\setlength{\parskip}{0pt}%" + "\n"
        for blockquote in details_notes:
            details_latex = format_details_blockquote(blockquote)
            if details_latex:
                latex_output += details_latex + "\n"
        latex_output += r"}}" + "\n"
        latex_output += r"\vspace{0.2em}" + "\n"

    return latex_output


def generate_latex_section_4_2a(text, details_blockquotes):
    """Convert extracted §4.2a markdown to LaTeX."""
    # Use unnumbered subsubsection with custom numbering and TOC entry
    latex_output = r"\subsubsection*{4.2a Reactive limits: PV to PQ switching}" + "\n"
    latex_output += (
        r"\addcontentsline{toc}{subsubsection}{4.2a Reactive limits: PV to PQ switching}"
        + "\n"
    )
    latex_output += r"\label{sec:reactive_limits_pv_pq}" + "\n\n"

    if text:
        latex_text = markdown_to_latex_text(text)
        latex_output += latex_text + "\n"

    if details_blockquotes and INCLUDE_FOOTNOTES:
        latex_output += r"\par\vspace{0.8em}" + "\n"
        latex_output += r"\noindent\fbox{\parbox[t]{0.98\textwidth}{%" + "\n"
        latex_output += r"\setlength{\parskip}{0pt}%" + "\n"
        for blockquote in details_blockquotes:
            details_latex = format_details_blockquote(blockquote)
            if details_latex:
                latex_output += details_latex + "\n"
        latex_output += r"}}" + "\n"
        latex_output += r"\vspace{0.2em}" + "\n"

    return latex_output


def generate_latex_section_4_2b(text, details_blockquotes):
    """Convert extracted §4.2b markdown to LaTeX."""
    # Use unnumbered subsubsection with custom numbering and TOC entry
    latex_output = r"\subsubsection*{4.2b AC OPF as a second Q-limit route}" + "\n"
    latex_output += (
        r"\addcontentsline{toc}{subsubsection}{4.2b AC OPF as a second Q-limit route}"
        + "\n"
    )
    latex_output += r"\label{sec:ac_opf_q_limits}" + "\n\n"

    if text:
        latex_text = markdown_to_latex_text(text)
        latex_output += latex_text + "\n"

    if details_blockquotes and INCLUDE_FOOTNOTES:
        latex_output += r"\par\vspace{0.8em}" + "\n"
        latex_output += r"\noindent\fbox{\parbox[t]{0.98\textwidth}{%" + "\n"
        latex_output += r"\setlength{\parskip}{0pt}%" + "\n"
        for blockquote in details_blockquotes:
            details_latex = format_details_blockquote(blockquote)
            if details_latex:
                latex_output += details_latex + "\n"
        latex_output += r"}}" + "\n"
        latex_output += r"\vspace{0.2em}" + "\n"

    return latex_output


# PHASE 5: Extract §4.3 (Wind becomes negative load)


def extract_section_4_3(md_content):
    """
    Extract section §4.3 "Wind becomes negative load" from markdown.
    Returns: (impl_note, text, details_blockquotes)
    """
    pattern_start = r"^### 4\.3\s+Wind"
    match_start = re.search(pattern_start, md_content, re.MULTILINE)
    if not match_start:
        return ("", "", [])

    start_pos = match_start.end()
    next_newline = md_content.find("\n", start_pos)
    if next_newline != -1:
        start_pos = next_newline + 1

    # Section ends at 4.4 or later
    pattern_end = r"^### 4\.4"
    match_end = re.search(pattern_end, md_content[start_pos:], re.MULTILINE)
    end_pos = start_pos + match_end.start() if match_end else len(md_content)

    section_text = md_content[start_pos:end_pos].strip()

    # Extract: implementation note, main text, and details blockquotes
    impl_note = ""
    text_lines = []
    details_blockquotes = []

    lines = section_text.split("\n")
    i = 0

    # Check if first blockquote is an implementation note
    if i < len(lines) and lines[i].strip().startswith(">"):
        if "note:" in lines[i].lower():
            blockquote_lines = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                blockquote_lines.append(lines[i])
                i += 1
            impl_note = "\n".join(blockquote_lines)

    # Collect main text and Details blockquotes
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if stripped.startswith(">") and (
            "Details:" in stripped or "note:" in stripped.lower()
        ):
            # Details blockquote
            blockquote_lines = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                blockquote_lines.append(lines[i])
                i += 1
            details_blockquotes.append("\n".join(blockquote_lines))
        else:
            text_lines.append(
                line
            )  # Keep all lines including blank ones for paragraph breaks
            i += 1

    text = "\n".join(text_lines).strip()
    return impl_note, text, details_blockquotes


def generate_latex_section_4_3(impl_note, text, details_blockquotes):
    """Convert extracted §4.3 markdown to LaTeX."""
    latex_output = r"\subsubsection*{4.3 Wind becomes negative load}" + "\n"
    latex_output += (
        r"\addcontentsline{toc}{subsubsection}{4.3 Wind becomes negative load}" + "\n"
    )
    latex_output += r"\label{sec:wind_negative_load}" + "\n\n"

    # Add implementation note at beginning if present
    if impl_note and INCLUDE_FOOTNOTES:
        latex_output += r"\noindent\fbox{\parbox[t]{0.98\textwidth}{%" + "\n"
        latex_output += r"\setlength{\parskip}{0pt}%" + "\n"
        note_latex = format_details_blockquote(impl_note)
        if note_latex:
            latex_output += note_latex + "\n"
        latex_output += r"}}" + "\n"
        latex_output += r"\vspace{0.4em}" + "\n\n"

    # Add main text
    if text:
        latex_text = markdown_to_latex_text(text)
        latex_output += latex_text + "\n"

    # Add Details blockquotes at the end
    if details_blockquotes and INCLUDE_FOOTNOTES:
        latex_output += r"\par\vspace{0.8em}" + "\n"
        latex_output += r"\noindent\fbox{\parbox[t]{0.98\textwidth}{%" + "\n"
        latex_output += r"\setlength{\parskip}{0pt}%" + "\n"
        for blockquote in details_blockquotes:
            details_latex = format_details_blockquote(blockquote)
            if details_latex:
                latex_output += details_latex + "\n"
        latex_output += r"}}" + "\n"
        latex_output += r"\vspace{0.2em}" + "\n"

    return latex_output


def extract_section_4_4(md_content):
    """
    Extract section §4.4 "The datacenter" from markdown.
    Returns: (impl_note, text) where text includes Details blockquotes inline.
    """
    pattern_start = r"^### 4\.4\s+The datacenter"
    match_start = re.search(pattern_start, md_content, re.MULTILINE)
    if not match_start:
        return ("", "")

    start_pos = match_start.end()
    next_newline = md_content.find("\n", start_pos)
    if next_newline != -1:
        start_pos = next_newline + 1

    # Section ends at §5 or later (next ## heading that's not ### subsection)
    pattern_end = r"^## 5\."
    match_end = re.search(pattern_end, md_content[start_pos:], re.MULTILINE)
    end_pos = start_pos + match_end.start() if match_end else len(md_content)

    section_text = md_content[start_pos:end_pos].strip()

    # Extract: implementation note at beginning, then all remaining text inline
    impl_note = ""
    text_lines = []

    lines = section_text.split("\n")
    i = 0

    # Check if first blockquote is an implementation note
    if i < len(lines) and lines[i].strip().startswith(">"):
        if "note:" in lines[i].lower():
            blockquote_lines = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                blockquote_lines.append(lines[i])
                i += 1
            impl_note = "\n".join(blockquote_lines)

    # Collect all remaining text (including Details blockquotes inline)
    while i < len(lines):
        text_lines.append(
            lines[i]
        )  # Keep all lines including blank ones and blockquotes
        i += 1

    text = "\n".join(text_lines).strip()
    return impl_note, text


def generate_latex_section_4_4(impl_note, text):
    """Convert extracted §4.4 markdown to LaTeX, with Details blockquotes inline."""
    latex_output = r"\subsubsection*{4.4 The datacenter}" + "\n"
    latex_output += r"\addcontentsline{toc}{subsubsection}{4.4 The datacenter}" + "\n"
    latex_output += r"\label{sec:datacenter}" + "\n\n"

    # Add implementation note at beginning if present
    if impl_note and INCLUDE_FOOTNOTES:
        latex_output += r"\noindent\fbox{\parbox[t]{0.98\textwidth}{%" + "\n"
        latex_output += r"\setlength{\parskip}{0pt}%" + "\n"
        note_latex = format_details_blockquote(impl_note)
        if note_latex:
            latex_output += note_latex + "\n"
        latex_output += r"}}" + "\n"
        latex_output += r"\vspace{0.4em}" + "\n\n"

    # Add main text with Details blockquotes inline (converted by markdown_to_latex_text)
    if text:
        latex_text = markdown_to_latex_text(text)
        latex_output += latex_text + "\n"

    return latex_output


# PHASE 7+: Extract §5-12 (Remaining sections)


def extract_section_n(md_content, section_num):
    """Generic extractor for sections 5-12."""
    pattern_start = rf"^## {section_num}\."
    match_start = re.search(pattern_start, md_content, re.MULTILINE)
    if not match_start:
        return ""

    start_pos = match_start.end()
    next_newline = md_content.find("\n", start_pos)
    if next_newline != -1:
        start_pos = next_newline + 1

    # Section ends at the next top-level "## " heading (next section, "## 8a",
    # or the standalone "## References" list, which the bibliography handles).
    pattern_end = r"^## "
    match_end = re.search(pattern_end, md_content[start_pos:], re.MULTILINE)
    end_pos = start_pos + match_end.start() if match_end else len(md_content)

    section_text = md_content[start_pos:end_pos].strip()
    return section_text


def _emit_mermaid_from_caption(caption, section_num):
    """Emit a mermaid figure (via its pre-generated PNG) from its caption line."""
    if not caption:
        return ""
    caption_text = caption.replace("*", "").strip()
    caption_text = re.sub(r"^Figure\s+[\d.]+\.\s*", "", caption_text)
    fig_label = generate_figure_label(section_num, caption_text)
    return generate_latex_figure_with_mermaid(None, caption, fig_label, width=0.72)


def markdown_section_body_to_latex(text, section_num):
    """Robust block-level converter for section bodies (§5-§12).

    Handles markdown tables, $$ display math, images, fenced code and mermaid
    blocks explicitly, and delegates prose, lists, blockquotes (incl. Details
    boxes) and headings to markdown_to_latex_text.
    """
    lines = text.split("\n")
    n = len(lines)
    out = []
    buf = []

    def flush():
        if buf:
            chunk = "\n".join(buf)
            if chunk.strip():
                out.append(markdown_to_latex_text(chunk))
            buf.clear()

    i = 0
    while i < n:
        line = lines[i]
        stripped = line.strip()

        # Fenced code / mermaid block
        if stripped.startswith("```"):
            flush()
            lang = stripped[3:].strip().lower()
            body = []
            j = i + 1
            while j < n and not lines[j].strip().startswith("```"):
                body.append(lines[j])
                j += 1
            if lang == "mermaid":
                k = j + 1
                while k < n and not lines[k].strip():
                    k += 1
                caption = None
                next_i = j + 1
                if k < n and lines[k].strip().startswith("*") and "Figure" in lines[k]:
                    caption = lines[k].strip()
                    next_i = k + 1
                fig = _emit_mermaid_from_caption(caption, section_num)
                if fig:
                    out.append(fig)
                i = next_i
            else:
                code = "\n".join(body)
                out.append("\\begin{verbatim}\n" + code + "\n\\end{verbatim}")
                i = j + 1
            continue

        # Display math $$...$$
        if stripped.startswith("$$"):
            flush()
            if stripped.endswith("$$") and len(stripped) > 3:
                math = stripped[2:-2]
                i += 1
            else:
                body = [stripped[2:]]
                j = i + 1
                while j < n and not lines[j].strip().endswith("$$"):
                    body.append(lines[j])
                    j += 1
                if j < n:
                    body.append(lines[j].strip()[:-2])
                math = "\n".join(body)
                i = j + 1
            out.append("\\[\n" + latexify_math(math.strip()) + "\n\\]")
            continue

        # Image ![alt](path), optionally followed by a *Figure ...* caption
        img = re.match(r"^\s*!\[([^\]]*)\]\(([^)]+)\)\s*$", line)
        if img:
            flush()
            path = img.group(2)
            k = i + 1
            while k < n and not lines[k].strip():
                k += 1
            caption = None
            next_i = i + 1
            if k < n and lines[k].strip().startswith("*") and "Figure" in lines[k]:
                caption = lines[k].strip()
                next_i = k + 1
            fig = generate_latex_figure_image(path, caption)
            if fig:
                out.append(fig)
            i = next_i
            continue

        # Markdown table (header row followed by |---| separator row)
        if (
            stripped.startswith("|")
            and i + 1 < n
            and re.match(r"^\s*\|[\s:|-]+\|?\s*$", lines[i + 1])
        ):
            flush()
            tbl = [line]
            j = i + 1
            while j < n and lines[j].strip().startswith("|"):
                tbl.append(lines[j])
                j += 1
            out.append(markdown_table_to_latex_tabular("\n".join(tbl), float=False))
            i = j
            continue

        buf.append(line)
        i += 1

    flush()
    return "\n\n".join(b for b in out if b and b.strip())


def generate_latex_section_n(section_num, section_title, text, label_name):
    """Generic generator for sections 5-12 using the robust block converter."""
    if not text:
        return ""

    heading = f"{section_num}. {section_title}"
    latex_output = rf"\subsection*{{{heading}}}" + "\n"
    latex_output += rf"\addcontentsline{{toc}}{{subsection}}{{{heading}}}" + "\n"
    latex_output += rf"\label{{sec:{label_name}}}" + "\n\n"
    latex_output += markdown_section_body_to_latex(text, section_num) + "\n\n"

    return latex_output


# ============================================================================
# POST-PROCESSOR: numbered headings + cross-reference resolution
# ============================================================================


def apply_numbered_headings_and_refs(latex):
    r"""Convert starred headings to numbered ones and wire up cross-references.

    The report loads titlesec, which (a) does not number \subsection*/\subsubsection*
    and (b) typesets their titles stretched to full justified width. The generators
    emit starred headings carrying a *manual* number in the title (e.g.
    "\subsection*{2. Glossary}"), which therefore lose the report's 9.1.x prefix and
    render stretched. This pass rewrites them to numbered commands so LaTeX supplies
    the correct 9.1.x / 9.1.x.y numbers automatically:

        \subsection*{2. Glossary}              -> \subsection{Glossary}
        \addcontentsline{toc}{subsection}{...}    \label{uqp:2}

    The manual number becomes a \label{uqp:<number>}, and every prose cross-reference
    ("{\textsection}7.6", i.e. the converted form of "§7.6") is rewritten to
    "{\textsection}\ref{uqp:7.6}" so it always prints the heading's real number.
    """
    heading_re = re.compile(r"^\\(subsection|subsubsection|paragraph)\*\{(.*)\}\s*$")
    # Leading manual number: "2", "4.2", "4.2a", "7.0" (optional trailing dot).
    num_re = re.compile(r"^(\d+(?:\.\d+)*[a-z]?)\.?\s+(.*)$")

    lines = latex.split("\n")
    out = []
    i = 0
    while i < len(lines):
        m = heading_re.match(lines[i])
        if not m:
            out.append(lines[i])
            i += 1
            continue

        cmd, title = m.group(1), m.group(2)
        # Numbered headings add themselves to the ToC, so drop a trailing
        # \addcontentsline line the generators emitted for the starred form.
        drop_acl = i + 1 < len(lines) and lines[i + 1].lstrip().startswith(
            r"\addcontentsline{toc}"
        )

        nm = num_re.match(title)
        if nm:
            number, clean_title = nm.group(1), nm.group(2)
            out.append(rf"\{cmd}{{{clean_title}}}")
            out.append(rf"\label{{uqp:{number}}}")
        else:
            # No manual number (e.g. a plain "One-time seed building"): just number it.
            out.append(rf"\{cmd}{{{title}}}")

        i += 2 if drop_acl else 1

    latex = "\n".join(out)

    # Clean up inbound references to skipped sections BEFORE resolving the rest.
    # The target heading no longer exists, so "\ref{uqp:N}" would be undefined; strip
    # the reference (and its scaffolding) so the sentence still reads cleanly.
    for n in sorted(SKIP_SECTIONS):
        sec = r"\{\\textsection\}" + str(n) + r"(?![0-9.])"
        before = latex
        # "(see §N)" / "(§N)" parentheticals, and a preceding space.
        latex = re.sub(r"\s*\(see " + sec + r"\)", "", latex)
        latex = re.sub(r"\s*\(" + sec + r"\)", "", latex)
        # Trailing "See §N." sentence.
        latex = re.sub(r"\s*See " + sec + r"\.", "", latex)
        # Possessive "§N's X" -> "the X".
        latex = re.sub(sec + r"'s\b", "the", latex)
        # Anything left: drop the bare reference and warn so it can be reviewed.
        if re.search(sec, latex):
            print(
                f"WARNING: leftover reference to skipped §{n}; removed bare '§{n}'.",
                file=sys.stderr,
            )
            latex = re.sub(r"\s*" + sec, "", latex)
        if latex != before:
            print(f"  ✓ Cleaned inbound references to skipped §{n}")

    # Prose "§N" / "§N.M" (already unicode-converted to "{\textsection}N...") -> \ref.
    latex = re.sub(
        r"\{\\textsection\}(\d+(?:\.\d+)*[a-z]?)",
        lambda mm: r"{\textsection}\ref{uqp:" + mm.group(1) + "}",
        latex,
    )
    return latex


def convert_pilot_report(phase="1"):
    """
    Main conversion function. Converts specified phases and all prior phases.

    Args:
        phase (str): Which phase to run. "1" = §1 only, "2" = §1 + §2, etc.
    """
    phase_int = int(phase)
    print(
        f"[Phase {phase}] Converting datacenter_pilot_report.md → pilot_uq_usecase.tex"
    )
    print(f"Source: {SOURCE_MD}")
    print(f"Target: {TARGET_PATH}")

    # Check source file exists
    if not Path(SOURCE_MD).exists():
        print(f"ERROR: Source markdown file not found: {SOURCE_MD}", file=sys.stderr)
        return False

    # Read markdown
    try:
        with open(SOURCE_MD, "r", encoding="utf-8") as f:
            md_content = f.read()
        print(f"✓ Read markdown file ({len(md_content)} bytes)")
    except Exception as e:
        print(f"ERROR: Could not read markdown file: {e}", file=sys.stderr)
        return False

    # Build LaTeX output progressively
    latex_sections = []

    # Phase 1: Extract §1
    try:
        section_1 = extract_section_1(md_content)
        print(f"✓ Extracted §1 'What this is' ({len(section_1)} bytes)")
        latex_1 = generate_latex_section_1(section_1)
        latex_sections.append(latex_1)
        print(f"✓ Generated §1 LaTeX ({len(latex_1)} bytes)")
    except Exception as e:
        print(f"ERROR: Could not process §1: {e}", file=sys.stderr)
        return False

    # Phase 2: Extract §2 (if requested)
    if phase_int >= 2:
        try:
            table_md, footnote_md = extract_section_2(md_content)
            print(
                f"✓ Extracted §2 Glossary (table: {len(table_md)} bytes, footnotes: {len(footnote_md)} bytes)"
            )
            latex_2 = generate_latex_section_2(table_md, footnote_md)
            latex_sections.append(latex_2)
            print(f"✓ Generated §2 LaTeX ({len(latex_2)} bytes)")
        except Exception as e:
            print(f"ERROR: Could not process §2: {e}", file=sys.stderr)
            return False

    # Phase 3: Extract §3 (The system) with mermaid diagram (if requested)
    if phase_int >= 3:
        try:
            text_before, mermaid_code, caption, text_after, details_md = (
                extract_section_3(md_content)
            )
            print(
                f"✓ Extracted §3 'The system' (text: {len(text_before) + len(text_after)} bytes, diagram: {'yes' if mermaid_code else 'no'})"
            )
            latex_3 = generate_latex_section_3(
                text_before, mermaid_code, caption, text_after, details_md
            )
            latex_sections.append(latex_3)
            print(f"✓ Generated §3 LaTeX ({len(latex_3)} bytes)")
        except Exception as e:
            print(f"ERROR: Could not process §3: {e}", file=sys.stderr)
            return False

    # Phase 4: Extract §4 (Building one operating condition with mermaid)
    # and §4.1 (Commitment subsection) (if requested)
    if phase_int >= 4:
        # First extract and generate §4 (main section with diagram)
        try:
            text_before, mermaid_code, caption, text_after = extract_section_4(
                md_content
            )
            print(
                f"✓ Extracted §4 'Building condition' (diagram: {'yes' if mermaid_code else 'no'})"
            )
            latex_4 = generate_latex_section_4(
                text_before, mermaid_code, caption, text_after
            )
            latex_sections.append(latex_4)
            print(f"✓ Generated §4 LaTeX ({len(latex_4)} bytes)")
        except Exception as e:
            print(f"ERROR: Could not process §4: {e}", file=sys.stderr)
            return False

        # Then extract and generate §4.1 (subsection)
        try:
            before_table, table_md, after_table, details_md = extract_section_4_1(
                md_content
            )
            print(
                f"✓ Extracted §4.1 'Commitment' (text: {len(before_table) + len(after_table)} bytes, table: {'yes' if table_md else 'no'})"
            )
            latex_4_1 = generate_latex_section_4_1(
                before_table, table_md, after_table, details_md
            )
            latex_sections.append(latex_4_1)
            print(f"✓ Generated §4.1 LaTeX ({len(latex_4_1)} bytes)")
        except Exception as e:
            print(f"ERROR: Could not process §4.1: {e}", file=sys.stderr)
            return False

        # Extract and generate §4.2 (Dispatch) (if requested and not excluded)
        if "4.2" not in EXCLUDE_SECTIONS:
            try:
                text_before, mermaid_code, caption, text_after, details_md = (
                    extract_section_4_2(md_content)
                )
                print(
                    f"✓ Extracted §4.2 'Dispatch' (diagram: {'yes' if mermaid_code else 'no'})"
                )
                latex_4_2 = generate_latex_section_4_2(
                    text_before, mermaid_code, caption, text_after, details_md
                )
                latex_sections.append(latex_4_2)
                print(f"✓ Generated §4.2 LaTeX ({len(latex_4_2)} bytes)")
            except Exception as e:
                print(f"ERROR: Could not process §4.2: {e}", file=sys.stderr)
                return False

        # Extract and generate §4.2a (Reactive limits) (if requested and not excluded)
        if "4.2a" not in EXCLUDE_SECTIONS:
            try:
                text, details_md = extract_section_4_2a(md_content)
                if text:
                    print(f"✓ Extracted §4.2a 'Reactive limits' ({len(text)} bytes)")
                    latex_4_2a = generate_latex_section_4_2a(text, details_md)
                    latex_sections.append(latex_4_2a)
                    print(f"✓ Generated §4.2a LaTeX ({len(latex_4_2a)} bytes)")
            except Exception as e:
                print(f"ERROR: Could not process §4.2a: {e}", file=sys.stderr)
                return False

        # Extract and generate §4.2b (AC OPF) (if requested and not excluded)
        if "4.2b" not in EXCLUDE_SECTIONS:
            try:
                text, details_md = extract_section_4_2b(md_content)
                if text:
                    print(f"✓ Extracted §4.2b 'AC OPF' ({len(text)} bytes)")
                    latex_4_2b = generate_latex_section_4_2b(text, details_md)
                    latex_sections.append(latex_4_2b)
                    print(f"✓ Generated §4.2b LaTeX ({len(latex_4_2b)} bytes)")
            except Exception as e:
                print(f"ERROR: Could not process §4.2b: {e}", file=sys.stderr)
                return False

        # Extract and generate §4.3 (Wind becomes negative load)
        if True:
            try:
                impl_note, text, details_md = extract_section_4_3(md_content)
                if text or impl_note:
                    print(
                        f"✓ Extracted §4.3 'Wind becomes negative load' ({len(text)} bytes)"
                    )
                    latex_4_3 = generate_latex_section_4_3(impl_note, text, details_md)
                    latex_sections.append(latex_4_3)
                    print(f"✓ Generated §4.3 LaTeX ({len(latex_4_3)} bytes)")
            except Exception as e:
                print(f"ERROR: Could not process §4.3: {e}", file=sys.stderr)
                return False

        if True:
            try:
                impl_note, text = extract_section_4_4(md_content)
                if text or impl_note:
                    print(f"✓ Extracted §4.4 'The datacenter' ({len(text)} bytes)")
                    latex_4_4 = generate_latex_section_4_4(impl_note, text)
                    latex_sections.append(latex_4_4)
                    print(f"✓ Generated §4.4 LaTeX ({len(latex_4_4)} bytes)")
            except Exception as e:
                print(f"ERROR: Could not process §4.4: {e}", file=sys.stderr)
                return False

        # Sections 5-12: phase N includes every section with number <= N
        section_info = [
            (5, "CCT by bisection", "cct_bisection"),
            (6, "Parallel execution on Kestrel", "parallel_kestrel"),
            (7, "Results", "results"),
            (8, "First look at ACTIVSg2000", "first_look_activsg2000"),
            (9, "Where this goes", "where_this_goes"),
            (10, "GridKit end-user feedback", "gridkit_feedback"),
            (11, "Limitations", "limitations"),
            (12, "Implementation", "implementation"),
        ]

        for sect_num, sect_title, sect_label in section_info:
            if sect_num in SKIP_SECTIONS:
                print(f"⏭  Skipping §{sect_num} '{sect_title}' (in SKIP_SECTIONS)")
                continue
            if phase_int >= sect_num:
                try:
                    text = extract_section_n(md_content, sect_num)
                    if text:
                        print(
                            f"✓ Extracted §{sect_num} '{sect_title}' ({len(text)} bytes)"
                        )
                        latex_section = generate_latex_section_n(
                            sect_num, sect_title, text, sect_label
                        )
                        latex_sections.append(latex_section)
                        print(
                            f"✓ Generated §{sect_num} LaTeX ({len(latex_section)} bytes)"
                        )
                except Exception as e:
                    print(f"ERROR: Could not process §{sect_num}: {e}", file=sys.stderr)
                    return False

    # Combine all sections
    latex_output = "".join(latex_sections)

    # Switch every inner heading to a numbered command and resolve "§N.M" refs so
    # the block inherits the report's 9.1.x numbering (titlesec stretches starred
    # headings and gives them no number — see apply_numbered_headings_and_refs).
    latex_output = apply_numbered_headings_and_refs(latex_output)

    # Number down to subsubsection for this block, then restore the report default
    # (report class default is 2) so later chapters are unaffected.
    latex_output = (
        "% Auto-generated from datacenter_pilot_report.md — do not edit by hand.\n"
        "\\setcounter{secnumdepth}{3}% number subsubsections within this block\n\n"
        + latex_output
        + "\n\\setcounter{secnumdepth}{2}% restore report default\n"
    )

    # Note: no \bibliography here. This file is \input into a parent document
    # (the SciDAC report, or test_main.tex) which owns the bibliography.

    # Write to both output locations (SciDAC report and repo)
    output_paths = [TARGET_PATH, REPO_OUTPUT_PATH]

    for output_path in output_paths:
        target_dir = Path(output_path).parent
        try:
            target_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            print(
                f"ERROR: Could not create directory {target_dir}: {e}", file=sys.stderr
            )
            return False

        try:
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(latex_output)
            print(f"✓ Wrote LaTeX file: {output_path}")
        except Exception as e:
            print(
                f"ERROR: Could not write LaTeX file {output_path}: {e}", file=sys.stderr
            )
            return False

    print(f"\n[SUCCESS] Phase {phase} conversion complete!")
    print(f"Next steps:")
    print(f"  1. Review output: {TARGET_PATH}")
    print(f"  2. Or repo version: {REPO_OUTPUT_PATH}")
    if phase_int < 2:
        print(f"  2. Run again: python convert_pilot_report.py phase=2 (add Glossary)")
    elif phase_int < 3:
        print(
            f"  2. Run again: python convert_pilot_report.py phase=3 (add System with diagram)"
        )
    elif phase_int < 4:
        print(
            f"  2. Run again: python convert_pilot_report.py phase=4 (add Section 4 with diagram and subsections 4.1-4.2b)"
        )
    elif phase_int < 5:
        print(
            f"  2. Run again: python convert_pilot_report.py phase=5 (add Section 4.3 Wind)"
        )
    elif phase_int < 6:
        print(
            f"  2. Run again: python convert_pilot_report.py phase=6 (add Section 4.4 Datacenter)"
        )
    elif phase_int < 7:
        print(
            f"  2. Run again: python convert_pilot_report.py phase=7 (add Section 5 CCT)"
        )
    elif phase_int < 8:
        print(
            f"  2. Run again: python convert_pilot_report.py phase=8 (add Section 6 Parallel)"
        )
    elif phase_int < 9:
        print(
            f"  2. Run again: python convert_pilot_report.py phase=9 (add Section 7 Results)"
        )
    elif phase_int < 10:
        print(
            f"  2. Run again: python convert_pilot_report.py phase=10 (add Section 8 ACTIVSg2000)"
        )
    elif phase_int < 11:
        print(
            f"  2. Run again: python convert_pilot_report.py phase=11 (add Section 9 Where)"
        )
    elif phase_int < 12:
        print(
            f"  2. Run again: python convert_pilot_report.py phase=12 (add Section 10 Feedback)"
        )
    elif phase_int < 13:
        print(
            f"  2. Run again: python convert_pilot_report.py phase=13 (add Sections 11-12 Limitations & Implementation)"
        )
    else:
        print(
            f"  2. Test: bash ../chapters/uq-subsections/build_test.sh  (-> _build/test_uq_pilot.pdf)"
        )
    print(
        f"  3. Integrate: Add to uq.tex: \\input{{chapters/uq-subsections/pilot_uq_usecase}}"
    )
    print(f"\nFor more info, see: convert_pilot_report.md")

    # Copy BibTeX file to repo location for version control
    import shutil

    bibtex_source = "/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/tex-report/uq-usecase_refs.bib"
    # Keep every copy next to a .tex that references it in sync with the source,
    # including the copy the faithful harness (test_uq_pilot.tex) reads from the
    # report directory.
    bibtex_targets = [
        REPO_BIBTEX_PATH,
        str(Path(TARGET_PATH).parent / "uq-usecase_refs.bib"),
        str(Path(REPO_OUTPUT_PATH).parent / "uq-usecase_refs.bib"),
    ]
    for bibtex_target in bibtex_targets:
        try:
            Path(bibtex_target).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(bibtex_source, bibtex_target)
            print(f"✓ Synced bibliography: {bibtex_target}")
        except Exception as e:
            print(
                f"WARNING: Could not copy BibTeX file to {bibtex_target}: {e}",
                file=sys.stderr,
            )

    return True


if __name__ == "__main__":
    # Parse phase from command line if provided
    import sys

    phase = "12"  # Default: convert the whole report (§1–§12)
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            if arg.startswith("phase="):
                phase = arg.split("=")[1]
                break

    success = convert_pilot_report(phase=phase)
    sys.exit(0 if success else 1)
