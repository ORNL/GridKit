#!/usr/bin/env python3
"""
Convert datacenter_pilot_report.md (Markdown) → pilot_uq_usecase.tex (LaTeX)

This script incrementally converts sections of the markdown pilot UQ report
into LaTeX format for inclusion in the SciDAC technical report.

Phase 1 (MVP): Converts §1 "What this is" (text only, no figures/tables).
Future phases will add tables (§2), figures/Mermaid (§3–7), math/footnotes (§5–10).

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
TARGET_PATH = "/Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections/pilot_uq_usecase.tex"

# Output figure directory (for final report)
FIGURES_OUTPUT_DIR = "/Users/isatkaus/projects/scidac/docs/report/SciDAC_OE_Technical_Report/chapters/uq-subsections/figs/"

# Phase 2+ options
INCLUDE_FOOTNOTES = True  # Set to False to omit cross-document references as footnotes
REPO_BRANCH = "isatkaus/uq-usecase"  # Git branch where details are located
REPO_BASE = "https://github.com/ORNL/GridKit"  # GitHub repository base URL

# Phase 3+ options (Mermaid diagrams)
MERMAID_FORMAT = "png"  # Options: "tikz" (inline, styled), "png" (requires mermaid-cli, creates placeholder if unavailable)

# Inline code formatting options
# INLINE_CODE_FORMAT = "math"   # `code` → $code$ (inline math mode, current behavior)
# INLINE_CODE_FORMAT = "texttt"  # `code` → \texttt{code} (monospace font, simpler, preserves spacing)
INLINE_CODE_FORMAT = "texttt"  # Choose preferred rendering: "math" or "texttt"

# Sections to exclude from conversion
# Useful for optional subsections that may not be needed in final report
# Example: EXCLUDE_SECTIONS = ["4.2a", "4.2b"]  # Skip reactive limits and AC OPF sections
EXCLUDE_SECTIONS = []  # Set to ["4.2a", "4.2b"] to exclude those sections

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

    def close_quote():
        nonlocal in_quote, quote_buffer
        if in_quote:
            processed_quote = "\n".join(quote_buffer)
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
            quote_buffer = []
            in_quote = False

    def close_list():
        nonlocal in_list, list_buffer
        if in_list:
            result_lines.append("\\begin{itemize}")
            for list_item in list_buffer:
                # Join multi-line items with spaces and apply formatting
                formatted_item = markdown_inline_formatting(list_item)
                result_lines.append(f"  \\item {formatted_item}")
            result_lines.append("\\end{itemize}")
            list_buffer = []
            in_list = False

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
                result_lines.append("\\begin{quote}")
            i += 1

        # Check for bullet list line (starts with * or -)
        elif (
            stripped
            and len(stripped) > 1
            and stripped[0] in ("*", "-")
            and stripped[1] == " "
        ):
            close_quote()
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
            # i is already incremented by while loop
            continue

        else:
            # Regular text or empty line
            close_quote()
            close_list()

            # Handle subsections: ### text → \subsubsection{text}
            if stripped.startswith("### "):
                subsection_text = stripped[4:].strip()
                result_lines.append(f"\\subsubsection{{{subsection_text}}}")
            elif stripped:  # Only add non-empty lines
                # Apply inline formatting to regular text lines
                formatted_line = markdown_inline_formatting(line)
                result_lines.append(formatted_line)
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

    # Wrap in LaTeX \subsection
    latex_output = r"""\subsection{Pilot Aleatoric UQ Usecase}
\label{sec:pilot_uq_usecase}

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
        # Simple compact format without borders or fancy styling
        latex = r"\vspace{0.15cm}" + "\n"
        latex += r"\centering\small" + "\n"
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
        latex += r"\vspace{0.15cm}" + "\n"
        latex += r"\normalsize\raggedright" + "\n"
        return latex

    # Determine column count and format
    n_cols = len(header)
    col_format = "p{0.25\\textwidth}p{0.7\\textwidth}"  # 2 columns, adjust if needed
    if n_cols != 2:
        col_format = " ".join(["l"] * n_cols)  # Fallback: left-aligned

    # Start with table wrapper (floating) or direct tabularx (inline)
    if float:
        latex = r"\begin{table}[htb]" + "\n"
        latex += r"  \small" + "\n"
        latex += r"  \renewcommand{\arraystretch}{1.4}" + "\n"
        latex += r"  \begin{tabularx}{\textwidth}{|p{0.35\textwidth}|X|}" + "\n"
    else:
        # Inline table: no float wrapper, just tabularx with spacing
        latex = r"\vspace{0.3cm}" + "\n"
        latex += r"\small" + "\n"
        latex += r"\renewcommand{\arraystretch}{1.4}" + "\n"
        latex += r"\begin{tabularx}{\textwidth}{|p{0.35\textwidth}|X|}" + "\n"

    latex += r"    \hline" + "\n"

    # Add header
    header_cells = []
    for cell in header:
        # Convert markdown formatting in header cells
        cell = markdown_inline_formatting(cell)
        header_cells.append(cell)
    latex += "    " + " & ".join(header_cells) + r" \\" + "\n"
    latex += r"    \hline" + "\n"

    # Add rows (no \hline after each row for cleaner spacing)
    for row in rows:
        row_cells = []
        for cell in row:
            # Convert markdown formatting in cells
            cell = markdown_inline_formatting(cell)
            row_cells.append(cell)
        latex += "    " + " & ".join(row_cells) + r" \\" + "\n"

    latex += r"    \hline" + "\n"
    latex += r"  \end{tabularx}" + "\n"

    if float:
        latex += r"\end{table}" + "\n"
    else:
        latex += r"\vspace{0.3cm}" + "\n"

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
    # Escape Unicode characters for LaTeX FIRST (before applying formatting)
    text = text.replace("§", r"{\textsection}")  # Section symbol
    text = text.replace("→", r"\textrightarrow{}")  # Right arrow
    text = text.replace("•", r"\textbullet{}")  # Bullet
    text = text.replace("–", "--")  # En-dash
    text = text.replace("—", "---")  # Em-dash

    # `code` → format based on INLINE_CODE_FORMAT setting
    # First, find all backtick-enclosed code blocks
    import re

    def convert_code(m):
        code = m.group(1)

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
            return r"\texttt{" + code + "}"

    text = re.sub(r"`([^`]+?)`", convert_code, text)

    # [text](url) → \url{url} (markdown links to LaTeX)
    # Handle both external URLs and file references
    text = re.sub(
        r"\[([^\]]+?)\]\(([^)]+?)\)",
        lambda m: r"\url{" + m.group(2) + "}",
        text,
    )

    # **bold** → \textbf{bold}
    text = re.sub(r"\*\*([^\*]+?)\*\*", r"\\textbf{\1}", text)

    # *italic* → \textit{italic}
    text = re.sub(r"(?<!\*)\*([^\*]+?)\*(?!\*)", r"\\textit{\1}", text)

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
        return f"\\textbf{{Implementation note:}} {text}"

    # Handle Details format (with links)
    if has_details:
        # Remove **Details:** prefix to extract just the links
        text = re.sub(r"\*\*Details:\*\*\s*", "", text)

        # Find all markdown links: [text](url)
        link_pattern = r"\[([^\]]+)\]\(([^)]+)\)"
        matches = list(re.finditer(link_pattern, text))

        if not matches:
            # No links found, just return the text with Details label
            return f"\\textbf{{Details:}} {text}"

        # Extract filename from FIRST link
        first_link_url = matches[0].group(2)
        filename = (
            first_link_url.split("#")[0] if "#" in first_link_url else first_link_url
        )

        # Build full GitHub URL (one URL for all sections from same file)
        github_url = f"{REPO_BASE}/blob/{REPO_BRANCH}/uq-usecase/work-notes/{filename}"

        # Extract all link texts (section names as they appear in markdown)
        section_names = [match.group(1) for match in matches]

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
        section=3, caption="One-time seed building (blue)..." → "fig:3_one_time_seed_building_blue_..."

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

        return f"fig:{section_num}_{slug}"
    else:
        return f"fig:{section_num}"


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


def generate_latex_figure_with_mermaid(
    mermaid_code, caption, fig_label="fig:system_overview", width=0.6
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

    # Copy to output directory
    png_filename = f"{fig_label}.png"
    Path(FIGURES_OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    png_output_path = Path(FIGURES_OUTPUT_DIR) / png_filename

    try:
        shutil.copy(str(png_source), str(png_output_path))
        print(f"  ✓ Copied figure: {png_output_path.name}")
    except Exception as e:
        print(f"WARNING: Could not copy figure: {e}", file=sys.stderr)
        return ""

    # Create LaTeX figure environment
    latex = r"\begin{figure}[htb]" + "\n"
    latex += r"  \centering" + "\n"
    latex += (
        r"  \includegraphics[width="
        + f"{width:.2f}"
        + r"\textwidth]{figs/"
        + png_filename
        + "}"
        + "\n"
    )

    # Add caption
    if caption:
        # Convert markdown *text* to LaTeX
        caption_text = caption.replace("*", "").strip()
        latex += r"  \caption{" + caption_text + "}" + "\n"

    latex += r"  \label{" + fig_label + "}" + "\n"
    latex += r"\end{figure}" + "\n"

    return latex


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
    latex_output = r"""\subsection{The system}
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
    latex_output = r"""\subsection{Glossary}
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
    latex_output = r"""\subsubsection{Commitment: a feasible band, not an economic optimum}
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
    latex_output = r"""\subsection{Building one operating condition}
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
    latex_output = r"\subsubsection*{0.0.4.2 Dispatch: DC OPF then AC PF}" + "\n"
    latex_output += (
        r"\addcontentsline{toc}{subsubsection}{0.0.4.2 Dispatch: DC OPF then AC PF}"
        + "\n"
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
    latex_output = (
        r"\subsubsection*{0.0.4.2a Reactive limits: PV to PQ switching}" + "\n"
    )
    latex_output += (
        r"\addcontentsline{toc}{subsubsection}{0.0.4.2a Reactive limits: PV to PQ switching}"
        + "\n"
    )
    latex_output += r"\label{sec:reactive_limits_pv_pq}" + "\n\n"
    latex_output += (
        r"\textcolor{red}{\emph{Optional section---may move to appendix or details later.}}"
        + "\n\n"
    )

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
    latex_output = r"\subsubsection*{0.0.4.2b AC OPF as a second Q-limit route}" + "\n"
    latex_output += (
        r"\addcontentsline{toc}{subsubsection}{0.0.4.2b AC OPF as a second Q-limit route}"
        + "\n"
    )
    latex_output += r"\label{sec:ac_opf_q_limits}" + "\n\n"
    latex_output += (
        r"\textcolor{red}{\emph{Optional section---may move to appendix or details later.}}"
        + "\n\n"
    )

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

    # Combine all sections
    latex_output = "".join(latex_sections)

    # Ensure target directory exists
    target_dir = Path(TARGET_PATH).parent
    try:
        target_dir.mkdir(parents=True, exist_ok=True)
        print(f"✓ Target directory ready: {target_dir}")
    except Exception as e:
        print(f"ERROR: Could not create target directory: {e}", file=sys.stderr)
        return False

    # Write LaTeX
    try:
        with open(TARGET_PATH, "w", encoding="utf-8") as f:
            f.write(latex_output)
        print(f"✓ Wrote LaTeX file: {TARGET_PATH}")
    except Exception as e:
        print(f"ERROR: Could not write LaTeX file: {e}", file=sys.stderr)
        return False

    print(f"\n[SUCCESS] Phase {phase} conversion complete!")
    print(f"Next steps:")
    print(f"  1. Review output: {TARGET_PATH}")
    if phase_int < 2:
        print(f"  2. Run again: python convert_pilot_report.py phase=2 (add Glossary)")
    elif phase_int < 3:
        print(
            f"  2. Run again: python convert_pilot_report.py phase=3 (add System with diagram)"
        )
    elif phase_int < 4:
        print(
            f"  2. Run again: python convert_pilot_report.py phase=4 (add Section 4 with diagram and subsection 4.1)"
        )
    else:
        print(f"  2. Test: cd ../chapters/uq-subsections && pdflatex test_main.tex")
    print(
        f"  3. Integrate: Add to uq.tex: \\input{{chapters/uq-subsections/pilot_uq_usecase}}"
    )
    print(f"\nFor more info, see: convert_pilot_report.md")

    return True


if __name__ == "__main__":
    # Parse phase from command line if provided
    import sys

    phase = "3"  # Default
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            if arg.startswith("phase="):
                phase = arg.split("=")[1]
                break

    success = convert_pilot_report(phase=phase)
    sys.exit(0 if success else 1)
