#!/usr/bin/env python3
"""
Extract and generate Mermaid diagrams from markdown report.

This script finds all mermaid diagram blocks in the markdown document,
converts them to PNG figures (with fallback placeholders), and saves them
to the figures directory with consistent naming conventions.

Design: Idempotent. Each run regenerates all diagrams from scratch.
No state files. Same markdown always produces identical figures.

Usage:
    python extract_and_generate_mermaid.py              # Generate all diagrams
    python extract_and_generate_mermaid.py format=tikz  # Use TikZ (for PDF inclusion)
"""

import re
import sys
import subprocess
from pathlib import Path

# ============================================================================
# CONFIGURATION
# ============================================================================

SOURCE_MD = "/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/datacenter_pilot_report.md"
FIGURES_OUTPUT_DIR = (
    "/Users/isatkaus/projects/scidac/repos/gridkit/uq-usecase/work-notes/figures/"
)

# Output format: "png" or "tikz"
# PNG: Requires mermaid-cli (fallback to placeholder if unavailable)
# TikZ: Inline LaTeX code (embedded in conversion script)
OUTPUT_FORMAT = "png"

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def extract_caption_from_text(text, mermaid_start_pos):
    """
    Extract caption from text following a mermaid block.
    Looks for: *Figure N.N. Caption text*
    Returns: (caption_text, position) or (None, None)
    """
    # Search after mermaid block for italic caption
    caption_pattern = r"\*Figure\s+[\d.]+\.\s+([^*]+)\*"
    match = re.search(caption_pattern, text[mermaid_start_pos:])
    if match:
        caption_full = match.group(0)  # Includes *...*
        caption_text = match.group(1)  # Just the text
        return caption_full, caption_text
    return None, None


def generate_diagram_label(section_num, caption_text=None):
    """
    Generate consistent diagram label from section and caption.
    Examples:
        section=3, caption="System Overview" → "system_overview"
        section=3, caption=None → "diagram_3"

    Returns: Label suitable for filename (no spaces, lowercase, shell-safe)
    """
    if caption_text:
        # Limit to first ~40 chars (concise filenames)
        short_caption = caption_text.strip()[:40].strip()

        # Convert to slug: remove special chars, spaces→underscores
        slug = re.sub(r"[^a-z0-9\s_-]", "", short_caption.lower())
        slug = re.sub(r"[\s-]+", "_", slug)  # Spaces/hyphens → underscores
        slug = re.sub(r"_+", "_", slug)  # Remove consecutive underscores
        slug = slug.strip("_")  # Remove leading/trailing underscores

        return f"fig:{section_num}_{slug}"
    else:
        return f"fig:{section_num}"


def find_all_mermaid_blocks(md_content):
    """
    Find all mermaid diagram blocks in markdown.

    Returns: List of tuples:
        (section_num, diagram_label, mermaid_code, caption_text, full_caption_markdown)
    """
    diagrams = []

    # Find all section headers and their positions
    section_pattern = r"^##\s+(\d+)\.\s+(.+?)$"
    section_matches = list(re.finditer(section_pattern, md_content, re.MULTILINE))

    # For each section, find mermaid blocks
    mermaid_pattern = r"```mermaid\n(.*?)\n```"

    for i, section_match in enumerate(section_matches):
        section_num = section_match.group(1)
        section_start = section_match.end()

        # Find where next section starts (or end of document)
        section_end = (
            section_matches[i + 1].start()
            if i + 1 < len(section_matches)
            else len(md_content)
        )
        section_content = md_content[section_start:section_end]

        # Find all mermaid blocks in this section
        for mermaid_match in re.finditer(mermaid_pattern, section_content, re.DOTALL):
            mermaid_code = mermaid_match.group(1).strip()
            mermaid_pos = section_start + mermaid_match.start()

            # Extract caption
            caption_full, caption_text = extract_caption_from_text(
                md_content, mermaid_pos
            )

            # Generate label
            label = generate_diagram_label(section_num, caption_text)

            diagrams.append(
                (section_num, label, mermaid_code, caption_text, caption_full)
            )

    return diagrams


def mermaid_to_png(mermaid_code, output_path):
    """
    Convert mermaid code to PNG using mermaid-cli.

    Returns: True if successful, False otherwise.
    On failure, caller should create placeholder.
    """
    try:
        # Try to use mermaid-cli (mmdc)
        # --scale 2: 2x resolution for better quality
        # --width/--height: Larger rendering area for better layout
        result = subprocess.run(
            [
                "mmdc",
                "-i",
                "/dev/stdin",
                "-o",
                output_path,
                "-s",
                "2",  # 2x scale for better quality
                "-w",
                "1200",
                "-H",
                "1200",
            ],  # Larger rendering area
            input=mermaid_code,
            capture_output=True,
            text=True,
            timeout=10,
        )
        return result.returncode == 0
    except FileNotFoundError:
        print(f"  ⚠ mermaid-cli not found; will create placeholder", file=sys.stderr)
        return False
    except subprocess.TimeoutExpired:
        print(f"  ⚠ mermaid-cli timed out; will create placeholder", file=sys.stderr)
        return False
    except Exception as e:
        print(f"  ⚠ mermaid-cli error: {e}; will create placeholder", file=sys.stderr)
        return False


def create_placeholder_png(output_path, diagram_label, section_num):
    """
    Create a placeholder PNG when mermaid-cli is unavailable.
    Uses PIL if available, otherwise creates minimal PNG.
    """
    try:
        from PIL import Image, ImageDraw

        # Create placeholder image
        img = Image.new("RGB", (600, 400), color=(240, 240, 240))
        draw = ImageDraw.Draw(img)

        # Draw border
        draw.rectangle([(10, 10), (590, 390)], outline=(100, 100, 100), width=2)

        # Add text
        text_lines = [
            "[Diagram Placeholder]",
            f"Section {section_num}: {diagram_label}",
            "",
            "Full diagram generation requires:",
            "  npm install -g @mermaid-js/mermaid-cli",
        ]
        y = 50
        for line in text_lines:
            draw.text((30, y), line, fill=(50, 50, 50))
            y += 50

        img.save(output_path)
        return True
    except ImportError:
        print(
            f"  ⚠ PIL not available; cannot create placeholder PNG",
            file=sys.stderr,
        )
        print(f"  → Install with: pip install Pillow", file=sys.stderr)
        return False


def generate_all_mermaid_diagrams():
    """
    Main function: Extract all mermaid blocks and generate PNG files.

    Returns: Dict mapping diagram_label → filename
    """
    # Read markdown
    if not Path(SOURCE_MD).exists():
        print(f"ERROR: Source markdown not found: {SOURCE_MD}", file=sys.stderr)
        return {}

    with open(SOURCE_MD, "r", encoding="utf-8") as f:
        md_content = f.read()

    print(f"[Mermaid Extraction] Scanning: {SOURCE_MD}")
    print(f"Output directory: {FIGURES_OUTPUT_DIR}\n")

    # Create output directory
    Path(FIGURES_OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    # Find all mermaid blocks
    diagrams = find_all_mermaid_blocks(md_content)
    if not diagrams:
        print("No mermaid diagrams found.")
        return {}

    print(f"Found {len(diagrams)} mermaid diagram(s):\n")

    # Generate each diagram
    mapping = {}
    for section_num, label, code, caption, caption_full in diagrams:
        filename = f"{label}.png"
        filepath = Path(FIGURES_OUTPUT_DIR) / filename

        print(f"  [{section_num}] {label}")
        if caption:
            print(f"      Caption: {caption}")
        print(f"      → {filepath}")

        # Generate PNG
        if mermaid_to_png(code, str(filepath)):
            print(f"      ✓ Generated from mermaid-cli")
            mapping[label] = filename
        elif create_placeholder_png(str(filepath), label, section_num):
            print(f"      ✓ Created placeholder PNG")
            mapping[label] = filename
        else:
            print(f"      ✗ Failed to generate PNG; skipping", file=sys.stderr)

    print(f"\n[SUCCESS] Generated {len(mapping)}/{len(diagrams)} diagram(s)")
    return mapping


# ============================================================================
# MAIN
# ============================================================================


if __name__ == "__main__":
    # Parse command-line arguments
    for arg in sys.argv[1:]:
        if arg.startswith("format="):
            OUTPUT_FORMAT = arg.split("=", 1)[1]

    # Generate all diagrams
    mapping = generate_all_mermaid_diagrams()

    # Print mapping for reference
    if mapping:
        print("\nDiagram mapping (for reference):")
        for label, filename in sorted(mapping.items()):
            print(f"  {label} → {filename}")
