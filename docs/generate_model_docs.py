#!/usr/bin/env python3
"""Generate Sphinx pages from GridKit Markdown documentation."""

from __future__ import annotations

import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MarkdownPage:
    source: Path
    output: Path
    title: str


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_DIR = REPO_ROOT / "docs"
MODEL_DIR = REPO_ROOT / "GridKit" / "Model"
EXAMPLES_DIR = REPO_ROOT / "examples"
GITHUB_TREE_URL = "https://github.com/ORNL/GridKit/tree/develop"

GENERAL_OUT_DIR = DOCS_DIR / "generated"
MODEL_OUT_DIR = DOCS_DIR / "models" / "generated"
EXAMPLE_OUT_DIR = DOCS_DIR / "examples" / "generated"
ROOT_INDEX = DOCS_DIR / "index.md"
ROOT_TOCTREE_ENTRIES = (
    "generated/install",
    "applications/index",
    "models/index",
    "examples/generated/index",
    "api",
    "development/index",
)

COMMON_MATH = REPO_ROOT / "GridKit" / "CommonMath.md"
PHASOR_INPUT_FORMAT = MODEL_DIR / "PhasorDynamics" / "INPUT_FORMAT.md"

PROJECT_DOCS = (
    MarkdownPage(REPO_ROOT / "README.md", ROOT_INDEX, "GridKit"),
    MarkdownPage(REPO_ROOT / "INSTALL.md", GENERAL_OUT_DIR / "install.md", "Installation"),
    MarkdownPage(
        REPO_ROOT / "CONTRIBUTING.md",
        GENERAL_OUT_DIR / "contributing.md",
        "Contributing",
    ),
    MarkdownPage(
        REPO_ROOT / "CHANGELOG.md",
        GENERAL_OUT_DIR / "changelog.md",
        "Changelog",
    ),
    MarkdownPage(
        DOCS_DIR / "README.md",
        GENERAL_OUT_DIR / "documentation-build.md",
        "Documentation Build",
    ),
    MarkdownPage(
        REPO_ROOT / "buildsystem" / "README.md",
        GENERAL_OUT_DIR / "buildsystem.md",
        "Buildsystem",
    ),
    MarkdownPage(
        REPO_ROOT / "application" / "PhasorDynamics" / "README.md",
        GENERAL_OUT_DIR / "application-input-format.md",
        "Application Input Format",
    ),
)

PROJECT_DOC_BY_SOURCE = {page.source: page for page in PROJECT_DOCS}
EXAMPLE_READMES = frozenset(EXAMPLES_DIR.glob("**/README.md"))


def example_doc_dirs() -> set[Path]:
    dirs = {EXAMPLES_DIR}
    for readme in EXAMPLE_READMES:
        directory = readme.parent
        while directory != EXAMPLES_DIR.parent:
            dirs.add(directory)
            if directory == EXAMPLES_DIR:
                break
            directory = directory.parent
    return dirs


EXAMPLE_DOC_DIRS = example_doc_dirs()


def slugify(value: str) -> str:
    value = value.replace("README.md", "").strip("/")
    value = re.sub(r"([a-z0-9])([A-Z])", r"\1-\2", value)
    value = value.replace("/", "-").replace("_", "-")
    value = re.sub(r"[^A-Za-z0-9.-]+", "-", value)
    value = re.sub(r"-+", "-", value)
    return value.strip("-").lower()


def title_for_dir(directory: Path, root: Path, root_title: str) -> str:
    if directory == root:
        return root_title
    return directory.name


def title_for(path: Path) -> str:
    if path in PROJECT_DOC_BY_SOURCE:
        return PROJECT_DOC_BY_SOURCE[path].title
    if path == COMMON_MATH:
        return "CommonMath"
    if path == PHASOR_INPUT_FORMAT:
        return "Input Format"
    if is_example_readme(path):
        return title_for_dir(path.parent, EXAMPLES_DIR, "Examples")
    return path.parent.name


def is_readme_under(path: Path, root: Path) -> bool:
    return path.name == "README.md" and root in path.parents


def model_parts_for_dir(model_dir: Path) -> list[str]:
    return [slugify(part) for part in model_dir.relative_to(MODEL_DIR).parts]


def example_parts_for_dir(example_dir: Path) -> list[str]:
    return [slugify(part) for part in example_dir.relative_to(EXAMPLES_DIR).parts]


def child_readme_dirs(model_dir: Path) -> list[Path]:
    return sorted(
        child
        for child in model_dir.iterdir()
        if child.is_dir() and (child / "README.md").exists()
    )


def example_child_dirs(directory: Path) -> list[Path]:
    return sorted(
        child
        for child in EXAMPLE_DOC_DIRS
        if child.parent == directory and child != directory
    )


def is_model_readme(path: Path) -> bool:
    return is_readme_under(path, MODEL_DIR)


def is_example_readme(path: Path) -> bool:
    return path in EXAMPLE_READMES


def is_model_container(path: Path) -> bool:
    return is_model_readme(path) and bool(child_readme_dirs(path.parent))


def is_example_container(path: Path) -> bool:
    return is_example_readme(path) and bool(example_child_dirs(path.parent))


def example_index_path_for_dir(directory: Path) -> Path:
    return EXAMPLE_OUT_DIR.joinpath(*example_parts_for_dir(directory), "index.md")


def example_readme_output_path(path: Path) -> Path:
    return example_index_path_for_dir(path.parent)


def output_path_for(path: Path) -> Path:
    if path in PROJECT_DOC_BY_SOURCE:
        return PROJECT_DOC_BY_SOURCE[path].output
    if path == COMMON_MATH:
        return MODEL_OUT_DIR / "common-math.md"
    if path == PHASOR_INPUT_FORMAT:
        return MODEL_OUT_DIR / "phasor-dynamics" / "input-format.md"
    if is_model_readme(path):
        return MODEL_OUT_DIR.joinpath(*model_parts_for_dir(path.parent), "index.md")
    if is_example_readme(path):
        return example_readme_output_path(path)
    raise ValueError(f"Unhandled documentation source: {path}")


def index_path_for(path: Path) -> Path:
    if is_model_readme(path):
        return MODEL_OUT_DIR.joinpath(*model_parts_for_dir(path.parent), "index.md")
    if is_example_readme(path):
        return example_index_path_for_dir(path.parent)
    return output_path_for(path)


def is_url(target: str) -> bool:
    return bool(re.match(r"^[A-Za-z][A-Za-z0-9+.-]*:", target))


def split_anchor(target: str) -> tuple[str, str]:
    path, sep, anchor = target.partition("#")
    return path, f"{sep}{anchor}" if sep else ""


def resolve_local(source_dir: Path, target: str) -> Path:
    return (source_dir / target).resolve()


def relpath_from_output(target: Path, current_output: Path) -> str:
    return Path(os.path.relpath(target, current_output.parent)).as_posix()


def repo_source_url(path: Path) -> str:
    relpath = path.relative_to(REPO_ROOT).as_posix()
    return f"{GITHUB_TREE_URL}/{relpath}"


def fallback_repo_path(source_dir: Path, path_part: str) -> Path | None:
    candidates: list[Path] = []
    if path_part.endswith("src/Model/PowerFlow/Gen/README.md"):
        candidates.append(MODEL_DIR / "PowerFlow" / "README.md")
    if path_part.startswith("src/Model/"):
        candidates.append(REPO_ROOT / "GridKit" / path_part.removeprefix("src/"))
    if path_part.startswith("Model/"):
        candidates.append(REPO_ROOT / "GridKit" / path_part)
    if "Model/" in path_part:
        candidates.append(
            REPO_ROOT / "GridKit" / "Model" / path_part.split("Model/", 1)[1]
        )

    source_parts = source_dir.relative_to(REPO_ROOT).parts
    if len(source_parts) >= 2 and source_parts[0] == "application":
        candidates.append(REPO_ROOT / "GridKit" / path_part.removeprefix("../../"))

    return next((candidate.resolve() for candidate in candidates if candidate.exists()), None)


def page_link(path: Path, anchor: str, current_output: Path) -> str | None:
    if path == COMMON_MATH and anchor == "#anti-windup-indicator":
        anchor = "#derived-functions"

    if path in PROJECT_DOC_BY_SOURCE or path == COMMON_MATH or path == PHASOR_INPUT_FORMAT:
        return f"{relpath_from_output(output_path_for(path), current_output)}{anchor}"
    if (is_model_readme(path) or is_example_readme(path)) and path.exists():
        return f"{relpath_from_output(output_path_for(path), current_output)}{anchor}"
    return None


def rewrite_target(source_dir: Path, target: str, current_output: Path) -> str:
    if not target or target.startswith("#") or is_url(target):
        return target

    path_part, anchor = split_anchor(target)
    resolved = resolve_local(source_dir, path_part)
    fallback_path = None
    if not resolved.exists():
        fallback_path = fallback_repo_path(source_dir, path_part)
        if fallback_path is not None:
            resolved = fallback_path

    generated_link = page_link(resolved, anchor, current_output)
    if generated_link is not None:
        return generated_link

    if resolved.exists() and not resolved.is_dir():
        return f"{relpath_from_output(resolved, current_output)}{anchor}"
    if resolved.exists() and resolved.is_dir():
        readme = resolved / "README.md"
        if readme.exists():
            generated_dir_link = page_link(readme, anchor, current_output)
            if generated_dir_link is not None:
                return generated_dir_link
        return f"{repo_source_url(resolved)}{anchor}"

    try:
        resolved.relative_to(REPO_ROOT)
    except ValueError:
        return target

    return f"{relpath_from_output(resolved, current_output)}{anchor}"


def rewrite_asset_target(source_dir: Path, target: str, current_output: Path) -> str:
    if not target or target.startswith("#") or is_url(target):
        return target

    path_part, anchor = split_anchor(target)
    resolved = resolve_local(source_dir, path_part)
    try:
        path = relpath_from_output(resolved, current_output)
    except ValueError:
        return target
    return f"{path}{anchor}"


def rewrite_markdown_links(text: str, source_dir: Path, current_output: Path) -> str:
    def image_repl(match: re.Match[str]) -> str:
        label, target = match.groups()
        return f"![{label}]({rewrite_asset_target(source_dir, target, current_output)})"

    def link_repl(match: re.Match[str]) -> str:
        label, target = match.groups()
        return f"[{label}]({rewrite_target(source_dir, target, current_output)})"

    text = re.sub(r"!\[([^\]]*)\]\(([^)\s]+)\)", image_repl, text)
    text = re.sub(r"(?<!!)\[([^\]]+)\]\(([^)\s]+)\)", link_repl, text)
    return text


def rewrite_html_paths(text: str, source_dir: Path, current_output: Path) -> str:
    text = re.sub(
        r"(?im)^\s*<div\b[^>]*\balign=[\"']center[\"'][^>]*>\s*$",
        "",
        text,
    )
    text = re.sub(r"(?im)^\s*</div>\s*$", "", text)

    def html_attr(attrs: str, name: str) -> str | None:
        match = re.search(rf"\b{name}\s*=\s*([\"'])(.*?)\1", attrs, re.IGNORECASE)
        return match.group(2) if match else None

    def repl(match: re.Match[str]) -> str:
        before, target, after = match.groups()
        attrs = f"{before} {after}"
        target = rewrite_asset_target(source_dir, target, current_output)
        options = []

        alt = html_attr(attrs, "alt")
        if alt:
            options.append(f":alt: {alt}")

        align = html_attr(attrs, "align")
        if align in {"left", "center", "right"}:
            options.append(f":align: {align}")

        option_block = "\n".join(options)
        if option_block:
            option_block += "\n"

        return f"\n```{{image}} {target}\n{option_block}```\n"

    return re.sub(
        r"<img\b([^>]*)\bsrc=[\"']([^\"']+)[\"']([^>]*)>",
        repl,
        text,
        flags=re.IGNORECASE,
    )


def normalize_heading_levels(text: str) -> str:
    levels = sorted(
        {len(match.group(1)) for match in re.finditer(r"(?m)^(#{1,6})(\s+)", text)}
    )
    mapping = {level: index + 2 for index, level in enumerate(levels)}

    def repl(match: re.Match[str]) -> str:
        level = len(match.group(1))
        return f"{'#' * mapping.get(level, level)}{match.group(2)}"

    return re.sub(r"(?m)^(#{1,6})(\s+)", repl, text)


def comparable_title(value: str) -> str:
    value = value.replace("™", "")
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def strip_duplicate_top_heading(text: str, title: str) -> str:
    match = re.match(r"\s*#\s+(.+?)\s*(?:\n+|$)", text)
    if match and comparable_title(match.group(1)) == comparable_title(title):
        return text[match.end() :]
    return text


def normalize_comment_fences(text: str, source: Path) -> str:
    if source != REPO_ROOT / "CONTRIBUTING.md":
        return text

    def repl(match: re.Match[str]) -> str:
        body = match.group(1)
        lines = [line.strip() for line in body.splitlines() if line.strip()]
        if lines and all(line.startswith(("/", "*")) for line in lines):
            return f"```text\n{body}```"
        return match.group(0)

    return re.sub(r"```c\+\+\n(.*?)```", repl, text, flags=re.DOTALL)


def normalize_markdown(text: str, source: Path, title: str) -> str:
    text = strip_duplicate_top_heading(text, title)
    text = normalize_comment_fences(text, source)
    text = re.sub(r"(?m)^```\s*math\s*$", "```{math}", text)
    text = re.sub(r"\$`([^`]+)`\$", r"$\1$", text)
    text = normalize_heading_levels(text)
    if source == COMMON_MATH:
        text = re.sub(
            r"(?m)^(#{2,6}\s+Derived Functions)",
            r"(anti-windup-indicator)=\n\1",
            text,
            count=1,
        )
    return text


def toctree_block(entries: list[str], *, hidden: bool = False, maxdepth: int = 2) -> str:
    options = [f":maxdepth: {maxdepth}", ":titlesonly:"]
    if hidden:
        options.append(":hidden:")

    return (
        "```{toctree}\n"
        + "\n".join(options)
        + "\n\n"
        + "\n".join(entries)
        + "\n```\n"
    )


def markdown_link(title: str, path: Path, current_output: Path) -> str:
    return f"[{title}]({relpath_from_output(path, current_output)})"


def tree_index_entries(source: Path, current_output: Path) -> list[str]:
    entries = []
    if source == MODEL_DIR / "PhasorDynamics" / "README.md":
        entries.append(
            relpath_from_output(
                output_path_for(PHASOR_INPUT_FORMAT).with_suffix(""),
                current_output,
            )
        )

    entries.extend(
        relpath_from_output(
            index_path_for(child_dir / "README.md").with_suffix(""),
            current_output,
        )
        for child_dir in child_readme_dirs(source.parent)
    )
    return entries


def example_index_entries(directory: Path, current_output: Path) -> list[str]:
    entries: list[str] = []

    entries.extend(
        relpath_from_output(example_index_path_for_dir(child).with_suffix(""), current_output)
        for child in example_child_dirs(directory)
    )
    return entries


def example_contents_block(directory: Path, current_output: Path) -> str:
    readme = directory / "README.md"
    child_dirs = example_child_dirs(directory)
    lines: list[str] = []

    if readme in EXAMPLE_READMES and child_dirs:
        lines.extend(
            [
                "## Pages",
                "",
                f"- {markdown_link('Overview', output_path_for(readme), current_output)}",
            ]
        )

    if child_dirs:
        if lines:
            lines.append("")
        lines.extend(["## Sections", ""])
        for child in child_dirs:
            child_link = markdown_link(
                title_for_dir(child, EXAMPLES_DIR, "Examples"),
                example_index_path_for_dir(child),
                current_output,
            )
            grandchild_links = [
                markdown_link(
                    title_for_dir(grandchild, EXAMPLES_DIR, "Examples"),
                    example_index_path_for_dir(grandchild),
                    current_output,
                )
                for grandchild in example_child_dirs(child)
            ]
            if grandchild_links:
                lines.append(f"- {child_link}: {', '.join(grandchild_links)}")
            else:
                lines.append(f"- {child_link}")

    return "\n".join(lines) + "\n"


def generated_page_title(source: Path) -> str:
    return title_for(source)


def generated_page_toctree(source: Path, current_output: Path) -> str:
    if source == REPO_ROOT / "README.md":
        return toctree_block(list(ROOT_TOCTREE_ENTRIES), hidden=True, maxdepth=4)
    if is_model_container(source):
        entries = tree_index_entries(source, current_output)
        return f"{toctree_block(entries, hidden=True)}\n" if entries else ""
    if is_example_container(source):
        entries = example_index_entries(source.parent, current_output)
        return f"{toctree_block(entries, hidden=True)}\n" if entries else ""
    return ""


def generate_page(source: Path) -> None:
    out = output_path_for(source)
    if source == REPO_ROOT / "README.md":
        out.write_text(
            "# GridKit\n\n"
            f"{generated_page_toctree(source, out)}"
            "```{include} ../README.md\n"
            ":relative-images:\n"
            "```\n",
            encoding="utf-8",
        )
        return

    source_dir = source.parent
    title = generated_page_title(source)
    body = source.read_text(encoding="utf-8")
    body = normalize_markdown(body, source, title)
    body = rewrite_markdown_links(body, source_dir, out)
    body = rewrite_html_paths(body, source_dir, out)

    out.parent.mkdir(parents=True, exist_ok=True)
    rel_source = source.relative_to(REPO_ROOT).as_posix()
    out.write_text(
        f"# {title}\n\n"
        f"{generated_page_toctree(source, out)}"
        f"_Source: `{rel_source}`_\n\n"
        f"{body.rstrip()}\n",
        encoding="utf-8",
    )


def generate_example_index_page(directory: Path) -> None:
    child_dirs = example_child_dirs(directory)
    readme = directory / "README.md"
    if readme in EXAMPLE_READMES:
        return
    if not child_dirs:
        return

    out = example_index_path_for_dir(directory)
    out.parent.mkdir(parents=True, exist_ok=True)
    entries = example_index_entries(directory, out)
    out.write_text(
        f"# {title_for_dir(directory, EXAMPLES_DIR, 'Examples')}\n\n"
        f"{toctree_block(entries, hidden=True)}\n"
        f"{example_contents_block(directory, out)}",
        encoding="utf-8",
    )


def main() -> None:
    for generated_dir in (GENERAL_OUT_DIR, MODEL_OUT_DIR, EXAMPLE_OUT_DIR):
        if generated_dir.exists():
            shutil.rmtree(generated_dir)
        generated_dir.mkdir(parents=True, exist_ok=True)

    project_sources = [page.source for page in PROJECT_DOCS]
    model_sources = [COMMON_MATH, PHASOR_INPUT_FORMAT]
    model_sources.extend(sorted(MODEL_DIR.glob("**/README.md")))
    example_sources = sorted(EXAMPLE_READMES)

    for source in project_sources + model_sources + example_sources:
        generate_page(source)
    for directory in sorted(EXAMPLE_DOC_DIRS):
        generate_example_index_page(directory)


if __name__ == "__main__":
    main()
