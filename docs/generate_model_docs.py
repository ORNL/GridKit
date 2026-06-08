#!/usr/bin/env python3
"""Generate Sphinx pages from selected GridKit Markdown files."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
MODEL = ROOT / "GridKit" / "Model"
EXAMPLES = ROOT / "examples"
GITHUB_REPO = "https://github.com/ORNL/GridKit"

GENERATED = DOCS / "generated"
GENERATED_MODELS = DOCS / "models" / "generated"
GENERATED_EXAMPLES = DOCS / "examples" / "generated"

COMMON_MATH = ROOT / "GridKit" / "CommonMath.md"
PHASOR_INPUT_FORMAT = MODEL / "PhasorDynamics" / "INPUT_FORMAT.md"

PROJECT_PAGES = {
    ROOT / "README.md": (GENERATED / "readme.md", "GridKit"),
    ROOT / "INSTALL.md": (GENERATED / "install.md", "Installation"),
    ROOT / "CONTRIBUTING.md": (GENERATED / "contributing.md", "Contributing"),
    ROOT / "CHANGELOG.md": (GENERATED / "changelog.md", "Changelog"),
    DOCS / "README.md": (GENERATED / "documentation-build.md", "Documentation Build"),
    ROOT / "buildsystem" / "README.md": (GENERATED / "buildsystem.md", "Buildsystem"),
    ROOT / "application" / "PhasorDynamics" / "README.md": (
        GENERATED / "application-input-format.md",
        "Application Input Format",
    ),
}

MODEL_READMES = tuple(path.resolve() for path in sorted(MODEL.glob("**/README.md")))
EXAMPLE_READMES = tuple(path.resolve() for path in sorted(EXAMPLES.glob("**/README.md")))


def git_output(*args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", *args],
            cwd=ROOT,
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def source_ref() -> str:
    if os.environ.get("READTHEDOCS_VERSION_TYPE") == "external":
        if ref := os.environ.get("READTHEDOCS_GIT_COMMIT_HASH"):
            return ref
    if ref := os.environ.get("READTHEDOCS_GIT_IDENTIFIER"):
        return ref
    if ref := git_output("rev-parse", "--abbrev-ref", "HEAD"):
        if ref != "HEAD":
            return ref
    return git_output("rev-parse", "HEAD") or "HEAD"


SOURCE_REF = source_ref()


def slug(text: str) -> str:
    text = re.sub(r"([a-z0-9])([A-Z])", r"\1-\2", text)
    text = re.sub(r"[^A-Za-z0-9.-]+", "-", text.replace("_", "-").replace("/", "-"))
    return re.sub(r"-+", "-", text).strip("-").lower()


def relative(target: Path, page: Path) -> str:
    return Path(os.path.relpath(target, page.parent)).as_posix()


def model_page(readme: Path) -> Path:
    parts = [slug(part) for part in readme.parent.relative_to(MODEL).parts]
    return GENERATED_MODELS.joinpath(*parts, "index.md")


def example_page(directory: Path) -> Path:
    parts = [slug(part) for part in directory.relative_to(EXAMPLES).parts]
    return GENERATED_EXAMPLES.joinpath(*parts, "index.md")


SOURCE_TO_PAGE = {source.resolve(): output for source, (output, _) in PROJECT_PAGES.items()}
SOURCE_TO_PAGE[COMMON_MATH.resolve()] = GENERATED_MODELS / "common-math.md"
SOURCE_TO_PAGE[PHASOR_INPUT_FORMAT.resolve()] = (
    GENERATED_MODELS / "phasor-dynamics" / "input-format.md"
)
SOURCE_TO_PAGE.update({readme: model_page(readme) for readme in MODEL_READMES})
SOURCE_TO_PAGE.update({readme: example_page(readme.parent) for readme in EXAMPLE_READMES})


def example_dirs() -> set[Path]:
    dirs = {EXAMPLES.resolve()}
    for readme in EXAMPLE_READMES:
        directory = readme.parent
        while directory != EXAMPLES.parent:
            dirs.add(directory.resolve())
            if directory == EXAMPLES:
                break
            directory = directory.parent
    return dirs


EXAMPLE_DOC_DIRS = example_dirs()


def title_for(source: Path) -> str:
    if source in PROJECT_PAGES:
        return PROJECT_PAGES[source][1]
    if source == COMMON_MATH:
        return "CommonMath"
    if source == PHASOR_INPUT_FORMAT:
        return "Input Format"
    if source in EXAMPLE_READMES and source.parent == EXAMPLES:
        return "Examples"
    return source.parent.name


def toctree(entries: list[str], maxdepth: int = 2) -> str:
    return (
        "```{toctree}\n"
        f":maxdepth: {maxdepth}\n"
        ":titlesonly:\n:hidden:\n\n"
        + "\n".join(entries)
        + "\n```\n"
    )


def children_of(directory: Path, readmes: tuple[Path, ...]) -> list[Path]:
    return sorted(readme.parent for readme in readmes if readme.parent.parent == directory)


def example_children(directory: Path) -> list[Path]:
    return sorted(path for path in EXAMPLE_DOC_DIRS if path.parent == directory and path != directory)


def generated_toctree(source: Path, page: Path) -> str:
    entries = []
    if source in MODEL_READMES:
        if source == MODEL / "PhasorDynamics" / "README.md":
            entries.append(relative(SOURCE_TO_PAGE[PHASOR_INPUT_FORMAT.resolve()].with_suffix(""), page))
        entries += [
            relative(model_page(child / "README.md").with_suffix(""), page)
            for child in children_of(source.parent, MODEL_READMES)
        ]
    elif source in EXAMPLE_READMES:
        entries = [
            relative(example_page(child).with_suffix(""), page)
            for child in example_children(source.parent)
        ]
    return f"{toctree(entries)}\n" if entries else ""


def split_anchor(target: str) -> tuple[str, str]:
    path, sep, anchor = target.partition("#")
    return path, f"{sep}{anchor}" if sep else ""


def is_external(target: str) -> bool:
    return bool(re.match(r"^[A-Za-z][A-Za-z0-9+.-]*:", target))


def fallback(source_dir: Path, target: str) -> Path | None:
    candidates = []
    if target.endswith("src/Model/PowerFlow/Gen/README.md"):
        candidates.append(MODEL / "PowerFlow" / "README.md")
    if target.startswith("src/Model/"):
        candidates.append(ROOT / "GridKit" / target.removeprefix("src/"))
    if target.startswith("Model/"):
        candidates.append(ROOT / "GridKit" / target)
    if "Model/" in target:
        candidates.append(MODEL / target.split("Model/", 1)[1])
    if source_dir.relative_to(ROOT).parts[:1] == ("application",):
        candidates.append(ROOT / "GridKit" / target.removeprefix("../../"))
    return next((path.resolve() for path in candidates if path.exists()), None)


def page_link(path: Path, anchor: str, current_page: Path) -> str | None:
    path = path.resolve()
    if path == COMMON_MATH.resolve() and anchor == "#anti-windup-indicator":
        anchor = "#derived-functions"
    if path in SOURCE_TO_PAGE:
        return f"{relative(SOURCE_TO_PAGE[path], current_page)}{anchor}"
    readme = (path / "README.md").resolve() if path.is_dir() else None
    if readme in SOURCE_TO_PAGE:
        return f"{relative(SOURCE_TO_PAGE[readme], current_page)}{anchor}"
    return None


def rewrite_link(source_dir: Path, target: str, current_page: Path) -> str:
    if not target or target.startswith("#") or is_external(target):
        return target
    path_text, anchor = split_anchor(target)
    path = (source_dir / path_text).resolve()
    if not path.exists():
        path = fallback(source_dir, path_text) or path
    link = page_link(path, anchor, current_page)
    if link:
        return link
    if path.exists() and path.is_dir():
        return f"{GITHUB_REPO}/tree/{SOURCE_REF}/{path.relative_to(ROOT).as_posix()}{anchor}"
    if path.exists() or ROOT in path.parents:
        return f"{relative(path, current_page)}{anchor}"
    return target


def rewrite_asset(source_dir: Path, target: str, current_page: Path) -> str:
    if not target or target.startswith("#") or is_external(target):
        return target
    path_text, anchor = split_anchor(target)
    return f"{relative((source_dir / path_text).resolve(), current_page)}{anchor}"


def rewrite_markdown_links(text: str, source_dir: Path, current_page: Path) -> str:
    text = re.sub(
        r"!\[([^\]]*)\]\(([^)\s]+)\)",
        lambda m: f"![{m.group(1)}]({rewrite_asset(source_dir, m.group(2), current_page)})",
        text,
    )
    return re.sub(
        r"(?<!!)\[([^\]]+)\]\(([^)\s]+)\)",
        lambda m: f"[{m.group(1)}]({rewrite_link(source_dir, m.group(2), current_page)})",
        text,
    )


def rewrite_html_images(text: str, source_dir: Path, current_page: Path) -> str:
    text = re.sub(r"(?im)^\s*<div\b[^>]*\balign=[\"']center[\"'][^>]*>\s*$", "", text)
    text = re.sub(r"(?im)^\s*</div>\s*$", "", text)

    def attr(attrs: str, name: str) -> str | None:
        match = re.search(rf"\b{name}\s*=\s*([\"'])(.*?)\1", attrs, re.IGNORECASE)
        return match.group(2) if match else None

    def image(match: re.Match[str]) -> str:
        before, target, after = match.groups()
        attrs = f"{before} {after}"
        options = []
        if alt := attr(attrs, "alt"):
            options.append(f":alt: {alt}")
        if (align := attr(attrs, "align")) in {"left", "center", "right"}:
            options.append(f":align: {align}")
        options_text = "\n".join(options)
        if options_text:
            options_text += "\n"
        target = rewrite_asset(source_dir, target, current_page)
        return f"\n```{{image}} {target}\n{options_text}```\n"

    return re.sub(
        r"<img\b([^>]*)\bsrc=[\"']([^\"']+)[\"']([^>]*)>",
        image,
        text,
        flags=re.IGNORECASE,
    )


def normalize_headings(text: str) -> str:
    levels = sorted({len(m.group(1)) for m in re.finditer(r"(?m)^(#{1,6})(\s+)", text)})
    levels = {level: index + 2 for index, level in enumerate(levels)}
    return re.sub(
        r"(?m)^(#{1,6})(\s+)",
        lambda m: f"{'#' * levels.get(len(m.group(1)), len(m.group(1)))}{m.group(2)}",
        text,
    )


def strip_title(text: str, title: str) -> str:
    def comparable(value: str) -> str:
        return re.sub(r"[^a-z0-9]+", "", value.replace("™", "").lower())

    match = re.match(r"\s*#\s+(.+?)\s*(?:\n+|$)", text)
    if match and comparable(match.group(1)) == comparable(title):
        return text[match.end() :]
    return text


def normalize_comment_fences(text: str, source: Path) -> str:
    if source != ROOT / "CONTRIBUTING.md":
        return text

    def fence(match: re.Match[str]) -> str:
        lines = [line.strip() for line in match.group(1).splitlines() if line.strip()]
        if lines and all(line.startswith(("/", "*")) for line in lines):
            return f"```text\n{match.group(1)}```"
        return match.group(0)

    return re.sub(r"```c\+\+\n(.*?)```", fence, text, flags=re.DOTALL)


def normalize_markdown(text: str, source: Path, title: str) -> str:
    text = strip_title(text, title)
    text = normalize_comment_fences(text, source)
    text = re.sub(r"(?m)^```\s*math\s*$", "```{math}", text)
    text = re.sub(r"\$`([^`]+)`\$", r"$\1$", text)
    text = normalize_headings(text)
    if source == COMMON_MATH:
        text = re.sub(
            r"(?m)^(#{2,6}\s+Derived Functions)",
            r"(anti-windup-indicator)=\n\1",
            text,
            count=1,
        )
    return text


def write_page(source: Path) -> None:
    source = source.resolve()
    page = SOURCE_TO_PAGE[source]
    page.parent.mkdir(parents=True, exist_ok=True)

    text = source.read_text(encoding="utf-8")
    if source == (ROOT / "README.md").resolve():
        text = rewrite_markdown_links(text, source.parent, page)
        page.write_text(text.rstrip() + "\n", encoding="utf-8")
        return

    title = title_for(source)
    text = normalize_markdown(text, source, title)
    text = rewrite_markdown_links(text, source.parent, page)
    text = rewrite_html_images(text, source.parent, page)
    source_name = source.relative_to(ROOT).as_posix()
    page.write_text(
        f"# {title}\n\n"
        f"{generated_toctree(source, page)}"
        f"_Source: `{source_name}`_\n\n"
        f"{text.rstrip()}\n",
        encoding="utf-8",
    )


def write_example_index(directory: Path) -> None:
    if (directory / "README.md").resolve() in EXAMPLE_READMES:
        return
    children = example_children(directory)
    if not children:
        return

    page = example_page(directory)
    page.parent.mkdir(parents=True, exist_ok=True)
    entries = [relative(example_page(child).with_suffix(""), page) for child in children]
    lines = ["## Sections", ""]
    for child in children:
        child_link = f"[{child.name}]({relative(example_page(child), page)})"
        grandchildren = [
            f"[{grandchild.name}]({relative(example_page(grandchild), page)})"
            for grandchild in example_children(child)
        ]
        lines.append(f"- {child_link}: {', '.join(grandchildren)}" if grandchildren else f"- {child_link}")

    page.write_text(
        f"# {'Examples' if directory == EXAMPLES else directory.name}\n\n"
        f"{toctree(entries)}\n"
        + "\n".join(lines)
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    for directory in (GENERATED, GENERATED_MODELS, GENERATED_EXAMPLES):
        if directory.exists():
            shutil.rmtree(directory)
        directory.mkdir(parents=True, exist_ok=True)

    sources = [
        *PROJECT_PAGES,
        COMMON_MATH,
        PHASOR_INPUT_FORMAT,
        *MODEL_READMES,
        *EXAMPLE_READMES,
    ]
    for source in sources:
        write_page(source)
    for directory in sorted(EXAMPLE_DOC_DIRS):
        write_example_index(directory)


if __name__ == "__main__":
    main()
