"""Sphinx extension providing GridKit model, case, and example documentation.

The C++ model documentation is the sole authority on what a model accepts, so
the JSON schema and every generated table are built from the same records.
"""

from __future__ import annotations

from pathlib import Path

from sphinx.application import Sphinx
from sphinx.config import Config
from sphinx.errors import ExtensionError
from sphinx.util.typing import ExtensionMetadata

from .directives import DIRECTIVES
from .doxygen import InventoryError
from .repository import RepositoryError, repository
from .schema import write_schema


def _generate(app: Sphinx, _config: Config) -> None:
    # Sphinx validates html_extra_path before the builder-inited event. Generate
    # the schema during configuration so a clean RTD checkout can copy it to the
    # stable /case.schema.json URL instead of only exposing a hashed download.
    models = Path(app.srcdir).parent / "GridKit" / "Model"
    try:
        write_schema(
            repository(app).models,
            models / "case.schema.base.json",
            models / "case.schema.json",
        )
    except (InventoryError, RepositoryError, OSError, ValueError, KeyError) as error:
        raise ExtensionError(f"cannot generate case schema: {error}") from error


def setup(app: Sphinx) -> ExtensionMetadata:
    app.connect("config-inited", _generate)
    for name, directive in DIRECTIVES.items():
        app.add_directive(name, directive)
    return {
        "version": "2.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
