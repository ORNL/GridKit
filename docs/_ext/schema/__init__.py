from pathlib import Path

from sphinx.application import Sphinx
from sphinx.config import Config
from sphinx.errors import ExtensionError
from sphinx.util.typing import ExtensionMetadata

from .cases import CaseModelsDirective
from .generate import write_schema
from .inventory import InventoryError


def _generate(app: Sphinx, _config: Config) -> None:
    repository = Path(app.srcdir).parent
    try:
        write_schema(
            Path(app.srcdir) / "xml",
            repository / "GridKit" / "Model" / "case.schema.base.json",
            repository / "GridKit" / "Model" / "case.schema.json",
        )
    except (InventoryError, OSError, ValueError, KeyError) as error:
        raise ExtensionError(f"cannot generate case schema: {error}") from error


def setup(app: Sphinx) -> ExtensionMetadata:
    # Sphinx validates html_extra_path before the builder-inited event. Generate
    # the schema during configuration so a clean RTD checkout can copy it to the
    # stable /case.schema.json URL instead of only exposing a hashed download.
    app.connect("config-inited", _generate)
    app.add_directive("case-models", CaseModelsDirective)
    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
