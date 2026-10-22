"""nanobuilder: crystal + ligand nanostructure builder for nanoHUB.

The scientific backend never imports ipywidgets; the GUI lives in
``nanobuilder.gui`` and only calls the functions re-exported here.
"""

__version__ = "0.1.0"

from .schema import BuilderError, RunConfig  # noqa: E402,F401
