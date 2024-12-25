from __future__ import annotations
try:
    import freecad
    # Some distribuitions of FreeCAD need this import to set up the paths
except ImportError:
    pass
from .math import *
from .optics import *
from .outputs import *
from .materials import *
from .ray import *
from .source import *
from .scene import *
from .experiments import *
from .movements import *
from .helpers import *

import importlib.metadata
from pathlib import Path
import toml

import logging
logger = logging.getLogger(__name__)

try:
    __version__ = importlib.metadata.version('otsun')
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"
    # adopt path to your pyproject.toml
    pyproject_toml_file = Path(__file__).parent.parent.parent / "pyproject.toml"
    print(pyproject_toml_file)
    if pyproject_toml_file.exists() and pyproject_toml_file.is_file():
        data = toml.load(pyproject_toml_file)
        # check project.version
        if "project" in data and "version" in data["project"]:
            __version__ = data["project"]["version"]
        # check tool.poetry.version
        elif "tool" in data and "poetry" in data["tool"] and "version" in data["tool"]["poetry"]:
            __version__ = data["tool"]["poetry"]["version"]
