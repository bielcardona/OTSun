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
from .logging_unit import *
from .movements import *
from .__about__ import __version__
