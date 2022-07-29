"""
Adjoint-based sensitivity-analysis module for pymeep.
Authors: Homer Reid <homer@homerreid.com>, Alec Hammond <alec.hammond@gatech.edu>, Ian Williamson <iwill@google.com>
"""

from . import utils
from .basis import BilinearInterpolationBasis
from .connectivity import *
from .filter_source import FilteredSource
from .filters import *
from .objective import *
from .optimization_problem import OptimizationProblem
from .utils import DesignRegion

try:
    from .wrapper import MeepJaxWrapper
except ModuleNotFoundError as _:
    pass
