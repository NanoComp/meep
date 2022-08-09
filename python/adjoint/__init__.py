"""
Adjoint-based sensitivity-analysis module for pymeep.
Authors: Homer Reid <homer@homerreid.com>, Alec Hammond <alec.hammond@gatech.edu>, Ian Williamson <iwill@google.com>
"""

from .objective import *

from . import utils
from .utils import DesignRegion

from .basis import BilinearInterpolationBasis

from .optimization_problem import OptimizationProblem

from .filter_source import FilteredSource

from .filters import *

from .connectivity import *

try:
    from .wrapper import MeepJaxWrapper
except ModuleNotFoundError as _:
    pass
