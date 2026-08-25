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

from .unfilter_design import *

# JAX is an optional dependency; everything that needs it lives in `wrapper`.
# Importing it also registers JAX as a way to differentiate objective functions,
# so objective functions written with `jax.numpy` need no special treatment.
try:
    from .wrapper import MeepJaxWrapper, value_and_jacobian
except ModuleNotFoundError as _:
    pass
