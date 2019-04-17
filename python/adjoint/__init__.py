# Adjoint-based sensitivity-analysis submodule for the MEEP python module.
#
#Documentation:
#http://https://meep.readthedocs.io/en/latest/Python_Tutorials/AdjointSolver.md
#
#The code is structured as follows:
#
#    __init__.py:   	     formal definitions
#
#    OptimizationProblem.py:  abstract base class from which user-supplied
#                             classes describing specific optimization
#                             geometries should inherit; implements a high-level
#                             notion of a differentiable objective function
#                             depending on one or more input variables, which we
#                             call "objective quantities."
#
#    Objective.py             lower-level routines for carrying out MEEP timestepping
#                             calculations to compute objective quantities and their
#                             adjoint derivatives.
#
#    ParallelDesignTester.py: simple module to spawn a multiprocessor server
#                             pool to run large numbers of objective-function
#                             calculations in parallel
#
#    Basis.py:                general support for spatially-varying permittivity
#                             functions described by expansions in user-defined sets
#                             of basis functions, plus implementations of simple basis sets
#                             for some common cases
#
#    Visualization.py:        routines for visualizing MEEP input geometries
#                             and computational results.
#
#
#__all__ = [ 'OptimizationProblem', 'DFTCell', 'adjoint_options', 'update_plot',
#            'EHTransverse', 'Exyz', 'Hxyz', 'EHxyz',
#            'xHat', 'yHat', 'zHat', 'origin', 'FluxLine',
#            'visualize_sim', 'plot_basis' ]

from .OptimizationProblem import OptimizationProblem

from .Objective import (adjoint_options, xHat, yHat, zHat, origin,
                        EHTransverse, Exyz, Hxyz, EHxyz, Grid, xyzw2grid,
                        abs2, unit_vector, rel_diff, FluxLine,
                        DFTCell, ObjectiveFunction, AdjointSolver,
                        get_dft_cell_names)

from .ParallelDesignTester import ParallelDesignTester

from .Basis import (ParameterizedDielectric, PlaneWaveBasis,
                    FourierLegendreBasis, FiniteElementBasis,
                    SimpleFiniteElementBasis)

from .Visualization import (set_plot_default, plot_basis, AFEClient,
                            visualize_sim, AdjointVisualizer)
