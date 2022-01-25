import meep as mp
import numpy as np
from autograd import grad, jacobian
from collections import namedtuple

from . import utils, DesignRegion

class OptimizationProblem(object):
    """Top-level class in the MEEP adjoint module.

    Intended to be instantiated from user scripts with mandatory constructor
    input arguments specifying the data required to define an adjoint-based
    optimization.

    The class knows how to do one basic thing: Given an input vector
    of design variables, compute the objective function value (forward
    calculation) and optionally its gradient (adjoint calculation).
    This is done by the __call__ method.

    """
    def __init__(
        self,
        simulation,
        objective_functions,
        objective_arguments,
        design_regions,
        frequencies=None,
        fcen=None,
        df=None,
        nf=None,
        decay_by=1e-11,
        decimation_factor=0,
        minimum_run_time=0,
        maximum_run_time=None,
        finite_difference_step=utils.FD_DEFAULT
    ):

        self.sim = simulation

        if isinstance(objective_functions, list):
            self.objective_functions = objective_functions
        else:
            self.objective_functions = [objective_functions]
        self.objective_arguments = objective_arguments
        self.f_bank = []  # objective function evaluation history

        if isinstance(design_regions, list):
            self.design_regions = design_regions
        else:
            self.design_regions = [design_regions]

        self.num_design_params = [
            ni.num_design_params for ni in self.design_regions
        ]
        self.num_design_regions = len(self.design_regions)

        # TODO typecheck frequency choices
        if frequencies is not None:
            self.frequencies = frequencies
            self.nf = np.array(frequencies).size
        else:
            if nf == 1:
                self.nf = nf
                self.frequencies = [fcen]
            else:
                fmax = fcen + 0.5 * df
                fmin = fcen - 0.5 * df
                dfreq = (fmax - fmin) / (nf - 1)
                self.frequencies = np.linspace(
                    fmin,
                    fmin + dfreq * nf,
                    num=nf,
                    endpoint=False,
                )
                self.nf = nf

        if self.nf == 1:
            self.fcen_idx = 0
        else:
            self.fcen_idx = int(
                np.argmin(
                    np.abs(
                        np.asarray(self.frequencies) -
                        np.mean(np.asarray(self.frequencies)))**
                    2))  # index of center frequency

        self.decay_by = decay_by
        self.decimation_factor = decimation_factor
        self.minimum_run_time = minimum_run_time
        self.maximum_run_time = maximum_run_time
        self.finite_difference_step = finite_difference_step # step size used in Aᵤ computation

        # store sources for finite difference estimations
        self.forward_sources = self.sim.sources

        # The optimizer has three allowable states : "INIT", "FWD", and "ADJ".
        #    INIT - The optimizer is initialized and ready to run a forward simulation
        #    FWD  - The optimizer has already run a forward simulation
        #    ADJ  - The optimizer has already run an adjoint simulation (but not yet calculated the gradient)
        self.current_state = "INIT"

        self.gradient = []

    def __call__(self, rho_vector=None, need_value=True, need_gradient=True, beta=None):
        """Evaluate value and/or gradient of objective function.
        """
        if rho_vector:
            self.update_design(rho_vector=rho_vector, beta=beta)

        # Run forward run if requested
        if need_value and self.current_state == "INIT":
            print("Starting forward run...")
            self.forward_run()

        # Run adjoint simulation and calculate gradient if requested
        if need_gradient:
            if self.current_state == "INIT":
                # we need to run a forward run before an adjoint run
                print("Starting forward run...")
                self.forward_run()
                print("Starting adjoint run...")
                self.D_a = []
                self.adjoint_run()
                print("Calculating gradient...")
                self.calculate_gradient()
            elif self.current_state == "FWD":
                print("Starting adjoint run...")
                self.D_a = []
                self.adjoint_run()
                print("Calculating gradient...")
                self.calculate_gradient()
            else:
                raise ValueError("Incorrect solver state detected: {}".format(
                    self.current_state))

        return self.f0, self.gradient

    def get_fdf_funcs(self):
        """construct callable functions for objective function value and gradient

        Returns
        -------
        2-tuple (f_func, df_func) of standalone (non-class-method) callables, where
            f_func(beta) = objective function value for design variables beta
           df_func(beta) = objective function gradient for design variables beta
        """
        def _f(x=None):
            (fq, _) = self.__call__(rho_vector=x, need_gradient=False)
            return fq

        def _df(x=None):
            (_, df) = self.__call__(need_value=False)
            return df

        return _f, _df

    def prepare_forward_run(self):
        # prepare forward run
        self.sim.reset_meep()

        # add forward sources
        self.sim.change_sources(self.forward_sources)

        # register user specified monitors
        self.forward_monitors = [
        ]  # save reference to monitors so we can track dft convergence
        for m in self.objective_arguments:
            self.forward_monitors.append(m.register_monitors(self.frequencies))

        # register design region
        self.design_region_monitors = utils.install_design_region_monitors(
            self.sim, self.design_regions, self.frequencies, self.decimation_factor
        )

    def forward_run(self):
        # set up monitors
        self.prepare_forward_run()

        # Forward run
        self.sim.run(until_after_sources=mp.stop_when_dft_decayed(
            self.decay_by,
            self.minimum_run_time,
            self.maximum_run_time
        ))

        # record objective quantities from user specified monitors
        self.results_list = []
        for m in self.objective_arguments:
            self.results_list.append(m())

        # evaluate objectives
        self.f0 = [fi(*self.results_list) for fi in self.objective_functions]
        if len(self.f0) == 1:
            self.f0 = self.f0[0]

        # Store forward fields for each set of design variables in array
        self.D_f = utils.gather_design_region_fields(self.sim,self.design_region_monitors,self.frequencies)

        # store objective function evaluation in memory
        self.f_bank.append(self.f0)

        # update solver's current state
        self.current_state = "FWD"

    def prepare_adjoint_run(self):
        # Compute adjoint sources
        self.adjoint_sources = [[]
                                for i in range(len(self.objective_functions))]
        for ar in range(len(self.objective_functions)):
            for mi, m in enumerate(self.objective_arguments):
                dJ = jacobian(self.objective_functions[ar],
                              mi)(*self.results_list)
                # get gradient of objective w.r.t. monitor
                if np.any(dJ):
                    self.adjoint_sources[ar] += m.place_adjoint_source(
                        dJ)  # place the appropriate adjoint sources

    def adjoint_run(self):
        # set up adjoint sources and monitors
        self.prepare_adjoint_run()

        # flip the m number
        if utils._check_if_cylindrical(self.sim):
            self.sim.m = -self.sim.m

        # flip the k point
        if self.sim.k_point:
            self.sim.k_point *= -1

        for ar in range(len(self.objective_functions)):
            # Reset the fields
            self.sim.reset_meep()

            # Update the sources
            self.sim.change_sources(self.adjoint_sources[ar])

            # register design flux
            self.design_region_monitors = utils.install_design_region_monitors(
            self.sim, self.design_regions, self.frequencies, self.decimation_factor
            )

            # Adjoint run
            self.sim.run(until_after_sources=mp.stop_when_dft_decayed(
                self.decay_by,
                self.minimum_run_time,
                self.maximum_run_time
            ))

            # Store adjoint fields for each design set of design variables
            self.D_a.append(utils.gather_design_region_fields(self.sim,self.design_region_monitors,self.frequencies))

        # reset the m number
        if utils._check_if_cylindrical(self.sim):
            self.sim.m = -self.sim.m

        # reset the k point
        if self.sim.k_point: self.sim.k_point *= -1

        # update optimizer's state
        self.current_state = "ADJ"

    def calculate_gradient(self):
        # Iterate through all design regions and calculate gradient
        self.gradient = [[
            dr.get_gradient(
                self.sim,
                self.D_a[ar][dri],
                self.D_f[dri],
                self.frequencies,
                self.finite_difference_step
            ) for dri, dr in enumerate(self.design_regions)
        ] for ar in range(len(self.objective_functions))]

        # Cleanup list of lists
        if len(self.gradient) == 1:
            self.gradient = self.gradient[0]  # only one objective function
            if len(self.gradient) == 1:
                self.gradient = self.gradient[
                    0]  # only one objective function and one design region
        else:
            if len(self.gradient[0]) == 1:
                self.gradient = [
                    g[0] for g in self.gradient
                ]  # multiple objective functions bu one design region
        # Return optimizer's state to initialization
        self.current_state = "INIT"

    def calculate_fd_gradient(
        self,
        num_gradients=1,
        db=1e-4,
        design_variables_idx=0,
        filter=None,
    ):
        '''
        Estimate central difference gradients.

        Parameters
        ----------
        num_gradients ... : scalar
            number of gradients to estimate. Randomly sampled from parameters.
        db ... : scalar
            finite difference step size
        design_variables_idx ... : scalar
            which design region to pull design variables from

        Returns
        -----------
        fd_gradient ... : lists
            [number of objective functions][number of gradients]

        '''
        if filter is None:
            filter = lambda x: x
        if num_gradients < self.num_design_params[design_variables_idx]:
            # randomly choose indices to loop estimate
            fd_gradient_idx = np.random.choice(
                self.num_design_params[design_variables_idx],
                num_gradients,
                replace=False,
            )
        elif num_gradients == self.num_design_params[design_variables_idx]:
            fd_gradient_idx = range(self.num_design_params[design_variables_idx])
        else:
            raise ValueError(
                "The requested number of gradients must be less than or equal to the total number of design parameters."
            )

        # cleanup simulation object
        self.sim.reset_meep()
        self.sim.change_sources(self.forward_sources)

        # preallocate result vector
        fd_gradient = []

        for k in fd_gradient_idx:

            b0 = np.ones((self.num_design_params[design_variables_idx], ))
            b0[:] = (self.design_regions[design_variables_idx].
                     design_parameters.weights)
            # -------------------------------------------- #
            # left function evaluation
            # -------------------------------------------- #
            self.sim.reset_meep()

            # assign new design vector
            b0[k] -= db
            self.design_regions[design_variables_idx].update_design_parameters(
                b0)

            # initialize design monitors
            self.forward_monitors = []
            for m in self.objective_arguments:
                self.forward_monitors.append(
                    m.register_monitors(self.frequencies))

            self.sim.run(until_after_sources=mp.stop_when_dft_decayed(
                self.decay_by,
                self.minimum_run_time,
                self.maximum_run_time
            ))

            # record final objective function value
            results_list = []
            for m in self.objective_arguments:
                results_list.append(m())
            fm = [fi(*results_list) for fi in self.objective_functions]

            # -------------------------------------------- #
            # right function evaluation
            # -------------------------------------------- #
            self.sim.reset_meep()

            # assign new design vector
            b0[k] += 2 * db  # central difference rule...
            self.design_regions[design_variables_idx].update_design_parameters(
                b0)

            # initialize design monitors
            self.forward_monitors = []
            for m in self.objective_arguments:
                self.forward_monitors.append(
                    m.register_monitors(self.frequencies))

            # add monitor used to track dft convergence
            self.sim.run(until_after_sources=mp.stop_when_dft_decayed(
                self.decay_by,
                self.minimum_run_time,
                self.maximum_run_time
            ))

            # record final objective function value
            results_list = []
            for m in self.objective_arguments:
                results_list.append(m())
            fp = [fi(*results_list) for fi in self.objective_functions]

            # -------------------------------------------- #
            # estimate derivative
            # -------------------------------------------- #
            fd_gradient.append([
                np.squeeze((fp[fi] - fm[fi]) / (2 * db))
                for fi in range(len(self.objective_functions))
            ])

        # Cleanup singleton dimensions
        if len(fd_gradient) == 1:
            fd_gradient = fd_gradient[0]

        return fd_gradient, fd_gradient_idx

    def update_design(self, rho_vector, beta=None):
        """Update the design permittivity function.

        rho_vector ....... a list of numpy arrays that maps to each design region
        """
        for bi, b in enumerate(self.design_regions):
            if np.array(rho_vector[bi]).ndim > 1:
                raise ValueError(
                    "Each vector of design variables must contain only one dimension."
                )
            b.update_design_parameters(rho_vector[bi])
            if beta:
                b.update_beta(beta)

        self.sim.reset_meep()
        self.current_state = "INIT"

    def get_objective_arguments(self):
        '''Return list of evaluated objective arguments.
        '''
        objective_args_evaluation = [
            m.get_evaluation() for m in self.objective_arguments
        ]
        return objective_args_evaluation

    def plot2D(self, init_opt=False, **kwargs):
        """Produce a graphical visualization of the geometry and/or fields,
           as appropriately autodetermined based on the current state of
           progress.
        """

        if init_opt:
            self.prepare_forward_run()

        self.sim.plot2D(**kwargs)

def atleast_3d(*arys):
    from numpy import array, asanyarray, newaxis
    '''
    Modified version of numpy's `atleast_3d`

    Keeps one dimensional array data in first dimension, as
    opposed to moving it to the second dimension as numpy's
    version does. Keeps the meep dimensionality convention.

    View inputs as arrays with at least three dimensions.
    Parameters
    ----------
    arys1, arys2, ... : array_like
        One or more array-like sequences.  Non-array inputs are converted to
        arrays.  Arrays that already have three or more dimensions are
        preserved.
    Returns
    -------
    res1, res2, ... : ndarray
        An array, or list of arrays, each with ``a.ndim >= 3``.  Copies are
        avoided where possible, and views with three or more dimensions are
        returned.  For example, a 1-D array of shape ``(N,)`` becomes a view
        of shape ``(N, 1, 1)``, and a 2-D array of shape ``(M, N)`` becomes a
        view of shape ``(M, N, 1)``.
    '''
    res = []
    for ary in arys:
        ary = asanyarray(ary)
        if ary.ndim == 0:
            result = ary.reshape(1, 1, 1)
        elif ary.ndim == 1:
            result = ary[:, newaxis, newaxis]
        elif ary.ndim == 2:
            result = ary[:, :, newaxis]
        else:
            result = ary
        res.append(result)
    if len(res) == 1:
        return res[0]
    else:
        return res
