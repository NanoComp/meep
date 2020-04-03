import meep as mp
import numpy as np
import autograd.numpy as npa
from autograd import grad, jacobian
from collections import namedtuple

Grid = namedtuple('Grid', ['x', 'y', 'z', 'w'])

class OptimizationProblem(object):
    """Top-level class in the MEEP adjoint module.

    Intended to be instantiated from user scripts with mandatory constructor
    input arguments specifying the data required to define an adjoint-based
    optimization.

    The class knows how to do one basic thing: Given an input vector
    of design variables, compute the objective function value (forward
    calculation) and optionally its gradient (adjoint calculation).
    This is done by the __call__ method. The actual computations
    are delegated to a hierarchy of lower-level classes, of which
    the uppermost is TimeStepper.

    """

    def __init__(self, 
                simulation,
                objective_function,
                objective_arguments,
                basis,
                fcen,
                df=0,
                nf=1,
                decay_dt=50,
                decay_fields=[mp.Ez],
                decay_by=1e-6
                 ):

        self.sim = simulation
        self.objective_function = objective_function
        self.objective_arguments = objective_arguments
        self.basis = basis
        self.design_regions = [dr.volume for dr in self.basis]
        self.num_bases = len(self.basis)
        self.f_bank = [] # objective function history
        
        self.fcen = fcen
        self.df = df
        self.nf = nf
        self.freq_min = self.fcen - self.df/2
        self.freq_max = self.fcen + self.df/2

        self.decay_by=decay_by
        self.decay_fields=decay_fields
        self.decay_dt=decay_dt

        self.num_design_params = [ni.num_design_params for ni in self.basis]
        
        # store sources for finite difference estimations
        self.forward_sources = self.sim.sources

        # The optimizer has three allowable states : "INIT", "FWD", and "ADJ".
        #    INIT - The optimizer is initialized and ready to run a forward simulation
        #    FWD  - The optimizer has already run a forward simulation
        #    ADJ  - The optimizer has already run an adjoint simulation (but not yet calculated the gradient)
        self.current_state = "INIT"

    def __call__(self, rho_vector=None, need_value=True, need_gradient=True):
        """Evaluate value and/or gradient of objective function.
        """
        if rho_vector is not None:
            self.update_design(rho_vector=rho_vector)

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
                self.adjoint_run()
                print("Calculating gradient...")
                self.calculate_gradient()
            elif self.current_state == "FWD":
                print("Starting adjoint run...")
                self.adjoint_run()
                print("Calculating gradient...")
                self.calculate_gradient()
            else:
                raise ValueError("Incorrect solver state detected: {}".format(self.current_state))
        
        return self.f0, self.gradient, self.design_grids

    def get_fdf_funcs(self):
        """construct callable functions for objective function value and gradient

        Returns
        -------
        2-tuple (f_func, df_func) of standalone (non-class-method) callables, where
            f_func(beta) = objective function value for design variables beta
           df_func(beta) = objective function gradient for design variables beta
        """

        def _f(x=None):
            (fq, _) = self.__call__(rho_vector = x, need_gradient = False)
            return fq

        def _df(x=None):
            (_, df) = self.__call__(need_value = False)
            return df

        return _f, _df
    
    def prepare_forward_run(self):
        # prepare forward run
        self.sim.reset_meep()

        # add forward sources
        self.sim.change_sources(self.forward_sources)

        # register user specified monitors
        for m in self.objective_arguments:
            m.register_monitors(self.fcen,self.df,self.nf)

        # register design region
        self.design_region_monitors = [self.sim.add_dft_fields([mp.Ex,mp.Ey,mp.Ez],self.fcen,self.df,self.nf,where=dr,yee_grid=False) for dr in self.design_regions]

        # store design region voxel parameters
        self.design_grids = [Grid(*self.sim.get_array_metadata(dft_cell=drm)) for drm in self.design_region_monitors]

    def forward_run(self):
        # set up monitors
        self.prepare_forward_run()

        # add monitor used to track dft convergence
        mdft = self.sim.add_dft_fields(self.decay_fields,self.fcen,self.df,1,center=self.design_regions[0].center,size=mp.Vector3(1/self.sim.resolution))

        # Forward run
        self.sim.run(until_after_sources=stop_when_dft_decayed(mdft, self.decay_dt, self.decay_fields, self.fcen, self.decay_by))

        # record objective quantities from user specified monitors
        self.results_list = []
        for m in self.objective_arguments:
            self.results_list.append(m())

        # evaluate objective
        self.f0 = self.objective_function(*self.results_list)

        # Store forward fields for each design basis in array (x,y,z,field_components,frequencies)
        # FIXME allow for multiple design regions
        self.d_E = [np.zeros((len(dg.x),len(dg.y),len(dg.z),3,self.nf),dtype=np.complex128) for dg in self.design_grids]
        for nb, dgm in enumerate(self.design_region_monitors):
            for f in range(self.nf):
                for ic, c in enumerate([mp.Ex,mp.Ey,mp.Ez]):
                    self.d_E[nb][:,:,:,ic,f] = np.atleast_3d(self.sim.get_dft_array(dgm,c,f))

        # store objective function evaluation in memory
        self.f_bank.append(self.f0)

        # update solver's current state
        self.current_state = "FWD"

    def adjoint_run(self):
        # Grab the simulation step size from the forward run
        self.dt = self.sim.fields.dt

        # Prepare adjoint run
        self.sim.reset_meep()

        # Replace sources with adjoint sources
        self.adjoint_sources = []
        for mi, m in enumerate(self.objective_arguments):
            dJ = jacobian(self.objective_function,mi)(*self.results_list) # get gradient of objective w.r.t. monitor
            self.adjoint_sources.append(m.place_adjoint_source(dJ,self.dt)) # place the appropriate adjoint sources
        self.sim.change_sources(self.adjoint_sources)

        # register design flux
        # TODO use yee grid directly 
        self.design_region_monitors = [self.sim.add_dft_fields([mp.Ex,mp.Ey,mp.Ez],self.fcen,self.df,self.nf,where=dr,yee_grid=False) for dr in self.design_regions]

        # add monitor used to track dft convergence
        mdft = self.sim.add_dft_fields(self.decay_fields,self.fcen,self.df,1,center=self.design_regions[0].center,size=mp.Vector3(1/self.sim.resolution))

        # Adjoint run
        self.sim.run(until_after_sources=stop_when_dft_decayed(mdft, self.decay_dt, self.decay_fields, self.fcen, self.decay_by))

        # Store adjoint fields for each design basis in array (x,y,z,field_components,frequencies)
        # FIXME allow for multiple design regions
        self.a_E = [np.zeros((len(dg.x),len(dg.y),len(dg.z),3,self.nf),dtype=np.complex128) for dg in self.design_grids]
        for nb, dgm in enumerate(self.design_region_monitors):
            for f in range(self.nf):
                for ic, c in enumerate([mp.Ex,mp.Ey,mp.Ez]):
                    self.a_E[nb][:,:,:,ic,f] = np.atleast_3d(self.sim.get_dft_array(dgm,c,f))
        
        # store frequencies (will be same for all monitors)
        self.frequencies = np.array(mp.get_flux_freqs(self.design_region_monitors[0]))

        # update optimizer's state
        self.current_state = "ADJ"

    def calculate_gradient(self):
        # Iterate through all design region bases and store gradient w.r.t. permittivity
        self.gradient = [2*np.sum(np.real(self.a_E[nb]*self.d_E[nb]),axis=(3)) for nb in range(self.num_bases)]
        # Return optimizer's state to initialization
        self.current_state = "INIT"
    
    def calculate_fd_gradient(self,num_gradients=1,db=1e-4,basis_idx=0):
        '''
        Estimate central difference gradients.
        '''

        if num_gradients > self.num_design_params[basis_idx]:
            raise ValueError("The requested number of gradients must be less than or equal to the total number of design parameters.")

        # cleanup simulation object
        self.sim.reset_meep()
        self.sim.change_sources(self.forward_sources)

        # preallocate result vector
        dummy_inputs = [np.zeros((self.nf,))]*len(self.objective_arguments)
        num_outputs = np.squeeze(self.objective_function(*dummy_inputs)).size
        fd_gradient = 0*np.ones((self.num_design_params[basis_idx],num_outputs))

        # randomly choose indices to loop estimate
        fd_gradient_idx = np.random.choice(self.num_design_params[basis_idx],num_gradients,replace=False)

        for k in fd_gradient_idx:
            
            b0 = np.ones((self.num_design_params[basis_idx],))
            b0[:] = self.basis[basis_idx].rho_vector
            # -------------------------------------------- #
            # left function evaluation
            # -------------------------------------------- #
            self.sim.reset_meep()
            
            # assign new design vector
            b0[k] -= db
            self.basis[basis_idx].set_rho_vector(b0)
            
            # initialize design monitors
            for m in self.objective_arguments:
                m.register_monitors(self.fcen,self.df,self.nf)
            
            # add monitor used to track dft convergence
            mdft = self.sim.add_dft_fields(self.decay_fields,self.fcen,self.df,1,center=self.design_regions[0].center,size=mp.Vector3(1/self.sim.resolution))
            self.sim.run(until_after_sources=stop_when_dft_decayed(mdft, self.decay_dt, self.decay_fields, self.fcen, self.decay_by))
            
            # record final objective function value
            results_list = []
            for m in self.objective_arguments:
                results_list.append(m())
            fm = self.objective_function(*results_list)

            # -------------------------------------------- #
            # right function evaluation
            # -------------------------------------------- #
            self.sim.reset_meep()

            # assign new design vector
            b0[k] += 2*db # central difference rule...
            self.basis[basis_idx].set_rho_vector(b0)

            # initialize design monitors
            for m in self.objective_arguments:
                m.register_monitors(self.fcen,self.df,self.nf)
            
            # add monitor used to track dft convergence
            mdft = self.sim.add_dft_fields(self.decay_fields,self.fcen,self.df,1,center=self.design_regions[0].center,size=mp.Vector3(1/self.sim.resolution))
            self.sim.run(until_after_sources=stop_when_dft_decayed(mdft, self.decay_dt, self.decay_fields, self.fcen, self.decay_by))
            
            # record final objective function value
            results_list = []
            for m in self.objective_arguments:
                results_list.append(m())
            fp = self.objective_function(*results_list)

            # -------------------------------------------- #
            # estimate derivative
            # -------------------------------------------- #
            fd_gradient[k,:] = (fp - fm) / (2*db)
        
        return np.squeeze(fd_gradient), fd_gradient_idx
    
    def update_design(self, rho_vector):
        """Update the design permittivity function.

        rho_vector ....... a list of numpy arrays that maps to each design region
        """
        for bi, b in enumerate(self.basis):
            b.set_rho_vector(rho_vector[bi])
        
        self.sim.reset_meep()
        self.current_state = "INIT"
    def get_objective_arguments(self):
        '''Return list of evaluated objective arguments.
        '''
        objective_args_evaluation = [m.get_evaluation() for m in self.objective_arguments]
        return objective_args_evaluation
        
    def plot2D(self,init_opt=False, **kwargs):
        """Produce a graphical visualization of the geometry and/or fields,
           as appropriately autodetermined based on the current state of
           progress.
        """

        if init_opt:
            self.prepare_forward_run()

        self.sim.plot2D(**kwargs)

def stop_when_dft_decayed(mon, dt, c, freq, decay_by):
    closure = {
        'previous_fields': np.array([1]*len(c)),
        't0': 0,
    }

    def _stop(sim):
        if sim.round_time() <= dt + closure['t0']:
            return False
        else:
            current_fields = np.array([np.abs(sim.get_dft_array(mon, ic, 0))[1] ** 2 for ic in c])
            previous_fields = closure['previous_fields']
            relative_change = np.abs(previous_fields - current_fields) / previous_fields
            closure['previous_fields'] = current_fields
            closure['t0'] = sim.round_time()
            idx = np.argmax(relative_change)
            if mp.cvar.verbosity > 0:
                fmt = "DFT decay(t = {0:1.1f}): abs({1:0.4e} - {2:0.4e}) / {3:0.4e} = {4:0.4e}"
                print(fmt.format(sim.meep_time(), previous_fields[idx], current_fields[idx], previous_fields[idx], relative_change[idx]))
            return relative_change[idx] <= decay_by
    return _stop