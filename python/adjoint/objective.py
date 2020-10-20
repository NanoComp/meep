"""Handling of objective functions and objective quantities."""

from abc import ABC, abstractmethod
import numpy as np
import meep as mp
from .filter_source import FilteredSource
from .optimization_problem import atleast_3d, Grid

# Transverse component definitions needed for flux adjoint calculations
EH_components = [[mp.Ey, mp.Ez, mp.Hy, mp.Hz], [mp.Ez, mp.Ex, mp.Hz, mp.Hx], [mp.Ex, mp.Ey, mp.Hx, mp.Hy]]

class ObjectiveQuantitiy(ABC):
    @abstractmethod
    def __init__(self):
        return
    @abstractmethod
    def register_monitors(self):
        return
    @abstractmethod
    def place_adjoint_source(self):
        return
    @abstractmethod
    def __call__(self):
        return
    @abstractmethod
    def get_evaluation(self):
        return

class EigenmodeCoefficient(ObjectiveQuantitiy):
    def __init__(self,sim,volume,mode,forward=True,kpoint_func=None,**kwargs):
        '''
        '''
        self.sim = sim
        self.volume=volume
        self.mode=mode
        self.forward = 0 if forward else 1
        self.normal_direction = None
        self.kpoint_func = kpoint_func
        self.eval = None
        self.EigenMode_kwargs = kwargs
        return

    def register_monitors(self,frequencies):
        self.frequencies = np.asarray(frequencies)
        self.monitor = self.sim.add_mode_monitor(frequencies,mp.ModeRegion(center=self.volume.center,size=self.volume.size),yee_grid=True)
        self.normal_direction = self.monitor.normal_direction
        return self.monitor

    def place_adjoint_source(self,dJ):
        '''Places an equivalent eigenmode monitor facing the opposite direction. Calculates the
        correct scaling/time profile.
        dJ ........ the user needs to pass the dJ/dMonitor evaluation
        '''
        dJ = np.atleast_1d(dJ)
        dt = self.sim.fields.dt # the timestep size from sim.fields.dt of the forward sim
        # determine starting kpoint for reverse mode eigenmode source
        direction_scalar = 1 if self.forward else -1
        if self.kpoint_func is None:
            if self.normal_direction == 0:
                k0 = direction_scalar * mp.Vector3(x=1)
            elif self.normal_direction == 1:
                k0 = direction_scalar * mp.Vector3(y=1)
            elif self.normal_direction == 2:
                k0 == direction_scalar * mp.Vector3(z=1)
        else:
            k0 = direction_scalar * self.kpoint_func(self.time_src.frequency,1)
        if dJ.ndim == 2:
            dJ = np.sum(dJ,axis=1)
        da_dE = 0.5 * self.cscale # scalar popping out of derivative

        scale = adj_src_scale(self, dt)

        if self.frequencies.size == 1:
            # Single frequency simulations. We need to drive it with a time profile.
            amp = da_dE * dJ * scale # final scale factor
            src = self.time_src
        else:
            # multi frequency simulations
            scale = da_dE * dJ * scale
            src = FilteredSource(self.time_src.frequency,self.frequencies,scale,dt) # generate source from broadband response
            amp = 1

        # generate source object
        self.source = [mp.EigenModeSource(src,
                    eig_band=self.mode,
                    direction=mp.NO_DIRECTION,
                    eig_kpoint=k0,
                    amplitude=amp,
                    eig_match_freq=True,
                    size=self.volume.size,
                    center=self.volume.center,
                    **self.EigenMode_kwargs)]
        
        return self.source

    def __call__(self):
        # We just need a workable time profile, so just grab the first available time profile and use that.
        self.time_src = self.sim.sources[0].src

        # Eigenmode data
        direction = mp.NO_DIRECTION if self.kpoint_func else mp.AUTOMATIC
        ob = self.sim.get_eigenmode_coefficients(self.monitor,[self.mode],direction=direction,kpoint_func=self.kpoint_func,**self.EigenMode_kwargs)
        self.eval = np.squeeze(ob.alpha[:,:,self.forward]) # record eigenmode coefficients for scaling
        self.cscale = ob.cscale # pull scaling factor

        return self.eval
    def get_evaluation(self):
        '''Returns the requested eigenmode coefficient.
        '''
        try:
            return self.eval
        except AttributeError:
            raise RuntimeError("You must first run a forward simulation before resquesting an eigenmode coefficient.")


class FourierFields(ObjectiveQuantitiy):
    def __init__(self,sim,volume, component):
        self.sim = sim
        self.volume=volume
        self.eval = None
        self.component = component
        return

    def register_monitors(self,frequencies):
        self.frequencies = np.asarray(frequencies)
        self.num_freq = len(self.frequencies)
        self.monitor = self.sim.add_dft_fields([self.component], self.frequencies, where=self.volume, yee_grid=False)
        return self.monitor

    def place_adjoint_source(self,dJ):
        # Correctly format the dJ matrix
        print(dJ.shape)
        dJ_shape = np.array(self.weights.shape)
        dJ_shape[-1] = self.num_freq
        dJ = np.ascontiguousarray((dJ).reshape(dJ_shape))
        
        dt = self.sim.fields.dt # the timestep size from sim.fields.dt of the forward sim

        mon_dv = self.sim.resolution ** self.volume.get_nonzero_dims()
        time_scale = adj_src_scale(self, dt)
        amp = dJ*time_scale*mon_dv
        self.sources = []
        if self.component in [mp.Ex, mp.Ey, mp.Ez]:
            amp = -amp
        if self.num_freq == 1:
            self.sources += [mp.Source(self.time_src,component=self.component,amp_data=amp[:,:,:,0],
                            center=self.volume.center,size=self.volume.size,amp_data_use_grid=True)]
        else:
            src = FilteredSource(self.time_src.frequency,self.frequencies,amp.reshape(-1, *amp.shape[-1:]),dt)
            (num_basis, num_pts) = src.nodes.shape
            fit_data = src.nodes
            fit_data = (fit_data.T)[:,np.newaxis,np.newaxis,:]#fit_data.reshape(dJ_shape)
            for basis_i in range(num_basis):
                self.sources += [mp.Source(src.time_src_bf[basis_i],component=self.component,amp_data=np.ascontiguousarray(fit_data[:,:,:,basis_i]),
                            center=self.volume.center,size=self.volume.size,amp_data_use_grid=True)]
        return self.sources 

    def __call__(self):
        self.time_src = self.sim.sources[0].src
        self.eval = np.array([self.sim.get_dft_array(self.monitor, self.component, i) for i in range(self.num_freq)]) #Shape = (num_freq, [pts])
        
        # Calculate, shape, and store weights
        self.dg = Grid(*self.sim.get_array_metadata(dft_cell=self.monitor))
        self.expected_dims = np.array([len(self.dg.x),len(self.dg.y),len(self.dg.z)])
        self.weights = self.dg.w.reshape(self.expected_dims)
        self.weights = self.weights[...,np.newaxis]
        
        return self.eval

    def get_evaluation(self):
        try:
            return self.eval
        except AttributeError:
            raise RuntimeError("You must first run a forward simulation.")


class Near2FarFields(ObjectiveQuantitiy):
    def __init__(self,sim,Near2FarRegions, far_pt):
        self.sim = sim
        self.Near2FarRegions=Near2FarRegions
        self.eval = None
        self.far_pt = far_pt
        return

    def register_monitors(self,frequencies):
        self.frequencies = np.asarray(frequencies)
        self.num_freq = len(self.frequencies)
        self.monitor = self.sim.add_near2far(self.frequencies, *self.Near2FarRegions, yee_grid=True)
        return self.monitor

    def place_adjoint_source(self,dJ):
        dt = self.sim.fields.dt # the timestep size from sim.fields.dt of the forward sim
        self.sources = []
        dJ = dJ.flatten()

        #TODO far_pts in 3d or cylindrical, perhaps py_v3_to_vec from simulation.py
        self.all_nearsrcdata = self.monitor.swigobj.near_sourcedata(mp.vec(self.far_pt.x, self.far_pt.y), dJ)
        for near_data in self.all_nearsrcdata:
            cur_comp = near_data.near_fd_comp
            amp_arr = np.array(near_data.amp_arr).reshape(-1, self.num_freq)
            scale = amp_arr * adj_src_scale(self, dt, include_resolution=False)

            if self.num_freq == 1:
                self.sources += [mp.IndexedSource(self.time_src, near_data, scale[:,0])]
            else:
                src = FilteredSource(self.time_src.frequency,self.frequencies,scale,dt)
                (num_basis, num_pts) = src.nodes.shape
                for basis_i in range(num_basis):
                    self.sources += [mp.IndexedSource(src.time_src_bf[basis_i], near_data, src.nodes[basis_i])]

        return self.sources

    def __call__(self):
        self.time_src = self.sim.sources[0].src
        self.eval = np.array(self.sim.get_farfield(self.monitor, self.far_pt))
        self.eval = self.eval.reshape((self.num_freq, 6))
        return self.eval

    def get_evaluation(self):
        try:
            return self.eval
        except AttributeError:
            raise RuntimeError("You must first run a forward simulation.")

class PoyntingFlux(ObjectiveQuantitiy):
    def __init__(self,sim,volume):
        self.sim = sim
        self.volume=volume
        self.eval = None
        self.normal_direction = None
        return

    def register_monitors(self,frequencies):
        self.frequencies = np.asarray(frequencies)
        self.num_freq = len(self.frequencies)
        self.monitor = self.sim.add_flux(self.frequencies, mp.FluxRegion(center=self.volume.center,size=self.volume.size))
        self.normal_direction = self.monitor.normal_direction
        return self.monitor

    def place_adjoint_source(self,dJ):        
        
        # Format jacobian (and use linearity)
        if dJ.ndim == 2:
            dJ = np.sum(dJ,axis=1)
        dJ = dJ.reshape(1,1,1,self.num_freq,1)
        
        # Adjust for time scaling
        dt = self.sim.fields.dt
        time_scale = adj_src_scale(self, dt).reshape(1,1,1,self.num_freq,1)
        
        # final source amplitude as a function of position
        amp = dJ * time_scale * np.conj(self.m_EH) * self.weights * (self.sim.resolution**self.volume.get_nonzero_dims())
        
        self.sources = []
        for ic,c in enumerate(EH_components[self.normal_direction]):
            # principle of equivalence
            amp_scale = -1 if c in [mp.Hx,mp.Hy,mp.Hz] else 1
            
            if np.any(amp[:,:,:,:,ic]):
                # single frequency case
                if self.num_freq == 1:
                    self.sources += [mp.Source(self.time_src,component=EH_components[self.normal_direction][3-ic],amp_data=np.ascontiguousarray(amp[:,:,:,0,ic]),
                                center=self.volume.center,size=self.volume.size,amp_data_use_grid=True,amplitude=amp_scale)]
                # multi frequency case
                else:
                    src = FilteredSource(self.time_src.frequency,self.frequencies,amp[...,ic].reshape(-1, self.num_freq),dt) # fit data to basis functions
                    fit_data = (src.nodes.T).reshape(amp[...,ic].shape) # format new amplitudes
                    for basis_i in range(self.num_freq):
                        self.sources += [mp.Source(src.time_src_bf[basis_i],EH_components[self.normal_direction][3-ic],amp_data=np.ascontiguousarray(fit_data[:,:,:,basis_i]),
                                    center=self.volume.center,size=self.volume.size,amp_data_use_grid=True,amplitude=amp_scale)]
        return self.sources 

    def __call__(self):
        self.time_src = self.sim.sources[0].src

        # Evaluate poynting flux
        self.eval = np.array(mp.get_fluxes(self.monitor))

        # Calculate, shape, and store weights
        self.dg = Grid(*self.sim.get_array_metadata(dft_cell=self.monitor))
        self.weights = self.dg.w.reshape((len(self.dg.x),len(self.dg.y),len(self.dg.z),1,1))

        # Store fields at monitor
        self.m_EH = np.zeros((len(self.dg.x), len(self.dg.y), len(self.dg.z), self.num_freq, 4), dtype=complex)
        for f in range(self.num_freq):
            for ic, c in enumerate(EH_components[self.normal_direction]):
                self.m_EH[..., f, ic] = self.sim.get_dft_array(self.monitor, c, f).reshape(len(self.dg.x), len(self.dg.y), len(self.dg.z))
        
        return self.eval

    def get_evaluation(self):
        try:
            return self.eval
        except AttributeError:
            raise RuntimeError("You must first run a forward simulation.")

def adj_src_scale(obj_quantity, dt, include_resolution=True):
    # -------------------------------------- #
    # Get scaling factor
    # -------------------------------------- #
    # leverage linearity and combine source for multiple frequencies
    T = obj_quantity.sim.meep_time()

    '''
    Integral-like adjoints (e.g. poynting flux, mode overlaps, etc)
    have a dV scale factor that needs to be included in the adjoint
    source. Other quantities (e.g. Near2Far, DFT fields, etc.) don't
    need the scale factor since no integral is involved.
    '''
    if not include_resolution:
        dV = 1
    elif obj_quantity.sim.cell_size.y == 0:
        dV = 1/obj_quantity.sim.resolution
    elif obj_quantity.sim.cell_size.z == 0:
        dV = 1/obj_quantity.sim.resolution * 1/obj_quantity.sim.resolution
    else:
        dV = 1/obj_quantity.sim.resolution * 1/obj_quantity.sim.resolution * 1/obj_quantity.sim.resolution

    iomega = (1.0 - np.exp(-1j * (2 * np.pi * obj_quantity.frequencies) * dt)) * (1.0 / dt) # scaled frequency factor with discrete time derivative fix

    src = obj_quantity.time_src

    # an ugly way to calcuate the scaled dtft of the forward source
    y = np.array([src.swigobj.current(t,dt) for t in np.arange(0,T,dt)]) # time domain signal
    fwd_dtft = np.matmul(np.exp(1j*2*np.pi*obj_quantity.frequencies[:,np.newaxis]*np.arange(y.size)*dt), y)*dt/np.sqrt(2*np.pi) # dtft

    # we need to compensate for the phase added by the time envelope at our freq of interest
    src_center_dtft = np.matmul(np.exp(1j*2*np.pi*np.array([src.frequency])[:,np.newaxis]*np.arange(y.size)*dt), y)*dt/np.sqrt(2*np.pi)
    adj_src_phase = np.exp(1j*np.angle(src_center_dtft))

    if obj_quantity.frequencies.size == 1:
        # Single frequency simulations. We need to drive it with a time profile.
        scale = dV * iomega / fwd_dtft / adj_src_phase # final scale factor
    else:
        # multi frequency simulations
        scale = dV * iomega / adj_src_phase
    return scale
