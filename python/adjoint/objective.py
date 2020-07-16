"""Handling of objective functions and objective quantities."""

from abc import ABC, abstractmethod
import numpy as np
import meep as mp
from .filter_source import FilteredSource
from .optimization_problem import YeeDims, atleast_3d

# ---------------------------------------------- #
# general-purpose constants and utility routines
# ---------------------------------------------- #
ORIGIN           = np.zeros(3)
XHAT, YHAT, ZHAT = [ mp.Vector3(a) for a in [[1.,0,0], [0,1.,0], [0,0,1.]] ]
E_CPTS           = [mp.Ex, mp.Ey, mp.Ez]
H_CPTS           = [mp.Hx, mp.Hy, mp.Hz]
EH_CPTS          = E_CPTS + H_CPTS
EH_TRANSVERSE    = [ [mp.Ey, mp.Ez, mp.Hy, mp.Hz],
                     [mp.Ez, mp.Ex, mp.Hz, mp.Hx],
                     [mp.Ex, mp.Ey, mp.Hx, mp.Hy] ]

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
        self.normal_direction = self.volume.swigobj.normal_direction()
        self.kpoint_func = kpoint_func
        self.eval = None
        self.EigenMode_kwargs = kwargs
        return

    def register_monitors(self,frequencies):
        self.frequencies = np.asarray(frequencies)
        self.monitor = self.sim.add_mode_monitor(frequencies,mp.FluxRegion(center=self.volume.center,size=self.volume.size),yee_grid=True)
        self.normal_direction = self.monitor.normal_direction
        return self.monitor

    def place_adjoint_source(self,dJ):
        '''Places an equivalent eigenmode monitor facing the opposite direction. Calculates the
        correct scaling/time profile.

        dJ ........ the user needs to pass the dJ/dMonitor evaluation
        dt ........ the timestep size from sim.fields.dt of the forward sim
        '''
        dt = self.sim.fields.dt
        T = self.sim.meep_time()

        dJ = np.atleast_1d(dJ)
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

        # -------------------------------------- #
        # Get scaling factor
        # -------------------------------------- #
        # leverage linearity and combine source for multiple frequencies
        if dJ.ndim == 2:
            dJ = np.sum(dJ,axis=1)

        # Determine the correct resolution scale factor
        if self.sim.cell_size.y == 0:
            dV = 1/self.sim.resolution
        elif self.sim.cell_size.z == 0:
            dV = 1/self.sim.resolution * 1/self.sim.resolution
        else:
            dV = 1/self.sim.resolution * 1/self.sim.resolution * 1/self.sim.resolution

        da_dE = 0.5 * self.cscale # scalar popping out of derivative
        iomega = (1.0 - np.exp(-1j * (2 * np.pi * self.frequencies) * dt)) * (1.0 / dt) # scaled frequency factor with discrete time derivative fix

        src = self.time_src

        # an ugly way to calcuate the scaled dtft of the forward source
        y = np.array([src.swigobj.current(t,dt) for t in np.arange(0,T,dt)]) # time domain signal
        fwd_dtft = np.matmul(np.exp(1j*2*np.pi*self.frequencies[:,np.newaxis]*np.arange(y.size)*dt), y)*dt/np.sqrt(2*np.pi) # dtft

        # we need to compensate for the phase added by the time envelope at our freq of interest
        src_center_dtft = np.matmul(np.exp(1j*2*np.pi*np.array([src.frequency])[:,np.newaxis]*np.arange(y.size)*dt), y)*dt/np.sqrt(2*np.pi)
        adj_src_phase = np.exp(1j*np.angle(src_center_dtft))

        if self.frequencies.size == 1:
            # Single frequency simulations. We need to drive it with a time profile.
            amp = da_dE * dV * dJ * iomega / fwd_dtft / adj_src_phase # final scale factor
        else:
            # multi frequency simulations
            scale = da_dE * dV * dJ * iomega / adj_src_phase
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

        return self.eval / self.dx
    def get_evaluation(self):
        '''Returns the requested eigenmode coefficient.
        '''
        try:
            return self.eval
        except AttributeError:
            raise RuntimeError("You must first run a forward simulation before resquesting an eigenmode coefficient.")
