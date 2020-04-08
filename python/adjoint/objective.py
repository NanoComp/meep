"""Handling of objective functions and objective quantities."""

from abc import ABC, abstractmethod
import numpy as np
import meep as mp
from .filter_source import FilteredSource
from matplotlib import pyplot as plt

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
    def __init__(self,sim,volume,mode,forward=True,k0=None,**kwargs):
        '''
        '''
        self.sim = sim
        self.volume=volume
        self.mode=mode
        self.forward = 0 if forward else 1
        self.normal_direction = None
        self.k0 = k0
        self.eval = None
        self.EigenMode_kwargs = kwargs
        return
    
    def register_monitors(self,frequencies):
        self.frequencies = np.asarray(frequencies)
        self.monitor = self.sim.add_mode_monitor(frequencies,mp.FluxRegion(center=self.volume.center,size=self.volume.size))
        self.normal_direction = self.monitor.normal_direction
        return self.monitor
    
    def place_adjoint_source(self,dJ,dt):
        '''Places an equivalent eigenmode monitor facing the opposite direction. Calculates the 
        correct scaling/time profile.

        dJ ........ the user needs to pass the dJ/dMonitor evaluation
        dt ........ the timestep size from sim.fields.dt of the forward sim
        '''
        dJ = np.atleast_1d(dJ)
        # determine starting kpoint for reverse mode eigenmode source
        direction_scalar = 1 if self.forward else -1
        if self.k0 is None:
            if self.normal_direction == 0:
                k0 = direction_scalar * mp.Vector3(x=1)
            elif self.normal_direction == 1:
                k0 = direction_scalar * mp.Vector3(y=1)
            elif self.normal_direction == 2:
                k0 == direction_scalar * mp.Vector3(z=1)
        else:
            k0 = direction_scalar * self.k0
        
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
        da_dE = 0.5*(dV * self.cscale)
        scale = da_dE * dJ * 1j * 2 * np.pi * self.frequencies / np.array([self.time_src.fourier_transform(f) for f in self.frequencies]) # final scale factor
        if self.frequencies.size == 1:
            # Single frequency simulations. We need to drive it with a time profile.
            src = self.time_src
            amp = scale
        else:
            # TODO: In theory we should be able drive the source without normalizing out the time profile.
            # But for some reason, there is a frequency dependent scaling discrepency. It works now for 
            # multiple monitors and multiple sources, but we should figure out why this is.
            src = FilteredSource(self.time_src.frequency,self.frequencies,scale,dt,self.time_src) # generate source from broadband response
            amp = 1
        # generate source object
        self.source = mp.EigenModeSource(src,
                    eig_band=self.mode,
                    direction=mp.NO_DIRECTION,
                    eig_kpoint=k0,
                    amplitude=amp,
                    size=self.volume.size,
                    center=self.volume.center,
                    **self.EigenMode_kwargs)
        
        return self.source

    def __call__(self):
        # We just need a workable time profile, so just grab the first available time profile and use that.
        self.time_src = self.sim.sources[0].src

        # Eigenmode data
        ob = self.sim.get_eigenmode_coefficients(self.monitor,[self.mode],**self.EigenMode_kwargs)
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