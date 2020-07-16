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
        self.components = EH_TRANSVERSE[self.normal_direction]
    
    def register_monitors(self,frequencies):
        self.frequencies = np.asarray(frequencies)
        self.nf = self.frequencies.size
        #self.monitor = self.sim.add_mode_monitor(frequencies,mp.FluxRegion(center=self.volume.center,size=self.volume.size))
        self.monitor = self.sim.add_dft_fields(self.components,frequencies,where=self.volume,yee_grid=True)
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
        
        iomega = ((1.0 - np.exp(-1j * (2 * np.pi * self.frequencies) * dt)) * (1.0 / dt))
        scale = self.dx * dV * dJ * iomega / np.array([self.time_src.fourier_transform(f) for f in self.frequencies]) # final scale factor
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
        
        sign  =  -1.0 if self.forward else 1.0
        signs = [ -1.0, +1.0, -1.0*sign, +1.0*sign ]
        # -------------------------------------- #
        # Place sources
        # -------------------------------------- #
        self.source = []
        for ci,c in enumerate(self.components):
            for ix,rx in enumerate(np.linspace(self.min_corners[ci].x,self.max_corners[ci].x,self.component_metadata[ci][0])):
                for iy,ry in enumerate(np.linspace(self.min_corners[ci].y,self.max_corners[ci].y,self.component_metadata[ci][1])):
                    for iz,rz in enumerate(np.linspace(self.min_corners[ci].z,self.max_corners[ci].z,self.component_metadata[ci][2])):
                        cur_amp = amp * self.EH[ci][ix,iy,iz] * signs[ci]
                        self.source += [mp.Source(src,component=self.components[3-ci],amplitude=cur_amp,size=mp.Vector3(),center=mp.Vector3(rx,ry,rz))]
        
        return self.source

    def __call__(self):
        # We just need a workable time profile, so just grab the first available time profile and use that.
        self.time_src = self.sim.sources[0].src

        # store the array metadata for each field component
        self.component_metadata = []; self.min_corners = []; self.max_corners = []
        for c in self.components:
            s, min_corner, max_corner = self.sim.get_array_slice_dimensions(c,vol=self.volume)
            self.component_metadata += [s]; self.min_corners += [min_corner]; self.max_corners += [max_corner]
        
        # pull the dft fields on the yee grid
        self.EH = [np.zeros((cm[0],cm[1],cm[2],self.nf),dtype=np.complex128) for cm in self.component_metadata]
        for ci,c in enumerate(self.components):
            for f in range(self.nf):
                self.EH[ci][:,:,:,f] = atleast_3d(self.sim.get_dft_array(self.monitor,c,f))

        # compute and store the eigenmode profiles on the yee grid
        direction = mp.NO_DIRECTION if self.kpoint_func else self.normal_direction
        eh = [np.zeros((cm[0],cm[1],cm[2],self.nf),dtype=np.complex128) for cm in self.component_metadata]
        for f in range(self.nf):
            kpoint = self.kpoint_func(self.frequencies[f],1) if self.kpoint_func else mp.Vector3()
            eig_data = self.sim.get_eigenmode(frequency=self.frequencies[f],direction=direction,
                kpoint=kpoint,where=self.volume,band_num=self.mode,resolution=2*self.sim.resolution,
                **self.EigenMode_kwargs)
            for ci,c in enumerate(self.components):
                for ix,rx in enumerate(np.linspace(self.min_corners[ci].x,self.max_corners[ci].x,self.component_metadata[ci][0])):
                    for iy,ry in enumerate(np.linspace(self.min_corners[ci].y,self.max_corners[ci].y,self.component_metadata[ci][1])):
                        for iz,rz in enumerate(np.linspace(self.min_corners[ci].z,self.max_corners[ci].z,self.component_metadata[ci][2])):
                            eh[ci][ix,iy,iz,f] = eig_data.amplitude(mp.Vector3(rx,ry,rz),c) 
        
        # compute eigenmode coefficient for each freq
        cp = np.sum(np.conj(eh[2]) * self.EH[1],axis=(0,1,2)) - np.sum(np.conj(eh[3]) * self.EH[0],axis=(0,1,2))
        cm = np.sum(np.conj(eh[0]) * self.EH[3],axis=(0,1,2)) - np.sum(np.conj(eh[1]) * self.EH[2],axis=(0,1,2))
        self.eval = -(cp-cm)/4 if self.forward else (cp+cm)/4

        if (((self.volume.size.y == 0) and (self.volume.size.z == 0)) or 
            ((self.volume.size.x == 0) and (self.volume.size.z == 0)) or
            ((self.volume.size.x == 0) and (self.volume.size.y == 0))):
            self.dx = self.sim.resolution ** (-1)
        elif (((self.volume.size.y != 0) and (self.volume.size.z != 0)) or 
            ((self.volume.size.x != 0) and (self.volume.size.z != 0)) or
            ((self.volume.size.x != 0) and (self.volume.size.y != 0))):
            self.dx = self.sim.resolution ** (-2)
        else:
            self.dx = self.sim.resolution ** (-3)

        return self.eval / self.dx
    def get_evaluation(self):
        '''Returns the requested eigenmode coefficient.
        '''
        try:
            return self.eval
        except AttributeError:
            raise RuntimeError("You must first run a forward simulation before resquesting an eigenmode coefficient.")