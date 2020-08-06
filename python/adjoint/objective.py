"""Handling of objective functions and objective quantities."""

from abc import ABC, abstractmethod
import numpy as np
import meep as mp
from .filter_source import FilteredSource

from collections import namedtuple

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

        return self.eval
    def get_evaluation(self):
        '''Returns the requested eigenmode coefficient.
        '''
        try:
            return self.eval
        except AttributeError:
            raise RuntimeError("You must first run a forward simulation before resquesting an eigenmode coefficient.")



class Fourier_Coefficients(ObjectiveQuantitiy):
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

        dt = self.sim.fields.dt

        if self.sim.cell_size.y == 0:
            dV = 1/self.sim.resolution
        elif self.sim.cell_size.z == 0:
            dV = 1/self.sim.resolution * 1/self.sim.resolution
        else:
            dV = 1/self.sim.resolution * 1/self.sim.resolution * 1/self.sim.resolution

        self.sources = []

        if self.num_freq == 1:
            scale = dV * 1j * 2 * np.pi * self.frequencies[0] / self.time_src.fourier_transform(self.frequencies[0])
            amp = -atleast_3d(dJ[0]) * scale
            if self.component in [mp.Hx, mp.Hy, mp.Hz]:
                amp = -amp
            for zi in range(len(self.dg.z)):
                for yi in range(len(self.dg.y)):
                    for xi in range(len(self.dg.x)):
                        if amp[xi, yi, zi] != 0:
                            self.sources += [mp.Source(self.time_src, component=self.component, amplitude= amp[xi, yi, zi],
                            center=mp.Vector3(self.dg.x[xi], self.dg.y[yi], self.dg.z[zi]))]
        else:
            dJ_4d = np.array([atleast_3d(dJ[f]) for f in range(self.num_freq)])
            if self.component in [mp.Hx, mp.Hy, mp.Hz]:
                dJ_4d = -dJ_4d
            for zi in range(len(self.dg.z)):
                for yi in range(len(self.dg.y)):
                    for xi in range(len(self.dg.x)):
                        scale = -dJ_4d[:,xi,yi,zi] * dV * 1j * 2 * np.pi * self.frequencies / np.array([self.time_src.fourier_transform(f) for f in self.frequencies])
                        src = FilteredSource(self.time_src.frequency,self.frequencies,scale,dt,self.time_src)
                        self.sources += [mp.Source(src, component=self.component, amplitude= 1,
                                center=mp.Vector3(self.dg.x[xi], self.dg.y[yi], self.dg.z[zi]))]

        return self.sources

    def __call__(self):
        self.time_src = self.sim.sources[0].src

        self.dg = Grid(*self.sim.get_array_metadata(dft_cell=self.monitor))
        self.eval = np.array([self.sim.get_dft_array(self.monitor, self.component, i) for i in range(self.num_freq)]) #Shape = (num_freq, [pts])
        return self.eval

    def get_evaluation(self):

        try:
            return self.eval
        except AttributeError:
            raise RuntimeError("You must first run a forward simulation.")


Grid = namedtuple('Grid', ['x', 'y', 'z', 'w'])

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





class Far_Coefficients(ObjectiveQuantitiy):
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

        dt = self.sim.fields.dt

        if self.sim.cell_size.y == 0:
            dV = 1/self.sim.resolution
        elif self.sim.cell_size.z == 0:
            dV = 1/self.sim.resolution * 1/self.sim.resolution
        else:
            dV = 1/self.sim.resolution * 1/self.sim.resolution * 1/self.sim.resolution

        self.sources = []


        freq_scale = dV * 1j * 2 * np.pi * self.frequencies / np.array([self.time_src.fourier_transform(f) for f in self.frequencies])

        #TODO far_pts in 3d or cylindrical, perhaps py_v3_to_vec from simulation.py
        self.all_near_data = self.monitor.swigobj.near_fds(mp.vec(self.far_pt.x, self.far_pt.y))


        for near_data in self.all_near_data:
            near_pt = near_data.near_x
            near_pt = mp.Vector3(near_pt.x(), near_pt.y(), near_pt.z())
            cur_comp = near_data.near_fd_comp
            cur_green_matrix = np.array(near_data.matrix_elt)

            if cur_comp in [mp.Ex, mp.Ey, mp.Ez]:
                scale = -freq_scale * np.sum(dJ*cur_green_matrix, axis=1)
            else:
                scale = freq_scale * np.sum(dJ*cur_green_matrix, axis=1)

            if self.num_freq == 1:
                self.sources += [mp.Source(self.time_src, component=cur_comp, amplitude=scale[0], center=near_pt)]
            else:
                src = FilteredSource(self.time_src.frequency,self.frequencies,scale,dt,self.time_src)
                self.sources += [mp.Source(src, component=cur_comp, amplitude=1, center=near_pt)]

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
