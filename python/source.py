from __future__ import division

import meep as mp
from meep.geom import Vector3, check_nonnegative
from scipy import signal
import numpy as np
from scipy.special import erf

def check_positive(prop, val):
    if val > 0:
        return val
    else:
        raise ValueError("{} must be positive. Got {}".format(prop, val))

class Source(object):

    def __init__(self, src, component, center, size=Vector3(), amplitude=1.0, amp_func=None,
                 amp_func_file='', amp_data=None):
        self.src = src
        self.component = component
        self.center = Vector3(*center)
        self.size = Vector3(*size)
        self.amplitude = complex(amplitude)
        self.amp_func = amp_func
        self.amp_func_file = amp_func_file
        self.amp_data = amp_data


class SourceTime(object):

    def __init__(self, is_integrated=False):
        self.is_integrated = is_integrated


class ContinuousSource(SourceTime):

    def __init__(self, frequency=None, start_time=0, end_time=1.0e20, width=0,
                 fwidth=float('inf'), cutoff=3.0, wavelength=None, **kwargs):

        if frequency is None and wavelength is None:
            raise ValueError("Must set either frequency or wavelength in {}.".format(self.__class__.__name__))

        super(ContinuousSource, self).__init__(**kwargs)
        self.frequency = 1 / wavelength if wavelength else float(frequency)
        self.start_time = start_time
        self.end_time = end_time
        self.width = max(width, 1 / fwidth)
        self.cutoff = cutoff
        self.swigobj = mp.continuous_src_time(self.frequency, self.width, self.start_time,
                                              self.end_time, self.cutoff)
        self.swigobj.is_integrated = self.is_integrated


class GaussianSource(SourceTime):

    def __init__(self, frequency=None, width=0, fwidth=float('inf'), start_time=0, cutoff=5.0, wavelength=None,
                 **kwargs):
        if frequency is None and wavelength is None:
            raise ValueError("Must set either frequency or wavelength in {}.".format(self.__class__.__name__))

        super(GaussianSource, self).__init__(**kwargs)
        self.frequency = 1 / wavelength if wavelength else float(frequency)
        self.width = max(width, 1 / fwidth)
        self.start_time = start_time
        self.cutoff = cutoff

        self.swigobj = mp.gaussian_src_time(self.frequency, self.width, self.start_time,
                                            self.start_time + 2 * self.width * self.cutoff)
        self.swigobj.is_integrated = self.is_integrated

    def fourier_transform(self, freq):
        return self.swigobj.fourier_transform(freq)

class CustomSource(SourceTime):

    def __init__(self, src_func, start_time=-1.0e20, end_time=1.0e20, center_frequency=0, **kwargs):
        super(CustomSource, self).__init__(**kwargs)
        self.src_func = src_func
        self.start_time = start_time
        self.end_time = end_time
        self.center_frequency = center_frequency
        self.swigobj = mp.custom_src_time(src_func, start_time, end_time, center_frequency)
        self.swigobj.is_integrated = self.is_integrated

class FilteredSource(CustomSource):
    def __init__(self,center_frequency,frequencies,frequency_response,dt,T,time_src,min_err=1e-6):
        dt = dt/2
        self.dt = dt
        self.center_frequency=center_frequency
        self.frequencies=frequencies
        self.time_src=time_src
        self.min_err = min_err
        f = self.func()

        # calculate dtft of input signal
        signal_t = np.array([time_src.swigobj.current(t,dt) for t in np.arange(0,T,dt)]) # time domain signal
        signal_dtft = np.exp(1j*2*np.pi*frequencies[:,np.newaxis]*np.arange(0,signal_t.size)[np.newaxis,:]*dt)@signal_t # vectorize dtft for speed

        # multiply sampled dft of input signal with filter transfer function
        H = signal_dtft#* frequency_response

        self.estimate_impulse_response(H)

        # initialize super
        super(FilteredSource, self).__init__(src_func=f,center_frequency=self.center_frequency,is_integrated=False)
        
    def gaussian(self,f,f0,fwidth):
        return np.exp(-0.5*((f-f0)/fwidth)**2)
    def antiderivative(self,f,n,f0,fwidth,T):
        a = np.sqrt(np.pi/2)*fwidth
        phase = np.exp(-2*np.pi*n*T*(np.pi*n*T*fwidth*fwidth + 1j*f0))
        kernel = erf(f/(np.sqrt(2)*fwidth)-f0/(np.sqrt(2)*fwidth)+1j*np.sqrt(2)*np.pi*n*fwidth*T)
        return a*phase*kernel
    def dtft_gaussian(self,n,f0,fwidth,T):
        f_start = 0
        f_end = 1/T
        return T*(self.antiderivative(f_end,n,f0,fwidth,T) - self.antiderivative(f_start,n,f0,fwidth,T))
    def __call__(self,t):
        n = int(np.round(t/self.dt))
        #print(t/self.dt,n)
        vec = self.dtft_gaussian(n,self.frequencies,self.gaus_widths,self.dt)
        # simple RBF with gaussian kernel reduces to inner product at time step
        return np.dot(vec,self.nodes)
    
    def func(self):
        def _f(t): 
            return self(t)
        return _f
    
    def estimate_impulse_response(self,H):
        '''
        find gaussian weighting coefficients.

        TODO use optimizer to find optimal gaussian widths
        '''
        # Use vandermonde matrix to calculate weights of each gaussian.
        # Each gaussian is centered at each frequency point
        def rbf_l2(fwidth):
            vandermonde = self.gaussian(self.frequencies[:,np.newaxis],self.frequencies[np.newaxis,:],fwidth)
            nodes = np.matmul(np.linalg.pinv(vandermonde),H)
            H_hat = np.matmul(vandermonde,nodes)
            l2_err = np.sum(np.abs(H-H_hat)**2)
            return nodes, l2_err
        
        df = self.frequencies[2] - self.frequencies[1]
        err_high = True
        fwidth = 1/self.time_src.width
        
        # Iterate through smaller and smaller widths until error is small enough or width is distance between frequency points
        while err_high:
            nodes, l2_err = rbf_l2(fwidth)
            if l2_err < self.min_err or fwidth < df:
                err_high = False
            else:
                fwidth = 0.5 * fwidth
            print(l2_err)
        self.gaus_widths = fwidth
        self.nodes = nodes

        from matplotlib import pyplot as plt

        temp = self.gaussian(self.frequencies[:,np.newaxis],self.frequencies,fwidth)
        i_hat = np.inner(self.nodes,temp)

class EigenModeSource(Source):

    def __init__(self,
                 src,
                 center,
                 eig_lattice_size=None,
                 eig_lattice_center=None,
                 component=mp.ALL_COMPONENTS,
                 direction=mp.AUTOMATIC,
                 eig_band=1,
                 eig_kpoint=Vector3(),
                 eig_match_freq=True,
                 eig_parity=mp.NO_PARITY,
                 eig_resolution=0,
                 eig_tolerance=1e-12,
                 **kwargs):

        super(EigenModeSource, self).__init__(src, component, center, **kwargs)
        self.eig_lattice_size = eig_lattice_size
        self.eig_lattice_center = eig_lattice_center
        self.component = component
        self.direction = direction
        self.eig_band = eig_band
        self.eig_kpoint = mp.Vector3(*eig_kpoint)
        self.eig_match_freq = eig_match_freq
        self.eig_parity = eig_parity
        self.eig_resolution = eig_resolution
        self.eig_tolerance = eig_tolerance

    @property
    def eig_lattice_size(self):
        return self._eig_lattice_size

    @eig_lattice_size.setter
    def eig_lattice_size(self, val):
        if val is None:
            self._eig_lattice_size = self.size
        else:
            self._eig_lattice_size = val

    @property
    def eig_lattice_center(self):
        return self._eig_lattice_center

    @eig_lattice_center.setter
    def eig_lattice_center(self, val):
        if val is None:
            self._eig_lattice_center = self.center
        else:
            self._eig_lattice_center = val

    @property
    def eig_band(self):
        return self._eig_band

    @eig_band.setter
    def eig_band(self, val):
        self._eig_band = check_positive('EigenModeSource.eig_band', val)

    @property
    def eig_resolution(self):
        return self._eig_resolution

    @eig_resolution.setter
    def eig_resolution(self, val):
        self._eig_resolution = check_nonnegative('EigenModeSource.eig_resolution', val)

    @property
    def eig_tolerance(self):
        return self._eig_tolerance

    @eig_tolerance.setter
    def eig_tolerance(self, val):
        self._eig_tolerance = check_positive('EigenModeSource.eig_tolerance', val)

    def eig_power(self,freq):
        amp = self.amplitude
        if callable(getattr(self.src, "fourier_transform", None)):
           amp *= self.src.fourier_transform(freq)
        return abs(amp)**2
