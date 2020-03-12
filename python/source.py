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
        self.peak_time = 0.5 * (self.start_time + self.start_time + 2 * self.width * self.cutoff)

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
    def __init__(self,center_frequency,frequencies,frequency_response,dt,T,time_src):
        dt = dt/2 # divide by two to compensate for staggered E,H time interval
        self.dt = dt
        self.frequencies=frequencies
        self.N = np.round(T/dt)
        f = self.func()

        # calculate dtft of input signal
        signal_t = np.array([time_src.swigobj.current(t,dt) for t in np.arange(0,T,dt)]) # time domain signal
        signal_dtft = np.exp(1j*2*np.pi*frequencies[:,np.newaxis]*np.arange(0,signal_t.size)[np.newaxis,:]*dt)@signal_t # vectorize dtft for speed        

        # multiply sampled dft of input signal with filter transfer function
        H = signal_dtft * frequency_response

        # estimate the impulse response using a sinc function RBN
        self.nodes, self.err = self.estimate_impulse_response(H)

        # initialize super
        super(FilteredSource, self).__init__(src_func=f,center_frequency=center_frequency,is_integrated=False)

    def sinc(self,f,f0):
        omega = 2*np.pi*f
        omega0 = 2*np.pi*f0
        num = np.where(f == f0, self.N, (1-np.exp(1j*(self.N+1)*(omega-omega0)*self.dt)))
        den = np.where(f == f0, 1, (1-np.exp(1j*(omega-omega0)*self.dt)))
        return num/den

    def rect(self,n,f0):
        omega0 = 2*np.pi*f0
        #return np.exp(-1j*omega0*(n-self.N/2)*self.dt)
        return np.where(n < 0 or n > self.N,0,np.exp(-1j*omega0*(n)*self.dt))
    
    def __call__(self,t):
        n = int(np.round((t)/self.dt))
        vec = self.rect(n,self.frequencies)
        return np.dot(vec,self.nodes)
    
    def func(self):
        def _f(t): 
            return self(t)
        return _f
    
    def estimate_impulse_response(self,H):
        # Use vandermonde matrix to calculate weights of each gaussian.
        # Each sinc is centered at each frequency point
        vandermonde = self.sinc(self.frequencies[:,np.newaxis],self.frequencies[np.newaxis,:])
        nodes = np.matmul(np.linalg.pinv(vandermonde),H)
        H_hat = np.matmul(vandermonde,nodes)
        l2_err = np.sum(np.abs(H-H_hat)**2/np.abs(H)**2)
        return nodes, l2_err

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
