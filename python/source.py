from __future__ import division

import meep as mp
from meep.geom import Vector3, check_nonnegative
from scipy import signal
import numpy as np
from scipy.interpolate import interp1d

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
    def __init__(self,center_frequency,frequencies,frequency_response,num_taps,dt,time_src):
        self.center_frequency=center_frequency
        self.frequencies=frequencies
        self.frequency_response=frequency_response
        self.num_taps=num_taps
        self.dt=dt/2 # the "effective" dt needs double resolution for staggered yee grid
        self.time_src=time_src
        self.current_time = None

        f = self.func()

        # initialize super
        super(FilteredSource, self).__init__(src_func=f,center_frequency=self.center_frequency)

        # calculate equivalent sample rate
        self.fs = 1/self.dt

        # estimate impulse response from frequency response
        self.estimate_impulse_response()

    def filter(self,t):
        idx = int(np.round(t/self.dt))
        if idx >= self.num_taps:
            return 0
        else:
            return self.taps[idx]
        
        
        '''print(idx)
        if self.current_time is None or self.current_time != t:
            self.current_time = t # increase current time (we went through all the current components)
            np.roll(self.memory, -1) # shift feedforward memory
            self.memory[0] = self.time_src.swigobj.dipole(t) # update current memory slot
            self.current_y = np.dot(self.memory,self.taps) # calculate filter response as inner product of taps and memory
        return self.current_y'''
    
    def estimate_impulse_response(self):
        # calculate band edges from target frequencies
        w = self.frequencies/(self.fs/2) * np.pi
        D = self.frequency_response
        self.taps = self.spline_fit(self.num_taps,w,D)

        # allocate filter memory taps
        self.memory = np.zeros(self.taps.shape,dtype=np.complex128)
    
    def func(self):
        def _f(t): 
            return self.filter(t)
        return _f
    
    def spline_fit(self,num_taps,freqs,h_desired):
        num_taps = 2000
        # fit real part
        real_x = np.concatenate(([0],freqs,[np.pi]))
        real_y = np.concatenate(([0],np.real(h_desired),[0]))
        fr = interp1d(real_x, real_y, kind='cubic')

        # fit imaginary part
        imag_x = np.concatenate(([0],freqs,[np.pi]))
        imag_y = np.concatenate(([0],np.imag(h_desired),[0]))
        fi = interp1d(imag_x, imag_y, kind='cubic')

        # formulate hermitian filter response with specified number of taps
        freqs_filter = numpy.fft.fftfreq#np.linspace(-np.pi,np.pi,num_taps)
        zero_freq = np.argmin(np.abs(freqs_filter))+1
        print(freqs_filter[zero_freq])

        filter_pos = fr(freqs_filter[zero_freq:]) + 1j*fi(freqs_filter[zero_freq:])
        filter_neg = np.flipud(np.conjugate(filter_pos))

        if num_taps %2 == 0: #even
            filter_both = np.concatenate((filter_pos,filter_neg))
        else:
            filter_both = np.concatenate((filter_pos,filter_neg[:-1]))

        from matplotlib import pyplot as plt
        plt.figure()
        #plt.plot(freqs_filter,np.abs(filter_both))
        plt.plot(freqs_filter[zero_freq:],fr(freqs_filter[zero_freq:]))
        plt.plot(freqs,np.real(h_desired),'o')
        plt.show()

        print(num_taps)
        print(filter_both.size)
        quit()

        # ifft to get sampled impulse response

        return
    def lstsqrs(self,num_taps,freqs,h_desired):
        n_freqs = freqs.size
        vandermonde_left = np.zeros((n_freqs,num_taps),dtype=np.complex128)
        vandermonde_right = np.zeros((n_freqs,num_taps),dtype=np.complex128)
        for iom, om in enumerate(freqs):
            for it in range(num_taps):
                vandermonde_left[iom,it] = np.exp(-1j*it*om)
                vandermonde_right[iom,it] = np.exp(1j*it*om)
        vandermonde = np.vstack((vandermonde_left,vandermonde_right))
        h_desired_full = np.hstack((h_desired,np.conj(h_desired)))
        
        a = np.matmul(np.linalg.pinv(vandermonde), h_desired_full)
        _, h_hat = signal.freqz(a,worN=freqs)

        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(freqs,np.abs(h_desired))
        plt.plot(freqs,np.abs(h_hat),'--')
        plt.show()
        quit()
        
        self.l2_error = np.sqrt(np.sum(np.abs(h_hat - h_desired)**2))
        print(self.l2_error)

        # account for dtft scaling
        a = a*self.dt*np.sqrt(2)
        return a

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
