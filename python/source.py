from __future__ import division

import meep as mp
from meep.geom import Vector3, check_nonnegative


def check_positive(prop, val):
    if val > 0:
        return val
    else:
        raise ValueError("{} must be positive. Got {}".format(prop, val))


class Source(object):

    def __init__(self, src, component, center, size=Vector3(), amplitude=1.0, amp_func=None):
        self.src = src
        self.component = component
        self.center = center
        self.size = size
        self.amplitude = amplitude
        self.amp_func = amp_func


class SourceTime(object):

    def __init__(self, is_integrated=False):
        self.is_integrated = is_integrated


class ContinuousSource(SourceTime):

    def __init__(self, frequency, start_time=0, end_time=float('inf'), width=0, cutoff=3.0):
        super(ContinuousSource, self).__init__()
        self.frequency = frequency
        self.start_time = start_time
        self.end_time = end_time
        self.width = width
        self.cutoff = cutoff
        self.swigobj = mp.continuous_src_time(frequency, width, start_time, end_time, cutoff)
        self.swigobj.is_integrated = self.is_integrated


class GaussianSource(SourceTime):

    def __init__(self, frequency, width=0, fwidth=float('inf'), start_time=0, cutoff=5.0):
        super(GaussianSource, self).__init__()
        self.frequency = frequency
        self.width = max(width, 1 / fwidth)
        self.start_time = start_time
        self.cutoff = cutoff
        self.swigobj = mp.gaussian_src_time(frequency, self.width, start_time, start_time + 2 * self.width * cutoff)
        self.swigobj.is_integrated = self.is_integrated


class CustomSource(SourceTime):

    def __init__(self, src_func, start_time=float('-inf'), end_time=float('inf')):
        super(CustomSource, self).__init__()
        self.src_func = src_func
        self.start_time = start_time
        self.end_time = end_time
        self.swigobj = mp.custom_src_time(src_func, start_time, end_time)
        self.swigobj.is_integrated = self.is_integrated


class EigenModeSource(Source):

    def __init__(self, src, center,
                 eig_lattice_size=None,
                 eig_lattice_center=None,
                 component=mp.Dielectric,
                 direction=-1,
                 eig_band=1,
                 eig_kpoint=Vector3(),
                 eig_match_freq=True,
                 eig_parity=0,
                 eig_resolution=0,
                 eig_tolerence=1e-7,
                 **kwargs):

        super(EigenModeSource, self).__init__(src, component, center, **kwargs)
        self.eig_lattice_size = eig_lattice_size
        self.eig_lattice_center = eig_lattice_center
        self.component = component
        self.direction = direction
        self.eig_band = eig_band
        self.eig_kpoint = eig_kpoint
        self.eig_match_freq = eig_match_freq
        self.eig_parity = eig_parity
        self.eig_resolution = eig_resolution
        self.eig_tolerence = eig_tolerence

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
    def eig_tolerence(self):
        return self._eig_tolerence

    @eig_tolerence.setter
    def eig_tolerence(self, val):
        self._eig_tolerence = check_positive('EigenModeSource.eig_tolerence', val)
