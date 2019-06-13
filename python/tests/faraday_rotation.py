from __future__ import division

import unittest
import numpy as np
import meep as mp

## Farady rotation rate for gyrotropic Lorentzian medium
def kgyro_lorentzian(freq, epsn, f0, gamma, sigma, b0):
    dfsq = (f0**2 - 1j*freq*gamma - freq**2)
    eperp = epsn + sigma * f0**2 * dfsq / (dfsq**2 - (freq*b0)**2)
    eta = sigma * f0**2 * freq * b0 / (dfsq**2 - (freq*b0)**2)
    return 2*np.pi*freq * np.sqrt(0.5*(eperp - np.sqrt(eperp**2 - eta**2)))

## Farady rotation rate for gyrotropic Drude medium
def kgyro_drude(freq, epsn, f0, gamma, sigma, b0):
    dfsq = - 1j*freq*gamma - freq**2
    eperp = epsn + sigma * f0**2 * dfsq / (dfsq**2 - (freq*b0)**2)
    eta = sigma * f0**2 * freq * b0 / (dfsq**2 - (freq*b0)**2)
    return 2*np.pi*freq * np.sqrt(0.5*(eperp - np.sqrt(eperp**2 - eta**2)))

## Farady rotation rate for Landau-Lifshitz-Gilbert medium
def kgyro_llg(freq, epsn, f0, gamma, sigma, alpha):
    df1 = f0 - 1j*freq*alpha
    df2 = freq + 1j*gamma
    eperp = epsn + sigma * df1/(df1**2 - df2**2)
    eta = sigma * df2 / (df1**2 - df2**2)
    return 2*np.pi*freq * np.sqrt(0.5*(eperp - np.sqrt(eperp**2 - eta**2)))

class TestFaradayRotation(unittest.TestCase):
    ## Simulate a linearly polarized plane wave traveling along the gyrotropy axis.
    ## Extract Faraday rotation angle by comparing the Ex and Ey amplitudes, and
    ## compare to the theoretical result predicted by rotation rate KPRED
    ## up to relative tolerance RTOL.
    def check_rotation(self, mat, L, fsrc, zsrc, resolution, tmax, zout, kpred, rtol):
        cell = mp.Vector3(0, 0, L)
        pml_layers = [mp.PML(thickness=1.0, direction=mp.Z)]
        sources = [mp.Source(mp.ContinuousSource(frequency=fsrc),
                             component=mp.Ex, center=mp.Vector3(0, 0, zsrc))]

        self.sim = mp.Simulation(cell_size=cell, geometry=[], sources=sources,
                                 boundary_layers=pml_layers,
                                 default_material=mat, resolution=resolution)

        record_vol = mp.Volume(center=mp.Vector3(0, 0, zout))
        record_Ex, record_Ey, record_t = [], [], []

        def record_ex_ey(sim):
            record_Ex.append(sim.get_array(vol=record_vol, component=mp.Ex))
            record_Ey.append(sim.get_array(vol=record_vol, component=mp.Ey))
            record_t.append(sim.meep_time())

        self.sim.run(mp.after_time(0.5*tmax, mp.at_every(1e-6, record_ex_ey)), until=tmax)

        ex_rel = np.amax(abs(np.fft.fft(record_Ex)))
        ey_rel = np.amax(abs(np.fft.fft(record_Ey)))
        result = np.arctan2(ey_rel, ex_rel)

        Ex_theory = np.abs(np.cos(kpred * (zout - zsrc)).real)
        Ey_theory = np.abs(np.sin(kpred * (zout - zsrc)).real)
        expected = np.arctan2(Ey_theory, Ex_theory)

        np.testing.assert_allclose(expected, result, rtol=rtol)

    def test_faraday_rotation(self):
        L, zsrc, zout = 12.0, -4.5, 4.0
        freq, tmax = 0.8, 300.0
        resolution = 24

        ## Test gyrotropic lorentzian (2% tolerance)
        epsn, f0, gamma, sn, b0  = 1.5, 1.0, 1e-3, 0.1, 0.15
        susc = [mp.GyrotropicLorentzianSusceptibility(frequency=f0, gamma=gamma, sigma=sn,
                                                      bias=mp.Vector3(0, 0, b0))]
        mat = mp.Medium(epsilon=epsn, mu=1, E_susceptibilities=susc)
        k = kgyro_lorentzian(freq, epsn, f0, gamma, sn, b0)
        self.check_rotation(mat, L, freq, zsrc, resolution, tmax, zout, k, 0.02)

        ## Test gyrotropic Drude medium (2% tolerance)
        susc = [mp.GyrotropicDrudeSusceptibility(frequency=f0, gamma=gamma, sigma=sn,
                                                 bias=mp.Vector3(0, 0, b0))]
        mat = mp.Medium(epsilon=epsn, mu=1, E_susceptibilities=susc)
        k = kgyro_drude(freq, epsn, f0, gamma, sn, b0)
        self.check_rotation(mat, L, freq, zsrc, resolution, tmax, zout, k, 0.02)

        ## Test Landau-Lifshitz-Gilbert medium (5% tolerance)
        alpha = 1e-5
        susc = [mp.GyrotropicSaturatedSusceptibility(frequency=f0, gamma=gamma, sigma=sn,
                                                     alpha=alpha,
                                                     bias=mp.Vector3(0, 0, 1.0))]
        mat = mp.Medium(epsilon=epsn, mu=1, E_susceptibilities=susc)
        k = kgyro_llg(freq, epsn, f0, gamma, sn, alpha)
        self.check_rotation(mat, L, freq, zsrc, resolution, tmax, zout, k, 0.05)

if __name__ == '__main__':
    unittest.main()
