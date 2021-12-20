import unittest
import numpy as np
import meep as mp

dB_cm_to_dB_um = 1e-4

class TestConductivity(unittest.TestCase):

    def wvg_flux(self, res, att_coeff):
        """
        Computes the Poynting flux in a single-mode waveguide at two
        locations (5 and 10 μm) downstream from the source given the
        grid resolution res (pixels/μm) and material attenuation
        coefficient att_coeff (dB/cm).
        """

        cell_size = mp.Vector3(14.,14.)

        pml_layers = [mp.PML(thickness=2.)]

        w = 1. # width of waveguide

        fsrc = 0.15 # frequency (in vacuum)

        # note: MPB can only compute modes of lossless material systems.
        #       The imaginary part of ε is ignored and the launched
        #       waveguide mode is therefore inaccurate. For small values
        #       of imag(ε) (which is proportional to att_coeff), this
        #       inaccuracy tends to be insignificant.
        sources = [mp.EigenModeSource(src=mp.GaussianSource(fsrc,fwidth=0.2*fsrc),
                                      center=mp.Vector3(-5.),
                                      size=mp.Vector3(y=10.),
                                      eig_parity=mp.EVEN_Y+mp.ODD_Z)]

        # ref: https://en.wikipedia.org/wiki/Mathematical_descriptions_of_opacity
        # Note that this is the loss of a planewave, which is only approximately
        # the loss of a waveguide mode.  In principle, we could compute the latter
        # semi-analytically if we wanted to run this unit test to greater accuracy
        # (e.g. to test convergence with resolution).
        n_eff = np.sqrt(12.) + 1j * (1/fsrc) * (dB_cm_to_dB_um * att_coeff) / (4 * np.pi)
        eps_eff = n_eff * n_eff
        # ref: https://meep.readthedocs.io/en/latest/Materials/#conductivity-and-complex
        sigma_D = 2 * np.pi * fsrc * np.imag(eps_eff) / np.real(eps_eff)

        geometry = [mp.Block(center=mp.Vector3(),
                             size=mp.Vector3(mp.inf,w,mp.inf),
                             material=mp.Medium(epsilon=np.real(eps_eff),
                                                D_conductivity=sigma_D))]

        sim = mp.Simulation(cell_size=cell_size,
                            resolution=res,
                            boundary_layers=pml_layers,
                            sources=sources,
                            geometry=geometry,
                            symmetries=[mp.Mirror(mp.Y)])

        tran1 = sim.add_flux(fsrc,
                             0,
                             1,
                             mp.FluxRegion(center=mp.Vector3(x=0.), size=mp.Vector3(y=10.)))

        tran2 = sim.add_flux(fsrc,
                             0,
                             1,
                             mp.FluxRegion(center=mp.Vector3(x=5.), size=mp.Vector3(y=10.)))

        sim.run(until_after_sources=20)

        return mp.get_fluxes(tran1)[0], mp.get_fluxes(tran2)[0]


    def test_conductivity(self):
        res = 25. # pixels/μm

        # compute the incident flux for a lossless waveguide
        incident_flux1, incident_flux2 = self.wvg_flux(res, 0.)
        self.assertAlmostEqual(incident_flux1/incident_flux2, 1., places=2)
        print("incident_flux:, {} (measured), 1.0 (expected)".format(incident_flux2/incident_flux1))

        # compute the flux for a lossy waveguide
        att_coeff = 37.46 # dB/cm
        attenuated_flux1, attenuated_flux2 = self.wvg_flux(res, att_coeff)

        L1 = 5.
        expected_att1 = np.exp(-att_coeff * dB_cm_to_dB_um * L1)
        self.assertAlmostEqual(attenuated_flux1/incident_flux2, expected_att1, places=2)
        print("flux:, {}, {:.6f} (measured), {:.6f} (expected)".format(L1,
                                                                       attenuated_flux1/incident_flux2,
                                                                       expected_att1))

        L2 = 10.
        expected_att2 = np.exp(-att_coeff * dB_cm_to_dB_um * L2)
        self.assertAlmostEqual(attenuated_flux2/incident_flux2, expected_att2, places=2)
        print("flux:, {}, {:.6f} (measured), {:.6f} (expected)".format(L2,
                                                                       attenuated_flux2/incident_flux2,
                                                                       expected_att2))


if __name__ == '__main__':
    unittest.main()
