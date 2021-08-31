import unittest
import numpy as np
import meep as mp

class TestModeDecomposition(unittest.TestCase):

    def test_mode_decomposition(self):
        resolution = 10
        w1 = 1
        w2 = 2
        Lw = 2
        dair = 3.0
        dpml = 5.0
        sy = dpml + dair + w2 + dair + dpml
        half_w1 = 0.5 * w1
        half_w2 = 0.5 * w2
        Si = mp.Medium(epsilon=12.0)
        boundary_layers = [mp.PML(dpml)]
        lcen = 6.67
        fcen = 1 / lcen
        symmetries = [mp.Mirror(mp.Y)]
        Lt = 2
        sx = dpml + Lw + Lt + Lw + dpml
        cell_size = mp.Vector3(sx, sy, 0)
        prism_x = sx + 1
        half_Lt = 0.5 * Lt
        src_pt = mp.Vector3(-0.5 * sx + dpml + 0.2 * Lw, 0, 0)
        sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
                                      center=src_pt,
                                      size=mp.Vector3(0, sy - 2 * dpml, 0),
                                      eig_match_freq=True,
                                      eig_parity=mp.ODD_Z + mp.EVEN_Y)]

        vertices = [mp.Vector3(-prism_x, half_w1),
                    mp.Vector3(prism_x, half_w1),
                    mp.Vector3(prism_x, -half_w1),
                    mp.Vector3(-prism_x, -half_w1)]

        sim = mp.Simulation(resolution=resolution,
                            cell_size=cell_size,
                            boundary_layers=boundary_layers,
                            geometry=[mp.Prism(vertices, height=mp.inf, material=Si)],
                            sources=sources,
                            symmetries=symmetries)

        mon_pt = mp.Vector3(-0.5 * sx + dpml + 0.5 * Lw, 0, 0)
        flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mon_pt, size=mp.Vector3(0, sy - 2 * dpml, 0)))

        sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, src_pt, 1e-9))

        res = sim.get_eigenmode_coefficients(flux, [1], eig_parity=mp.ODD_Z + mp.EVEN_Y)
        incident_coeffs = res.alpha
        incident_flux = mp.get_fluxes(flux)
        incident_flux_data = sim.get_flux_data(flux)

        sim.reset_meep()

        vertices = [mp.Vector3(-prism_x, half_w1),
                    mp.Vector3(-half_Lt, half_w1),
                    mp.Vector3(half_Lt, half_w2),
                    mp.Vector3(prism_x, half_w2),
                    mp.Vector3(prism_x, -half_w2),
                    mp.Vector3(half_Lt, -half_w2),
                    mp.Vector3(-half_Lt, -half_w1),
                    mp.Vector3(-prism_x, -half_w1)]

        sim = mp.Simulation(resolution=resolution,
                            cell_size=cell_size,
                            boundary_layers=boundary_layers,
                            geometry=[mp.Prism(vertices, height=mp.inf, material=Si)],
                            sources=sources,
                            symmetries=symmetries)

        refl_flux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mon_pt, size=mp.Vector3(0, sy - 2 * dpml, 0)))
        sim.load_minus_flux_data(refl_flux, incident_flux_data)

        sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, src_pt, 1e-9))

        res = sim.get_eigenmode_coefficients(refl_flux, [1], eig_parity=mp.ODD_Z + mp.EVEN_Y)
        coeffs = res.alpha
        taper_flux = mp.get_fluxes(refl_flux)

        self.assertAlmostEqual(abs(coeffs[0, 0, 1])**2 / abs(incident_coeffs[0, 0, 0])**2,
                               -taper_flux[0] / incident_flux[0], places=4)

    def test_oblique_waveguide_backward_mode(self):
        sxy = 12.0
        cell_size = mp.Vector3(sxy,sxy,0)

        dpml = 0.6
        pml_layers = [mp.PML(thickness=dpml)]

        fcen = 1/1.55
        rot_angle = np.radians(35.0)
        kpoint = mp.Vector3(1,0,0).rotate(mp.Vector3(0,0,1), rot_angle) * -1.0
        sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen,fwidth=0.1),
                                      center=mp.Vector3(0.5*sxy-3.4,0,0),
                                      size=mp.Vector3(0,sxy,0),
                                      direction=mp.NO_DIRECTION,
                                      eig_kpoint=kpoint,
                                      eig_band=1,
                                      eig_parity=mp.ODD_Z,
                                      eig_match_freq=True)]

        geometry = [mp.Block(center=mp.Vector3(),
                             size=mp.Vector3(mp.inf,1,mp.inf),
                             e1 = mp.Vector3(1,0,0).rotate(mp.Vector3(0,0,1), rot_angle),
                             e2 = mp.Vector3(0,1,0).rotate(mp.Vector3(0,0,1), rot_angle),
                             material=mp.Medium(index=3.5))]

        sim = mp.Simulation(cell_size=cell_size,
                            resolution=20,
                            boundary_layers=pml_layers,
                            sources=sources,
                            geometry=geometry)

        mode = sim.add_mode_monitor(fcen, 0, 1,
                                    mp.FluxRegion(center=mp.Vector3(-0.5*sxy+dpml,0,0),
                                                  size=mp.Vector3(0,sxy,0)),
                                    decimation_factor=1)
        mode_decimated = sim.add_mode_monitor(fcen, 0, 1,
                                              mp.FluxRegion(center=mp.Vector3(-0.5*sxy+dpml,0,0),
                                                            size=mp.Vector3(0,sxy,0)),
                                              decimation_factor=10)

        sim.run(until_after_sources=30)

        flux = mp.get_fluxes(mode)[0]
        coeff = sim.get_eigenmode_coefficients(mode,[1],
                                               direction=mp.NO_DIRECTION,
                                               kpoint_func=lambda f,n: kpoint).alpha[0,0,0]
        flux_decimated = mp.get_fluxes(mode_decimated)[0]
        coeff_decimated = sim.get_eigenmode_coefficients(mode_decimated,[1],
                                                         direction=mp.NO_DIRECTION,
                                                         kpoint_func=lambda f,n: kpoint).alpha[0,0,0]

        print("oblique-waveguide-flux:, {:.6f}, {:.6f}".format(-flux, abs(coeff)**2))
        print("oblique-waveguide-flux (decimated):, {:.6f}, {:.6f}".format(-flux_decimated,
                                                                           abs(coeff_decimated)**2))
        ## the magnitude of |flux| is 100.008731 and so we check two significant digits of accuracy
        self.assertAlmostEqual(-1,abs(coeff)**2/flux,places=2)
        self.assertAlmostEqual(flux,flux_decimated,places=3)
        self.assertAlmostEqual(coeff,coeff_decimated,places=3)


if __name__ == '__main__':
    unittest.main()
