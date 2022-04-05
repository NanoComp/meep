import unittest
import numpy as np
import meep as mp

class TestModeDecomposition(unittest.TestCase):

    def test_linear_taper_2d(self):
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
        cell_size = mp.Vector3(sx,sy,0)
        prism_x = sx + 1
        half_Lt = 0.5 * Lt
        src_pt = mp.Vector3(-0.5*sx+dpml+0.2*Lw,0,0)
        sources = [mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.2*fcen),
                                      center=src_pt,
                                      size=mp.Vector3(0,sy-2*dpml,0),
                                      eig_match_freq=True,
                                      eig_parity=mp.ODD_Z+mp.EVEN_Y)]

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

        mon_pt = mp.Vector3(-0.5*sx+dpml+0.5*Lw,0,0)
        flux = sim.add_flux(fcen,
                            0,
                            1,
                            mp.FluxRegion(center=mon_pt,
                                          size=mp.Vector3(0,sy-2*dpml,0)))

        sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, src_pt, 1e-9))

        res = sim.get_eigenmode_coefficients(flux,
                                             [1],
                                             eig_parity=mp.ODD_Z+mp.EVEN_Y)
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

        refl_flux = sim.add_flux(fcen,
                                 0,
                                 1,
                                 mp.FluxRegion(center=mon_pt,
                                               size=mp.Vector3(0,sy-2*dpml,0)))
        sim.load_minus_flux_data(refl_flux, incident_flux_data)

        sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, src_pt, 1e-9))

        res = sim.get_eigenmode_coefficients(refl_flux,
                                             [1],
                                             eig_parity=mp.ODD_Z+mp.EVEN_Y)
        coeffs = res.alpha
        taper_flux = mp.get_fluxes(refl_flux)

        self.assertAlmostEqual(abs(coeffs[0,0,1])**2 / abs(incident_coeffs[0,0,0])**2,
                               -taper_flux[0] / incident_flux[0],
                               places=4)

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
        mode_decimated = sim.add_mode_monitor(fcen,
                                              0,
                                              1,
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


    def test_grating_3d(self):
        """Unit test for mode decomposition in 3d.

        Verifies that the reflectance and transmittance in the z
        direction at a single wavelength for a unit cell of a
        3d grating using a normally incident planewave is equivalent
        to the sum of the Poynting flux (normalized by the flux
        of the input source) for all the individual reflected
        and transmitted diffracted orders.
        """
        resolution = 25  # pixels/μm

        nSi = 3.45
        Si = mp.Medium(index=nSi)
        nSiO2 = 1.45
        SiO2 = mp.Medium(index=nSiO2)

        wvl = 0.5  # wavelength
        fcen = 1/wvl

        dpml = 2.0  # PML thickness
        dsub = 3.0  # substrate thickness
        dair = 3.0  # air padding
        hcyl = 0.5  # cylinder height
        rcyl = 0.2  # cylinder radius

        sx = 1.1
        sy = 0.8
        sz = dpml+dsub+hcyl+dair+dpml

        cell_size = mp.Vector3(sx,sy,sz)

        boundary_layers = [mp.PML(thickness=dpml,direction=mp.Z)]

        # periodic boundary conditions
        k_point = mp.Vector3()

        src_cmpt = mp.Ex
        sources = [mp.Source(src=mp.GaussianSource(fcen,fwidth=0.2*fcen),
                             size=mp.Vector3(sx,sy,0),
                             center=mp.Vector3(0,0,-0.5*sz+dpml),
                             component=src_cmpt)]

        sim = mp.Simulation(resolution=resolution,
                            cell_size=cell_size,
                            sources=sources,
                            default_material=SiO2,
                            boundary_layers=boundary_layers,
                            k_point=k_point)

        flux = sim.add_mode_monitor(fcen,
                                    0,
                                    1,
                                    mp.ModeRegion(center=mp.Vector3(0,0,-0.5*sz+dpml+0.5*dsub),
                                                  size=mp.Vector3(sx,sy,0)))

        stop_cond = mp.stop_when_fields_decayed(20,
                                                src_cmpt,
                                                mp.Vector3(0,0,0.5*sz-dpml-0.5*dair),
                                                1e-6)
        sim.run(until_after_sources=stop_cond)

        input_flux = mp.get_fluxes(flux)
        input_flux_data = sim.get_flux_data(flux)

        sim.reset_meep()

        geometry = [mp.Block(size=mp.Vector3(mp.inf,mp.inf,dpml+dsub),
                             center=mp.Vector3(0,0,-0.5*sz+0.5*(dpml+dsub)),
                             material=SiO2),
                    mp.Cylinder(height=hcyl,
                                radius=rcyl,
                                center=mp.Vector3(0,0,-0.5*sz+dpml+dsub+0.5*hcyl),
                                material=Si)]

        sim = mp.Simulation(resolution=resolution,
                            cell_size=cell_size,
                            sources=sources,
                            geometry=geometry,
                            boundary_layers=boundary_layers,
                            k_point=k_point)

        refl_flux = sim.add_mode_monitor(fcen,
                                         0,
                                         1,
                                         mp.ModeRegion(center=mp.Vector3(0,0,-0.5*sz+dpml+0.5*dsub),
                                                       size=mp.Vector3(sx,sy,0)))
        sim.load_minus_flux_data(refl_flux,input_flux_data)

        tran_flux = sim.add_mode_monitor(fcen,
                                         0,
                                         1,
                                         mp.ModeRegion(center=mp.Vector3(0,0,0.5*sz-dpml),
                                                       size=mp.Vector3(sx,sy,0)))

        sim.run(until_after_sources=stop_cond)

        # sum the Poynting flux in z direction for all reflected orders
        Rsum = 0

        # number of reflected modes/orders in SiO2 in x and y directions (upper bound)
        nm_x = int(fcen*nSiO2*sx) + 1
        nm_y = int(fcen*nSiO2*sy) + 1
        for m_x in range(nm_x):
            for m_y in range(nm_y):
                for S_pol in [False,True]:
                    res = sim.get_eigenmode_coefficients(refl_flux,
                                                         mp.DiffractedPlanewave([m_x,m_y,0],
                                                                                mp.Vector3(1,0,0),
                                                                                1 if S_pol else 0,
                                                                                0 if S_pol else 1))
                    r_coeffs = res.alpha
                    Rmode = abs(r_coeffs[0,0,1])**2/input_flux[0]
                    print("refl-order:, {}, {}, {}, {:.6f}".format("s" if S_pol else "p",m_x,m_y,Rmode))
                    if m_x == 0 and m_y == 0:
                        Rsum += Rmode
                    elif (m_x != 0 and m_y == 0) or (m_x == 0 and m_y != 0):
                        Rsum += 2*Rmode
                    else:
                        Rsum += 4*Rmode


        # sum the Poynting flux in z direction for all transmitted orders
        Tsum = 0

        # number of transmitted modes/orders in air in x and y directions (upper bound)
        nm_x = int(fcen*sx) + 1
        nm_y = int(fcen*sy) + 1
        for m_x in range(nm_x):
            for m_y in range(nm_y):
                for S_pol in [False,True]:
                    res = sim.get_eigenmode_coefficients(tran_flux,
                                                         mp.DiffractedPlanewave([m_x,m_y,0],
                                                                                mp.Vector3(1,0,0),
                                                                                1 if S_pol else 0,
                                                                                0 if S_pol else 1))
                    t_coeffs = res.alpha
                    Tmode = abs(t_coeffs[0,0,0])**2/input_flux[0]
                    print("tran-order:, {}, {}, {}, {:.6f}".format("s" if S_pol else "p",m_x,m_y,Tmode))
                    if m_x == 0 and m_y == 0:
                        Tsum += Tmode
                    elif (m_x != 0 and m_y == 0) or (m_x == 0 and m_y != 0):
                        Tsum += 2*Tmode
                    else:
                        Tsum += 4*Tmode


        r_flux = mp.get_fluxes(refl_flux)
        t_flux = mp.get_fluxes(tran_flux)
        Rflux = -r_flux[0]/input_flux[0]
        Tflux =  t_flux[0]/input_flux[0]

        print("refl:, {}, {}".format(Rsum,Rflux))
        print("tran:, {}, {}".format(Tsum,Tflux))
        print("sum:,  {}, {}".format(Rsum+Tsum,Rflux+Tflux))

        self.assertAlmostEqual(Rsum,Rflux,places=2)
        self.assertAlmostEqual(Tsum,Tflux,places=2)
        self.assertAlmostEqual(Rsum+Tsum,1.00,places=2)


if __name__ == '__main__':
    unittest.main()
