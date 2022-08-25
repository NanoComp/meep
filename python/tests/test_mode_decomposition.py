import cmath
import math
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
        cell_size = mp.Vector3(sx, sy, 0)
        prism_x = sx + 1
        half_Lt = 0.5 * Lt
        src_pt = mp.Vector3(-0.5 * sx + dpml + 0.2 * Lw, 0, 0)
        sources = [
            mp.EigenModeSource(
                src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
                center=src_pt,
                size=mp.Vector3(0, sy - 2 * dpml, 0),
                eig_match_freq=True,
                eig_parity=mp.ODD_Z + mp.EVEN_Y,
            )
        ]

        vertices = [
            mp.Vector3(-prism_x, half_w1),
            mp.Vector3(prism_x, half_w1),
            mp.Vector3(prism_x, -half_w1),
            mp.Vector3(-prism_x, -half_w1),
        ]

        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cell_size,
            boundary_layers=boundary_layers,
            geometry=[mp.Prism(vertices, height=mp.inf, material=Si)],
            sources=sources,
            symmetries=symmetries,
        )

        mon_pt = mp.Vector3(-0.5 * sx + dpml + 0.5 * Lw, 0, 0)
        flux = sim.add_flux(
            fcen,
            0,
            1,
            mp.FluxRegion(center=mon_pt, size=mp.Vector3(0, sy - 2 * dpml, 0)),
        )

        sim.run(
            until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, src_pt, 1e-9)
        )

        res = sim.get_eigenmode_coefficients(flux, [1], eig_parity=mp.ODD_Z + mp.EVEN_Y)
        incident_coeffs = res.alpha
        incident_flux = mp.get_fluxes(flux)
        incident_flux_data = sim.get_flux_data(flux)

        sim.reset_meep()

        vertices = [
            mp.Vector3(-prism_x, half_w1),
            mp.Vector3(-half_Lt, half_w1),
            mp.Vector3(half_Lt, half_w2),
            mp.Vector3(prism_x, half_w2),
            mp.Vector3(prism_x, -half_w2),
            mp.Vector3(half_Lt, -half_w2),
            mp.Vector3(-half_Lt, -half_w1),
            mp.Vector3(-prism_x, -half_w1),
        ]

        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cell_size,
            boundary_layers=boundary_layers,
            geometry=[mp.Prism(vertices, height=mp.inf, material=Si)],
            sources=sources,
            symmetries=symmetries,
        )

        refl_flux = sim.add_flux(
            fcen,
            0,
            1,
            mp.FluxRegion(center=mon_pt, size=mp.Vector3(0, sy - 2 * dpml, 0)),
        )
        sim.load_minus_flux_data(refl_flux, incident_flux_data)

        sim.run(
            until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, src_pt, 1e-9)
        )

        res = sim.get_eigenmode_coefficients(
            refl_flux, [1], eig_parity=mp.ODD_Z + mp.EVEN_Y
        )
        coeffs = res.alpha
        taper_flux = mp.get_fluxes(refl_flux)

        self.assertAlmostEqual(
            abs(coeffs[0, 0, 1]) ** 2 / abs(incident_coeffs[0, 0, 0]) ** 2,
            -taper_flux[0] / incident_flux[0],
            places=4,
        )

    def test_oblique_waveguide_backward_mode(self):
        sxy = 12.0
        cell_size = mp.Vector3(sxy, sxy, 0)

        dpml = 0.6
        pml_layers = [mp.PML(thickness=dpml)]

        fcen = 1 / 1.55
        rot_angle = np.radians(35.0)
        kpoint = mp.Vector3(1, 0, 0).rotate(mp.Vector3(0, 0, 1), rot_angle) * -1.0
        sources = [
            mp.EigenModeSource(
                src=mp.GaussianSource(fcen, fwidth=0.1),
                center=mp.Vector3(0.5 * sxy - 3.4, 0, 0),
                size=mp.Vector3(0, sxy, 0),
                direction=mp.NO_DIRECTION,
                eig_kpoint=kpoint,
                eig_band=1,
                eig_parity=mp.ODD_Z,
                eig_match_freq=True,
            )
        ]

        geometry = [
            mp.Block(
                center=mp.Vector3(),
                size=mp.Vector3(mp.inf, 1, mp.inf),
                e1=mp.Vector3(1, 0, 0).rotate(mp.Vector3(0, 0, 1), rot_angle),
                e2=mp.Vector3(0, 1, 0).rotate(mp.Vector3(0, 0, 1), rot_angle),
                material=mp.Medium(index=3.5),
            )
        ]

        sim = mp.Simulation(
            cell_size=cell_size,
            resolution=20,
            boundary_layers=pml_layers,
            sources=sources,
            geometry=geometry,
        )

        mode = sim.add_mode_monitor(
            fcen,
            0,
            1,
            mp.FluxRegion(
                center=mp.Vector3(-0.5 * sxy + dpml, 0, 0), size=mp.Vector3(0, sxy, 0)
            ),
            decimation_factor=1,
        )
        mode_decimated = sim.add_mode_monitor(
            fcen,
            0,
            1,
            mp.FluxRegion(
                center=mp.Vector3(-0.5 * sxy + dpml, 0, 0), size=mp.Vector3(0, sxy, 0)
            ),
            decimation_factor=10,
        )

        sim.run(until_after_sources=30)

        flux = mp.get_fluxes(mode)[0]
        coeff = sim.get_eigenmode_coefficients(
            mode, [1], direction=mp.NO_DIRECTION, kpoint_func=lambda f, n: kpoint
        ).alpha[0, 0, 0]
        flux_decimated = mp.get_fluxes(mode_decimated)[0]
        coeff_decimated = sim.get_eigenmode_coefficients(
            mode_decimated,
            [1],
            direction=mp.NO_DIRECTION,
            kpoint_func=lambda f, n: kpoint,
        ).alpha[0, 0, 0]

        print(f"oblique-waveguide-flux:, {-flux:.6f}, {abs(coeff) ** 2:.6f}")
        print(
            "oblique-waveguide-flux (decimated):, {:.6f}, {:.6f}".format(
                -flux_decimated, abs(coeff_decimated) ** 2
            )
        )
        ## the magnitude of |flux| is 100.008731 and so we check two significant digits of accuracy
        self.assertAlmostEqual(-1, abs(coeff) ** 2 / flux, places=2)
        self.assertAlmostEqual(flux, flux_decimated, places=3)
        self.assertAlmostEqual(coeff, coeff_decimated, places=3)

    def test_grating_3d(self):
        """Unit test for mode decomposition in 3d with zero k_point.

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
        fcen = 1 / wvl

        dpml = 1.0  # PML thickness
        dsub = 3.0  # substrate thickness
        dair = 3.0  # air padding
        hcyl = 0.5  # cylinder height
        rcyl = 0.2  # cylinder radius

        sx = 1.1
        sy = 0.8
        sz = dpml + dsub + hcyl + dair + dpml

        cell_size = mp.Vector3(sx, sy, sz)

        boundary_layers = [mp.PML(thickness=dpml, direction=mp.Z)]

        # periodic boundary conditions
        k_point = mp.Vector3()

        src_cmpt = mp.Ex
        sources = [
            mp.Source(
                src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
                size=mp.Vector3(sx, sy, 0),
                center=mp.Vector3(0, 0, -0.5 * sz + dpml),
                component=src_cmpt,
            )
        ]

        symmetries = [
            mp.Mirror(direction=mp.X, phase=-1),
            mp.Mirror(direction=mp.Y, phase=+1),
        ]

        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cell_size,
            sources=sources,
            default_material=SiO2,
            boundary_layers=boundary_layers,
            k_point=k_point,
            symmetries=symmetries,
        )

        refl_pt = mp.Vector3(0, 0, -0.5 * sz + dpml + 0.5 * dsub)
        refl_flux = sim.add_mode_monitor(
            fcen, 0, 1, mp.ModeRegion(center=refl_pt, size=mp.Vector3(sx, sy, 0))
        )

        stop_cond = mp.stop_when_energy_decayed(20, 1e-6)
        sim.run(until_after_sources=stop_cond)

        input_flux = mp.get_fluxes(refl_flux)
        input_flux_data = sim.get_flux_data(refl_flux)

        sim.reset_meep()

        geometry = [
            mp.Block(
                size=mp.Vector3(mp.inf, mp.inf, dpml + dsub),
                center=mp.Vector3(0, 0, -0.5 * sz + 0.5 * (dpml + dsub)),
                material=SiO2,
            ),
            mp.Cylinder(
                height=hcyl,
                radius=rcyl,
                center=mp.Vector3(0, 0, -0.5 * sz + dpml + dsub + 0.5 * hcyl),
                material=Si,
            ),
        ]

        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cell_size,
            sources=sources,
            geometry=geometry,
            boundary_layers=boundary_layers,
            k_point=k_point,
            symmetries=symmetries,
        )

        refl_flux = sim.add_mode_monitor(
            fcen, 0, 1, mp.ModeRegion(center=refl_pt, size=mp.Vector3(sx, sy, 0))
        )
        sim.load_minus_flux_data(refl_flux, input_flux_data)

        tran_flux = sim.add_mode_monitor(
            fcen,
            0,
            1,
            mp.ModeRegion(
                center=mp.Vector3(0, 0, 0.5 * sz - dpml), size=mp.Vector3(sx, sy, 0)
            ),
        )

        sim.run(until_after_sources=stop_cond)

        # sum the Poynting flux in z direction for all reflected orders
        Rsum = 0

        # number of reflected modes/orders in SiO2 in x and y directions (upper bound)
        nm_x = int(fcen * nSiO2 * sx) + 1
        nm_y = int(fcen * nSiO2 * sy) + 1
        for m_x in range(nm_x):
            for m_y in range(nm_y):
                for S_pol in [False, True]:
                    res = sim.get_eigenmode_coefficients(
                        refl_flux,
                        mp.DiffractedPlanewave(
                            [m_x, m_y, 0],
                            mp.Vector3(1, 0, 0),
                            1 if S_pol else 0,
                            0 if S_pol else 1,
                        ),
                    )
                    r_coeffs = res.alpha
                    Rmode = abs(r_coeffs[0, 0, 1]) ** 2 / input_flux[0]
                    print(
                        "refl-order:, {}, {}, {}, {:.6f}".format(
                            "s" if S_pol else "p", m_x, m_y, Rmode
                        )
                    )
                    if m_x == 0 and m_y == 0:
                        Rsum += Rmode
                    elif (m_x != 0 and m_y == 0) or (m_x == 0 and m_y != 0):
                        Rsum += 2 * Rmode
                    else:
                        Rsum += 4 * Rmode

        # sum the Poynting flux in z direction for all transmitted orders
        Tsum = 0

        # number of transmitted modes/orders in air in x and y directions (upper bound)
        nm_x = int(fcen * sx) + 1
        nm_y = int(fcen * sy) + 1
        for m_x in range(nm_x):
            for m_y in range(nm_y):
                for S_pol in [False, True]:
                    res = sim.get_eigenmode_coefficients(
                        tran_flux,
                        mp.DiffractedPlanewave(
                            [m_x, m_y, 0],
                            mp.Vector3(1, 0, 0),
                            1 if S_pol else 0,
                            0 if S_pol else 1,
                        ),
                    )
                    t_coeffs = res.alpha
                    Tmode = abs(t_coeffs[0, 0, 0]) ** 2 / input_flux[0]
                    print(
                        "tran-order:, {}, {}, {}, {:.6f}".format(
                            "s" if S_pol else "p", m_x, m_y, Tmode
                        )
                    )
                    if m_x == 0 and m_y == 0:
                        Tsum += Tmode
                    elif (m_x != 0 and m_y == 0) or (m_x == 0 and m_y != 0):
                        Tsum += 2 * Tmode
                    else:
                        Tsum += 4 * Tmode

        r_flux = mp.get_fluxes(refl_flux)
        t_flux = mp.get_fluxes(tran_flux)
        Rflux = -r_flux[0] / input_flux[0]
        Tflux = t_flux[0] / input_flux[0]

        print(f"refl:, {Rsum}, {Rflux}")
        print(f"tran:, {Tsum}, {Tflux}")
        print(f"sum:,  {Rsum + Tsum}, {Rflux + Tflux}")

        ## to obtain agreement for two decimal digits,
        ## the resolution must be increased to 200
        self.assertAlmostEqual(Rsum, Rflux, places=1)
        self.assertAlmostEqual(Tsum, Tflux, places=2)
        self.assertAlmostEqual(Rsum + Tsum, 1.00, places=1)

    def test_triangular_lattice_oblique(self):
        """Unit test for mode decomposition in 3d with nonzero k_point.

        Verifies that the sum of the diffraction efficiencies of all
        the reflected and transmitted orders of a binary grating with
        triangular lattice given an oblique planewave incident from
        within the high-index medium is equivalent to the reflectance and
        transmittance, respectively, obtained using the Poynting flux.
        """
        resolution = 30

        ng = 1.5
        glass = mp.Medium(index=ng)

        wvl = 0.5
        fcen = 1 / wvl

        dpml = 1.0
        dsub = 2.0
        dair = 2.0
        rcyl = 0.1
        hcyl = 0.3

        a = 0.6

        sx = a
        sy = a * np.sqrt(3)

        sz = dpml + dsub + hcyl + dair + dpml

        cell_size = mp.Vector3(sx, sy, sz)

        boundary_layers = [mp.PML(thickness=dpml, direction=mp.Z)]

        # plane of incidence is yz
        # 0° is +z with CCW rotation about x
        theta = math.radians(34.6)

        if theta == 0:
            k = mp.Vector3()
        else:
            # The planewave source is incident from within the high-index
            # medium which means ω = c|k|/n where n is the index of medium.
            # In Meep units (c=1), this implies |k| = nω.
            k = mp.Vector3(0, 0, ng * fcen).rotate(mp.Vector3(1, 0, 0), theta)

        def pw_amp(k, x0):
            def _pw_amp(x):
                return cmath.exp(1j * 2 * math.pi * k.dot(x + x0))

            return _pw_amp

        src_pt = mp.Vector3(0, 0, -0.5 * sz + dpml)
        src_cmpt = mp.Ex  # S-pol: Ex / P-pol: Ey
        sources = [
            mp.Source(
                src=mp.GaussianSource(fcen, fwidth=0.1 * fcen),
                size=mp.Vector3(sx, sy, 0),
                center=src_pt,
                component=src_cmpt,
                amp_func=pw_amp(k, src_pt),
            )
        ]

        symmetries = [mp.Mirror(direction=mp.X, phase=-1 if src_cmpt == mp.Ex else +1)]

        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cell_size,
            sources=sources,
            default_material=glass,
            boundary_layers=boundary_layers,
            k_point=k,
            symmetries=symmetries,
        )

        refl_pt = mp.Vector3(0, 0, -0.5 * sz + dpml + 0.5 * dsub)
        refl_flux = sim.add_mode_monitor(
            fcen, 0, 1, mp.ModeRegion(center=refl_pt, size=mp.Vector3(sx, sy, 0))
        )

        stop_cond = mp.stop_when_fields_decayed(25, src_cmpt, src_pt, 1e-6)
        sim.run(until_after_sources=stop_cond)

        input_flux = mp.get_fluxes(refl_flux)[0]
        input_flux_data = sim.get_flux_data(refl_flux)

        sim.reset_meep()

        substrate = [
            mp.Block(
                size=mp.Vector3(mp.inf, mp.inf, dpml + dsub),
                center=mp.Vector3(0, 0, -0.5 * sz + 0.5 * (dpml + dsub)),
                material=glass,
            )
        ]

        grating = [
            mp.Cylinder(
                center=mp.Vector3(0, 0, -0.5 * sz + dpml + dsub + 0.5 * hcyl),
                radius=rcyl,
                height=hcyl,
                material=glass,
            ),
            mp.Cylinder(
                center=mp.Vector3(
                    0.5 * sx, 0.5 * sy, -0.5 * sz + dpml + dsub + 0.5 * hcyl
                ),
                radius=rcyl,
                height=hcyl,
                material=glass,
            ),
            mp.Cylinder(
                center=mp.Vector3(
                    -0.5 * sx, 0.5 * sy, -0.5 * sz + dpml + dsub + 0.5 * hcyl
                ),
                radius=rcyl,
                height=hcyl,
                material=glass,
            ),
            mp.Cylinder(
                center=mp.Vector3(
                    0.5 * sx, -0.5 * sy, -0.5 * sz + dpml + dsub + 0.5 * hcyl
                ),
                radius=rcyl,
                height=hcyl,
                material=glass,
            ),
            mp.Cylinder(
                center=mp.Vector3(
                    -0.5 * sx, -0.5 * sy, -0.5 * sz + dpml + dsub + 0.5 * hcyl
                ),
                radius=rcyl,
                height=hcyl,
                material=glass,
            ),
        ]

        geometry = substrate + grating

        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cell_size,
            sources=sources,
            geometry=geometry,
            boundary_layers=boundary_layers,
            k_point=k,
            symmetries=symmetries,
        )

        refl_flux = sim.add_mode_monitor(
            fcen, 0, 1, mp.ModeRegion(center=refl_pt, size=mp.Vector3(sx, sy, 0))
        )

        sim.load_minus_flux_data(refl_flux, input_flux_data)

        tran_pt = mp.Vector3(0, 0, 0.5 * sz - dpml)
        tran_flux = sim.add_mode_monitor(
            fcen, 0, 1, mp.ModeRegion(center=tran_pt, size=mp.Vector3(sx, sy, 0))
        )

        sim.run(until_after_sources=stop_cond)

        Rsum = 0
        Tsum = 0
        m = 5
        tol = 1e-6
        for nx in range(-m, m + 1):
            for ny in range(-m, m + 1):
                # convert supercell order to unit cell order
                mx = nx
                my = (nx + ny) // 2

                # consider only propagating modes in high-index medium
                kz2 = (ng * fcen) ** 2 - (k.x + nx / sx) ** 2 - (k.y + ny / sy) ** 2
                if kz2 > 0:
                    Rpol = 0
                    for S_pol in [True, False]:
                        res = sim.get_eigenmode_coefficients(
                            refl_flux,
                            mp.DiffractedPlanewave(
                                (nx, ny, 0),
                                mp.Vector3(0, 1, 0),
                                1 if S_pol else 0,
                                0 if S_pol else 1,
                            ),
                        )

                        coeffs = res.alpha
                        refl = abs(coeffs[0, 0, 1]) ** 2 / input_flux

                        pol_str = "S" if S_pol else "P"

                        if refl > tol:
                            # determine whether diffracted order is for the unit cell or super cell
                            if (nx + ny) % 2 == 0:
                                Rpol += refl
                                print(
                                    "refl:, {}, {:2d}, {:2d}, {:.5f}, (unit cell)".format(
                                        pol_str, mx, my, refl
                                    )
                                )
                            else:
                                print(
                                    "refl:, {}, {:2d}, {:2d}, {:.7f}, (super cell)".format(
                                        pol_str, nx, ny, refl
                                    )
                                )

                    Rsum += Rpol

                # consider only propagating modes in air
                kz2 = fcen**2 - (k.x + nx / sx) ** 2 - (k.y + ny / sy) ** 2
                if kz2 > 0:
                    Tpol = 0
                    for S_pol in [True, False]:
                        res = sim.get_eigenmode_coefficients(
                            tran_flux,
                            mp.DiffractedPlanewave(
                                (nx, ny, 0),
                                mp.Vector3(0, 1, 0),
                                1 if S_pol else 0,
                                0 if S_pol else 1,
                            ),
                        )
                        coeffs = res.alpha
                        tran = abs(coeffs[0, 0, 0]) ** 2 / input_flux

                        pol_str = "S" if S_pol else "P"

                        if tran > tol:
                            # determine whether diffracted order is for the unit cell or super cell
                            if (nx + ny) % 2 == 0:
                                Tpol += tran
                                print(
                                    "tran:, {}, {:2d}, {:2d}, {:.5f}, (unit cell)".format(
                                        pol_str, mx, my, tran
                                    )
                                )
                            else:
                                print(
                                    "tran:, {}, {:2d}, {:2d}, {:.7f}, (super cell)".format(
                                        pol_str, nx, ny, tran
                                    )
                                )

                    Tsum += Tpol

        Rflux = -mp.get_fluxes(refl_flux)[0] / input_flux
        err = abs(Rflux - Rsum) / Rflux
        print(
            "refl:, {:.6f} (flux), {:.6f} (orders), {:.6f} (error)".format(
                Rflux, Rsum, err
            )
        )

        Tflux = mp.get_fluxes(tran_flux)[0] / input_flux
        err = abs(Tflux - Tsum) / Tflux
        print(
            "tran:, {:.6f} (flux), {:.6f} (orders), {:.6f} (error)".format(
                Tflux, Tsum, err
            )
        )

        self.assertAlmostEqual(Rsum, Rflux, places=3)
        self.assertAlmostEqual(Tsum, Tflux, places=3)
        self.assertAlmostEqual(Rsum + Tsum, 1.00, places=2)


if __name__ == "__main__":
    unittest.main()
