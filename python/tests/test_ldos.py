import unittest
import numpy as np
import meep as mp


class TestLDOS(unittest.TestCase):
    @classmethod
    def setUp(cls):
        cls.resolution = 25  # pixels/μm
        cls.dpml = 0.5  # thickness of PML
        cls.dair = 1.0  # thickness of air padding
        cls.L = 6.0  # length of non-PML region
        cls.n = 2.4  # refractive index of surrounding medium
        cls.wvl = 1.0  # wavelength (in vacuum)

        cls.fcen = 1 / cls.wvl

        # termination criteria
        cls.tol = 1e-8

    def bulk_ldos_cyl(self):
        """Computes the LDOS of a point dipole in a homogeneous dielectric
        medium in cylindrical coordinates.
        """
        sr = self.L + self.dpml
        sz = self.L + 2 * self.dpml
        cell_size = mp.Vector3(sr, 0, sz)

        pml_layers = [mp.PML(self.dpml)]

        sources = [
            mp.Source(
                src=mp.GaussianSource(self.fcen, fwidth=0.1 * self.fcen),
                component=mp.Er,
                center=mp.Vector3(),
            )
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=cell_size,
            boundary_layers=pml_layers,
            sources=sources,
            dimensions=mp.CYLINDRICAL,
            m=-1,
            default_material=mp.Medium(index=self.n),
        )

        sim.run(
            mp.dft_ldos(self.fcen, 0, 1),
            until_after_sources=mp.stop_when_fields_decayed(
                20, mp.Er, mp.Vector3(), self.tol
            ),
        )

        return sim.ldos_data[0]

    def cavity_ldos_cyl(self, sz):
        """Computes the LDOS of a point dipole in a planar cavity with
        lossless metallic walls in cylindrical coordinates.

        Args:
          sz: thickness of cavity.
        """
        sr = self.L + self.dpml
        cell_size = mp.Vector3(sr, 0, sz)

        pml_layers = [mp.PML(self.dpml, direction=mp.R)]

        sources = [
            mp.Source(
                src=mp.GaussianSource(self.fcen, fwidth=0.1 * self.fcen),
                component=mp.Er,
                center=mp.Vector3(),
            )
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=cell_size,
            boundary_layers=pml_layers,
            sources=sources,
            dimensions=mp.CYLINDRICAL,
            m=-1,
            default_material=mp.Medium(index=self.n),
        )

        sim.run(
            mp.dft_ldos(self.fcen, 0, 1),
            until_after_sources=mp.stop_when_fields_decayed(
                20, mp.Er, mp.Vector3(), self.tol
            ),
        )

        return sim.ldos_data[0]

    def bulk_ldos_3D(self):
        """Computes the LDOS of a point dipole in a homogeneous dielectric
        medium in 3D Cartesian coordinates.
        """
        s = self.L + 2 * self.dpml
        cell_size = mp.Vector3(s, s, s)

        pml_layers = [mp.PML(self.dpml)]

        sources = [
            mp.Source(
                src=mp.GaussianSource(self.fcen, fwidth=0.2 * self.fcen),
                component=mp.Ex,
                center=mp.Vector3(),
            )
        ]

        symmetries = [
            mp.Mirror(direction=mp.X, phase=-1),
            mp.Mirror(direction=mp.Y),
            mp.Mirror(direction=mp.Z),
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=cell_size,
            boundary_layers=pml_layers,
            sources=sources,
            symmetries=symmetries,
            default_material=mp.Medium(index=self.n),
        )

        sim.run(
            mp.dft_ldos(self.fcen, 0, 1),
            until_after_sources=mp.stop_when_fields_decayed(
                20, mp.Ex, mp.Vector3(), self.tol
            ),
        )

        return sim.ldos_data[0]

    def cavity_ldos_3D(self, sz):
        """Computes the LDOS of a point dipole in a planar cavity with
        lossless metallic walls in 3D Cartesian coordinates.

        Args:
          sz: thickness of cavity.
        """
        sxy = self.L + 2 * self.dpml
        cell_size = mp.Vector3(sxy, sxy, sz)

        boundary_layers = [
            mp.PML(self.dpml, direction=mp.X),
            mp.PML(self.dpml, direction=mp.Y),
        ]

        sources = [
            mp.Source(
                src=mp.GaussianSource(self.fcen, fwidth=0.2 * self.fcen),
                component=mp.Ex,
                center=mp.Vector3(),
            )
        ]

        symmetries = [
            mp.Mirror(direction=mp.X, phase=-1),
            mp.Mirror(direction=mp.Y),
            mp.Mirror(direction=mp.Z),
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=cell_size,
            boundary_layers=boundary_layers,
            sources=sources,
            symmetries=symmetries,
            default_material=mp.Medium(index=self.n),
        )

        sim.run(
            mp.dft_ldos(ldos=mp.Ldos(self.fcen, 0, 1)),
            until_after_sources=mp.stop_when_fields_decayed(
                20, mp.Ex, mp.Vector3(), self.tol
            ),
        )

        return sim.ldos_data[0]

    def purcell_enh_theory(self, c):
        """Computes the Purcell enhancement factor of a point dipole in a planar
        dielectric cavity with lossless metallic walls using equation 7 of:
        I. Abram et al., IEEE J. Quantum Electronics, Vol. 34, pp. 71-76 (1998).

        Args:
          c: cavity thickness in units of wavelength in the cavity medium.
        """
        return 3 * np.fix(c + 0.5) / (4 * c) + (
            4 * np.power(np.fix(c + 0.5), 3) - np.fix(c + 0.5)
        ) / (16 * np.power(c, 3))

    def ext_eff_cyl(self, dmat, h, m):
        """Computes the extraction efficiency of a point dipole embedded
        within a dielectric layer above a lossless ground plane in
        cylindrical coordinates.

        Args:
          dmat: thickness of dielectric layer.
          h: height of dipole above ground plane as fraction of dmat.
          m: angular dependence of the fields, exp(imφ).
        """
        sr = self.L + self.dpml
        sz = dmat + self.dair + self.dpml
        cell_size = mp.Vector3(sr, 0, sz)

        boundary_layers = [
            mp.PML(self.dpml, direction=mp.R),
            mp.PML(self.dpml, direction=mp.Z, side=mp.High),
        ]

        src_cmpt = mp.Er
        src_pt = mp.Vector3(0, 0, -0.5 * sz + h * dmat)
        sources = [
            mp.Source(
                src=mp.GaussianSource(self.fcen, fwidth=0.1 * self.fcen),
                component=src_cmpt,
                center=src_pt,
            ),
        ]

        geometry = [
            mp.Block(
                material=mp.Medium(index=self.n),
                center=mp.Vector3(0, 0, -0.5 * sz + 0.5 * dmat),
                size=mp.Vector3(mp.inf, mp.inf, dmat),
            )
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=cell_size,
            dimensions=mp.CYLINDRICAL,
            m=m,
            boundary_layers=boundary_layers,
            sources=sources,
            geometry=geometry,
        )

        flux_air = sim.add_flux(
            self.fcen,
            0,
            1,
            mp.FluxRegion(
                center=mp.Vector3(0.5 * self.L, 0, 0.5 * sz - self.dpml),
                size=mp.Vector3(self.L, 0, 0),
            ),
            mp.FluxRegion(
                center=mp.Vector3(self.L, 0, 0.5 * sz - self.dpml - 0.5 * self.dair),
                size=mp.Vector3(0, 0, self.dair),
            ),
        )

        sim.run(
            mp.dft_ldos(self.fcen, 0, 1),
            until_after_sources=mp.stop_when_fields_decayed(
                20, src_cmpt, src_pt, self.tol
            ),
        )

        out_flux = mp.get_fluxes(flux_air)[0]
        dV = np.pi / (self.resolution**3)
        total_flux = -np.real(sim.ldos_Fdata[0] * np.conj(sim.ldos_Jdata[0])) * dV
        ext_eff = out_flux / total_flux
        print(f"extraction efficiency (cyl):, " f"{dmat:.4f}, {h:.4f}, {ext_eff:.6f}")

        return ext_eff

    def ext_eff_3D(self, dmat, h):
        """Computes the extraction efficiency of a point dipole embedded
        within a dielectric layer above a lossless ground plane in
        3D Cartesian coordinates.

        Args:
          dmat: thickness of dielectric layer.
          h: height of dipole above ground plane as fraction of dmat.
        """
        sxy = self.L + 2 * self.dpml
        sz = dmat + self.dair + self.dpml
        cell_size = mp.Vector3(sxy, sxy, sz)

        symmetries = [mp.Mirror(direction=mp.X, phase=-1), mp.Mirror(direction=mp.Y)]

        boundary_layers = [
            mp.PML(self.dpml, direction=mp.X),
            mp.PML(self.dpml, direction=mp.Y),
            mp.PML(self.dpml, direction=mp.Z, side=mp.High),
        ]

        src_cmpt = mp.Ex
        src_pt = mp.Vector3(0, 0, -0.5 * sz + h * dmat)
        sources = [
            mp.Source(
                src=mp.GaussianSource(self.fcen, fwidth=0.1 * self.fcen),
                component=src_cmpt,
                center=src_pt,
            )
        ]

        geometry = [
            mp.Block(
                material=mp.Medium(index=self.n),
                center=mp.Vector3(0, 0, -0.5 * sz + 0.5 * dmat),
                size=mp.Vector3(mp.inf, mp.inf, dmat),
            )
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=cell_size,
            boundary_layers=boundary_layers,
            sources=sources,
            geometry=geometry,
            symmetries=symmetries,
        )

        flux_air = sim.add_flux(
            self.fcen,
            0,
            1,
            mp.FluxRegion(
                center=mp.Vector3(0, 0, 0.5 * sz - self.dpml),
                size=mp.Vector3(self.L, self.L, 0),
            ),
            mp.FluxRegion(
                center=mp.Vector3(
                    0.5 * self.L, 0, 0.5 * sz - self.dpml - 0.5 * self.dair
                ),
                size=mp.Vector3(0, self.L, self.dair),
            ),
            mp.FluxRegion(
                center=mp.Vector3(
                    -0.5 * self.L, 0, 0.5 * sz - self.dpml - 0.5 * self.dair
                ),
                size=mp.Vector3(0, self.L, self.dair),
                weight=-1.0,
            ),
            mp.FluxRegion(
                center=mp.Vector3(
                    0, 0.5 * self.L, 0.5 * sz - self.dpml - 0.5 * self.dair
                ),
                size=mp.Vector3(self.L, 0, self.dair),
            ),
            mp.FluxRegion(
                center=mp.Vector3(
                    0, -0.5 * self.L, 0.5 * sz - self.dpml - 0.5 * self.dair
                ),
                size=mp.Vector3(self.L, 0, self.dair),
                weight=-1.0,
            ),
        )

        sim.run(
            mp.dft_ldos(self.fcen, 0, 1),
            until_after_sources=mp.stop_when_fields_decayed(
                20, src_cmpt, src_pt, self.tol
            ),
        )

        out_flux = mp.get_fluxes(flux_air)[0]
        dV = 1 / (self.resolution**3)
        total_flux = -np.real(sim.ldos_Fdata[0] * np.conj(sim.ldos_Jdata[0])) * dV
        ext_eff = out_flux / total_flux
        print(f"extraction efficiency (3D):, " f"{dmat:.4f}, {h:.4f}, {ext_eff:.6f}")

        return ext_eff

    def test_ldos_cyl(self):
        """Verifies that the Purcell enhancement factor of a parallel dipole
        in a planar dielectric cavity with lossless metallic walls computed in
        cylindrical coordinates agrees with the analytic result.
        """

        ldos_bulk = self.bulk_ldos_cyl()

        # not a Van Hove singularity
        cavity_thickness = 1.63
        gap = cavity_thickness * self.wvl / self.n

        ldos_cavity = self.cavity_ldos_cyl(gap)

        # Purcell enhancement factor (relative to bulk medium)
        pe_meep = ldos_cavity / ldos_bulk

        pe_theory = self.purcell_enh_theory(cavity_thickness)

        rel_err = abs(pe_meep - pe_theory) / pe_theory

        print(
            "ldos-cyl:, {:.6f} (Meep), {:.6f} (theory), "
            "{:.6f} (error)".format(pe_meep, pe_theory, rel_err)
        )

        self.assertAlmostEqual(pe_meep, pe_theory, delta=0.1)

    def test_ldos_3D(self):
        """Verifies that the Purcell enhancement factor of a parallel dipole
        in a planar dielectric cavity with lossless metallic walls computed in
        3D Cartesian coordinates agrees with the analytic result.
        """

        ldos_bulk = self.bulk_ldos_3D()

        # not a Van Hove singularity
        cavity_thickness = 0.75
        gap = cavity_thickness * self.wvl / self.n

        ldos_cavity = self.cavity_ldos_3D(gap)

        # Purcell enhancement factor (relative to bulk medium)
        pe_meep = ldos_cavity / ldos_bulk

        pe_theory = self.purcell_enh_theory(cavity_thickness)

        rel_err = abs(pe_meep - pe_theory) / pe_theory

        print(
            "ldos-3D:, {:.6f} (Meep), {:.6f} (theory), "
            "{:.6f} (error)".format(pe_meep, pe_theory, rel_err)
        )

        self.assertAlmostEqual(pe_meep, pe_theory, delta=0.1)

    def test_ldos_ext_eff(self):
        """Verifies that the extraction efficiency of a point dipole in a
        dielectric layer above a lossless ground plane computed in cylindrical
        coordinates (for m=±1, separately) and 3D Cartesian agree.
        """
        layer_thickness = 0.5 * self.wvl / self.n
        dipole_height = 0.5

        ext_eff_cyl = self.ext_eff_cyl(layer_thickness, dipole_height, -1.0)
        ext_eff_3D = self.ext_eff_3D(layer_thickness, dipole_height)

        self.assertAlmostEqual(ext_eff_cyl, ext_eff_3D, delta=0.01)

        ext_eff_cyl_m_plus = self.ext_eff_cyl(
            layer_thickness,
            dipole_height,
            +1.0,
        )
        self.assertEqual(ext_eff_cyl, ext_eff_cyl_m_plus)


if __name__ == "__main__":
    unittest.main()
