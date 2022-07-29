import unittest
import numpy as np
import meep as mp

# Computes the Purcell enhancement factor of a parallel dipole in a planar
# dielectric cavity with lossless metallic walls. The result is computed in
# cylindrical and 3D coordinates and validated using analytic theory from:
# I. Abram et al., IEEE J. Quantum Electronics, Vol. 34, pp. 71-76 (1998).


class TestLDOS(unittest.TestCase):
    @classmethod
    def setUp(cls):
        cls.resolution = 25  # pixels/μm
        cls.dpml = 0.5  # thickness of PML
        cls.L = 6.0  # length of non-PML region
        cls.n = 2.4  # refractive index of surrounding medium
        cls.wvl = 1.0  # wavelength (in vacuum)

        cls.fcen = 1 / cls.wvl

    def bulk_ldos_cyl(self):
        sr = self.L + self.dpml
        sz = self.L + 2 * self.dpml
        cell_size = mp.Vector3(sr, 0, sz)

        pml_layers = [mp.PML(self.dpml)]

        sources = [
            mp.Source(
                src=mp.GaussianSource(self.fcen, fwidth=0.2 * self.fcen),
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
                20, mp.Er, mp.Vector3(), 1e-6
            ),
        )

        return sim.ldos_data[0]

    def cavity_ldos_cyl(self, sz):
        sr = self.L + self.dpml
        cell_size = mp.Vector3(sr, 0, sz)

        pml_layers = [mp.PML(self.dpml, direction=mp.R)]

        sources = [
            mp.Source(
                src=mp.GaussianSource(self.fcen, fwidth=0.2 * self.fcen),
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
                20, mp.Er, mp.Vector3(), 1e-6
            ),
        )

        return sim.ldos_data[0]

    def bulk_ldos_3D(self):
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
                20, mp.Ex, mp.Vector3(), 1e-6
            ),
        )

        return sim.ldos_data[0]

    def cavity_ldos_3D(self, sz):
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
                20, mp.Ex, mp.Vector3(), 1e-6
            ),
        )

        return sim.ldos_data[0]

    def purcell_enh_theory(self, c):
        # equation 7 of reference
        return 3 * np.fix(c + 0.5) / (4 * c) + (
            4 * np.power(np.fix(c + 0.5), 3) - np.fix(c + 0.5)
        ) / (16 * np.power(c, 3))

    def test_ldos_cyl(self):
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


if __name__ == "__main__":
    unittest.main()
