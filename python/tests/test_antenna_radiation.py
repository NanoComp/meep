import math
import unittest

import numpy as np
from utils import ApproxComparisonTestCase

import meep as mp


class TestAntennaRadiation(ApproxComparisonTestCase):
    @classmethod
    def setUp(cls):
        cls.resolution = 100  # pixels/μm

        cls.h = 1.125  # height of point source above ground plane
        cls.n = 1.2  # refractive index of surrounding medium

        cls.src_cmpt = mp.Ez
        cls.wvl = 0.65

        cls.npts = 50  # number of points in [0,pi/2) range of angles
        cls.angles = 0.5 * math.pi / cls.npts * np.arange(cls.npts)
        cls.r = 1000 * cls.wvl  # radius of far-field hemicircle

    def radial_flux(self, sim, nearfield_box, r):
        E = np.zeros((self.npts, 3), dtype=np.complex128)
        H = np.zeros((self.npts, 3), dtype=np.complex128)
        for n in range(self.npts):
            ff = sim.get_farfield(
                nearfield_box,
                mp.Vector3(r * math.sin(self.angles[n]), r * math.cos(self.angles[n])),
            )
            E[n, :] = [np.conj(ff[j]) for j in range(3)]
            H[n, :] = [ff[j + 3] for j in range(3)]

        Px = np.real(E[:, 1] * H[:, 2] - E[:, 2] * H[:, 1])  # Ey*Hz-Ez*Hy
        Py = np.real(E[:, 2] * H[:, 0] - E[:, 0] * H[:, 2])  # Ez*Hx-Ex*Hz
        return np.sqrt(np.square(Px) + np.square(Py))

    def free_space_radiation(self):
        sxy = 4
        dpml = 1
        cell_size = mp.Vector3(sxy + 2 * dpml, sxy + 2 * dpml)

        pml_layers = [mp.PML(dpml)]

        fcen = 1 / self.wvl

        sources = [
            mp.Source(
                src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
                center=mp.Vector3(),
                component=self.src_cmpt,
            )
        ]

        if self.src_cmpt == mp.Hz:
            symmetries = [mp.Mirror(mp.X, phase=-1), mp.Mirror(mp.Y, phase=-1)]
        elif self.src_cmpt == mp.Ez:
            symmetries = [mp.Mirror(mp.X, phase=+1), mp.Mirror(mp.Y, phase=+1)]
        else:
            symmetries = []

        sim = mp.Simulation(
            cell_size=cell_size,
            resolution=self.resolution,
            sources=sources,
            symmetries=symmetries,
            boundary_layers=pml_layers,
            default_material=mp.Medium(index=self.n),
        )

        nearfield_box = sim.add_near2far(
            fcen,
            0,
            1,
            mp.Near2FarRegion(
                center=mp.Vector3(0, +0.5 * sxy), size=mp.Vector3(sxy, 0)
            ),
            mp.Near2FarRegion(
                center=mp.Vector3(0, -0.5 * sxy), size=mp.Vector3(sxy, 0), weight=-1
            ),
            mp.Near2FarRegion(
                center=mp.Vector3(+0.5 * sxy, 0), size=mp.Vector3(0, sxy)
            ),
            mp.Near2FarRegion(
                center=mp.Vector3(-0.5 * sxy, 0), size=mp.Vector3(0, sxy), weight=-1
            ),
        )

        sim.run(until_after_sources=mp.stop_when_dft_decayed())

        return self.radial_flux(sim, nearfield_box, self.r)

    def pec_ground_plane_radiation(self):
        L = 8.0  # length of non-PML region
        dpml = 1.0  # thickness of PML
        sxy = dpml + L + dpml
        cell_size = mp.Vector3(sxy, sxy, 0)

        boundary_layers = [mp.PML(dpml)]

        fcen = 1 / self.wvl

        # The near-to-far field transformation feature only supports
        # homogeneous media which means it cannot explicitly take into
        # account the ground plane. As a workaround, we use two antennas
        # of _opposite_ sign surrounded by a single near2far box which
        # encloses both antennas. We then use an odd mirror symmetry to
        # cut the computational cell in half which is effectively
        # equivalent to a PEC boundary condition on one side.
        # Note: This setup means that the radiation pattern can only
        # be measured in the top half above the dipole.
        sources = [
            mp.Source(
                src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
                component=self.src_cmpt,
                center=mp.Vector3(0, +self.h),
            ),
            mp.Source(
                src=mp.GaussianSource(fcen, fwidth=0.2 * fcen),
                component=self.src_cmpt,
                center=mp.Vector3(0, -self.h),
                amplitude=-1 if self.src_cmpt == mp.Ez else +1,
            ),
        ]

        symmetries = [
            mp.Mirror(direction=mp.X, phase=+1 if self.src_cmpt == mp.Ez else -1),
            mp.Mirror(direction=mp.Y, phase=-1 if self.src_cmpt == mp.Ez else +1),
        ]

        sim = mp.Simulation(
            resolution=self.resolution,
            cell_size=cell_size,
            boundary_layers=boundary_layers,
            sources=sources,
            symmetries=symmetries,
            default_material=mp.Medium(index=self.n),
        )

        nearfield_box = sim.add_near2far(
            fcen,
            0,
            1,
            mp.Near2FarRegion(
                center=mp.Vector3(0, 2 * self.h), size=mp.Vector3(2 * self.h, 0)
            ),
            mp.Near2FarRegion(
                center=mp.Vector3(0, -2 * self.h),
                size=mp.Vector3(2 * self.h, 0),
                weight=-1,
            ),
            mp.Near2FarRegion(
                center=mp.Vector3(self.h, 0), size=mp.Vector3(0, 4 * self.h)
            ),
            mp.Near2FarRegion(
                center=mp.Vector3(-self.h, 0), size=mp.Vector3(0, 4 * self.h), weight=-1
            ),
        )

        sim.run(until_after_sources=mp.stop_when_dft_decayed())

        return self.radial_flux(sim, nearfield_box, self.r)

    def test_pec_ground_plane(self):
        """Unit test for near-to-far field transformation and symmetries.

        Verifies that the radiation pattern for a point dipole source a
        given height above a perfect-electric conductor (PEC) ground plane
        agrees with the theoretical result.

        The radiation pattern of a two-element antenna array is equivalent
        to the radiation pattern of a single antenna multiplied by its array
        factor as derived in Section 6.2 "Two-Element Array" of Antenna Theory:
        Analysis and Design, Fourth Edition (2016) by C.A. Balanis.
        """
        Pr_fsp = self.free_space_radiation()
        Pr_pec = self.pec_ground_plane_radiation()

        k = 2 * np.pi / (self.wvl / self.n)  # wavevector in medium
        Pr_theory = np.zeros(
            self.npts,
        )
        for i, ang in enumerate(self.angles):
            Pr_theory[i] = Pr_fsp[i] * 2 * np.sin(k * self.h * np.cos(ang))

        Pr_pec_norm = Pr_pec / np.max(Pr_pec)
        Pr_theory_norm = (Pr_theory / max(Pr_theory)) ** 2

        tol = 0.02
        self.assertClose(Pr_pec_norm, Pr_theory_norm, epsilon=tol)

    def test_poynting_theorem(self):
        """Unit test for near-to-far field transformation in 2d.

        Verifies that the Poynting flux of an Ez-polarized point
        dipole source in vacuum is independent of the shape of the
        enclosing measurement box due to Poynting's theorem by
        considering three arrangements:
             (1) bounding box of thenear fields,
             (2) bounding circle of the far fields, and
             (3) bounding box of the far fields.
        """
        resolution = 50
        sxy = 4
        dpml = 1
        cell = mp.Vector3(sxy + 2 * dpml, sxy + 2 * dpml)

        pml_layers = mp.PML(dpml)

        fcen = 1.0
        df = 0.4

        sources = [
            mp.Source(
                src=mp.GaussianSource(fcen, fwidth=df),
                center=mp.Vector3(),
                component=mp.Ez,
            )
        ]

        symmetries = [mp.Mirror(mp.X), mp.Mirror(mp.Y)]

        sim = mp.Simulation(
            cell_size=cell,
            resolution=resolution,
            sources=sources,
            symmetries=symmetries,
            boundary_layers=[pml_layers],
        )

        nearfield_box = sim.add_near2far(
            fcen,
            0,
            1,
            mp.Near2FarRegion(mp.Vector3(y=0.5 * sxy), size=mp.Vector3(sxy)),
            mp.Near2FarRegion(
                mp.Vector3(y=-0.5 * sxy), size=mp.Vector3(sxy), weight=-1
            ),
            mp.Near2FarRegion(mp.Vector3(0.5 * sxy), size=mp.Vector3(y=sxy)),
            mp.Near2FarRegion(
                mp.Vector3(-0.5 * sxy), size=mp.Vector3(y=sxy), weight=-1
            ),
        )

        flux_box = sim.add_flux(
            fcen,
            0,
            1,
            mp.FluxRegion(mp.Vector3(y=0.5 * sxy), size=mp.Vector3(sxy)),
            mp.FluxRegion(mp.Vector3(y=-0.5 * sxy), size=mp.Vector3(sxy), weight=-1),
            mp.FluxRegion(mp.Vector3(0.5 * sxy), size=mp.Vector3(y=sxy)),
            mp.FluxRegion(mp.Vector3(-0.5 * sxy), size=mp.Vector3(y=sxy), weight=-1),
        )

        sim.run(
            until_after_sources=mp.stop_when_fields_decayed(
                50, mp.Ez, mp.Vector3(), 1e-8
            )
        )

        near_flux = mp.get_fluxes(flux_box)[0]

        r = 1000 / fcen  # radius of far field circle
        Pr = self.radial_flux(sim, nearfield_box, r)
        far_flux_circle = 4 * np.sum(Pr) * 0.5 * np.pi * r / len(Pr)

        rr = 20 / fcen  # length of far-field square box
        res_far = 20  # resolution of far-field square box
        far_flux_square = (
            nearfield_box.flux(
                mp.Y,
                mp.Volume(center=mp.Vector3(y=0.5 * rr), size=mp.Vector3(rr)),
                res_far,
            )[0]
            - nearfield_box.flux(
                mp.Y,
                mp.Volume(center=mp.Vector3(y=-0.5 * rr), size=mp.Vector3(rr)),
                res_far,
            )[0]
            + nearfield_box.flux(
                mp.X,
                mp.Volume(center=mp.Vector3(0.5 * rr), size=mp.Vector3(y=rr)),
                res_far,
            )[0]
            - nearfield_box.flux(
                mp.X,
                mp.Volume(center=mp.Vector3(-0.5 * rr), size=mp.Vector3(y=rr)),
                res_far,
            )[0]
        )

        print(
            "flux:, {:.6f}, {:.6f}, {:.6f}".format(
                near_flux, far_flux_circle, far_flux_square
            )
        )

        self.assertAlmostEqual(near_flux, far_flux_circle, places=2)
        self.assertAlmostEqual(far_flux_circle, far_flux_square, places=2)
        self.assertAlmostEqual(far_flux_square, near_flux, places=2)


if __name__ == "__main__":
    unittest.main()
