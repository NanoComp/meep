import unittest

import meep as mp

# compute group velocity of a waveguide mode using two different methods
# (1) ratio of Poynting flux to energy density
# (2) via MPB from get-eigenmode-coefficients


class TestDftEnergy(unittest.TestCase):
    def test_dft_energy(self):
        resolution = 20
        cell = mp.Vector3(10, 5)
        geom = [
            mp.Block(size=mp.Vector3(mp.inf, 1, mp.inf), material=mp.Medium(epsilon=12))
        ]
        pml = [mp.PML(1)]
        fsrc = 0.15
        sources = [
            mp.EigenModeSource(
                src=mp.GaussianSource(frequency=fsrc, fwidth=0.2 * fsrc),
                center=mp.Vector3(-3),
                size=mp.Vector3(y=5),
                eig_band=1,
                eig_parity=mp.ODD_Z + mp.EVEN_Y,
                eig_match_freq=True,
            )
        ]

        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cell,
            geometry=geom,
            boundary_layers=pml,
            sources=sources,
            symmetries=[mp.Mirror(direction=mp.Y)],
        )

        flux = sim.add_flux(
            fsrc,
            0,
            1,
            mp.FluxRegion(center=mp.Vector3(3), size=mp.Vector3(y=5)),
            decimation_factor=1,
        )
        energy = sim.add_energy(
            fsrc,
            0,
            1,
            mp.EnergyRegion(center=mp.Vector3(3), size=mp.Vector3(y=5)),
            decimation_factor=1,
        )
        energy_decimated = sim.add_energy(
            fsrc,
            0,
            1,
            mp.EnergyRegion(center=mp.Vector3(3), size=mp.Vector3(y=5)),
            decimation_factor=10,
        )
        sim.run(until_after_sources=100)

        res = sim.get_eigenmode_coefficients(flux, [1], eig_parity=mp.ODD_Z + mp.EVEN_Y)
        mode_vg = res.vgrp[0]
        poynting_flux = mp.get_fluxes(flux)[0]
        e_energy = mp.get_electric_energy(energy)[0]
        ratio_vg = (0.5 * poynting_flux) / e_energy
        m_energy = mp.get_magnetic_energy(energy)[0]
        t_energy = mp.get_total_energy(energy)[0]

        self.assertAlmostEqual(m_energy + e_energy, t_energy)
        self.assertAlmostEqual(ratio_vg, mode_vg, places=3)

        e_energy_decimated = mp.get_electric_energy(energy_decimated)[0]
        m_energy_decimated = mp.get_magnetic_energy(energy_decimated)[0]
        self.assertAlmostEqual(e_energy, e_energy_decimated, places=1)
        self.assertAlmostEqual(m_energy, m_energy_decimated, places=1)


if __name__ == "__main__":
    unittest.main()
