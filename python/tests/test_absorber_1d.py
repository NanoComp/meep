import unittest

from meep.materials import Al
from utils import ApproxComparisonTestCase

import meep as mp


def track_max_abs(component, pt, record):
    """Returns a step function recording the maximum of |Re(component)| at pt."""
    record["max_abs"] = 0.0

    def _track(sim):
        field = sim.get_field_point(component, pt).real
        record["max_abs"] = max(record["max_abs"], abs(field))

    return _track


def record_field(component, pt, record, key):
    """Returns a step function storing Re(component) at pt in record[key]."""

    def _record(sim):
        record[key] = sim.get_field_point(component, pt).real

    return _record


class TestAbsorber(ApproxComparisonTestCase):
    def setUp(self):

        resolution = 40
        cell_size = mp.Vector3(z=10)

        absorber_layers = [mp.Absorber(1, direction=mp.Z)]

        sources = [
            mp.Source(
                src=mp.GaussianSource(1 / 0.803, fwidth=0.1),
                center=mp.Vector3(),
                component=mp.Ex,
            )
        ]

        self.sim = mp.Simulation(
            cell_size=cell_size,
            resolution=resolution,
            dimensions=1,
            default_material=Al,
            boundary_layers=absorber_layers,
            sources=sources,
        )

    def test_absorber(self):
        """An absorber overlapping a dispersive (Drude-Lorentz) medium.

        Al is a metal at this wavelength, so the fields are confined near the
        source and this does not test the absorption itself (that is done by
        test_absorber_1d_vacuum and test_absorber_2d). Rather, it checks that
        the fields are reproducible and decay, i.e. remain stable.
        """
        pt = mp.Vector3()
        record = {}
        self.sim.run(
            track_max_abs(mp.Ex, pt, record),
            mp.at_time(50, record_field(mp.Ex, pt, record, "f50")),
            until_after_sources=mp.stop_when_fields_decayed(50, mp.Ex, pt, 1e-6),
        )

        tol = 1e-4 if mp.is_single_precision() else 1e-6
        self.assertClose(record["f50"], -1.710379541e-02, epsilon=tol)
        # The source turns off at t=100, so the fields must have decayed by
        # the first check after that.
        self.assertAlmostEqual(self.sim.meep_time(), 150.0375, places=6)
        f_end = self.sim.get_field_point(mp.Ex, pt).real
        self.assertLess(abs(f_end), 1e-6 * record["max_abs"])

    def test_absorber_1d_vacuum(self):
        """A 1d absorber in vacuum, checked against reference values.

        The fields are stepped until (after t=50) the field at the source has
        decayed to 1e-6 of its maximum. The field at t=50 and the final time
        and field are then compaerd with reference values, which makes this
        a sensitive regression test of the absorber.
        """
        pt = mp.Vector3()
        sim = mp.Simulation(
            cell_size=mp.Vector3(z=10),
            resolution=20,
            dimensions=1,
            boundary_layers=[mp.Absorber(1, direction=mp.Z)],
            sources=[
                mp.Source(
                    # is_integrated=True for consistency with the reference
                    # values, which were computed using the default source
                    # time convention of the C++ interface.
                    src=mp.GaussianSource(1 / 0.803, fwidth=0.1, is_integrated=True),
                    center=pt,
                    component=mp.Ex,
                )
            ],
        )

        record = {"max_abs": 0.0, "f50": None}

        def fields_decayed(sim):
            t = sim.meep_time()
            field = sim.get_field_point(mp.Ex, pt).real
            if t >= 50 and record["f50"] is None:
                record["f50"] = field
            record["max_abs"] = max(record["max_abs"], abs(field))
            record["field"] = field
            return t > 50 and abs(field) <= 1e-6 * record["max_abs"]

        sim.run(until=fields_decayed)

        f_final_ref = 1.413336e-07 if mp.is_single_precision() else 1.624782e-07
        self.assertClose(record["f50"], -5.090114e-01, epsilon=1e-6)
        self.assertClose(sim.meep_time(), 9.195000e01, epsilon=1e-6)
        self.assertClose(record["field"], f_final_ref, epsilon=1e-6)

    def test_absorber_2d(self):
        source = mp.Source(
            src=mp.GaussianSource(frequency=0.1, fwidth=0.1),
            component=mp.Hz,
            center=mp.Vector3(),
        )

        sim = mp.Simulation(
            cell_size=mp.Vector3(20, 20, 0),
            resolution=10,
            sources=[source],
            boundary_layers=[mp.Absorber(5)],
        )

        v = mp.Vector3(4.13, 3.75, 0)
        record = {}
        sim.run(
            track_max_abs(mp.Hz, v, record),
            mp.at_time(50, record_field(mp.Hz, v, record, "f50")),
            until_after_sources=1000,
        )

        tol = 1e-5 if mp.is_single_precision() else 1e-6
        self.assertClose(record["f50"], 5.334524937e-02, epsilon=tol)
        # Long after the source has turned off, the absorber must have
        # removed (nearly) all of the energy.
        p = sim.get_field_point(mp.Hz, v).real
        self.assertAlmostEqual(abs(p), 1e-7 * record["max_abs"])


if __name__ == "__main__":
    unittest.main()
