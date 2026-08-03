import unittest

import meep as mp
from meep import materials
from meep.materials import (
    Ag,
    Cr,
    Ge,
    InP,
    LiNbO3,
    HfO2,
    MgF2,
    Si,
    SiO2_aniso,
    Ta2O5,
    TiO2,
    get_material,
)


class TestMaterialsLibrary(unittest.TestCase):
    def test_materials_library(self):
        self.assertAlmostEqual(InP.epsilon(1 / 3.3)[0][0], (3.1031) ** 2, places=2)

        self.assertAlmostEqual(Ge.epsilon(1 / 6.8)[0][0], (4.0091) ** 2, places=2)

        self.assertAlmostEqual(Si.epsilon(1 / 1.55)[0][0], (3.4777) ** 2, places=2)

        self.assertAlmostEqual(LiNbO3.epsilon(1 / 1.55)[0][0], (2.2111) ** 2, places=2)
        self.assertAlmostEqual(LiNbO3.epsilon(1 / 1.55)[1][1], (2.2111) ** 2, places=2)
        self.assertAlmostEqual(LiNbO3.epsilon(1 / 1.55)[2][2], (2.1376) ** 2, places=2)

        self.assertAlmostEqual(
            SiO2_aniso.epsilon(1 / 1.55)[0][0], (1.5277) ** 2, places=2
        )
        self.assertEqual(SiO2_aniso.epsilon(1 / 1.55)[1][0], 0)
        self.assertAlmostEqual(
            SiO2_aniso.epsilon(1 / 1.55)[1][1], (1.5277) ** 2, places=2
        )
        self.assertAlmostEqual(
            SiO2_aniso.epsilon(1 / 1.55)[2][2], (1.5362) ** 2, places=2
        )

        self.assertAlmostEqual(
            Ag.epsilon(1 / 0.65)[0][0], (0.14623 + 1j * 3.9367) ** 2, places=2
        )

        self.assertAlmostEqual(
            Cr.epsilon(1 / 0.71)[0][0], (3.8275 + 1j * 4.3457) ** 2, places=2
        )

        # amorphous TiO2 (Sarkar 2019), lossless over the visible; n(532 nm) ≈ 2.1747
        self.assertAlmostEqual(TiO2.epsilon(1 / 0.532)[0][0], (2.1747) ** 2, places=2)

        # amorphous Ta2O5 (Rodriguez-de Marcos 2016), high-index coating material
        # with a weak absorption tail; at 600 nm n ≈ 2.102, k ≈ 0.0083.
        self.assertAlmostEqual(
            Ta2O5.epsilon(1 / 0.6)[0][0].real, (2.102) ** 2, places=2
        )
        self.assertAlmostEqual(
            Ta2O5.epsilon(1 / 0.6)[0][0].imag, 2 * 2.102 * 0.0083, places=2
        )

        # amorphous HfO2 (Franta 2015), high-index UV/DUV coating materials:
        # transparent in the visible (n(633 nm) ≈ 2.088) with a deep-UV
        # absorption edge (k(220 nm) ≈ 0.0385).
        self.assertAlmostEqual(
            HfO2.epsilon(1 / 0.633)[0][0].real, (2.0877) ** 2, places=2
        )
        self.assertAlmostEqual(
            HfO2.epsilon(1 / 0.22)[0][0].imag, 2 * 2.6141 * 0.03853, places=2
        )

        # MgF2 film (Franta 2017), the standard low-index AR-coating material;
        # lossless across the visible with n(550 nm) ≈ 1.3837.
        self.assertAlmostEqual(MgF2.epsilon(1 / 0.55)[0][0], (1.3837) ** 2, places=2)

        try:
            Ag.epsilon(1 / 0.2)[0][0]
        except ValueError:
            pass
        else:
            raise AssertionError("Ag is not defined at a wavelength of 0.2 μm")

    def test_registry(self):
        # every entry of the registry is an mp.Medium, and lookups by name
        # return the same object exposed as a module attribute.
        self.assertGreater(len(materials.materials_library), 0)
        for name, medium in materials.materials_library.items():
            self.assertIsInstance(medium, mp.Medium)
            self.assertIs(getattr(materials, name), medium)
            self.assertIs(get_material(name), medium)

        # spot-check a few representative names.
        for name in ("Si", "Ag", "SiO2", "SiO2_aniso", "LiNbO3"):
            self.assertIn(name, materials.materials_library)


if __name__ == "__main__":
    unittest.main()
