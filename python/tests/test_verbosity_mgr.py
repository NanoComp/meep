import unittest

from meep.verbosity_mgr import Verbosity


class VerbosityForTest(Verbosity):
    """Allows for testing of Verbosity without interfering with the singleton."""

    _instance = None


class MyCvar:
    def __init__(self):
        self.verbosity = 1


class TestVerbosity(unittest.TestCase):
    def setUp(self):
        VerbosityForTest.reset()
        self.v1 = VerbosityForTest(name="foo")
        self.v2 = VerbosityForTest(MyCvar(), "bar")

    def test_identity(self):
        # Ensure each verbosity is really the same singleton instance
        v1, v2 = self.v1, self.v2
        self.assertTrue(v1 is v2)
        self.assertEqual(id(v1), id(v2))
        self.assertEqual(v1.get_all(), [1, 1])

    def test_initial_value(self):
        v1, v2 = self.v1, self.v2
        self.assertEqual(v1.get(), 1)
        v2.set(2)
        self.assertEqual(v1.get(), 2)

    def test_properties(self):
        v1, v2 = self.v1, self.v2
        self.assertEqual(v1.foo, 1)
        self.assertEqual(v1.bar, 1)
        v1.foo = 2
        v2.bar = 3
        self.assertEqual(v2.foo, 2)
        self.assertEqual(v2.bar, 3)

    def test_operators(self):
        v1, v2 = self.v1, self.v2

        self.assertTrue(v1 == 1)
        self.assertFalse(v1 == 2)
        self.assertFalse(v1 > 1)
        self.assertTrue(v1 < 3)
        self.assertFalse(v1 >= 2)
        self.assertTrue(v1 <= 1)

        v1(3)
        self.assertFalse(v2 == 1)
        self.assertFalse(v2 == 2)
        self.assertTrue(v2 == 3)

    def test_out_of_range(self):
        v1, v2 = self.v1, self.v2

        with self.assertRaises(ValueError):
            v1.set(5)
        with self.assertRaises(ValueError):
            v1.set(-5)
        with self.assertRaises(ValueError):
            v2.foo = 5
        with self.assertRaises(ValueError):
            v2.bar = -5


if __name__ == "__main__":
    unittest.main()
