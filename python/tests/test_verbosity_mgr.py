import copy
import pickle
import unittest

from meep.verbosity_mgr import Verbosity, VerbosityLevel


class VerbosityForTest(Verbosity):
    """Allows for testing of Verbosity without interfering with the singleton."""

    _instance = None


class MyCvar:
    def __init__(self, verbosity=1):
        self.verbosity = verbosity


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

        # The dedicated helpers invoke the operator they are named for, so they
        # stand in wherever the expected result is True. Where it is False they
        # cannot: assertLessEqual(v1, 1) would exercise __le__ rather than a
        # false __gt__. Those keep the explicit form, with a message that
        # reports the level the way the helpers do.
        self.assertEqual(v1, 1)
        self.assertFalse(v1 == 2, f"{v1!r} should not equal 2")
        self.assertFalse(v1 > 1, f"f{v1!r} should not be greater than 1")
        self.assertLess(v1, 3)
        self.assertFalse(v1 >= 2, f"{v1!r} should not be >= 2")
        self.assertTrue(v1 <= 1)

        v1(3)
        self.assertFalse(v2 == 1, f"{v2!r} should not equal 1")
        self.assertFalse(v2 == 2, f"{v2!r} should not equal 2")
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

    def test_not_equal(self):
        # __ne__ is not defined; Python derives it from __eq__. Both operators
        # are probed in both directions, so assertNotEqual can only stand in
        # where the operator under test really is `!=`.
        v1 = self.v1
        self.assertNotEqual(v1, 2)
        self.assertFalse(v1 != 1, f"{v1!r} != 1 should be False")
        # Falls back to identity for things that aren't numbers, rather than
        # raising from inside the comparison.
        self.assertFalse(v1 == "one", f"{v1!r} == 'one' should be False")
        self.assertTrue(v1 != "one")

    def test_ordering_with_non_number(self):
        with self.assertRaises(TypeError):
            self.v1 > "one"

    def test_int_and_repr(self):
        v1 = self.v1
        self.assertEqual(int(v1), 1)
        self.assertEqual(repr(v1), "Verbosity: level=1")
        v1(2)
        self.assertEqual(int(v1), 2)
        self.assertEqual(repr(v1), "Verbosity: level=2")

    def test_set_returns_former_value(self):
        v1 = self.v1
        self.assertEqual(v1.set(3), 1)
        self.assertEqual(v1.set(0), 3)
        self.assertEqual(v1(2), 0)

    def test_get_is_not_stale_after_individual_set(self):
        # Setting a single flag used to leave get() reporting the old level.
        v1 = self.v1
        v1.foo = 3
        self.assertEqual(v1.get(), 3)
        self.assertEqual(int(v1), 3)
        self.assertEqual(v1, 3)
        self.assertEqual(repr(v1), "Verbosity: level=3")

    def test_get_reads_the_first_registered_flag(self):
        # 'foo' was registered first, so it is the one get() reports.
        v1 = self.v1
        v1.bar = 3
        self.assertEqual(v1.foo, 1)
        self.assertEqual(v1.get(), 1)

    def test_index(self):
        self.v1(2)
        self.assertEqual(list(range(self.v1)), [0, 1])

    def test_compares_against_floats_and_bools(self):
        v1 = self.v1
        self.assertGreater(v1, 0.5)
        self.assertLess(v1, 1.5)
        self.assertEqual(v1, True)
        self.assertEqual(v1, v1)


if __name__ == "__main__":
    unittest.main()
