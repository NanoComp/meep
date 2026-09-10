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

    def test_copy_returns_the_singleton(self):
        # Without __copy__/__deepcopy__, copying detaches the singleton from
        # its cvars: __new__ hands back the live instance and then copy
        # overwrites its __dict__ with copies of them.
        v1 = self.v1
        cvar = v1._cvars["bar"]
        self.assertTrue(copy.copy(v1) is v1)
        self.assertTrue(copy.deepcopy(v1) is v1)
        self.assertTrue(v1._cvars["bar"] is cvar)
        v1(0)
        self.assertEqual(cvar.verbosity, 0)

    def test_pickle_returns_the_singleton(self):
        # Uses the real Verbosity, since pickle needs an importable class.
        self.assertTrue(pickle.loads(pickle.dumps(Verbosity())) is Verbosity())

    def test_delete_flag_is_rejected(self):
        with self.assertRaises(AttributeError):
            del self.v1.foo

    def test_dir_includes_flags(self):
        listing = dir(self.v1)
        self.assertIn("foo", listing)
        self.assertIn("bar", listing)

    def test_hashable(self):
        # Defining __eq__ without __hash__ used to make instances unhashable.
        v1 = self.v1
        self.assertIsInstance(hash(v1), int)
        self.assertEqual(len({v1, self.v2}), 1)

    def test_no_args_does_not_add_a_flag(self):
        # A bare Verbosity() should just hand back the singleton, not register
        # another throwaway dummy flag.
        before = self.v1.get_all()
        v3 = VerbosityForTest()
        self.assertTrue(v3 is self.v1)
        self.assertEqual(v3.get_all(), before)

    def test_late_flag_adopts_current_level(self):
        # A flag registered after the level has been chosen should adopt it,
        # instead of keeping whatever defaults its C library was built with.
        v1 = self.v1
        v1(0)
        late = MyCvar(verbosity=2)
        v1.add_verbosity_var(late, "baz", 1)
        self.assertEqual(late.verbosity, 0)
        self.assertEqual(v1.baz, 0)
        self.assertEqual(v1.get(), 0)

    def test_first_flag_uses_initial_level(self):
        VerbosityForTest.reset()
        cvar = MyCvar(verbosity=3)
        v = VerbosityForTest(cvar, "solo", 2)
        self.assertEqual(cvar.verbosity, 2)
        self.assertEqual(v.get(), 2)

    def test_subclass_does_not_leak_flags_onto_base(self):
        # The flags registered on VerbosityForTest must not become attributes of
        # Verbosity itself, or this test file would corrupt meep.verbosity.
        self.assertNotIn("foo", Verbosity.__dict__)
        self.assertNotIn("bar", Verbosity.__dict__)

    def test_unknown_flag_raises_attribute_error(self):
        with self.assertRaises(AttributeError):
            self.v1.nonexistent

    def test_ordinary_attributes_still_work(self):
        self.v1.not_a_flag = 42
        self.assertEqual(self.v1.not_a_flag, 42)

    def test_reserved_names_are_rejected(self):
        with self.assertRaises(ValueError):
            self.v1.add_verbosity_var(MyCvar(), "set")
        with self.assertRaises(ValueError):
            self.v1.add_verbosity_var(MyCvar(), "_cvars")

    def test_non_integer_level(self):
        v1 = self.v1
        with self.assertRaises(TypeError):
            v1.set(1.5)
        with self.assertRaises(TypeError):
            v1.set("2")
        with self.assertRaises(TypeError):
            v1.foo = 2.5

    def test_verbosity_level_enum(self):
        v1 = self.v1
        v1(VerbosityLevel.DEBUG)
        self.assertEqual(v1.get(), 3)
        v1(VerbosityLevel.SILENT)
        self.assertEqual(v1.get(), 0)

    def test_temporary(self):
        v1 = self.v1
        v1(1)
        v1.bar = 2
        with v1.temporary(0) as v:
            self.assertTrue(v is v1)
            self.assertEqual(v1.get_all(), [0, 0])
        # Each flag goes back to its own former level, not just the global one.
        self.assertEqual(v1.get_all(), [1, 2])

    def test_temporary_with_flag_added_inside(self):
        # A flag registered inside the block should end up at the restored
        # global level, not stuck at the temporary one.
        v1 = self.v1
        v1(2)
        with v1.temporary(0):
            v1.add_verbosity_var(MyCvar(), "late")
            self.assertEqual(v1.late, 0)
        self.assertEqual(v1.late, 2)
        self.assertEqual(v1.get_all(), [2, 2, 2])

    def test_temporary_rejects_bad_level(self):
        v1 = self.v1
        with self.assertRaises(ValueError):
            with v1.temporary(9):
                pass
        self.assertEqual(v1.get_all(), [1, 1])

    def test_temporary_restores_on_exception(self):
        v1 = self.v1
        v1(1)
        with self.assertRaises(RuntimeError):
            with v1.temporary(3):
                raise RuntimeError("boom")
        self.assertEqual(v1.get_all(), [1, 1])

    def test_reset_drops_the_instance(self):
        v1 = self.v1
        VerbosityForTest.reset()
        v2 = VerbosityForTest(MyCvar(), "qux")
        self.assertFalse(v2 is v1)
        self.assertEqual(v2.get_all(), [1])
        # The dropped instance is left alone rather than being emptied out.
        self.assertEqual(v1.get_all(), [1, 1])


if __name__ == "__main__":
    unittest.main()
