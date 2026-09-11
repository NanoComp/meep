import enum
import numbers
from contextlib import contextmanager
from typing import Any, Dict, Iterator, List, Optional


class VerbosityLevel(enum.IntEnum):
    """
    Symbolic names for the verbosity levels understood by Meep. These are plain
    integers, so `mp.verbosity(mp.VerbosityLevel.SILENT)` and `mp.verbosity(0)`
    are equivalent.
    """

    SILENT = 0
    NORMAL = 1
    VERBOSE = 2
    DEBUG = 3


def _check_level(level: Any) -> int:
    """
    Validate a verbosity level and return it as a plain `int`. Raises `TypeError`
    if it is not an integer, or `ValueError` if it is out of range.
    """
    if not isinstance(level, numbers.Integral):
        raise TypeError(
            f"Verbosity level must be an integer, got {type(level).__name__}"
        )
    level = int(level)
    if level < VerbosityLevel.SILENT or level > VerbosityLevel.DEBUG:
        raise ValueError("Only verbosity levels 0-3 are supported")
    return level


class Verbosity:
    """
    A class to help make accessing and setting the global verbosity level a bit
    more Pythonic. It manages one or more verbosity flags that are located in
    the C/C++ libraries used by Meep.

    The verbosity levels are:

    * 0: minimal output (`VerbosityLevel.SILENT`)
    * 1: a little (`VerbosityLevel.NORMAL`)
    * 2: a lot (`VerbosityLevel.VERBOSE`)
    * 3: debugging (`VerbosityLevel.DEBUG`)

    An instance of `Verbosity` is created when meep is imported, and is
    accessible as `meep.verbosity`. The `meep.mpb` package also has a verbosity
    flag in its C library, and it can also be managed via the `Verbosity` class
    after `meep.mpb` is imported.

    Note that this class is a Singleton: every `Verbosity()` gives you the same
    instance, the one already available as `meep.verbosity`.

    The `Verbosity` instance can be used as a global verbosity controller, and
    assignments to any instance of `Verbosity` will set the global verbosity
    level for all library components. For example, this:

    ```python
    meep.verbosity(2)
    # or meep.verbosity.set(2) if you prefer being more explicit
    ```

    will set all of the managed verbosity flags to level 2.

    Each managed verbosity flag can also be accessed individually if desired,
    using the name it was registered under. Currently the names that are
    available are simply `meep` and `mpb`. This means that you can set two
    different verbosity levels like this:

    ```python
    verbosity = meep.verbosity # not required, it's just to save some typing
    verbosity.meep = 2
    verbosity.mpb = 1
    ```

    Note that while Meep is calling MPB internally the `mpb` flag is temporarily
    overridden to one level quieter than the `meep` flag (see the RAII class
    `meep::adjust_mpb_verbosity` in `src/adjust_verbosity.hpp`). That adjustment
    is undone as soon as the call returns, so it is not observable from Python;
    it only affects how chatty MPB is during Meep's own eigenmode calculations.
    """

    # Each verbosity flag is reached through a `cvar` object, which is what SWIG
    # calls the proxy it generates for a wrapped library's global variables.
    # Reading or assigning `cvar.verbosity` reads or writes the corresponding
    # C/C++ global directly, which is how this class controls the output of the
    # compiled libraries from Python:
    #
    # * `meep.cvar.verbosity` is the C++ global `meep::verbosity`, declared in
    #   `src/meep.hpp` and defined in `src/mympi.cpp`.
    # * `meep.mpb.cvar.verbosity` is MPB's C global `mpb_verbosity`, which is
    #    exposed under the shorter name by `%rename(verbosity) mpb_verbosity` in
    #   `python/mpb.i`.
    #
    # Those `cvar` objects are exactly what gets passed to `Verbosity()` or to
    # `add_verbosity_var()`.  That is, passing a `cvar` or a `name` to `Verbosity`
    # adds a new C `verbosity` flag to the list of flags managed by this class.
    # New flags can also be added explicitly with `add_verbosity_var()`.

    _instance: Optional["Verbosity"] = None

    def __new__(cls, *args, **kw) -> "Verbosity":
        # Create the real instance only the first time, and return the same each
        # time Verbosity is called thereafter.
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._init()
        return cls._instance

    def _init(self) -> None:
        """
        Set up the initial state of the singleton. Called only when the first
        instance is created.
        """
        # Note that these plain assignments are safe: __setattr__ falls through
        # to object.__setattr__ for anything that is not a registered flag name.
        self._cvars: Dict[str, Any] = {}
        # The name of the first flag registered, which is the one get() reports.
        self._primary: Optional[str] = None
        # Only used to answer get() before any flag has been registered.
        self._default_level = int(VerbosityLevel.NORMAL)

    @classmethod
    def reset(cls) -> None:
        # Probably just for testing. Drops the existing singleton instance so a
        # new one will be created the next time a new Verbosity is instantiated.
        # The dropped instance is deliberately left untouched, since other code
        # may still be holding a reference to it.
        cls._instance = None

    def __init__(
        self, cvar: Any = None, name: Optional[str] = None, initial_level: int = 1
    ) -> None:
        """See `add_verbosity_var()`"""
        # __init__ runs on every `Verbosity(...)` call, even though __new__ hands
        # back the existing singleton. Registering a flag only when the caller
        # actually asked for one keeps a bare `Verbosity()` from accumulating
        # throwaway dummy flags.
        if cvar is None and name is None:
            return
        self.add_verbosity_var(cvar, name, initial_level)

    def add_verbosity_var(
        self, cvar: Any = None, name: Optional[str] = None, initial_level: int = 1
    ) -> None:
        """
        Add a new verbosity flag to be managed. `cvar` is a SWIG `cvar` proxy
        for a wrapped library's global variables, such as `meep.cvar` or
        `meep.mpb.cvar` (see the class docstring for what those map onto in
        C/C++). Any object with a mutable `verbosity` attribute will do.

        The new flag is set to `initial_level` if it is the first one to be
        registered, and otherwise to the current global level, so that a flag
        registered late (for example by a deferred `import meep.mpb`) does not
        quietly ignore a verbosity level that has already been chosen.
        """
        if cvar is None or not hasattr(cvar, "verbosity"):
            # If we're not given a module.cvar (e.g., while testing) or if the
            # cvar does not have a verbosity member (e.g. the lib hasn't been
            # updated to have a verbosity flag yet) then use a dummy object so
            # things can still run without it.
            class _dummy:
                def __init__(self):
                    self.verbosity = 1

            cvar = _dummy()

        # If a name is not given then manufacture one
        if name is None:
            name = f"cvar_{len(self._cvars)}"
        # __setattr__ intercepts every assignment while __getattr__ only runs
        # when normal lookup fails, so a flag that shadows a class attribute
        # would be settable but not gettable. Leading underscores are reserved
        # for this class's own state.
        if not isinstance(name, str) or name.startswith("_"):
            raise ValueError(f"{name!r} is not a valid verbosity flag name")
        if name not in self._cvars and hasattr(type(self), name):
            raise ValueError(
                f"{name!r} can not be used as a verbosity flag name; it collides "
                f"with an existing attribute of {type(self).__name__}"
            )

        first = not self._cvars
        level = _check_level(initial_level) if first else self.get()
        self._cvars[name] = cvar
        if first:
            self._primary = name
        cvar.verbosity = level
        self._default_level = level

    def get(self) -> int:
        """
        Returns the current global verbosity level. This reads the first flag
        that was registered — the `meep` flag in a normal Meep session — so it
        reflects the live value of the C variable rather than a cached copy.
        """
        primary = self._primary
        if primary is None or primary not in self._cvars:
            return self._default_level
        return int(self._cvars[primary].verbosity)

    def get_all(self) -> List[int]:
        """
        Return a list of the values of all verbosity flags being managed. This
        is mostly intended for debugging this class and won't likely be useful
        otherwise.
        """
        return [int(cvar.verbosity) for cvar in self._cvars.values()]

    def set(self, level: int) -> int:
        """
        Validates the range, and sets the global verbosity level. Returns the
        former value.
        """
        level = _check_level(level)
        old = self.get()
        for cvar in self._cvars.values():
            cvar.verbosity = level
        self._default_level = level
        return old

    def __call__(self, level: int) -> int:
        """
        Convenience for setting the verbosity level. This lets you set the
        global level by calling the instance like a function. For example, if
        `verbosity` is an instance of this class, then its value can be changed
        like this:

        ```
        verbosity(0)
        ```
        """
        return self.set(level)

    @contextmanager
    def temporary(self, level: int) -> Iterator["Verbosity"]:
        """
        A context manager that sets the global verbosity level for the duration
        of a `with` block and then restores every managed flag to the level it
        had before, including when the block exists via an exception.

        ```python
        with meep.verbosity.temporary(0):
            sim.run(until=200)  # this part is quiet
        ```
        """
        saved = {name: int(cvar.verbosity) for name, cvar in self._cvars.items()}
        saved_default = self._default_level
        # set() validates before mutating anything, so a bad level raises here
        # and there is nothing to unwind.
        self.set(level)
        try:
            yield self
        finally:
            # Iterate the current flags rather than the snapshot, so a flag
            # registered inside the block (say by a deferred `import meep.mpb`)
            # is left at the restored global level instead of the temporary one.
            for name, cvar in self._cvars.items():
                cvar.verbosity = saved.get(name, saved_default)
            self._default_level = saved_default

    def __int__(self) -> int:
        """
        A convenience for getting the global verbosity level anywhere an integer
        is expected.
        """
        return self.get()

    def __index__(self) -> int:
        """
        Lets the global verbosity level be used directly anywhere Python wants
        an index, such as `range(verbosity)` or a numpy subscript.
        """
        return self.get()

    def __repr__(self) -> str:
        return f"Verbosity: level={self.get()}"

    # Verbosity is a mutable process-wide singleton, so a value-based hash would
    # not be stable. Using the identity hash keeps instances usable as dict keys
    # and set members, which defining __eq__ would otherwise take away. Note the
    # tradeoff: `verbosity == 1` is True but `hash(verbosity) != hash(1)`.
    __hash__ = object.__hash__

    # Copying has to hand back the same instance. Without these, `__new__`
    # returns the singleton and then `copy` overwrites its `__dict__` with
    # copies of the cvars, silently detaching it from the real C variables.
    def __copy__(self) -> "Verbosity":
        return self

    def __deep_copy__(self, memo: Dict[int, Any]) -> "Verboisty":
        return self

    def __reduce__(self):
        # Round-trips correctly only because Verbosity() with no arguments
        # registers nothing and just returns the existing singleton.
        return (type(self), ())

    def _comparable(self, other: Any) -> Optional[Any]:
        """
        Return `other` as something the global level can be compared against, or
        `None` if no meaningful comparison exists. `numbers.Real` covers `int`,
        `bool`, `float`, and the NumPy scalar types.
        """
        if isinstance(other, Verbosity):
            return other.get()
        if isinstance(other, numbers.Real):
            return other
        return None

    # Some comparison operators. These return NotImplemented rather than raising
    # for operands that are not numbers, so that Python can fall back to its own
    # handling (a clear TypeError for ordering, identity-based ==/!=).
    def __gt__(self, o: Any) -> Any:
        other = self._comparable(o)
        return NotImplemented if other is None else self.get() > other

    def __lt__(self, o: Any) -> Any:
        other = self._comparable(o)
        return NotImplemented if other is None else self.get() < other

    def __eq__(self, o: Any) -> Any:
        other = self._comparable(o)
        return NotImplemented if other is None else self.get() == other

    def __le__(self, o: Any) -> Any:
        other = self._comparable(o)
        return NotImplemented if other is None else self.get() <= other

    def __ge__(self, o: Any) -> Any:
        other = self._comparable(o)
        return NotImplemented if other is None else self.get() >= other

    def __getattr__(self, name: str) -> int:
        """
        Look up a managed verbosity flag by the name it was registered under, so
        that `verbosity.meep` returns the current level of the `meep` flag. Only
        consulted when normal attribute lookup fails.
        """
        if name.startswith("_"):
            # Covers this class's own state before _init() has run, plus every
            # dunder probe (__deepcopy__, __setstate__, _repr_html_, ...) that
            # would otherwise have to fall through the lookup below.
            raise AttributeError(
                f"{type(self).__name__!r} object has no attribute {name!r}"
            )
        cvars = self.__dict__.get("_cvars")
        if cvars is not None and name in cvars:
            return int(cvars[name].verbosity)
        raise AttributeError(
            f"{type(self).__name__!r} object has no verbosity flag {name!r}"
        )

    def __setattr__(self, name: str, value: Any) -> None:
        """
        Assign to a managed verbosity flag by name, validating the level, and
        otherwise behave like an ordinary attribute assignment.
        """
        cvars = self.__dict__.get("_cvars")
        if cvars is not None and name in cvars:
            cvars[name].verbosity = _check_level(value)
        else:
            object.__setattr__(self, name, value)

    def __delattr__(self, name: str) -> None:
        """
        Managed verbosity flags are removed by `reset()`, not by `del`.
        """
        cvars = self.__dict__.get("_cvars")
        if cvars is not None and name in cvars:
            raise AttributeError(f"cannot delete managed verbosity flag {name!r}")
        object.__delattr__(self, name)

    def __dir__(self) -> List[str]:
        return sorted(set(super().__dir__()) | set(self.__dict__.get("_cvars", {})))
