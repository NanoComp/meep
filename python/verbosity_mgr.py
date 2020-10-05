# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function


class Verbosity(object):
    """
    A class to help make accessing and setting the global verbosity level a bit
    more pythonic. It manages one or more verbosity flags that are located in
    the C/C++ libraries used by Meep.

    The verbosity levels are:

    * 0: minimal output
    * 1: a little
    * 2: a lot
    * 3: debugging

    An instace of `Verbosity` is created when meep is imported, and is accessible
    as `meep.verbosity`.

    Note that this class is a Singleton, meaning that each call to
    `Verbosity(cvar)` gives you the same instance. The new `cvar` will be added
    to a list of verbosity flags managed by this class.
    """
    _instance = None

    def __new__(cls, *args, **kw):
        # Create the real instance only the first time, and return the same each
        # time Verbosity is called thereafter.
        if cls._instance is None:
            cls._instance = super(Verbosity, cls).__new__(cls)
            cls._instance._init()
        return cls._instance

    def _init(self):
        """
        Set up the initial state of the singleton. Called only when the first
        instance is created.
        """
        self._master_verbosity = -1
        self._cvars = dict()

    @classmethod
    def reset(cls):
        # Probably just for testing. Drops the existing the singleton instance
        # so a new one will be created the next time a new Verbosity is
        # instantiated.
        if cls._instance is not None:
            cls._instance._init()
            cls._instance = None

    def __init__(self, cvar=None, name=None, initial_level=1):
        """See `add_verbosity_var()`"""
        self.add_verbosity_var(cvar, name, initial_level)

    def add_verbosity_var(self, cvar=None, name=None, initial_level=1):
        """
        Add a new verbosity flag to be managed. `cvar` should be some object
        that has a `verbosity` attribute, such as `meep.cvar` or `mpb.cvar`.
        """
        if cvar is None or not hasattr(cvar, 'verbosity'):
            # If we're not given a module.cvar (e.g., while testing) or if the
            # cvar does not have a verbosity member (e.g. the lib hasn't been
            # updated to have a verbosity flag yet) then use a dummy object so
            # things can still run without it.
            class _dummy():
                def __init__(self): self.verbosity = 1
            cvar = _dummy()

        # If a name is not given then manufacture one
        if name is None:
            name = 'cvar_{}'.format(len(self._cvars))
        self._cvars[name] = cvar

        # And create a property for so it can be accessed like `verbosity.mpb`
        self.make_property(name)

        # The initial master (global) level is determined by the first instance created
        if self._master_verbosity == -1:
            self.set(initial_level)

    def get(self):
        """
        Returns the current global verbosity level.
        """
        return self._master_verbosity

    def get_all(self):
        """
        Return a list of the values of all verbosity flags being managed. This
        is mostly intended for debugging this class and won't likely be useful
        otherwise.
        """
        return [value.verbosity for value in self._cvars.values() ]

    def set(self, level):
        """
        Validates the range, and sets the global verbosity level. Returns the
        former value.
        """
        if level < 0 or level > 3:
            raise ValueError('Only verbosity levels 0-3 are supported')
        old = self._master_verbosity
        for cvar in self._cvars.values():
            cvar.verbosity = level
        self._master_verbosity = level
        return old

    def __call__(self, level):
        """
        Convenience for setting the verbosity level. This lets you set the
        global level by calling the instance like a function. For example, if
        `verbosity` is an instance of this class, then it's value can be changed
        like this:

        ```
        verbosity(0)
        ```
        """
        return self.set(level)

    def __int__(self):
        """
        A convenience for getting the global verbosity level anywhere an integer
        is expected.
        """
        return self.get()

    def __repr__(self):
        return "Verbosity: level={}".format(self.get())

    # Some comparison operators
    def __gt__(self, o): return self.get() > o
    def __lt__(self, o): return self.get() < o
    def __eq__(self, o): return self.get() == o
    def __ne__(self, o): return self.get() != o
    def __le__(self, o): return self.get() <= o
    def __ge__(self, o): return self.get() >= o

    def make_property(self, name):
        """
        Add a property to the class with the given name that gets or sets the
        verbosity in the cvar with that name in self._cvars.
        """
        def _getter(self, name=name):
            return self._cvars[name].verbosity
        def _setter(self, level, name=name):
            if level < 0 or level > 3:
                raise ValueError('Only verbosity levels 0-3 are supported')
            self._cvars[name].verbosity = level

        setattr(Verbosity, name, property(_getter, _setter))

