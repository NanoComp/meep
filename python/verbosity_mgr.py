# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function


class Verbosity(object):
    """
    A class to help make accessing and setting the global verbosity level a bit
    more pythonic. It manages one or more verbosity flags that are located in
    the C/C++ libraries used by Meep.

    Verbosity levels are:
        0: minimal output
        1: a little
        2: a lot (default)
        3: debugging

    Note that this class is a Singleton, meaning that each call to
    Verbosity(cvar) gives you the same instance. The new cvar will be added to a
    list of verbosity flags managed by this class.
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
        self._master_verbosity = 1
        self._cvars = list()

    def __init__(self, cvar=None, initial_level=None):
        """See add_verbosity_var()"""
        self.add_verbosity_var(cvar, initial_level)

    def add_verbosity_var(self, cvar=None, initial_level=None):
        """
        Add a new verbosity flag to be managed. `cvar` should be some object
        that has a `verbosity` attribute, such as meep.cvar or mpb.cvar.
        """
        if cvar is None or not hasattr(cvar, 'verbosity'):
            # If we're not given a module.cvar (e.g., while testing) or if the
            # cvar does not have a verbosity member (e.g. the lib hasn't been
            # updated to have a verbosity flag yet) then use a dummy object so
            # things can still run without it.
            class _dummy():
                def __init__(self): self.verbosity = 1
            cvar = _dummy()
        self._cvars.append(cvar)
        if initial_level is not None:
            self.set(initial_level)

    def get(self):
        """
        Returns the current verbosity level.
        """
        return self._master_verbosity

    def get_all(self):
        return [cvar.verbosity for cvar in self._cvars ]

    def set(self, level):
        """
        Validates the range, and sets the global verbosity level. Returns the
        former value.
        """
        if level < 0 or level > 3:
            raise ValueError('Only verbosity levels 0-3 are supported')
        old = self._master_verbosity
        for cvar in self._cvars:
            cvar.verbosity = level
        self._master_verbosity = level
        return old

    def __call__(self, level):
        """
        Convenience for setting the verbosity level. This lets you set the level
        by calling the instance like a function. For example, if `verbosity` is
        an instance of this class, then it's value can be changed like this:

        ```
        verbosity(0)
        ```
        """
        return self.set(level)

    def __int__(self):
        """
        A convenience for getting the verbosity level anywhere an integer is
        expected.
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




def main():
    # simple test code
    import sys

    class MyCvar:
        def __init__(self): self.verbosity = 1

    verbosity = Verbosity()
    v2 = Verbosity(MyCvar())
    print(verbosity is v2)
    print(id(verbosity), id(v2))
    print(verbosity.get_all())

    print('initial value: {}'.format(verbosity.get()))

    print('v==1: {}'.format(verbosity == 1))
    print('v==2: {}'.format(verbosity == 2))

    print('v>1: {}'.format(verbosity > 1))
    print('v<3: {}'.format(verbosity < 3))
    print('v>=2: {}'.format(verbosity >= 2))
    print('v<=1: {}'.format(verbosity <= 1))

    verbosity(3)
    print('v==1: {}'.format(verbosity == 1))
    print('v==2: {}'.format(verbosity == 2))
    print('v==3: {}'.format(verbosity == 3))

    print(verbosity.get_all())

    # should raise ValueError
    verbosity(5)


if __name__ == '__main__':
    main()
