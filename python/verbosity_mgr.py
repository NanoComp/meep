# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function


class Verbosity(object):
    """
    A class to help make accessing and setting the global verbosity level a bit
    more pythonic.

    Verbosity levels are:
        0: minimal output
        1: a little
        2: a lot (default)
        3: debugging
    """
    def __init__(self, cvar=None, initial_level=None):
        """
        Initialize the Verbosity manager. `cvar` should be some object that has
        a `verbosity` attribute, such as meep.cvar or mpb.cvar.
        """
        super(Verbosity, self).__init__()
        self.cvar = cvar
        if cvar is None or not hasattr(cvar, 'verbosity'):
            # If we're not given a module.cvar (e.g., while testing) or if the
            # cvar does not have a verbosity member (the lib hasn't been updated
            # yet) then use a dummy object so things can still run without it.
            class _dummy():
                verbosity = 1
            self.cvar = _dummy()
        if initial_level is not None:
            self.set(initial_level)

    def get(self):
        """
        Returns the current verbosity level.
        """
        return self.cvar.verbosity

    def set(self, level):
        """
        Validates the range, and sets the global verbosity level. Returns the
        former value.
        """
        if level < 0 or level > 3:
            raise ValueError('Only verbosity levels 0-3 are supported')
        old = self.cvar.verbosity
        self.cvar.verbosity = level
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
    verbosity = Verbosity()
    print('initial value: {}'.format(verbosity.get()))

    print('v==1: {}'.format(verbosity == 1))
    print('v==2: {}'.format(verbosity == 2))

    print('v>1: {}'.format(verbosity > 1))
    print('v<3: {}'.format(verbosity < 3))
    print('v>=2: {}'.format(verbosity >= 2))
    print('v<=1: {}'.format(verbosity <= 1))

    verbosity(1)
    print('v==1: {}'.format(verbosity == 1))
    print('v==2: {}'.format(verbosity == 2))

    # should raise ValueError
    verbosity(5)


if __name__ == '__main__':
    main()
