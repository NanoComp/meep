# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function

from ._mpb import cvar

class Verbosity(int):
    """
    A class to help make accessing and setting the global verbosity level a bit
    more pythonic.

    Verbosity levels are:
        0: minimal output
        1: a little
        2: a lot (default)
        3: debugging
    """

    def get(self):
        """
        Returns the current verbosity level.
        """
        return cvar.verbosity

    def set(self, value):
        """
        Validates the range, and sets the global verbosity level.
        """
        if value < 0 or value > 3:
            raise ValueError('Only verbosity levels 0-3 are supported')
        cvar.verbosity = value

    def __call__(self, value):
        """
        Convenience for setting the verbosity level. This lets you set the level
        by calling the instance like a function. For example, if `verbosity` is
        an instance of this class, then it's value can be changed like this:

        ```
        verbosity(0)
        ```
        """
        self.set(value)

    def __int__(self):
        """
        A convenience for getting the verbosity level anywhere an integer is
        expected.
        """
        #print('__int__')
        return self.get()

    # Some comparrison operators
    def __gt__(self, o): return self.get() > o
    def __lt__(self, o): return self.get() < o
    def __eq__(self, o): return self.get() == o
    def __ne__(self, o): return self.get() != o
    def __le__(self, o): return self.get() <= o
    def __ge__(self, o): return self.get() >= o



def main():
    # simple test code
    verbosity = Verbosity()
    print(f'initial value: {verbosity.get()}')

    print(f'v==1: {verbosity == 1}')
    print(f'v==2: {verbosity == 2}')

    print(f'v>1: {verbosity > 1}')
    print(f'v<3: {verbosity < 3}')
    print(f'v>=2: {verbosity >= 2}')
    print(f'v<=1: {verbosity <= 1}')

    verbosity(1)
    print(f'v==1: {verbosity == 1}')
    print(f'v==2: {verbosity == 2}')

    # should raise ValueError
    verbosity(5)


if __name__ == '__main__':
    main()
