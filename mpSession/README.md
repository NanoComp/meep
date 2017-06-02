Second attempt at implementing `mpSession.`

-- `libctlgeom-noscheme` is a scheme-independent refactoring of libctlgeom.

-- General idea: `mpSession` will provide all functionality currently
   available in the libctl-enabled meep binary, so that any meep
   calculation currently driven by a `.ctl` file can be ported
   basically line-by-line to a `.cpp` program.

-- As examples,
   `pw-source.cpp` and `holey-wvg-bands.cpp` are C++ ports of
   `pw-source.ctl` and `holey-wvg-bands.ctl` in `meep/examples.`

-- Then the python module can live atop `mpSession,` making use
   of its C++ implementations of low-level functionality for
   e.g. geometry operations, while potentially replacing other
   functionality (callbacks, etc.) with pythonic alternatives.
