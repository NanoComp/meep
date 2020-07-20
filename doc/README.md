This is the documentation tree for MEEP.

markdown (.md) files are in the `doc/docs` folder.

To build and visualize the HTML documentation locally using the mkdocs package
(useful for verifying changes on your local machine before committing), first
install `mkdocs` as well as two auxiliary packages via e.g.:

```
% pip3 install --user mkdocs python-markdown-math mkdocs-material
```

The main Python API document (`Python_User_Interface.md`) is generated from a
template file (named the same, with a `.in` extension), and the docstrings
located in the Python source code. Changes or additions for any text that is not
part of a docstring should be added or edited in the
`Python_User_Interface.md.in` file. Anything that is specifically about a
function, class or method should be documented in the docstring for that
function, class or method. The docstrings are expected to use MarkDown
formatting. To control where the docstrings are inserted into the documentation
a simple tagging system is used. See the documentation in the
`doc/generate_py_api.py` file for details.

If a docstring contains an "alternate function signature", usually to make
functions that accept a variety of positional or keyword arguments more clear to
the reader. Functions that fall into this category typically use Python's `*` or
`**` syntax. These lines can be tagged to indicate to the tool that they should
be moved or copied to the code block in the documentation where they function
signature is declared.  For example:

```python
    def add_flux(self, *args):
        """
        `add_flux(fcen, df, nfreq, freq, FluxRegions...)` ##sig

        ...
        """
```

Notice the use of  the `##sig` tag at the end of the line. If the alternate
signature is located in the middle of the docstring, and needs to remain there,
then the `##sig-keep` tag can be used instead. For an example of this usage, see
the `Simulation.run` method.

To (re)generate the Python API documentation (extracted from the docstrings)
run the following command in the project root folder. Note that this requires
that the Python extension for the meep library has been built and is importable:

```
% PYTHONPATH=./python python doc/generate_py_api.py
```

Note that this command should be rerun after making any changes to main template
file or the docstrings in the source, and rebuilding the project, in order to
update the documentation.

The view an auto-updating version of the documentation, run the following
command from the top-level MEEP repository tree:

```
% mkdocs serve
```

Finally, open the following address in a browser window:

http://127.0.0.1:8000

This launches a web server on your local machine plus a filesystem hook for
rebuilding the documentation tree automatically whenever any .md file is
modified. This enables viewing the HTML documentation in real time as the
source files are edited.

To build the static HTML version of the documentation, run:

```
% mkdocs build
```

The results will be put into the `./doc/site` folder.

Alternatively, the project `Makefile` can also be used to build the
documentation and to build a distribution archive containing the
documentation. To do so, run one of the following commands from the
top-level MEEP repository tree:

```
# Only regenerate the Python API document, useful if you're using `mkdocs serve`
% make python_api_doc

# Regenerate the API doc and run `mkdocs build`
% make docs

# Bundle the HTML version of the document files into a tarball archive
% make docs-dist
```
