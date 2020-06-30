This is the documentation tree for MEEP.

markdown (.md) files are in the `doc/docs` folder.

To build and visualize the HTML documentation locally using the
mkdocs package (useful for verifying changes on your local machine
before committing), first install `mkdocs` (version 0.17.5) as
well as two auxiliary packages via e.g.:

```
% pip3 install --user mkdocs==0.17.5 python-markdown-math mkdocs-material
```

To (re)generate the Python API documentation (extracted from the docstrings)
run the following command in the project root folder. Note that this requires
that the Python extension for the meep library has been built:

```
% PYTHONPATH=./python python generate_py_api.py
```

Rerun this after making any changes to the docstrings in the source and
rebuilding the project, in order to update the documentation.

Next, run the following command from the top-level MEEP repository tree:

```
% mkdocs serve
```

Finally, open the following address in a browser window:

http://127.0.0.1:8000

This launches a web server on your local machine plus a filesystem hook for
rebuilding the documentation tree automatically whenever any .md file is
modified. This enables viewing the HTML documentation in real time as the
source files are edited.

To build the HTML version of the documentation, run:

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

# Bundle the document files into a tarball archive
% make docs-dist
```
