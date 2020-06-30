ACLOCAL_AMFLAGS=-I m4

SUBDIRS = src tests

if WITH_LIBCTL
if WITH_SCHEME
  SUBDIRS += scheme
endif
endif

if WITH_PYTHON
if WITH_MPB
if WITH_LIBCTLGEOM
  SUBDIRS += libpympb
endif
endif
  SUBDIRS += python
endif

EXTRA_DIST = NEWS.md LICENSE COPYRIGHT m4 meep-pkgconfig.in \
             doc generate_py_api.py mkdocs.yml
DIST_TARGETS = docs dist-gzip

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = meep.pc

meep.pc: meep-pkgconfig
	cp -f $(top_builddir)/meep-pkgconfig $@

libmeep:
	cd src && $(MAKE)

# Generate the Python API Document
# Set PYTHONPATH to use the in-place build of meep
python_api_doc: all
	test -d doc || cp -rf $(top_srcdir)/doc . && chmod -R +rw doc
	PYTHONPATH=./python $(PYTHON) generate_py_api.py

# Build the HTML version of the docs
docs: python_api_doc
	mkdir -p doc/site
	if test -z $$SKIP_BUILD_DOCS; then $(PYTHON) -m mkdocs build; else true; fi

# Build a tarball containing just the docs
docs-dist: docs
	tar czvf meep-documentation-$(PACKAGE_VERSION).tgz doc/site

DISTCLEANFILES = meep.pc
