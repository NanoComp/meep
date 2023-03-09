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

EXTRA_DIST = NEWS.md LICENSE COPYRIGHT m4 meep-pkgconfig.in version.sh

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = meep.pc

meep.pc: meep-pkgconfig
	cp -f $(top_builddir)/meep-pkgconfig $@

libmeep:
	cd src && $(MAKE)

# Generate the Python API Document
# PYTHONPATH is set to use the in-place build of meep
python_api_doc: all
	if test '$(top_srcdir)' != . ; then \
		cp -rf $(top_srcdir)/doc .; \
		chmod -R +rw doc; \
	fi
	if test -z $$SKIP_BUILD_DOCS; then \
		PYTHONPATH=./python $(PYTHON) doc/generate_py_api.py; \
	fi

# Build the HTML version of the docs
docs: python_api_doc
	if test -z $$SKIP_BUILD_DOCS; then \
		$(PYTHON) -m mkdocs build -f $(top_srcdir)/mkdocs.yml; \
	fi

# Build a tarball containing just the html docs
docs-dist: docs
	tar czvf meep-documentation-$(PACKAGE_VERSION).tgz doc/site

clean-local:
	rm -f meep-documentation-$(PACKAGE_VERSION).tgz
	if test '$(top_srcdir)' != . ; then \
		rm -rf doc; \
	fi

DISTCLEANFILES = meep.pc
