noinst_PROGRAMS = example

TEX = complicated.tex energy_cons.tex energy_cons_1d.tex epsilon_polariton_1d.tex gpl.tex lossgain_epsilon.tex meep.tex omniguide.tex polaritonbands.tex simple.tex simplebands.tex

EPS = polaritonbands-out/bands.eps epsilon_polariton_1d-out/eps.eps omniguide-out/bands.eps lossgain_epsilon-out/eps.eps energy_cons_1d-out/energy.eps energy_cons-out/energy.eps energy_cons-out/energy.eps simple-out/ez-000200-00.eps Yee_bulk.eps

PDF = $(EPS:.eps=.pdf)

FREQS = complicated-out/freqs

EXTRA_DIST = html meep.pdf $(TEX) Yee_bulk.eps preproc.pl

AM_CPPFLAGS = -I$(top_srcdir)/src

LIBMEEP = $(top_builddir)/src/libmeep.la

example_SOURCES = example.cpp
example_LDADD = $(LIBMEEP)

###########################################################################

SUFFIXES = .dac .tex .dvi .ps .pdf .eps

.o.dac: $(LIBMEEP)
	$(CXXLINK) $< $(LIBMEEP) $(LIBS)

.tex.cpp:
	perl ./preproc.pl $< > $@

.eps.pdf:
	epstopdf --outfile=$@ $<

html: $(TEX) $(EPS) $(FREQS)
	rm -rf $@
	mkdir $@
	@LATEX2HTML@ -split +1 -local_icons -dir $@ meep

meep.pdf: $(TEX) $(PDF) $(FREQS)
	pdflatex meep

meep.ps: meep.dvi
	dvips -o meep.ps meep.dvi

meep.dvi: $(TEX) $(EPS) $(FREQS)
	latex meep.tex

###########################################################################

simple-out/ez-000200-00.eps: simple.dac
	rm -rf simple-out && $(RUNCODE) ./simple.dac
	mv simple-out/ez-000200.00.eps $@

complicated-out/freqs: complicated.dac
	rm -rf complicated-out/ && $(RUNCODE) ./complicated.dac

energy_cons_1d-out/energy.eps : energy_cons_1d.dac
	rm -rf energy_cons_1d-out/ && $(RUNCODE) ./energy_cons_1d.dac
	mv $@ $@-tmp
	sed -e 's/fill//' $@-tmp > $@
energy_cons-out/energy.eps : energy_cons.dac
	rm -rf energy_cons-out/ && $(RUNCODE) ./energy_cons.dac
	mv $@ $@-tmp
	sed -e 's/fill//' $@-tmp > $@
polaritonbands-out/bands.eps : polaritonbands.dac
	rm -rf polaritonbands-out/ && $(RUNCODE) ./polaritonbands.dac
	mv $@ $@-tmp
	sed -e 's/fill//' $@-tmp > $@

epsilon_polariton_1d-out/eps.eps : epsilon_polariton_1d.dac
	rm -rf epsilon_polariton_1d-out/
	$(RUNCODE) ./epsilon_polariton_1d.dac
	mv $@ $@-tmp
	sed -e 's/fill//' $@-tmp > $@

omniguide-out/bands.eps : omniguide.dac
	rm -rf omniguide-out/
	$(RUNCODE) ./omniguide.dac
	mv $@ $@-tmp
	sed -e 's/fill//' $@-tmp > $@

lossgain_epsilon-out/eps.eps : lossgain_epsilon.dac
	rm -rf lossgain_epsilon-out/
	$(RUNCODE) ./lossgain_epsilon.dac
	mv $@ $@-tmp
	sed -e 's/fill//' $@-tmp > $@

check-local: simplebands.dac complicated.dac polaritonbands.dac \
	epsilon_polariton_1d.dac omniguide.dac lossgain_epsilon.dac \
	energy_cons_1d.dac energy_cons.dac simple.dac


clean-local:
	rm -rf *.o *.hi *.dac *-out *-out-* *.ps *.dvi *.doc *.log *.hi *.dac $(PDF) $(EPS)

maintainer-clean-local:
	rm -rf html *.pdf

