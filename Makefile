STEM = qtl_in_mpp

FIGS = Figs/intercross.pdf \
	   Figs/lodcurve_insulin.pdf \
	   Figs/data_fig.png \
	   Figs/ail.pdf \
	   Figs/risib.pdf \
	   Figs/riself.pdf \
	   Figs/hs.pdf \
	   Figs/ri8.pdf \
	   Figs/ri8self.pdf \
	   Figs/do_genome.pdf \
	   Figs/congenic.pdf \
	   Figs/permtest.pdf
#	   Figs/genome_reconstr.pdf \
#	   Figs/qtl_scan.png \
#	   Figs/hmm.pdf \
#	   Figs/founder_pop.pdf \
#	   Figs/cc_xchr_reconstr.pdf \
#	   Figs/do_qtl.pdf

R_OPTS=--no-save --no-restore --no-init-file --no-site-file

all: docs/$(STEM).pdf

docs/%.pdf: %.pdf
	cp $^ $@

$(STEM).pdf: $(STEM).tex header.tex $(FIGS)
	xelatex $^

Figs/%.pdf: R/%.R
	cd R;R CMD BATCH $(R_OPTS) $(<F)

Figs/%.png: R/%.R
	cd R;R CMD BATCH $(R_OPTS) $(<F)

clean:
	rm *.aux *.log *.nav *.out *.snm *.toc *.vrb
