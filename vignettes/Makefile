LATEXFILES = *.aux\
        *.bbl\
        *.blg\
        *.ilg\
        *.log\
        *.nlo\
        *.nls\
        *.toc\
        *.aux\
	.Rhistory\
	Rplots.pdf\
	synapter.tex\
	*.dvi\
	*.map\
        *.out\
	*.tikz\

all:
	make synapter
	make clean

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif

synapter: 
	## mkdir -p sweave-cache/figs/
	"$(R_HOME)/bin/R" --vanilla -e "library(pgfSweave); pgfSweave('synapter.Rnw');"
	bibtex synapter
	"$(R_HOME)/bin/R" CMD pdflatex synapter.tex
	"$(R_HOME)/bin/R" CMD pdflatex synapter.tex


clean:
	rm -f $(LATEXFILES)
	rm -f *~
	rm -rf sweave-cache

pdfclean:
	make clean
	rm synapter.pdf
