LATEXFILES = *.aux\
        *.bbl\
        *.blg\
        *.ilg\
        *.log\
        *.nlo\
        *.nls\
        *.toc\
        *.aux\
        synapter.tex\
        fragmentmatching.tex\
	Rplots.pdf\
	*.dvi\
	*.map\
        *.out\
	*.figlist\
	*.dep\
	*.dpth\
	*.makefile\
	*-figure*.pdf\
	*tikzDictionary

all: synapter fragmentmatching clean

setvars:
ifeq (${R_HOME},)
R_HOME=$(shell R RHOME)
endif

synapter:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('synapter.Rnw');"
	"$(R_HOME)/bin/R" CMD pdflatex synapter.tex
	bibtex synapter
	"$(R_HOME)/bin/R" CMD pdflatex synapter.tex
	"$(R_HOME)/bin/R" CMD pdflatex synapter.tex

fragmentmatching:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('fragmentmatching.Rnw');"
	"$(R_HOME)/bin/R" CMD pdflatex fragmentmatching.tex
	bibtex fragmentmatching
	"$(R_HOME)/bin/R" CMD pdflatex fragmentmatching.tex
	"$(R_HOME)/bin/R" CMD pdflatex fragmentmatching.tex

clean:
	rm -f $(LATEXFILES)
	rm -f *~
	rm -f */*~
	rm -f .build.timestamp
	rm -f .Rhistory

clean2: clean
	rm -rf figure
	rm -rf cache

pdfclean:
	rm -f fragmentmatching.pdf
	rm -f synapter.pdf

all-clean: clean2 pdfclean
