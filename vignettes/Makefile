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
        ionmobility.tex\
        crossmatching.tex\
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

all:
	make synapter
	make crossmatching
	make ionmobility
	make clean

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif

synapter: 
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('synapter.Rnw');"
	"$(R_HOME)/bin/R" CMD pdflatex synapter.tex
	bibtex synapter
	"$(R_HOME)/bin/R" CMD pdflatex synapter.tex
	"$(R_HOME)/bin/R" CMD pdflatex synapter.tex

ionmobility: 
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('ionmobility.Rnw');"
	"$(R_HOME)/bin/R" CMD pdflatex ionmobility.tex
	bibtex ionmobility
	"$(R_HOME)/bin/R" CMD pdflatex ionmobility.tex
	"$(R_HOME)/bin/R" CMD pdflatex ionmobility.tex

crossmatching: 
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('crossmatching.Rnw');"
	"$(R_HOME)/bin/R" CMD pdflatex crossmatching.tex
	bibtex crossmatching
	"$(R_HOME)/bin/R" CMD pdflatex crossmatching.tex
	"$(R_HOME)/bin/R" CMD pdflatex crossmatching.tex

clean:
	rm -f $(LATEXFILES)
	rm -f *~
	rm -f */*~
	rm -f .build.timestamp
	rm -f .Rhistory

clean2: 
	make clean
	rm -rf figure
	rm -rf cache

pdfclean: 
	rm -f crossmatching.pdf
	rm -f synapter.pdf
	rm -f ionmobility.pdf

all-clean:
	make clean2
	make pdfclean
