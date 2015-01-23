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
	*tikzDictionary

all:
	make synapter
	make clean

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi


synapter: 
	"$(R_HOME)/bin/R" CMD pdflatex synapter.tex
	bibtex synapter
	"$(R_HOME)/bin/R" CMD pdflatex synapter.tex
	"$(R_HOME)/bin/R" CMD pdflatex synapter.tex

## manual
synapter2: 
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('synapter.Rnw');"
	bibtex synapter
	"$(R_HOME)/bin/R" CMD pdflatex synapter.tex
	"$(R_HOME)/bin/R" CMD pdflatex synapter.tex


clean:
	rm -f $(LATEXFILES)
	rm -f *~
	rm -rf figure cache

pdfclean:
	make clean
	rm synapter.pdf
