### Introduction

The `synapter` package provides functionality to re-analyse MSe
label-free proteomics data acquired on a Waters Synapt Series mass
spectrometer (and probably any Waters instrument). It allows to
combine acquisitions that have been optimised for better
identification (typically using ion mobility separation - HDMSe) and
quantitation accuracy. It also allows to transfer identifications
across multiple runs to reduce missing data across an experiment.

The official release is the Bioconductor version, available
[here](http://bioconductor.org/packages/devel/bioc/html/synapter.html). The
[github page](https://lgatto.github.io/synapter/) is a useful resource
that gives access to all vignettes and manuals.

### Installation

`synapter` is available from the
[Bioconductor](http://www.bioconductor.org) repository. The package
and its dependencies can be installed with

     source("http://www.bioconductor.org/biocLite.R")
     biocLite("synapter")

### Help

`synapter` comes with plenty of documentation. Have a start with the
package documentation page `?synapter` and the vignette

     vignette("synapter", package="synapter")


See also the `synapter` [Bioconductor
page](http://bioconductor.org/packages/devel/bioc/html/synapter.html)
for on-line access to the vignette and the reference manual.

### GitHub build status

Current: [![Build
Status](https://travis-ci.org/lgatto/synapter.svg?branch=master)](https://travis-ci.org/lgatto/synapter)

### PLGS processing

The raw data files produced must first be processed by Water's PLGS
software to produce `synapter` input files. This is described in
details in the vignette. Additional information with lots of
screenshots can be found in [these
slides](http://proteome.sysbiol.cam.ac.uk/lgatto/synapter/PLGS_Data_Processing.pdf).
