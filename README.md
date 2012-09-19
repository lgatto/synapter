### Introduction

The `synapter` package provides functionality to re-analyse MSe label-free proteomics data acquired on a Waters Synapt Series mass spectrometer (and probably any Waters instrument). It allows to combine acquisitions that have been optimised for better identification (typically using ion mobility separation - HDMSe) and quantitation accuracy. It also allows to transfer identifications across multiple runs to reduce missing data across an experiment. 

The analysis pipeline can be executed using a simple graphical user interface (started with `synapterGUI()`) or a high-level function `synergise` to produce a [html report](http://proteome.sysbiol.cam.ac.uk/lgatto/synapter/Report/index.html). Alternatively, or low-level interface is available (see `?Synapter`).

### Help

`synapter` comes with plenty of documentation. Have a start with the package documentation page `?synapter` and the vignette 

     vignette("synapter", package="synapter")

Do not hesitate to contact [me](http://proteome.sysbiol.cam.ac.uk/lgatto/) for questions/comments/suggestions.

### PLGS processing

The raw data files produced must first be processed by Water's PLGS software to produce `synapter` input files. This is described in details in the vignette. Additional information with lots of screenshots can be found in [these slides](http://proteome.sysbiol.cam.ac.uk/lgatto/synapter/PLGS_Data_Processing.pdf).

### Installation
`synapter` is available from the [Bioconductor](http://www.bioconductor.org) repository. The package and its dependencies can be installed

     source("http://www.bioconductor.org/biocLite.R")
     biocLite("synapter")

See also the `synapter` [Bioconductor page](http://bioconductor.org/packages/devel/bioc/html/synapter.html) for on-line access to the vignette and the reference manual.

**NB:** Installation of the package requires `R` version 2.15.1 (the latest stable release). It is currently in the Bioconductor development branch (next release is scheduled for the 2nd of October 2012). To do so, on a clean `R` installation, type the following:

     source("http://www.bioconductor.org/biocLite.R") ## to get up-to-speed with Bioc
     useDevel(TRUE) ## to use devel branch
     biocLite("synapter") ## to install synapter

Next time a package needs to be installed, first `load("BiocInstaller")` and then install the package with `biocLite("packageName")`.

### Source code

The code on [github](https://github.com/lgatto/synapter/) is for sharing, testing, issue tracking and forking/pulling purposes. Although it should be in sync with the code on the [Bioconductor svn server](https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/), the latter is the official repository for the working source code. Get is with

    svn co https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/rols

(user name: `readonly`, password: `readonly`)
