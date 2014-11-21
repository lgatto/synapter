## useful refs for tcltk/bwidgets
## http://www.activestate.com/activetcl
## http://www2.warwick.ac.uk/fac/sci/systemsbiology/facilities/computing/support/limmagui/
## https://stat.ethz.ch/pipermail/bioconductor/2007-April/016878.html
## https://stat.ethz.ch/pipermail/r-sig-mac/2006-October/003301.html
## http://stackoverflow.com/questions/6234690/resizeable-tk-windows-done-right

##' Note: The tcltk interface is being deprecated. Please consider
##' using its underlying \code{synergise()} function instead. 
##' 
##' This function starts the synapter graphical user interface
##' that allows to (1) load \code{n} sets of files to be analyses,
##' (2) parameters used for data filtering and retention time
##' modelling and (3) parameters used for grid search.
##' The \code{n} file sets and the parameters recorded by the GUI
##' are used to call \code{\link{synergise}} that performs the
##' analysis and generates \code{n} html reports.
##' See the \code{synapter} vignette for a description of the
##' interface and \code{\link{synergize}} for a description of the
##' parameters and data processing.
##'
##' The GUI code requires the Tcl package 'BWidget' to display the tree
##' drill-down file selection menu. The package is available on Windows
##' (in case of issues, please read
##' \url{http://www.stat.berkeley.edu/users/spector/s133/Bwidget.html}).
##' On GNU/Linux, the package needs to be installed separately;
##' on Ubuntu, the package is called 'bwidget'.
##'
##' @title Launch the synapter GUI
##' @param n Number of analysis to be run.
##' @return Invisibly returns \code{null}. Used for its side effects.
##' @author Laurent Gatto
synapterGUI <- function(n = 1) {
    .Deprecated("synergise", package = "synapter",
                msg = c("The synapterGUI() interface is deprecated.\n",
                    "Please consider using synergise() instead."))
    if (!require("tcltk") | !require("tcltk2"))
        stop("The gui requires the 'tcltk' and 'txltk2' packages.")
    if (n < 1)
        stop("You must provide at least one set of files.")
    gui_output <- NULL
    .GUI(n, e = environment())
    if (!is.null(gui_output)) {
        for (i in 1:n) {
            sel <- grep(paste0(i, "$"), names(gui_output$filenames))    
            .filenames <- gui_output$filenames[sel]
            out <- grep(paste0("Out", i), names(.filenames))
            .outputdir <- .filenames[out]
            .filenames <- c(.filenames[-out],
                            gui_output$filenames["fasta"])
            .filenames <- as.list(.filenames)
            names(.filenames) <- c("quantpep3d", "identpeptide",                           
                                   "quantpeptide", "fasta")
            message("(", i, ") Processing ", basename(.filenames[[2]]))
            synergise(filenames = .filenames,
                      master = gui_output$master,
                      outputdir = .outputdir,
                      fdr = gui_output$pepfdr,
                      fpr = gui_output$prtfpr,
                      identppm = gui_output$identppm,
                      quantppm = gui_output$quantppm,
                      uniquepep = gui_output$uniquepep,
                      span = gui_output$span,
                      grid.ppm.from = gui_output$grid.ppm.from,
                      grid.ppm.to = gui_output$grid.ppm.to,
                      grid.ppm.by = gui_output$grid.ppm.by,
                      grid.nsd.from = gui_output$grid.nsd.from,
                      grid.nsd.to = gui_output$grid.nsd.to,
                      grid.nsd.by = gui_output$grid.nsd.by,
                      grid.subset = gui_output$grid.subset,
                      grid.n = gui_output$grid.n)
        }
    }
    invisible(NULL)
}

## Rnews 2001-3.pdf for packer/grid intro
## use grid to align Tab2/3 elements

.GUI <- function(n = 1, e) {
    ## Create a variable to keep track of the state of the dialog window:
    ##  If the window is active,                                           done = 0
    ##  If the window has been closed using the Run button,                done = 1
    ##  If the window has been closed using the Close button or destroyed, done = 2
    done <- tclVar(0) 
    fls <- new.env()
    tcltk::tclRequire("BWidget")
    tt <- tktoplevel() ## toplevel windowx
    tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
    tkwm.title(tt,"synapter")
    ## ===============================================
    upperframe <- tkframe(tt) ## notebook with 3 tabs
    lowerframe <- tkframe(tt) ## run/close buttons
    ## ===============================================
    ## .tt.upperframe - notebook
    tkpack(upperframe, lowerframe, side = "top", expand = TRUE, fill = "both")
    notebook <- tk2notebook(upperframe,
                            tabs = c("Files", "Model","Grid"))
    tb1 <- tk2notetab(notebook, "Files")
    tb2 <- tk2notetab(notebook, "Model")
    tb3 <- tk2notetab(notebook, "Grid")
    tkpack(notebook, fill = "both", expand = TRUE)
    ## -----------------------------------------------
    ## .tt.upperframe.tb1 - Files
    tb1.upperframe <- tkframe(tb1)
    add <- tk2button(tb1.upperframe,
                     text = "Add",
                     command = function() {  
                         x <- tclvalue(tcl(treeWidget,"selection","get"))
                         if (x != "" &
                             (length(grep("Quant", x)) == 1 |
                                  length(grep("Ident", x)) == 1 |
                                      length(grep("fasta", x)) == 1 ) &
                             length(grep("Val", x)) == 0) {
                             f <- tk_choose.files("",
                                                  caption = "Select the corresponding file",
                                                  multi = FALSE)
                             val <- paste0(x, "Val")
                             tkinsert(treeWidget, "end", x, val, text = f)
                             assign(x, f, envir = fls)
                         } else if (x != "" &
                                    (length(grep("Out", x)) == 1 ) &
                                    length(grep("Val", x)) == 0) {
                             f <- tk_choose.dir("", caption = "Select the output folder")
                             val <- paste0(x, "Val")
                             tkinsert(treeWidget, "end", x, val, text = f)
                             assign(x, f, envir = fls)
                         } 
                     })
    remove <- tk2button(tb1.upperframe,
                        text = "Remove",
                        command = function() {
                            x <- tclvalue(tcl(treeWidget,"selection","get"))
                            if (x != "" & length(grep("Val$", x)) == 1) {
                                tkdelete(treeWidget, x)
                            }
                        })
    tkpack(tb1.upperframe)
    tkpack(add, side = "left", expand = TRUE)
    tkpack(remove, side = "right", expand = TRUE)
    ## -----------------------------------------------
    tb1.lowerframe <- tkframe(tb1)
    tkpack(tb1.lowerframe, expand = TRUE, fill = "both")
    xScr <- tkscrollbar(tb1.lowerframe,
                        command = function(...) tkxview(treeWidget,...),
                        orient = "horizontal")
    yScr <- tkscrollbar(tb1.lowerframe,
                        command = function(...) tkyview(treeWidget,...))
    treeWidget <- tkwidget(tb1.lowerframe,
                           "Tree",
                           xscrollcommand = function(...) tkset(xScr,...),
                           yscrollcommand = function(...) tkset(yScr,...))
    tkgrid(treeWidget, yScr)
    tkgrid.configure(treeWidget, sticky = "nsew")
    tkgrid.configure(yScr, sticky = "ns")
    tkgrid.columnconfigure(tb1.lowerframe, 0, weight = 1)
    tkgrid.rowconfigure(tb1.lowerframe, 0, weight = 1) 
    tkgrid.configure(xScr, sticky = "wse")
    tkinsert(treeWidget, "end", "root",
             "fasta", text = "fasta")
    tmp <- sapply(1:n, function(i) {
        rec <- paste0("Node", i)
        rectxt <- paste0("File set ", i)
        ident <- paste0("Ident", i)
        quant <- paste0("Quant", i)
        pep3d <- paste0("QuantPep3D", i)
        out <- paste0("Out", i)
        tkinsert(treeWidget,"end", "root", rec, text = rectxt)
        tkinsert(treeWidget,"end", rec, ident, text = "Ident")
        tkinsert(treeWidget,"end", rec, quant,  text = "Quant")
        tkinsert(treeWidget,"end", rec, pep3d, text = "Quant Pep3D")
        tkinsert(treeWidget,"end", rec, out, text = "Output folder")
    })    
    ## -----------------------------------------------
    tb1.bottomframe <- tkframe(tb1, relief = "groove", borderwidth = 1)
    tkpack(tb1.bottomframe, fill = "both", expand = TRUE, anchor = "s")
    master <- tclVar(0)
    cb.master <- tkcheckbutton(tb1.bottomframe)
    tkconfigure(cb.master, variable = master)
    masterlabel <- tklabel(tb1.bottomframe, text = "Master ")
    tkpack(masterlabel, cb.master, side = "left", expand = TRUE, padx = 5)
    ## ===============================================
    ## .tt.upperframe.tb2 - Model
    tb2.frame1 <- tkframe(tb2, relief = "groove", borderwidth = 1)
    tb2.frame2 <- tkframe(tb2, relief = "groove", borderwidth = 1)
    tkpack(tb2.frame1, expand = TRUE, fill = "both")
    tkpack(tb2.frame2, expand = TRUE, fill = "both")
    ## -----------------------------------------------  
    tkpack(tklabel(tb2.frame1, text = "Filtering"), side = "top",
           anchor = "w", pady = 1)
    .f1 <- tkframe(tb2.frame1)
    pepfdr <- tclVar(0.01)
    entry.pepfdr <- tkentry(.f1, textvariable = pepfdr, width = 10)
    tkpack(tklabel(.f1, text = "Peptide FDR", width = 10),
           entry.pepfdr, padx = 5, side = "left")
    .f2 <- tkframe(tb2.frame1)
    prtfpr <- tclVar(0.01)
    entry.prtfpr <- tkentry(.f2, textvariable = prtfpr, width = 10)
    tkpack(tklabel(.f2, text="Protein FPR", width = 10),
           entry.prtfpr, padx = 5, side = "left")
    .f3 <- tkframe(tb2.frame1)
    identppm <- tclVar(20)
    entry.identppm <- tkentry(.f3, textvariable = identppm, width = 10)
    tkpack(tklabel(.f3, text = "Ident ppm ", width = 10),
           entry.identppm, padx = 5, side = "left")
    .f4 <- tkframe(tb2.frame1)
    quantppm <- tclVar(20)
    entry.quantppm <- tkentry(.f4, textvariable = quantppm, width = 10)
    tkpack(tklabel(.f4,text="Quant ppm ", width = 10),
           entry.quantppm, padx = 5, side = "left")
    .f5 <- tkframe(tb2.frame1)
    uniquepep <- tclVar(1)
    cb.uniquepep <- tkcheckbutton(.f5)
    tkconfigure(cb.uniquepep, variable = uniquepep)
    tkpack(tklabel(.f5, text = "Unique peptides ", width = 15),
           cb.uniquepep, padx = 2, side = "left")
    tkpack(.f1, .f2, .f3, .f4, .f5, side = "top",
           padx = 1, pady = 1)
    tkpack(tklabel(tb2.frame2, text = "Retention model"), side = "top",
           anchor = "w", pady = 1)
    .f6 <- tkframe(tb2.frame2)
    span <- tclVar(0.05)
    entry.span <- tkentry(.f6, textvariable = span, width = 10) 
    tkpack(tklabel(.f6, text = "Loess span ", width = 15),
           entry.span, padx = 5, side = "left")
    tkpack(.f6, side = "top", padx = 1, pady = 1)
    ## ===============================================
    ## .tt.upperframe.tb3 - Grid
    tb3.frame1 <- tkframe(tb3, relief = "groove", borderwidth = 1)
    tb3.frame2 <- tkframe(tb3, relief = "groove", borderwidth = 1)
    tb3.frame3 <- tkframe(tb3, relief = "groove", borderwidth = 1)
    tkpack(tb3.frame1, expand = TRUE, fill = "both")
    tkpack(tb3.frame2, expand = TRUE, fill = "both")
    tkpack(tb3.frame3, expand = TRUE, fill = "both")
    ## -----------------------------------------------  
    tkpack(tklabel(tb3.frame1, text = "Mass tolerance (ppm)"), side = "top",
           anchor = "w", pady = 1)
    .f7 <- tkframe(tb3.frame1)
    ppm.from <- tclVar(5)
    entry.ppm.from <- tkentry(.f7, textvariable = ppm.from, width = 10)
    tkpack(tklabel(.f7, text = "from ", width = 5),
           entry.ppm.from, padx = 5, side = "left")  
    .f8 <- tkframe(tb3.frame1)
    ppm.to <- tclVar(20)
    entry.ppm.to <- tkentry(.f8, textvariable = ppm.to, width = 10)
    tkpack(tklabel(.f8, text = "to ", width = 5),
           entry.ppm.to, padx = 5, side = "left")
    .f9 <- tkframe(tb3.frame1)
    ppm.by <- tclVar(2)
    entry.ppm.by <- tkentry(.f9, textvariable = ppm.by, width = 10)
    tkpack(tklabel(.f9, text = "by ", width = 5),
           entry.ppm.by, padx = 5, side = "left")  
    tkpack(.f7, .f8, .f9, side = "top",
           padx = 1, pady = 1)
    ## -----------------------------------------------  
    tkpack(tklabel(tb3.frame2, text = "Retention time (nsd)"), side = "top",
           anchor = "w", pady = 1)
    .f10 <- tkframe(tb3.frame2)
    nsd.from <- tclVar(0.5)
    entry.nsd.from <- tkentry(.f10, textvariable = nsd.from, width = 10)
    tkpack(tklabel(.f10, text = "from ", width = 5),
           entry.nsd.from, padx = 5, side = "left")  
    .f11 <- tkframe(tb3.frame2)
    nsd.to <- tclVar(5)
    entry.nsd.to <- tkentry(.f11, textvariable = nsd.to, width = 10)
    tkpack(tklabel(.f11, text = "to ", width = 5),
           entry.nsd.to, padx = 5, side = "left")  
    .f12 <- tkframe(tb3.frame2)
    nsd.by <- tclVar(0.5)
    entry.nsd.by <- tkentry(.f12, textvariable = nsd.by, width = 10)
    tkpack(tklabel(.f12, text = "by ", width = 5),
           entry.nsd.by, padx = 5, side = "left")  
    tkpack(.f10, .f11, .f12, side = "top",
           padx = 1, pady = 1)
    ## -----------------------------------------------  
    tkpack(tklabel(tb3.frame3, text = "Data subset"),
           side = "top", anchor = "w", pady = 1)
    .f13 <- tkframe(tb3.frame3)
    subset <- tclVar(1)
    entry.subset <- tkentry(.f13, textvariable = subset, width = 10)
    tkpack(tklabel(.f13, text = "percentage ", width = 10),
           entry.subset, padx = 5, side = "left")  
    .f14 <- tkframe(tb3.frame3)
    grid.n <- tclVar(0)
    entry.grid.n <- tkentry(.f14, textvariable = grid.n, width = 10)
    tkpack(tklabel(.f14, text = "absolute  ", width = 10),
           entry.grid.n, padx = 5, side = "left")  
    tkpack(.f13, .f14, side = "top",
           padx = 1, pady = 1)
    ## ===============================================
    ## .tt.lowerframe - run/close buttons
    run <- tk2button(lowerframe,
                     text = "Run",                   
                     command = function() {
                         if (length(fls) != ((4 * n) + 1)) {
                             message("Not all input files have been provided")
                             tclvalue(done) <- 2
                             tkdestroy(tt)
                         }
                         if (as.numeric(tclvalue(grid.n)) != 0 &
                             as.numeric(tclvalue(subset)) !=0) {
                             message("Data subsetting via percentage or absolute - using percentage")
                             tclvalue(subset) <- 0
                         }
                         .gui_output <- list(filenames = unlist(as.list(fls)),
                                             master = as.logical(as.numeric(tclvalue(master))),
                                             pepfdr = as.numeric(tclvalue(pepfdr)),
                                             prtfpr = as.numeric(tclvalue(prtfpr)),
                                             identppm = as.numeric(tclvalue(identppm)),
                                             quantppm = as.numeric(tclvalue(quantppm)),
                                             uniquepep = as.logical(as.numeric(tclvalue(uniquepep))),
                                             span = as.numeric(tclvalue(span)),
                                             grid.ppm.from = as.numeric(tclvalue(ppm.from)),
                                             grid.ppm.to = as.numeric(tclvalue(ppm.to)),
                                             grid.ppm.by = as.numeric(tclvalue(ppm.by)),
                                             grid.nsd.from = as.numeric(tclvalue(nsd.from)),
                                             grid.nsd.to = as.numeric(tclvalue(nsd.to)),
                                             grid.nsd.by = as.numeric(tclvalue(nsd.by)),
                                             grid.subset = as.numeric(tclvalue(subset)),
                                             grid.n = as.numeric(tclvalue(grid.n)))
                         assign("gui_output", .gui_output, envir = e)
                         tclvalue(done) <- 1
                         tkdestroy(tt)
                     })
    close <- tk2button(lowerframe,
                       text = "Close",
                       command = function() {
                           tclvalue(done) <- 2
                           tkdestroy(tt)                       
                       })
    tkpack(run, side = "left", expand = TRUE)
    tkpack(close, side = "right", expand = TRUE)
    tkfocus(tt)
    tkwait.variable(done) ## wait for done != 0
}

