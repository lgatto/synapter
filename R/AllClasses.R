setClass("MasterPeptides",
         representation = representation(
           fdr = "numeric",
           method = "character",
           ## fpr = "numeric",
           masters = "list",
           pepfiles = "character",
           orders = "list"),
         validity = function(object) {
           msg <- NULL           
           if (length(object@masters) != length(object@orders))
             msg <- c(msg, "Master dataframes and orders have different length.")
           if (length(object@orders[[1]]) != length(object@orders[[2]]))
             msg <- c(msg, "Master orders have different length.")
           if (length(object@orders[[1]]) != length(object@pepfiles))
             msg <- c(msg, "Peptides files and orders have different length.")
           if (is.null(msg)) TRUE           
           else msg
         })

setClass("MasterFdrResults",
         representation = representation(
           all = "data.frame",
           files = "character",
           best = "data.frame",
           masterFdr = "numeric"))

