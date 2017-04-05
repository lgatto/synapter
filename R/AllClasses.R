setClass("MasterPeptides",
         contains = c("Versioned"),
         representation = representation(
           fdr = "numeric",
           method = "character",
           ## fpr = "numeric",
           masters = "list",
           pepfiles = "character",
           orders = "list",
           fragmentfiles = "character",
           fragments = "data.frame",
           fragmentlibrary = "MSnExp"),
         prototype = prototype(
          new("Versioned",
              versions = c(MasterPeptides = "2.0.0"))),
         validity = function(object) {
           msg <- NULL
           if (!all(isCurrent(object)))
             msg <- validMsg(msg, paste0("MasterPeptides object is out of date. ",
                                         "Please run ",
                                         sQuote("object <- updateObject(object)"), "."))
           if (length(object@masters) != length(object@orders))
             msg <- validMsg(msg, "Master dataframes and orders have different length.")
           if (length(object@orders[[1]]) != length(object@orders[[2]]))
             msg <- validMsg(msg, "Master orders have different length.")
           if (length(object@orders[[1]]) != length(object@pepfiles))
             msg <- validMsg(msg, "Peptides files and orders have different length.")
           if (length(object@fragmentfiles) &&
               length(object@fragmentfiles) != length(object@pepfiles))
             msg <- validMsg(msg, "Peptides and fragment files have different length.")
           if (is.null(msg)) TRUE
           else msg
         })

setClass("MasterFdrResults",
         contains = c("Versioned"),
         representation = representation(
           all = "data.frame",
           files = "character",
           best = "data.frame",
           masterFdr = "numeric"),
         prototype = prototype(
          new("Versioned",
              versions = c(MasterFdrResults = "2.0.0"))),
         validity = function(object) {
           msg <- NULL
           if (!all(isCurrent(object)))
             msg <- validMsg(msg, paste0("MasterFdrResults object is out of date. ",
                                         "Please run ",
                                         sQuote("object <- updateObject(object)"), "."))
           if (is.null(msg)) TRUE
           else msg
         })

