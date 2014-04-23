##' Updates an old synapter object
##'
##' This function updates an old synapter object. Please see the details section
##' for more information about the changes in different versions of synapter
##' objects.
##'
##' \describe{
##'   \item{synapter 1.7.1}{
##'     Introduce \code{updateSynapterObject} function.}
##' }
##'
##' @param obj synapter object
##' @param verbose print verbose output
##' @return new synapter object
##' @seealso \code{\link{Synapter}}
##' @examples
##' synapterTinyData()
##' synapterTiny
##' newSynapterObj <- updateSynapterObject(synapterTiny)
##' newSynapterObj
##'
updateSynapterObject <- function(obj, verbose=TRUE) {

  current <- as.package_version(obj$Version)

  minimalVersion <- as.package_version("0.0.1")

  isDeprecated <- current < minimalVersion

  newobj <- obj$copy()

  if (isDeprecated) {
    newobj$Version <- as.character(packageVersion("synapter"))
    newobj$SynapterLog <- c(newobj$SynapterLog,
                            paste("Instance updated to synapter",
                                  newobj$Version, "on", date()))
    if (verbose) {
      message("You are using an old synapter object (", obj$Version, "). ",
              "There are some internal changes in the definition of synapter ",
              "objects. Your object is updated to synapter ", newobj$Version,
              ". Please see ", sQuote("?updateSynapterObject"), " for details.")
    }
  }

  return(newobj)
}
