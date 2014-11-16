##' Checks version of an Synapter object
##'
##' IMPORTANT: update the minimalVersion value everytime you change
##' synapter-class.R (regardless if you just change the code in a method or the
##' complete API).
##'
##' @param object synapter object
##' @return TRUE if the object needs to be updated, FALSE otherwise
##' @noRd
.isSynapterObjectOutOfDate <- function(object) {
  minimalVersion <- as.package_version("1.7.2")

  version <- as.package_version(object$Version)

  isTRUE(version < minimalVersion)
}

##' Updates an old synapter object
##'
##' This function updates an old synapter object. Please see the details section
##' for more information about the changes in different versions of synapter
##' objects.
##'
##' \describe{
##'   \item{synapter 1.7.2}{
##'     Introduce \code{updateObject} method.\cr}
##' }
##'
##' @param object old synapter object
##' @param verbose print verbose output
##' @return new synapter object
##' @seealso \code{\link{Synapter}}
##' @examples
##' synapterTinyData()
##' synapterTiny
##' newSynapterObj <- updateObject(synapterTiny)
##' newSynapterObj
##'
##' @noRd
.updateSynapterObject <- function(object, verbose=TRUE) {

  if (.isSynapterObjectOutOfDate(object)) {
    newObject <- object$copy()

    newObject$Version <- as.character(packageVersion("synapter"))
    newObject$SynapterLog <- c(newObject$SynapterLog,
                               paste("Instance updated to synapter",
                                      newObject$Version, "on", date()))
    if (verbose) {
      message("You are using an old synapter object (", object$Version, "). ",
              "There are some internal changes in the definition of synapter ",
              "objects. Your object is updated to synapter ",
              newObject$Version, ". ")
    }

    return(newObject)
  }

  object
}

