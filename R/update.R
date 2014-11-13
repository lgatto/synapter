##' Checks version of an Synapter object
##' @param object synapter object
##' @param throwError if TRUE throw an error if version is too old
##' @return new synapter object
##' @noRd
.isSynapterObjectOutOfDate <- function(object, throwError=TRUE) {
  minimalVersion <- as.package_version("1.7.2")

  version <- as.package_version(object$Version)

  isOutOfDate <- version < minimalVersion

  if (isOutOfDate && throwError) {
    stop("Your Synapter object is out of date. Please run ",
         sQuote("object <- updateObject(object)"), ".")
  }

  isOutOfDate
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

  if (.isSynapterObjectOutOfDate(object, throwError=FALSE)) {
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

