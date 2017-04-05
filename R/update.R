##' Checks version of an Synapter object
##'
##' @param object synapter object
##' @return FALSE if the object needs to be updated, TRUE otherwise
##' @noRd
.isCurrent <- function(object) {
  ## synpater class version < 2.0.0 does not have a ClassVersion field
  if (!"ClassVersion" %in% ls(object)) {
    return(FALSE)
  }

  ## synapter class version >= 2.0.0
  version <- as.package_version(object$ClassVersion)

  isTRUE(version == .synapterClassVersion)
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
.updateSynapterObject <- function(object, verbose=interactive()) {

  if (!isCurrent(object)) {
    newObject <- object$copy()

    ## 1.99 was Function before: see issue #67
    colnames(newObject$QuantPep3DData)[
      colnames(newObject$QuantPep3DData) == "Function"] <- "matchedEMRTs"

    newObject$ClassVersion <- .synapterClassVersion
    newObject$Version <- as.character(packageVersion("synapter"))
    newObject$SynapterLog <- c(newObject$SynapterLog,
                               paste("Instance updated to synapter",
                                      newObject$Version, "on", date()))
    if (verbose) {
      message("You are using an old synapter object. ",
              "There are some internal changes in the definition of synapter ",
              "objects. Your object is updated to synapter ",
              newObject$Version, ". ")
    }

    return(newObject)
  }

  object
}
