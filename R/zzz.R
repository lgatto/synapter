.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("This is the 'synapter' version ",
                               packageVersion("synapter"), ".\n",
                              "Read '?synapter' to get started.\n", sep=""))
  addVigs2WinMenu("synapter")
}

