.onLoad <- function(libname = find.package("grattan"), pkgname = "grattan") {
  # to avoid CRAN "Note"
  utils::globalVariables(c("run_EntryWise", "run_MultiEntry", "run_preproc"))
}