#' Compute Multi Entry Perturbation Expectation
#'
#' This function takes a jacobian matrix and computes the multi-entry 
#' perturbation expectation.
#' @param input_file Input comma separated file for the jacobian matrix.
#' @param num_iterates Number of iterates in the Monte Carlo sampling to perform.
#' Default: 10000
#' @param interval_length Interval length over which to make the perturbations.
#' Default: 0.01
#' @param threads Number of threads to use. Default: 1
#' @return returns a scalar
#' @export
#' @examples
#' \dontrun{
#' infile <- system.file("extdata", "Modules", "IGP.csv", 
#'     package = "PressPurtCoreAlg") 
#' ComputeMultiEntryPerturbationExpectation(input_file = infile)
#' }

ComputeMultiEntryPerturbationExpectation <- function(
  input_file, num_iterates=1000, 
  interval_length=0.01, 
  threads=1){
  NaiveSS <- reticulate::import_from_path(
    "NaiveSS", 
    system.file("python", package = "PressPurtCoreAlg"), 
    convert = T)
  NumSwitch <- reticulate::import_from_path(
    "NumSwitch", 
    system.file("python", package = "PressPurtCoreAlg"), 
    convert = T)
  reticulate::source_python(system.file("python", 
                                        "ComputeMultiEntryPerturbationExpectation.py", 
                                        package = "PressPurtCoreAlg"), convert = F)
  MultiEntry <- py_to_r(run_MultiEntry(input_file=input_file, 
                                       num_iterates=num_iterates, 
                                       interval_length=interval_length, 
                                       threads=threads))
  return(MultiEntry)
}
