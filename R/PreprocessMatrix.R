#' Preprocess Matrix
#'
#' This script pre-processes a matrix by figuring out what the intervals of asymptotic
#' stability are, as well as finding which perturbation values lead to a sign switch.
#' @param input_file Input comma separated file for the jacobian matrix.
#' @param output_folder Optional output folder to save python objects to disk. 
#' A number of files will be created in the form ‘output_folder/<prefix>_*.npy’.
#' Default is NULL.
#' @param prefix Prefix of output files, if you so choose.
#' @param max_bound some of the matrices are unbounded stable towards one end, this 
#' is the limit the user imposes. Default: 10
#' @param zero_perturb Flag to indicate you want to perturb the zero entries. 
#' Default: FALSE
#' @param threads Number of threads to use. Default: 1
#' @param verbose Default: FALSE
#' @return A list of with the following objects: matrix_size, column_names, row_names,
#' non_zero, num_switch_functions, asymptotic_stability_start,
#' asymptotic_stability_end, num_switch_funcs_r
#' @export
#' @examples
#' \dontrun{
#' infile <- system.file("extdata", "Modules", "IGP.csv", 
#'     package = "PressPurtCoreAlg") 
#' PreProsMatrix <- PreprocessMatrix(input_file = infile, 
#'     output_folder = NULL, max_bound = 10, threads = 2)
#' }


PreprocessMatrix <- function(input_file, output_folder=NULL, prefix=NULL, max_bound=10, 
                             zero_perturb=FALSE, threads=1, verbose=FALSE){
  # import NumSwitch module and source PreprocessMatrix.py
  reticulate::import_from_path("NumSwitch", 
                               system.file("python", package = "PressPurtCoreAlg"), 
                               convert = F)
  reticulate::source_python(system.file("python", "PreprocessMatrix.py", 
                                        package = "PressPurtCoreAlg"), 
                            convert = F)
  # If output folder is specified, don't output R object
  if(!is.null(output_folder)){
    tt <- py_to_r(run_preproc(input_file, output_folder, prefix, max_bound, 
                              zero_perturb, threads, verbose))
    return(paste0("Output saved to output_folder: ", output_folder)) 
  } else{
    py_temp <- run_preproc(input_file, output_folder, prefix, max_bound, 
                              zero_perturb, threads, verbose)
    tt <- list("original_matrix" = py_to_r(py_temp[[0]]),
               "matrix_size" = py_to_r(py_temp[[1]]),
               "column_names" = py_to_r(py_temp[[2]]),
               "row_names" = py_to_r(py_temp[[3]]),
               "non_zero" = py_to_r(py_temp[[4]]),
               "num_switch_functions_py" = py_temp[[5]],
               "num_switch_functions" = py_to_r(py_temp[[5]]),
               "asymptotic_stability" = py_to_r(py_temp[[6]]))
    #names(tt) <- c("original_matrix", "matrix_size", "column_names", "row_names", "non_zero", 
    #               "num_switch_functions", "asymptotic_stability")
    # separate AS matrix into start and end
    tt$asymptotic_stability_start <- tt$asymptotic_stability[,,1]
    tt$asymptotic_stability_end <- tt$asymptotic_stability[,,2]
    tt$asymptotic_stability <- NULL
    # r based names
    names(tt$num_switch_functions) <- .r_index(names = tt$num_switch_functions, to_r = T)
    # unlist num_switch_funcs
    tt$num_switch_funcs_r <- lapply(names(tt$num_switch_functions), function(x) 
      .NS_func_r(num_switch_funcs = tt$num_switch_functions, name = x))
    names(tt$num_switch_funcs_r) <- names(tt$num_switch_functions)
    tt$num_switch_functions <- NULL
    return(tt)
  }
}
