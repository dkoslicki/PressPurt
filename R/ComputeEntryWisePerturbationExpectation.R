#' Compute Entry Wise Perturbation Expectation
#'
#' This function computes the expected number of sign switches from perturbing 
#' each entry individually. Run after PreprocessMatrix(). 
#' @param input_folder Input folder. The location of the files created by PreprocessMatrix
#' if you specified an output_folder. If this option is specified, this
#' is also where the num switch array will be saved. Must specify an 
#' input_folder OR PreProsMatrix. Default: NULL
#' @param PreProsMatrix Object where the PreprocessMatrix output was saved.
#' Must specify an input_folder OR PreProsMatrix. Default: NULL
#' @param prefix Prefix of output files, if you so choose.
#' @param distribution_type Kind of distribution to use. Valid choices are: 
#' truncnorm, uniform, trunc_lognorm, beta. Default: “truncnorm”
#' @param a First parameter to the distribution you choose. For truncnorm, this is the mean.
#' Default: 0
#' @param b First parameter to the distribution you choose. For truncnorm, this is the variance.
#' Using a negative value indicates you want the standard deviation to be the length of the 
#' interval divided by the absolute value of the input parameter. Default: -2
#' @param threads Number of threads to use. Default: 1
#' @return If an input folder is specified the objects will be saved to that folder.
#' If the PreProsMatrix object is specified, an R list object with the following: 
#' original_matrix, matrix_size, column_names, row_names,
#' non_zero, num_switch_functions, asymptotic_stability_start, 
#' asymptotic_stability_end, num_switch_funcs_r, distributions,
#' expected_num_switch, distributions_object
#' @export
#' @examples Entrywise <- ComputeEntryWisePerturbationExpectation(PreProsMatrix = PreProsMatrix,
#' distribution_type = "truncnorm", input_a = 0, input_b = -2, threads = 1)


ComputeEntryWisePerturbationExpectation <- function(input_folder=NULL, 
                                                    PreProsMatrix = NULL,
                                                    prefix=NULL, 
                                                    distribution_type="truncnorm", 
                                                    input_a=0, input_b=-2, threads=1
){
  if(is.null(input_folder) & is.null(PreProsMatrix)){
    stop('No input provided. Please set input_folder path or 
         provide PreProsMatrix object')
  }
  # import NumSwitch, MRS & SS modules and source ComputeEntryWisePerturbationExpectation.py
  reticulate::import_from_path("NumSwitch", 
                               system.file("python", package = "PressPurtCoreAlg"), 
                               convert = F)
  reticulate::import_from_path("MRS", 
                               system.file("python", package = "PressPurtCoreAlg"), 
                               convert = F)
  reticulate::import_from_path("SS", 
                               system.file("python", package = "PressPurtCoreAlg"), 
                               convert = F)
  reticulate::source_python(system.file("python", 
                                        "ComputeEntryWisePerturbationExpectation.py", 
                                        package = "PressPurtCoreAlg"), convert = F)
  # If input/output folder is specified, don't output R object
  if(!is.null(input_folder)){
    entrywise <- py_to_r(run_EntryWise(input_folder, prefix, distribution_type, 
                                       input_a, input_b, threads))
    return(paste0("Output saved to input_folder: ", input_folder)) 
  } else{
    # back to single array for python
    asymptotic_stability <- .back_to_array(
      as_start = PreProsMatrix$asymptotic_stability_start,
      as_stop = PreProsMatrix$asymptotic_stability_end)
    # go back to python based
    #names(PreProsMatrix$num_switch_functions) <- .r_index(
    #  names = PreProsMatrix$num_switch_functions, to_r = F)
    # run main function
    entrywise <- py_to_r(run_EntryWise(NULL, prefix, distribution_type, 
                                       input_a, input_b, threads, 
                                       matrix_size=PreProsMatrix$matrix_size, 
                                       col_names=PreProsMatrix$column_names, 
                                       row_names=PreProsMatrix$row_names, 
                                       num_switch=PreProsMatrix$num_switch_functions_py, 
                                       asymp_stab=asymptotic_stability))
    names(entrywise) <- c("distributions", "expected_num_switch")
    combined <- c(PreProsMatrix, entrywise)
    #names(combined$num_switch_functions) <- .r_index(
    #  names = combined$num_switch_functions, to_r = T)
    names(combined$distributions) <- .r_index(
      names = combined$distributions, to_r = T)
    distributions_object <- lapply(names(combined$distributions), function(x){
      split_name <- unlist(lapply(
        strsplit(gsub("\\(", "", gsub(")", "", x)), ","), as.numeric))
      k <- split_name[1]
      l <- split_name[2]
      single_dist <- get_distributions_single(
        matrix_entry = c(k,l), 
        distribution_list = combined$distributions, 
        asymp_stab = c(combined$asymptotic_stability_start[k,l], 
                       combined$asymptotic_stability_end[k,l]))
      return(single_dist)
    })
    names(distributions_object) <- names(combined$distributions)
    combined$distributions_object <- distributions_object
    combined$distributions <- NULL
    ns_object <- lapply(names(combined$num_switch_funcs_r), function(x){
      split_name <- unlist(lapply(
        strsplit(gsub("\\(", "", gsub(")", "", x)), ","), as.numeric))
      k <- split_name[1]
      l <- split_name[2]
      num_switch_func <- combined$num_switch_funcs_r[
        paste("(", k, ", ", l, ")", sep = '')][[1]]
      ns_step <- ns_to_step(asymp_stab_start = combined$asymptotic_stability_start[k,l],
                            asymp_stab_end = combined$asymptotic_stability_end[k,l],
                            num_switch_func = num_switch_func)
      return(ns_step)
    })
    names(ns_object) <- names(combined$num_switch_funcs_r)
    combined$ns_object_plot <- ns_object
    return(combined)
  }
}
