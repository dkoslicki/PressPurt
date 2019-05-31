#' Find Python versions, Conda, & Virtual Environments
#'
#' This function lists available python versions, conda
#' environments, and virtual environments. One may show
#' all three or just one.
#' @param python If TRUE will list available python versions.
#' Default: TRUE
#' @param conda If TRUE will list available conda environments.
#' Default: TRUE
#' @param virtualenv If TRUE will list available virtual environments.
#' Default: TRUE
#' @return Returns a list of python versions, conda envs and 
#' virtual envs.
#' @export
#' @examples find_python()
#' @import reticulate

find_python <- function(python=TRUE, conda=TRUE, virtualenv=TRUE){
  f_py <- py_discover_config()
  condalist <- conda_list()
  virtlist <- virtualenv_list()
  if(python == TRUE & conda == TRUE & virtualenv == TRUE){
    cat("Default Python:\n", f_py$python, "\n\n",
        "Python versions found:\n", f_py$python_versions, "\n\n",
        "List of condaenvironments:\n",
        paste0(condalist$name, ": ", condalist$python, ","),
        "\n\n",
        "List of virtualenvs:\n", virtlist)
  } else if(python == TRUE & conda == FALSE & virtualenv == FALSE){
    cat("Default Python:\n", f_py$python, "\n\n",
        "Python versions found:\n", f_py$python_versions, "\n\n")
  } else if(python == FALSE & conda == TRUE & virtualenv == FALSE){
    cat("List of condaenvironments:\n", 
        paste0(condalist$name, ": ", condalist$python, ","))
  } else if(python == FALSE & conda == FALSE & virtualenv == TRUE){
    cat("List of virtualenvs:\n", virtlist)
  } else {
    cat("Please set python, virtualevl, and/or conda to TRUE")
  }
}

#' Set Up Python Conda Configuration
#'
#' This function sets your python version and conda environment.
#' Run this command before PreprocessMatrix. Install python dependencies in the
#' same conda environment that you set here.
#' If conda environment does not exist, this function will make a new one.
#' If the conda environment already exists, no need to set the python version.
#' When making a new conda environment, if the python version isn't set, then
#' your default one will be used.
#' @param condaenv Specify conda environment name
#' @param version Set path to specific version of python.
#' @param verbose TRUE or FALSE. When TRUE, shows python and conda configuration.
#' Default: TRUE
#' @export
#' @examples set_python(version = "~/anaconda3/bin/python", condaenv = "r-reticulate", verbose = TRUE)
#' @import reticulate

set_python <- function(condaenv, version=NULL, verbose = TRUE){
  if(!is.null(version)){
    # set python version
    use_python(version, required = T)
  } else {
    f_py <- py_discover_config()
    cat("No python version set. Default python is:\n", f_py$python, "\n")
  }
  # check if conda environment exists, if not make it
  condalist <- conda_list()
  if(!(condaenv %in% condalist$name)){
    cat("\n Your specified conda environment, ", condaenv, 
        " was not found.\n Making new conda environment. \n")
    # create conda env
    conda_create(condaenv, packages = "python", conda = "auto")
    # set conda environment
    use_condaenv(condaenv = condaenv, required = T) 
  } else {
    cat('\n Setting condaenvironment \n')
    # set conda environment
    use_condaenv(condaenv = condaenv, required = T)
  }
  if(verbose == TRUE){
    cat("\n Python/conda environment in use: \n")
    return(py_config())
  }
}


#' Set Up Virtual Environment for Python Configuration
#'
#' This function sets your python version and Virtual environment.
#' Run this command before PreprocessMatrix. Install python dependencies in the
#' same Virtual environment that you set here.
#' If the virtual environment does not exist, this function will make a new one.
#' If the virtual environment already exists, no need to set the python version.
#' When making a new virtual environment, if the python version isn't set, then
#' your default one will be used.
#' @param virtualenv Specify conda environment name
#' @param version Set path to specific version of python.
#' @param verbose TRUE or FALSE. When TRUE, shows python and conda configuration.
#' Default: TRUE
#' @export
#' @examples set_python_virtual(version = "/usr/bin/python3", virtualenv = "r-reticulate", verbose = TRUE)
#' @import reticulate

set_python_virtual <- function(virtualenv, version=NULL, verbose = TRUE){
  if(!is.null(version)){
    # set python version
    use_python(version, required = T)
  } else {
    f_py <- py_discover_config()
    cat("No python version set. Default python is:\n", f_py$python, "\n")
  }
  # check if virtual environment exists, if not make it
  virtlist <- virtualenv_list()
  if(!(virtualenv %in% virtlist)){
    cat("\n Your specified virtual environment, ", virtualenv, 
        " was not found.\n Making new virtual environment. \n")
    # create virtual env
    virtualenv_create(envname = virtualenv)
    # set virtual environment
    use_virtualenv(virtualenv = virtualenv, required = T) 
  } else {
    cat('\n Setting virtual environment \n')
    use_virtualenv(virtualenv = virtualenv, required = T)
  }
  if(verbose == TRUE){
    cat("\n Python/virtual environment in use: \n")
    return(py_config())
  }
}



#' Install Python Dependencies
#'
#' This function installs needed python libraries into the specified conda
#' environment OR virtual environment. Should be the same as the one 
#' specified in set_python.
#' Required python libraries: matplotlib, numpy, pandas, pathos,
#' scipy and sympy
#' On CentOS 7 pandas & scipy may need to be installed with pip install 
#' from the command line. Will get the error: /lib/libstdc++.so.6: version
#' `CXXABI_1.3.9' not found
#' See vingette for more information.
#' @param condaenv Name of conda environment to install python libraries to.
#' Default: NULL
#' @param virtualenv Name of virtual environment to install python libraries to.
#' Default: NULL
#' @export
#' @examples py_depend(condaenv = "r-reticulate", virtualenv = NULL)
#' @examples py_depend(virtualenv = "r-reticulate", condaenv = NULL)
#' @import reticulate

py_depend <- function(condaenv=NULL, virtualenv=NULL){
  # symengine not installing on windows, will add later
  # issues with pandas & scipy on CentOS7
  required_py <- c("matplotlib", "numpy", "pandas", "pathos", "scipy", "sympy")
  if(!is.null(condaenv)){
    conda_install(envname = condaenv, packages = required_py, pip = FALSE)
  } else if(!is.null(virtualenv)){
    virtualenv_install(envname = virtualenv, 
                       packages = required_py, ignore_installed = TRUE)
  } else {
    # add spedify conda or virt
  }
  
}


#' Convert data to R format if saved to disk
#'
#' This function will convert objects saved to disk to
#' R friendly objects, or the same output as 
#' ComputeEntryWisePerturbationExpectation.
#' If you used the "save to disk" option or ran via python
#' directly, run this function to read the data into R.
#' Files read in: asymptotic_stability.npy, column_names.txt,
#' distributions.pkl, expected_num_switch.csv,
#' num_non_zero.npy, num_switch_funcs.pkl, 
#' row_names.txt and size.npy.
#' Note how most of these objects are python based objects-
#' numpy or pickle objects. 
#' @param matrix path to the original matrix.
#' @param type csv or tab. Is the oringal matrix comma separated
#' or tab separated?
#' @param folder path to the folder where output data was saved.
#' @return object formatted in the same way the output of 
#' ComputeEntryWisePerturbationExpectation
#' @export
#' @examples data <- process_data(matrix = infile, type = "csv", folder = "test_r/test3")
#' @import reticulate

process_data <- function(matrix, type, folder){
  np <- reticulate::import("numpy")
  if(type == "csv"){
    original_matrix <- as.matrix(read.csv(matrix, header = F))
    colnames(original_matrix) <- 1:ncol(original_matrix)
  } else if (type == "tab"){
    original_matrix <- as.matrix(read.table(matrix, sep = "\t"))
    colnames(original_matrix) <- 1:ncol(original_matrix)
  } else {
    cat("Please specify csv file or tab separated file in type")
    stop()
  }
  out_files <- list.files(folder, full.names = T)
  # matrix size
  matrix_size <- as.integer(np$load(out_files[
    grep("size.npy", out_files)]))
  # col names
  column_names <- as.character(read.table(out_files[
    grep("column_names.txt", out_files)])$V1)
  # row names
  row_names <- as.character(read.table(out_files[
    grep("row_names.txt", out_files)])$V1)
  # non zero
  non_zero <- as.integer(np$load(out_files[
    grep("num_non_zero.npy", out_files)]))
  # get AS and split into start and stop
  asymptotic_stability <- np$load(out_files[
    grep("asymptotic_stability.npy", out_files)])
  asymptotic_stability_start <- asymptotic_stability[,,1]
  asymptotic_stability_end <- asymptotic_stability[,,2]
  # load num switch funcs
  num_switch_functions <- reticulate::py_load_object(out_files[
    grep("num_switch_funcs.pkl", out_files)])
  # convert to R format
  names(num_switch_functions) <- .r_index(
    names = num_switch_functions, to_r = T)
  num_switch_funcs_r <- lapply(names(num_switch_functions), function(x) 
    .NS_func_r(num_switch_funcs = num_switch_functions, name = x))
  names(num_switch_funcs_r) <- names(num_switch_functions)
  # expected num switch
  expected_num_switch <- read.csv(out_files[
    grep("expected_num_switch.csv", out_files)], 
    header = T, row.names = 1)
  colnames(expected_num_switch) <- 1:ncol(expected_num_switch)
  rownames(expected_num_switch) <- 1:nrow(expected_num_switch)
  # distributions
  distributions <- reticulate::py_load_object(out_files[
    grep("distributions.pkl", out_files)])
  names(distributions) <- .r_index(
    names = distributions, to_r = T)
  # get in format
  distributions_object <- lapply(names(distributions), function(x){
    split_name <- unlist(lapply(
      strsplit(gsub("\\(", "", gsub(")", "", x)), ","), as.numeric))
    k <- split_name[1]
    l <- split_name[2]
    single_dist <- get_distributions_single(
      matrix_entry = c(k,l), 
      distribution_list = distributions, 
      asymp_stab = c(asymptotic_stability_start[k,l], 
                     asymptotic_stability_end[k,l]))
    return(single_dist)
  })
  names(distributions_object) <- names(distributions)
  # place everything in list
  output <- list(original_matrix, matrix_size, 
                 column_names, row_names, non_zero,
                 asymptotic_stability_start,
                 asymptotic_stability_end,
                 num_switch_funcs_r, distributions_object)
  names(output) <- c("original_matrix", "matrix_size",
                     "column_names", "row_names",
                     "non_zero", 
                     "asymptotic_stability_start",
                     "asymptotic_stability_end",
                     "num_switch_funcs_r",
                     "distributions_object")
  return(output)
}

############ Helper functions

#' Get PDF distribution
#'
#' This function retrieves the PDF (Probablity Distribution Function)
#' object from the scipy method 
#' <scipy.stats._distn_infrastructure.rv_frozen>.
#' @param matrix_entry Position in the matrix. Example: c(0, 0) # matrix entry 
#' @param distribution_list list of scipy distribution
#' @param asymp_stab asymptotic stability interval
#' @param points the number of values in x range
#' @export
#' @examples ax2 <- get_distributions(matrix_entry = c(0, 0), distribution_list, asymp_stab)
#' @import reticulate

get_distributions <- function(matrix_entry, 
                              distribution_list, asymp_stab,
                              points=250){
  np <- reticulate::import("numpy")
  k <- matrix_entry[1]
  l <- matrix_entry[2]
  interval <- asymp_stab[(k+1), (l+1),] # change back to R 1 based
  padding <- (interval[2] - interval[1])/100
  x_range <- np$linspace((interval[1] - padding), (interval[2] + padding), points)
  dist.py <- distribution_list[paste("(", k, ", ", l, ")", sep = '')]
  dist_vals <- sapply(x_range, function(x){dist.py[[1]]$pdf(x)})
  ax2 <- data.frame(x_range = as.numeric(x_range), dist_vals = dist_vals)
  return(ax2)
}

# new version of get_distribution?

get_distributions_single <- function(matrix_entry, 
                                     distribution_list, asymp_stab,
                                     points=250){
  np <- reticulate::import("numpy")
  k <- matrix_entry[1]
  l <- matrix_entry[2]
  interval <- asymp_stab
  padding <- (interval[2] - interval[1])/points
  x_range <- np$linspace((interval[1] - padding), (interval[2] + padding), points)
  dist.py <- distribution_list[paste("(", k, ", ", l, ")", sep = '')]
  dist_vals <- sapply(x_range, function(x){dist.py[[1]]$pdf(x)})
  ax2 <- data.frame(x_range = as.numeric(x_range), dist_vals = dist_vals)
  out_list <- list(position = c(k, l), interval = interval,
                   distribution = ax2)
  return(out_list)
}

# check if file exists

.check_read_file <- function(file_path){
  if(!file.exists(file_path)){
    # check if the file exists
    stop(paste("The file, ", file_path, ", does not appear to exist.", sep = ''))
  }
  tr <- data.table::fread(file_path)
  return(tr)
}

# Helper function for plotting

.grid_helper <- function(n, m, plots){
  # info for grid
  matrix.grid <- data.frame(n=n, m=m, total=n*m)
  # empty matrix for gridExtra::grid.arrange
  empty.grid <- matrix(as.numeric(NA), nrow = matrix.grid$n, ncol = matrix.grid$m)
  # now need to set up matrix position to list name
  layout.matrix <- data.frame(do.call(rbind, strsplit(names(plots), "[.]")))
  colnames(layout.matrix) <- c("rows", "cols")
  layout.matrix$index <- rownames(layout.matrix)
  # change to 1 based index
  layout.matrix$rows.r <- as.numeric(as.character(layout.matrix$rows)) + 1
  layout.matrix$cols.r <- as.numeric(as.character(layout.matrix$cols)) + 1
  # populate empty.grid
  for(i in 1:nrow(layout.matrix)){
    empty.grid[layout.matrix$rows.r[i], layout.matrix$cols.r[i]] <- as.numeric(layout.matrix$index[i])
    #print(layout.matrix$index[i])
  }
  return(empty.grid)
}

# recombine asymptotic_stability into 1 object

.back_to_array <- function(as_start, as_stop){
  ar <- array(dim = c(nrow(as_start), ncol(as_start), 2))
  ar[,,1] <- as_start
  ar[,,2] <- as_stop
  return(ar)
}

# add 1 (r-based 1) or subtract 1 (python 0 based)

.r_index <- function(names, to_r=TRUE){
  if (to_r == T) {
    index <- 1
  } else if (to_r == F) {
    index <- -1
  } else {
    print("Set to_r to T or F")
  }
  # remove par, split str, convert to num
  remove_par <- gsub("\\(", "", gsub(")", "", names(names)))
  remove_par <- lapply(strsplit(remove_par, ","), as.numeric)
  converted <- lapply(remove_par, function(x){return(x + index)})
  # change back to character vector
  back_ch <- lapply(converted, function(x){
    tel <- paste0("(", x[1], ", ", x[2], ")")
    return(tel)})
  return(unlist(back_ch))
}

# unlist num_switch_funcs

.NS_func_r <- function(num_switch_funcs, name){
  num_test <- do.call(rbind, lapply(num_switch_funcs[[name]], unlist))
  rownames(num_test) <- num_test[,1]
  colnames(num_test) <- c("num_switch_val", "start", "end")
  return(num_test)
}
