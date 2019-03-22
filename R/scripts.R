#' Find Python versions and Conda Environments
#'
#' This function lists available python versions and conda
#' environments.
#' @param python If TRUE will list available python versions.
#' Default: TRUE
#' @param conda If TRUE will list available conda environments.
#' Default: TRUE
#' @return
#' @export
#' @examples find_python()
#' @import reticulate

find_python <- function(python=TRUE, conda=TRUE){
  f_py <- py_discover_config()
  condalist <- conda_list()
  if(python == TRUE & conda == TRUE){
    cat("Default Python:\n", f_py$python, "\n\n",
                "Python versions found:\n", f_py$python_versions, "\n\n",
                        "List of condaenvironments:\n",
                        paste0(condalist$name, ": ", condalist$python, ","))
  } else if(python == TRUE & conda == FALSE){
    cat("Default Python:\n", f_py$python, "\n\n",
                "Python versions found:\n", f_py$python_versions, "\n\n")
  } else if(python == FALSE & conda == TRUE){
    cat("List of condaenvironments:\n", 
                paste0(condalist$name, ": ", condalist$python, ","))
  } else {
    cat("Please set python and/or conda to TRUE")
  }
}

#' Set Up Python Configuration
#'
#' This function sets your python version and conda environment.
#' Run this command before PreprocessMatrix. Install python dependencies in the
#' same conda environment that you set here.
#' If conda environment does not exist, this function will make a new one.
#' @param condaenv Specify conda environment name
#' @param version Set path to specific version of python.
#' @param verbose TRUE or FALSE. When TRUE, shows python and conda configuration.
#' Default: TRUE
#' @return
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
    cat("\n Your specified condaenvironment, ", condaenv, 
                " was not found.\n Making new condaenvironment. \n")
    # create conda env
    conda_create(condaenv, packages = "python", conda = "auto")
    use_condaenv(condaenv = condaenv, required = T) # set conda environment
  } else {
    cat('\n Setting condaenvironment \n')
    use_condaenv(condaenv = condaenv, required = T) # set conda environment
  }
  if(verbose == TRUE){
    cat("\n Python/condaenvironment in use: \n")
    return(py_config())
  }
}


#' Install Python Dependencies
#'
#' This function installs needed python libraries into the specified conda
#' environment. Should be the same as the one specified in set_python.
#' On CentOS 7 pandas & scipy must be installed with pip install from 
#' the command line. Will get the error: /lib/libstdc++.so.6: version
#' `CXXABI_1.3.9' not found
#' @param condaenv Name of conda environment to install python libraries to.
#' @export
#' @examples py_depend(condaenv = "r-reticulate")
#' @import reticulate

py_depend <- function(condaenv){
  # symengine not installing on windows, will add later
  # issues with pandas & scipy on CentOS7
  required_py <- c("matplotlib", "numpy", "pandas", "scipy")
  conda_install(envname = condaenv, packages = required_py, pip = FALSE)
}



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
#' @param zero_perturb Flag to indicate you want to pertub the zero entries. 
#' Default: FALSE
#' @param threads Number of threads to use. Default: 1
#' @param verbose TODO
#' @return A list of with the following objects: matrix_size, column_names, row_names,
#' non_zero, num_switch_functions, asymptotic_stability
#' @export
#' @examples


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
    tt <- py_to_r(run_preproc(input_file, output_folder, prefix, max_bound, 
                              zero_perturb, threads, verbose))
    names(tt) <- c("matrix_size", "column_names", "row_names", "non_zero", 
                   "num_switch_functions", "asymptotic_stability")
    return(tt)
  }
}



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
#' @return If and input folder is specified the objects will be saved to that folder.
#' If the PreProsMatrix object is specified, an R list object with the following: 
#' matrix_size, column_names, row_names,
#' non_zero, num_switch_functions, asymptotic_stability, distributions,
#' expected_num_switch
#' @export
#' @examples


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
    entrywise <- py_to_r(run_EntryWise(NULL, prefix, distribution_type, 
                                       input_a, input_b, threads, 
                                       matrix_size=PreProsMatrix$matrix_size, 
                                       col_names=PreProsMatrix$column_names, 
                                       row_names=PreProsMatrix$row_names, 
                                       num_switch=PreProsMatrix$num_switch_functions, 
                                       asymp_stab=PreProsMatrix$asymptotic_stability))
    names(entrywise) <- c("distributions", "expected_num_switch")
    combined <- c(PreProsMatrix, entrywise)
    return(combined) 
  }
}


#' Generate Entry Wise Figures
#' 
#' This function plots the number of mis-predictions versus perturbation
#' value, overlaid with distribution over stable pertrubation values.
#' Run after ComputeEntryWisePerturbationExpectation()
#' @param input_folder Input folder. The location of the files created by PreprocessMatrix
#' if you specified an output_folder. This is also where the num switch array was saved. 
#' Must specify an input_folder OR EntryWise object. Default: NULL
#' @param EntryWise Object where the ComputeEntryWisePerturbationExpectation output was saved.
#' @param prefix Prefix of output files, if you so choose.
#' @param all_numswitch_plots set to TRUE if you ant to plot all num switch plots 
#' (potentially very large). Default: FALSE
#' @param list_of_numswitch_to_plot List of entries you want visualized with num switch. 
#' Should be a list of vectors. Example: list(c(0, 0), c(0, 1))
#' @return 
#' @export
#' @examples
#' @import ggplot2

GenerateEntryWiseFigures <- function(input_folder=NULL, EntryWise=NULL, prefix=NULL, 
                                     all_numswitch_plots = FALSE,
                                     list_of_numswitch_to_plot){
  # check parameters
  if((all_numswitch_plots == FALSE) & length(list_of_numswitch_to_plot) < 1){
    stop('No plots selected. Please set all_numswitch_plots or list_of_numswitch_to_plot')
  }
  if(is.null(input_folder) & is.null(EntryWise)){
    stop('No input provided. Please set input_folder path or 
         provide EntryWise object')
  }
  NumSwitch <- reticulate::import_from_path(
    "NumSwitch", 
    system.file("python", package = "PressPurtCoreAlg"), 
    convert = T)
  if(!is.null(input_folder)){
    # import NumSwitch and other python modules
    np <- import("numpy")
    os <- import("os")
    pickle <- import("pickle")
    sys <- import("sys")
    scipy.stats <- import("scipy.stats")
    # input files
    output_folder <- input_folder
    asymp_stab_file <- paste(output_folder, "asymptotic_stability.npy", sep = '')
    num_switch_file <- paste(output_folder, "num_switch_funcs.pkl", sep = '')
    exp_num_switch_file <- paste(output_folder, "expected_num_switch.csv", sep = '')
    distribution_file <- paste(output_folder, "distributions.pkl", sep = '')
    matrix_size_file <- paste(output_folder, "size.npy", sep = '')
    row_names_file <- paste(output_folder, "row_names.txt", sep = '')
    column_names_file <- paste(output_folder, "column_names.txt", sep = '')
    num_nonzero_file <- paste(output_folder, "num_non_zero.npy", sep = '')
    # load matrix_size_file
    n <- as.integer(np$load(matrix_size_file))
    m <- n
    # read in rownames
    row_names <- .check_read_file(row_names_file)
    # read in colnames
    column_names <- .check_read_file(column_names_file)
    # If appropriate, get the indices to plot
    # skipping this part for now
    # load the intervals of stability
    intervals <- np$load(asymp_stab_file)
    # load the num_switch functions
    num_switch_funcs <- r_to_py(py_load_object(num_switch_file))
    # load the distributions
    dists <- py_load_object(distribution_file)
  } else {
    n <- EntryWise$matrix_size
    m <- n
    column_names <- EntryWise$column_names
    row_names <- EntryWise$row_names 
    num_switch_funcs <- EntryWise$num_switch_functions 
    intervals <- EntryWise$asymptotic_stability
    dists <- EntryWise$distributions
  }
  if(length(list_of_numswitch_to_plot) > 0){
    sub_list <- do.call(rbind, list_of_numswitch_to_plot)
    rownames(sub_list) <- paste("(", sub_list[,1], ", ", sub_list[,2], ")", sep = '')
    dists.2 <- dists[names(dists) %in% rownames(sub_list)]
    nf.df.l <- list()
    plots <- list()
    for(k in 0:(m-1)){
      for(l in 0:(n-1)){
        if(paste("(", k, ", ", l, ")", sep = '') %in% names(dists.2)){
          ax2 <- get_distributions(matrix_interval = c(k, l), 
                                   distribution_list = dists, 
                                   asymp_stab = intervals)
          dist_vals <- ax2$dist_vals
          x_range <- ax2$x_range
          ns <- NumSwitch$num_switch_to_step(num_switch_funcs, intervals, k, l)
          ns.df <-  data.frame(k = k, l = l, nsx = do.call(c, ns[[1]]), 
                               nsy = do.call(c, ns[[2]]))
          ax2$dist_vals.trans <- ax2$dist_vals/(max(ax2$dist_vals)/max(ns.df$nsy))
          nf.df.l[[paste0(k, ".", l)]] <- ns.df
          # second y axis: max(dist_vals) = max(nsy)*x
          # where x is the fraction to get axis to match up
          p <- ggplot() + geom_step(data = ns.df, aes(nsx, nsy), size = 1.25, direction = 'vh') + 
            geom_area(data = ax2, aes(x=x_range, y= dist_vals.trans),
                      alpha = .2, color = "grey") + 
            theme_bw() +
            ggtitle(paste("(", ns.df$k[1], ", ", ns.df$l, ")", " entry", sep = '')) +
            scale_x_continuous("Epsilon value") +
            scale_y_continuous("Number of incorrect predictions",
                               sec.axis = sec_axis(~ . * (max(ax2$dist_vals)/max(ns.df$nsy)), 
                                                   name = "Probability density")) +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.line.y.right = element_line(color = "grey"), 
                  axis.ticks.y.right = element_line(colour = "grey"),
                  axis.text.y.right = element_text(color = "grey"),
                  axis.title.y.right = element_text(color = "grey"))
          plots[[paste0(k, ".", l)]] <- p
          print(p)
        }
      }
    }
  }
  if(all_numswitch_plots == TRUE){
    nf.df.l <- list()
    plots <- list()
    for(k in 0:(m-1)){
      for(l in 0:(n-1)){
        if(paste("(", k, ", ", l, ")", sep = '') %in% names(dists)){
          ax2 <- get_distributions(matrix_interval = c(k, l), 
                                         distribution_list = dists, 
                                         asymp_stab = intervals)
          dist_vals <- ax2$dist_vals
          x_range <- ax2$x_range
          ns <- NumSwitch$num_switch_to_step(num_switch_funcs, intervals, k, l)
          ns.df <-  data.frame(k = k, l = l, nsx = do.call(c, ns[[1]]), 
                               nsy = do.call(c, ns[[2]]))
          ax2$dist_vals.trans <- ax2$dist_vals/(max(ax2$dist_vals)/max(ns.df$nsy))
          nf.df.l[[paste0(k, ".", l)]] <- ns.df
          # second y axis: max(dist_vals) = max(nsy)*x
          # where x is the fraction to get axis to match up
          p <- ggplot() + geom_step(data = ns.df, aes(nsx, nsy), size = 1.25, direction = 'vh') + 
            geom_area(data = ax2, aes(x=x_range, y= dist_vals.trans),
                      alpha = .2, color = "grey") + 
            theme_bw() +
            ggtitle(paste("(", ns.df$k[1], ", ", ns.df$l, ")", " entry", sep = '')) +
            scale_y_continuous(sec.axis = sec_axis(~ . * (max(ax2$dist_vals)/max(ns.df$nsy)))) +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.line.y.right = element_line(color = "grey"), 
                  axis.ticks.y.right = element_line(colour = "grey"),
                  axis.text.y.right = element_text(color = "grey"),
                  axis.title.y.right = element_text(color = "grey"),
                  # remove for single plot
                  axis.title.x = element_blank(), 
                  axis.title.y.left = element_blank())
          plots[[paste0(k, ".", l)]] <- p
        }
      }
    }
    filled.grid <- .grid_helper(n,m, plots)
    gridExtra::grid.arrange(
      grobs = plots,
      top = "Number of mis-predictions versus perturbation value,\noverlaid with distribution over stable perturbation values",
      left = "Number of incorrect predictions",
      right = grid::textGrob("Probability density", rot =270, gp=grid::gpar(col="grey")),
      layout_matrix = filled.grid
    )  
  }
}

############ Helper functions

#' Get PDF distribution
#'
#' This function retrieves the PDF (Probablity Distribution Function)
#' object from the scipy method 
#' <scipy.stats._distn_infrastructure.rv_frozen>.
#' @param matrix_interval Position in the matrix. Example: c(0, 0)
#' @param distribution_list list of scipy distribution
#' @param asymp_stab asymptotic stability matrix
#' @export
#' @examples ax2 <- get_distributions(matrix_interval = c(0, 0), distribution_list, asymp_stab)
#' @import reticulate

get_distributions <- function(matrix_interval, distribution_list, asymp_stab){
  np <- reticulate::import("numpy")
  k <- matrix_interval[1]
  l <- matrix_interval[2]
  interval <- asymp_stab[(k+1), (l+1),] # change back to R 1 based
  padding <- (interval[2] - interval[1])/100
  x_range <- np$linspace((interval[1] - padding), (interval[2] + padding), 250)
  dist.py <- distribution_list[paste("(", k, ", ", l, ")", sep = '')]
  dist_vals <- sapply(x_range, function(x){dist.py[[1]]$pdf(x)})
  ax2 <- data.frame(x_range = as.numeric(x_range), dist_vals = dist_vals)
  return(ax2)
}


.check_read_file <- function(file_path){
  if(!file.exists(file_path)){
    # check if the file exists
    stop(paste("The file, ", file_path, ", does not appear to exist.", sep = ''))
  }
  tr <- data.table::fread(file_path)
  return(tr)
}

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
