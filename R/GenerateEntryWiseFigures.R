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
                                     list_of_numswitch_to_plot = NULL){
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
    asymptotic_stability <- .back_to_array(as_start = EntryWise$asymptotic_stability_start,
                                           as_stop = EntryWise$asymptotic_stability_end)
    intervals <- asymptotic_stability
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
          ax2 <- get_distributions(matrix_entry = c(k, l), 
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
          ax2 <- get_distributions(matrix_entry = c(k, l), 
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
