#' Generate Entry Wise Figures
#' 
#' This function plots the number of mis-predictions versus perturbation
#' value, overlaid with distribution over stable perturbation values.
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
#' @export
#' @examples
#' \dontrun{
#' # Set input file
#' infile <- system.file("extdata", "Modules", "IGP.csv", 
#'     package = "PressPurt")
#' # Preprocess the matrix
#' PreProsMatrix <- PreprocessMatrix(input_file = infile, 
#'     output_folder = NULL, max_bound = 10, threads = 2)
#'
#' # Run ComputeEntryWisePerturbationExpectation
#' Entrywise <- ComputeEntryWisePerturbationExpectation(PreProsMatrix = PreProsMatrix,
#'     distribution_type = "truncnorm", 
#'     input_a = 0, input_b = -2, threads = 1)
#'
#' # Plot specific entries using entrywise object
#' list_of_numswitch_to_plot <- list(c(1, 1), c(1, 2))
#' GenerateEntryWiseFigures(EntryWise=Entrywise, 
#'     all_numswitch_plots = FALSE, 
#'     list_of_numswitch_to_plot=list_of_numswitch_to_plot)
#'      
#' 
#' # Plot specific entries from folder
#' GenerateEntryWiseFigures(input_folder = "test_r/test3", 
#'     all_numswitch_plots = FALSE, 
#'     list_of_numswitch_to_plot=list_of_numswitch_to_plot)
#'
#' # Plot all numswitch plots
#' GenerateEntryWiseFigures(EntryWise=Entrywise, 
#'     all_numswitch_plots = TRUE)
#' }
#' @import ggplot2
#' @import utils


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
  if(!is.null(input_folder)){
    # Convert to R object
    # Need distributions and Num Switch
    EntryWise <- process_data(matrix = NULL, type = "csv", 
                              folder = input_folder, prefix = prefix)
  }
  n <- EntryWise$matrix_size
  m <- n
  if(length(list_of_numswitch_to_plot) > 0){
    sub_list <- do.call(rbind, list_of_numswitch_to_plot)
    rownames(sub_list) <- paste("(", sub_list[,1], ", ", sub_list[,2], ")", sep = '')
    if(length(which(rownames(sub_list) %in% names(EntryWise$ns_object_plot))) < 1){
      stop("No Matrix entries input found in list_of_numswitch_to_plot were found 
           in the available matrix entries. Perhaps you are thinking in 0-based python?")
    }
    dists <- EntryWise$distributions_object[rownames(sub_list)]
    ns <- EntryWise$ns_object_plot[rownames(sub_list)]
    plots <- list()
    for(i in 1:length(dists)){
      ns.df <- ns[[i]] # num switch
      entry <- sub_list[i,] # matrix entry
      dists.2 <- dists[[i]]$distribution # distributions
      # transform axis 2 to work with ggplot2
      # second y axis: max(dist_vals) = max(nsy)*x
      # where x is the fraction to get axis to match up
      dists.2$dist_vals.trans <- dists.2$dist_vals/
        (max(dists.2$dist_vals)/max(ns.df$nsy))
      p <- ggplot() + geom_step(data = ns.df, aes_string("nsx", "nsy"), size = 1.25, direction = 'vh') + 
        geom_area(data = dists.2, aes_string(x="x_range", y="dist_vals.trans"),
                  alpha = .2, color = "grey") + 
        theme_bw() +
        ggtitle(paste("(", entry[1], ", ", entry[2], ")", " entry", sep = '')) +
        scale_x_continuous("Epsilon value") +
        scale_y_continuous("Number of incorrect predictions",
                           sec.axis = sec_axis(~ . * (max(dists.2$dist_vals)/max(ns.df$nsy)), 
                                               name = "Probability density")) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.line.y.right = element_line(color = "grey"), 
              axis.ticks.y.right = element_line(colour = "grey"),
              axis.text.y.right = element_text(color = "grey"),
              axis.title.y.right = element_text(color = "grey"))
      plots[[paste0(entry[1], ".", entry[2])]] <- p
      print(p)
    }
    }
  if(all_numswitch_plots == TRUE){
    # get k and l from names
    remove_par <- gsub("\\(", "", gsub(")", "", names(EntryWise$distributions_object)))
    remove_par <- lapply(strsplit(remove_par, ","), as.numeric)
    sub_list <- do.call(rbind, remove_par)
    rownames(sub_list) <- names(EntryWise$distributions_object)
    dists <- EntryWise$distributions_object
    ns <- EntryWise$ns_object_plot
    plots <- list()
    for(i in 1:length(dists)){
      ns.df <- ns[[i]] # num switch
      entry <- sub_list[i,] # matrix entry
      dists.2 <- dists[[i]]$distribution # distributions
      # transform axis 2 to work with ggplot2
      # second y axis: max(dist_vals) = max(nsy)*x
      # where x is the fraction to get axis to match up
      dists.2$dist_vals.trans <- dists.2$dist_vals/
        (max(dists.2$dist_vals)/max(ns.df$nsy))
      p <- ggplot() + geom_step(data = ns.df, aes_string("nsx", "nsy"), size = 1.25, direction = 'vh') + 
        geom_area(data = dists.2, aes_string(x="x_range", y="dist_vals.trans"),
                  alpha = .2, color = "grey") + 
        theme_bw() +
        ggtitle(paste("(", entry[1], ", ", entry[2], ")", " entry", sep = '')) +
        #scale_x_continuous("Epsilon value") +
        #scale_y_continuous("Number of incorrect predictions",
                           #sec.axis = sec_axis(~ . * (max(dists.2$dist_vals)/max(ns.df$nsy)), 
                                               #name = "Probability density")) +
        scale_y_continuous(sec.axis = 
                             sec_axis(~ . * (max(dists.2$dist_vals)/max(ns.df$nsy)))) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              axis.line.y.right = element_line(color = "grey"), 
              axis.ticks.y.right = element_line(colour = "grey"),
              axis.text.y.right = element_text(color = "grey"),
              axis.title.y.right = element_text(color = "grey"))
      plots[[paste0(entry[1], ".", entry[2])]] <- p
    }
    filled.grid <- .grid_helper(n,m, plots)
    gridExtra::grid.arrange(
      grobs = plots,
      top = "Number of mis-predictions versus perturbation value,\noverlaid with distribution over stable perturbation values",
      left = "Number of incorrect predictions",
      right = grid::textGrob("Probability density", rot =270, gp=grid::gpar(col="grey")),
      bottom = "Epsilon value",
      layout_matrix = filled.grid
    )  
  }
}
