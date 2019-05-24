############ Helper functions

#' Get PDF distribution
#'
#' This function retrieves the PDF (Probablity Distribution Function)
#' object from the scipy method 
#' <scipy.stats._distn_infrastructure.rv_frozen>.
#' @param matrix_interval Position in the matrix. Example: c(0, 0)
#' @param distribution_list list of scipy distribution
#' @param asymp_stab asymptotic stability interval
#' @param points
#' @export
#' @examples ax2 <- get_distributions(matrix_interval = c(0, 0), distribution_list, asymp_stab)
#' @import reticulate

get_distributions <- function(matrix_interval, 
                              distribution_list, asymp_stab,
                              points=250){
  np <- reticulate::import("numpy")
  k <- matrix_interval[1]
  l <- matrix_interval[2]
  interval <- asymp_stab[(k+1), (l+1),] # change back to R 1 based
  padding <- (interval[2] - interval[1])/100
  x_range <- np$linspace((interval[1] - padding), (interval[2] + padding), points)
  dist.py <- distribution_list[paste("(", k, ", ", l, ")", sep = '')]
  dist_vals <- sapply(x_range, function(x){dist.py[[1]]$pdf(x)})
  ax2 <- data.frame(x_range = as.numeric(x_range), dist_vals = dist_vals)
  return(ax2)
}

# new version of get_distribution?

get_distributions_single <- function(matrix_interval, 
                                     distribution_list, asymp_stab,
                                     points=250){
  np <- reticulate::import("numpy")
  k <- matrix_interval[1]
  l <- matrix_interval[2]
  interval <- asymp_stab
  padding <- (interval[2] - interval[1])/100
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
