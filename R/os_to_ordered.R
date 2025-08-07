#' Output OSPCA to ordered factors
#'
#' translate the results of the os_pca function to ordered factors, to use as
#' input for `rvinecopulib`'s `vine` function
#'
#' @param os_object Object with os_scores, os_loadings, os_data, os_catquants
#' from the `os_pca` function
#' @param data A data.frame input data for os_pca, if ommitted data from os_pca
#' function can be used (only when `keep_data` is TRUE)
#'
#'
#' @returns Object with oredered categorical data (cat_data) and the levels of
#' those ordered factors
#'
#' @export os_to_ordered
os_to_ordered <- function(os_object, data) {
  if (missing(data)) {
    if (is.null(os_object$data)) {
      stop("Add data to data object!")
    }
    data <- os_object$data
  }
  cat_data <- data
  os_cols <- names(cat_data)[os_object$level != "numerical"]
  cat_data[, os_cols] <- os_object$os_data[, os_cols]
  trans_list <- list()
  for (column in os_cols) {
    t_df <- unique(data.frame(os = cat_data[, column], orig = data[, column]))
    levels_os <- as.character(t_df[order(t_df$os), "orig"])
    trans_list <- c(trans_list, list(levels_os))
    cat_data[, column] <- factor(cat_data[, column], ordered = TRUE)
    levels(cat_data[, column]) <- levels_os
  }
  names(trans_list) <- os_cols

  attr(cat_data, "os_levels") <- trans_list

  return(cat_data)
}
