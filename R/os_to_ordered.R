#' Output OSPCA to ordered factors
#'
#' translate the results of the os_pca function to ordered factors, to use as
#' input for `rvinecopulib`'s `vine` function
#'
#' @param os_object Object with os_scores, os_loadings, os_data, os_catquants
#' from the `os_pca` function
#' @param data A data.frame input data for os_pca, if omitted data from os_pca
#' function can be used (only when `keep_data` is TRUE)
#' @param trans_list A named list of ordered character vectors (as returned in
#' `attr(cat_data, "os_levels")`), one per categorical column. If provided,
#' `os_object$data` can be NULL.
#'
#' @returns Object with ordered categorical data (cat_data) and the levels of
#' those ordered factors
#'
#' @export os_to_ordered
os_to_ordered <- function (os_object, data, trans_list = NULL)
{
  if (missing(data)) {
    if (is.null(os_object$data)) {
      stop("Add original data to os_object or provide data to function")
    }
    data <- os_object$data
    message("Changing original data to ordered factor")
  }
  os_cols <- rownames(os_object$os_loadings)[os_object$level !=
                                               "numerical"]
  if (is.null(trans_list)) {
    if (is.null(os_object$data)) {
      stop("Add original data to os_object or add provide trans_list to function")
    }
    trans_list <- list()
    for (column in os_cols) {
      orig_data <- os_object$data[, column]
      os_data <- os_object$os_data[, column]
      lookup <- unique(data.frame(orig = orig_data, os = os_data))
      if (os_object$level[which(rownames(os_object$os_loadings) == column)] == "ordinal") {
        lookup <- merge(data.frame(orig = ordered(levels(orig_data), levels(orig_data))),
                        lookup, all.x = TRUE)
        lookup$os <- approxfun(lookup$orig, lookup$os, rule = 2)(lookup$orig)
      }
      trans_list[[column]] <- as.character(lookup$orig[order(lookup$os)])
    }
  }
  cat_data <- data
  for (column in os_cols) {
    levels_os <- trans_list[[column]]
    if (anyNA(match(data[, column], levels_os))) {
      unseen <- unique(data[is.na(match(data[, column],
                                        levels_os)), column])
      warning(sprintf("Column '%s': %d value(s) not found in levels and set to NA: %s",
                      column, length(unseen), paste(unseen, collapse = ", ")))
    }
    cat_data[, column] <- factor(data[, column], levels = levels_os,
                                 ordered = TRUE)
  }
  attr(cat_data, "os_levels") <- trans_list
  return(cat_data)
}
