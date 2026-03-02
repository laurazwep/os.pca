# OS_PCA, based on homals in R. Errors in ordinal optimal scaling have been corrected.
# eigenvalues, centroids, component loadings and component scores are appropriately rescaled
# matrices for homals results have been corrected
# output: os_scores, os_loadings, os_data, os_catquants

#' Optimal Scaling PCA
#'
#' OS_PCA, based on homals in R. Errors in ordinal optimal scaling have been corrected.
#' eigenvalues, centroids, component loadings and component scores are appropriately rescaled
#' matrices for homals results have been corrected
#' output: os_scores, os_loadings, os_data, os_catquants
#'
#' @param data A data.frame input data for os_pca
#' @param level Which quantification levels. Possible values are `"nominal"`, `"ordinal"`, and `"numerical"` which can be defined as single character (if all variable are of the same level) or as vector which length corresponds to the number of variables.
#' @param ndim Number of dimensions to be extracted.
#' @param reflec1 Boolean, if TRUE, let the largest loading point to positive side in the first dimension
#' @param reflec2 Boolean, if TRUE, the direction of the original variable is preserved in the transformation function
#' @param homals_only Output only homals output or also other objects
#'
#' @returns Object with os_scores, os_loadings, os_data, os_catquants
#'
#' @import homals
#' @export os_pca
os_pca <- function(data, level, ndim = 2, reflec1 = 1, reflec2 = 1, homals_only = FALSE, as.is = TRUE, keep_data = FALSE) {
  data_b <- recode(data)
  if (missing(level)) {
    level <- rep("nominal", ncol(data))
  }
  #scale using the homals algorithm
  output <- homals(data_b, ndim = ndim, rank = 1, level = level, sets = 0,
                   active = TRUE, eps = 1E-8, itermax = 10000, verbose = 0)
  if (homals_only) {
    return(output)
  }

  nvar <- dim(data)[2]
  n <- dim(data)[1]

  # compute the number of categosies (ncat) and cumulative ncat
  ncat <- numeric(length = nvar)
  for (j in 1:nvar) {
    ncat[j] <- length(unique(data[, j]))
  }
  c_ncat <- c(0, cumsum(ncat))

  # rescale eigenvalues and loadings; eigenvalues should equal sum of squares of the loadings;
  # loadings should be correlations between transformed variables and principal components
  # (colSums((loadings*sqrt(nvar))^2)) = SSQ(loadings) from homals * sqrt(nvar); number of variables = nvar
  # principal components are called objscores; rescale to variance=1
  # construct matrix with loadings

  os_loadings <- sqrt(nvar)*matrix(unlist(output$loadings), ncol = ndim, byrow = TRUE)
  rownames(os_loadings) <- colnames(data)

  os_scores <- sqrt(n*nvar) * output$objscores
  os_centroids <- sqrt(n*nvar) * do.call("rbind", output$cat.centroids)
  catscores <- sqrt(n*nvar) * do.call("rbind", output$catscores)
  os_catquants <- sqrt(n) * matrix(unlist(output$low.rank))

  #catquants per category
  os_cat_list <- lapply(output$low.rank, function(x) sqrt(n)*x)
  orig_values <- lapply(output$low.rank, function(x) row.names(x))

  #calculate eigenvalues and variance accounted for
  EV1 <- colSums(os_loadings^2)
  VAF <- rowSums(os_loadings^2)

  # reflect first dimension for loadings, object scores and centroids
  # when largest value of loadings in first dimension is negative
  if (reflec1) {
    if (max(abs(os_loadings[, 1])) > max(os_loadings[, 1])) {
      os_loadings[, 1] <- os_loadings[, 1] * -1
      os_scores[, 1] <- os_scores[, 1] * -1
      os_centroids[, 1] <- os_centroids[, 1] * -1
    }
  }

  # transformed data matrix; for homals, multiple transformed data matrices!
  # we don't use these. Normalization probably also wrong. Order also wrong?

  # Form transformed data matrix on the basis of indicator matrices and category
  # quantifications, possibly with reversed signs
  num_data <- data
  num_data[] <- lapply(data, function(x) as.numeric(as.factor(as.numeric(x))))

  # transformed data matrix
  GV <- data_b
  for (j in 1:nvar) {
    map <- setNames(os_cat_list[[j]], orig_values[[j]])
    GV[, j] <- map[GV[, j]]
  }
  # GV is the transformed data matrix
  os_data <- GV
  cc1 <- suppressWarnings(cor(cbind(num_data, os_data), use = "pairwise.complete.obs"))

  # check whether transformation is predominantly increasing or decreasing.
  # this is related to direction of vectors for loadings in plot
  # correlation between original and transformed should be positive. If not, reverse the signs. Also for loadings.
  # compute new transformed data matrix
  if (reflec2 == 1) {
    cor_quant_orig <- diag(cc1[1:nvar, (nvar + 1):(nvar*2)])
    if (any(cor_quant_orig < 0)) {
      col_flipped <- which(cor_quant_orig < 0)
      for (j in col_flipped) {
        inj <- (1 + c_ncat[j]):c_ncat[j + 1]
        if (cc1[j, nvar + j] < 0) {
          os_catquants[inj, ] <- os_catquants[inj, ] * -1
          os_loadings[j, ] <- os_loadings[j, ] * -1
          os_cat_list[[j]] <- os_cat_list[[j]] * -1
        }
      }
      GV <- data_b
      for (j in 1:nvar) {
        map <- setNames(os_cat_list[[j]], orig_values[[j]])
        GV[, j] <- map[GV[, j]]
      }
      os_data <- GV
      cc1 <- suppressWarnings(cor(cbind(num_data, os_data), use = "pairwise.complete.obs"))
    }
  }
  cor_matrix <- cc1[1:nvar, (nvar + 1):(nvar*2)]

  results <- list(
    eigenvalues = EV1, VAF = VAF, os_catquants = as.data.frame(os_catquants),
    os_loadings = as.data.frame(os_loadings), os_scores = as.data.frame(os_scores),
    os_centroids = as.data.frame(os_centroids), os_data = os_data,
    c_ncat = c_ncat, level = level, ndim = ndim, cc1 = cor_matrix
  )
  if (keep_data) {
    results$data <- data
  }

  return(results)
}
