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
#' @param output_only Output only homals output or also other objects
#'
#' @returns Object with os_scores, os_loadings, os_data, os_catquants
#'
#' @import homals
#' @export

os_pca <- function(data, level, ndim, output_only = FALSE) {
  data_b<-recode(data)
  if (missing(level))  level<-rep("nominal",ncol(data))
  if (missing(ndim)) ndim<-2
  output <- homals::homals(data_b, ndim = ndim, rank = 1, level = level, sets = 0,
                active = TRUE, eps = 1E-8, itermax = 10000, verbose = 0)
  if (output_only) {
    return(output)
  }
nvar<-dim(data)[2];n<-dim(data)[1]
# compute ncat and cumulative ncat
ncat <- numeric(length = nvar)
  for (j in 1:nvar) {
    ncat[j] <- length(unique(data[, j]))
  }
  c_ncat <- c(0, cumsum(ncat))

# rescale eigenvalues and loadings; eigenvalues should equal sum of squares of the loadings;
#loadings should be correlations between transformed variables and principal components
#(colSums((loadings*(nvar**.5))**2)) = SSQ(loadings) from homals * nvar**.5; number of variables = nvar
#principal components are called objscores; rescale to variance=1
# construct matrix with loadings


loadings<-matrix(unlist(output$loadings), ncol=ndim,byrow=TRUE);
var_names<-colnames(data)
rownames(loadings)<-var_names
#loadings <- t(sapply(output$loadings, function(xy) xy[1, c(1,2)]))
os_loadings<-loadings*(nvar**.5)


objscores<-output$objscores; n<-(dim(objscores))[1]
os_scores<-objscores*((n*nvar)**.5)
os_centroids<-((n*nvar)**.5)*matrix(unlist(output$cat.centroids),ncol=ndim,byrow=TRUE)
catscores<-((n*nvar)**.5)*matrix(unlist(output$catscores),ncol=ndim,byrow=TRUE)
os_catquants<-(n**.5)*matrix(unlist(output$low.rank),ncol=1,byrow=TRUE)

#matrices for centroids and catscores are constructed wrongly. This error is corrected here.
for (j in 1:nvar) {
                  invar<-j;inj<-(c_ncat[invar]+1):c_ncat[invar+1];
                  os_centroids[inj,]<-matrix(as.vector(t(os_centroids[inj,])),nrow=ncat[j])
                  catscores[inj,]<-matrix(as.vector(t(catscores[inj,])),nrow=ncat[j])
                  }


EV1<-colSums(os_loadings**2)
VAF<-rowSums(os_loadings**2)
#Check eigenvalues:
EV2<-(2*nvar**2)*(output$eigenvalues)
# still to check: for 3 or more dimensions. Same factor? or 3*?
# colSums(os_scores**2);colSums(os_loadings**2)

# reflect first dimension for loadings, object scores and centroids if necessary.
# for example, when largest value of loadings in first dimension is negative
if (max (abs(os_loadings[,1])) > max(os_loadings[,1]) )
   {os_loadings[,1]<-os_loadings[,1]*-1
    os_scores[,1]<-os_scores[,1]*-1
    os_centroids[,1]<-os_centroids[,1]*-1
   }

#GV1<-output$scoremat[,,1];GV2<-output$scoremat[,,2]
#transformed data matrix; for homals, multiple transformed data matrices!
# we don't use these. Normalization probably also wrong. Order also wrong?


# check whether transformation is predominantly increasing or decreasing.
# this is related to direction of vectors for loadings in plot
# first category quantificaton should be smaller than last one. If not, reverse the signs. Also for loadings.
for (j in 1:nvar) {inj<-(1+c_ncat[j]):c_ncat[j+1];
                  if (mean(os_catquants[(1:3)+c_ncat[j],]) > os_catquants[c_ncat[j+1],])
                     {os_catquants[inj,]<-os_catquants[inj,]*-1;
                    os_loadings[j,]<-os_loadings[j,]*-1}
                  }

# Form transformed data matrix on the basis of indicator matrices and category quantifications, possibly with reversed signs
Vtot<-os_catquants;Gtot<-output$ind.mat;GV<-data*0;
for (j in 1:nvar) {inj<-((1+c_ncat[j]):c_ncat[j+1]);GV[,j]<-Gtot[,inj]%*%Vtot[inj,]}
# GV is the transformed data matrix
os_data<-GV
results <- list(eigenvalues = EV1, VAF = VAF, os_catquants = as.data.frame(os_catquants), os_loadings = as.data.frame(os_loadings), os_scores = as.data.frame(os_scores),
                 os_data = os_data,c_ncat = c_ncat)

  return(results)
}
# end function -------------------------------------------------------
















