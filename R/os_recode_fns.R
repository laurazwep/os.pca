#' @noRd
del_missing<-function(dat) {
d<-dat[1,]*0
for (i in 1:nrow(dat)) {
    if (min(dat[i,])!=0) d<-rbind(d,dat[i,])
    }
    d<-(d[2:nrow(d),])
    return(d)
 }

#' @noRd
recode <- function(data) {
  for (j in 1:ncol(data)) {
    var_j <- as.factor(data[, j])
    levels(var_j) <- sort(levels(var_j))
    data[, j] <- var_j
  }
  return(data)
}
