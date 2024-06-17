#' Plot loadings from os_pca
#'
#'
#' @param output Output from os_pca
#' @param plot_vars Select which variables should be plotted by giving numeric indices
#' @param var_sel if TRUE select only large variables for loading plot (Variance Accounted For (VAF) + sd(VAF) > mean(VAF))
#' @param pldim number of dimensions for which the loadings should be calculated on the basis of the optimally scaled data
#'
#' @returns Loading plot from optimal scaling PCA object
#'
#' @export plot.os_loadings
plot.os_loadings<-function(output,plot_vars,var_sel,pldim){
if (missing(pldim)) pldim<-2
os_loadings<-output$os_loadings
os_scores<-output$os_scores
ndim<-dim(os_loadings)[2]
if (ndim == 1) {
   os_data<-output$os_data; nobj<-nrow(os_data); nvar<-ncol(os_data);os_svd<-svd(os_data,nu = pldim,nv = pldim)
   sval<-os_svd[[1]][1:pldim]; sval<-diag(sval/nobj**.5)
   os_loadings<-os_svd[[3]]%*%sval; varnames<-colnames(os_data); rownames(os_loadings)<-varnames
   os_scores<-os_svd[[2]]*nobj**.5

# reflect first dimension for loadings, object scores
# when largest value of loadings in first dimension is negative
   if (max (abs(os_loadings[,1])) > max(os_loadings[,1]) )
      {os_loadings[,1]<-os_loadings[,1]*-1
       os_scores[,1]<-os_scores[,1]*-1
   }
}
nvar<-dim(os_loadings)[1]
if (missing(plot_vars)) plot_vars<-1:nvar
if (missing(var_sel)) var_sel<-0
if (var_sel > 0) {VAF<-output$VAF;invar<-(VAF+var_sel*sd(VAF))>mean(VAF)
    plot_vars<-plot_vars[invar]}
par(mfrow=c(1,1),pty='s')
minx<-min(os_loadings[,1]*1.5);maxx<-max(os_loadings[,1]*1.5)
miny<-min(os_loadings[,2]*1.5);maxy<-max(os_loadings[,2]*1.5)
xlab="dimension 1"; ylab="dimension 2"; main="os component loadings"
#dev.new(width=10, height=10, unit="cm")
plot(rbind(os_loadings*1.05,os_loadings*-.50), type = "n", pch = 20,
xlab = xlab, ylab = ylab, main = main, cex = 0.5,asp=1,axes=FALSE)
r<-c(1:4)/4;axis(1,c(rev(r*-1),0,r)); axis(2,c(rev(r*-1),0,r),las=2);box()
lines(matrix(c(min(minx,miny),0,maxx,0),nrow=2,byrow=T),lty=2)
lines(matrix(c(0,min(minx,miny),0,maxx),nrow=2,byrow=T),lty=2)
for (i in plot_vars) {
    arrows(0, 0, os_loadings[i, 1], os_loadings[i, 2], length = 0.08)
    posvec <- apply(os_loadings, 1, sign)[2, ] + 2
    text(os_loadings[plot_vars,], labels = rownames(os_loadings)[plot_vars], pos = posvec, cex = 0.7)
    }
# points(os_scores[,1:2]) ; # perhaps later in biplot option
# if (ndim > 1) points(output$os_centroids)[,1:2]
}
# end plot loadings

#' Plot quantifications from os_pca
#'
#'
#' @param output Output from os_pca
#' @param plot_vars Select which variables should be plotted by giving numeric indices
#' @param rcplot select the number of rows and column to be displayed
#'
#' @returns Loading plot from optimal scaling PCA object
#'
#' @export plot.os_catquants
plot.os_catquants<-function(output,plot_vars,rcplot){
os_loadings<-output$os_loadings; Title<-rownames(os_loadings)
os_catquants<-output$os_catquants;c_ncat<-output$c_ncat
level<-output$level
#plot original categories and quantifications plots for selected variables
if (missing(plot_vars)) plot_vars<-1:(length(c_ncat)-1)
if (missing(rcplot)) rcplot <- c(3,ceiling(length(plot_vars)/3))
par(mfrow = rcplot);
for (j in plot_vars) {
            invar<-j;inj<-(c_ncat[invar]+1):c_ncat[invar+1];
            pl1<-cbind(as.vector(c(1:length(inj))),os_catquants[inj,])
            title<-c(Title[j],level[j])
            plot(type='n',pl1,xlab="original categories",ylab="optimal quantifications",main=title,axes=FALSE);
            clm<-ceiling(max(abs(os_catquants[inj,])));r<-c(1:(2*clm))/2
            axis(1,1:length(inj));axis(2,c(rev(r*-1),0,r),las=2);box()
            lines(type='s',pl1,col='black',lty=1,lwd=2)
            points(pl1,col='black',lty=1,lwd=2)
            }
}
#end plot catquants --------------------------------------------

# Scree plot
#' Scree plot from os_pca
#'
#'
#' @param output Output from os_pca
#'
#' @returns Scree plot from optimal scaling PCA object
#'
#' @export plot.eigval
# Scree plot
plot.eigval<-function(output) {
Evalues<-output[[1]]
ndim<-output$ndim
par(mfrow=c(1,1),pty='m')
GV<-output$os_data
Evalues<-eigen(cor(GV));Evalues<-Evalues[[1]]
pl1<-Evalues
title<-"Eigenvalues of correlation matrix in decreasing order"
plot(type='n',pl1,xlab="number",ylab="eigenvalues",main=title,axes=FALSE);
clm<-ceiling(max(pl1));r<-c(1:(2*clm))/2
axis(1,1:length(pl1));axis(2,c(rev(r*-1),0,r),las=2);box()
lines(type='S',pl1,col='black',lty=1,lwd=2)
points(pl1,col='black',lty=1,lwd=2)
lines(matrix(c(0,pl1[ndim+1],length(pl1)+1,pl1[ndim+1]),nrow=2,byrow=T),lty=2)
}

#end screeplot-------------------------------------------

#' Ordering plot
#'
#'
#' @param output Output from os_pca
#' @param plot_var Select which variable should be plotted by giving numeric index
#'
#' @returns plot of the optimal ordering of the categories the selected variable
#'
#' @export plot.ordering
plot.ordering<-function(output,plot_var){
c_ncat<-output$c_ncat;os_catquants<-output$os_catquants
invar<-plot_var;inj<-(c_ncat[invar]+1):c_ncat[invar+1];
pl1<-cbind(as.vector(c(1:length(inj))),os_catquants[inj,])
Title<-rownames(output$os_loadings);title<-c("transformation",Title[invar])
par(mfrow=c(1,2))
plot(type='n',pl1,xlab="original categories",ylab="optimal quantifications",main=title,axes=FALSE);
clm<-ceiling(max(abs(os_catquants[inj,])));r<-c(1:(2*clm))/2
axis(1,1:length(inj));axis(2,c(rev(r*-1),0,r),las=2);box()
lines(type='s',pl1,col='black',lty=1,lwd=2); points(pl1,col='black',lty=1,lwd=2)

# ordering of categories using optimal quantifications
title<-c("optimal ordering",Title[invar])
pl1[,2]<-pl1[or<-order(pl1[,2]),2]
plot(type='n',pl1,xlab="optimally ordered categories",ylab="optimal quantifications",main=title,axes=FALSE);
clm<-ceiling(max(abs(os_catquants[inj,])));r<-c(1:(2*clm))/2
axis(1,c(1:length(inj)), label=or)
axis(2,c(rev(r*-1),0,r),las=2);box()
lines(type='s',pl1,col='black',lty=1,lwd=2); points(pl1,col='black',lty=1,lwd=2)
}
