#' Select the clustering size for sparse tensor clustering via BIC
#' 
#' Select the clustering size for three-way clustering. The function searches over a range of clustering sizes and outputs the one that minimizes BIC. The clustering size (\eqn{d_1}, \eqn{d_2}, \eqn{d_3}) is a length-3 vector consisting of the number of clusters in each mode. 
#' @param x a three-dimensional array
#' @param k a vector, the possible numbers of clusters at mode 1
#' @param r a vector, the possible numbers of clusters at mode 2
#' @param l a vector, the possible numbers of clusters at mode 3
#' @param lambda a numeric value, regularization coefficient
#' @param sim.times the number of simulation replicates when performing clustering
#' @param method two options: "L0", "L1". "L0" indicates L0 penalty, and "L1" indicates Lasso penalty
#' @param n.cores the number of cores in parallel implementation
#' @return a list   
#' 
#' \code{estimated_krl} a 1*3 matrix consisting of the estimated clustering size  
#' 
#' \code{BIC} a vector consisting of the BIC value for all combinations of clustering sizes 
#' 
#' @export

choosekrl_bic = function (x,k,r,l,lambda=0,sim.times=1,method="L0",n.cores=NULL){
  #k = 2:5;r=2:5;l=2:5;lambda=0;sim.times=1;method="L0"
    ## x = x - mean(x) ## commented out
  if (sum(diff(k) <= 0) > 0 || sum(diff(r) <= 0) > 0 || sum(diff(l) <= 0) > 0) 
    stop("k and r has to be an increasing sequence.  Please sort k and r before using the function")
  n=dim(x)[1];p=dim(x)[2];q=dim(x)[3]
  krl = matrix(c(rep(1:length(k),each=length(r)*length(l)),
                 rep(1:length(r),times=length(k)*length(l)),
                 rep(rep(1:length(l),each=length(r)),times=length(k))),byrow=TRUE,
               nrow=3)
  krl_list = as.list(as.data.frame(krl))
  if (.Platform$OS.type == "windows") {
    bires = apply(krl,MARGIN=2,label_for_krl,k,r,l,sim.times=sim.times,lambda=lambda,xmiss=x,method=method,crossvalidation=FALSE)
    CBIC = lapply(bires,tensor_calculateBIC, x=x,method=method)
    BIC = unlist(CBIC)
  } else {
    bires = mclapply(krl_list, label_for_krl,k,r,l,sim.times=sim.times,lambda=lambda,xmiss=x,method=method,crossvalidation=FALSE,mc.cores=n.cores)
    CBIC = mclapply(bires,tensor_calculateBIC, x=x,method=method,mc.cores = n.cores)
    BIC = unlist(CBIC)
  }
  names(BIC) = apply(krl,MARGIN=2,FUN=function(x)paste(k[x[1]],r[x[2]],l[x[3]]))
  best = krl_list[[which(BIC == min(BIC))[1]]]
  return(list(estimated_krl = t(as.matrix(c(k[best[1]],r[best[2]],l[best[3]]))), BIC = BIC))
}
