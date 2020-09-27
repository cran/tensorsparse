#' Perform tensor clustering via tensor block model (TBM)
#' 
#' Perform tensor clustering via tensor block model (TBM) method.
#' @param x an order-3 data tensor
#' @param k an positive integer, the numbers of clusters at mode 1
#' @param r an positive integer, the numbers of clusters at mode 2
#' @param l an positive integer, the numbers of clusters at mode 3
#' @param lambda a numeric value, regularization coefficient
#' @param max.iter a positive integer, the maximum numbers of iteration
#' @param threshold a positive small numeric value for convergence threshold
#' @param sim.times the number of simulation replicates when performing clustering
#' @param trace logic value, print result per each iteration if TRUE
#' @param Cs.init vector or NULL, initial cluster label assignment at mode 1
#' @param Ds.init vector or NULL, initial cluster label assignment at mode 2
#' @param Es.init vector or NULL, initial cluster label assignment at mode 3
#' @param method two options: "L0", "L1". "L0" indicates L0 penalty, and "L1" indicates Lasso penalty
#' @return a list  
#' 
#' 
#' \code{judgeX} estimated underlying signal tensor
#' 
#'                \code{Cs} clustering result at mode 1
#'                
#'                \code{Ds} clustering result at mode 2
#'                
#'                \code{Es} clustering result at mode 3 
#'                
#'                \code{mus} estimated block means 
#'                
#' @export
#' @examples
#' x = getOrder3Tensor(20,20,20,2,2,2)$x
#' tbmClustering(x,2,2,2)
#' 
#' @references {M. Wang and Y. Zeng, "Multiway clustering via tensoe block models". Advances in Neural Information Processing System 32 (NeurIPS), 715-725, 2019.}
#' @author Yuchen Zeng \email{yzeng58@@wisc.edu}
#' 
tbmClustering = function(x,k,r,l,lambda=0,max.iter=1000,threshold = 1e-10,sim.times=1,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,method="L0"){
  #x=test;lambda=1e-3;max.iter=200;threshold = 5e-3;sim.times=10
  if (sim.times == 1) return(classify2(x,k,r,l,lambda=lambda,max.iter = max.iter,threshold = threshold,Cs.init = Cs.init,Ds.init = Ds.init,Es.init = Es.init,method=method))
  if (.Platform$OS.type == "windows") {
    result = lapply(rep(list(x),sim.times), classify2, k,r,l,lambda,max.iter,threshold,trace,Cs.init,Ds.init,Es.init,method=method)
    objs = unlist(lapply(result, function(result){result$objs}))
  } else {
    result = mclapply(rep(list(x),sim.times), classify2, k,r,l,lambda,max.iter,threshold,trace,Cs.init,Ds.init,Es.init,nstart = sample(1:1000,1),method=method,mc.cores = n.cores)
    objs = unlist(lapply(result, function(result){result$objs}))
  }
  result = result[[which(objs == min(objs))[1]]]
  return(result)
}
