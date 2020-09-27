#' Select the regularization coefficient for order-3 sparse tensor clustering via BIC
#' 
#' Select the regularization coefficient for three-way clustering. The clustering size is assumed to be known. The function searches over a range of regularization sizes and outputs the one that minimizes the BIC.
#' @param x a three-dimensional array
#' @param k an positive integer, the numbers of clusters at mode 1
#' @param r an positive integer, the numbers of clusters at mode 2
#' @param l an positive integer, the numbers of clusters at mode 3
#' @param lambda a vector of possible lambda, eg: lambda = c(0,50,100,200)
#' @param method two options: "L0", "L1". "L0" indicates L0 penalty, and "L1" indicates Lasso penalty
#' @return a list   
#' 
#' 
#'                \code{lambda} the lambda with lowest BIC
#' 
#'                \code{BIC} the BIC for each lambda in the given range
#'                
#'                \code{nonzeromus} the number of clusters with non-zero means
#'   
#' @export
chooseLambda = function (x, k, r, l, lambda=NULL,method="L0") {
  ##x = x - mean(x) commented out
  if (is.null(lambda)){
    n=dim(x)[1];p=dim(x)[2];q=dim(x)[3]
    if (method == "L0") lambda = sqrt((n*p*q)/(k*r*l))*seq(0,2,by=0.1)
    if (method == "L1") lambda = (n*p*q)/(k*r*l)*seq(0,2,by=0.1)
    if (is.null(lambda)) stop("No such kind of method:", method, ".\n")
  } 
  if (.Platform$OS.type == "windows") {
    bires = lapply(lambda,FUN=classify2,x=x,k=k,r=r,l=l,method=method)
    CBIC = lapply(bires,tensor_calculateBIC, x=x,method=method)
    BIC = unlist(CBIC)
    nonzero = unlist(lapply(bires, FUN=function(bires){return(sum(bires$mus!=0))}))
  } else {
    bires = mclapply(lambda,FUN=classify2,x=x,k=k,r=r,l=l,method=method,mc.cores = n.cores)
    CBIC = mclapply(bires,tensor_calculateBIC, x=x,method=method,mc.cores = n.cores)
    BIC = unlist(CBIC)
    nonzero = unlist(mclapply(bires, FUN=function(bires){return(sum(bires$mus!=0))},mc.cores = n.cores))
  }
  return(list(lambda = lambda[which(BIC == min(BIC))[1]], BIC = BIC, 
              nonzeromus = nonzero))
}
