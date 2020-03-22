#' Generate a random order-4 tensor
#' 
#' Generate a random order-4 tensor based on tensor block model.
#' @param n the dimension at mode 1
#' @param p the dimension at mode 2
#' @param q the dimension at mode 3
#' @param s the dimension at mode 4
#' @param k an positive integer, the numbers of clusters at mode 1
#' @param r an positive integer, the numbers of clusters at mode 2
#' @param l an positive integer, the numbers of clusters at mode 3
#' @param m an positive integer, the numbers of clusters at mode 4
#' @param error a positive numeric value, noise level
#' @param sort if TRUE, the tensor entries belonging to the same cluster would be assumed together
#' @param sparse.percent the proportion of zero entries based on the Gaussian tensor block model
#' @param center if True, the data tensor would be centered to zero-mean before clustering
#' @param seed a positive integer, used to specify the random seed
#' @param mumin a numeric value, the lower bound of the block mean
#' @param mumax a numeric value, the upper bound of the block mean
#' @return a list 
#' 
#'                \code{x} the tensor   
#' 
#'                \code{truthX} the underlying signal tensor following block model
#'                
#'                \code{truthCs} true cluster label assignment at mode 1   
#'                
#'                \code{truthDs} true cluster label assignment at mode 2   
#'                
#'                \code{truthEs} true cluster label assignment at mode 3 
#'                
#'                \code{truthFs} true cluster label assignment at mode 4
#'                
#'                \code{mus} the block means
#'                
#'                \code{binaryX} the 0-1 tensor (0:the mean signal = 0; 1:the mean signal != 0)   
#'            
#' @export
#' 
#' get_data4(10,10,10,10,2,2,2,2)
get_data4 = function(n,p,q,s,k=NULL,r=NULL,l=NULL,m=NULL,error=3,sort=TRUE,sparse.percent=0,center=FALSE,seed=NULL,mumin = -3, mumax = 3){
  if(!is.null(seed)) set.seed(seed)
  #print(m)
  mus = runif(k*r*l*m,mumin,mumax)#take the mean of k*r*l biclusters/cubes
  if(sparse.percent!=0) mus[sample(k*r*l*m,floor(k*r*l*m*sparse.percent),replace=F)]= 0
  mus = array(mus,c(k,r,l,m))
  
  if(is.null(mus)) stop("multiplicative must be a positive integer!") 
  if (k!=1) truthCs = ReNumber(sample(1:k,n,replace=TRUE)) else truthCs = rep(1,n)
  if (r!=1) truthDs = ReNumber(sample(1:r,p,replace=TRUE)) else truthDs = rep(1,p)
  if (l!=1) truthEs = ReNumber(sample(1:l,q,replace=TRUE)) else truthEs = rep(1,q)
  if (m!=1) truthFs = ReNumber(sample(1:m,s,replace=TRUE)) else truthFs = rep(1,s)
  
  ##### added
  if(sort==TRUE){
    truthCs=sort(truthCs)
    truthDs=sort(truthDs)
    truthEs=sort(truthEs)
    truthFs=sort(truthFs)
  }
  ######
  
  x = array(rnorm(n*p*q*s,mean=0,sd=error),dim = c(n,p,q,s))
  truthX = array(rep(0,n*p*q*s),c(n,p,q,s))
  for(i1 in 1:max(truthCs)){
    for(i2 in 1:max(truthDs)){
      for(i3 in 1:max(truthEs)){
        for(i4 in 1:max(truthFs)){
          x[truthCs==i1, truthDs==i2, truthEs==i3, truthFs==i4] = x[truthCs==i1, truthDs==i2, truthEs==i3, truthFs==i4] + mus[i1,i2,i3,i4]
          truthX[truthCs==i1, truthDs==i2, truthEs==i3, truthFs==i4] =  mus[i1,i2,i3,i4]
        }
      }
    }
  }
  if (center == TRUE) x = x - mean(x)
  binaryX = (truthX!=0)*1
  result = list("x"=x,"truthX"=truthX,"truthCs"=truthCs,"truthDs"=truthDs,"truthEs"=truthEs,"truthFs"=truthFs,"mus"=mus,"binaryX"=binaryX)
  #}
  return(result)
}
