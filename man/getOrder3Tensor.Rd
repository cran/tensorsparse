% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getOrder3Tensor.R
\name{getOrder3Tensor}
\alias{getOrder3Tensor}
\title{Generate a random order-3 tensor}
\usage{
getOrder3Tensor(
  n,
  p,
  q,
  k = NULL,
  r = NULL,
  l = NULL,
  error = 3,
  sort = TRUE,
  sparse.percent = 0,
  center = FALSE,
  seed = NULL,
  mumin = -3,
  mumax = 3
)
}
\arguments{
\item{n}{the dimension at mode 1}

\item{p}{the dimension at mode 2}

\item{q}{the dimension at mode 3}

\item{k}{an positive integer, the numbers of clusters at mode 1}

\item{r}{an positive integer, the numbers of clusters at mode 2}

\item{l}{an positive integer, the numbers of clusters at mode 3}

\item{error}{a positive numeric value, noise level}

\item{sort}{if TRUE, the tensor entries belonging to the same cluster would be assumed together}

\item{sparse.percent}{the proportion of zero entries based on the Gaussian tensor block model}

\item{center}{if True, the data tensor would be centered to zero-mean before clustering}

\item{seed}{a positive integer, used to specify the random seed}

\item{mumin}{a numeric value, the lower bound of the block mean}

\item{mumax}{a numeric value, the upper bound of the block mean}
}
\value{
a list 

               \code{x} the tensor   
  
               \code{truthX} the underlying signal tensor following block model
               
               \code{truthCs} true cluster label assignment at mode 1   
               
               \code{truthDs} true cluster label assignment at mode 2   
               
               \code{truthEs} true cluster label assignment at mode 3   
               
               \code{mus} the block means
               
               \code{binaryX} the 0-1 tensor (0:the mean signal = 0; 1:the mean signal != 0)
}
\description{
Generate an order-3 random tensor based on tensor block model.
}
\examples{

getOrder3Tensor(20,20,20,2,2,2)$x
}
