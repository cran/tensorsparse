% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chooseLambda.R
\name{chooseLambda}
\alias{chooseLambda}
\title{Perform tuning parameter (lambda) selection for sparse tensor clustering via BIC criterion}
\usage{
chooseLambda(x, k, r, l, lambda = NULL, method = "L0")
}
\arguments{
\item{x}{a three-dimensional array}

\item{k}{\eqn{d_1}: the clusters number of mode 1}

\item{r}{\eqn{d_2}: the clusters number of mode 2}

\item{l}{\eqn{d_3}: the clusters number of mode 3}

\item{lambda}{a vector of possible lambda, eg: lambda = c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200).}

\item{method}{two options: "L0", "L1". Two methods use different penalties, where "L0" indicating L0 penalty and "L1" indicating Lasso penalty.}
}
\value{
a list   
               \code{lambda} the lambda with lowest BIC;  

               \code{BIC} the corresponding BIC for each lambda in the given range;  
               
               \code{nonzeromus} the number clusters with non-zero mean.
}
\description{
We assume that \eqn{d_1}, \eqn{d_2}, \eqn{d_3} are known. A range of values of lambda is usually considered - value that results in the lowest BIC is selected.
}