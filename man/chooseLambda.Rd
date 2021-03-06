% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chooseLambda.R
\name{chooseLambda}
\alias{chooseLambda}
\title{Select the regularization coefficient for order-3 sparse tensor clustering via BIC}
\usage{
chooseLambda(x, k, r, l, lambda = NULL, method = "L0")
}
\arguments{
\item{x}{a three-dimensional array}

\item{k}{an positive integer, the numbers of clusters at mode 1}

\item{r}{an positive integer, the numbers of clusters at mode 2}

\item{l}{an positive integer, the numbers of clusters at mode 3}

\item{lambda}{a vector of possible lambda, eg: lambda = c(0,50,100,200)}

\item{method}{two options: "L0", "L1". "L0" indicates L0 penalty, and "L1" indicates Lasso penalty}
}
\value{
a list   


               \code{lambda} the lambda with lowest BIC

               \code{BIC} the BIC for each lambda in the given range
               
               \code{nonzeromus} the number of clusters with non-zero means
}
\description{
Select the regularization coefficient for three-way clustering. The clustering size is assumed to be known. The function searches over a range of regularization sizes and outputs the one that minimizes the BIC.
}
