% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmi510chengding.R
\name{rep_mat}
\alias{rep_mat}
\title{Replicate a matrix or data frame by repeating rows or columns}
\usage{
rep_mat(x, M = 1, N = 1)
}
\arguments{
\item{x}{A matrix or data frame to replicate}

\item{M}{The number of times to repeat the rows of \code{x} (default = 1)}

\item{N}{The number of times to repeat the columns of \code{x} (default = 1)}
}
\value{
A matrix with rows and/or columns replicated M and/or N times
}
\description{
This function replicates a matrix or data frame \code{x} by repeating its rows \code{M} times or its columns \code{N} times, depending on the parameters passed.
}
