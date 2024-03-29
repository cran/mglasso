% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_clusterpath}
\alias{plot_clusterpath}
\title{Plot MGLasso Clusterpath}
\usage{
plot_clusterpath(X, mglasso_res, colnames_ = NULL)
}
\arguments{
\item{X}{numeric matrix}

\item{mglasso_res}{object of class \code{mglasso}}

\item{colnames_}{columns labels}
}
\value{
no return value.
}
\description{
Plot MGLasso Clusterpath
}
\details{
This function plot the clustering path of mglasso method on the 2 principal components
axis of X. As the centroids matrices are not of the same dimension as X, we choose to plot the
predicted X matrix path.
}
