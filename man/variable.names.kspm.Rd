% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variableNamesKspm.R
\name{variable.names.kspm}
\alias{variable.names.kspm}
\title{Variable names of fitted models}
\usage{
\method{variable.names}{kspm}(object, ...)
}
\arguments{
\item{object}{an object of class "kspm", usually, a result of a call to \code{kspm}.}

\item{...}{additional optional argument (currently unused).}
}
\value{
a list of character vectors. The first element correspond to the names of variables included in the linear part of the model. Then, a vector containing names of variables including in kernel part is provided for each kernel.
}
\description{
Simple utility returning names of variables involved in a kernel semi parametric model.
}
\seealso{
\link{kspm}, \link{summary.kspm}, \link{case.names.kspm}.
}
\author{
Catherine Schramm, Aurelie Labbe, Celia Greenwood
}
