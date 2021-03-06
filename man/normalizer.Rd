% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalizer.R
\name{normalizer}
\alias{normalizer}
\title{Function to normalize distributions}
\usage{
normalizer(x, plotit = FALSE, return_all = FALSE, ...)
}
\arguments{
\item{x}{a numeric vector without absent data.}

\item{plotit}{Logical If TRUE histograms of each transformation will be ploted.}

\item{return_all}{Logical. If TRUE it will return all values.
Default is FALSE, and only the selected vector is returned.}

\item{...}{Other parameters passed to the \code{hist} function.}
}
\value{
The function will return a list with transformed vector
(or the raw data) that is closer to a normal distribution and the D value.
If \code{return_all} equals TRUE, than a list with all transformations
is returned.
}
\description{
The function will apply several transformations,
and indicate the transformation that best approximates a normal distribution.
}
\details{
The function will apply eigth different transformations to the data:
\enumerate{
\item Square root
\item Natural logarithm
\item Cube root
\item Box Cox
\item Logarithm in the base 10
\item Inverse (1 / x)
\item Inverse Natural logarithm (revert the data, apply a natural logarithm
and then invert it back)
\item Tukey Ladder of Powers (if the vector is between 3 and 5000 long)
}
The \code{normalizer} function uses the standardized D value calculated by a
Kolmogorov-Smirnov test using the function \code{stats::ks.test} to determine
the transformation that most approximate a normal distribution.
}
\examples{
# Example 1:
library(letsR)
x <- as.vector(na.exclude(values(temp)))
x <- x / 100
x_t <- normalizer(x, plotit = TRUE)

# Example 2:
library(datasets)
data(mtcars)
mpg_t <- normalizer(mtcars$mpg, plotit = TRUE, return_all = TRUE)
}
\author{
Bruno Vilela (email: \email{bvilela@wustl.edu})
}
