% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmargint-fn.R
\name{psi.tukey}
\alias{psi.tukey}
\title{Derivative of Tukey's bi-square loss function.}
\usage{
psi.tukey(r, k = 4.685)
}
\arguments{
\item{r}{A vector of real numbers}

\item{k}{A positive tuning constant.}
}
\value{
A vector of the same length as \code{r}.
}
\description{
This function evaluates the first derivative of Tukey's bi-square loss function.
}
\details{
This function evaluates the first derivative of Tukey's bi-square loss function.
}
\examples{
x <- seq(-2, 2, length=10)
psi.tukey(r=x, k = 1.5)

}
\author{
Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
}
