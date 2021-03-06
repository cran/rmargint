% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmargint-fn.R
\name{kernel6}
\alias{kernel6}
\title{Order 6 kernel}
\usage{
kernel6(x)
}
\arguments{
\item{x}{A vector of real numbers.}
}
\value{
A vector of the same length as \code{x} where each entry is \code{( 105/256 ) * ( 1 - x^2 ) * ( 5 - 30 * x^2 + 33 * x^4 )} if \code{abs(x) < 1} and 0 otherwise.
}
\description{
This function evaluates a kernel of order 6.
}
\details{
This function evaluates a kernel of order 6. A kernel L is a kernel of order 6 if it integrates 1, the integrals of u^j L(u) are 0 for 1 <= j < 6 (j integer) and the integral of u^6 L(u) is different from 0.
}
\examples{
x <- seq(-2,2,length=10)
kernel6(x)

}
\author{
Alejandra Martinez, \email{ale_m_martinez@hotmail.com}, Matias Salibian-Barrera
}
