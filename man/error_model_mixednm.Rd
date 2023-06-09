% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residual_error.R
\name{error_model_mixednm}
\alias{error_model_mixednm}
\title{Residual error model mixed (idem NONMEM)}
\usage{
error_model_mixednm(f, sigma)
}
\arguments{
\item{f}{Numeric vector, output of a pharmacokinetic model}

\item{sigma}{Matrix of the coefficients for the
residual error model}
}
\value{
Numeric vector, residual error
}
\description{
Mixed residual error model, similar to NONMEM implementation.
}
