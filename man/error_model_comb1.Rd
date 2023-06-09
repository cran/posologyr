% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residual_error.R
\name{error_model_comb1}
\alias{error_model_comb1}
\title{Residual error model combined 1}
\usage{
error_model_comb1(f, sigma)
}
\arguments{
\item{f}{Numeric vector, output of a pharmacokinetic model}

\item{sigma}{Numeric vector of the coefficients for the
residual error model}
}
\value{
Numeric vector, residual error
}
\description{
Residual error model combined 1. Constant error model
if no proportional coefficient is provided. Proportional
error model if no constant (or additive) error
coefficient is provided.
}
\details{
Implements the following function:
\code{g <- sigma[1] + sigma[2]*f}
}
