% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parametric-hrf-fit-class.R
\name{new_parametric_hrf_fit}
\alias{new_parametric_hrf_fit}
\title{Construct a parametric_hrf_fit object (Sprint 2 enhanced version)}
\usage{
new_parametric_hrf_fit(
  estimated_parameters,
  amplitudes,
  parameter_names,
  hrf_model = "lwu",
  r_squared = NULL,
  residuals = NULL,
  parameter_ses = NULL,
  convergence_info = list(),
  metadata = list()
)
}
\arguments{
\item{estimated_parameters}{numeric matrix of parameter estimates (voxels x parameters)}

\item{amplitudes}{numeric vector of fitted amplitudes}

\item{parameter_names}{character vector naming the parameters}

\item{hrf_model}{character string identifying the HRF model}

\item{r_squared}{numeric vector of R-squared values for each voxel}

\item{residuals}{numeric matrix of residuals (timepoints x voxels) or NULL}

\item{parameter_ses}{numeric matrix of standard errors (voxels x parameters) or NULL}

\item{convergence_info}{list of convergence diagnostics including trajectory}

\item{metadata}{list containing additional metadata such as the call,
number of voxels and time points, the parameter seed and bounds}
}
\value{
An object of class \code{parametric_hrf_fit}
}
\description{
Creates a new S3 object storing results from parametric HRF estimation.
This enhanced constructor includes additional fields for diagnostics and
standard errors introduced in Sprint 2.
}
\keyword{internal}
