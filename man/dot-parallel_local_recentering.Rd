% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parallel-processing.R
\name{.parallel_local_recentering}
\alias{.parallel_local_recentering}
\title{Parallel local re-centering}
\usage{
.parallel_local_recentering(
  moderate_indices,
  fit_data,
  prepared_data,
  hrf_interface,
  parallel_config,
  theta_bounds,
  lambda_ridge = 0.01,
  verbose = FALSE
)
}
\arguments{
\item{moderate_indices}{Indices of moderate voxels}

\item{fit_data}{Current fit data}

\item{prepared_data}{Prepared input data}

\item{hrf_interface}{HRF interface functions}

\item{parallel_config}{Parallel configuration}

\item{...}{Additional parameters}
}
\value{
Updated fit results for moderate voxels
}
\description{
Perform local re-centering for moderate voxels in parallel.
}
\keyword{internal}
