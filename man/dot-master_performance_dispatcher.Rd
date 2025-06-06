% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smart_performance_dispatcher.R
\name{.master_performance_dispatcher}
\alias{.master_performance_dispatcher}
\title{MASTER PERFORMANCE DISPATCHER}
\usage{
.master_performance_dispatcher(
  n_voxels,
  n_timepoints,
  n_kernels = 4,
  kernel_length = 61,
  is_iterative = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{n_voxels}{Number of voxels to process}

\item{n_timepoints}{Number of timepoints}

\item{n_kernels}{Number of basis functions}

\item{kernel_length}{Length of each kernel}

\item{is_iterative}{Is this part of an iterative algorithm?}

\item{verbose}{Print performance decisions?}
}
\value{
Comprehensive performance configuration. The returned list
includes the smart decision objects, a human-readable strategy
summary, and a \code{functions} element referencing the smart
implementations.
}
\description{
This is the central intelligence that makes all performance decisions.
It analyzes the problem and automatically configures optimal algorithms.
}
