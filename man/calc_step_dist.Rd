% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_step_dist.R
\name{calc_step_dist}
\alias{calc_step_dist}
\title{Step distribution of a path section (independent steps).}
\usage{
calc_step_dist(fixed_constant, times, behavs, move_params)
}
\arguments{
\item{fixed_constant}{List with components:
\describe{
\item{\code{num_states}}{Numeric, number of behavioural states.}
\item{\code{q.speed_mean}}{Vector, index values in \code{move_params} which correspond to mean speed parameters.}
\item{\code{q.speed_var}}{Vector, index values in \code{move_params} which correspond to speed variance parameters.}
}}

\item{times}{Vector of refined times.}

\item{behavs}{Vector of refined behaviours.}

\item{move_params}{Vector of movement parameters.}
}
\value{
List with the components:
  \describe{
  \item{\code{mean}}{Vector of mean step lengths.}
  \item{\code{sd}}{Vector of step length standard deviations.}
  \item{\code{var}}{Vector of step length variances.}
  }
}
\description{
\code{calc_step_dist} returns the distribution of a set of steps when the independent step model is assumed.
}
\seealso{
Other Step distribution: \code{\link{calc_dist_OU_backward}},
  \code{\link{calc_dist_OU_bridge}},
  \code{\link{calc_dist_OU_forward}},
  \code{\link{calc_dist_step_backward}},
  \code{\link{calc_dist_step_bridge}},
  \code{\link{calc_dist_step_forward}}
}
\concept{Step distribution}
