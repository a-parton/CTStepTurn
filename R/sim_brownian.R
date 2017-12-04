#' @title Simulate Brownian bridge.
#'
#' @description \code{sim_brownian_bridge} returns a simulated realisation of Brownian motion between two known values and times.
#'
#' @family Brownian motion simulations
#'
#' @param num_states Numeric, number of (behavioural) states that the Brownian bridge will simulate over.
#' @param fixed_values List with components:
#'   \describe{
#'   \item{\code{times}}{Vector of length two, start and end times of the simulation.}
#'   \item{\code{behav}}{Vector of length one, behaviour/state at the start time of the simulation.}
#'   \item{\code{bearings}}{Vector of length two, start and end known values of the simulation.}
#'   }
#' @param times Vector (length variable) of the times at which the brownian bridge is to be simulated over.
#' @param behavs Vector (length variable but equal to \code{times}) of the behaviour/state at each simulation time.
#' @param variance Vector of length \code{num_states}, giving the variance of the Brownian motion for each behaviour/state.
#'
#' @return Vector (length variable but equal to \code{times}) of the simulated Brownian bridge (not including known, fixed end points).
#'
#' @example examples/sim_brownian_bridge.R
#'
#' @export sim_brownian_bridge
sim_brownian_bridge <- function(num_states, fixed_values, times, behavs, variance) {

  fixed_times <- fixed_values$times; fixed_behav <- fixed_values$behav; fixed_bearings <- fixed_values$bearings

  if (num_states == 1) {

    variance_vec <- variance
  } else {

    # checks on the lengths of the input variables
    if (length(variance) != num_states) stop('Should have a turn variance for each behavioural state')
    if (length(behavs) != length(times)) stop('Behav and time length mis-match')
    # create a vector of the variance at each time the simulation is drawn at, based on the behaviours/states
    variance_vec <- variance[c(fixed_behav, behavs)]
  }

  # transform the times, so that the starting time is 0 and the process has a constant variance of 1 (includes the end time in this transformation)
  transformed_times <- cumsum(c(0, variance_vec * diff(c(fixed_times[1], times, fixed_times[2]))))

  # draw brownian motion simulation (including the final point)
  W <- cumsum(sqrt(diff(transformed_times)) * rnorm(length(times) + 1))

  # transform the brownian motion simulation to get a bridge, using a standard formula
  return(fixed_bearings[1] + W[-length(W)] - (transformed_times[-c(1,length(transformed_times))] / transformed_times[length(transformed_times)]) * (W[length(W)] - fixed_bearings[2] + fixed_bearings[1]))
}


#' @title Simulate reverse-time Brownian motion.
#'
#' @description \code{sim_brownian_backward} returns a simulated realisation of Brownian motion with known end value and start/end times.
#'
#' @family Brownian motion simulations
#'
#' @param num_states Numeric, number of (behavioural) states that the Brownian motion will simulate over.
#' @param fixed_values List with components:
#'   \describe{
#'   \item{\code{time}}{Vector of length one, end time of the simulation.}
#'   \item{\code{bearing}}{Vector of length one, end known value of the simulation.}
#'   }
#' @param times Vector (length variable) of the times at which the brownian motion is to be simulated over (including the start time).
#' @param behavs Vector (length variable but equal to \code{times}) of the behaviour/state at each simulation time.
#' @param variance Vector of length \code{num_states}, giving the variance of the Brownian motion for each behaviour/state.
#'
#' @return Vector (length variable but equal to \code{times}) of the simulated Brownian motion (not including known end point).
#'
#' @example examples/sim_brownian_backward.R
#'
#' @export sim_brownian_backward
sim_brownian_backward <- function(num_states, fixed_values, times, behavs, variance) {

  fixed_time <- fixed_values$time; fixed_bearing <- fixed_values$bearing

  if (num_states == 1) {

    variance_vec <- variance
  } else {

    # checks on the lengths of the input variables
    if (length(variance) != num_states) stop('Should have a turn variance for each behavioural state')
    if (length(behavs) != length(times)) stop('Behav and time length mis-match')
    # create a vector of the variance at each time the simulation is drawn at, based on the behaviours/states
    variance_vec <- variance[behavs]
  }

  # transform the times, so that the starting time is 0 and the process has a constant variance of 1
  # (includes the end time in this transformation)
  transformed_times <- cumsum(c(0, variance_vec * diff(c(times, fixed_time))))

  # draw reverse-time brownian motion simulation, reverse to be on the same ordering as `times'
  return(rev(fixed_bearing - cumsum(rev(sqrt(diff(transformed_times)) * rnorm(length(times))))))
}


#' @title Simulate Brownian motion.
#'
#' @description \code{sim_brownian_forward} returns a simulated realisation of Brownian motion with known start value and start/end times.
#'
#' @family Brownian motion simulations
#'
#' @param num_states Numeric, number of (behavioural) states that the Brownian motion will simulate over.
#' @param fixed_values List with components:
#'   \describe{
#'   \item{\code{time}}{Vector of length one, start time of the simulation.}
#'   \item{\code{behav}}{Vector of length one, start state of the simulation.}
#'   \item{\code{bearing}}{Vector of length one, start known value of the simulation.}
#'   }
#' @param times Vector (length variable) of the times at which the brownian motion is to be simulated over (including the end time).
#' @param behavs Vector (length variable but equal to \code{times}) of the behaviour/state at each simulation time.
#' @param variance Vector of length \code{num_states}, giving the variance of the Brownian motion for each behaviour/state.
#'
#' @return Vector (length variable but equal to \code{times}) of the simulated Brownian motion (not including known start point).
#'
#' @example examples/sim_brownian_forward.R
#'
#' @export sim_brownian_forward
sim_brownian_forward <- function(num_states, fixed_values, times, behavs, variance) {

  fixed_time <- fixed_values$time; fixed_behav <- fixed_values$behav; fixed_bearing <- fixed_values$bearing

  if (num_states == 1) {

    variance_vec <- variance
  } else {

    # checks on the lengths of the input variables
    if (length(variance) != num_states) stop('Should have a turn variance for each behavioural state')
    if (length(behavs) != length(times)) stop('Behav and time length mis-match')
    # create a vector of the variance at each time the simulation is drawn at, based on the behaviours/states
    variance_vec <- variance[c(fixed_behav, behavs[-length(behavs)])]
  }

  # transform the times, so that the starting time is 0 and the process has a constant variance of 1 (includes the start time in this transformation)
  transformed_times <- cumsum(c(0, variance_vec * diff(c(fixed_time, times))))

  # draw Brownian motion simulation
  return(fixed_bearing + cumsum(sqrt(diff(transformed_times)) * rnorm(length(times))))
}
