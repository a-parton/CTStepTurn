#' @title Simulate Markov chain bridge.
#'
#' @description \code{sim_markov_bridge} returns a simulated realisation of a Markov chain between two known values and times.
#'
#' @details Simulation is by rejection. Markov chain simulated forwards in time from the starting point, and then compared with the
#'   known end point. If the end point is correct, the simulation is accepted, otherwise it is rejected.
#'
#' @family Markov chain simulations
#'
#' @param num_states Numeric, number of (behavioural) states in the Markov chain.
#' @param fixed_values List with components:
#'   \describe{
#'   \item{\code{inc_times}}{Vector of length two, start and end times of the simulation.}
#'   \item{\code{inc_behavs}}{Vector of length two, start and end known values (behaviours) of the simulation.}
#'   }
#' @param behav_params List with components:
#'   \describe{
#'   \item{\code{lambda}}{Vector of length \code{num_states}, switching rates out of each state.}
#'   \item{\code{q}}{Sqaure matrix size \code{num_states}, probability of switching from (row) each state to (column) another. Diagonal elements are NA.}
#'   }
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{accept}}{Binary element that the rejection simulation has been accepted/rejected. If rejected, this is the only return value.}
#'   \item{\code{times}}{Vector, switching times of the simulation. Includes the fixed end points.}
#'   \item{\code{states}}{Vector, states at switching times of the simulation. Includes the fixed end points.}
#'   }
#'
#' @example examples/sim_markov_bridge.R
#'
#' @export sim_markov_bridge
#'
sim_markov_bridge <- function(num_states, fixed_values, behav_params) {

  fixed_times <- fixed_values$inc_times; fixed_states <- fixed_values$inc_behavs
  # intialise switching times and states
  switch_times <- curr_time <- fixed_times[1]
  switch_states <- curr_state <- fixed_states[1]

  repeat {

    # find out when the time of the next switch will occur
    curr_time <- curr_time + rexp(1, behav_params$lambda[curr_state])
    # if this is past the end time, then we stop and do not consider it
    if (curr_time > fixed_times[2]) break

    # otherwise we choose what state we are switching into
    # handle 2 state case separately otherwise sample() will behave unexpectedly
    curr_state <- if (num_states == 2) (1:2)[-curr_state] else sample((1:num_states)[-curr_state], 1, prob = behav_params$q[curr_state, -curr_state])

    # then add this information to a running vector of switches
    switch_times <- c(switch_times, curr_time)
    switch_states <- c(switch_states, curr_state)
  }

    # if the final state matches the fixed end state, then the simulation is accepted, otherwise it is rejected
    if (curr_state == fixed_states[2]) return(list(accept = TRUE, times = c(switch_times, fixed_times[2]), states = c(switch_states, curr_state)))
    return(list(accept = FALSE))
}



#' @title Simulate Markov chain backwards in time.
#'
#' @description \code{sim_markov_backward} returns a simulated realisation of a Markov chain backwards in time from a known value and time.
#'
#' @family Markov chain simulations
#' @seealso \code{\link{calc_markov_steady_state}} and \code{\link{calc_reverse_markov_gen}} for calculating
#'   Markov steady state/equilibrium distribution and reverse-time generator matrices, respectively.
#'
#' @param num_states Numeric, number of (behavioural) states in the Markov chain.
#' @param fixed_values List with components:
#'   \describe{
#'   \item{\code{inc_times}}{Vector of length two, start and end times of the simulation.}
#'   \item{\code{inc_behav}}{Vector of length one, end known value (behaviour) of the simulation.}
#'   }
#' @param behav_params List with components:
#'   \describe{
#'   \item{\code{lambda}}{Vector of length \code{num_states}, switching rates out of each state.}
#'   \item{\code{q}}{Sqaure matrix size \code{num_states}, probability of switching from (row) each state to (column) another. Diagonal elements are NA.}
#'   }
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{times}}{Vector, switching times of the simulation. Includes the fixed end points.}
#'   \item{\code{states}}{Vector, states at switching times of the simulation. Includes the fixed end points.}
#'   }
#'
#' @example examples/sim_markov_backward.R
#'
#' @export sim_markov_backward
#'
sim_markov_backward <- function(num_states, fixed_values, behav_params) {

  fixed_times <- fixed_values$inc_times
  # compute the time-reversal generator matrix
  if (num_states == 2) { # reverse process is the same when there are only 2 states
    reverse_behav_param <- behav_params
  } else {
    reverse_behav_param <- calc_reverse_markov_gen(num_states, behav_params)
  }

  # set up values (these are the fixed end values)
  switch_times <- curr_time <- fixed_times[2]
  switch_states <- curr_state <- fixed_values$inc_behav

  repeat {

    # find out when the time of the next switch will occur, we are taking the time away because we're going backwards in time
    curr_time <- curr_time - rexp(1, reverse_behav_param$lambda[curr_state])

    # if this is past the end time (which is actually the start time because we are going backwards), then we stop and do not consider it
    if (curr_time < fixed_times[1]) break

    # then add this information to a running vector of switches.
    # do this now because of the backwards simulation, but want vectors to appear like the forwards time case
    switch_times <- c(curr_time, switch_times)
    switch_states <- c(curr_state, switch_states)

    # choose what state we are switching into
    # handle 2 state case separately otherwise sample() will behave unexpectedly
    curr_state <- if (num_states == 2) (1:2)[-curr_state] else sample((1:num_states)[-curr_state], 1, prob = reverse_behav_param$q[curr_state, -curr_state])
  }

  return(list(times = c(fixed_times[1], switch_times), states = c(curr_state, switch_states)))
}


#' @title Calculate the generator matrix of a Markov chain in reverse time.
#'
#' @description \code{calc_reverse_markov_gen} returns the generator matrix of the Markov chain which is the time-reversed process of a given Markov chain.
#'
#' @seealso \code{\link{calc_markov_steady_state}} for calculating a Markov steady state/equilibrium distribution and
#'  \code{\link{sim_markov_backward}} for simulating from this reverse-time process.
#'
#' @param num_states Numeric, number of (behavioural) states in the Markov chain.
#' @param behav_params List with components:
#'   \describe{
#'   \item{\code{lambda}}{Vector of length \code{num_states}, switching rates out of each state.}
#'   \item{\code{q}}{Sqaure matrix size \code{num_states}, probability of switching from (row) each state to (column) another. Diagonal elements are NA.}
#'   }
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{times}}{Vector, switching times of the simulation. Includes the fixed end points.}
#'   \item{\code{states}}{Vector, states at switching times of the simulation. Includes the fixed end points.}
#'   }
#'
#' @example examples/calc_reverse_markov_gen.R
#'
#' @export calc_reverse_markov_gen
#'
calc_reverse_markov_gen <- function(num_states, behav_params) {

  # calculate the steady state solution of the markov chain
  ss <- calc_markov_steady_state(num_states, behav_params)
  steady_state <- ss$steady_state

  # create the reverse generator matrix as the transpose of the original multiplied by steady state ratio
  generator <- ss$generator
  rev_generator <- t(generator * outer(X = steady_state, Y = steady_state, FUN = function(X, Y) {X / Y}))

  # calculate the reverse switching probabilities based on the rates
  rev_q <- matrix(NA, nrow = num_states, ncol = num_states)
  diag(rev_q) <- NA
  for (i in 1:num_states) rev_q[i, -i] <- rev_generator[i, -i] / rep(behav_params$lambda[i], num_states - 1)

  # switching rates remain the same as in the forward case
  return(list(lambda = behav_params$lambda, q = rev_q))
}


#' @title Simulate Markov chain.
#'
#' @description \code{sim_markov_forward} returns a simulated realisation of a Markov chain given known start point.
#'
#' @family Markov chain simulations
#'
#' @param num_states Numeric, number of (behavioural) states in the Markov chain.
#' @param fixed_values List with components:
#'   \describe{
#'   \item{\code{inc_times}}{Vector of length two, start and end times of the simulation.}
#'   \item{\code{inc_behav}}{Vector of length one, start known value (behaviour) of the simulation.}
#'   }
#' @param behav_params List with components:
#'   \describe{
#'   \item{\code{lambda}}{Vector of length \code{num_states}, switching rates out of each state.}
#'   \item{\code{q}}{Sqaure matrix size \code{num_states}, probability of switching from (row) each state to (column) another. Diagonal elements are NA.}
#'   }
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{times}}{Vector, switching times of the simulation. Includes the fixed end points.}
#'   \item{\code{states}}{Vector, states at switching times of the simulation. Includes the fixed end points.}
#'   }
#'
#' @example examples/sim_markov_forward.R
#'
#' @export sim_markov_forward
#'
sim_markov_forward <- function(num_states, fixed_values, behav_params) {

  fixed_times <- fixed_values$inc_times
  # set up values
  switch_times <- curr_time <- fixed_times[1]
  switch_states <- curr_state <- fixed_values$inc_behav

  repeat {

    # find out when the time of the next switch will occur
    curr_time <- curr_time + rexp(1, behav_params$lambda[curr_state])

    # if this is past the end time, then we stop and do not consider it
    if (curr_time > fixed_times[2]) break

    # otherwise we choose what state we are switching into
    # handle 2 state case separately otherwise sample() will behave unexpectedly
    curr_state <- if (num_states == 2) (1:2)[-curr_state] else sample((1:num_states)[-curr_state], 1, prob = behav_params$q[curr_state, -curr_state])

    # then add this information to a running vector of switches
    switch_times <- c(switch_times, curr_time)
    switch_states <- c(switch_states, curr_state)
  }

  return(list(times = c(switch_times, fixed_times[2]), states = c(switch_states, curr_state)))
}
