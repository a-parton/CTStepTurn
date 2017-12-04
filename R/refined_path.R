#' @title Update a given section of the refined path.
#'
#' @description \code{update_refined_path} returns a MH update of a given section of the refined path.
#'
#' @family Refined path
#'
#' @param fixed_constant List with components:
#'   \describe{
#'   \item{\code{num_states}}{Numeric, the number of behavioural states.}
#'   \item{\code{indep_step}}{Boolean, True if independent step model, False if correlated.}
#'   \item{\code{obs_error}}{Boolean, True if observation error is to be included (independent or correlated).}
#'   \item{\code{corr_obs_error}}{Boolean, True if the observation error is correlated rather than independent.}
#'   \item{\code{q.bearings}}{Vector, index positions in \code{move_params} corresponding to bearing parameters.}
#'   \item{\code{q.speed_mean}}{Vector, index positions in \code{move_params} corresponding to mean speed parameters.}
#'   \item{\code{q.speed_corr}}{Vector, index positions in \code{move_params} corresponding to speed correlation parameters (only if correlated step model).}
#'   \item{\code{q.speed_var}}{Vector, index positions in \code{move_params} corresponding to speed variance parameters.}
#'   \item{\code{q.obs_error}}{Vector, index positions in \code{move_params} corresponding to error parameters.}
#'   \item{\code{ideal_refined_time_diff}}{Numeric, ideal time scale for the refined path (only if multistate case).}
#'   }
#' @param type Character, the type of section being updated (start, middle, or end) from within the full path.
#' @param fixed_values List with components set by \code{set_fixed_values}.
#' @param curr_values List with componenets set by \code{set_current_values}.
#' @param move_params Vector, movement parameters.
#' @param behav_params List with components:
#'   \describe{
#'   \item{\code{lambda}}{Vector, switching rates.}
#'   \item{\code{q}}{Square matrix (diagonal \code{NA}), switching probabilities.}
#'   }
#' @param obs List with components:
#'   \describe{
#'   \item{\code{x}}{Vector, observed x locations.}
#'   \item{\code{y}}{Vector, observed y locations.}
#'   \item{\code{times}}{Vector, observation times.}
#'   \item{\code{time_diffs}}{Vector, differences in observation times.}
#'   \item{\code{index_on_refined_path}}{Vector, index positions on refined path of observation times.}
#'   \item{\code{errors}}{List with components:
#'     \describe{
#'     \item{\code{x}}{Vector, errors in x direction.}
#'     \item{\code{y}}{Vector, errors in y direction.}
#'     }}
#'   }
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{accept}}{Boolean, True if proposal is accepted.}
#'   \item{\code{behav_fail}}{Boolean, True if proposal failed because of the rejection step in the behavioural proposal.}
#'   \item{\code{bearings}}{Vector of bearings.}
#'   \item{\code{steps}}{Vector of steps.}
#'   \item{\code{times}}{Vector of times (only in multi-behaviour).}
#'   \item{\code{behavs}}{Vector of behaviours (only in multi-behaviour).}
#'   }
#'
#' @export update_refined_path
#'
update_refined_path <- function(fixed_constant, type, fixed_values, curr_values, move_params, behav_params, obs) {

  # the three types of path section are handled separately (starting with the middle as this is the most common)
  if (type == "middle") {

    # handle single and multi behaviour cases separately
    if (fixed_constant$num_states == 1) {

      # propose a new section of bearings by simulating a brownian bridge
      prop_bearings <- sim_brownian_bridge(fixed_constant$num_states, fixed_values, times = curr_values$times, behavs = NULL, variance = move_params[fixed_constant$q.bearings])

      # construct the distribution of the steps, handle two step models separately
      if (fixed_constant$indep_step) {

        curr_step_dist <- calc_step_dist(fixed_constant, times = c(curr_values$times, fixed_values$times[2]), behavs = NULL, move_params)
        curr_step_dist$covar <- diag(curr_step_dist$var)
        prop_step_dist <- curr_step_dist
      } else {

        curr_step_dist <- prop_step_dist <- calc_dist_step_bridge(fixed_constant, fixed_values, times = curr_values$times, behavs = NULL, move_params)
      }

      # construct the distribution of the fixed locations, given the bearings (for both the new and the old)
      prop_true_locs_dist <- calc_dist_true_locs_forward(start_loc = fixed_values$start_loc, fixed_index = curr_values$loc_index, bearings = prop_bearings, step_dist = prop_step_dist)
      curr_true_locs_dist <- calc_dist_true_locs_forward(start_loc = fixed_values$start_loc, fixed_index = curr_values$loc_index, bearings = curr_values$bearings, step_dist = curr_step_dist)
      if (fixed_constant$obs_error) {
        if (fixed_constant$corr_obs_error) {
          prop_noisy_locs_dist <- calc_dist_noisy_locs_bridge(fixed_values, error_params = move_params[fixed_constant$q.obs_error], loc_dist = prop_true_locs_dist)
          curr_noisy_locs_dist <- calc_dist_noisy_locs_bridge(fixed_values, error_params = move_params[fixed_constant$q.obs_error], loc_dist = curr_true_locs_dist)
        } else {

          error_covar <- diag(c(rep(move_params[fixed_constant$q.obs_error], length(prop_true_locs_dist$loc_mean) - 2), 0, 0))
          prop_noisy_locs_dist <- list(loc_mean_error = prop_true_locs_dist$loc_mean, loc_covar_error = prop_true_locs_dist$loc_covar + error_covar, error_covar = error_covar)
          curr_noisy_locs_dist <- list(loc_mean_error = curr_true_locs_dist$loc_mean, loc_covar_error = curr_true_locs_dist$loc_covar + error_covar, error_covar = error_covar)
        }
      }
    } else {

      # propose a new section of behaviours
      prop_behav <- sim_markov_bridge(fixed_constant$num_states, fixed_values, behav_params)
      # if the rejection proposal is instantly rejected, break out of the path proposal altogether
      if (prop_behav$accept == 0) return(list(accept = 0, behav_fail = 1))

      # using the new behaviours, create the refined time scale and behaviours
      prop_times <- create_refined_times(ideal_refined_time_diff = fixed_constant$ideal_refined_time_diff, behav_proc_times = prop_behav$times, contained_obs_times = fixed_values$obs_times)
      # dimensionality issues can occur if a short section then has a switch removed, making it even shorter, so we skip out these
      if (length(prop_times$refined_times) < 3) return(list(accept = 0, behav_fail = 1))
      prop_behav_refined <- create_refined_behavs(behav_proc = prop_behav, refined_times = prop_times$refined_times)

      # propose a new section of bearings by simulating a brownian bridge
      prop_bearings <- sim_brownian_bridge(fixed_constant$num_states, fixed_values, times = prop_times$refined_times, behavs = prop_behav_refined, variance = move_params[fixed_constant$q.bearings])

      # construct the distribution of the steps (for both the new and the old as they have different times and behaviours)
      if (fixed_constant$indep_step) {

        prop_step_dist <- calc_step_dist(fixed_constant, times = c(prop_times$refined_times, fixed_values$times[2]), behavs = prop_behav_refined, move_params)
        curr_step_dist <- calc_step_dist(fixed_constant, times = c(curr_values$times, fixed_values$times[2]), behavs = curr_values$behavs, move_params)
        prop_step_dist$covar <- diag(prop_step_dist$var)
        curr_step_dist$covar <- diag(curr_step_dist$var)
      } else {

        prop_step_dist <- calc_dist_step_bridge(fixed_constant, fixed_values, times = prop_times$refined_times, behavs = prop_behav_refined, move_params)
        curr_step_dist <- calc_dist_step_bridge(fixed_constant, fixed_values, times = curr_values$times, behavs = curr_values$behavs, move_params)
      }

      # construct the distribution of the fixed locations, given the bearings (for both the new and the old)
      prop_true_locs_dist <- calc_dist_true_locs_forward(start_loc = fixed_values$start_loc, fixed_index = prop_times$loc_index, bearings = prop_bearings, step_dist = prop_step_dist)
      curr_true_locs_dist <- calc_dist_true_locs_forward(start_loc = fixed_values$start_loc, fixed_index = curr_values$loc_index, bearings = curr_values$bearings, step_dist = curr_step_dist)
      if (fixed_constant$obs_error) {
        if (fixed_constant$corr_obs_error) {

          prop_noisy_locs_dist <- calc_dist_noisy_locs_bridge(fixed_values, error_params = move_params[fixed_constant$q.obs_error], loc_dist = prop_true_locs_dist)
          curr_noisy_locs_dist <- calc_dist_noisy_locs_bridge(fixed_values, error_params = move_params[fixed_constant$q.obs_error], loc_dist = curr_true_locs_dist)
        } else {

          error_covar <- diag(c(rep(move_params[fixed_constant$q.obs_error], length(prop_true_locs_dist$loc_mean) - 2), 0, 0))
          prop_noisy_locs_dist <- list(loc_mean_error = prop_true_locs_dist$loc_mean, loc_covar_error = prop_true_locs_dist$loc_covar + error_covar, error_covar = error_covar)
          curr_noisy_locs_dist <- list(loc_mean_error = curr_true_locs_dist$loc_mean, loc_covar_error = curr_true_locs_dist$loc_covar + error_covar, error_covar = error_covar)
        }
      }
    }
  } else if (type == "start") {

    # handle single and multi behaviour cases separately
    if (fixed_constant$num_states == 1) {

      # propose a new section of bearings by simulating brownian motion backwards in time
      prop_bearings <- sim_brownian_backward(fixed_constant$num_states, fixed_values, times = curr_values$times, behavs = NULL, variance = move_params[fixed_constant$q.bearings])

      # construct the distribution of the steps, handle two step models separately
      if (fixed_constant$indep_step) {

        curr_step_dist <- calc_step_dist(fixed_constant, times = c(curr_values$times, fixed_values$time), behavs = NULL, move_params)
        curr_step_dist$covar <- diag(curr_step_dist$var)
        prop_step_dist <- curr_step_dist
      } else {

        curr_step_dist <- prop_step_dist <- calc_dist_step_backward(fixed_constant, fixed_values, times = curr_values$times, behavs = NULL, move_params)
      }

      # construct the distribution of the fixed locations, given the bearings (for both the new and the old)
      prop_true_locs_dist <- calc_dist_true_locs_backward(start_loc = fixed_values$start_loc, fixed_index = curr_values$loc_index, bearings = prop_bearings, step_dist = prop_step_dist)
      curr_true_locs_dist <- calc_dist_true_locs_backward(start_loc = fixed_values$start_loc, fixed_index = curr_values$loc_index, bearings = curr_values$bearings, step_dist = curr_step_dist)
      if (fixed_constant$obs_error) {
        if (fixed_constant$corr_obs_error) {
          prop_noisy_locs_dist <- calc_dist_noisy_locs_backward(fixed_values, error_params = move_params[fixed_constant$fixed_constant$q.obs_error], loc_dist = prop_true_locs_dist, obs)
          curr_noisy_locs_dist <- calc_dist_noisy_locs_backward(fixed_values, error_params = move_params[fixed_constant$fixed_constant$q.obs_error], loc_dist = curr_true_locs_dist, obs)
        } else {

          error_covar <- diag(rep(move_params[fixed_constant$q.obs_error], length(prop_true_locs_dist$loc_mean)))
          prop_noisy_locs_dist <- list(loc_mean_error = prop_true_locs_dist$loc_mean, loc_covar_error = prop_true_locs_dist$loc_covar + error_covar, error_covar = error_covar)
          curr_noisy_locs_dist <- list(loc_mean_error = curr_true_locs_dist$loc_mean, loc_covar_error = curr_true_locs_dist$loc_covar + error_covar, error_covar = error_covar)
        }
      }
    } else {

      # propose a new section of behaviours by simulating a markov process backwards in time
      prop_behav <- sim_markov_backward(fixed_constant$num_states, fixed_values, behav_params)

      # using the new behaviours, create the refined time scale for this section and a refined set of behaviours, as this will be more useful than switching times
      prop_times <- create_refined_times(ideal_refined_time_diff = fixed_constant$ideal_refined_time_diff, behav_proc_times = prop_behav$times, contained_obs_times = fixed_values$obs_times)
      if (length(prop_times$refined_times) < 3) return(list(accept = 0, behav_fail = 1))
      prop_behav_refined <- create_refined_behavs(behav_proc = prop_behav, refined_times = prop_times$refined_times)

      # propose a new section of bearings by simulating brownian motion backwards in time
      prop_bearings <- sim_brownian_backward(fixed_constant$num_states, fixed_values, times = prop_times$refined_times, behavs = prop_behav_refined, variance = move_params[fixed_constant$q.bearings])

      # construct the distribution of the steps (for both the new and the old as they have different times and behaviours)
      if (fixed_constant$indep_step) {

    	  prop_step_dist <- calc_step_dist(fixed_constant, times = c(prop_times$refined_times, fixed_values$time), behavs = prop_behav_refined, move_params)
        curr_step_dist <- calc_step_dist(fixed_constant, times = c(curr_values$times, fixed_values$time), behavs = curr_values$behavs, move_params)
        prop_step_dist$covar <- diag(prop_step_dist$var)
        curr_step_dist$covar <- diag(curr_step_dist$var)

    	} else {

        prop_step_dist <- calc_dist_step_backward(fixed_constant, fixed_values, times = prop_times$refined_times, behavs = prop_behav_refined, move_params)
        curr_step_dist <- calc_dist_step_backward(fixed_constant, fixed_values, times = curr_values$times, behavs = curr_values$behavs, move_params)
    	}

      # construct the distribution of the fixed locations, given the bearings (for both the new and the old)
      prop_true_locs_dist <- calc_dist_true_locs_backward(start_loc = fixed_values$start_loc, fixed_index = prop_times$loc_index_start, bearings = prop_bearings, step_dist = prop_step_dist)
      curr_true_locs_dist <- calc_dist_true_locs_backward(start_loc = fixed_values$start_loc, fixed_index = curr_values$loc_index, bearings = curr_values$bearings, step_dist = curr_step_dist)
      if (fixed_constant$obs_error) {
        if (fixed_constant$corr_obs_error) {
          prop_noisy_locs_dist <- calc_dist_noisy_locs_backward(fixed_values, error_params = move_params[fixed_constant$q.obs_error], loc_dist = prop_true_locs_dist, obs)
          curr_noisy_locs_dist <- calc_dist_noisy_locs_backward(fixed_values, error_params = move_params[fixed_constant$q.obs_error], loc_dist = curr_true_locs_dist, obs)
        } else {

          error_covar <- diag(rep(move_params[fixed_constant$q.obs_error], length(prop_true_locs_dist$loc_mean)))
          prop_noisy_locs_dist <- list(loc_mean_error = prop_true_locs_dist$loc_mean, loc_covar_error = prop_true_locs_dist$loc_covar + error_covar, error_covar = error_covar)
          curr_noisy_locs_dist <- list(loc_mean_error = curr_true_locs_dist$loc_mean, loc_covar_error = curr_true_locs_dist$loc_covar + error_covar, error_covar = error_covar)
        }
      }
    }
  } else if (type == "end") {

    # handle single and multi behaviour cases separately
    if (fixed_constant$num_states == 1) {

      # propose a new section of bearings by simulating brownian motion forwards in time
      prop_bearings <- sim_brownian_forward(fixed_constant$num_states, fixed_values, times = curr_values$times, behavs = NULL, variance = move_params[fixed_constant$q.bearings])

      # construct the distribution of the steps, handle two step models separately
      if (fixed_constant$indep_step) {

        curr_step_dist <- calc_step_dist(fixed_constant, times = c(curr_values$times, fixed_values$after_time), behavs = NULL, move_params)
        curr_step_dist$covar <- diag(curr_step_dist$var)
        prop_step_dist <- curr_step_dist
      } else {

        curr_step_dist <- prop_step_dist <- calc_dist_step_forward(fixed_constant, fixed_values, times = curr_values$times, behavs = NULL, move_params)
      }

      # construct the distribution of the fixed locations, given the bearings (for both the new and the old)
      prop_true_locs_dist <- calc_dist_true_locs_forward(start_loc = fixed_values$start_loc, fixed_index = curr_values$loc_index, bearings = prop_bearings, step_dist = prop_step_dist)
      curr_true_locs_dist <- calc_dist_true_locs_forward(start_loc = fixed_values$start_loc, fixed_index = curr_values$loc_index, bearings = curr_values$bearings, step_dist = curr_step_dist)
      if (fixed_constant$obs_error) {
        if (fixed_constant$corr_obs_error) {
          prop_noisy_locs_dist <- calc_dist_noisy_locs_forward(fixed_values, error_params = move_params[fixed_constant$q.obs_error], loc_dist = prop_true_locs_dist, obs)
          curr_noisy_locs_dist <- calc_dist_noisy_locs_forward(fixed_values, error_params = move_params[fixed_constant$q.obs_error], loc_dist = curr_true_locs_dist, obs)
        } else {

          error_covar <- diag(rep(move_params[fixed_constant$q.obs_error], length(prop_true_locs_dist$loc_mean)))
          prop_noisy_locs_dist <- list(loc_mean_error = prop_true_locs_dist$loc_mean, loc_covar_error = prop_true_locs_dist$loc_covar + error_covar, error_covar = error_covar)
          curr_noisy_locs_dist <- list(loc_mean_error = curr_true_locs_dist$loc_mean, loc_covar_error = curr_true_locs_dist$loc_covar + error_covar, error_covar = error_covar)
        }
      }
    } else {

      # propose a new section of behaviours by simulating a markov process forwards in time
      prop_behav <- sim_markov_forward(fixed_constant$num_states, fixed_values, behav_params)

      # using the new behaviours, create the refined time scale for this section and a refined set of behaviours, as this will be more useful than switching times
      prop_times <- create_refined_times(ideal_refined_time_diff = fixed_constant$ideal_refined_time_diff, behav_proc_times = prop_behav$times, contained_obs_times = fixed_values$obs_times)
      if (length(prop_times$refined_times) < 3) return(list(accept = 0, behav_fail = 1))
      prop_behav_refined <- create_refined_behavs(behav_proc = prop_behav, refined_times = prop_times$refined_times)

      # propose a new section of bearings by simulating brownian motion forwards in time
      prop_bearings <- sim_brownian_forward(fixed_constant$num_states, fixed_values, times = prop_times$refined_times, behavs = prop_behav_refined, variance = move_params[fixed_constant$q.bearings])

      # construct the distribution of the steps (for both the new and the old as they have different times and behaviours)
      if (fixed_constant$indep_step) {

    	  prop_step_dist <- calc_step_dist(fixed_constant, times = c(prop_times$refined_times, fixed_values$after_time), behavs = prop_behav_refined, move_params)
        curr_step_dist <- calc_step_dist(fixed_constant, times = c(curr_values$times, fixed_values$after_time), behavs = curr_values$behavs, move_params)
        prop_step_dist$covar <- diag(prop_step_dist$var)
        curr_step_dist$covar <- diag(curr_step_dist$var)

    	} else {

        prop_step_dist <- calc_dist_step_forward(fixed_constant, fixed_values, times = prop_times$refined_times, behavs = prop_behav_refined, move_params)
        curr_step_dist <- calc_dist_step_forward(fixed_constant, fixed_values, times = curr_values$times, behavs = curr_values$behavs, move_params)
    	}

      # construct the distribution of the fixed locations, given the bearings (for both the new and the old)
      prop_true_locs_dist <- calc_dist_true_locs_forward(start_loc = fixed_values$start_loc, fixed_index = prop_times$loc_index, bearings = prop_bearings, step_dist = prop_step_dist)
      curr_true_locs_dist <- calc_dist_true_locs_forward(start_loc = fixed_values$start_loc, fixed_index = curr_values$loc_index, bearings = curr_values$bearings, step_dist = curr_step_dist)
      if (fixed_constant$obs_error) {
        if (fixed_constant$corr_obs_error) {
          prop_noisy_locs_dist <- calc_dist_noisy_locs_forward(fixed_values, error_params = move_params[fixed_constant$q.obs_error], loc_dist = prop_true_locs_dist, obs)
          curr_noisy_locs_dist <- calc_dist_noisy_locs_forward(fixed_values, error_params = move_params[fixed_constant$q.obs_error], loc_dist = curr_true_locs_dist, obs)
        } else {

          error_covar <- diag(rep(move_params[fixed_constant$q.obs_error], length(prop_true_locs_dist$loc_mean)))
          prop_noisy_locs_dist <- list(loc_mean_error = prop_true_locs_dist$loc_mean, loc_covar_error = prop_true_locs_dist$loc_covar + error_covar, error_covar = error_covar)
          curr_noisy_locs_dist <- list(loc_mean_error = curr_true_locs_dist$loc_mean, loc_covar_error = curr_true_locs_dist$loc_covar + error_covar, error_covar = error_covar)
        }
      }
    }
  }

  # acceptance probability for the proposed path over the current
  if (fixed_constant$obs_error) {
    accept_prob <- exp(dmvn(fixed_values$locs, mu = prop_noisy_locs_dist$loc_mean_error, sigma = prop_noisy_locs_dist$loc_covar_error, log = TRUE) -
                       dmvn(fixed_values$locs, mu = curr_noisy_locs_dist$loc_mean_error, sigma = curr_noisy_locs_dist$loc_covar_error, log = TRUE))
  } else {
    accept_prob <- exp(dmvn(fixed_values$locs, mu = prop_true_locs_dist$loc_mean, sigma = prop_true_locs_dist$loc_covar, log = TRUE) -
                       dmvn(fixed_values$locs, mu = curr_true_locs_dist$loc_mean, sigma = curr_true_locs_dist$loc_covar, log = TRUE))
  }


  if (runif(1) < accept_prob) {

    if (fixed_constant$corr_obs_error) {
      # propose a new section of steps by drawing from the conditioned distribution
      prop_steps <- sim_steps(step_mean = prop_step_dist$mean, step_covar = prop_step_dist$covar,
                              A = prop_true_locs_dist$A,
                              loc_covar_error = prop_noisy_locs_dist$loc_covar_error, error_mean = prop_noisy_locs_dist$error_mean, error_covar = prop_noisy_locs_dist$error_covar,
                              loc_step_covar = prop_true_locs_dist$loc_step_covar,
                              start_loc = fixed_values$start_loc, fixed_locs = fixed_values$locs, type)
    } else if (fixed_constant$obs_error){
      prop_steps <- sim_steps(step_mean = prop_step_dist$mean, step_covar = prop_step_dist$covar,
                              A = prop_true_locs_dist$A,
                              loc_covar_error = prop_noisy_locs_dist$loc_covar_error, error_mean = rep(0, nrow(prop_true_locs_dist$A)), error_covar = prop_noisy_locs_dist$error_covar,
                              loc_step_covar = prop_true_locs_dist$loc_step_covar,
                              start_loc = fixed_values$start_loc, fixed_locs = fixed_values$locs, type)
    } else {
      prop_steps <- sim_steps(step_mean = prop_step_dist$mean, step_covar = prop_step_dist$covar,
                              A = prop_true_locs_dist$A,
                              loc_covar_error = prop_true_locs_dist$loc_covar, error_mean = NA, error_covar = NA,
                              loc_step_covar = prop_true_locs_dist$loc_step_covar,
                              start_loc = fixed_values$start_loc, fixed_locs = fixed_values$locs, type)
    }


    if(fixed_constant$num_states == 1) {

      return(list(accept = 1, bearings = prop_bearings, steps = prop_steps))
    } else {

      return(list(accept = 1, bearings = prop_bearings, steps = prop_steps, times = prop_times$refined_times, behavs = prop_behav_refined))
    }
  } else {

    return(list(accept = 0, behav_fail = 0))
  }
}



#' @title Simulate conditional steps.
#'
#' @description \code{sim_steps} simulates a set of conditional steps.
#'
#' @family Refined path
#'
#' @param step_mean Vector, mean of unconditional steps.
#' @param step_covar Square matrix, covariance of unconditional steps.
#' @param A linear constraint matrix.
#' @param loc_covar_error Square matrix, covariance of observed locations (including error if applicable).
#' @param error_mean Vector, mean of the errors.
#' @param error_covar Square matrix, covariance of the errors.
#' @param loc_step_covar Matrix, covariance between unconditional steps and the observed locations.
#' @param start_loc Vector, length 2, giving the location at the start of the section.
#' @param fixed_locs Vector, giving the observed locations along the section.
#' @param type Character ("start", "middle", "end") describing where on the full path the section is from.
#'
#' @return Vector of simulated steps.
#'
#' @export sim_steps
#'
sim_steps <- function(step_mean, step_covar, A, loc_covar_error, error_mean, error_covar, loc_step_covar, start_loc, fixed_locs, type) {

  inv_loc_covar <- solve(loc_covar_error)
  uncond_steps <- as.vector(rmvn(1, mu = step_mean, sigma = step_covar))

  if (length(error_covar) > 1) {

    if (type == "middle") {

      if (length(fixed_locs) == 2) {

        true_loc <- fixed_locs
      } else {

        true_loc <- fixed_locs - error_mean + c(rmvn(1, rep(0, length(fixed_locs)-2), error_covar[1:(length(fixed_locs)-2),1:(length(fixed_locs)-2)]), 0, 0)
      }
    } else {

      true_loc <- c(rmvn(1, fixed_locs - error_mean, error_covar))
    }
    uncond_steps - as.vector(loc_step_covar %*% inv_loc_covar %*% ((start_loc + A %*% uncond_steps) - true_loc))

  } else {
    uncond_steps - as.vector(loc_step_covar %*% inv_loc_covar %*% ((start_loc + A %*% uncond_steps) - fixed_locs))
  }
}
