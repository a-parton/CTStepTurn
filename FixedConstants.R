#' @title Set fixed inference values.
#'
#' @description \code{FixedConstants} returns a list of values controlling the MCMC sampler.
#'
#' @param pert_sd Vector of perturbation variances for the speed parameters (one for each state, laid out as mean[1],mean[2],var[1],var[2]).
#' @param error_pert_sd Vector of perturbation variances for the error parameters (only applicable if correlated errors are involved).
#' @param num_states Numeric value, number of behavioural states.
#' @param indep_step Boolean, True if independent step model is to be used, False for correlated steps.
#' @param obs_error Boolean, True if any observation error is to be assumed.
#' @param corr_obs_error Boolean, True if the error is correlated, False if no error or it is independent.
#' @param turn_prec_prior_shape Vector of length \code{num_states} with shape prior for turn variance.
#' @param turn_prec_prior_rate Vector of length \code{num_states} with rate prior for turn variance.
#' @param error_prec_prior_shape Only applicable for independent errors, numeric shape prior for error variance.
#' @param error_prec_prior_rate Only applicable for correlated errors, numeric rate prior for error variance.
#'
#' @return List of fixed constants, to be used in the MCMC sampler.

FixedConstants <- function(pert_sd, error_pert_sd = NULL,
                           num_states = 1,
                           indep_step = FALSE, obs_error = FALSE, corr_obs_error = FALSE,
                           turn_prec_prior_shape = 0.001, turn_prec_prior_rate = 0.001,
                           error_prec_prior_shape = NULL, error_prec_prior_rate = NULL)
{

  if (indep_step) {
    num_move_params <- 3 * num_states + obs_error + corr_obs_error
    num_speed_param <- 2 * num_states
    q.bearings <- 1:num_states
    q.speed_mean <- (num_states+1):(2*num_states)
    q.speed_var <- (2*num_states+1):(2*num_states)
    if(corr_obs_error) q.obs_error <- (3*num_states+1):(3*num_states+2)
    else q.obs_error <- 3*num_states+1

    fc <- list(
      pert_sd = pert_sd,
      error_pert_sd = error_pert_sd,
      num_states = 1,
      indep_step = indep_step, obs_error = obs_error, corr_obs_error = corr_obs_error,
      num_move_params = num_move_params, num_speed_param = num_speed_param,
      q.bearings = q.bearings, q.speed_mean = q.speed_mean, q.speed_var = q.speed_var,
      q.obs_error = q.obs_error,
      turn_prec_prior = list(rate = turn_prec_prior_rate,
                             shape = turn_prec_prior_shape),
      error_prec_prior = list(rate = error_prec_prior$rate,
                              shape = error_prec_prior$shape)
    )
  } else {
    num_move_params <- 4 * num_states + obs_error + corr_obs_error
    num_speed_param <- 3 * num_states
    q.bearings <- 1:num_states
    q.speed_mean <- (num_states+1):(2*num_states)
    q.speed_corr <- (2*num_states+1):(3*num_states)
    q.speed_var <- (3*num_states+1):(4*num_states)
    if(corr_obs_error) q.obs_error <- (4*num_states+1):(4*num_states+2)
    else q.obs_error <- 4*num_states+1

    fc <- list(
      pert_sd = pert_sd,
      error_pert_sd = error_pert_sd,
      num_states = 1,
      indep_step = indep_step, obs_error = obs_error, corr_obs_error = corr_obs_error,
      num_move_params = num_move_params, num_speed_param = num_speed_param,
      q.bearings = q.bearings, q.speed_mean = q.speed_mean, q.speed_corr = q.speed_corr, q.speed_var = q.speed_var,
      q.obs_error = q.obs_error,
      turn_prec_prior = list(rate = turn_prec_prior_rate,
                             shape = turn_prec_prior_shape),
      error_prec_prior = list(rate = error_prec_prior_rate,
                              shape = error_prec_prior_shape)
    )
  }

  ## Set the name for the class
  class(fc) <- append(class(fc),"FixedConstants")
  return(fc)
}
