#' @title Simulate behaviour parameters.
#'
#' @description \code{update_behav_param} returns a conditional Gibbs sample of the behaviour parameters.
#'
#' @family Behaviour parameters
#'
#' @param fixed_constant List with components:
#'   \describe{
#'   \item{\code{num_states}}{Numeric, number of behavioural states.}
#'   \item{\code{prior_param}}{List with components:
#'     \describe{
#'     \item{\code{lambda_rate}}{Vector (length \code{num_states}) with prior rate of the switching rates.}
#'     \item{\code{lambda_shape}}{Vector (length \code{num_states}) with prior shape of the switching rates.}
#'     \item{\code{q_conc}}{Matrix (square, size \code{num_states} with \code{NA} diagonal elements) with prior concentration of the switching probabilities.}
#'     }}
#'   }
#' @param behav_proc List with components:
#'   \describe{
#'   \item{\code{times}}{Vector, switch times of the behaviour process.}
#'   \item{\code{states}}{Vector, behavioural states at switching times.}
#'   }
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{lambda}}{Vector (length \code{num_states}) of sampled switching rates.}
#'   \item{\code{q}}{Matrix (square size \code{num_states} with \code{NA} diagonal elements) of sampled switching probabilities. Note: row = from, column = to.}
#'   }
#'
#' @example examples/update_behav_param.R
#'
#' @export update_behav_param
#'
update_behav_param <- function(fixed_constant, behav_proc) {

  switch_times <- behav_proc$times; switch_states <- behav_proc$states
  num_states <- fixed_constant$num_states; prior_param <- fixed_constant$prior_param

  # check that the input is a valid Markov chain
  check_valid_markov(num_states, switch_times, switch_states)

  # calculate the sufficient statistics for this Markov chain
  suff_stat <- calc_markov_suff_stats(num_states, switch_times, switch_states)

  # calculate the posterior parameters given sufficient statistics
  post_param <- calc_markov_post_param(prior_param, num_states, suff_stat)

  # draw sample from posterior distribution
  lambda <- rgamma(num_states, shape = post_param$lambda_shape, rate = post_param$lambda_rate)
  if (num_states == 2) q <- matrix(c(NA,1,1,NA), nrow = 2, ncol = 2) else {
    q <- matrix(NA, nrow = num_states, ncol = num_states)
    for (i in 1:num_states) q[i, -i] <- rdirichlet(1, alpha = post_param$q_conc[i, -i]) # sample q in rows, then they sum to one as required
  }

  return(list(lambda = lambda, q = q))
}


#' @title Check Markov chain is valid.
#'
#' @description \code{check_valid_markov} is a testing function to check that a Markov chain is valid.
#'
#' @family Behaviour parameters
#'
#' @param num_states Numeric, number of behavioural states.
#' @param switch_times Vector, times of behavioural switches.
#' @param switch_states Vector, states at switching times.
#'
#' @return \code{NULL}
#'
#' @export check_valid_markov
#'
check_valid_markov <- function(num_states, switch_times, switch_states) {

  # check if times same length as states, throw exception if not
  if (length(switch_times) != length(switch_states)) stop('Time and state vectors not of equal length')
  # check that times is strictly increasing vector
  if (!identical(switch_times, sort(switch_times))) stop('Time vector is not strictly increasing')
  # check there is only valid switch_states in the state vector
  if (!all(switch_states %in% 1:num_states)) stop('State vector contains invalid switch_states')
  # check that no state switches to itself, but that the last state does (to give an end point to the process)
  if (any(diff(switch_states[-length(switch_states)]) == 0) || !identical(switch_states[length(switch_states)], switch_states[length(switch_states) - 1])) stop('State vector is invalid process')
}


#' @title Calculate Markov chain sufficient statistics.
#'
#' @description \code{calc_markov_suff_stats} calculates the sufficient statistics (time in each state and the number of transitions between state pairs) for a Markov chain.
#'
#' @family Behaviour parameters
#'
#' @param num_states Numeric, number of behavioural states.
#' @param switch_times Vector, times of behavioural switches.
#' @param switch_states Vector, states at switching times.
#'
#' @return List with components:
#'   \describe{
#'   \item{\code{time_in_state}}{Vector (length \code{num_states}) with total time spent in each behavioural state.}
#'   \item{\code{num_trans}}{Matrix (square, size \code{num_states} with \code{NA} diagonal elements) with number of transitions from each state (row) to every other (column).}
#'   }
#'
#' @example examples/calc_markov_suff_stats.R
#'
#' @export calc_markov_suff_stats
#'
calc_markov_suff_stats <- function(num_states, switch_times, switch_states) {

  # initialise sufficient statistics at 0
  time_in_state <- rep(0, times = num_states)
  num_trans <- matrix(0, nrow = num_states, ncol = num_states) # row is "from", column is "to"
  diag(num_trans) <- NA # can't switch to the same state so diagonal is always NA

  # loop over each switching time and increment the sufficient statistics
  for (i in 2:length(switch_times)) {

    time_in_state[switch_states[i - 1]] <- time_in_state[switch_states[i - 1]] + switch_times[i] - switch_times[i - 1]
    num_trans[switch_states[i - 1], switch_states[i]] <- num_trans[switch_states[i - 1], switch_states[i]] + 1
  }

  return(list(time_in_state = time_in_state, num_trans = num_trans))
}


##' @title Conditional posterior distribution of Markov chain.
#'
#' @description \code{calc_markov_post_param} returns the conditional posterior distribution of the behaviour parameters (given prior and observation of the chain).
#'
#' @family Behaviour parameters
#'
#' @param prior_param List with components:
#'   \describe{
#'   \item{\code{lambda_rate}}{Vector (length \code{num_states}) with prior rate of the switching rates.}
#'   \item{\code{lambda_shape}}{Vector (length \code{num_states}) with prior shape of the switching rates.}
#'   \item{\code{q_conc}}{Matrix (square, size \code{num_states} with \code{NA} diagonal elements) with prior concentration of the switching probabilities.}
#'   }
#' @param num_states Numeric, number of behavioural states.
#' @param suff_stat List with components:
#'   \describe{
#'   \item{\code{time_in_state}}{Vector (length \code{num_states}) with total time spent in each behavioural state.}
#'   \item{\code{num_trans}}{Matrix (square, size \code{num_states} with \code{NA} diagonal elements) with number of transitions from each state (row) to every other (column).}
#'   }
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{lambda_rate}}{Vector (length \code{num_states}) with posterior rate of the switching rates.}
#'   \item{\code{lambda_shape}}{Vector (length \code{num_states}) with posterior shape of the switching rates.}
#'   \item{\code{q_conc}}{Matrix (square, size \code{num_states} with \code{NA} diagonal elements) with posterior concentration of the switching probabilities.}
#'   }
#'
#' @example examples/calc_markov_post_param.R
#'
#' @export calc_markov_post_param
#'
calc_markov_post_param <- function(prior_param, num_states, suff_stat) {

  # initialise the posterior parameters
  lambda_shape <- rep(NA, times = num_states)
  lambda_rate <- rep(NA, times = num_states)

  if (num_states == 2) { # ignore q in this case at it doesn't play a part

    for (i in 1:2) {

      # posterior lambda shape is prior plus the number of transitions out of that state (row of matrix)
      lambda_shape[i] <- prior_param$lambda_shape[i] + sum(suff_stat$num_trans[i, -i])
      # posterior lambda rate is prior plus time in that state
      lambda_rate[i] <- prior_param$lambda_rate[i] + suff_stat$time_in_state[i]
      q_conc <- matrix(c(NA, 1, 1, NA), nrow = num_states, ncol = num_states)
    }
  } else {

    q_conc <- matrix(NA, nrow = num_states, ncol = num_states)
    for (i in 1:num_states) {

      lambda_shape[i] <- prior_param$lambda_shape[i] + sum(suff_stat$num_trans[i, -i])
      lambda_rate[i] <- prior_param$lambda_rate[i] + suff_stat$time_in_state[i]
      # posterior q concentration is prior plus transitions between the pairs (complete a row at a time)
      # NA's take care of themselves and remain NA's
      q_conc[i, ] <- prior_param$q_conc[i, ] + suff_stat$num_trans[i, ]
    }
  }

  return(list(lambda_shape = lambda_shape, lambda_rate = lambda_rate, q_conc = q_conc))
}


#' @title Simulate speed parameters.
#'
#' @description \code{update_speed_param} returns a conditional MH proposal (with acceptance probability) of the speed parameters.
#'
#' @family Speed parameters
#'
#' @param fixed_constant List with components:
#'   \describe{
#'   \item{\code{pert_sd}}{Vector, perturbation standard deviations for the speed parameter proposals.}
#'   \item{\code{num_states}}{Numeric, number of behavioural states.}
#'   \item{\code{indep_step}}{Boolean, True if independent step model, False if correlated step model.}
#'   \item{\code{num_speed_param}}{Numeric, number of speed parameters.}
#'   \item{\code{q.bearings}}{Vector, index positions within \code{curr_move_params} relating to bearing parameters.}
#'   \item{\code{q.speed_mean}}{Vector, index positions within \code{move_params} relating to speed mean parameters.}
#'   \item{\code{q.speed_corr}}{Vector, index positions within \code{move_params} relating to speed correlation parameters.}
#'   \item{\code{q.speed_var}}{(Will only apply if \code{indep_step} is False) Vector, index positions within \code{move_params} relating to speed variance parameters.}
#'   \item{\code{q.obs_error}}{Vector, index positions within \code{curr_move_params} relating to observation error parameters.}
#'   }
#' @param curr_move_params Vector of movement parameters (as turn_var_1, turn_var_2, mean_speed_1, .., speed_var_2).
#' @param behav_params List with components:
#'   \describe{
#'   \item{\code{lambda}}{Vector (length \code{num_states}) of switching rates.}
#'   \item{\code{q}}{Matrix (square size \code{num_states} with \code{NA} diagonal elements) of switching probabilities. Note: row = from, column = to.}
#'   }
#' @param refined_path List with components:
#'   \describe{
#'   \item{\code{time_diffs}}{Vector, differences in the refined time scale.}
#'   \item{\code{behavs}}{Vector, refined behavioural/state sequence.}
#'   \item{\code{steps}}{Vector, refined step process.}
#'   }
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{accept_prob}}{Numeric, probability of accepting the MH proposal.}
#'   \item{\code{prop_speed_params}}{Vector, proposed speed parameters.}
#'   }
#'
#' @export update_speed_param
#'
update_speed_param <- function(fixed_constant, curr_move_params, behav_params, refined_path) {

  pert_sd <- fixed_constant$pert_sd
  # propose the new speed parameters from independent truncated normals
  curr_speed_params <- curr_move_params[-c(fixed_constant$q.bearings, fixed_constant$q.obs_error)]
  prop_speed_params <- rtnorm(fixed_constant$num_speed_param, mean = curr_speed_params, sd = pert_sd, lower = 0, upper = Inf)

  # ratio of the log proposal likelihoods (not symmetric because truncated distribution)
  prop_lik_ratio <- sum(dtnorm(curr_speed_params, mean = prop_speed_params, sd = pert_sd, lower = 0, upper = Inf, log = TRUE) -
                        dtnorm(prop_speed_params, mean = curr_speed_params, sd = pert_sd, lower = 0, upper = Inf, log = TRUE))

  # ratio of the log likelihood of the prior distribution given current and proposed params
  prop_move_params <- c(rep(NA, fixed_constant$num_states), prop_speed_params, NA)
  prior_lik_ratio <- calc_prior_speed_lik(move_params = prop_move_params) - calc_prior_speed_lik(move_params = curr_move_params)

  # calculate the log likelihood of the steps, given current and proposed parameters
  prop_step_lik <- calc_step_lik(fixed_constant, refined_path, move_params = prop_move_params, behav_params)
  curr_step_lik <- calc_step_lik(fixed_constant, refined_path, move_params = curr_move_params, behav_params)

  # calculate the standard MH acceptance ratio, taking exponential as all likelihoods are log likelihoods
  list(accept_prob = exp(prop_lik_ratio +  prior_lik_ratio + prop_step_lik - curr_step_lik), prop_speed_params = prop_speed_params)
}


#' @title Step likelihood.
#'
#' @description \code{calc_step_lik} returns the likelihood of a set of steps given speed parameters.
#'
#' @family Speed parameters
#'
#' @param fixed_constant List with components:
#'   \describe{
#'   \item{\code{num_states}}{Numeric, number of behavioural states.}
#'   \item{\code{indep_step}}{Boolean, True if model is independent step, False if it is correlated step.}
#'   \item{\code{q.speed_mean}}{Vector, index positions within \code{move_params} relating to speed mean parameters.}
#'   \item{\code{q.speed_corr}}{Vector, index positions within \code{move_params} relating to speed correlation parameters.}
#'   \item{\code{q.speed_var}}{Vector, index positions within \code{move_params} relating to speed variance parameters.}
#'   }
#' @param refined_path List with components:
#'   \describe{
#'   \item{\code{time_diffs}}{Vector, differences in the refined time scale.}
#'   \item{\code{behavs}}{Vector, refined behavioural/state sequence.}
#'   \item{\code{steps}}{Vector, refined step process.}
#'   }
#' @param move_params Vector of movement parameters (as turn_var_1, turn_var_2, mean_speed_1, .., speed_var_2).
#' @param behav_params List with components:
#'   \describe{
#'   \item{\code{lambda}}{Vector (length \code{num_states}) of switching rates.}
#'   \item{\code{q}}{Matrix (square size \code{num_states} with \code{NA} diagonal elements) of switching probabilities. Note: row = from, column = to.}
#'   }
#'
#' @return Numeric, likelihood of the set of steps.
#'
#' @example examples/calc_step_lik.R
#'
#' @export calc_step_lik
#'
calc_step_lik <- function(fixed_constant, refined_path, move_params, behav_params) {

  time_diffs <- refined_path$time_diffs; steps <- refined_path$steps; behavs <- refined_path$behavs # in single behav case "behavs" will just be null
  num_states <- fixed_constant$num_states; q.speed_mean <- fixed_constant$q.speed_mean; q.speed_var <- fixed_constant$q.speed_var

  # two different step models, independent and correlated, are handled differently
  if (fixed_constant$indep_step) {

    # likelihood of all the steps (initial has same distribution as the rest)
    # handle single and multi behaviour separately
    if (num_states == 1) {

      step_mean <- time_diffs * move_params[q.speed_mean]
      step_sd <- sqrt(time_diffs * move_params[q.speed_var])
    } else {

      step_mean <- time_diffs * move_params[q.speed_mean[behavs]]
      step_sd <- sqrt(time_diffs * move_params[q.speed_var[behavs]])
    }

    # all likelihoods are log likelihoods so sum
    return(sum(dnorm(steps, mean = step_mean, sd = step_sd, log = TRUE)))
  } else {

    # handle single and multi behaviour separately for the intial step likelihood
    if (num_states == 1) {

      initial_step_lik <- dnorm(steps[1], mean = time_diffs[1] * move_params[q.speed_mean], sd = time_diffs[1] * sqrt(move_params[q.speed_var]), log = TRUE)
    } else {

      # use equilibrium distribution of the behavioural process to calculate initial step likelihood as we don't know the initial behaviour
      steady_state <- calc_markov_steady_state(num_states, behav_params)$steady_state
      initial_step_per_state <- dnorm(steps[1], mean = time_diffs[1] * move_params[q.speed_mean], sd = time_diffs[1] * sqrt(move_params[q.speed_var]))
      initial_step_lik <- log(sum(initial_step_per_state * steady_state))
    }

    # likelihood of the rest of the steps can be gained by a conditional argument
    cond_step_dist <- calc_cond_step_dist(fixed_constant, steps, time_diffs, behavs, move_params)
    # check that there will not be a length mis-match leading to an incorrect recycling of values in the dnorm call
    if (length(steps[-1]) != length(cond_step_dist$mean) || length(steps[-1]) != length(cond_step_dist$sd)) stop('Conditional step length mis-match')
    cond_step_lik <- sum(dnorm(steps[-1], mean = cond_step_dist$mean, sd = cond_step_dist$sd, log = TRUE))

    return(initial_step_lik + cond_step_lik)
  }
}


#' @title Markov chain equilibrium distribution.
#'
#' @description \code{calc_markov_steady_state} returns the equilibrium distribution of a Markov chain.
#'
#' @family Speed parameters
#' @family Markov chain simulations
#'
#' @param num_states Numeric, number of behavioural states.
#' @param behav_params List with components:
#'   \describe{
#'   \item{\code{lambda}}{Vector (length \code{num_states}) of switching rates.}
#'   \item{\code{q}}{Matrix (square size \code{num_states} with \code{NA} diagonal elements) of switching probabilities. Note: row = from, column = to.}
#'   }
#'
#' @return List with components:
#'   \describe{
#'   \item{\code{steady_state}}{Vector (length \code{num_states}) giving the equilibrium distribution.}
#'   \item{\code{generator}}{Matrix (square, size \code{num_states}) giving the generator matrix of the Markov chain (just an alternative presentation of behav params).}
#'   }
#'
#' @example examples/calc_markov_steady_state.R
#'
#' @export calc_markov_steady_state
#'
calc_markov_steady_state <- function(num_states, behav_params) {

  # checks dimensions of the input parameters
  if (length(behav_params$lambda) != num_states) stop('Should have switch rates for each behavioural state')
  if (length(behav_params$q) != num_states^2) stop('Should have switch probs for each behavioural state')

  # create the generator matrix Q from the switching rates and probabilities
  Q <- matrix(NA, nrow = num_states, ncol = num_states)
  diag(Q) <- -behav_params$lambda
  for (i in 1:num_states) Q[i, -i] <- behav_params$q[i, -i] * rep(behav_params$lambda[i], num_states - 1)
  generator <- Q

  # solve the following to get the stationary distribution
  Q[ , num_states] <- rep(1, times = num_states)
  return(list(steady_state = as.vector(solve(t(Q), matrix(c(rep(0, times = num_states - 1), 1), ncol = 1))),
              generator = generator))
}


#' @title OU Conditional distribution.
#'
#' @description \code{calc_cond_step_dist} returns the mean and sd of a set of steps, conditional on each step before it.
#'
#' @family Speed parameters
#'
#' @param fixed_constant List with components:
#'   \describe{
#'   \item{\code{num_states}}{Numeric, number of behavioural states.}
#'   \item{\code{q.speed_mean}}{Vector, index positions within \code{move_params} relating to speed mean parameters.}
#'   \item{\code{q.speed_corr}}{Vector, index positions within \code{move_params} relating to speed correlation parameters.}
#'   \item{\code{q.speed_var}}{Vector, index positions within \code{move_params} relating to speed variance parameters.}
#'   }
#' @param steps Vector, refined step lengths.
#' @param time_diffs Vector, time period each step is over.
#' @param behavs Vector, refined behavioural states.
#' @param move_params Vector, parameters (arranged turn_var_1,turn_var_2,mean_speed_1,...,speed_corr_1,...,speed_var_1,...).
#'
#' @return List with componenets:
#'   \describe{
#'   \item{\code{mean}}{Vector of means.}
#'   \item{\code{sd}}{Vector of standard deviations.}
#'   }
#'
#' @example examples/calc_cond_step_dist.R
#'
#' @export calc_cond_step_dist
#'
calc_cond_step_dist <- function(fixed_constant, steps, time_diffs, behavs, move_params) {

  speed_mean <- move_params[fixed_constant$q.speed_mean]; speed_corr <- move_params[fixed_constant$q.speed_corr]; speed_var <- move_params[fixed_constant$q.speed_var]
  num_states <- fixed_constant$num_states

  if (num_states == 1) {

    speed_mean_vec <- speed_mean
    speed_corr_vec <- speed_corr
    speed_var_vec <- speed_var
  } else {

    # checks on the sizes of the input parameters
    if (length(steps) != length(behavs)) stop('Step and behaviour vectors not of equal length')
    if (length(speed_mean) != num_states || length(speed_corr) != num_states || length(speed_var) != num_states) stop('Should have parameter values for each behavioural state')
    # get the parameters for each step using the behaviour at each point
    speed_mean_vec <- speed_mean[behavs[-length(behavs)]]
    speed_corr_vec <- speed_corr[behavs[-length(behavs)]]
    speed_var_vec <- speed_var[behavs[-length(behavs)]]
  }

  # mean and variance are the conditional OU process, multiplied by time difference because we want steps not speeds
  expon_term <- exp(-speed_corr_vec * time_diffs[-length(time_diffs)])
  step_mean <- time_diffs[-1] * (speed_mean_vec + expon_term * (steps[-length(steps)] / time_diffs[-length(time_diffs)] - speed_mean_vec))
  step_sd <- time_diffs[-1] * sqrt(speed_var_vec * (1 - expon_term^2))

  return(list(mean = step_mean, sd = step_sd))
}


#' @title Simulate bearing parameters.
#'
#' @description \code{update_bearing_param} returns a conditional Gibbs sample of the bearing parameters.
#'
#' @family Bearing parameters
#'
#' @param fixed_constant List with components:
#'   \describe{
#'   \item{\code{num_states}}{Numeric, number of behavioural states.}
#'   \item{\code{turn_prec_prior}}{List with components:
#'     \describe{
#'     \item{\code{rate}}{Vector (length \code{num_states}) with prior rate of the turn precision (not variance).}
#'     \item{\code{shape}}{Vector (length \code{num_states}) with prior shape of the turn precision (not variance).}
#'     }}
#'   }
#' @param refined_path List with components:
#'   \describe{
#'   \item{\code{times}}{Vector of the refined times on the path.}
#'   \item{\code{behavs}}{Vector of refined behaviours.}
#'   \item{\code{bearings}}{Vector of refined bearings}
#'   }
#'
#' @return Vector, sampled bearing variances for each behavioural state.
#'
#' @example examples/update_bearing_param.R
#'
#' @export update_bearing_param
#'
update_bearing_param <- function(fixed_constant, refined_path) {

  bearings <- refined_path$bearings; times <- refined_path$times; behavs <- refined_path$behavs # `behavs' will be null in single behav case
  num_states <- fixed_constant$num_states

  # calculate the turns of the refined path and standardise them by time
  stand_turns <- diff(bearings) / sqrt(diff(times[-length(times)]))

  # handle single and multi behaviour cases separately
  if (num_states == 1) {

    # calculate posterior, note that prior is a global variable as it does not change
    shape_posterior <- fixed_constant$turn_prec_prior$shape + length(stand_turns) / 2
    rate_posterior <- fixed_constant$turn_prec_prior$rate + sum(stand_turns^2) / 2
    # sample from the posterior
    new_precision <- rgamma(1, shape = shape_posterior, rate = rate_posterior)
  } else {

    # checks on lengths of inputs to avoid mis-match recycling
    if (length(bearings) != length(behavs)) stop('Bearing and behaviour vectors not of equal length')
    if (length(bearings) != length(times[-length(times)])) stop('Times vector should have one more element that bearings')

    new_precision <- rep(NA, times = num_states)
    for (i in 1:num_states) { # for each state

      # take only the standardised turns corresponding to that state
      red_stand_turns <- stand_turns[behavs[-length(behavs)] == i]
      # calculate posterior, note that prior is a global variable as it does not change
      shape_posterior <- fixed_constant$turn_prec_prior$shape[i] + length(red_stand_turns) / 2
      rate_posterior <- fixed_constant$turn_prec_prior$rate[i] + sum(red_stand_turns^2) / 2
      # sample from the posterior
      new_precision[i] <- rgamma(1, shape = shape_posterior, rate = rate_posterior)
    }
  }

  # sampled precision, but want variance
  return(1 / new_precision)
}


#' @title Simulate observation error (independent) parameter.
#'
#' @description \code{Gibbs_update_obs_error_param} returns a conditional Gibbs sample of the observation error parameter when the errors are assumed independent.
#'
#' @family Observation error parameters
#'
#' @param error_prec_prior List with components:
#'   \describe{
#'     \item{\code{rate}}{Numeric with prior rate of the error precision (not variance).}
#'     \item{\code{shape}}{Numeric with prior shape of the error precision (not variance).}
#'   }
#' @param refined_path List with components:
#'   \describe{
#'   \item{\code{X}}{Vector of the refined x locations on the path.}
#'   \item{\code{Y}}{Vector of refined y locations on the path.}
#'   }
#' @param obs List with componenets:
#'   \describe{
#'   \item{\code{x}}{Observed x locations.}
#'   \item{\code{y}}{Observed y locations.}
#'   \item{\code{index_on_refined_path}}{Index values on the refined path where observations occur.}
#'   }
#'
#' @return Numeric, sampled error variance.
#'
#' @export Gibbs_update_obs_error_param
#'
Gibbs_update_obs_error_param <- function(error_prec_prior, refined_path, obs) {

  # calculate the actual error based on observations and refined path (put x and y all in one vector as it doesn't matter)
  obs_error <- c(refined_path$X[obs$index_on_refined_path] - obs$x, refined_path$Y[obs$index_on_refined_path] - obs$y)

  # calculate posterior, note that prior is a global variable as it does not change
  shape_posterior <- error_prec_prior$shape + length(obs_error) / 2
  rate_posterior <- error_prec_prior$rate + sum(obs_error^2) / 2
  # sample from the posterior
  new_precision <- rgamma(1, shape = shape_posterior, rate = rate_posterior)

  # sampled precision, but want variance
  return(1 / new_precision)
}

#' @title Simulate observation error (correlated) parameters.
#'
#' @description \code{MH_update_obs_error_param} returns an MH sample of the observation error parameters when the errors are assumed correlated.
#'
#' @family Observation error parameters
#'
#' @param fixed_constant List with components:
#'   \describe{
#'     \item{\code{q.obs_error}}{Vector with index positions in \code{move_params} that correspond to the error parameters.}
#'     \item{\code{error_pert_sd}}{Vector (length 2) of perturbation standard deviations of the error parameters.}
#'   }
#' @param obs List with componenets:
#'   \describe{
#'   \item{\code{errors}}{List with componenets:
#'     \describe{
#'     \item{\code{x}}{Vector of errors in x direction.}
#'     \item{\code{y}}{Vector of errors in y direction.}
#'     }}
#'   \item{\code{time_diffs}}{Vector of differences in the observation time scale.}
#'   }
#' @param curr_move_params Vector of parameters.
#'
#' @return Numeric, sampled error variance.
#'
#' @export MH_update_obs_error_param
#'
MH_update_obs_error_param <- function(fixed_constant, obs, curr_move_params) {

  error_pert_sd <- fixed_constant$error_pert_sd

  # propose the new obs error parameters from independent truncated normals
  curr_error_params <- curr_move_params[fixed_constant$q.obs_error]
  prop_error_params <- rtnorm(2, mean = curr_error_params, sd = error_pert_sd, lower = 0, upper = Inf)

  # ratio of the log proposal likelihoods (not symmetric because truncated distribution)
  prop_lik_ratio <- sum(dtnorm(curr_error_params, mean = prop_error_params, sd = error_pert_sd, lower = 0, upper = Inf, log = TRUE) -
                        dtnorm(prop_error_params, mean = curr_error_params, sd = error_pert_sd, lower = 0, upper = Inf, log = TRUE))

  # ratio of the log likelihood of the prior distribution given current and proposed params
  prior_lik_ratio <- calc_prior_error_lik(error_params = prop_error_params) - calc_prior_error_lik(error_params = curr_error_params)

  # calculate the log likelihood of the observation errors, given current and proposed parameters
  prop_error_lik <- calc_error_lik(obs, error_params = prop_error_params)
  curr_error_lik <- calc_error_lik(obs, error_params = curr_error_params)

  # calculate the standard MH acceptance ratio, taking exponential as all likelihoods are log likelihoods
  return(list(accept_prob = exp(prop_lik_ratio +  prior_lik_ratio + prop_error_lik - curr_error_lik), prop_error_params = prop_error_params))
}

#' @title Calculate observation error (correlated) likelihood
#'
#' @description \code{calc_error_lik} returns the likelihood of a set of errors (correlated).
#'
#' @family Observation error parameters
#'
#' @param obs List with componenets:
#'   \describe{
#'   \item{\code{errors}}{List with componenets:
#'     \describe{
#'     \item{\code{x}}{Vector of errors in x direction.}
#'     \item{\code{y}}{Vector of errors in y direction.}
#'     }}
#'   \item{\code{time_diffs}}{Vector of differences in the observation time scale.}
#'   }
#' @param error_params Vector (length 2) of the error variance and correlation.
#'
#' @return Numeric, likelihood.
#'
#' @export calc_error_lik
#'
calc_error_lik <- function(obs, error_params) {

  error_x <- obs$errors$x
  error_y <- obs$errors$y
  time_diffs <- obs$time_diffs

  # initial error likelihood, equilibrium OU process
  initial_error_lik <- sum(dnorm(c(error_x[1],error_y[1]), mean = 0, sd = error_params[1], log = TRUE))

  # rest of errors likelihood, conditional OU process
  expon_term <- exp(-error_params[2] * time_diffs)
  error_sd <- sqrt(error_params[1] * (1 - expon_term^2))
  cond_error_lik <- sum(dnorm(c(error_x[-1], error_y[-1]), c(error_x[-length(error_x)], error_y[-length(error_y)]) * expon_term, error_sd, log = TRUE))

  return(initial_error_lik + cond_error_lik)
}
