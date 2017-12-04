#' @title Step distribution of a path section (independent steps).
#'
#' @description \code{calc_step_dist} returns the distribution of a set of steps when the independent step model is assumed.
#'
#' @family Step distribution
#'
#' @param fixed_constant List with components:
#'   \describe{
#'   \item{\code{num_states}}{Numeric, number of behavioural states.}
#'   \item{\code{q.speed_mean}}{Vector, index values in \code{move_params} which correspond to mean speed parameters.}
#'   \item{\code{q.speed_var}}{Vector, index values in \code{move_params} which correspond to speed variance parameters.}
#'   }
#' @param times Vector of refined times.
#' @param behavs Vector of refined behaviours.
#' @param move_params Vector of movement parameters.
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{mean}}{Vector of mean step lengths.}
#'   \item{\code{sd}}{Vector of step length standard deviations.}
#'   \item{\code{var}}{Vector of step length variances.}
#'   }
#'
#' @export calc_step_dist
#'
calc_step_dist <- function(fixed_constant, times, behavs, move_params) {

  speed_mean <- move_params[fixed_constant$q.speed_mean]; speed_var <- move_params[fixed_constant$q.speed_var]
  num_states <- fixed_constant$num_states

  # get the time intervals of each step
  time_diffs <- diff(times)

  # handle single and mulit behaviour cases separately
  if (num_states == 1) {

    # mean and variance are BM with drift
    step_mean <- time_diffs * speed_mean
    step_var <- time_diffs * speed_var
    step_sd <- sqrt(step_var)
  } else {

    # checks on the sizes of the input parameters
    if (length(times) != length(behavs) + 1) stop('Times vector should have one more element that behavs')
    if (length(speed_mean) != num_states || length(speed_var) != num_states) stop('Should have parameter values for each behavioural state')

    step_mean <- time_diffs * speed_mean[behavs]
    step_var <- time_diffs * speed_var[behavs]
    step_sd <- sqrt(step_var)
  }

  return(list(mean = step_mean, sd = step_sd, var = step_var))
}


#' @title Step distribution of a path section (correlated steps, middle section).
#'
#' @description \code{calc_dist_step_bridge} returns the distribution of a set of steps when the correlated model is assumed and the section in the middle of the full path.
#'
#' @family Step distribution
#'
#' @param fixed_constant List with components:
#'   \describe{
#'   \item{\code{num_states}}{Numeric, number of behavioural states.}
#'   \item{\code{q.speed_mean}}{Vector, index values in \code{move_params} which correspond to mean speed parameters.}
#'   \item{\code{q.speed_corr}}{Vector, index values in \code{move_params} which correspond to speed correlation parameters.}
#'   \item{\code{q.speed_var}}{Vector, index values in \code{move_params} which correspond to speed variance parameters.}
#'   }
#' @param fixed_values List with components:
#'   \describe{
#'   \item{\code{times}}{Vector, length 2, of the times at the ends of the section.}
#'   \item{\code{after_time}}{Numeric, the time after the end point (used to work out the end speed).}
#'   \item{\code{behav}}{Numeric, fixed behaviour at the start.}
#'   \item{\code{steps}}{Vector, length 2, of the steps at the ends of the section.}
#'   }
#' @param times Vector of refined times.
#' @param behavs Vector of refined behaviours.
#' @param move_params Vector of movement parameters.
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{mean}}{Vector of mean step lengths.}
#'   \item{\code{covar}}{Square matrix of step length covariances.}
#'   }
#'
#' @export calc_dist_step_bridge
#'
calc_dist_step_bridge <- function(fixed_constant, fixed_values, times, behavs, move_params) {

  fixed_times <- fixed_values$times; after_end_time <- fixed_values$after_time; fixed_behav <- fixed_values$behav;  fixed_steps <- fixed_values$steps
  mean_para <- move_params[fixed_constant$q.speed_mean]; corr_para <- move_params[fixed_constant$q.speed_corr]; var_para <- move_params[fixed_constant$q.speed_var]

  # convert the fixed value steps to speeds
  fixed_speeds <- fixed_steps / c(times[1] - fixed_times[1], after_end_time - fixed_times[2])

  # call the OU bridge function with the fixed_speeds
  speed_bridge_dist <- calc_dist_OU_bridge(fixed_constant$num_states, fixed_times, fixed_behav, fixed_speeds, times, behavs, mean_para, corr_para, var_para)

  # convert the speed bridge to a step bridge
  time_diffs <- diff(c(times, fixed_times[2]))
  step_bridge_mean <- speed_bridge_dist$mean * time_diffs
  step_bridge_cov <- diag(time_diffs) %*% speed_bridge_dist$covar %*% diag(time_diffs)

  return(list(mean = step_bridge_mean, covar = step_bridge_cov))
}


#' @title Step distribution of a path section (correlated steps, start section).
#'
#' @description \code{calc_dist_step_backward} returns the distribution of a set of steps when the correlated model is assumed and the section at the start of the full path.
#'
#' @family Step distribution
#'
#' @param fixed_constant List with components:
#'   \describe{
#'   \item{\code{num_states}}{Numeric, number of behavioural states.}
#'   \item{\code{q.speed_mean}}{Vector, index values in \code{move_params} which correspond to mean speed parameters.}
#'   \item{\code{q.speed_corr}}{Vector, index values in \code{move_params} which correspond to speed correlation parameters.}
#'   \item{\code{q.speed_var}}{Vector, index values in \code{move_params} which correspond to speed variance parameters.}
#'   }
#' @param fixed_values List with components:
#'   \describe{
#'   \item{\code{time}}{Numeric, the time at the end of the section.}
#'   \item{\code{after_time}}{Numeric, the time after the end point (used to work out the end speed).}
#'   \item{\code{step}}{Numeric, the step at the end of the section.}
#'   }
#' @param times Vector of refined times.
#' @param behavs Vector of refined behaviours.
#' @param move_params Vector of movement parameters.
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{mean}}{Vector of mean step lengths.}
#'   \item{\code{covar}}{Square matrix of step length covariances.}
#'   }
#'
#' @export calc_dist_step_backward
#'
calc_dist_step_backward <- function(fixed_constant, fixed_values, times, behavs, move_params) {

  fixed_time <- fixed_values$time; after_end_time <- fixed_values$after_time; fixed_step <- fixed_values$step
  mean_para <- move_params[fixed_constant$q.speed_mean]; corr_para <- move_params[fixed_constant$q.speed_corr]; var_para <- move_params[fixed_constant$q.speed_var]

  # convert the fixed value step to speed
  fixed_speed <- fixed_step / (after_end_time - fixed_time)

  # call the OU backwards function with the fixed_speed
  speed_dist <- calc_dist_OU_backward(fixed_constant$num_states, fixed_time, fixed_speed, times, behavs, mean_para, corr_para, var_para)

  # convert the speeds to steps
  time_diffs <- diff(c(times, fixed_time))
  step_mean <- speed_dist$mean * time_diffs
  step_cov <- diag(time_diffs) %*% speed_dist$covar %*% diag(time_diffs)

  return(list(mean = step_mean, covar = step_cov))
}


#' @title Step distribution of a path section (correlated steps, end section).
#'
#' @description \code{calc_dist_step_forward} returns the distribution of a set of steps when the correlated model is assumed and the section at the end of the full path.
#'
#' @family Step distribution
#'
#' @param fixed_constant List with components:
#'   \describe{
#'   \item{\code{num_states}}{Numeric, number of behavioural states.}
#'   \item{\code{q.speed_mean}}{Vector, index values in \code{move_params} which correspond to mean speed parameters.}
#'   \item{\code{q.speed_corr}}{Vector, index values in \code{move_params} which correspond to speed correlation parameters.}
#'   \item{\code{q.speed_var}}{Vector, index values in \code{move_params} which correspond to speed variance parameters.}
#'   }
#' @param fixed_values List with components:
#'   \describe{
#'   \item{\code{time}}{Numeric, the time at the start of the section.}
#'   \item{\code{after_time}}{Numeric, the time after the end point (used to work out the end speed).}
#'   \item{\code{behav}}{Numeric, fixed behaviour at the start.}
#'   \item{\code{step}}{Numeric, the step at the start of the section.}
#'   }
#' @param times Vector of refined times.
#' @param behavs Vector of refined behaviours.
#' @param move_params Vector of movement parameters.
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{mean}}{Vector of mean step lengths.}
#'   \item{\code{covar}}{Square matrix of step length covariances.}
#'   }
#'
#' @export calc_dist_step_forward
#'
calc_dist_step_forward <- function(fixed_constant, fixed_values, times, behavs, move_params) {

  fixed_time <- fixed_values$time; after_end_time <- fixed_values$after_time; fixed_behav <- fixed_values$behav; fixed_step <- fixed_values$step
  mean_para <- move_params[fixed_constant$q.speed_mean]; corr_para <- move_params[fixed_constant$q.speed_corr]; var_para <- move_params[fixed_constant$q.speed_var]

  # convert the fixed value step to speed
  fixed_speed <- fixed_step / (times[1] - fixed_time)

  # call the OU forwards function with the fixed_speed
  speed_dist <- calc_dist_OU_forward(fixed_constant$num_states, fixed_time, fixed_behav, fixed_speed, times, behavs, mean_para, corr_para, var_para)

  # convert the speed bridge to a step bridge
  time_diffs <- diff(c(times, after_end_time))
  step_mean <- speed_dist$mean * time_diffs
  step_cov <- diag(time_diffs) %*% speed_dist$covar %*% diag(time_diffs)

  return(list(mean = step_mean, covar = step_cov))
}


#' @title Calculate OU bridge distribution.
#'
#' @description \code{calc_dist_OU_bridge} returns the distribution of an OU bridge.
#'
#' @family Step distribution
#'
#' @param num_states Numeric, number of behavioural states.
#' @param fixed_times Vector, length 2, of the times at the ends of the section.
#' @param fixed_behav Numeric, fixed behaviour at the start.
#' @param fixed_speeds Vector, length 2, of the speeds at the ends of the section.
#' @param times Vector of refined times.
#' @param behavs Vector of refined behaviours.
#' @param mean_para Vector of mean speed parameters.
#' @param corr_para Vector of speed correlation parameters.
#' @param var_para Vector of speed variance parameters.
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{mean}}{Vector of mean.}
#'   \item{\code{covar}}{Square matrix of covariances.}
#'   }
#'
#' @export calc_dist_OU_bridge
#'
calc_dist_OU_bridge <- function(num_states, fixed_times, fixed_behav, fixed_speeds, times, behavs, mean_para, corr_para, var_para) {

  # set up things
  num_points <- length(times)
  OU_mean <- rep(NA, times = num_points + 2)
  OU_covar <- matrix(NA, nrow = num_points + 2, ncol = num_points + 2)
  time_diffs <- diff(c(fixed_times[1], times, fixed_times[2]))

  # initialise start values (which are the known fixed values, hence 0 variance)
  OU_mean[1] <- fixed_speeds[1]
  OU_covar[1, ] <- OU_covar[ , 1] <- 0

  if (num_states == 1) {

    expon <- exp(-corr_para * time_diffs)
    for(j in 2:(num_points + 2)) {

      expon_t <- expon[j - 1]
      # iterative mean and variance, standard formula being followed
      OU_mean[j] <- OU_mean[j - 1] * expon_t + mean_para * (1 - expon_t)
      OU_covar[j, j] <- OU_covar[j - 1, j - 1] * expon_t^2 + var_para * (1 - expon_t^2)

      OU_covar[j, -(1:j)] <- OU_covar[-(1:j), j] <- cumprod(expon[-(1:(j-1))]) * OU_covar[j, j]
    }
  } else {

    all_behavs <- c(fixed_behav, behavs)
    expon <- exp(-corr_para[all_behavs] * time_diffs)
    # iterate through, filling in the mean and covariance of each entry based on the last, following an OU process (not bridge yet)
    for(j in 2:(num_points + 2)) {

      # set current params, use the behaviour from the previous time and the time difference between now and the previous point
      c_behav <- all_behavs[j - 1]
      expon_t <- expon[j - 1]

      # iterative mean and variance, standard formula being followed
      OU_mean[j] <- OU_mean[j - 1] * expon_t + mean_para[c_behav] * (1 - expon_t)
      OU_covar[j, j] <- OU_covar[j - 1, j - 1] * expon_t^2 + var_para[c_behav] * (1 - expon_t^2)

      OU_covar[j, -(1:j)] <- OU_covar[-(1:j), j] <- cumprod(expon[-(1:(j-1))]) * OU_covar[j, j]
    }
  }

  # separate the covariance matrix into the bit for the bridge (2,..,num_points+1) and the end bit to condition upon
  # want to keep sigma_12 as a matrix (it would be a single column so would automatically collapse to a vector),
  # but as sigma_22 is a single value, we want it to turn into a vector, so we can avoid some matrix multiplication
  sigma_12 <- OU_covar[2:(num_points + 1), num_points + 2, drop = FALSE]
  sigma_22 <- OU_covar[num_points + 2, num_points + 2]
  # calculate the bridge mean and covariance by conditioning upon the final value, this is by a standard formula
  bridge_mean <- OU_mean[2:(num_points + 1)] + (fixed_speeds[2] - OU_mean[num_points + 2]) / sigma_22 * sigma_12
  bridge_covar <- OU_covar[2:(num_points + 1), 2:(num_points + 1)] - (sigma_12 / sigma_22) %*% t(sigma_12)

  list(mean = bridge_mean, covar = bridge_covar)
}



#' @title Calculate reverse-time OU distribution.
#'
#' @description \code{calc_dist_OU_backward} returns the distribution of a reverse time OU process.
#'
#' @family Step distribution
#'
#' @param num_states Numeric, number of behavioural states.
#' @param fixed_time Numeric, of the time at the end of the section.
#' @param fixed_speed Numeric, of the speed at the end of the section.
#' @param times Vector of refined times.
#' @param behavs Vector of refined behaviours.
#' @param mean_para Vector of mean speed parameters.
#' @param corr_para Vector of speed correlation parameters.
#' @param var_para Vector of speed variance parameters.
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{mean}}{Vector of mean.}
#'   \item{\code{covar}}{Square matrix of covariances.}
#'   }
#'
#' @export calc_dist_OU_backward
#'
calc_dist_OU_backward <- function(num_states, fixed_time, fixed_speed, times, behavs, mean_para, corr_para, var_para) {

  # set up things
  num_points <- length(times)
  OU_mean <- rep(NA, times = num_points + 1)
  OU_covar <- matrix(NA, nrow = num_points + 1, ncol = num_points + 1)
  time_diffs <- diff(c(times, fixed_time))


  # initialise start values (which are the known fixed values, hence 0 variance)
  OU_mean[num_points + 1] <- fixed_speed
  OU_covar[num_points + 1, ] <- OU_covar[ , num_points + 1] <- 0

  if (num_states == 1) {

    expon <- exp(-corr_para * time_diffs)
    # iterate through, backwards, filling in the mean and covariance of each entry based on the last, following a time reversed OU process
    for(j in num_points:1) {

      expon_t <- expon[j]
      # iterative mean and variance, standard formula being followed
      OU_mean[j] <- OU_mean[j + 1] * expon_t + mean_para * (1 - expon_t)
      OU_covar[j, j] <- OU_covar[j + 1, j + 1] * expon_t^2 + var_para * (1 - expon_t^2)

      if (j != 1) {
        OU_covar[j, (j-1):1] <- OU_covar[(j-1):1, j] <- cumprod(expon[(j-1):1]) * OU_covar[j, j]
      }
    }
  } else {

    expon <- exp(-corr_para[behavs] * time_diffs)
    # iterate through, backwards, filling in the mean and covariance of each entry based on the last, following a time reversed OU process
    for(j in num_points:1) {

      # set current params, use the behaviour from the current point, and the time period between now and the next point (because backwards)
      c_behav <- behavs[j]
      expon_t <- expon[j]

      # iterative mean and variance, standard formula being followed
      OU_mean[j] <- OU_mean[j + 1] * expon_t + mean_para[c_behav] * (1 - expon_t)
      OU_covar[j, j] <- OU_covar[j + 1, j + 1] * expon_t^2 + var_para[c_behav] * (1 - expon_t^2)

      if (j != 1) {
        OU_covar[j, (j-1):1] <- OU_covar[(j-1):1, j] <- cumprod(expon[(j-1):1]) * OU_covar[j, j]
      }
    }
  }

  list(mean = OU_mean[-(num_points + 1)], covar = OU_covar[-(num_points + 1), -(num_points + 1)])
}



#' @title Calculate OU distribution.
#'
#' @description \code{calc_dist_OU_forward} returns the distribution of an OU prcoess
#'
#' @family Step distribution
#'
#' @param num_states Numeric, number of behavioural states.
#' @param fixed_time Numeric, of the time at the start of the section.
#' @param fixed_behav Numeric, fixed behaviour at the start.
#' @param fixed_speed Numeric, of the speed at the start of the section.
#' @param times Vector of refined times.
#' @param behavs Vector of refined behaviours.
#' @param mean_para Vector of mean speed parameters.
#' @param corr_para Vector of speed correlation parameters.
#' @param var_para Vector of speed variance parameters.
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{mean}}{Vector of mean.}
#'   \item{\code{covar}}{Square matrix of covariances.}
#'   }
#'
#' @export calc_dist_OU_forward
#'
calc_dist_OU_forward <- function(num_states, fixed_time, fixed_behav, fixed_speed, times, behavs, mean_para, corr_para, var_para) {

  # set up things
  num_points <- length(times)
  OU_mean <- rep(NA, times = num_points + 1)
  OU_covar <- matrix(NA, nrow = num_points + 1, ncol = num_points + 1)
  time_diffs <- diff(c(fixed_time, times))

  # initialise start values (which are the known fixed values, hence 0 variance)
  OU_mean[1] <- fixed_speed
  OU_covar[1, ] <- OU_covar[ , 1] <- 0

  if (num_states == 1) {

    expon <- exp(-corr_para * time_diffs)
    # iterate through, filling in the mean and covariance of each entry based on the last, following an OU process
    for(j in 2:(num_points + 1)) {

      expon_t <- expon[j-1]
      # iterative mean and variance, standard formula being followed
      OU_mean[j] <- OU_mean[j - 1] * expon_t + mean_para * (1 - expon_t)
      OU_covar[j, j] <- OU_covar[j - 1, j - 1] * expon_t^2 + var_para * (1 - expon_t^2)

      OU_covar[j, -(1:j)] <- OU_covar[-(1:j), j] <- cumprod(expon[-(1:(j-1))]) * OU_covar[j, j]
    }
  } else {

    all_behavs <- c(fixed_behav, behavs[-length(behavs)])
    expon <- exp(-corr_para[all_behavs] * time_diffs)
    # iterate through, filling in the mean and covariance of each entry based on the last, following an OU process
    for(j in 2:(num_points + 1)) {

      # set current params, use the behaviour from the previous time and the time difference between now and the previous point
      c_behav <- all_behavs[j - 1]
      expon_t <- expon[j-1]

      # iterative mean and variance, standard formula being followed
      OU_mean[j] <- OU_mean[j - 1] * expon_t + mean_para[c_behav] * (1 - expon_t)
      OU_covar[j, j] <- OU_covar[j - 1, j - 1] * expon_t^2 + var_para[c_behav] * (1 - expon_t^2)

      OU_covar[j, -(1:j)] <- OU_covar[-(1:j), j] <- cumprod(expon[-(1:(j-1))]) * OU_covar[j, j]
    }
  }

  return(list(mean = OU_mean[-1], covar = OU_covar[-1, -1]))
}
