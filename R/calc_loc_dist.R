#' @title Calculate distribution of fixed locations.
#'
#' @description \code{calc_dist_true_locs_forward} returns the distribution of a set of locations, given step distribution and actual bearings.
#'
#' @family Location distribution
#'
#' @param start_loc Vector (length 2) of the x,y starting position of the section to update.
#' @param fixed_index Vector (variable length, at least 1) of the index points on the section where a fixed location is (will always include the last point).
#' @param bearings Vector (variable length) of the actual bearings on the section.
#' @param step_dist List with components:
#'   \describe{
#'   \item{\code{mean}}{Vector of the mean step lengths on the section.}
#'   \item{\code{covar}}{Square matrix of the covariances between the steps on the section.}
#'   }
#'
#' @return List with components:
#'   \describe{
#'   \item{\code{A}}{Linear constraint matrix.}
#'   \item{\code{loc_mean}}{Vector of location means (arranged as x_1,y_1,x_2,...).}
#'   \item{\code{loc_covar}}{Square matrix of location covariances.}
#'   \item{\code{loc_step_covar}}{Matrix of covariances between the steps and the locations.}
#'   }
#'
#' @export calc_dist_true_locs_forward
#'
calc_dist_true_locs_forward <- function(start_loc, fixed_index, bearings, step_dist) {

  step_mean <- step_dist$mean; step_covar <- step_dist$covar
  num_points <- length(fixed_index)
  sec_length <- length(bearings)

  # create the linear constraint matrix
  A <- matrix(0, nrow = 2 * num_points, ncol = sec_length)
  for (i in 1:num_points) {
    A[2 * i - 1, 1:fixed_index[i]] <- cos(bearings[1:fixed_index[i]])
    A[2 * i,     1:fixed_index[i]] <- sin(bearings[1:fixed_index[i]])
  }

  # mean of the locs
  loc_mean <- start_loc + A %*% step_mean
  # cov of steps and locs
  loc_step_covar <- step_covar %*% t(A)
  # cov of locs
  loc_covar <- A %*% loc_step_covar

  return(list(A = A, loc_mean = loc_mean, loc_covar = loc_covar, loc_step_covar = loc_step_covar))
}


#' @title Calculate distribution of fixed locations (but in reverse time).
#'
#' @description \code{calc_dist_true_locs_backward} returns the distribution of a set of locations, given step distribution and actual bearings, in reverse-time.
#'
#' @family Location distribution
#'
#' @param start_loc Vector (length 2) of the x,y ending position of the section to update.
#' @param fixed_index Vector (variable length, at least 1) of the index points on the section where a fixed location is (will always include the initial location).
#' @param bearings Vector (variable length) of the actual bearings on the section.
#' @param step_dist List with components:
#'   \describe{
#'   \item{\code{mean}}{Vector of the mean step lengths on the section.}
#'   \item{\code{covar}}{Square matrix of the covariances between the steps on the section.}
#'   }
#'
#' @return List with components:
#'   \describe{
#'   \item{\code{A}}{Linear constraint matrix.}
#'   \item{\code{loc_mean}}{Vector of location means (arranged as x_1,y_1,x_2,...).}
#'   \item{\code{loc_covar}}{Square matrix of location covariances.}
#'   \item{\code{loc_step_covar}}{Matrix of covariances between the steps and the locations.}
#'   }
#'
#' @export calc_dist_true_locs_backward
#'
calc_dist_true_locs_backward <- function(start_loc, fixed_index, bearings, step_dist) {

  step_mean <- step_dist$mean; step_covar <- step_dist$covar
  num_points <- length(fixed_index)
  sec_length <- length(bearings)

  # create the linear constraint matrix
  A <- matrix(0, nrow = 2 * num_points, ncol = sec_length)
  for (i in 1:num_points) {
    A[2 * i - 1, fixed_index[i]:sec_length] <- -cos(bearings[fixed_index[i]:sec_length])
    A[2 * i,     fixed_index[i]:sec_length] <- -sin(bearings[fixed_index[i]:sec_length])
  }

  # mean of the locs
  loc_mean <- start_loc + A %*% step_mean
  # cov of steps and locs
  loc_step_covar <- step_covar %*% t(A)
  # cov of locs
  loc_covar <- A %*% loc_step_covar

  return(list(A = A, loc_mean = loc_mean, loc_covar = loc_covar, loc_step_covar = loc_step_covar))
}

calc_dist_noisy_locs_bridge <- function(fixed_values, error_params, loc_dist) {

  contained_obs_times <- fixed_values$obs_times; fixed_error_times <- fixed_values$error_times; fixed_errors <- fixed_values$errors
  loc_mean <- loc_dist$loc_mean; loc_covar <- loc_dist$loc_covar

  # set up things
  num_points <- length(contained_obs_times)
  if (num_points == 0) {return(list(loc_mean_error = loc_mean, loc_covar_error = loc_covar))} # no internal observations to add error to
  OU_mean_x <- OU_mean_y <- rep(NA, times = num_points + 2)
  OU_covar <- matrix(NA, nrow = num_points + 2, ncol = num_points + 2)
  time_diffs <- diff(c(fixed_error_times[1], contained_obs_times, fixed_error_times[2]))

  # initialise start values (which are the known fixed values, hence 0 variance)
  OU_mean_x[1] <- fixed_errors$x[1]
  OU_mean_y[1] <- fixed_errors$y[1]
  OU_covar[1, ] <- OU_covar[ , 1] <- 0

  expon <- exp(-error_params[2] * time_diffs)
  for(j in 2:(num_points + 2)) {

    expon_t <- expon[j - 1]
    # iterative mean and variance, standard formula being followed
    OU_mean_x[j] <- OU_mean_x[j - 1] * expon_t
    OU_mean_y[j] <- OU_mean_y[j - 1] * expon_t
    OU_covar[j, j] <- OU_covar[j - 1, j - 1] * expon_t^2 + error_params[1] * (1 - expon_t^2)
    OU_covar[j, -(1:j)] <- OU_covar[-(1:j), j] <- cumprod(expon[-(1:(j-1))]) * OU_covar[j, j]
  }

  # separate the covariance matrix into the bit for the bridge (2,..,num_points+1) and the end bit to condition upon
  # want to keep sigma_12 as a matrix (it would be a single column so would automatically collapse to a vector),
  # but as sigma_22 is a single value, we want it to turn into a vector, so we can avoid some matrix multiplication
  sigma_12 <- OU_covar[2:(num_points + 1), num_points + 2, drop = FALSE]
  sigma_22 <- OU_covar[num_points + 2, num_points + 2, drop = TRUE]
  # calculate the bridge mean and covariance by conditioning upon the final value, this is by a standard formula
  bridge_mean_x <- OU_mean_x[2:(num_points + 1)] + (fixed_errors$x[2] - OU_mean_x[num_points + 2]) / sigma_22 * sigma_12
  bridge_mean_y <- OU_mean_y[2:(num_points + 1)] + (fixed_errors$y[2] - OU_mean_y[num_points + 2]) / sigma_22 * sigma_12
  bridge_covar <- OU_covar[2:(num_points + 1), 2:(num_points + 1)] - (sigma_12 / sigma_22) %*% t(sigma_12)

  # end location has to meet the end point exactly, so there are zero errors there
  error_mean <- c(rbind(as.numeric(bridge_mean_x), as.numeric(bridge_mean_y)), 0, 0)
  loc_mean_error <- loc_mean + error_mean
  covar_error <- matrix(0, nrow = 2 * nrow(bridge_covar) + 2, ncol = 2 * nrow(bridge_covar) + 2)
  for (i in 1:(num_points)) {
    covar_error[2*i-1, 1:(2*nrow(bridge_covar))] <- c(rbind(bridge_covar[i, ], 0))
    covar_error[2*i, 1:(2*nrow(bridge_covar))]  <- c(rbind(0, bridge_covar[i, ]))
  }
  loc_covar_error <- loc_covar + covar_error

  return(list(loc_mean_error = loc_mean_error, loc_covar_error = loc_covar_error, error_mean = error_mean, error_covar = covar_error))
}


calc_dist_noisy_locs_forward <- function(fixed_values, error_params, loc_dist, obs) {

  contained_obs_times <- c(fixed_values$obs_times, obs$times[length(obs$times)])
  fixed_error_times <- fixed_values$error_times; fixed_errors <- fixed_values$errors
  loc_mean <- loc_dist$loc_mean; loc_covar <- loc_dist$loc_covar

  # set up things
  num_points <- length(contained_obs_times)
  OU_mean_x <- OU_mean_y <- rep(NA, times = num_points + 1)
  OU_covar <- matrix(NA, nrow = num_points + 1, ncol = num_points + 1)
  time_diffs <- diff(c(fixed_error_times, contained_obs_times))

  # initialise start values (which are the known fixed values, hence 0 variance)
  OU_mean_x[1] <- fixed_errors$x
  OU_mean_y[1] <- fixed_errors$y
  OU_covar[1, ] <- OU_covar[ , 1] <- 0

  expon <- exp(-error_params[2] * time_diffs)
  for(j in 2:(num_points + 1)) {

    expon_t <- expon[j - 1]
    # iterative mean and variance, standard formula being followed
    OU_mean_x[j] <- OU_mean_x[j - 1] * expon_t
    OU_mean_y[j] <- OU_mean_y[j - 1] * expon_t
    OU_covar[j, j] <- OU_covar[j - 1, j - 1] * expon_t^2 + error_params[1] * (1 - expon_t^2)
    OU_covar[j, -(1:j)] <- OU_covar[-(1:j), j] <- cumprod(expon[-(1:(j-1))]) * OU_covar[j, j]
  }

  mean_x <- OU_mean_x[-1]
  mean_y <- OU_mean_y[-1]
  covar <- OU_covar[-1, -1, drop = FALSE]

  error_mean <- c(rbind(as.numeric(mean_x), as.numeric(mean_y)))
  loc_mean_error <- loc_mean + error_mean
  covar_error <- matrix(0, nrow = 2 * nrow(covar), ncol = 2 * nrow(covar))
  for (i in 1:(num_points)) {
    covar_error[2*i-1, 1:(2*nrow(covar))] <- c(rbind(covar[i, ], 0))
    covar_error[2*i, 1:(2*nrow(covar))]  <- c(rbind(0, covar[i, ]))
  }
  loc_covar_error <- loc_covar + covar_error

  return(list(loc_mean_error = loc_mean_error, loc_covar_error = loc_covar_error, error_mean = error_mean, error_covar = covar_error))
}


calc_dist_noisy_locs_backward <- function(fixed_values, error_params, loc_dist, obs) {

  contained_obs_times <- c(obs$times[1], fixed_values$obs_times)
  fixed_error_times <- fixed_values$error_times; fixed_errors <- fixed_values$errors
  loc_mean <- loc_dist$loc_mean; loc_covar <- loc_dist$loc_covar

  # set up things
  num_points <- length(contained_obs_times)
  OU_mean_x <- OU_mean_y <- rep(NA, times = num_points + 1)
  OU_covar <- matrix(NA, nrow = num_points + 1, ncol = num_points + 1)
  time_diffs <- diff(c(contained_obs_times, fixed_error_times))

  # initialise start values (which are the known fixed values, hence 0 variance)
  OU_mean_x[num_points + 1] <- fixed_errors$x
  OU_mean_y[num_points + 1] <- fixed_errors$y
  OU_covar[num_points + 1, ] <- OU_covar[ , num_points + 1] <- 0

  expon <- exp(-error_params[2] * time_diffs)
  # iterate through, backwards, filling in the mean and covariance of each entry based on the last, following a time reversed OU process
  for(j in num_points:1) {

    expon_t <- expon[j]
    # iterative mean and variance, standard formula being followed
    OU_mean_x[j] <- OU_mean_x[j + 1] * expon_t
    OU_mean_y[j] <- OU_mean_y[j + 1] * expon_t
    OU_covar[j, j] <- OU_covar[j + 1, j + 1] * expon_t^2 + error_params[1] * (1 - expon_t^2)

    if (j != 1) {
      OU_covar[j, (j-1):1] <- OU_covar[(j-1):1, j] <- cumprod(expon[(j-1):1]) * OU_covar[j, j]
    }
  }

  mean_x <- OU_mean_x[-(num_points + 1)]
  mean_y <- OU_mean_y[-(num_points + 1)]
  covar <- OU_covar[-(num_points + 1), -(num_points + 1), drop = FALSE]

  error_mean <- c(rbind(as.numeric(mean_x), as.numeric(mean_y)))
  loc_mean_error <- loc_mean + error_mean
  covar_error <- matrix(0, nrow = 2 * nrow(covar), ncol = 2 * nrow(covar))
  for (i in 1:(num_points)) {
    covar_error[2*i-1, 1:(2*nrow(covar))] <- c(rbind(covar[i, ], 0))
    covar_error[2*i, 1:(2*nrow(covar))]  <- c(rbind(0, covar[i, ]))
  }
  loc_covar_error <- loc_covar + covar_error

  return(list(loc_mean_error = loc_mean_error, loc_covar_error = loc_covar_error, error_mean = error_mean, error_covar = covar_error))
}
