# this needs changing to handle multiple behaviours and independent error and no error cases need adding
InitialPath <- function(obs, time_scale, error = FALSE, corr_error = FALSE, error_para = NULL) {

  refined_path <- list(times = seq(obs$times[1],obs$times[length(obs$times)],time_scale))
  refined_path$time_diffs <- diff(refined_path$times)

  if(error) {
    if(corr_error) {
      error_x <- error_y <- rep(NA, length(obs$x))
      error_x[1] <- rnorm(1, 0, sqrt(error_para[1]))
      error_y[1] <- rnorm(1, 0, sqrt(error_para[1]))
      for(i in 2:length(error_x)) {
        error_x[i] <- rnorm(1, error_x[i-1] * exp(-error_para[2] * (obs$times[i]-obs$times[i-1])), sqrt(error_para[1]*(1-exp(-2*error_para[2]*(obs$times[i]-obs$times[i-1])))))
        error_y[i] <- rnorm(1, error_y[i-1] * exp(-error_para[2] * (obs$times[i]-obs$times[i-1])), sqrt(error_para[1]*(1-exp(-2*error_para[2]*(obs$times[i]-obs$times[i-1])))))
      }
      refined_path$X <- spline(obs$times, obs$x + error_x, n = length(refined_path$times), method="natural")$y
      refined_path$Y <- spline(obs$times, obs$y + error_y, n = length(refined_path$times), method="natural")$y
    } else {
      error_x <- error_y <- rep(NA, length(obs$x))
      error_x <- rnorm(1, 0, sqrt(error_para))
      error_y <- rnorm(1, 0, sqrt(error_para))
      refined_path$X <- spline(obs$times, obs$x + error_x, n = length(refined_path$times), method="natural")$y
      refined_path$Y <- spline(obs$times, obs$y + error_y, n = length(refined_path$times), method="natural")$y
    }
  } else {
    refined_path$X <- spline(obs$times, obs$x, n = length(refined_path$times), method="natural")$y
    refined_path$Y <- spline(obs$times, obs$y, n = length(refined_path$times), method="natural")$y
  }

  refined_path$steps <- sqrt(diff(refined_path$Y)^2 + diff(refined_path$X)^2)
  refined_path$bearings <- atan2(diff(refined_path$Y), diff(refined_path$X))
  for(i in which(diff(refined_path$bearings) > pi)){refined_path$bearings <- c(refined_path$bearings[1:i], refined_path$bearings[(i+1):length(refined_path$bearings)] - 2 * pi)}
  for(i in which(diff(refined_path$bearings) < (-pi))){refined_path$bearings <- c(refined_path$bearings[1:i], refined_path$bearings[(i+1):length(refined_path$bearings)] + 2 * pi)}

  return(refined_path)
}
