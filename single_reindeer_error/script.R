# Set a seed
set.seed(1234)
# Packages required
library(gtools); library(msm); library(mvnfast); library(fdrtool)
library(CTStepTurn)
source("../FixedConstants.R"); source("../InitialPath.R")

# Set observations (x, y, times)
obs <- as.list(read.table("observations.txt", header = T))
obs$time_diffs <- diff(obs$times)

# Set initial refined path
refined_path <- InitialPath(obs, time_scale = 0.25, error = TRUE, corr_error = FALSE, error_para = 100)
write.table(list(times = refined_path$times,
                 bearings = c(refined_path$bearings, NA),
                 steps = c(refined_path$steps, NA),
                 X = refined_path$X,
                 Y = refined_path$Y), "initial_path.txt", row.names = F)
obs$index_on_refined_path <- which(refined_path$times %in% obs$times)
obs$errors <- list(x = obs$x - refined_path$X[obs$index_on_refined_path],
                   y = obs$y - refined_path$Y[obs$index_on_refined_path])

# Set initial parameters - also need behaviour parameters if applicable
move_params <-  c(var(diff(refined_path$bearings))/0.25,
                  mean(refined_path$steps)/0.25,
                  log(acf(refined_path$steps,plot=F)$acf[2])/(-2*0.25),
                  var(refined_path$steps)/0.25^2,
                  100)

# Assign the fixed constants
fixed_constant <- FixedConstants(pert_sd = sqrt(c(4, 0.3, 100)), error_pert_sd = NULL,
                                 num_states = 1,
                                 turn_prec_prior_shape = 0.5, turn_prec_prior_rate = 0.5,
                                 error_prec_prior_shape = 2, error_prec_prior_rate = 200,
                                 indep_step = FALSE, obs_error = TRUE, corr_obs_error = FALSE)

# Assign fixed prior distributions - also need behaviour priors if applicable
calc_prior_speed_lik <- function(move_params) {
  neg_step <- ifelse(pnorm(0,
                           mean = 5 * exp(-move_params[3]*0.25) + move_params[2]*(1-exp(-move_params[3]*0.25)),
                           sd = sqrt(move_params[4]*(1-exp(-2*move_params[3]*0.25)))) < 0.1, log(1), log(0))
  beta_prior <- dhalfnorm(move_params[3], theta = sqrt(0.7), log = TRUE)
  neg_step + beta_prior
}

# Assign variables relating to the number of MCMC iterations, thinning, etc.
section_length_limits <- c(3, 16)
speed_param_count <- 0; refined_path_count <- 0
num_iterations <- 5 * 10^8; num_extra_path_updates <- 50; thin <- 10^5
stored_move_params <- matrix(NA, nrow = num_iterations / thin, ncol = fixed_constant$num_move_params)
stored_refined_path <- vector("list", num_iterations / thin)

####### Carry out the full MCMC loop ##############################

# loop over the number of path updates
for (i in 1:num_iterations) {

  # but update parameters less often than the path
  if (i %% num_extra_path_updates == 0) {

    # propose an update to the speed parameters in a random walk MH
    speed_param_proposal <- update_speed_param(fixed_constant, curr_move_params = move_params, behav_params = NULL, refined_path)
    # decide whether to accept the MH proposal
    if (runif(1) < speed_param_proposal$accept_prob) {

      speed_param_count <- speed_param_count + 1 # keep track of acceptance rate
      move_params[-c(fixed_constant$q.bearings, fixed_constant$q.obs_error)] <- speed_param_proposal$prop_speed_params
    }

    # update the bearing parameters by a Gibbs sampler
    move_params[fixed_constant$q.bearings] <- update_bearing_param(fixed_constant, refined_path)

    # carry out a Gibbs update to the observation error parameter
    move_params[fixed_constant$q.obs_error] <- Gibbs_update_obs_error_param(fixed_constant$error_prec_prior, refined_path, obs)

  }

  # update the refined path (in sections, chosen randomly)
  # flag if it is a special end case and pick all relevent info out from the full path
  drawn_section <- draw_random_path_section(section_length_limits, obs_index = obs$index_on_refined_path)
  type <- if (drawn_section$start == 1) "start" else if (drawn_section$end == length(refined_path$bearings)) "end"  else  "middle"
  fixed_values <- set_fixed_values(type, drawn_section, refined_path, obs)
  curr_values <- set_current_values(type, drawn_section, refined_path)

  # update the path section
  refined_path_proposal <- update_refined_path(fixed_constant, type, fixed_values, curr_values, move_params, behav_params = NULL, obs)
  # update values in the full path if it was accepted
  if (refined_path_proposal$accept == 1) {

    refined_path_count <- refined_path_count + 1 # keep track of acceptance rate

    start_to_end <- drawn_section$start:drawn_section$end
    refined_path$bearings[start_to_end] <- refined_path_proposal$bearings
    refined_path$steps[start_to_end] <- refined_path_proposal$steps
    obs$errors <- list(x = obs$x - refined_path$X[obs$index_on_refined_path], y = obs$y - refined_path$Y[obs$index_on_refined_path])
    if (type == "middle") {
        
      one_to_start <- 1:(drawn_section$start)
      one_to_end <- 1:(drawn_section$end + 1)
      refined_path$X <- c(refined_path$X[one_to_start], refined_path$X[drawn_section$start] + cumsum(refined_path_proposal$steps * cos(refined_path_proposal$bearings)), refined_path$X[-one_to_end])
      refined_path$Y <- c(refined_path$Y[one_to_start], refined_path$Y[drawn_section$start] + cumsum(refined_path_proposal$steps * sin(refined_path_proposal$bearings)), refined_path$Y[-one_to_end])
    } else if (type == "start") {
        
      one_to_end <- 1:(drawn_section$end)
      refined_path$X <- c(rev(refined_path$X[drawn_section$end + 1] - cumsum(rev(refined_path_proposal$steps * cos(refined_path_proposal$bearings)))), refined_path$X[-one_to_end])
      refined_path$Y <- c(rev(refined_path$Y[drawn_section$end + 1] - cumsum(rev(refined_path_proposal$steps * sin(refined_path_proposal$bearings)))), refined_path$Y[-one_to_end])
      # initial bearing is uniform on [-pi,pi] so need to thether it to that interval
      while(refined_path$bearings[1] > pi) {refined_path$bearings <- refined_path$bearings - 2 * pi}
      while(refined_path$bearings[1] < (-pi)) {refined_path$bearings <- refined_path$bearings + 2 * pi}
    } else if (type == "end") {
        
      one_to_start <- 1:(drawn_section$start)
      refined_path$X <- c(refined_path$X[one_to_start], refined_path$X[drawn_section$start] + cumsum(refined_path_proposal$steps * cos(refined_path_proposal$bearings)))
      refined_path$Y <- c(refined_path$Y[one_to_start], refined_path$Y[drawn_section$start] + cumsum(refined_path_proposal$steps * sin(refined_path_proposal$bearings)))
    }
    
  # Store a thinned sample of the full system
  if (i %% thin == 0) {

    stored_move_params[i / thin, ] <- move_params
    stored_refined_path[[i / thin]] <- refined_path
  }
}
############################################

save.image("results.RData")
write.table(stored_move_params, file = "move_param.txt", row.names = F)
write.table(c(speed_param_count, speed_param_count / (num_iterations / num_extra_path_updates),
              refined_path_count, refined_path_count / num_iterations), file = "acceptance.txt", row.names = F)
n <- num_iterations / thin
l <- length(refined_path$X)
samp_bearings <- matrix(NA, nrow = n, ncol = length(refined_path$bearings))
samp_steps<- matrix(NA, nrow = n, ncol = length(refined_path$bearings))
samp_loc <- matrix(NA, nrow = n*l, ncol = 3)
for(i in 1:n){
  samp_bearings[i, ] <- stored_refined_path[[i]]$bearings
  samp_steps[i, ] <- stored_refined_path[[i]]$steps
  samp_loc[((i-1)*l+1):(i*l),1] <- stored_refined_path[[i]]$X
  samp_loc[((i-1)*l+1):(i*l),2] <- stored_refined_path[[i]]$Y
  samp_loc[((i-1)*l+1):(i*l),3] <- i
}
write.table(samp_bearings, "bearings.txt", row.names = F)
write.table(samp_steps, "steps.txt", row.names = F)
write.table(samp_loc, "loc.txt", row.names = F)
