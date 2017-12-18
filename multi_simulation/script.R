# Set a seed
set.seed(123456)
# Packages required
library(gtools); library(msm); library(mvnfast); library(fdrtool)
library(CTStepTurn)
source("../FixedConstants.R"); source("../InitialPath.R")

# Set observations (x, y, times)
obs <- as.list(read.table("observations.txt", header = T))
obs$time_diffs <- diff(obs$times)

# Set initial refined path
refined_path <- InitialPath(obs, time_scale = 1, error = FALSE, corr_error = FALSE, error_para = 100)
refined_path$behavs <- rep(1, length(refined_path$bearings))
refined_path$behavs[refined_path$steps > 20] <- 2
obs$index_on_refined_path <- which(refined_path$times %in% obs$times)
behav_proc <- create_behav_proc_from_refined(refined_path$times,refined_path$behavs)

# Set initial parameters - also need behaviour parameters if applicable
move_params <- c(var(diff(refined_path$bearings)[refined_path$behavs==1], na.rm = T)/1,
                 var(diff(refined_path$bearings)[refined_path$behavs==2], na.rm = T)/1,
                 mean(refined_path$steps[refined_path$behavs==1], na.rm = T)/1,
                 mean(refined_path$steps[refined_path$behavs==2], na.rm = T)/1,
                 log(acf(refined_path$steps[refined_path$behavs==1],plot=F)$acf[2])/(-2*1),
                 log(acf(refined_path$steps[refined_path$behavs==2],plot=F)$acf[2])/(-2*1),
                 var(refined_path$steps[refined_path$behavs==1], na.rm = T)/1^2,
                 var(refined_path$steps[refined_path$behavs==2], na.rm = T)/1^2)
behav_params <- list(lambda = c(0.01, 0.02),
                     q = matrix(c(NA,1,1,NA),nrow=2,ncol=2))

# Assign the fixed constants
fixed_constant <- FixedConstants(pert_sd = c(0.13,0.18,0.017,0.014,2.7,2.1), error_pert_sd = NULL,
                                 num_states = 2,
                                 turn_prec_prior_shape = c(0.001,0.001), turn_prec_prior_rate = c(0.001,0.001),
                                 indep_step = FALSE, obs_error = FALSE, corr_obs_error = FALSE,
                                 time_scale = 1,
                                 prior_param = list(lambda_shape = c(0.001,0.001),
                                                    lambda_rate = c(0.001,0.001),
                                                    q_conc = matrix(c(NA,1,1,NA),nrow=2,ncol=2)))

# Assign fixed prior distributions - also need behaviour priors if applicable
calc_prior_speed_lik <- function(move_params) {
  ifelse(move_params[4] > move_params[3], log(1), log(0))
}

# Assign variables relating to the number of MCMC iterations, thinning, etc.
section_length_limits <- c(3, 13)
speed_param_count <- 0; refined_path_count <- 0; behav_fail_count <- 0
num_iterations <- 10^9; num_extra_path_updates <- 100; thin <- 10^5
stored_move_params <- matrix(NA, nrow = num_iterations / thin, ncol = fixed_constant$num_move_params)
stored_refined_path <- vector("list", num_iterations / thin)
stored_behav_params <- vector("list", num_iterations / thin)

####### Carry out the full MCMC loop ##############################

# loop over the number of path updates
for (i in 1:num_iterations) {

  # but update parameters less often than the path
  if (i %% num_extra_path_updates == 0) {

    # update the behaviour parameters by a Gibbs sampler
    behav_params <- update_behav_param(fixed_constant, behav_proc)

    # propose an update to the speed parameters in a random walk MH step
    speed_param_proposal <- update_speed_param(fixed_constant, curr_move_params = move_params, behav_params, refined_path)
    # decide whether to accept the MH proposal
    if (runif(1) < speed_param_proposal$accept_prob) {

      speed_param_count <- speed_param_count + 1 # keep track of acceptance rate
      move_params[-c(fixed_constant$q.bearings, fixed_constant$q.obs_error)] <- speed_param_proposal$prop_speed_params
    }

    # update the bearing parameters by a Gibbs sampler
    move_params[fixed_constant$q.bearings] <- update_bearing_param(fixed_constant, refined_path)
  }

  # update the refined path (in sections, chosen randomly)
  # flag if it is a special end case and pick all relevent info out from the full path
  drawn_section <- draw_random_path_section(section_length_limits, obs_index = obs$index_on_refined_path)
  type <- if (drawn_section$start == 1) "start" else if (drawn_section$end == length(refined_path$bearings)) "end"  else  "middle"
  fixed_values <- set_fixed_values(type, drawn_section, refined_path, obs)
  curr_values <- set_current_values(type, drawn_section, refined_path)

  # update the path section
  refined_path_proposal <- update_refined_path(fixed_constant, type, fixed_values, curr_values, move_params, behav_params, obs)

  # update values in the full path if it was accepted
  if (refined_path_proposal$accept == 1) {

    refined_path_count <- refined_path_count + 1 # keep track of acceptance rate

    time_l <- length(refined_path$times)
    path_l <- time_l - 1
    start_to_pathend <- drawn_section$start:time_l
    one_to_end <- 1:drawn_section$end
    refined_path$times <- c(refined_path$times[-start_to_pathend], refined_path_proposal$times, refined_path$times[-one_to_end])
    refined_path$behavs <- c(refined_path$behavs[-start_to_pathend], refined_path_proposal$behavs, refined_path$behavs[-one_to_end])
    refined_path$bearings <- c(refined_path$bearings[-start_to_pathend], refined_path_proposal$bearings, refined_path$bearings[-one_to_end])
    refined_path$steps <- c(refined_path$steps[-start_to_pathend], refined_path_proposal$steps, refined_path$steps[-one_to_end])
    refined_path$time_diffs <- diff(refined_path$times)
    obs$index_on_refined_path <- which(refined_path$times %in% obs$times)
    obs$errors <- list(x = obs$x - refined_path$X[obs$index_on_refined_path], y = obs$y - refined_path$Y[obs$index_on_refined_path])
    behav_proc <- create_behav_proc_from_refined(times = refined_path$times, states = refined_path$behavs)

    if (type == "middle") {

      one_to_start <- 1:drawn_section$start
      one_to_end <- 1:(drawn_section$end + 1)
      refined_path$X <- c(refined_path$X[one_to_start], refined_path$X[drawn_section$start] + cumsum(refined_path_proposal$steps * cos(refined_path_proposal$bearings)), refined_path$X[-one_to_end])
      refined_path$Y <- c(refined_path$Y[one_to_start], refined_path$Y[drawn_section$start] + cumsum(refined_path_proposal$steps * sin(refined_path_proposal$bearings)), refined_path$Y[-one_to_end])
    } else if (type == "start") {

      one_to_end <- 1:(drawn_section$end)
      refined_path$X <- c(rev(refined_path$X[drawn_section$end + 1] - cumsum(rev(refined_path_proposal$steps * cos(refined_path_proposal$bearings)))), refined_path$X[-one_to_end])
      refined_path$Y <- c(rev(refined_path$Y[drawn_section$end + 1] - cumsum(rev(refined_path_proposal$steps * sin(refined_path_proposal$bearings)))), refined_path$Y[-one_to_end])
      # keep the first bearing between -pi,pi to "thether" the bearing process
      while(refined_path$bearings[1] > pi) {refined_path$bearings <- refined_path$bearings - 2 * pi}
      while(refined_path$bearings[1] < (-pi)) {refined_path$bearings <- refined_path$bearings + 2 * pi}
    } else if (type == "end") {

      one_to_start <- 1:drawn_section$start
      refined_path$X <- c(refined_path$X[one_to_start], refined_path$X[drawn_section$start] + cumsum(refined_path_proposal$steps * cos(refined_path_proposal$bearings)))
      refined_path$Y <- c(refined_path$Y[one_to_start], refined_path$Y[drawn_section$start] + cumsum(refined_path_proposal$steps * sin(refined_path_proposal$bearings)))
    }
  } else if (refined_path_proposal$behav_fail == 1) {

    behav_fail_count <- behav_fail_count + 1 # keep track of how the behavioural bridge rejection proposal is going
  }

  # Store a thinned sample of the full system
  if (i %% thin == 0) {

    stored_move_params[i / thin, ] <- move_params
    stored_refined_path[[i / thin]] <- refined_path
    stored_behav_params[[i / thin]] <- behav_params
  }
}
############################################

save.image("results.RData")
write.table(stored_move_params, file = "move_param.txt", row.names = F)
stored_lambda <- matrix(NA, nrow = 10000, ncol = 2)
for(i in 1:10000){
  stored_lambda[i, ] <- stored_behav_params[[i]]$lambda
}
write.table(stored_lambda, file = "behav_param.txt", row.names = F)
write.table(c(speed_param_count, speed_param_count / (num_iterations / num_extra_path_updates),
              refined_path_count, refined_path_count / num_iterations), file = "acceptance.txt", row.names = F)
n <- num_iterations / thin
l <- length(refined_path$X)
samp_bearings <- matrix(NA, nrow = 100, ncol = 1100)
samp_steps <- matrix(NA, nrow = 100, ncol = 1100)
samp_times <- matrix(NA, nrow = 100, ncol = 1100)
samp_behavs <- matrix(NA, nrow = 100, ncol = 1100)
samp_loc <- matrix(NA, nrow = 100*1100, ncol = 3)
for(i in 1:100){
  j <- seq(4060,10000,60)[i]
  samp_bearings[i, ] <- c(stored_refined_path[[j]]$bearings, rep(NA,1100-length(stored_refined_path[[j]]$bearings)))
  samp_steps[i, ] <- c(stored_refined_path[[j]]$steps, rep(NA,1100-length(stored_refined_path[[j]]$steps)))
  samp_behavs[i, ] <- c(stored_refined_path[[j]]$behavs, rep(NA,1100-length(stored_refined_path[[j]]$behavs)))
  samp_times[i, ] <- c(stored_refined_path[[j]]$times, rep(NA,1100-length(stored_refined_path[[j]]$times)))
  samp_loc[((i-1)*1100+1):(i*1100),1] <- c(stored_refined_path[[j]]$X,rep(NA,1100-length(stored_refined_path[[j]]$X)))
  samp_loc[((i-1)*1100+1):(i*1100),2] <- c(stored_refined_path[[j]]$Y,rep(NA,1100-length(stored_refined_path[[j]]$Y)))
  samp_loc[((i-1)*1100+1):(i*1100),3] <- i
}
write.table(samp_bearings, "bearings.txt", row.names = F)
write.table(samp_steps, "steps.txt", row.names = F)
write.table(samp_behavs, "behavs.txt", row.names = F)
write.table(samp_times, "times.txt", row.names = F)
write.table(samp_loc, "locs.txt", row.names = F)
