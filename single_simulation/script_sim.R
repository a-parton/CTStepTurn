# Set a seed
set.seed(1234)
# Packages required
library(gtools); library(msm); library(mvnfast); #library(fdrtool)
library(CTStepTurn)
source("../FixedConstants.R"); source("../InitialPath.R")

# Set observations (x, y, times)
obs <- as.list(read.table("observations.txt", header = T))
obs$time_diffs <- diff(obs$times)

# Set initial refined path
refined_path <- as.list(read.table("simulated_path.txt", header = T))
refined_path$bearings <- refined_path$bearings[-length(refined_path$bearings)]
refined_path$steps <- refined_path$steps[-length(refined_path$steps)]
refined_path$time_diffs <- diff(refined_path$times)
obs$index_on_refined_path <- which(refined_path$times %in% obs$times)

# Set initial parameters - also need behaviour parameters if applicable
move_params <-  c(var(diff(refined_path$bearings))/0.1,
                  mean(refined_path$steps)/0.1,
                  var(refined_path$steps)/0.1)

s# Assign the fixed constants
fixed_constant <- FixedConstants(pert_sd = sqrt(c(0.02, 0.5)), error_pert_sd = NULL,
                                 num_states = 1,
                                 indep_step = TRUE, obs_error = FALSE, corr_obs_error = FALSE)

# Assign fixed prior distributions - also need behaviour priors if applicable
calc_prior_speed_lik <- function(move_params) {
  0
}

# Assign variables relating to the number of MCMC iterations, thinning, etc.
speed_param_count <- 0
num_iterations <- 10^5; thin <- 100
stored_move_params <- matrix(NA, nrow = num_iterations / thin, ncol = fixed_constant$num_move_params)

####### Carry out the full MCMC loop ##############################

# loop over the number of path updates
for (i in 1:num_iterations) {

  speed_param_proposal <- update_speed_param(fixed_constant, curr_move_params = move_params, behav_params = NULL, refined_path)
  # decide whether to accept the MH proposal
  if (runif(1) < speed_param_proposal$accept_prob) {

    speed_param_count <- speed_param_count + 1 # keep track of acceptance rate
    move_params[-c(fixed_constant$q.bearings, fixed_constant$q.obs_error)] <- speed_param_proposal$prop_speed_params
  }

  # update the bearing parameters by a Gibbs sampler
  move_params[fixed_constant$q.bearings] <- update_bearing_param(fixed_constant, refined_path)

  # Store a thinned sample of the full system
  if (i %% thin == 0) {

    stored_move_params[i / thin, ] <- move_params
  }
}
############################################

save.image("results_sim.RData")
write.table(stored_move_params, file = "move_param_sim.txt", row.names = F)
write.table(c(speed_param_count, speed_param_count / (num_iterations)), file = "acceptance_sim.txt", row.names = F)
