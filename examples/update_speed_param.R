const <- list(num_states = 3,
              indep_step = FALSE,
              q.bearings = 1:3,
              q.speed_mean = 4:6,
              q.speed_corr = 7:9,
              q.speed_var = 10:12,
              q.obs_error = 13,
              pert_sd = c(0.5, 1, 0.3, 0.01, 0.02, 0.03, 0.4, 0.6, 0.2),
              num_speed_param = 9)
path <- list(steps = c(5.3, 6.8, 15.3, 14.5, 9.4, 3.4, 8.5, 14.6, 17.8, 22.4, 16.7),
             time_diffs = c(1, 1, 1.5, 1.5, 1, 0.5, 1, 1, 1, 2, 1),
             behavs = c(1, 1, 1, 2, 2, 3, 3, 2, 2, 2, 1))
curr_move_params <- c(NA, NA, NA, 5, 10, 3, 0.1, 0.2, 0.3, 4, 6, 2, NA)
behav_params <- list(lambda = c(0.1, 0.2, 0.3),
                     q = matrix(c(NA, 0.5, 0.2, 0.3, NA, 0.8, 0.7, 0.5, NA), nrow = 3))
calc_prior_speed_lik <- function(move_params) {0}

set.seed(1234)
update_speed_param(const, curr_move_params, behav_params, path)
# $accept_prob
# [1] 0.01626418
#
# $prop_speed_params
# [1]  4.39646713 10.27742924  3.32533235  0.07654302  0.20858249  0.31518168  3.77010402  5.67202089  1.88710960

set.seed(543)
update_speed_param(const, curr_move_params, behav_params, path)
# $accept_prob
# [1] 138.744
#
# $prop_speed_params
# [1]  5.67570559 10.18547949  3.12945796  0.09809393  0.18056898  0.32304201  3.92861090  5.99740665  2.33194589
