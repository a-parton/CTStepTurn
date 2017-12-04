const <- list(num_states = 3,
              indep_step = FALSE,
              q.speed_mean = 4:6,
              q.speed_corr = 7:9,
              q.speed_var = 10:12)
path <- list(steps = c(5.3, 6.8, 15.3, 14.5, 9.4, 3.4, 8.5, 14.6, 17.8, 22.4, 16.7),
             time_diffs = c(1, 1, 1.5, 1.5, 1, 0.5, 1, 1, 1, 2, 1),
             behavs = c(1, 1, 1, 2, 2, 3, 3, 2, 2, 2, 1))
move_params <- c(NA, NA, NA, 5, 10, 3, 0.1, 0.2, 0.3, 4, 6, 2)
behav_params <- list(lambda = c(0.1, 0.2, 0.3),
                     q = matrix(c(NA, 0.5, 0.2, 0.3, NA, 0.8, 0.7, 0.5, NA), nrow = 3))

calc_step_lik(const, path, move_params, behav_params)
# [1] -78.24704

path <- list(steps = c(5.3, 6.8, 15.3, 14.5, 9.4, 3.4, 8.5, 14.6, 17.8, 22.4, 16.7),
             time_diffs = c(2, 2, 0.5, 1, 2, 0.25, 2, 1.5, 1, 1.2, 0.6),
             behavs = c(1, 1, 1, 2, 2, 3, 3, 2, 2, 2, 1))
calc_step_lik(const, path, move_params, behav_params)
# [1] -791.1856
