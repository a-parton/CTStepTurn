const <- list(prior_param = list(lambda_rate = c(0.1, 0.2),
                                 lambda_shape = c(2, 1),
                                 q_conc = matrix(c(NA, 1, 1, NA), nrow = 2)),
              num_states = 2)
behav_proc <- list(times = c(0, 12, 15, 20, 26, 43, 50),
                   states = c(1, 2, 1, 2, 1, 2, 2))

set.seed(1234)
update_behav_param(const, behav_proc)
# $lambda
# [1] 0.06755656 0.18652884
#
# $q
# [,1] [,2]
# [1,]   NA    1
# [2,]    1   NA

const <- list(prior_param = list(lambda_rate = c(0.1, 0.2, 0.3),
                                 lambda_shape = c(4, 3, 2),
                                 q_conc = matrix(c(NA, 0.5, 0.3, 0.2, NA, 0.7, 0.8, 0.5, NA), nrow = 3)),
              num_states = 3)
behav_proc <- list(times = c(100, 140, 167, 189, 241, 264, 279, 304, 365, 397, 400),
                   states = c(1, 2, 1, 3, 2, 3, 1, 2, 3, 2, 2))

set.seed(7654)
update_behav_param(const, behav_proc)
# $lambda
# [1] 0.04658404 0.05334004 0.08813678
#
# $q
# [,1]      [,2]      [,3]
# [1,]          NA 0.6686524 0.3313476
# [2,] 0.008780271        NA 0.9912197
# [3,] 0.366922353 0.6330776        NA
