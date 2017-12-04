prior <- list(lambda_rate = c(0.1, 0.2),
              lambda_shape = c(2, 1),
              q_conc = matrix(c(NA, 1, 1, NA), nrow = 2))
suff <- list(time_in_state = c(34, 16),
             num_trans = matrix(c(NA, 2, 3, NA), nrow = 2))

calc_markov_post_param(prior, 2, suff)
# $lambda_shape
# [1] 5 3
#
# $lambda_rate
# [1] 34.1 16.2
#
# $q_conc
# [,1] [,2]
# [1,]   NA    1
# [2,]    1   NA

prior <- list(lambda_rate = c(0.1, 0.2, 0.3),
              lambda_shape = c(4, 3, 2),
              q_conc = matrix(c(NA, 0.5, 0.3, 0.2, NA, 0.7, 0.8, 0.5, NA), nrow = 3))
suff <- list(time_in_state = c(87, 114, 99),
             num_trans = matrix(c(NA, 1, 1, 2, NA, 2, 1, 2, NA), nrow = 3))

calc_markov_post_param(prior, 3, suff)
# $lambda_shape
# [1] 7 6 5
#
# $lambda_rate
# [1]  87.1 114.2  99.3
#
# $q_conc
# [,1] [,2] [,3]
# [1,]   NA  2.2  1.8
# [2,]  1.5   NA  2.5
# [3,]  1.3  2.7   NA
