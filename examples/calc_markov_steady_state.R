para <- list(lambda = c(0.1, 0.2),
             q = matrix(c(NA, 1, 1, NA), nrow = 2))
calc_markov_steady_state(2, para)
# $steady_state
# [1] 0.6666667 0.3333333
#
# $generator
# [,1] [,2]
# [1,] -0.1  0.1
# [2,]  0.2 -0.2

para <- list(lambda = c(0.1, 0.2, 0.3),
             q = matrix(c(NA, 0.5, 0.2, 0.3, NA, 0.8, 0.7, 0.5, NA), nrow = 3))
calc_markov_steady_state(3, para)
$steady_state
# [1] 0.4568528 0.3274112 0.2157360
#
# $generator
# [,1]  [,2]  [,3]
# [1,] -0.10  0.03  0.07
# [2,]  0.10 -0.20  0.10
# [3,]  0.06  0.24 -0.30
