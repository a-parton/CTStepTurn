end_points <- list(inc_times = c(0,10),
                   inc_behav = c(3))
param <- list(lambda = c(0.1,0.2,0.3),
              q = matrix(c(NA,0.5,0.3,0.35,NA,0.7,0.65,0.5,NA),nrow=3))

set.seed(123)
sim_markov_backward(3, end_points, param)
# $times
# [1]  0.000000  7.188476 10.000000
#
# $states
# [1] 1 3 3

set.seed(2712)
sim_markov_backward(3, end_points, param)
# $times
# [1]  0.000000  4.467302  4.691276  8.101915 10.000000
#
# $states
# [1] 2 3 1 3 3
