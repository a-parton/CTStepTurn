end_points <- list(inc_times = c(0,10),
                   inc_behavs = c(1,3))
param <- list(lambda = c(0.1,0.2,0.3),
              q = matrix(c(NA,0.5,0.3,0.35,NA,0.7,0.65,0.5,NA),nrow=3))

set.seed(123) # simulation does not meet correct end point and is rejected
sim_markov_bridge(3, end_points, param)
# $accept
# [1] FALSE

set.seed(2712) # simultion is accepted
sim_markov_bridge(3, end_points, param)
# $accept
# [1] TRUE
#
# $times
# [1]  0.000000  5.694254  6.831134  7.503054 10.000000
#
# $states
# [1] 1 3 1 3 3
