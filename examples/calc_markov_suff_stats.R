times <- c(0, 12, 15, 20, 26, 43, 50)
states <- c(1, 2, 1, 2, 1, 2, 2)

calc_markov_suff_stats(2, times, states)
# $time_in_state
# [1] 34 16
#
# $num_trans
# [,1] [,2]
# [1,]   NA    3
# [2,]    2   NA

times <- c(100, 140, 167, 189, 241, 264, 279, 304, 365, 397, 400)
states <- c(1, 2, 1, 3, 2, 3, 1, 2, 3, 2, 2)
calc_markov_suff_stats(3, times, states)
# time_in_state
# [1]  87 114  99
#
# $num_trans
# [,1] [,2] [,3]
# [1,]   NA    2    1
# [2,]    1   NA    2
# [3,]    1    2   NA
