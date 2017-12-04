# In a two state process, the reverse chain is equal to the forward chain
num_states <- 2
behav_params <- list(lambda = c(0.1,0.2),
                     q = matrix(c(NA,1,1,NA), nrow = 2))

calc_reverse_markov_gen(num_states, behav_params)
# $lambda
# [1] 0.1 0.2
# 
# $q
# [,1] [,2]
# [1,]   NA    1
# [2,]    1   NA


num_states <- 3
behav_params <- list(lambda = c(0.1,0.2,0.3),
                     q = matrix(c(NA,0.15/0.2,0.05/0.3,0.05/0.1,NA,0.25/0.3,0.05/0.1,0.05/0.2,NA),nrow=3))

calc_reverse_markov_gen(num_states, behav_params)
# $lambda
# [1] 0.1 0.2 0.3
# 
# $q
# [,1]      [,2]      [,3]
# [1,]        NA 0.8684211 0.1315789
# [2,] 0.4318182        NA 0.5681818
# [3,] 0.6333333 0.3666667        NA