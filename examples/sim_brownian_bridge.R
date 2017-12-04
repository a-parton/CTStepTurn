num_states <- 2
fixed_values <- list(times = c(0, 10),
                     behav = 1,
                     bearings = c(0.1, 2.3))
times <- seq(0, 10, length.out = 22)[-c(1, 22)]
behavs <- c(rep(1, 5), rep(2, 9), rep(1, 3), rep(2, 3))
variance <- c(0.3, 1.5)

set.seed(123)
sim_brownian_bridge(num_states, fixed_values, times, behavs, variance)
# [1] -0.08125057 -0.13766018  0.48206550  0.53930448  0.61875996  1.29758290  1.84007475
# [8]  0.92384942  0.49629936  0.27259280  1.46007729  1.91712203  2.40878228  2.65527280
# [15]  2.33844785  3.04442684  3.26318595  2.55046384  3.29616431  3.04952920

num_states <- 2
fixed_values <- list(times = c(0, 10),
                     behav = 1,
                     bearings = c(0.1, 2.3))
times <- seq(0, 10, length.out = 22)[-c(1, 22)]
behavs <- rep(1, 20)
variance <- c(0.03, 1.5)

set.seed(2712)
sim_brownian_bridge(num_states, fixed_values, times, behavs, variance)
# [1] 0.3056246 0.4697955 0.5912883 0.3732920 0.5674697 0.9203932 0.9556434 1.0832605 1.3823015
# [10] 1.4719894 1.4580763 1.6156492 1.7283444 1.9258414 2.0658545 2.1636051 2.1708882 2.1778333
# [19] 2.3330903 2.3130439
