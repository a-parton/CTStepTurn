num_states <- 2
fixed_values <- list(time = 10,
                     bearing = 2.3)
times <- seq(0, 10, length.out = 22)[-22]
behavs <- c(rep(1, 6), rep(2, 9), rep(1, 3), rep(2, 3))
variance <- c(0.3, 1.5)

set.seed(123)
sim_brownian_backward(num_states, fixed_values, times, behavs, variance)
# [1] 2.210662 1.998823 1.911824 2.500960 2.527610 2.576476 3.224709 3.614255 2.545083 1.964586
# [11] 1.587933 2.622471 2.926569 3.265283 3.358827 2.889055 3.564445 3.752615 3.009303 3.602057
# [21] 3.202476

num_states <- 2
fixed_values <- list(time = 10,
                     bearing = 2.3)
times <- seq(0, 10, length.out = 22)[-22]
behavs <- rep(1, 21)
variance <- c(0.03, 1.5)

set.seed(2712)
sim_brownian_backward(num_states, fixed_values, times, behavs, variance)
# [1] 2.439727 2.533936 2.586692 2.596769 2.267357 2.350119 2.591627 2.515462 2.531663 2.719289
# [11] 2.697561 2.572232 2.618390 2.619669 2.705751 2.734348 2.720683 2.616551 2.512080 2.555922
# [21] 2.424459
