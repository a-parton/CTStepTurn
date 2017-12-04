const <- list(num_states = 3,
              turn_prec_prior = list(rate = c(1, 2, 3),
                                     shape = c(0.1, 0.2, 0.3)))
path <- list(times = 1:20,
             behavs = c(1, 1, 1, 1, 2, 2, 2, 3, 1, 1, 3, 3, 3, 2, 2, 2, 2, 2, 2),
             bearings = cumsum(rnorm(19, 0, 1)))

set.seed(1234)
update_bearing_param(const, path)
# [1] 4.8920033 0.8789057 1.7115489

set.seed(765)
update_bearing_param(const, path)
# [1] 2.0455226 0.7477887 2.6614361
