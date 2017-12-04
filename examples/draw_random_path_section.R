limits <- c(3, 10)
obs_index <- c(0,5,10,18,24,28,40,45,51,56,70)

set.seed(1234)
draw_random_path_section(limits, obs_index)
# $start
# [1] 10
#
# $end
# [1] 15
#
# $contained_obs_index
# numeric(0)
#
# $contained_obs
# integer(0)

set.seed(7654)
draw_random_path_section(limits, obs_index)
# $start
# [1] 46
#
# $end
# [1] 54
#
# $contained_obs_index
# [1] 51
#
# $contained_obs
# [1] 9

limits <- c(8, 15)
set.seed(678)
draw_random_path_section(limits, obs_index)
# $start
# [1] 14
#
# $end
# [1] 27
#
# $contained_obs_index
# [1] 18 24
#
# $contained_obs
# [1] 4 5
