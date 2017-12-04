#' @title Draw a random path section.
#'
#' @description \code{draw_random_path_section} returns details on a randomly selected section of path to update.
#'
#' @family Choose path
#'
#' @param section_length_limits Vector of length two, lower and upper length limits for the section to choose.
#' @param obs_index Vector of length number of observations, giving the index on the refined path where observation times coincide.
#'
#' @return List with components:
#'   \describe{
#'   \item{\code{start}}{Numeric, the index of the start of chosen section.}
#'   \item{\code{end}}{Numeric, the index of the end of chosen section.}
#'   \item{\code{contained_obs_index}}{Vector (length variable) of index which are also observations (within the section).}
#'   \item{\code{contained_obs}}{Vector (length variable) of observation indexs included in the section.}
#'   }
#'
#' @example examples/draw_random_path_section.R
#'
#' @export draw_random_path_section
draw_random_path_section <- function(section_length_limits, obs_index) {

  # how many points are on the refined path to choose from?
  total_length <- obs_index[length(obs_index)] - 1
  # choose an initial random starting point for section
  start <- sample(1:(total_length - section_length_limits[1]), 1)
  # choose an initial random length of section
  if (start == (total_length - section_length_limits[1])) {
    length <- section_length_limits[1]
  } else {
    length <- sample(section_length_limits[1]:min(section_length_limits[2], total_length - start), 1)
  }
  # which observations are contained within the chosen section?
  contained_obs_index <- obs_index[obs_index > start & obs_index <= (start + length)]
  contained_obs <- which(obs_index > start & obs_index <= (start + length))

  # if there are no observations inside the chosen section then we can end straight away
  if (length(contained_obs) == 0) {

    return(list(start = start, end = start + length, contained_obs_index = contained_obs_index, contained_obs = contained_obs))
  } else {

    # otherwise we have to deal with possible difficulties with where the observation points are

    # if there is an observation 1 or 2 points into the chosen section it causes problems due to dimensionality
    if(obs_index[contained_obs[1]] == (start + 1) || obs_index[contained_obs[1]] == (start + 2)) {

      # change the start of the section to be at the observation point instead and reduce the length to still end at the same point
      length <- length - obs_index[contained_obs[1]] + start
      start <- obs_index[contained_obs[1]]
      contained_obs_index <- obs_index[obs_index > start & obs_index <= (start + length)]
      contained_obs <- which(obs_index > start & obs_index <= (start + length))

      # if this leads to no observations being contained inside the section anymore, end
      if (length(contained_obs) == 0) {

        # but if this has reduced the section to being too small to realistically update over (because of dimensionality), we start again
        if(length < 3) return(draw_random_path_section(section_length_limits, obs_index))
        return(list(start = start, end = start + length, contained_obs_index = contained_obs_index, contained_obs = contained_obs))
      }
    }

    # the location at index start + length + 1 is fixed (because the bearing and step at index start + length get us to that location)
    # so if there is an observation at time start + length or start + length - 1 then we again have dimensionality issues
    if(obs_index[contained_obs[length(contained_obs)]] == (start + length) || obs_index[contained_obs[length(contained_obs)]] == (start + length - 1)) {

      # the start point stays the same, but the length is reduced so that the observation is at start + length + 1
      length <- obs_index[contained_obs[length(contained_obs)]] - start - 1
      contained_obs_index <- obs_index[obs_index > start & obs_index <= (start + length)]
      contained_obs <- which(obs_index > start & obs_index <= (start + length))
    }

    # path section is too short after the changes made, that we have to start again
    if(length < 3) return(draw_random_path_section(section_length_limits, obs_index = obs_index))
    return(list(start = start, end = start + length, contained_obs_index = contained_obs_index, contained_obs = contained_obs))
  }
}


#' @title Identify fixed values relating to a chosen path section.
#'
#' @description \code{set_fixed_values} returns details on the fixed values needed to update a path section.
#'
#' @family Choose path
#'
#' @param type Character, whether the path section is at the "start, "end", or "middle" of the full path.
#' @param drawn_section List, output from \code{draw_random_path_section}.
#' @param refined_path List with components:
#'   \describe{
#'   \item{\code{times}}{Vector, refined time scale.}
#'   \item{\code{behavs}}{Vector, refined behavioural/state sequence.}
#'   \item{\code{bearings}}{Vector, refined bearing process.}
#'   \item{\code{steps}}{Vector, refined step process.}
#'   \item{\code{X}}{Vector, refined x location.}
#'   \item{\code{Y}}{Vector, refined y location.}
#'   }
#' @param obs List with components:
#'   \describe{
#'   \item{\code{times}}{Vector, observation time scale.}
#'   \item{\code{x}}{Vector, observed x location.}
#'   \item{\code{y}}{Vector, observed y location.}
#'   }
#'
#' @return List of fixed values.
#'
#' @export set_fixed_values
set_fixed_values <- function(type, drawn_section, refined_path, obs) {

  start <- drawn_section$start; end <- drawn_section$end; contained_obs <- drawn_section$contained_obs

  # handle the three cases of path sections separately
  if (type == "middle") {

    index_errors <- c(max(which(obs$times <= refined_path$times[start])), min(which(obs$times >= refined_path$times[end])))
    list(inc_behavs = refined_path$behavs[c(start, end)], # null in single behav case
         inc_times = refined_path$times[c(start, end)],
         obs_times = obs$times[contained_obs],
         contained_obs = contained_obs,
         errors = list(x = obs$errors$x[index_errors],
                       y = obs$errors$y[index_errors]),
         error_times = obs$times[index_errors],
         times = refined_path$times[c(start - 1, end + 1)],
         behav = refined_path$behavs[start - 1], # null in single behav case
         bearings = refined_path$bearings[c(start - 1, end + 1)],
         steps = refined_path$steps[c(start - 1, end + 1)],
         after_time = refined_path$times[end + 2],
         locs = as.vector(rbind(c(obs$x[contained_obs], refined_path$X[end + 1]), c(obs$y[contained_obs], refined_path$Y[end + 1]))),
         start_loc = c(refined_path$X[start], refined_path$Y[start])
        )
  } else if (type == "start") {

    index_errors <- min(which(obs$times >= refined_path$times[end]))
    list(inc_behav = refined_path$behavs[end], # null in single behav case
         inc_times = refined_path$times[c(1, end)],
         obs_times = obs$times[contained_obs],
         contained_obs = contained_obs,
         errors = list(x = obs$errors$x[index_errors],
                       y = obs$errors$y[index_errors]),
         error_times = obs$times[index_errors],
         time = refined_path$times[end + 1],
         after_time = refined_path$times[end + 2],
         behav = refined_path$behavs[end], # null in single behav case
         bearing = refined_path$bearings[end + 1],
         step = refined_path$steps[end + 1],
         locs = as.vector(rbind(c(obs$x[1], obs$x[contained_obs]), c(obs$y[1], obs$y[contained_obs]))),
         start_loc = c(refined_path$X[end + 1], refined_path$Y[end + 1])
        )
  } else if (type == "end") {

    index_errors <- max(which(obs$times <= refined_path$times[start]))
    list(inc_behav = refined_path$behavs[start], # null in single behav case
         inc_times = refined_path$times[c(start, end)],
         obs_times = obs$times[contained_obs],
         contained_obs = contained_obs,
         errors = list(x = obs$errors$x[index_errors],
                       y = obs$errors$y[index_errors]),
         error_times = obs$times[index_errors],
         time = refined_path$times[start - 1],
         behav = refined_path$behavs[start - 1], # null in single behav case
         bearing = refined_path$bearings[start - 1],
         step = refined_path$steps[start - 1],
         after_time = refined_path$times[end + 1],
         locs = as.vector(rbind(c(obs$x[contained_obs], obs$x[length(obs$x)]), c(obs$y[contained_obs], obs$y[length(obs$y)]))),
         start_loc = c(refined_path$X[start], refined_path$Y[start])
        )
  }
}


#' @title Identify current values relating to a chosen path section.
#'
#' @description \code{set_current_values} returns details on the current values to be simulated.
#'
#' @family Choose path
#'
#' @param type Character, whether the path section is at the "start, "end", or "middle" of the full path.
#' @param drawn_section List, output from \code{draw_random_path_section}.
#' @param refined_path List with components:
#'   \describe{
#'   \item{\code{times}}{Vector, refined time scale.}
#'   \item{\code{behavs}}{Vector, refined behavioural/state sequence.}
#'   \item{\code{bearings}}{Vector, refined bearing process.}
#'   \item{\code{steps}}{Vector, refined step process.}
#'   \item{\code{X}}{Vector, refined x location.}
#'   \item{\code{Y}}{Vector, refined y location.}
#'   }
#'
#' @return List of current values.
#'
#' @export set_current_values
set_current_values <- function(type, drawn_section, refined_path) {

  start <- drawn_section$start; end <- drawn_section$end; contained_obs_index <- drawn_section$contained_obs_index

  # handle the start" case separately
  if (type == "start") {

    list(times = refined_path$times[1:end],
         behavs = refined_path$behavs[1:end], # null in single behav case
         bearings = refined_path$bearings[1:end],
         loc_index = c(1, contained_obs_index)
        )
  } else {

    list(times = refined_path$times[start:end],
         behavs = refined_path$behavs[start:end], # null in single behav case
         bearings = refined_path$bearings[start:end],
         loc_index = c(contained_obs_index - start, end - start + 1) # this is the index on the refined path that will "reach" an observation, i.e. the time directly before the obs
        )
  }
}


#' @title Create refined time scale.
#'
#' @description \code{create_refined_times} creates a set of refined times, based on an ideal time scale, including all observation and switch times.
#'
#' @family Refined scales
#'
#' @param ideal_refined_time_diff Numeric, ideal time scale.
#' @param behav_proc_times Vector (variable length), times of behavioural switches.
#' @param contained_obs_times Vector (variable length), times of observations.
#'
#' @return List with elements:
#'   \describe{
#'   \item{\code{refined_times}}{Vector of refined times.}
#'   \item{\code{loc_index}}{Index value within section when observations occur.}
#'   \item{\code{loc_index_start}}{For sections at the start of the path the loc index must be treated differently}
#'   }
#'
#' @example examples/create_refined_times.R
#'
#' @export create_refined_times

create_refined_times <- function(ideal_refined_time_diff, behav_proc_times, contained_obs_times) {

  # the refined times must at least include the start and end times (included in behav times),
  # the behavioural switching times and observation times
  base_times <- sort(c(behav_proc_times, contained_obs_times))

  # then it can be filled in with more times, that fit as best as possible to the desired time scale interval
  points_between <- ceiling(diff(base_times) / ideal_refined_time_diff) + 1
  refined_times <- base_times[1]
  for (i in 1:length(points_between)) {

    refined_times <- c(refined_times, seq(base_times[i], base_times[i + 1], length.out = points_between[i])[-1])
  }

  # But then considerations must be made to whether observation points will cause dimensionality issues
  # in the same way as when the section of path was originally chosen
  if (length(contained_obs_times) > 0) {

    # If the second point in the new refined timescale is an observation
    if(refined_times[2] == contained_obs_times[1]) {

      # Then we add some more points into the refined time scale to push this observation away from the start
      refined_times <- c(refined_times[1], seq(refined_times[1], contained_obs_times[1], length.out = 4)[2:3], refined_times[-1])

    # If the third point in the new refined timescale is an observation
    } else if(refined_times[3] == contained_obs_times[1]) {

      refined_times <- c(refined_times[1], seq(refined_times[1], contained_obs_times[1], length.out = 4)[2:3], refined_times[-(1:2)])
    }

    # If the second to last refined time is an observation
    if(refined_times[length(refined_times) - 1] == contained_obs_times[length(contained_obs_times)]) {

      refined_times <- c(refined_times[-length(refined_times)], seq(contained_obs_times[length(contained_obs_times)], refined_times[length(refined_times)], length.out = 4)[2:3], refined_times[length(refined_times)])
    }
    # We don't need to check whether the last refined time is an observation time, because we would have already not allowed this to happen when we chose the section
  }

  # Record where observation times are within this section of times
  # These are the index on the section that gets you to an observed location, so are the times directly before an observation time
  loc_index <- c(which(refined_times %in% contained_obs_times) - 1, length(refined_times))
  loc_index_start <- c(1, which(refined_times %in% contained_obs_times))

  list(refined_times = refined_times, loc_index = loc_index, loc_index_start = loc_index_start)
}


#' @title Create refined behaviours.
#'
#' @description \code{create_refined_behavs} creates behavioural process at the refined time scale given info about switch times and states.
#'
#' @family Refined scales
#'
#' @param behav_proc List with elements:
#'   \describe{
#'   \item{\code{times}}{Vector of switch times.}
#'   \item{\code{states}}{Vector of switch states.}
#'   }
#' @param refined_times Vector, refined time scale.
#'
#' @return Vector, refined behaviours
#'
#' @example examples/create_refined_behavs.R
#'
#' @export create_refined_behavs
create_refined_behavs <- function(behav_proc, refined_times) {

  # identify where on the refined scale there are behavioural switches
  which_switches <- which(refined_times %in% behav_proc$times)
  refined_behavs <- rep(NA, length(refined_times))

  # Then fill in all points between switches with the state at that switch
  for(i in 1:(length(which_switches) - 1)) {
    refined_behavs[which_switches[i]:which_switches[i + 1]] <- behav_proc$states[i]
  }

  return(refined_behavs)
}


#' @title Create behavioural switch times and states.
#'
#' @description \code{create_behav_proc_from_refined} creates behavioural process as switch times and states given refined behaviours.
#'
#' @family Refined scales
#'
#' @param times Vector of switch times.
#' @param states Vector of switch states.
#'
#' @return List with the components:
#'   \describe{
#'   \item{\code{times}}{Vector, switching times of the simulation.}
#'   \item{\code{states}}{Vector, states at switching times of the simulation.}
#'   }
#'
#' @example examples/create_behav_proc_from_refined.R
#'
#' @export create_behav_proc_from_refined
create_behav_proc_from_refined <- function(times, states) {

  switch_times <- times[1]
  switch_states <- states[1]

  for (i in 2:length(states)) {
    if (states[i] != states[i - 1]) {
      switch_times <- c(switch_times, times[i])
      switch_states <- c(switch_states, states[i])
    }
  }

  switch_times <- c(switch_times, times[length(times)])
  switch_states <- c(switch_states, states[length(states)])

  list(times = switch_times, states = switch_states)
}
