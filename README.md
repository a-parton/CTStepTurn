# CTStepTurn

CTStepTurn is an R package containing the functions needed to carry out inference for multistate continuous-step step-and-turn movement modelling, as detailed in Parton and Blackwell (2017), available at https://link.springer.com/article/10.1007/s13253-017-0286-5. 

A number of example implementation scripts are provided, detailed below.

For more information, please contact aparton2@sheffield.ac.uk.

## Installation
Using the `devtools' package:

```
> install.packages("devtools")
> library(devtools)
> devtools::install_github("a-parton/CTStepTurn")
```

### Prerequisites

A number of existing R packages are also required:

```
> install.packages("gtools", "msm", "mvnfast")
```

In some of the example implementations, the additional package is also required:

```
> install.packages("fdrtool")
```

## Usage

This package only provides the functions necessary to implement inference for continuous-time step-and-turn models. The full inference script will need to be created and run for a given instance. A number of example implementations are included that can be modified for use. In all cases, the necessary script involves the following elements:
1) The observed locations are read and initial values for the inference algorithm are assigned. Fixed constants relating to the specific scenario are set (e.g. number of behavioural states, prior distributions, pertubation variances).
2) The MCMC algorithm is implemented. This involves a loop sampling parameters and refined path, and uses the functions from the package `CTStepTurn`.
3) The results of the algorithm are stored as text files.

### Example single_simulation

This example implements single state movement for a set of simulated observations assuming the independent step model (all other examples below use correlated steps). Four implementations are included that use a different refined time scale. Refined time scales at the simulation time (0.01) and observation time (5) only implement parameter updates because the refined path is known. Refined time scales of 1 and 0.5 implement the full algorithm. To run this example, source the `single_simulation/script_sim.R`, `single_simulation/script_half.R`, `single_simulation/script_one.R`, `single_simulation/script_obs.R` files. 

1) In this case, the observations are given by `single_simulation/observations.txt`. The initial refined path is created using `InitialPath.r`. For this example, the independent step model is assumed (`indep_step = TRUE`) and no observation error is assumed (`obs_error = FALSE`, `corr_obs_error = FALSE`). Initial parameters are set as estimates from the initial path. Perturbation variances and prior distributions are assigned using `FixedConstants.R`. The speed parameter prior is set as a function `calc_prior_speed_lik`. Variables detailing the run length of the MCMC algorithm are set.

2) The MCMC algorithm is carried out as a for loop over the pre-set run length. Parameters are updated using 
```
update_speed_param()
update_bearing_param()
```
Note that `Gibbs_update_obs_error_param()` is not implemented as no error is assumed here.
Refined path (only in the `half` and `one` cases) is updated using
```
draw_random_path_section()
set_fixed_values()
set_current_values()
update_refined_path()
```
Samples of the parameters and path are stored.

3) The results of the algorithm are written to text files to store.

Figures of the results of this example can be produced using the `single_simulation/plot_results/plot.Rnw` file. This is an R Sweave document that can be compiled in R to produce individual pdf's of each of the results figures. Note that this relies on having implemented all included scripts to plot and compare the results.

### Example single_reindeer_error

This example implements single state movement for a set of reindeer observations, with unknown, independent observation errors. To run this example, source the `single_reindeer_error/script.R` file. 

1) In this case, the observations are given by `single_reindeer_error/observations.txt`. The initial refined path is created using `InitialPath.r` at a time scale of 0.25. For this example, the correlated step model is assumed (`indep_step = FALSE`) and independent observation error is assumed (`obs_error = TRUE`, `corr_obs_error = FALSE`). Initial parameters are set as estimates from the initial path. Perturbation variances and prior distributions are assigned using `FixedConstants.R`. The speed parameter prior is set as a function `calc_prior_speed_lik`. Variables detailing the run length of the MCMC algorithm are set.

2) The MCMC algorithm is carried out as a for loop over the pre-set run length. Parameters are updated using 
```
update_speed_param()
update_bearing_param()
Gibbs_update_obs_error_param()
```
Refined path is updated using
```
draw_random_path_section()
set_fixed_values()
set_current_values()
update_refined_path()
```
Samples of the parameters and path are stored.

3) The results of the algorithm are written to text files to store.

Figures of the results of this example can be produced using the `single_reindeer_error/plot_results/plot.Rnw` file. This is an R Sweave document that can be compiled in R to produce individual pdf's of each of the results figures. 

### Example single_simulation_error

This example implements single state movement for a set of simulated observations. The observations were simulated with independent errors, with variance 25. Scripts are provided to carry out a number of implementation in which the error variance parameter is assumed known and fixed for this set of observations (ranging from 0-50). To run this example, source the `single_simulation_error/script_X.R` file, where `X` is the fixed error variance. 

1) In this case, the observations are given by `single_simulation_error/noisy_observations.txt`. The initial refined path is created using `InitialPath.r` at a time scale of 0.25. For this example, the correlated step model is assumed (`indep_step = FALSE`) and independent observation error is assumed (`obs_error = TRUE`, `corr_obs_error = FALSE`). Initial parameters are set as estimates from the initial path. Perturbation variances and prior distributions are assigned using `FixedConstants.R`. The speed parameter prior is set as a function `calc_prior_speed_lik`. Variables detailing the run length of the MCMC algorithm are set.

2) The MCMC algorithm is carried out as a for loop over the pre-set run length. Parameters are updated using 
```
update_speed_param()
update_bearing_param()
```
Note that `Gibbs_update_obs_error_param()` is not implemented in this example because the error variance parameter is assumed known and fixed, so is not updated.
Refined path is updated using
```
draw_random_path_section()
set_fixed_values()
set_current_values()
update_refined_path()
```
Samples of the parameters and path are stored.

3) The results of the algorithm are written to text files to store.

Figures of the results of this example can be produced using the `single_simulation_error/plot_results/plot.Rnw` file. This is an R Sweave document that can be compiled in R to produce individual pdf's of each of the results figures. Note that this relies on having implemented all included scripts to plot and compare the results.

## Authors

* **Alison Parton** - [personal webpage](http://alisonparton.co.uk/).

## License

This project is licensed under the GNU License - see the [LICENSE.md](LICENSE.md) file for details.


