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

This package only provides the functions necessary to implement inference for continuous-time step-and-turn models. The full inference script will need to be created and run for a given instance. A number of example implementations are included that can be modified for use.

### Example single_reindeer_error

This example implements single state movement for a set of reindeer observations, with unknown, independent observation errors. To run this example, source the `single_reindeer_error/script.R` file. This file is broken up into a number of stages:
1) The observed locations are read and initial values for the inference algorithm are assigned. Fixed constants relating to the specific scenario are set (e.g. number of behavioural states, prior distributions, pertubation variances).
2) The MCMC algorithm is implemented. This involves a loop sampling parameters and refined path, and uses the functions from the package `CTStepTurn`.
3) The results of the algorithm are stored as text files.

Figures of the results of this example can be produced using the `single_reindeer_error/plot_results/plot.Rnw` file. This is an R Sweave document that can be compiled in R to produce individual pdf's of each of the results figures. 

## Authors

* **Alison Parton** - [personal webpage](http://alisonparton.co.uk/).

## License

This project is licensed under the GNU License - see the [LICENSE.md](LICENSE.md) file for details.


