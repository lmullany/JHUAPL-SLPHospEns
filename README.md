## Pre Smoothed Simulated Linear Pooling for COVID-19 Hospitalization Ensembles

This repo provides two R packages; see individual README.md files for more information. Briefly,

-   psSLP: This is a general use package that provides functions aiding in combining quantile functions or modeling/forecasts outputs. Includes:

    -   functions to pre-smoothing longitudinal sets of quantile estimations

    -   functions for generating both simulated linear pooling and vincentization ensemble estimates

    -   is intended for pre-smoothing and simulating tens of thousands of forecasts, and by default is
        parallelized

    -   Installation: devtools::install_github("lmullany/psSLP.git)

-   covid19hosp: This is a context/task -specific, non-general wrapper package for preparing data inputs to **psSLP** from data available at the Covid 19 Forecasting Hub (<https://github.com/reichlab/covid19-forecast-hub>).

    -   Requires local access to the hub data, preferably through a clone of the repo

    -   Provides functions for pulling modeling group specific quantile functions

    -   Includes checks to examine data for (longitudinal) consistencies

    -   Built-in functions for pulling empirical data from healthdata.gov

    -   Aids in estimating coverage of ensembles

    -   Includes plotting functions to aid in displaying results
